"""
Last updated 10/21/21 by Akansha Gupta

Follows MIT License:
Copyright (c) 2012-2021 Scott Chacon and others

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import re
import pandas as pd
from openpyxl import load_workbook
import itertools
import sys

clusters = {}
origLengths = {}
UTSLengths = {}
all_start = {}
all_end = {}
ranges = {}
min_start = {}
max_end = {}
GTF_length = {}
multiple_exons = []
exon_intervals = {}
orig_ranges = {}
UTS_min_start = {}
UTS_min_end = {}
splitMergedNames = {}
mergedNameClusters = {}

"""
Setup via Pandas Dataframes
"""
def setup(dataframe, mean_fragment_length = 200):
    dataframe = dataframe.drop(columns = ['seqname', 'source','score','frame'])
    dataframe['attribute'] = dataframe['attribute'].str.split(';').str[0]
    dataframe['attribute'] = dataframe['attribute'].str.split(' ').str[1]
    dataframe = dataframe[dataframe['feature'].str.contains("exon")]
    dataframe = dataframe.drop(columns=['feature'])
    strand_dfs = [x for _, x in dataframe.groupby('strand')]

    # sort and parse first strand of info
    strand1_df = strand_dfs[0]
    strand1_gene_dfs = strand1_df.sort_values('start', kind = "mergesort").groupby('attribute', sort = False)
    parseGroups(strand1_gene_dfs)

    # sort and parse second strand of info (in case both forward and reverse strands are used)
    try:
        strand2_df = strand_dfs[1]
        strand2_gene_dfs = strand2_df.sort_values('start', kind = "mergesort").groupby('attribute', sort = False)
        parseGroups(strand2_gene_dfs)
    except IndexError:
        pass

    orig_lengths(mean_fragment_length)
"""
Parse forward and reverse strand data for performing operations
"""
def parseGroups(proc_group):
    transcripts = []
    global all_start
    global all_end
    global min_start
    global max_end
    global GTF_length
    global ranges
    for group in proc_group.groups:
        curGroup = proc_group.get_group(group)
        exonNumber = 2
        for index, row in curGroup.iterrows():
            gene = row['attribute'].strip('\"')
            if gene not in all_start:
                exonNumber = 2
                all_start[gene] = row['start']
                all_end[gene] = row['end']
                min_start[gene] = row['start']
                max_end[gene] = row['end']
                GTF_length[gene] = row['end'] - row['start'] + 1
                ranges[gene] = [[row['start'], row['end']]]
                transcripts.append(gene)
            else:
                if gene not in multiple_exons:
                    multiple_exons.append(gene)
                all_start[gene + " " + str(exonNumber)] = row['start']
                all_end[gene + " " + str(exonNumber)] = row['end']

                # Isoform handling
                # Complete overlap
                if row['start'] in range(min_start[gene], max_end[gene] + 1) and row['end'] in range(min_start[gene], max_end[gene] + 1):
                    continue
                # Partial overlap
                elif row['start'] in range(min_start[gene], max_end[gene] + 1):
                    ranges[gene].append([max_end[gene], row['end']])
                    min_start[gene] = min(max_start[gene], row['start'])
                    max_end[gene] = max(max_end[gene], row['end'])
                    GTF_length[gene] += (row['end'] - row['start'] + 1) - (max_end[gene] - row['start'] + 1)
                # No overlap
                else:
                    ranges[gene].append([row['start'], row['end']])
                    min_start[gene] = min(min_start[gene], row['start'])
                    max_end[gene] = max(max_end[gene], row['end'])
                    GTF_length[gene] += row['end'] - row['start'] + 1
                exonNumber += 1
    global orig_ranges
    orig_ranges = ranges.copy()
    define_clusters(transcripts)

"""
Define clusters
"""
def define_clusters(transcripts):
    global clusters
    if len(clusters) == 0:
        for gene in transcripts:
            key = "+" + str(max_end[gene])
            if key not in clusters.keys():
                clusters[key] = [gene]
            else:
                clusters[key].append(gene)
    else:
        for gene in transcripts:
            key = "-" + str(min_start[gene])
            if key not in clusters.keys():
                clusters[key] = [gene]
            else:
                clusters[key].append(gene)
    return clusters

"""Intersections
Source: https://stackoverflow.com/questions/40367461/intersection-of-two-lists-of-ranges-in-python/40368603 
"""
def intersections(a,b):
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]
        if a_right < b_right:
            i += 1
        else:
            j += 1
        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)
    ri = 0
    while ri < len(ranges)-1:
        if ranges[ri][1] == ranges[ri+1][0]:
            ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]
        ri += 1
    return ranges

"""
Original lengths - calculate within clusters
"""
def orig_lengths(mean_fragment_length = 200):
    pastTranscript = None
    global GTF_length
    for clus in clusters:
        if len(clusters[clus]) == 1:
            origLengths[clusters[clus][0]] = GTF_length[clusters[clus][0]]
        else:
            # Check if any transcript in cluster has multiple exons
            multi = [True for transcript in clusters[clus] if transcript in multiple_exons]
            # No multi-exon genes
            if len(multi) == 0:
                if clus[0] == "-":
                    clusters[clus].reverse()
                for transcript in clusters[clus]:
                    if transcript == clusters[clus][0]:
                        pastTranscript = transcript
                        continue
                    elif transcript == clusters[clus][-1]:
                        origLengths[pastTranscript] = GTF_length[pastTranscript] - GTF_length[transcript]
                        origLengths[transcript] = GTF_length[transcript]
                    else:
                        origLengths[pastTranscript] = GTF_length[pastTranscript] - GTF_length[transcript]
                        pastTranscript = transcript
            # Multi-exon genes
            else:
                if clus[0] == "-":
                    clusters[clus].reverse()
                if len(clusters[clus]) == 2:
                    overlaps = intersections(ranges[clusters[clus][0]], ranges[clusters[clus][1]])
                    overlapLength = 0
                    origLengths[clusters[clus][0]] = GTF_length[clusters[clus][0]]
                    origLengths[clusters[clus][1]] = GTF_length[clusters[clus][1]]
                    # One gene completely encapsulated as overlap region
                    if overlaps == ranges[clusters[clus][0]]:
                        # Other gene truncated to remove overlap
                        for ar in overlaps:
                            overlapLength += ar[1] - ar[0] + 1
                        origLengths[clusters[clus][1]] -= overlapLength
                    elif overlaps == ranges[clusters[clus][1]]:
                        for ar in overlaps:
                            overlapLength += ar[1] - ar[0] + 1
                        origLengths[clusters[clus][0]] -= overlapLength
                    else:
                        for ar in overlaps:
                            # Overlap is a shared exon between two clusters --> not unique --> remove entirely
                            if ar in ranges[clusters[clus][0]] and ar in ranges[clusters[clus][1]]:
                                origLengths[clusters[clus][0]] -= (ar[1] - ar[0] + 1)
                                origLengths[clusters[clus][1]] -= (ar[1] - ar[0] + 1)
                                overlaps.remove(ar)
                        # Otherwise, overlap is between exons but only segment of them
                        # If the transcript has any unique ranges of its own other than these overlaps, get rid of overlaps from there
                        for ar in overlaps:
                            overlapLength += ar[1] - ar[0] + 1
                        origLengths[clusters[clus][0]] = GTF_length[clusters[clus][0]] - overlapLength
                        origLengths[clusters[clus][1]] = GTF_length[clusters[clus][1]] - overlapLength
    UTS_lengths(mean_fragment_length)
    return clusters

"""
Truncate genes 
"""
def truncate_genes(clusters):
    truncated_lengths = {}

    # for every combination of clusters (other than comparing with self)
    for clus, compclus in itertools.combinations(clusters, 2):
        if clus[0] == compclus[0]:
            for transcript in itertools.product(clusters[clus], clusters[compclus]):
                # Overlap detected
                if min_start[transcript[1]] in range(min_start[transcript[0]], max_end[transcript[0]]):
                    overlaps = intersections(ranges[transcript[0]], ranges[transcript[1]])
                    res = []
                    [res.append(x) for x in overlaps if x not in res]
                    if len(res) == 0:
                        continue

                    for trans in clusters[clus]:
                        newArray = []
                        for exon in ranges[trans]:
                            for intersec in res:
                                # complete overlap --> skip completely
                                if exon[0] >= intersec[0] and exon[1] <= intersec[1]:
                                    break
                                # no overlap if start[intersec] > end[exon] --> add original range to new array
                                elif intersec[0] > exon[1]:
                                    if exon not in newArray:
                                        newArray.append(exon)
                                    continue      # check next intersection point
                                # complete overlap if start[intersec] >= start[exon] AND end[intersec] <= end[exon] --> don't add to new array
                                elif intersec[0] >= exon[0] and intersec[1] <= exon[1]:
                                    if intersec[0] == exon[0]:
                                        newArray.append([intersec[1] + 1, exon[1]])
                                    elif intersec[1] == exon[1]:
                                        newArray.append([exon[0], intersec[0] - 1])
                                    else: # in middle --> append portion before and after intersection
                                        newArray.append([exon[0], intersec[0] - 1])
                                        newArray.append([intersec[1] + 1, exon[1]])
                                # partial overlap if start[intersec] >= start[exon] (intersection extends right) --> truncate and add to new array
                                elif intersec[0] <= exon[1] and intersec[0] >= exon[0]:
                                    newArray.append([exon[0], intersec[0] - 1])
                                # partial overlap if end[intersec] <= end[exon] (intersection extends left) --> truncate and add to new array
                                elif intersec[1] >= exon[0] and intersec[1] <= exon[1]:
                                    newArray.append([intersec[1] + 1, exon[1]])
                                # Reach last one and exons remaining are greater than intersections --> no overlap --> append exon
                                elif intersec == res[-1] and exon[0] > intersec[1]:
                                    newArray.append(exon)
                        ranges[trans] = newArray

                    for trans in clusters[compclus]:
                        newArray = []
                        for exon in ranges[trans]:
                            for intersec in res:
                                # complete overlap --> skip completely
                                if exon[0] >= intersec[0] and exon[1] <= intersec[1]:
                                    # if exon == intersec:
                                    break
                                # no overlap if start[intersec] > end[exon] --> add original range to new array
                                elif intersec[0] > exon[1]:
                                    if exon not in newArray:
                                        newArray.append(exon)
                                    continue      # check next intersection point
                                # complete overlap if start[intersec] >= start[exon] AND end[intersec] <= end[exon] --> don't add to new array
                                elif intersec[0] >= exon[0] and intersec[1] <= exon[1]:
                                    if intersec[0] == exon[0]:
                                        newArray.append([intersec[1] + 1, exon[1]])
                                    elif intersec[1] == exon[1]:
                                        newArray.append([exon[0], intersec[0] - 1])
                                    else: # in middle --> append portion before and after intersection
                                        newArray.append([exon[0], intersec[0] - 1])
                                        newArray.append([intersec[1] + 1, exon[1]])
                                # partial overlap if start[intersec] >= start[exon] (intersection extends right) --> truncate and add to new array
                                elif intersec[0] <= exon[1] and intersec[0] >= exon[0]:
                                    newArray.append([exon[0], intersec[0] - 1])
                                # partial overlap if end[intersec] <= end[exon] (intersection extends left) --> truncate and add to new array
                                elif intersec[1] >= exon[0] and intersec[1] <= exon[1]:
                                    newArray.append([intersec[1] + 1, exon[1]])
                                # Reach last one and exons remaining are greater than intersections --> no overlap --> append exon
                                elif intersec == res[-1] and exon[0] > intersec[1]:
                                    newArray.append(exon)
                        ranges[trans] = newArray

    # Repeat to get rid of any remaining overlaps between original copy and current copy
    # for every combination of clusters (other than comparing with self)
    for clus, compclus in itertools.combinations(clusters, 2):
        if clus[0] == compclus[0]:
            for transcript in itertools.product(clusters[clus], clusters[compclus]):
                # Overlap detected
                if min_start[transcript[1]] in range(min_start[transcript[0]], max_end[transcript[0]]):
                    overlaps = intersections(orig_ranges[transcript[0]], ranges[transcript[1]])
                    overlaps += intersections(orig_ranges[transcript[1]], ranges[transcript[0]])
                    res = []
                    [res.append(x) for x in overlaps if x not in res]
                    if len(res) == 0:
                        continue

                    for trans in clusters[clus]:
                        newArray = []
                        for exon in ranges[trans]:
                            for intersec in res:
                                # complete overlap --> skip completely
                                if exon[0] >= intersec[0] and exon[1] <= intersec[1]:
                                    break
                                # no overlap if start[intersec] > end[exon] --> add original range to new array
                                elif intersec[0] > exon[1]:
                                    if exon not in newArray:
                                        newArray.append(exon)
                                    continue      # check next intersection point
                                # complete overlap if start[intersec] >= start[exon] AND end[intersec] <= end[exon] --> don't add to new array
                                elif intersec[0] >= exon[0] and intersec[1] <= exon[1]:
                                    if intersec[0] == exon[0]:
                                        newArray.append([intersec[1] + 1, exon[1]])
                                    elif intersec[1] == exon[1]:
                                        newArray.append([exon[0], intersec[0] - 1])
                                    else: # in middle --> append portion before and after intersection
                                        newArray.append([exon[0], intersec[0] - 1])
                                        newArray.append([intersec[1] + 1, exon[1]])
                                # partial overlap if start[intersec] >= start[exon] (intersection extends right) --> truncate and add to new array
                                elif intersec[0] <= exon[1] and intersec[0] >= exon[0]:
                                    newArray.append([exon[0], intersec[0] - 1])
                                # partial overlap if end[intersec] <= end[exon] (intersection extends left) --> truncate and add to new array
                                elif intersec[1] >= exon[0] and intersec[1] <= exon[1]:
                                    newArray.append([intersec[1] + 1, exon[1]])
                                # Reach last one and exons remaining are greater than intersections --> no overlap --> append exon
                                elif intersec == res[-1] and exon[0] > intersec[1]:
                                    newArray.append(exon)
                        ranges[trans] = newArray

                    for trans in clusters[compclus]:
                        newArray = []
                        for exon in ranges[trans]:
                            for intersec in res:
                                # complete overlap --> skip completely
                                if exon[0] >= intersec[0] and exon[1] <= intersec[1]:
                                    # if exon == intersec:
                                    break
                                # no overlap if start[intersec] > end[exon] --> add original range to new array
                                elif intersec[0] > exon[1]:
                                    if exon not in newArray:
                                        newArray.append(exon)
                                    continue      # check next intersection point
                                # complete overlap if start[intersec] >= start[exon] AND end[intersec] <= end[exon] --> don't add to new array
                                elif intersec[0] >= exon[0] and intersec[1] <= exon[1]:
                                    if intersec[0] == exon[0]:
                                        newArray.append([intersec[1] + 1, exon[1]])
                                    elif intersec[1] == exon[1]:
                                        newArray.append([exon[0], intersec[0] - 1])
                                    else: # in middle --> append portion before and after intersection
                                        newArray.append([exon[0], intersec[0] - 1])
                                        newArray.append([intersec[1] + 1, exon[1]])
                                # partial overlap if start[intersec] >= start[exon] (intersection extends right) --> truncate and add to new array
                                elif intersec[0] <= exon[1] and intersec[0] >= exon[0]:
                                    newArray.append([exon[0], intersec[0] - 1])
                                # partial overlap if end[intersec] <= end[exon] (intersection extends left) --> truncate and add to new array
                                elif intersec[1] >= exon[0] and intersec[1] <= exon[1]:
                                    newArray.append([intersec[1] + 1, exon[1]])
                                # Reach last one and exons remaining are greater than intersections --> no overlap --> append exon
                                elif intersec == res[-1] and exon[0] > intersec[1]:
                                    newArray.append(exon)
                        ranges[trans] = newArray

    for gene in ranges:
        for ran in ranges[gene]:
            if gene in truncated_lengths:
                truncated_lengths[gene] += ran[1] - ran[0] + 1
            else:
                truncated_lengths[gene] = ran[1] - ran[0] + 1
    return truncated_lengths

"""
UTS lengths - take truncated lengths calculate within clusters
"""

def UTS_lengths(mean_fragment_length = 200):
    truncated_length = truncate_genes(clusters)
    pastTranscript = None
    for clus in clusters:
        if len(clusters[clus]) == 1:
            UTSLengths[clusters[clus][0]] = truncated_length[clusters[clus][0]]
        else:
            # Check if any transcript in cluster has multiple exons
            multi = [True for transcript in clusters[clus] if transcript in multiple_exons]
            # No multi-exon genes
            if len(multi) == 0:
                for transcript in clusters[clus]:
                    if transcript == clusters[clus][0]:
                        pastTranscript = transcript
                        continue
                    elif transcript == clusters[clus][-1]:
                        UTSLengths[pastTranscript] = truncated_length[pastTranscript] - truncated_length[transcript]
                        UTSLengths[transcript] = truncated_length[transcript]
                    else:
                        UTSLengths[pastTranscript] = truncated_length[pastTranscript] - truncated_length[transcript]
                        pastTranscript = transcript
            # Multi-exon genes
            else:
                if len(clusters[clus]) == 2:
                    overlaps = intersections(ranges[clusters[clus][0]], ranges[clusters[clus][1]])
                    overlapLength = 0
                    UTSLengths[clusters[clus][0]] = truncated_length[clusters[clus][0]]
                    UTSLengths[clusters[clus][1]] = truncated_length[clusters[clus][1]]
                    # One gene completely encapsulated as overlap region
                    if overlaps == ranges[clusters[clus][0]]:
                        # Other gene truncated to remove overlap
                        for ar in overlaps:
                            overlapLength += ar[1] - ar[0] + 1
                        UTSLengths[clusters[clus][1]] -= overlapLength
                    elif overlaps == ranges[clusters[clus][1]]:
                        for ar in overlaps:
                            overlapLength += ar[1] - ar[0] + 1
                        UTSLengths[clusters[clus][0]] -= overlapLength
                    else:
                        for ar in overlaps:
                            # Overlap is a shared exon between two clusters --> not unique --> remove entirely
                            if ar in ranges[clusters[clus][0]] and ar in ranges[clusters[clus][1]]:
                                UTSLengths[clusters[clus][0]] -= (ar[1] - ar[0] + 1)
                                UTSLengths[clusters[clus][1]] -= (ar[1] - ar[0] + 1)
                                overlaps.remove(ar)
                        # Otherwise, overlap is between exons but only segment of them
                        # If the transcript has any unique ranges of its own other than these overlaps, get rid of overlaps from there
                        for ar in overlaps:
                            overlapLength += ar[1] - ar[0] + 1
                        UTSLengths[clusters[clus][0]] = truncated_length[clusters[clus][0]] - overlapLength
                        UTSLengths[clusters[clus][1]] = truncated_length[clusters[clus][1]] - overlapLength
    effectiveLengthCalc(mean_fragment_length)

"""
Effective length calculation
"""
def effectiveLengthCalc(mean_fragment_length = 200):
    global effectiveLengthA
    mean_fragment_length = int(mean_fragment_length)
    if mean_fragment_length == 0:
        effectiveLength = UTSLengths.copy()
        effectiveLengthA = effectiveLength.copy()
        return

    half_mean_fragment_length = mean_fragment_length / 2
    effectiveLength = UTSLengths.copy()

    UTS_min_start = {}
    UTS_max_end = {}
    for transcript in ranges:
        UTS_min_start[transcript] = ranges[transcript][0][0]
        UTS_max_end[transcript] = ranges[transcript][-1][1]

    # Loop through samples and corresponding clusters
    for clus in clusters:
        # Correction to be completed if cluster only contains one transcript
        if len(clusters[clus]) == 1:
            transcript = clusters[clus][0]
            if UTSLengths[transcript] < 300:
                effectiveLength[transcript] = UTSLengths[transcript]
            # if truncated flag = true (i.e. this overlapped and was truncated):
            elif ranges[transcript] != orig_ranges[transcript]:
                # check if 5' end was truncated
                if min_start[transcript] != UTS_min_start[transcript]:
                    diff = UTS_min_start[transcript] - min_start[transcript]
                    if diff < half_mean_fragment_length:
                        trunc = half_mean_fragment_length - diff
                        effectiveLength[transcript] -= (trunc + 1)
                else:
                    effectiveLength[transcript] -= half_mean_fragment_length

                # check if 3' end was truncated
                if max_end[transcript] != UTS_max_end[transcript]:
                    diff = max_end[transcript] - UTS_max_end[transcript]
                    if diff < half_mean_fragment_length:
                        trunc = half_mean_fragment_length - diff
                        effectiveLength[transcript] -= (trunc + 1)
                else:
                    effectiveLength[transcript] -= half_mean_fragment_length
            else:
                effectiveLength[transcript] -= (mean_fragment_length + 1)
        # Else, for each transcript in cluster, correct for overlaps
        else:
            for transcript in clusters[clus]:
                # first gene in cluster
                if transcript == clusters[clus][0]:
                    if UTSLengths[transcript] < 300:
                        effectiveLength[transcript] = UTSLengths[transcript]
                    # check if truncated at all
                    elif (ranges[transcript] != orig_ranges[transcript]):
                        # check if 5' end was truncated
                        if (min_start[transcript] != UTS_min_start[transcript]) and clus[0] == "+":
                            diff = UTS_min_start[transcript] - min_start[transcript]
                            if diff < half_mean_fragment_length:
                                trunc = half_mean_fragment_length - diff
                                effectiveLength[transcript] -= (trunc + 1)
                        elif (max_end[transcript] != UTS_max_end[transcript]) and clus[0] == "-":
                            diff = max_end[transcript] - UTS_max_end[transcript]
                            if diff < half_mean_fragment_length:
                                trunc = half_mean_fragment_length - diff
                                effectiveLength[transcript] -= (trunc + 1)
                        else:
                            effectiveLength[transcript] -= (half_mean_fragment_length + 1)
                    else:
                        effectiveLength[transcript] -= (half_mean_fragment_length + 1)
                elif transcript == clusters[clus][-1]:
                    if UTSLengths[transcript] < 300:
                        effectiveLength[transcript] = UTSLengths[transcript]
                    # check if truncated at all
                    elif (ranges[transcript] != orig_ranges[transcript]):
                        # check if 3' end was truncated
                        if (min_start[transcript] != UTS_min_start[transcript]) and clus[0] == "-":
                            diff = UTS_min_start[transcript] - min_start[transcript]
                            if diff < half_mean_fragment_length:
                                trunc = half_mean_fragment_length - diff
                                effectiveLength[transcript] -= (trunc + 1)
                        elif (max_end[transcript] != UTS_max_end[transcript]) and clus[0] == "+":
                            diff = max_end[transcript] - UTS_max_end[transcript]
                            if diff < half_mean_fragment_length:
                                trunc = half_mean_fragment_length - diff
                                effectiveLength[transcript] -= (trunc + 1)
                        else:
                            effectiveLength[transcript] -= (half_mean_fragment_length + 1)
                    else:
                        effectiveLength[transcript] -= (half_mean_fragment_length + 1)
                else:
                    effectiveLength[transcript] = UTSLengths[transcript]
    effectiveLengthA = effectiveLength.copy()

"""
Determine merged names from GTF file --> allows for comparison with mmquant file
"""
def mergedGenes():
    pastTranscript = None
    for clus in clusters:
        # If only one gene in cluster, then assign appropriate length, save same name, and save associated reads
        if len(clusters[clus]) == 1:
            splitMergedNames[clusters[clus][0]] = [clusters[clus][0]]
            mergedNameClusters[clus] = [clusters[clus][0]]
        # Else, break the merged name into its components and calculate associated length and reads
        else:
            # Forward stranded case
            if clus[0] == "+":
                mergedName = ""
                mergedNameClusters[clus] = []
                for transcript in clusters[clus]:
                    # Case for first transcript
                    if transcript == clusters[clus][0]:
                        mergedName += transcript
                        pastTranscript = transcript
                    # Append to create merged name
                    else:
                        overlaps = intersections(ranges[transcript], ranges[pastTranscript])
                        if ranges[transcript] == overlaps:
                            mergedName = transcript + "--" + mergedName
                        else:
                            splitMergedNames[mergedName] = [mergedName]
                            mergedNameClusters[clus].append(mergedName)
                            mergedName = transcript
                            splitMergedNames[mergedName] = [mergedName]
                            mergedNameClusters[clus].append(mergedName)
                            continue
                    # Store information about length associated with this name
                    sort_split = sorted(mergedName.split("--"))
                    splitMergedNames.update({mergedName: sort_split})
                    mergedNameClusters[clus].append(mergedName)
                mergedName = ""
                pastTranscript = None
            # Reverse stranded case
            elif clus[0] == "-":
                mergedName = ""
                mergedNameClusters[clus] = []
                for transcript in clusters[clus]:
                    # Case for first transcript
                    if transcript == clusters[clus][0]:
                        mergedName += transcript
                    # Append to create merged name
                    else:
                        mergedName = transcript + "--" + mergedName
                    # Store information about length associated with this name
                    sort_split = sorted(mergedName.split("--"))
                    splitMergedNames.update({mergedName: sort_split})
                    mergedNameClusters[clus].append(mergedName)
                mergedName = ""
                pastTranscript = None

"""
Read in necessary read count values from mmquant file
Use merged names to determine which are relevant vs not
"""
def readOrigReads(mmquantfile):
    mergedGenes()
    origReads = {}
    # Read mmquant file line by line
    with open(mmquantfile, 'r') as f:
        for line in f:
            line = re.sub("\s+", " ", line)
            curRow = line.strip().split(' ')
            rowGenes = curRow[0].split('--')
            # Check if transcript at line exists within sMergedReads and thus needs to be evaluated
            if sorted(rowGenes) in splitMergedNames.values():
                # If exists, then find location in sMergedReads and store reads associated with this value
                name = list(splitMergedNames.keys())[list(splitMergedNames.values()).index(sorted(rowGenes))]
                if name in origReads.keys():
                    origReads[name].append(curRow[1])
                else:
                    origReads[name] = [curRow[1]]

    # Go through and assign oReads according to the name of the transcript
    asgnOReads = {}
    newNameClus = {}

    # After storing associated reads, determine UTS reads and lengths for each transcript
    for clus in mergedNameClusters:
        try:
            # if cluster has only one transcript, then directly assign same values
            if len(mergedNameClusters[clus]) == 1:
                asgnOReads[mergedNameClusters[clus][0]] = origReads[mergedNameClusters[clus][0]]
                newNameClus[clus] = [mergedNameClusters[clus][0]]
            else:
                # Otherwise, loop through each transcript in cluster
                for transcript in mergedNameClusters[clus]:
                    # Directly save information for first transcript
                    if transcript == mergedNameClusters[clus][0]:
                        if transcript in origReads.keys():
                            asgnOReads[transcript] = origReads[transcript]
                            newNameClus[clus] = [transcript]
                    else:
                        # Truncate mergedName to get next transcript name and assign
                        try:
                            corNameIndex = transcript.index("--")
                            corName = transcript[: corNameIndex]
                        except ValueError:
                            corName = transcript
                        if transcript in origReads.keys():
                            asgnOReads[corName] = origReads[transcript]
                        else:
                            asgnOReads[corName] = '0'
                        if clus in newNameClus.keys():
                            newNameClus[clus].append(corName)
                        else:
                            newNameClus[clus] = [corName]
        except KeyError:
            asgnOReads[mergedNameClusters[clus][0]] = '0'

    global assignedOReads
    assignedOReads = asgnOReads.copy()
    global newNameClusters
    newNameClusters = newNameClus.copy()

"""
Calculate RPM and RPKM values (before normalizing)
"""
def readspermil(GTFfile, mmquantfile, samtools):
    readOrigReads(mmquantfile)
    # variable setup
    transcriptRPM = {}
    transcriptRPKM = {}
    effectiveLength = {}

    scaleF = int(samtools) / (1000000)

    for transcript in assignedOReads:
        transcriptRPM[transcript] = int(assignedOReads[transcript][0]) / scaleF
        transcriptRPKM[transcript] = transcriptRPM[transcript] / (effectiveLengthA[transcript] / 1000)

    global RPKM
    RPKM = transcriptRPKM.copy()

"""
Correct RPKMs within cluster and report values for final output
"""
def cNormReads(GTFfile, mmquantfile, samtools, mean_frag_length):
    readspermil(GTFfile, mmquantfile, samtools)
    # set up variables and attain calculated reads
    cNorms = {}
    pastTranscript = None

    # Loop through samples and corresponding clusters
    for clus in newNameClusters:
        # No correction to be completed if cluster only contains one transcript
        if len(newNameClusters[clus]) == 1:
            cNorms[newNameClusters[clus][0]] = round(RPKM[newNameClusters[clus][0]], 2)
        # Else, for each transcript in cluster, correct for overlaps
        else:
            for transcript in newNameClusters[clus]:
                if transcript in RPKM:
                    if transcript == newNameClusters[clus][0]:
                        cNorms[transcript] = round(RPKM[transcript], 2)
                        pastTranscript = transcript
                        continue
                    else:
                        nReads = RPKM[transcript] - RPKM[pastTranscript]
                        # If the final number of reads in a corrected value are negative, set to 0
                        if nReads < 0:
                            nReads = 0
                        cNorms[transcript] = round(nReads, 2)
                        pastTranscript = transcript
    # Calculate expected counts (reverse of RPKM calculation but after normalization step is completed)
    expRPM = {}
    expCounts = {}

    for transcript in RPKM:
        expRPM[transcript] = RPKM[transcript] * (GTF_length[transcript] / 1000)
        expCounts[transcript] = int(expRPM[transcript] * (int(samtools) / (1000000)))
    global expected_counts
    expected_counts = expCounts.copy()

    # GTF length - mean fragment length
    GTF_eff = {}

    for transcript in GTF_length:
        if int(GTF_length[transcript]) > int(mean_frag_length):
            GTF_eff[transcript] = int(GTF_length[transcript]) - int(mean_frag_length) - 1

        else:
            GTF_eff[transcript] = int(GTF_length[transcript])
    global GTF_effective_length
    GTF_effective_length = GTF_eff.copy()
    return cNorms

"""
Calls for calculations and then writes output as tab-delimited text file saved to name requested by user
    args: <input-GTF-file> <input-mmquant-file> <log_final_bam> <output-file-name>
    returns: tab-delimited text file containing output organized with transcript name followed by corresponding
             calculation for RPKM value
"""

def writeFile(annotationGTF, mmquantFile, output, samtools, mean_length = 200):
    dataframe = pd.read_csv(annotationGTF, sep="\t", lineterminator='\n', names=['seqname','source','feature','start','end','score','strand','frame','attribute'])
    dataframe['seqname'] = dataframe['seqname'].apply(lambda s: s.replace('\n',''))
    grouped = dataframe.groupby('seqname')
    keys = grouped.groups.keys()
    run = 1
    for key in keys:
        splitdf = grouped.get_group(key)
        global clusters
        clusters = {}
        global origLengths
        origLengths = {}
        global UTSLengths
        UTSLengths = {}
        global all_start
        all_start = {}
        global all_end
        all_end = {}
        global ranges
        ranges = {}
        global min_start
        min_start = {}
        global max_end
        max_end = {}
        global GTF_length
        GTF_length = {}
        global multiple_exons
        multiple_exons = []
        global exon_intervals
        exon_intervals = {}
        global orig_ranges
        orig_ranges = {}
        global UTS_min_start
        UTS_min_start = {}
        global UTS_min_end
        UTS_min_end = {}
        global splitMergedNames
        splitMergedNames = {}
        global mergedNameClusters
        mergedNameClusters = {}
        global effectiveLengthA
        effectiveLengthA = {}
        global assignedOReads
        assignedOReads = {}
        global newNameClusters
        newNameClusters = {}
        global RPKM
        RPKM = {}
        global expected_counts
        expected_counts = {}
        global GTF_effective_length
        GTF_effective_length = {}
        setup(splitdf, mean_length)
        cNorms = cNormReads(annotationGTF, mmquantFile, samtools, mean_length)
        out_file = open(output, "a")
        if run == 1:
            out_file.write("gene_id\tlength\teffective_length\texpected_counts\tRPKM\n")
            run = 2

        for transcript in cNorms:
            out_file.write("%s\t%s\t%s\t%s\t%s\n" % (transcript, GTF_length[transcript], GTF_effective_length[transcript], expected_counts[transcript], cNorms[transcript]))
    out_file.close()

"""
Beginning of program: manages correct argument entry and begins execution of calculations
    args: UTS_GTF.py <input-GTF-file> <input-mmquant-file> <output-file-name> <optional: mean fragment length>
    returns: writes a tab-delimited text file containing output organized with transcript name followed by corresponding
             calculation for RPKM value
"""

def main():
    n = len(sys.argv)
    if n == 5:
        writeFile(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
        print("Successfully ran script. Output saved in given output file name")
    elif n==6:
        writeFile(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        print("Successfully ran script. Output saved in given output file name")
    else:
        print("Error: syntax to use script is python3 <python-script-name> <input-GTF-file> <input-mmquant-file> <output-file-name> <samtools value> <optional: mean fragment length>")

if __name__ == "__main__":
    main()