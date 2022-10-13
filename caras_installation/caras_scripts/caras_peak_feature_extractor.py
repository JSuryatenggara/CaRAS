#!/usr/bin/env python3
#pyright: reportUnboundVariable=false

print('Importing required modules')
import pandas as pd
import subprocess
import multiprocessing
import os
import argparse
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from sklearn.manifold import MDS
import rbo
from prince import MCA
from prince import MFA

# import numpy as np
# import Levenshtein as lev
# from scipy.spatial import distance
# from sklearn.cluster import AgglomerativeClustering
# from sklearn.feature_selection import VarianceThreshold
# from kmodes.kmodes import KModes
# from functools import partial
# import sys
# import csv
# from sklearn.metrics import jaccard_score
# from scipy import stats



print('Defining required functions')

def file_basename_parsing(file_relative_path_list):
    absolute_path = [os.path.abspath(relative_path_element) for relative_path_element in file_relative_path_list]
        
    original_name = []
    original_extension = []

    for absolute_path_element in absolute_path:
        current_extension = '.bam'
        current_file_name = absolute_path_element.split('/')[-1].strip(current_extension)
        original_name.append(current_file_name)
        original_extension.append(current_extension)
        
    return absolute_path, original_name, original_extension



# Function to generate Peak_ID out of chromosomal coordinates of the current peak
# INPUT     - chr_value     - the chromosome number at which the current peak is located
#           - start_value   - the smallest base number at which the current peak is located
#           - end_value     - the largest base number at which the current peak is located
# OUTPUT    - the current Peak_ID, made out of the current peak coordinates, in the format of chr:start-end
def peak_ID_generator_function(chr_value, start_value, end_value):
    return '{}:{}-{}'.format(chr_value, start_value, end_value)



# Function to split a pipe-separated string, converting it into a list, then count the number of elements
# INPUT     - peak_caller_combination - pipe-separated string containing name of the peak callers that called for the current peak
# OUTPUT    - the number of the peak callers names contained in the peak_caller_combination string
def overlap_number_calculator_function(peak_caller_combination):
    return len(peak_caller_combination.split('|'))



# Function to only output the input values before the brackets
# INPUT     - annotation - the annotation value of the current peak in the default HOMER format
# OUTPUT    - simplified annotation: without the brackets and bracketed values that comes after
def annotation_simplifier_function(annotation):
    if isinstance(annotation, float):
        return 'None'
    else:
        return annotation.split('(')[0]

def deltaIfRemoved(row, counts, countsAll, numClusters, multiplier, power):
    ret = 0.0

    for i in range(len(row)):
        if row[i] == 1:
            ret -= pow((counts[i] * multiplier + 1.0) / (countsAll[i] * multiplier + numClusters + 0.0), power)
            ret += pow(((counts[i] - 1.0) * multiplier + 1.0) / (countsAll[i] * multiplier + numClusters + 0.0), power)

    return ret

def deltaIfAdded(row, counts, countsAll, numClusters, multiplier, power):
    ret = 0.0

    for i in range(len(row)):
        if row[i] == 1:
            ret += pow(((counts[i] + 1.0) * multiplier + 1.0) / (countsAll[i] * multiplier + numClusters + 0.0), power)
            ret -= pow((counts[i] * multiplier + 1.0) / (countsAll[i] * multiplier + numClusters + 0.0), power)

    return ret

def popc(samples, multiplier = 1000.0, power = 10.0):
    kmeans = KMeans(n_clusters=int(len(samples)/2), random_state=0).fit(samples)

    labels = kmeans.labels_

    clusters = []
    clusters_counts = []

    counts_all = [0 for i in range(len(samples[0]))]

    for i in range(max(labels) + 1):
        clusters.append([])
        clusters_counts.append([0 for j in range(len(samples[0]))])

    for i in range(len(samples)):
        cs = samples[i]
        cl = labels[i]
        clusters[cl].append(i)

        for j in range(len(cs)):
            if cs[j] == 1:
                clusters_counts[cl][j] += 1
                counts_all[j] += 1

    changed = True

    while changed:
        changed = False
        i = 0 
        while i < len(clusters):
            j = 0
            while j < len(clusters[i]):
                largestGainWhere = -1
                largestGain = -1.0
                deltaBase = deltaIfRemoved(samples[clusters[i][j]], clusters_counts[i], counts_all, len(clusters), multiplier, power)
                for k in range(len(clusters)):
                    if i != k:
                        delta = deltaBase + deltaIfAdded(samples[clusters[i][j]], clusters_counts[k], counts_all, len(clusters), multiplier, power)
                        if delta > largestGain:
                            largestGain = delta
                            largestGainWhere = k
                if largestGain > 0:
                    changed = True
                    cs = clusters[i][j]
                    del clusters[i][j]
                    clusters[largestGainWhere].append(cs)
                    for k in range(len(samples[0])):
                        if samples[cs][k] == 1:
                            clusters_counts[i][k] -= 1
                            clusters_counts[largestGainWhere][k] += 1
                    j -= 1
                j += 1

            if len(clusters[i]) == 0:
                del clusters[i]
                del clusters_counts[i]
                i -= 1
            i += 1

    for i in range(len(clusters)):
        for j in range(len(clusters[i])):
            labels[clusters[i][j]] = i
    
    return labels



def popc_plot(labelsWithSamples, title):
        fig = plt.gcf()
        fig.set_size_inches(12, 7)

        colors = ['b','g','r','c','m','y','k']
        labelsWithSamples = sorted(labelsWithSamples, key=lambda s: s[0])

        for i in range(len(labelsWithSamples)):
                for j in range(len(labelsWithSamples[i][1])):
                        if labelsWithSamples[i][1][j] == 1:
                                plt.plot(j, i, colors[labelsWithSamples[i][0] % 7] + 'o', markersize=1)

        plt.title(title)
        plt.savefig('{}/{}_POPC_{}.png'.format(output_folder, dataset_name, title), dpi = 300)
        plt.clf()
        print('POPC plot saved as: {}/{}_POPC.png'.format(output_folder, dataset_name))



def called_peak_feature_extraction_function(interval_array_row):
    feature_replicate_list = []
    feature_value_list_of_lists = []
    for input_array_row in input_array:
        if interval_array_row[0] == input_array_row[0]:
            if int(interval_array_row[1]) - interpeaks_distance < int(input_array_row[2]):
                if int(interval_array_row[2]) + interpeaks_distance > int(input_array_row[1]):
                    peak_caller_list = input_array_row[3].split('|')
                    # print(peak_caller_list)
                    feature_value_list = []
                    for feature in feature_list:      
                        if 'overlaps' in feature:
                            feature_replicate = feature.split('_')[-1]
                            feature_overlap_number = feature.split('_')[0].strip('overlaps')
                            if len(peak_caller_list) >= int(feature_overlap_number):
                                if all(peak_caller.endswith(str(feature_replicate)) for peak_caller in peak_caller_list):
                                    if feature_replicate not in feature_replicate_list:
                                        feature_replicate_list.append(feature_replicate)
                                    feature_value_list.append(1)
                                else:
                                    feature_value_list.append(0)
                            else:
                                feature_value_list.append(0)

                        elif '-' in feature:
                            feature_replicate = feature.split('_')[-1]
                            feature_element_list = feature.split('_')[0].split('-')
                            feature_element_list = ['{}_{}'.format(feature_element, feature_replicate) for feature_element in feature_element_list]
                            if all(feature_element in peak_caller_list for feature_element in feature_element_list):
                                if feature_replicate not in feature_replicate_list:
                                    feature_replicate_list.append(feature_replicate)                                
                                feature_value_list.append(1)
                            else:
                                feature_value_list.append(0)

                        elif '+' in feature:
                            feature_replicate = feature.split('_')[-1]
                            feature_element_list = feature.split('_')[0].split('+')
                            feature_element_list = ['{}_{}'.format(feature_element, feature_replicate) for feature_element in feature_element_list]
                            if any(feature_element in peak_caller_list for feature_element in feature_element_list):
                                if feature_replicate not in feature_replicate_list:
                                    feature_replicate_list.append(feature_replicate)                                
                                feature_value_list.append(1)
                            else:
                                feature_value_list.append(0)

                        else:
                            feature_replicate = feature.split('_')[-1]
                            if feature in peak_caller_list:
                                if feature_replicate not in feature_replicate_list:
                                    feature_replicate_list.append(feature_replicate)                                
                                feature_value_list.append(1)
                            else:
                                feature_value_list.append(0)

                    # print(feature_value_list)

                    feature_value_list_of_lists.append(feature_value_list)

    summed_feature_value_list = [sum(x) for x in zip(*feature_value_list_of_lists)]

    binarized_feature_value_list = [0 if x == 0 else 1 for x in summed_feature_value_list]

    if len(feature_replicate_list) >= 1:

        interval_array_row = interval_array_row + binarized_feature_value_list
        return interval_array_row



def called_peak_ranking_function(peaklist_df):

    chr_list = peaklist_df['Chr'].tolist()
    start_list = peaklist_df['Start'].tolist()
    end_list = peaklist_df['End'].tolist()

    peak_rank_list = []

    for row_counter in range(peaklist_df.shape[0]):
        for interval_array_row in interval_array:
            if interval_array_row[0] == chr_list[row_counter]:
                if int(interval_array_row[1]) - interpeaks_distance < int(end_list[row_counter]):
                    if int(interval_array_row[2]) + interpeaks_distance > int(start_list[row_counter]):
                        current_peak_ID = '{}:{}-{}'.format(interval_array_row[0], interval_array_row[1], interval_array_row[2])
                        if current_peak_ID not in peak_rank_list:
                            peak_rank_list.append(current_peak_ID)

                        # peak_caller_list = combination_list[row_counter].split('|')

                        # if 'overlaps' in current_filter_string:
                        #     feature_overlap_number = current_filter_string.strip('overlaps')
                        #     if len(peak_caller_list) >= int(feature_overlap_number):
                        #         peak_rank_list.append('{}:{}-{}'.format(interval_array_row[0], interval_array_row[1], interval_array_row[2]))

                        # elif '-' in current_filter_string:
                        #     feature_element_list = current_filter_string.split('-')
                        #     if all(feature_element in peak_caller_list for feature_element in feature_element_list):                            
                        #         peak_rank_list.append('{}:{}-{}'.format(interval_array_row[0], interval_array_row[1], interval_array_row[2]))

                        # elif '+' in current_filter_string:
                        #     feature_element_list = current_filter_string.split('+')
                        #     if any(feature_element in peak_caller_list for feature_element in feature_element_list):                           
                        #         peak_rank_list.append('{}:{}-{}'.format(interval_array_row[0], interval_array_row[1], interval_array_row[2]))

                        # else:
                        #     if current_filter_string in peak_caller_list:                           
                        #         peak_rank_list.append('{}:{}-{}'.format(interval_array_row[0], interval_array_row[1], interval_array_row[2]))

    return peak_rank_list



def rbo_calculator_function(peak_ranking_list_a):
    
    rbo_distance_list = []

    for counter, peak_ranking_list_b in enumerate(peak_ranking_array):
        rbo_similarity = rbo.RankingSimilarity(peak_ranking_list_a, peak_ranking_list_b).rbo()
        rbo_distance = (rbo_similarity - 1) * -1
        rbo_distance_list.append(rbo_distance)

    return rbo_distance_list



print('Setting argument parser')

parser = argparse.ArgumentParser()

parser.add_argument('--thread', 
                    help = '<Optional> Maximum number of processes to use. Default is half the maximum available.', 
                    type = int, 
                    choices = range(1, (multiprocessing.cpu_count() + 1), 1),
                    metavar = "[1-{}]".format(multiprocessing.cpu_count()),
                    default = int(multiprocessing.cpu_count() / 2))

parser.add_argument('--interval_bed', 
                    help = '<Required> Input', 
                    required = True)

parser.add_argument('--input_tsv', 
                    help = '<Required> Input', 
                    required = True)

parser.add_argument('--peaklist_tsv', 
                    help = '<Required> All peaks calculated', 
                    nargs = '+')

parser.add_argument('--output', 
                    help = '<Required> Your desired output folder.', 
                    required = True)

parser.add_argument('--setname', 
                    help = '<Required> The prefix that will be used to label output and intermediate files.', 
                    required = True)

parser.add_argument('--motif', 
                    help = '<Optional> Your predicted/known motif file, in HOMER matrix format, .motif extension')

parser.add_argument('--annotate_arg', 
                    help = '<Optional> Flag arguments for annotatePeaks.pl run parameterizations.',
                    nargs = '+')

parser.add_argument('--ref', 
                    help = '<Optional> Your sample organism genome reference build. Default is hg38 (human).', 
                    choices = ['hg19', 'hg38', 'mm9', 'mm10', 'dm6', 'sacCer3'], 
                    default = 'hg38')

parser.add_argument('--filter', 
                    help = '<Required> The peakset to process. Accepted values are the (space-separated) names of the peak caller (case-sensitive) to select all peaks called by the given peak caller (e.g., MACS2, Genrich, etc). Use - for AND combination (e.g. GEM-HOMER). Use + for OR combination (e.g. SICER2+SEACR)',
                    nargs = '+',
                    required = True)

parser.add_argument('--repnum', 
                    help = '<Required> Number of replicates or samples', 
                    type = int,
                    required = True)

parser.add_argument('--distance', 
                    help = '<Required> Maximum distance between peaks in different replicates to be considered as to be in the same location.', 
                    default = 0)

parser.add_argument('--min_pos', 
                    help = '<Optional> Minimum fraction of positive peaks to be accounted as an informative variable.', 
                    type = float,
                    default = 0.001)

parser.add_argument('--max_pos', 
                    help = '<Optional> Maximum fraction of positive peaks to be accounted as an informative variable.', 
                    type = float,
                    default = 0.8)


parser.add_argument('--top_rank', 
                    help = '<Required> Top k rank of terms to be considered for plot datapoints colormap.', 
                    type = int)

parser.add_argument('--database', 
                    help = '<Required> HOMER GO or pathway database to use.', 
                    choices = ['biological_process', 'molecular_function', 'cellular_component', 'interactions', 'cosmic', 'kegg', 'biocyc', 'pathwayInteractionDB', 'reactome', 'smpdb', 'wikipathways'], 
                    required = True)

args = parser.parse_args()



subprocess.run('ulimit -n 2000', shell = True)

print('Parsing arguments')

cpu_count = args.thread

interval_bed_full_path = os.path.abspath(args.interval_bed)

input_tsv_full_path = os.path.abspath(args.input_tsv)

peaklist_tsv_full_path_list = [os.path.abspath(peaklist) for peaklist in args.peaklist_tsv]

dataset_name = args.setname # Resulting files will be named based in this

output_folder = os.path.abspath(args.output) # Absolute path to the folder [dataset_name] that contains all the outputs of the suite

if args.motif:
    motif_file_full_path = os.path.abspath(args.motif) # Absolute path to the .motif file (HOMER motif matrix file)
if not args.motif:
    motif_file_full_path = ''
    
genome_ref = args.ref

if args.annotate_arg:
    homer_annotatePeaks_arg = ' '.join(args.annotate_arg)
if not args.annotate_arg:
    homer_annotatePeaks_arg = ''

minimum_positive_fraction = args.min_pos

maximum_positive_fraction = args.max_pos

interpeaks_distance = int(args.distance)

if args.top_rank:
    top_rankers = args.top_rank

database = args.database


peakset_filter_list = [peakset_filter for peakset_filter in args.filter]

new_peakset_filter_list = []

for peakset_filter in peakset_filter_list:   
    if peakset_filter.isdigit():
        peakset_filter = '{}overlaps'.format(peakset_filter)

    new_peakset_filter_list.append(peakset_filter)

peakset_filter_list = new_peakset_filter_list


n_reps = args.repnum

feature_list = []

for replicate_counter in range(n_reps):
    for peakset_filter in peakset_filter_list:
        feature_list.append('{}_{}'.format(peakset_filter, replicate_counter + 1))

peakset_filter_string = '_'.join([peakset_filter for peakset_filter in peakset_filter_list])

# cell_line_list = ['H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 
#     'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562']
# color_list = ['red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red',
#     'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue']

cell_line_list = ['CTCF', 'CTCF', 'CTCF', 'CTCF', 'CTCF', 'CTCF', 'CTCF',
    'H3K27me3', 'H3K27me3', 'H3K27me3', 'H3K27me3', 'H3K27me3', 'H3K27me3', 'H3K27me3']
color_list = ['red', 'red', 'red', 'red', 'red', 'red', 'red',
    'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue']

# cell_line_list = ['H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 
#     'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562',
#     'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1']
# color_list = ['red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red',
#     'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue',
#     'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red']

selected_cell_line_list = cell_line_list

peakwise_feature_df = None
replicatewise_feature_df = None
filtered_replicatewise_feature_df = None



if args.peaklist_tsv:

    print('Opening file: {}'.format(interval_bed_full_path)) # Reading the HOMER annotated interval file
    interval_df = pd.read_csv(interval_bed_full_path, delimiter = '\t', header = None)
    interval_df = interval_df.iloc[:,[0,1,2]]
    interval_df.sort_values(by = [0,1], inplace = True) # Sorting the intervals list based on chromosomal position
    interval_array = interval_df.values.tolist()

    peaklist_df_list = []

    for peaklist_tsv_full_path in peaklist_tsv_full_path_list:
        print('Opening file: {}'.format(peaklist_tsv_full_path))
        peaklist_df = pd.read_csv(peaklist_tsv_full_path, delimiter = '\t')
        peaklist_df.sort_values(by = ['Peak Caller Overlaps', 'Fold Change', 'ChIP Tag Count'], ascending = False, inplace = True)
        # peaklist_df[peaklist_df.columns[:(list(peaklist_df.columns).index('IDR'))+1]]
        peaklist_df = peaklist_df[['Chr', 'Start', 'End', 'Peak Caller Overlaps', 'Peak Caller Combination']]
        peaklist_df_list.append(peaklist_df)

    for filter_counter in range(len(peakset_filter_list)):        
        current_filter_string = peakset_filter_list[filter_counter]
        print('Processing peaklists based on filter: {}'.format(current_filter_string))

        filtered_peaklist_df_list = []

        if 'overlaps' in current_filter_string:
            for peaklist_df in peaklist_df_list:
                filtered_peaklist_df_list.append(peaklist_df[peaklist_df['Peak Caller Overlaps'] >= int(current_filter_string.strip('overlaps'))])
                
        elif '-' in current_filter_string:
            for peaklist_df in peaklist_df_list:
                filtered_peaklist_df_list.append(peaklist_df[peaklist_df['Peak Caller Combination'].apply(lambda x: all(item[0] in x.split('|') for item in current_filter_string.split('-')))])

        elif '+' in current_filter_string:
            for peaklist_df in peaklist_df_list:
                filtered_peaklist_df_list.append(peaklist_df[peaklist_df['Peak Caller Combination'].apply(lambda x: any(item[0] in x.split('|') for item in current_filter_string.split('+')))])

        else:
            for peaklist_df in peaklist_df_list:
                filtered_peaklist_df_list.append(peaklist_df[peaklist_df['Peak Caller Combination'].apply(lambda x: current_filter_string in x)])

        print('Calculating peak ranks')
        pool_called_peak_ranking = multiprocessing.Pool(processes = cpu_count)
        peak_ranking_array = pool_called_peak_ranking.map(called_peak_ranking_function, filtered_peaklist_df_list)
        pool_called_peak_ranking.close()
        pool_called_peak_ranking.join()

        # set_a = set()
        # dupe_a_list = [x for x in peak_ranking_list_a if x in set_a or set_a.add(x)]
        # for dupe_a in dupe_a_list:
        #     dupe_a_index_list = [i for i, x in enumerate(peak_ranking_list_a) if x == dupe_a]
        #     for dupe_a_index_list_counter, dupe_a_index in enumerate(dupe_a_index_list):
        #         if dupe_a_index_list_counter != 0:
        #             del peak_ranking_list_a[dupe_a_index]

        print('Comparing peak ranks')
        pool_rbo_calculation = multiprocessing.Pool(processes = cpu_count)
        rbo_distance_matrix = pool_rbo_calculation.map(rbo_calculator_function, peak_ranking_array)
        pool_rbo_calculation.close()
        pool_rbo_calculation.join()



        print('Writing RBO distance matrix file')
        pd.DataFrame(rbo_distance_matrix).to_csv('{}/{}_{}_d{}_RBO_distance_matrix.csv'.format(output_folder, dataset_name, current_filter_string, interpeaks_distance), sep = '\t', index = False)

        print('Separating peak profiles with Multi Dimensional Scaling')
        mds = MDS(dissimilarity = 'precomputed', random_state = 0)

        transformed_rbo_distance_matrix = mds.fit_transform(rbo_distance_matrix)
        plt.scatter(transformed_rbo_distance_matrix[:,0], transformed_rbo_distance_matrix[:,1], c = color_list, cmap = 'rainbow')
        plt.title('MDS RBO Distance')
        plt.tight_layout()
        plt.savefig('{}/{}_{}_d{}_MDS_RBO.png'.format(output_folder, dataset_name, current_filter_string, interpeaks_distance), dpi = 300)
        plt.clf()

        fig = plt.figure()
        fig.set_figheight(15)
        fig.set_figwidth(15)

        print('Condensing Distance Matrices')
        condensed_rbo_distance_matrix = squareform(rbo_distance_matrix)

        print('Generating dendrograms')

        linkage_matrix = linkage(condensed_rbo_distance_matrix, 'average')
        dendrogram(linkage_matrix, orientation = 'left', labels = cell_line_list, leaf_font_size = 10)
        plt.title("RBO Distance Dendrogram (Average Linkage)")
        plt.tight_layout()
        plt.savefig('{}/{}_{}_d{}_dendrogram_average_rbo.png'.format(output_folder, dataset_name, current_filter_string, interpeaks_distance), dpi = 300)
        plt.clf()

        linkage_matrix = linkage(condensed_rbo_distance_matrix, 'single')
        dendrogram(linkage_matrix, orientation = 'left', labels = cell_line_list, leaf_font_size = 10)
        plt.title("RBO Distance Dendrogram (Single Linkage)")
        plt.tight_layout()
        plt.savefig('{}/{}_{}_d{}_dendrogram_single_rbo.png'.format(output_folder, dataset_name, current_filter_string, interpeaks_distance), dpi = 300)
        plt.clf()

        linkage_matrix = linkage(condensed_rbo_distance_matrix, 'complete')
        dendrogram(linkage_matrix, orientation = 'left', labels = cell_line_list, leaf_font_size = 10)
        plt.title("RBO Distance Dendrogram (Complete Linkage)")
        plt.tight_layout()
        plt.savefig('{}/{}_{}_d{}_dendrogram_complete_rbo.png'.format(output_folder, dataset_name, current_filter_string, interpeaks_distance), dpi = 300)
        plt.clf()



if not os.path.isfile('{}/{}_{}_d{}_peakwise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance)):

    print('Opening file: {}'.format(interval_bed_full_path)) # Reading the HOMER annotated interval file
    interval_df = pd.read_csv(interval_bed_full_path, delimiter = '\t', header = None)
    interval_df = interval_df.iloc[:,[0,1,2]]
    interval_df.sort_values(by = [0,1], inplace = True) # Sorting the intervals list based on chromosomal position
    interval_array = interval_df.values.tolist()

    print('Opening file: {}'.format(input_tsv_full_path)) # Reading the HOMER annotated input file
    input_df = pd.read_csv(input_tsv_full_path, delimiter = '\t', header = None)
    input_df.sort_values(by = [0,1], inplace = True) # Sorting the inputs list based on chromosomal position
    input_array = input_df.values.tolist()

    print('Extracting features by each detected peak')
    pool_called_peak_feature = multiprocessing.Pool(processes = cpu_count)
    peakwise_feature_array = pool_called_peak_feature.map(called_peak_feature_extraction_function, interval_array)
    pool_called_peak_feature.close()
    pool_called_peak_feature.join()

    peakwise_feature_array = [peakwise_feature_array_row for peakwise_feature_array_row in peakwise_feature_array if peakwise_feature_array_row != None]

    peakwise_feature_df = pd.DataFrame(data = peakwise_feature_array, columns = None)
    peakwise_feature_df.to_csv('{}/{}_{}_d{}_peakwise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance), sep = '\t', index = False, header = None)
    print('Peak-wise features are saved in: {}/{}_{}_d{}_peakwise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance))



if not os.path.isfile('{}/{}_{}_d{}_replicatewise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance)):

    if peakwise_feature_df is None:

        print('Opening file: {}/{}_{}_d{}_peakwise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance))
        peakwise_feature_df = pd.read_csv('{}/{}_{}_d{}_peakwise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance), delimiter = '\t', header = None)
        peakwise_feature_array = peakwise_feature_df.values.tolist()

    print('Converting features from peak-wise to replicate-wise')

    replicatewise_feature_array = []

    for replicate_counter in range(n_reps):
        # print('Replicate counter at: ' + str(replicate_counter))
        replicatewise_feature_array_row = []

        counter_start = 3 + (replicate_counter * len(peakset_filter_list))


        for peakwise_feature_array_row in peakwise_feature_array:
            # print('Counter start: ' + str(counter_start) + ' Peak {}:{}-{}'.format(peakwise_feature_array_row[0],peakwise_feature_array_row[1],peakwise_feature_array_row[2]))

            for feature_counter in range(counter_start, counter_start + len(peakset_filter_list)):
                # print('Feature counter at: ' + str(feature_counter) + ' Peak {}:{}-{}'.format(peakwise_feature_array_row[0],peakwise_feature_array_row[1],peakwise_feature_array_row[2]) + ' Value {}'.format(peakwise_feature_array_row[feature_counter]))
                replicatewise_feature_array_row.append(peakwise_feature_array_row[feature_counter])


        if not replicatewise_feature_array_row:
            print('Replicate-wise List is empty. Something is wrong with the input files. Exiting prematurely.')
            exit()
        
        else:
            replicatewise_feature_array.append(replicatewise_feature_array_row)


    replicatewise_feature_df = pd.DataFrame(data = replicatewise_feature_array, columns = None)
    replicatewise_feature_df.to_csv('{}/{}_{}_d{}_replicatewise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance), sep = '\t', index = False, header = None)
    print('Replicate-wise features are saved in: {}/{}_{}_d{}_replicatewise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance))



if not os.path.isfile('{}/{}_{}_d{}_replicatewise_filtered.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance)):

    if replicatewise_feature_df is None:

        print('Opening file: {}/{}_{}_d{}_replicatewise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance))
        replicatewise_feature_df = pd.read_csv('{}/{}_{}_d{}_replicatewise.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance), delimiter = '\t', header = None)


    rep_label_list = ['rep{}'.format(replicate_counter + 1) for replicate_counter in range(n_reps)]

    peak_label_list = []

    for peak_counter in range(int(replicatewise_feature_df.shape[1] / len(peakset_filter_list))):
        for peakset_filter in peakset_filter_list:
            peak_label_list.append('{}_#{}'.format(peakset_filter, peak_counter + 1))

    n_peaks = replicatewise_feature_df.shape[1]

    replicatewise_feature_df.columns = peak_label_list

    replicatewise_feature_df.index = rep_label_list



    print('Number of variables before filtering: {}'.format(replicatewise_feature_df.shape[1]))

    print('Removing non-informative variables')

    included_column_index_list = []

    for replicatewise_column in range(replicatewise_feature_df.shape[1]):
        column_value_df = replicatewise_feature_df.iloc[:, replicatewise_column]
        
        try:
            pos_occurence = column_value_df.value_counts()[1]
        except:
            pos_occurence = 0
        
        # try:
        #     neg_occurence = column_value_df.value_counts()[0]
        # except:
        #     neg_occurence = 0

        if minimum_positive_fraction <= (pos_occurence / column_value_df.shape[0]) <= maximum_positive_fraction:
            # if (neg_occurence / column_value_df.shape[0]) > 0.01:
            included_column_index_list.append(replicatewise_column)

    filtered_replicatewise_feature_df = replicatewise_feature_df.iloc[:, included_column_index_list]

    print('Number of variables after filtering: {}'.format(filtered_replicatewise_feature_df.shape[1]))

    filtered_replicatewise_feature_df.to_csv('{}/{}_{}_d{}_replicatewise_filtered.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance), sep = '\t')
    print('Filtered replicate-wise features are saved in: {}/{}_{}_d{}_replicatewise_filtered.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance))



if filtered_replicatewise_feature_df is None:

    print('Opening file: {}/{}_{}_d{}_replicatewise_filtered.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance))
    filtered_replicatewise_feature_df = pd.read_csv('{}/{}_{}_d{}_replicatewise_filtered.tsv'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance), delimiter = '\t')
    filtered_replicatewise_feature_df.set_index(filtered_replicatewise_feature_df.columns[0], inplace = True, drop = True)

data = filtered_replicatewise_feature_df

data.replace(0, 'N', inplace = True)
data.replace(1, 'Y', inplace = True)

# print(data.head())
# print(data.columns.tolist()[0:10])
# print(data.index.tolist()[0:10])
# print(data.values.tolist()[0][0:10])

peakset_group_list = {}

for peakset_filter in peakset_filter_list:
    peakset_group_list[peakset_filter] = [peak_label for peak_label in data.columns if peak_label.startswith(peakset_filter)]

# print(peakset_group_list)

# print(data.shape[0])
# print(data.shape[1])

# print(data.values.tolist())
# exit()


########## MCA ##########

print('Performing dimensionality reduction with MCA')

mca = MCA(n_components = 2, n_iter = 3, random_state = 101)
mca.fit(data)
data_mca = mca.transform(data)
print(data_mca)
mca.plot_coordinates(X = data)

# print('MCA - Variance of each columns:')
# print(mca.explained_inertia)

plt.title('MCA')
plt.tight_layout()
plt.savefig('{}/{}_{}_d{}_MCA.png'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance), dpi = 300)
plt.clf()
print('Plot saved as: {}/{}_{}_d{}_MCA.png'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance))

########## MCA ##########



########## MFA ##########

print('Performing dimensionality reduction with MFA')

mfa = MFA(groups = peakset_group_list, n_components = 2, n_iter = 3, random_state = 101)
mfa.fit(data)
mfa.transform(data)
mfa.plot_row_coordinates(X = data, labels = data.index, color_labels = ['{}'.format(cell_line) for cell_line in selected_cell_line_list])

# print('MFA (Full) - Variance of each columns:')
# print(mfa.explained_inertia)

plt.title('MFA (Full)')
plt.tight_layout()
plt.savefig('{}/{}_{}_d{}_MFA_full.png'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance), dpi = 300)
plt.clf()
print('Plot saved as: {}/{}_{}_d{}_MFA_full.png'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance))

mfa.partial_row_coordinates(data)
mfa.plot_partial_row_coordinates(X = data, color_labels = ['{}'.format(cell_line) for cell_line in selected_cell_line_list])

# print('MFA (Partial) - Variance of each columns:')
# print(mfa.explained_inertia)

plt.title('MFA (Partial)')
plt.tight_layout()
plt.savefig('{}/{}_{}_d{}_MFA_partial.png'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance), dpi = 300)
plt.clf()
print('Plot saved as: {}/{}_{}_d{}_MFA_partial.png'.format(output_folder, dataset_name, peakset_filter_string, interpeaks_distance))

########## MFA ##########



# fig = plt.figure()
# fig.set_figheight(15)
# fig.set_figwidth(15)



# cluster_label_file = open('{}/{}_cluster_labels.txt'.format(output_folder, dataset_name), 'w')
# with cluster_label_file:    
#     cluster_label_file_writer = csv.writer(cluster_label_file)
#     cluster_label_file_writer.writerows(map(lambda x : [x], cao_clusters))
#     print('Cluster labels saved as: {}/{}_cluster_labels.txt'.format(output_folder, dataset_name))



# cluster_label_file = open('{}/{}_cluster_labels.txt'.format(output_folder, dataset_name), 'r')
# with cluster_label_file:
#     cluster_label_list = cluster_label_file.readlines()

# concatenated_peak_list = ['{}/{}_rep{}_all_peaks_concatenated.tsv'.format(output_folder, dataset_name, n_reps) for n_reps in range(1, len(cluster_label_list) + 1)]



# clustered_files_nested_list = []

# for num_clusters in range(cao_optimal_n_clusters):

#     clustered_files_list = []
    
#     for cluster_label_counter in range(len(cao_clusters)):
#         if cao_clusters[cluster_label_counter] ==  num_clusters:
#             clustered_files_list.append(concatenated_peak_list[cluster_label_counter])
    
#     clustered_files_nested_list.append(clustered_files_list)



# print('Generating concatenated peak list for each cluster')
# for list_counter in range(cao_optimal_n_clusters):
#     clustered_files_string = ' '.join(clustered_files_nested_list[list_counter])
#     print('Cluster {}: {}'.format(list_counter + 1, clustered_files_string))
#     subprocess.run('cat {} > {}/{}_cluster{}_all_peaks_concatenated.tsv'.format(clustered_files_string,
#         output_folder,
#         dataset_name,
#         list_counter + 1),
#         shell = True)



# if not os.path.isfile('{}/{}_cluster{}_all_peaks_concatenated.tsv'.format(output_folder, dataset_name, cao_optimal_n_clusters)):

#     process_list = []

#     for list_counter in range(cao_optimal_n_clusters):
#         # Bash commands to call HOMER annotatePeaks.pl to append gene annotations to each peak
#         process = 'annotatePeaks.pl {} {}/{}_cluster{}_all_peaks_concatenated.tsv {} -m {} -nmotifs -go {}/{}_cluster{}_gene_ontology > {}/{}_cluster{}_all_peaks_annotated.tsv\n\n'.format(
#             homer_annotatePeaks_arg,
#             output_folder, 
#             dataset_name, 
#             (list_counter + 1),
#             genome_ref, 
#             motif_file_full_path, 
#             output_folder, 
#             dataset_name, 
#             (list_counter + 1),
#             output_folder, 
#             dataset_name, 
#             (list_counter + 1))

#         process_list.append(process)

#     popen_list = [subprocess.Popen(process, shell = True) for process in process_list]

#     for popen in popen_list:
#         popen.wait()




# cluster_go_tsv_list = ['{}/{}_cluster{}_gene_ontology/{}.txt'.format(output_folder, dataset_name, (list_counter + 1), database) for list_counter in range(cao_optimal_n_clusters)]
# rep_go_tsv_list = ['{}/{}_rep{}_gene_ontology/{}.txt'.format(output_folder, dataset_name, (list_counter + 1), database) for list_counter in range(n_reps)]

# match_score_data_array = []

# cluster_counter = 0



# for cluster_tsv in cluster_go_tsv_list:
#     print('Opening file: {}'.format(cluster_tsv)) # Reading the HOMER annotated cluster file
#     cluster_df = pd.read_csv(cluster_tsv, delimiter = '\t')
#     cluster_array = cluster_df.values.tolist()

#     if not args.top_rank:
#         top_rankers = len(cluster_array)

#     cluster_top_terms = cluster_df['Term'].unique()[0:top_rankers]

#     top_rankers_match_score_list = []

#     rep_counter = 0

#     cluster_counter += 1



#     for rep_tsv in rep_go_tsv_list:
#         print('Opening file: {}'.format(rep_tsv)) # Reading the HOMER annotated rep file
#         rep_df = pd.read_csv(rep_tsv, delimiter = '\t')
#         rep_array = rep_df.values.tolist()

#         rep_top_terms = rep_df['Term'].unique()[0:top_rankers]

#         # spearman_rho_score = stats.spearmanr(cluster_top_terms, rep_top_terms)
#         # jaccard_micro_score = jaccard_score(cluster_top_terms, rep_top_terms, average = 'micro')
#         # jaccard_macro_score = jaccard_score(cluster_top_terms, rep_top_terms, average = 'macro')
#         # jaccard_weighted_score = jaccard_score(cluster_top_terms, rep_top_terms, average = 'weighted')
#         top_rankers_match_score = round((len(list(set(cluster_top_terms).intersection(rep_top_terms))) / len(rep_top_terms)), 2)

#         top_rankers_match_score_list.append(top_rankers_match_score)

#         # print(cluster_top_terms)
#         # print(rep_top_terms)

#         # print('Matches: ' + str(top_rankers_match_score))
#         # print('Spearman Rho Score: ' + str(spearman_rho_score.correlation))
#         # print('Jaccard Micro Score: ' + str(jaccard_micro_score))
#         # print('Jaccard Macro Score: ' + str(jaccard_macro_score))
#         # print('Jaccard Weighted Score: ' + str(jaccard_weighted_score))

#         rep_counter += 1

#         match_score_data_list = ['cluster{}'.format(cluster_counter), 'rep{}'.format(rep_counter), top_rankers_match_score]

#         match_score_data_array.append(match_score_data_list)



#     datapoint_color_list = np.array([[top_rankers_match_score, 0, 0] for top_rankers_match_score in top_rankers_match_score_list])

#     print('Separating each replicate peakset profile with Multi Dimensional Scaling')
#     transformed_hamming_distance_matrix = mds.fit_transform(hamming_distance_matrix)
#     plt.scatter(transformed_hamming_distance_matrix[:,0], transformed_hamming_distance_matrix[:,1], c = datapoint_color_list)

#     for i, txt in enumerate(n):
#         plt.annotate(txt, (transformed_hamming_distance_matrix[i,0], transformed_hamming_distance_matrix[i,1]))

#     plt.title('MDS Hamming Distance ({}) - Cluster {}'.format(database, cluster_counter))
#     plt.tight_layout()
#     plt.savefig('{}/{}_MDS_hamming_{}_cluster{}.png'.format(output_folder, dataset_name, database, cluster_counter), dpi = 300)
#     plt.clf()
#     print('Plot saved as: {}/{}_MDS_hamming_{}_cluster{}.png'.format(output_folder, dataset_name, database, cluster_counter))



# match_score_data_df = pd.DataFrame(data = match_score_data_array, columns = None)
# match_score_data_df.to_csv('{}/{}_match_score.tsv'.format(output_folder, dataset_name), sep = '\t', index = False, header = None)
# print('Match score data is saved in: {}/{}_match_score.tsv'.format(output_folder, dataset_name))



# rep_top_terms_list = []

# for rep_tsv in rep_go_tsv_list:
#     print('Opening file: {}'.format(rep_tsv)) # Reading the HOMER annotated rep file
#     rep_df = pd.read_csv(rep_tsv, delimiter = '\t')
#     # rep_array = rep_df.values.tolist()

#     rep_top_terms = rep_df['Term'].unique()
#     rep_top_terms_list.append(rep_top_terms)


# def spearmanr_calculator(list_a):
        
#     spearman_rho_score_list = []

#     for list_b in rep_top_terms_list:
#         spearman_rho_score = stats.spearmanr(list_a, list_b)
#         spearman_rho_score_list.append(spearman_rho_score.correlation)

#     return spearman_rho_score_list



# print('Generating Spearman Rho Score matrix')
# pool_spearman_rho_score = multiprocessing.Pool(processes = cpu_count)
# spearman_rho_score_matrix = pool_spearman_rho_score.map(spearmanr_calculator, rep_top_terms_list)
# pool_spearman_rho_score.close()
# pool_spearman_rho_score.join()

# print('Writing Spearman Rho Score matrix file')
# pd.DataFrame(spearman_rho_score_matrix).to_csv('{}/{}_spearman_rho_score_matrix.csv'.format(output_folder, dataset_name), sep = '\t', index = False, header = None)



# transformed_spearman_rho_score_matrix = mds.fit_transform(spearman_rho_score_matrix)

# plt.scatter(transformed_spearman_rho_score_matrix[:,0], transformed_spearman_rho_score_matrix[:,1], c = datapoint_color_list, cmap = 'rainbow')

# for i, txt in enumerate(n):
#     plt.annotate(txt, (transformed_spearman_rho_score_matrix[i,0], transformed_spearman_rho_score_matrix[i,1]))

# plt.title('MDS Spearman Rho Score ({})'.format(database))
# plt.tight_layout()
# plt.savefig('{}/{}_MDS_spearman_rho_{}.png'.format(output_folder, dataset_name, database), dpi = 300)
# plt.clf()
# print('Plot saved as: {}/{}_MDS_spearman_rho_{}.png'.format(output_folder, dataset_name, database))




# def extended_tau(list_a, list_b):
#     """ Calculate the extended Kendall tau from two lists. """
#     ranks = join_ranks(create_rank(list_a), create_rank(list_b)).fillna(len(list_a))
#     dummy_df = pd.DataFrame([{'rank_a': len(list_a), 'rank_b': len(list_b)} for i in range(len(list_a)*2-len(ranks))])
#     total_df = ranks.append(dummy_df)
#     return scale_tau(len(list_a), stats.kendalltau(total_df['rank_a'], total_df['rank_b'])[0])

# def scale_tau(length, value):
#     """ Scale an extended tau correlation such that it falls in [-1, +1]. """
#     n_0 = 2*length*(2*length-1)
#     n_a = length*(length-1)
#     n_d = n_0 - n_a
#     min_tau = (2.*n_a - n_0) / (n_d)
#     return 2*(value-min_tau)/(1-min_tau) - 1

# def create_rank(a):
#     """ Convert an ordered list to a DataFrame with ranks. """
#     return pd.DataFrame(
#                   zip(a, range(len(a))),
#                   columns=['key', 'rank'])\
#              .set_index('key')

# def join_ranks(rank_a, rank_b):
#     """ Join two rank DataFrames. """
#     return rank_a.join(rank_b, lsuffix='_a', rsuffix='_b', how='outer')



# def kendall_tau_calculator(list_a):
        
#     kendall_tau_score_list = []

#     for list_b in rep_top_terms_list:
#         kendall_tau_score = extended_tau(list_a, list_b)
#         kendall_tau_score_list.append(kendall_tau_score)

#     return kendall_tau_score_list



# print('Generating Kendall Tau Score matrix')
# pool_kendall_tau_score = multiprocessing.Pool(processes = cpu_count)
# kendall_tau_score_matrix = pool_kendall_tau_score.map(kendall_tau_calculator, rep_top_terms_list)
# pool_kendall_tau_score.close()
# pool_kendall_tau_score.join()

# print('Writing Kendall Tau Score matrix file')
# pd.DataFrame(kendall_tau_score_matrix).to_csv('{}/{}_kendall_tau_score_matrix.csv'.format(output_folder, dataset_name), sep = '\t', index = False, header = None)



# transformed_kendall_tau_score_matrix = mds.fit_transform(kendall_tau_score_matrix)

# plt.scatter(transformed_kendall_tau_score_matrix[:,0], transformed_kendall_tau_score_matrix[:,1], c = datapoint_color_list, cmap = 'rainbow')

# for i, txt in enumerate(n):
#     plt.annotate(txt, (transformed_kendall_tau_score_matrix[i,0], transformed_kendall_tau_score_matrix[i,1]))

# plt.title('MDS Kendall Tau Score ({})'.format(database))
# plt.tight_layout()
# plt.savefig('{}/{}_MDS_kendall_tau_{}.png'.format(output_folder, dataset_name, database), dpi = 300)
# plt.clf()
# print('Plot saved as: {}/{}_MDS_kendall_tau_{}.png'.format(output_folder, dataset_name, database))



print('Done! Exiting script!')