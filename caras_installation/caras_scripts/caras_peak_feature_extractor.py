#!/usr/bin/env python3
#pyright: reportUnboundVariable=false

print('Importing required modules')

import pandas as pd
import subprocess
import multiprocessing
import os
import sys
import argparse
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.spatial.distance import hamming
from scipy.spatial.distance import jaccard
from sklearn.manifold import MDS
import rbo
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csgraph
from numpy import linalg as LA
from sklearn.cluster import SpectralClustering, AffinityPropagation


def removesuffix(s, suf):
    if suf and s.endswith(suf):
        return s[:-len(suf)]
    return s


def remove_constant_prefix_suffix(string_list):
    
    unique_substring_list = []

    character_counter = 0
    prefix_index_found = 0

    while True:
        current_character = string_list[0][character_counter]

        for string in string_list:
            if string[character_counter] != current_character:
                prefix_index = character_counter
                prefix_index_found = 1
                break

        if prefix_index_found == 1:
            break

        character_counter += 1

    character_counter = -1
    suffix_index_found = 0

    while True:
        current_character = string_list[0][character_counter]

        for string in string_list:
            if string[character_counter] != current_character:
                suffix_index = character_counter
                suffix_index_found = 1
                break

        if suffix_index_found == 1:
            break

        character_counter -= 1

    for string in string_list:
        unique_substring_list.append(string[prefix_index : suffix_index + 1])

    return unique_substring_list


def getAffinityMatrix(coordinates, k = 10):
    """
    Calculate affinity matrix based on input coordinates matrix and the number
    of nearest neighbours.
    
    Apply local scaling based on the k nearest neighbour
        References:
    https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf
    """
    # calculate euclidian distance matrix
    dists = squareform(pdist(coordinates)) 
    
    # for each row, sort the distances ascendingly and take the index of the 
    #k-th position (nearest neighbour)
    knn_distances = np.sort(dists, axis=0)[k]
    knn_distances = knn_distances[np.newaxis].T
    
    # calculate sigma_i * sigma_j
    local_scale = knn_distances.dot(knn_distances.T)

    affinity_matrix = dists * dists
    affinity_matrix = -affinity_matrix / local_scale
    # divide square distance matrix by local scale
    affinity_matrix[np.where(np.isnan(affinity_matrix))] = 0.0
    # apply exponential
    affinity_matrix = np.exp(affinity_matrix)
    np.fill_diagonal(affinity_matrix, 0)
    return affinity_matrix


def eigenDecomposition(A, plot = True, topK = 5):
    """
    :param A: Affinity matrix
    :param plot: plots the sorted eigen values for visual inspection
    :return A tuple containing:
    - the optimal number of clusters by eigengap heuristic
    - all eigen values
    - all eigen vectors
    
    This method performs the eigen decomposition on a given affinity matrix,
    following the steps recommended in the paper:
    1. Construct the normalized affinity matrix: L = D-1/2AD^ -1/2.
    2. Find the eigenvalues and their associated eigen vectors
    3. Identify the maximum gap which corresponds to the number of clusters
    by eigengap heuristic
    
    References:
    https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf
    http://www.kyb.mpg.de/fileadmin/user_upload/files/publications/attachments/Luxburg07_tutorial_4488%5b0%5d.pdf
    """
    L = csgraph.laplacian(A, normed=True)
    n_components = A.shape[0]
    
    # LM parameter : Eigenvalues with largest magnitude (eigs, eigsh), that is, largest eigenvalues in 
    # the euclidean norm of complex numbers.
    # eigenvalues, eigenvectors = eigsh(L, k=n_components, which="LM", sigma=1.0, maxiter=5000)
    eigenvalues, eigenvectors = LA.eig(L)
    
    if plot:
        plt.title('Largest eigen values of input matrix')
        plt.scatter(np.arange(len(eigenvalues)), eigenvalues)
        plt.grid()
        
    # Identify the optimal number of clusters as the index corresponding
    # to the larger gap between eigen values
    index_largest_gap = np.argsort(np.diff(eigenvalues))[::-1][:topK]
    nb_clusters = index_largest_gap + 1
        
    return nb_clusters, eigenvalues, eigenvectors


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

parser.add_argument('--annotate_skip', 
                    help = '<Optional> Use to skip annotatePeaks.pl run (when rerunning this script, optimizing clusters with the same existing --filter)',
                    action = 'store_true',
                    required = False)

parser.add_argument('--ref', 
                    help = '<Optional> Your sample organism genome reference build. Default is hg38 (human).', 
                    choices = ['hg19', 'hg38', 'mm9', 'mm10', 'dm6', 'sacCer3'], 
                    default = 'hg38')

parser.add_argument('--filter', 
                    help = '<Required> The peakset to process. Accepted values are the (space-separated) names of the peak caller (case-sensitive) to select all peaks called by the given peak caller (e.g., MACS2, Genrich, etc). Use - for AND combination (e.g. GEM-HOMER). Use + for OR combination (e.g. SICER2+SEACR). Can also select peaks based on number of peak callers by using "Xcaller" where X is the minimum number of peak callers.',
                    nargs = '+',
                    required = True)

parser.add_argument('--repnum', 
                    help = '<Required> Number of replicates or samples', 
                    type = int,
                    required = True)

parser.add_argument('--clustnum', 
                    help = '<Optional> User-determined expected number of resulting clusters', 
                    type = int,
                    default = 0)

parser.add_argument('--distance', 
                    help = '<Required> Minimum distance between peaks in different replicates to be considered as peaks in different locations.', 
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
                    help = '<Required> HOMER GO or pathway database to use. Multiple databases can be used in one go (space-separated).',
                    nargs = '+',
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

minimum_positive_fraction = args.min_pos

maximum_positive_fraction = args.max_pos

interpeaks_distance = int(args.distance)

if args.top_rank:
    top_rankers = args.top_rank

database_list = [database for database in args.database]

peakset_filter_list = [peakset_filter for peakset_filter in args.filter]

new_peakset_filter_list = []

for peakset_filter in peakset_filter_list:   
    if peakset_filter.isdigit():
        peakset_filter = '{}caller'.format(peakset_filter)

    new_peakset_filter_list.append(peakset_filter)

peakset_filter_list = new_peakset_filter_list


n_reps = args.repnum

expected_cluster =  args.clustnum

feature_list = []

for replicate_counter in range(n_reps):
    for peakset_filter in peakset_filter_list:
        feature_list.append('{}_{}'.format(peakset_filter, replicate_counter + 1))

peakset_filter_string = '_'.join([peakset_filter for peakset_filter in peakset_filter_list])



# cell_line_list = ['H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 
#     'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562']
# color_list = ['red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red',
#     'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue']
# true_label_list = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

# cell_line_list = ['CTCF', 'CTCF', 'CTCF', 'CTCF', 'CTCF', 'CTCF', 'CTCF',
#     'H3K27me3', 'H3K27me3', 'H3K27me3', 'H3K27me3', 'H3K27me3', 'H3K27me3', 'H3K27me3']
# color_list = ['red', 'red', 'red', 'red', 'red', 'red', 'red',
#     'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue']
# true_label_list = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]

# cell_line_list = ['H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 
#     'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562', 'K562',
#     'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1', 'H1']
# color_list = ['red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red',
#     'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue',
#     'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red']

# selected_cell_line_list = cell_line_list

peakwise_feature_df = None
replicatewise_feature_df = None
filtered_replicatewise_feature_df = None



# ########## TERM RANKING ##########


def called_term_ranking_function(term_ranking_file):
    
    print(term_ranking_file)
    print('Reading {} term ranking from sample {}'.format(term_ranking_file.split('/')[-1], term_ranking_file.split('/')[-2])) # Reading the GO term ranking file
    
    term_ranking_df = pd.read_csv(term_ranking_file, delimiter = '\t')
    sample_top_terms = term_ranking_df['Term'].unique()[0:top_rankers]
    
    return sample_top_terms


if args.database:

    try:
        os.mkdir('{}/clustering_files'.format(output_folder))
    except FileExistsError:
        pass

    replicate_number_list = remove_constant_prefix_suffix(peaklist_tsv_full_path_list)

    for peaklist_counter, peaklist_tsv_full_path in enumerate(peaklist_tsv_full_path_list):

        print('Opening file: {}'.format(peaklist_tsv_full_path)) # Reading the HOMER annotated interval file
        replicate_number = replicate_number_list[peaklist_counter]
        peaklist_df = pd.read_csv(peaklist_tsv_full_path, delimiter = '\t')
        peaklist_header = peaklist_df.columns.tolist()
        peaklist_array = peaklist_df.values.tolist()

        for peakset_filter in peakset_filter_list:
            
            filtered_peaklist_array = []

            for peaklist_array_row in peaklist_array:
                peak_caller_list = peaklist_array_row[5].split('|')

                if 'caller' in peakset_filter:
                    feature_caller_number = peakset_filter.rstrip('caller')
                    if len(peak_caller_list) >= int(feature_caller_number):
                        filtered_peaklist_array.append(peaklist_array_row[:5])
                elif '-' in peakset_filter:
                    feature_element_list = peakset_filter.split('-')
                    if all(feature_element in peak_caller_list for feature_element in feature_element_list):
                        filtered_peaklist_array.append(peaklist_array_row[:5])
                elif '+' in peakset_filter:
                    feature_element_list = peakset_filter.split('+')
                    if any(feature_element in peak_caller_list for feature_element in feature_element_list):
                        filtered_peaklist_array.append(peaklist_array_row[:5])
                else:
                    if peakset_filter in peak_caller_list:
                        filtered_peaklist_array.append(peaklist_array_row[:5])

            filtered_peaklist_df = pd.DataFrame(filtered_peaklist_array, columns = peaklist_header[:5])

            filtered_peaklist_df.to_csv('{}/clustering_files/{}_rep{}_{}.tsv'.format(output_folder, dataset_name, replicate_number, peakset_filter), sep = '\t', index = None)

            if not args.annotate_skip:
                subprocess.run('annotatePeaks.pl {}/clustering_files/{}_rep{}_{}.tsv hg38 -go {}/clustering_files/{}_rep{}_{}_gene_ontology > /dev/null'.format(
                    output_folder, 
                    dataset_name, 
                    replicate_number, 
                    peakset_filter, 
                    output_folder, 
                    dataset_name, 
                    replicate_number, 
                    peakset_filter), 
                    shell = True)
        
    for peakset_filter in peakset_filter_list:
        for current_database in database_list:
            term_ranking_file_list = ['{}/clustering_files/{}_rep{}_{}_gene_ontology/{}.txt'.format(output_folder, dataset_name, replicate_number, peakset_filter, current_database) for replicate_number in replicate_number_list]

            # Rank the term according to each sample's enrichment-based term list
            print('Calculating peak ranks based on {} database'.format(current_database))
            pool_enriched_term_ranking = multiprocessing.Pool(processes = cpu_count)
            ranked_term_array = pool_enriched_term_ranking.map(called_term_ranking_function, term_ranking_file_list)
            pool_enriched_term_ranking.close()
            pool_enriched_term_ranking.join()
        

            def term_ranking_rbo_calculator_function(ranked_term_list_a):
                
                term_ranking_rbo_distance_list = []

                for counter, ranked_term_list_b in enumerate(ranked_term_array): # Iterate through all other samples' ranked term list
                    term_ranking_rbo_similarity = rbo.RankingSimilarity(ranked_term_list_a, ranked_term_list_b).rbo() # Calculate term ranking RBO similarity
                    term_ranking_rbo_distance = (term_ranking_rbo_similarity - 1) * -1 # Calculate term ranking RBO distance
                    term_ranking_rbo_distance_list.append(term_ranking_rbo_distance) # Append to the list

                return term_ranking_rbo_distance_list # Return the distance matrix row
            

            # Generate RBO ranking distance matrix based on term ranking of each replicate
            print('Generating RBO term ranking distance matrix. Peakset: {}. Database: {}'.format(peakset_filter, current_database))
            pool_term_ranking_rbo_distance = multiprocessing.Pool(processes = cpu_count)
            term_ranking_rbo_distance_matrix = pool_term_ranking_rbo_distance.map(term_ranking_rbo_calculator_function, ranked_term_array)
            pool_term_ranking_rbo_distance.close()
            pool_term_ranking_rbo_distance.join()


            print('Transforming term ranking RBO distance matrix with Multi Dimensional Scaling. Peakset: {}. Database: {}'.format(peakset_filter, current_database))
            mds = MDS(dissimilarity = 'precomputed', random_state = 0)
            transformed_distance_matrix = mds.fit_transform(term_ranking_rbo_distance_matrix)

            if expected_cluster != 0:
                n_clust = expected_cluster
                print('Peakset: {}. Database: {}. User-determined expected number of clusters: {}'.format(peakset_filter, current_database, n_clust))

            else:
                affinity_matrix = getAffinityMatrix(transformed_distance_matrix, k = 5)
                optimal_cluster_list, _,  _ = eigenDecomposition(affinity_matrix)
                print('Peakset: {}. Database: {}. Number of clusters with largest EigenGap: {}'.format(peakset_filter, current_database, optimal_cluster_list))

                for optimal_cluster in optimal_cluster_list:
                    if optimal_cluster <= 10:
                        n_clust = optimal_cluster
                        print('Peakset: {}. Database: {}. Optimal number of clusters: {}'.format(peakset_filter, current_database, n_clust))
                        break

            clustering = SpectralClustering(n_clusters = n_clust, assign_labels = "discretize", random_state = 0).fit(transformed_distance_matrix)
            predicted_label_list = clustering.labels_
            predicted_label_list = [x + 1 for x in predicted_label_list]
            plt.figure(figsize = (6,6))
            # plt.subplot(121)
            plt.title('Spectral Clustering Prediction (Peakset: {}. Database: {})'.format(peakset_filter, current_database))
            plt.scatter(transformed_distance_matrix[:, 0], transformed_distance_matrix[:, 1], s = 50, c = predicted_label_list);

            # plt.subplot(122)
            # plt.title('Ground Truth (Peakset: {}. Database: {})'.format(peakset_filter, current_database))
            # plt.scatter(transformed_distance_matrix[:, 0], transformed_distance_matrix[:, 1], s = 50, c = true_label_list);
            
            plt.tight_layout()
            plt.savefig('{}/{}_{}_{}_Term_Ranking_RBO_Distance_MDS_Spectral_Clustering.png'.format(output_folder, dataset_name, peakset_filter, current_database), dpi = 300)
            plt.clf()

            cluster_group_tsv = open('{}/{}_{}_{}_Term_Ranking_RBO_Distance_MDS_Spectral_Clustering.tsv'.format(output_folder, dataset_name, peakset_filter, current_database), 'w')
            for label_counter, predicted_label in enumerate(predicted_label_list):
                cluster_group_tsv.write(peaklist_tsv_full_path_list[label_counter].split('/')[-1] + '\t' + str(predicted_label) + '\n')
            cluster_group_tsv.close()

            unique_label_list = np.unique(predicted_label_list)

            for unique_label in unique_label_list:
                group_file_list = []
                for label_counter, predicted_label in enumerate(predicted_label_list):
                    if predicted_label == unique_label:
                        group_file_list.append(peaklist_tsv_full_path_list[label_counter])
                        print('Concatenating {} into cluster {} peaklist'.format(peaklist_tsv_full_path_list[label_counter].split('/')[-1], unique_label))
                group_file_string = ' '.join(group_file_list)

                subprocess.run('cat {} | cut -f 1-7 | grep -v "Peak ID" | sort -k 1 -V > {}/{}_{}_{}_group{}_all_peaks_clustered.tsv'.format(group_file_string, output_folder, dataset_name, peakset_filter, current_database, unique_label), shell = True)

                if sys.platform == 'linux' or sys.platform == 'linux2':
                    # subprocess.run("sed -i '1i Peak ID\tChr\tStart\tEnd\tStrand\tPeak Caller Combination\tPeak Caller Overlaps\tChIP Tag Count\tControl Tag Count\tFold Change\tPeak Center\tNumber of Motifs\tAnnotation\tDetailed Annotation\tDistance to TSS\tNearest PromoterID\tEntrez ID\tNearest Unigene\tNearest Refseq\tNearest Ensembl\tGene Name\tGene Alias\tGene Description\tGene Type\tCpG%\tGC%' {}/{}_{}_{}_group{}_all_peaks_clustered.tsv".format(output_folder, dataset_name, peakset_filter, current_database, unique_label), shell = True)
                    subprocess.run("sed -i '1i Peak ID\tChr\tStart\tEnd\tStrand\tPeak Caller Combination\tPeak Caller Overlaps' {}/{}_{}_{}_group{}_all_peaks_clustered.tsv".format(output_folder, dataset_name, peakset_filter, current_database, unique_label), shell = True)

                if sys.platform == 'darwin':
                    # subprocess.run("gsed -i '1i Peak ID\tChr\tStart\tEnd\tStrand\tPeak Caller Combination\tPeak Caller Overlaps\tChIP Tag Count\tControl Tag Count\tFold Change\tPeak Center\tNumber of Motifs\tAnnotation\tDetailed Annotation\tDistance to TSS\tNearest PromoterID\tEntrez ID\tNearest Unigene\tNearest Refseq\tNearest Ensembl\tGene Name\tGene Alias\tGene Description\tGene Type\tCpG%\tGC%' {}/{}_{}_{}_group{}_all_peaks_clustered.tsv".format(output_folder, dataset_name, peakset_filter, current_database, unique_label), shell = True)
                    subprocess.run("gsed -i '1i Peak ID\tChr\tStart\tEnd\tStrand\tPeak Caller Combination\tPeak Caller Overlaps' {}/{}_{}_{}_group{}_all_peaks_clustered.tsv".format(output_folder, dataset_name, peakset_filter, current_database, unique_label), shell = True)

                subprocess.run('annotatePeaks.pl {}/{}_{}_{}_group{}_all_peaks_clustered.tsv hg38 -go {}/{}_{}_{}_group{}_gene_ontology > /dev/null'.format(
                    output_folder, 
                    dataset_name, 
                    peakset_filter, 
                    current_database, 
                    unique_label,
                    output_folder, 
                    dataset_name, 
                    peakset_filter, 
                    current_database, 
                    unique_label), 
                    shell = True)
                    
print('Done! Exiting script!')



# ########## PEAK RANKING ##########


# if args.peaklist_tsv:

#     print('Opening file: {}'.format(interval_bed_full_path)) # Reading the HOMER annotated interval file
#     interval_df = pd.read_csv(interval_bed_full_path, delimiter = '\t', header = None)
#     interval_df = interval_df.iloc[:,[0,1,2]] # Get only the chr, start, end
#     interval_df.sort_values(by = [0,1], inplace = True) # Sorting the intervals list based on chromosomal position
#     interval_array = interval_df.values.tolist()

#     peaklist_df_list = []

#     for peaklist_tsv_full_path in peaklist_tsv_full_path_list:
#         print('Opening file: {}'.format(peaklist_tsv_full_path))
#         peaklist_df = pd.read_csv(peaklist_tsv_full_path, delimiter = '\t')
#         peaklist_df.sort_values(by = ['Fold Change', 'ChIP Tag Count'], ascending = False, inplace = True) # Sort the peaklist. Most confident on top.
#         # peaklist_df[peaklist_df.columns[:(list(peaklist_df.columns).index('IDR'))+1]]
#         peaklist_df = peaklist_df[['Chr', 'Start', 'End', 'Peak Caller Overlaps', 'Peak Caller Combination']] # Remove unused columns to save space
#         peaklist_df_list.append(peaklist_df) # Dataframe list of all peaklists (all peaks calculated from each replicate)

#     for filter_counter in range(len(peakset_filter_list)):        
#         current_filter_string = peakset_filter_list[filter_counter]
#         print('Processing peaklists based on filter: {}'.format(current_filter_string))

#         filtered_peaklist_df_list = []

#         if 'caller' in current_filter_string: # Select only peaks called by minimum certain number of peak callers
#             for peaklist_df in peaklist_df_list:
#                 filtered_peaklist_df_list.append(peaklist_df[peaklist_df['Peak Caller Overlaps'] >= int(current_filter_string.rstrip('caller'))])
                
#         elif '-' in current_filter_string: # Select only the consensus peaks between peak caller sets of user choice
#             for peaklist_df in peaklist_df_list:
#                 filtered_peaklist_df_list.append(peaklist_df[peaklist_df['Peak Caller Combination'].apply(lambda x: all(item[0] in x.split('|') for item in current_filter_string.split('-')))])

#         elif '+' in current_filter_string: # Select only the union peaks between peak caller sets of user choice
#             for peaklist_df in peaklist_df_list:
#                 filtered_peaklist_df_list.append(peaklist_df[peaklist_df['Peak Caller Combination'].apply(lambda x: any(item[0] in x.split('|') for item in current_filter_string.split('+')))])

#         else:
#             for peaklist_df in peaklist_df_list: # Select only one peak set of user choice
#                 filtered_peaklist_df_list.append(peaklist_df[peaklist_df['Peak Caller Combination'].apply(lambda x: current_filter_string in x)])


#         def called_peak_ranking_function(peaklist_df):

#             chr_list = peaklist_df['Chr'].tolist()
#             start_list = peaklist_df['Start'].tolist()
#             end_list = peaklist_df['End'].tolist()

#             ranked_peak_list = []

#             for row_counter in range(peaklist_df.shape[0]): # Iterate through each replicate's peak calculated list
#                 for interval_array_row in interval_array: # Compare to the combined master peak interval list, looking for overlaps (including interpeaks distance)
#                     if interval_array_row[0] == chr_list[row_counter]:
#                         if int(interval_array_row[1]) - interpeaks_distance < int(end_list[row_counter]):
#                             if int(interval_array_row[2]) + interpeaks_distance > int(start_list[row_counter]):
#                                 current_peak_ID = '{}:{}-{}'.format(interval_array_row[0], interval_array_row[1], interval_array_row[2]) # For every overlap, record the master peak ID
#                                 if current_peak_ID not in ranked_peak_list: # Multiple peaks in the peak calculated list may overlap the same master peak interval
#                                     ranked_peak_list.append(current_peak_ID) # Append only if not pre-existing within the list (only one record per peak ID)

#                                 if len(ranked_peak_list) >= 10000: # Record at maximum 100 top peaks only
#                                     return ranked_peak_list

#             return ranked_peak_list # Return anyway even if the list does not reach 100 peaks


#         # Rank the peak ID according to each sample's confidence-ranked peak list
#         print('Calculating peak ranks')
#         pool_called_peak_ranking = multiprocessing.Pool(processes = cpu_count)
#         ranked_peak_array = pool_called_peak_ranking.map(called_peak_ranking_function, filtered_peaklist_df_list)
#         pool_called_peak_ranking.close()
#         pool_called_peak_ranking.join()


#         def peak_ranking_rbo_calculator_function(ranked_peak_list_a):
            
#             peak_ranking_rbo_distance_list = []

#             for counter, ranked_peak_list_b in enumerate(ranked_peak_array): # Iterate through all other samples' ranked peak list
#                 peak_ranking_rbo_similarity = rbo.RankingSimilarity(ranked_peak_list_a, ranked_peak_list_b).rbo() # Calculate peak ranking RBO similarity
#                 peak_ranking_rbo_distance = (peak_ranking_rbo_similarity - 1) * -1 # Calculate peak ranking RBO distance
#                 peak_ranking_rbo_distance_list.append(peak_ranking_rbo_distance) # Append to the list

#             return peak_ranking_rbo_distance_list # Return the distance matrix row


#         # Generate peak ranking RBO distance matrix based on peak ID ranking of each replicate
#         print('Generating peak ranking RBO distance matrix')
#         pool_peak_ranking_rbo_distance = multiprocessing.Pool(processes = cpu_count)
#         peak_ranking_rbo_distance_matrix = pool_peak_ranking_rbo_distance.map(peak_ranking_rbo_calculator_function, ranked_peak_array)
#         pool_peak_ranking_rbo_distance.close()
#         pool_peak_ranking_rbo_distance.join()


#         print('Transforming peak ranking RBO distance matrix with Multi Dimensional Scaling')
#         mds = MDS(dissimilarity = 'precomputed', random_state = 0)
#         transformed_distance_matrix = mds.fit_transform(peak_ranking_rbo_distance_matrix)

#         affinity_matrix = getAffinityMatrix(transformed_distance_matrix, k = 5)
#         number_of_cluster_list, _,  _ = eigenDecomposition(affinity_matrix)
#         print('Number of clusters with largest EigenGap: {}'.format(number_of_cluster_list))

#         for number_of_cluster in number_of_cluster_list:
#             if number_of_cluster <= 10:
#                 optimal_cluster = number_of_cluster
#                 print('Optimal number of clusters: {}'.format(optimal_cluster))
#                 break

#         clustering = SpectralClustering(n_clusters = optimal_cluster, assign_labels = "discretize", random_state = 0).fit(transformed_distance_matrix)
#         predicted_label_list = clustering.labels_
#         plt.figure(figsize = (14,6))
#         plt.subplot(121)
#         plt.title('Spectral Clustering Prediction (Based on {} Peak Ranking at minimum interpeaks distance = {})'.format(current_filter_string, interpeaks_distance))
#         plt.scatter(transformed_distance_matrix[:, 0], transformed_distance_matrix[:, 1], s = 50, c = predicted_label_list);

#         plt.subplot(122)
#         plt.title('Ground Truth (Based on {} Peak Ranking at minimum interpeaks distance = {})'.format(current_filter_string, interpeaks_distance))
#         plt.scatter(transformed_distance_matrix[:, 0], transformed_distance_matrix[:, 1], s = 50, c = true_label_list);
        
#         plt.tight_layout()
#         plt.savefig('{}/{}_{}_d{}_Peak_Ranking_RBO_Distance_MDS_Spectral_Clustering.png'.format(output_folder, dataset_name, current_filter_string, interpeaks_distance), dpi = 300)
#         plt.clf()