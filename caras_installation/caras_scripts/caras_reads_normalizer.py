#!/usr/bin/env python3
#pyright: reportUnboundVariable=false

script_version = '1.0'

import subprocess
import os
import multiprocessing
import argparse
import sys
import math
import csv
import numpy as np
import pandas as pd
from datetime import datetime

subprocess.run('ulimit -n 2000', shell = True)
subprocess.run('ulimit -s 65536', shell = True)




# def normalized_bam_generator(arg_1, arg_2):

#     with open('{}/{}_{}.normalized.sam'.format(bam_file_dir_list[bam_file_counter], bam_file_name_list[bam_file_counter], arg_1), 'w') as normalized_sam_file:
#         normalized_sam_file = csv.writer(normalized_sam_file, delimiter = '\t')
#         normalized_sam_file.writerows(sam_header_array)

#         for normalized_sam_body_list in arg_2:
#             normalized_sam_file.writerows(normalized_sam_body_list)

#     subprocess.run('samtools view -h -bS {}/{}_{}.normalized.sam > {}/{}_{}.normalized.bam'.format(
#         bam_file_dir_list[bam_file_counter],
#         bam_file_name_list[bam_file_counter],
#         arg_1,
#         bam_file_dir_list[bam_file_counter],
#         bam_file_name_list[bam_file_counter],
#         arg_1),
#         shell = True)

#     subprocess.run('rm -r -f {}/{}_{}.normalized.sam'.format(
#         bam_file_dir_list[bam_file_counter],
#         bam_file_name_list[bam_file_counter],
#         arg_1),
#         shell = True)

    # print('Normalized partial bam file is saved as {}/{}_{}.normalized.bam'.format(bam_file_dir_list[bam_file_counter], bam_file_name_list[bam_file_counter], arg_1))



# def sam_numerator(sam_array_row):
#     normalized_sam_list = []
#     for counter in range(int(sam_array_row[0])):
#         normalized_sam_list.append(sam_array_row[1])

#     return normalized_sam_list



parser = argparse.ArgumentParser()

parser.add_argument('--thread', 
                    help = '<Optional> Maximum number of threads/processes to use. Default is half the maximum available.', 
                    type = int)

parser.add_argument('--norm', 
                    help = '<Optional> Normalize bedgraph signal based on the number of specific reads mapped to a spiked-in DNA genome. Default is no normalization.',
                    choices = ['mapped', 'paired', 'properly_paired', 'value', 'none'],
                    default = 'none')

parser.add_argument('--bam', 
                    nargs = '+', 
                    help = '<Required> <file> bam file of reads aligned to ChIP-seq sample organism genome, to be converted into bedgraph and bigwig files.',
                    required = True)

parser.add_argument('--norm_bam', 
                    nargs = '+', 
                    help = '<Optional> bam file of reads aligned to spiked-in or carry-over DNA genome, to be used as signal normalizer for the output files. OR the path to greenlist.')

parser.add_argument('--make_bam', 
                    help = '<Optional> Switch flag to generate bam file.', 
                    action = 'store_true')

parser.add_argument('--make_bdg', 
                    help = '<Optional> Switch flag to generate bedgraph file.', 
                    action = 'store_true')

parser.add_argument('--make_bw', 
                    help = '<Optional> Switch flag to generate bigwig file.', 
                    action = 'store_true')

parser.add_argument('--bin', 
                    help = '<Optional> Bin size of the output bedgraph and bigwig file. Default is 10 base pairs.', 
                    type = int,
                    default = 10)

parser.add_argument('--num', 
                    help = '<Optional> Numerator value for the scaling factor. Sometimes spiked-in DNA can be too high in number, causing the generated bedgraph signal amplitude to be extremely low. Increase the default value of this number to deal with such case, or use the default automatic mode (0).', 
                    type = int,
                    default = 0)

parser.add_argument('--log_dir', 
                    help = '<Optional> Directory to save all processes log files.')

args = parser.parse_args()



###=====================================================================================================================



if args.thread:
    cpu_count = args.thread
elif not args.thread:    
    cpu_count = multiprocessing.cpu_count() / 2

bam_file_list = [os.path.abspath(bam_file) for bam_file in args.bam]
# bam_file = os.path.abspath(args.bam)
bam_file_extension_list = [bam_file.split('.')[-1] for bam_file in bam_file_list]

if any(bam_file_extension != 'bam' for bam_file_extension in bam_file_extension_list):
    print('Please use only .bam file inputs. Exiting program.')
    exit()

bam_file_name_list = [bam_file.split('/')[-1].strip('.bam') for bam_file in bam_file_list]
# bam_file_name = bam_file.split('/')[-1].strip('.' + bam_file_extension)
bam_file_dir_list = ['/'.join(bam_file.split('/')[:-1]) for bam_file in bam_file_list]

if any(bam_file_dir != bam_file_dir_list[0] for bam_file_dir in bam_file_dir_list):
    print('Please make sure all bam files are located in the same directory. Exiting program.')
    exit()

norm_mode = args.norm

if args.norm_bam:
    if norm_mode in ['mapped', 'paired', 'properly_paired']:
        norm_bam_file_list = [os.path.abspath(norm_bam_file) for norm_bam_file in args.norm_bam]
        norm_bam_file_extension_list = [norm_bam_file.split('.')[-1] for norm_bam_file in norm_bam_file_list]

        if any(norm_bam_file_extension != 'bam' for norm_bam_file_extension in norm_bam_file_extension_list):
            print('Please use only file inputs with .bam extension. Exiting program.')
            exit()

        if len(bam_file_list) != len(norm_bam_file_list):
            print('Number of bam files need to be the same as normalizer bam files. Exiting program')
            exit()

    if norm_mode == 'value':
        norm_bam_file = os.path.abspath(args.norm_bam[0])

if args.log_dir:
    log_dir = args.log_dir

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

norm_scale_numerator = args.num

bdg_bin_size = args.bin



###=====================================================================================================================



if norm_mode == 'value' and len(bam_file_list) > 1 and args.norm_bam:

    norm_scale_factor_list = []

    dir_file_list = os.listdir(bam_file_dir_list[0])
    
    for dir_file in dir_file_list:
       if '_greenlist_size_factors.tsv' in dir_file:
           greenlist_file = dir_file
           break
    
    greenlist_table_df = pd.read_csv('{}/{}'.format(bam_file_dir_list[0], greenlist_file), delimiter = '\t')
    greenlist_table_df = greenlist_table_df.replace(np.nan, '')
    
    greenlist_table_norm_scale_factor_list = greenlist_table_df.iloc[:,-1].values.tolist()
    greenlist_table_bam_file_list = greenlist_table_df.index.tolist()

    for bam_file in bam_file_list:
        for greenlist_table_bam_file_counter in range(len(greenlist_table_bam_file_list)):
            if greenlist_table_bam_file_list[greenlist_table_bam_file_counter] in bam_file:
                norm_scale_factor_list.append(greenlist_table_norm_scale_factor_list[greenlist_table_bam_file_counter])

    norm_scale_factor_min = np.min(norm_scale_factor_list)

    for norm_scale_factor_counter, norm_scale_factor in enumerate(norm_scale_factor_list):
        norm_scale_factor_list[norm_scale_factor_counter] = norm_scale_factor / norm_scale_factor_min

    norm_scale_numerator = 'NA'
    norm_read_count_list = ['NA' for bam_file_counter in range(len(bam_file_list))]



elif norm_mode in ['mapped', 'paired', 'properly_paired'] and len(bam_file_list) > 1 and args.norm_bam:

    norm_read_count_list = []

    for bam_file_counter in range(len(bam_file_list)): 

        if norm_mode == 'mapped':
            popen_norm = subprocess.Popen('samtools view -F4 -c {} -@ {}'.format(norm_bam_file_list[bam_file_counter], cpu_count), shell = True, stdout = subprocess.PIPE)
            # Get the number of mapped reads in the ChIP .bam file
            norm_out = popen_norm.communicate()[0]
            norm_read_count = int(norm_out.strip())

        if norm_mode == 'paired':
            popen_norm = subprocess.Popen('samtools view -f1 -c {} -@ {}'.format(norm_bam_file_list[bam_file_counter], cpu_count), shell = True, stdout = subprocess.PIPE)
            # Get the number of paired reads in the ChIP .bam file
            norm_out = popen_norm.communicate()[0]
            norm_read_count = int(norm_out.strip())

        if norm_mode == 'properly_paired':
            popen_norm = subprocess.Popen('samtools view -f2 -c {} -@ {}'.format(norm_bam_file_list[bam_file_counter], cpu_count), shell = True, stdout = subprocess.PIPE)
            # Get the number of properly_paired reads in the ChIP .bam file
            norm_out = popen_norm.communicate()[0]
            norm_read_count = int(norm_out.strip())

        norm_read_count_list.append(norm_read_count)

        current_sample_name = norm_bam_file_list[bam_file_counter].split('/')[-1].rstrip('.normmapped.bam')

        print('{} reads are {} to normalizer (spiked-in DNA / carry-over DNA) genome in sample {}'.format(norm_read_count, norm_mode, current_sample_name))



    norm_read_count_max = np.unique(norm_read_count_list)[-1]
    try:
        norm_read_count_second_max = np.max()[-2]
    except:
        norm_read_count_second_max = norm_read_count_max



    if norm_scale_numerator == 0: # norm_scale_numerator = 0 is default mode (automatic), where the script will determine the numerator value based on mapped carry over DNA fragments.
        if norm_read_count_second_max != norm_read_count_max:
            norm_scale_multiplier = round(1 / ((norm_read_count_max / norm_read_count_second_max) - 1))
        else:
            norm_scale_multiplier = 1
            
        if norm_scale_multiplier > 10:
            norm_scale_multiplier = 10
            
        if norm_scale_multiplier < 1:
            norm_scale_multiplier = 1

        norm_scale_numerator = norm_scale_multiplier * norm_read_count_max



    norm_scale_factor_list = []

    for bam_file_counter, norm_read_count in enumerate(norm_read_count_list): 

        if norm_read_count != 0:
            norm_scale_factor = round((norm_scale_numerator/norm_read_count), 4)
        else:
            norm_scale_factor = round(norm_scale_numerator, 4)
        
        if norm_scale_factor < 1:
            norm_scale_factor = 1

        if norm_scale_factor > 100:
            norm_scale_factor = 100

        norm_scale_factor_list.append(norm_scale_factor)

else:
    norm_scale_numerator = 'NA'
    norm_read_count_list = ['NA' for bam_file_counter in range(len(bam_file_list))]
    norm_scale_factor_list = [1 for bam_file_counter in range(len(bam_file_list))]



print('bam_file_names: \n{}'.format(bam_file_name_list))
print('norm_scale_numerator: \n{}'.format(norm_scale_numerator))
print('norm_read_counts: \n{}'.format(norm_read_count_list))
print('norm_scale_factors: \n{}'.format(norm_scale_factor_list))



if norm_mode in ['mapped', 'paired', 'properly_paired']:
    now = datetime.now()
    dt_string = now.strftime('%d_%m_%Y_%H_%M_%S')

    with open('{}/scaling_numbers_{}.tsv'.format(bam_file_dir_list[0], dt_string), 'w') as scaling_numbers_file:
        scaling_numbers_file = csv.writer(scaling_numbers_file, delimiter = '\t')
        scaling_numbers_file.writerow(['bam_file_name', 'norm_scale_numerator', 'norm_read_count', 'norm_scale_factor'])

        for bam_file_counter in range(len(bam_file_list)):
            scaling_numbers_file.writerow([bam_file_name_list[bam_file_counter], norm_scale_numerator, norm_read_count_list[bam_file_counter], norm_scale_factor_list[bam_file_counter]])

    print('Scaling factors for all samples are recorded in: {}/scaling_numbers_{}.tsv'.format(bam_file_dir_list[0], dt_string))



###=====================================================================================================================



for bam_file_counter in range(len(bam_file_list)): 

    print('Processing bam file: {}'.format(bam_file_name_list[bam_file_counter]))

    if args.make_bdg:

        print('Generating normalized bedgraph file with normalization scaling factor of {}'.format(norm_scale_factor_list[bam_file_counter]))

        if args.log_dir:
            bedgraph_log_string = '1> {}/{}.bedGraphing.out 2> {}/{}.bedGraphing.err'.format(log_dir, bam_file_name_list[bam_file_counter], log_dir, bam_file_name_list[bam_file_counter])
        else:
            bedgraph_log_string = ''

        # Running bamCoverage to generate a bedgraph file, to be used later as input for SEACR
        subprocess.run('bamCoverage -p {} -b {} -o {}/{}.normalized.bdg --skipNAs -of bedgraph --binSize {} --scaleFactor {} {}'.format(
            cpu_count, 
            bam_file_list[bam_file_counter],
            bam_file_dir_list[bam_file_counter],
            bam_file_name_list[bam_file_counter],
            bdg_bin_size,
            norm_scale_factor_list[bam_file_counter],
            bedgraph_log_string),
            shell = True)

        print('Normalized bedgraph file succesfully generated: {}/{}.normalized.bdg'.format(bam_file_dir_list[bam_file_counter], bam_file_name_list[bam_file_counter]))


    if args.make_bw:

        print('Generating normalized bigwig file with normalization scaling factor of {}'.format(norm_scale_factor_list[bam_file_counter]))

        if args.log_dir:
            bigwig_log_string = '1> {}/{}.bigWiging.out 2> {}/{}.bigWiging.err'.format(log_dir, bam_file_name_list[bam_file_counter], log_dir, bam_file_name_list[bam_file_counter])
        else:
            bigwig_log_string = ''

        # Running bamCoverage to generate a bigwig file, to be used later for peak visualization
        subprocess.run('bamCoverage -p {} -b {} -o {}/{}.normalized.bw -of bigwig --binSize {} --scaleFactor {} {}'.format(
            cpu_count, 
            bam_file_list[bam_file_counter],
            bam_file_dir_list[bam_file_counter],
            bam_file_name_list[bam_file_counter],
            bdg_bin_size,
            norm_scale_factor_list[bam_file_counter],
            bigwig_log_string),
            shell = True)

        print('Normalized bigwig file succesfully generated: {}/{}.normalized.bw'.format(bam_file_dir_list[bam_file_counter], bam_file_name_list[bam_file_counter]))



###=====================================================================================================================



    if args.make_bam:

        if norm_mode != 'none' and len(bam_file_list) > 1 and args.norm_bam and round(norm_scale_factor_list[bam_file_counter]) > 1:
            
            print('Generating normalized bam file with rounded normalization scaling factor of {}'.format(round(norm_scale_factor_list[bam_file_counter])))

            current_bam_filename = bam_file_list[bam_file_counter].rstrip('.bam')

            merge_bam_list = ['{}.bam'.format(current_bam_filename) for norm_scale_factor_counter in range(round(norm_scale_factor_list[bam_file_counter]))]
            merge_bam_string = (' '.join(merge_bam_list))

            subprocess.run('samtools merge -@ {} -f -c -p {}.normalized.bam {}'.format(cpu_count, current_bam_filename, merge_bam_string), shell = True)

            print('Normalization done. Normalized bam file is saved as {}.normalized.bam'.format(current_bam_filename))

        else:

            subprocess.run('cp -r -f {} {}/{}.normalized.bam'.format(
                bam_file_list[bam_file_counter],
                bam_file_dir_list[bam_file_counter],
                bam_file_name_list[bam_file_counter]),
                shell = True)

            print('No normalization was done. The original bam file is copy-pasted to {}/{}.normalized.bam'.format(bam_file_dir_list[bam_file_counter], bam_file_name_list[bam_file_counter]))