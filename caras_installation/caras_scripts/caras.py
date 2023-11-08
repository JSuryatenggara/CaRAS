#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


# script_version = '5.0'


# A fully automated ChIP-seq pipeline, processes raw sequencing reads (or aligned reads),
#   identifies potential protein binding sites, augment them with related informations,
#   and finally compile everything in one comprehensive, spreadsheet-compatible, complete list.

# INPUT         - Chromatin IP (+ replicates if available):
#                   - FASTQ or gzip -fped FASTQ files if reads are raw sequencer output
#                   - BAM files if reads are already preprocessed and aligned
#               - Background control (+ replicates if available): 
#                   - FASTQ or gzip -fped FASTQ files if reads are raw sequencer output
#                   - BAM files if reads are already preprocessed and aligned

# GENERATED     - 00_raw_data_script.sh                             - copy and rename raw sequencing data
# SCRIPTS       - 01_raw_reads_quality_control_script.sh            - quality check on raw sequencing data
#               - 02_deduplicating_script.sh                        - removes duplicate reads
#               - 02_adapter_trimming_script.sh                     - trims away adapter sequences
#               - 03_quality_trimming_script.sh                     - trims away low quality base reads
#               - 05_preprocessed_reads_quality_control_script.sh   - quality check on preprocessed sequencing data
#               - 06_bwa_mem_aligning_script.sh                     - map and aligns the reads to the reference genome
#               - 07_MAPQ_filtering_script.sh                       - filter out reads with low alignment score
#               - 08_results_script.sh                              - generates the files of the aligned reads
#                                                                       - BAM (+ index) files for further analysis
#                                                                       - BigWig files for track visualization
#               - 09_aligned_reads_quality_control_script.sh        - quality check on aligned sequencing data
#               - 11_macs2_peak_calling_script.sh                   - identifies potential binding sites using MACS2
#               - 12_gem_peak_calling_script.sh                     - identifies potential binding sites using GEM
#               - 13_homer_peak_calling_script.sh                   - identifies potential binding sites using HOMER
#               - 14_genrich_peak_calling_script.sh                 - identifies potential binding sites using Genrich
#               - 21_peaks_merging_script.sh                        - merge together all combinations of called peaks
#               - 22_peaks_processing_script.sh                     - adds genetic annotations, calculate peaks stats
#               - 23_go_annotation_script.sh                        - augment related genes with GO terms from databases
#               - 23_pathway_annotation_script.sh                   - augment related genes with pathways from databases
#       
#               - MASTER_script.sh                                  - calls and runs all scripts above in ordered fashion

# OUTPUT        - dataset_name_all_peaks_calculated                 - standard output (genetic annotations)
#               - dataset_name_all_peaks_GO_annotated               - standard output + GO terms (optional)
#               - dataset_name_all_peaks_pathway_annotated          - standard output + pathways (optional)
#               - dataset_name_all_peaks_GO_pathway_annotated       - standard output + GO terms + pathways (optional)
#               - dataset_name_peak_caller_combinations_statistics  - statistical summary of all called peak sets


# PATCH NOTES
# 
# 
# 


# Import required modules
print('Importing required modules')

from os.path import dirname as up
import argparse
import os
import errno
import math
import multiprocessing
import subprocess
import numpy as np
from shutil import which
import pandas as pd
import requests
import pathlib


print('Defining functions')

# Function to parse the absolute path, basename, and extension of the assigned argument: 

# INPUT     - file_relative_path_list - a list of relative paths to n files

# OUTPUT    - absolute_path - a list of the absolute paths of file_relative_path_list elements
#           - original_name - a list of the file basenames of file_relative_path_list elements
#           - original_extension a list of the file extensions of file_relative_path_list elements

def file_basename_parsing(file_relative_path_list):
    absolute_path = [os.path.abspath(relative_path_element) for relative_path_element in file_relative_path_list]
        
    original_name = []
    original_extension = []

    for absolute_path_element in absolute_path:
        if absolute_path_element.lower().endswith('.fq.gz'):
            current_extension = '.fq.gz'
            current_file_name = absolute_path_element.split('/')[-1].replace(current_extension, '')
            original_name.append(current_file_name)
            original_extension.append(current_extension)
            
        if absolute_path_element.lower().endswith('.fastq.gz'):
            current_extension = '.fastq.gz'
            current_file_name = absolute_path_element.split('/')[-1].replace(current_extension, '')
            original_name.append(current_file_name)
            original_extension.append(current_extension)

        if absolute_path_element.lower().endswith('.fq'):
            current_extension = '.fq'
            current_file_name = absolute_path_element.split('/')[-1].replace(current_extension, '')
            original_name.append(current_file_name)
            original_extension.append(current_extension)

        if absolute_path_element.lower().endswith('.fastq'):
            current_extension = '.fastq'
            current_file_name = absolute_path_element.split('/')[-1].replace(current_extension, '')
            original_name.append(current_file_name)
            original_extension.append(current_extension)

        if absolute_path_element.lower().endswith('.bam'):
            current_extension = '.bam'
            current_file_name = absolute_path_element.split('/')[-1].replace(current_extension, '')
            original_name.append(current_file_name)
            original_extension.append(current_extension)
            
    return absolute_path, original_name, original_extension



# Function to search for the location of programs that are required to be called 
#   by their full path (e.g. .jar files of trimmomatic & GEM)

# INPUT     - program_name - the name of the program (as how it is called in bash terminal) to be searched for
#           - search_path - the directory under which the program will be searched in

# OUTPUT    - True/False - True if the program exists. False if the program does not exist.
#           - PATH/None - the path to the program searched, if it exists. None if the program does not exist.

def find_program(program_name, search_path):
    for program_root, program_dirs, program_files in os.walk(search_path):
        if program_name in program_files:
            return True, os.path.join(program_root, program_name)

# program_dir is not used nor assigned in any of downstream processes and operations, 
#   yet necessary due to the structure of os.walk; thus, may cause mild alerts.

    for program_root, program_dirs, program_files in os.walk('/'):
        if program_name in program_files:
            return True, os.path.join(program_root, program_name)

    return False, None



# Function to actually try to call the programs by python subprocess module, 
#   to test if they can at least be started up, during suite requirement checks

# INPUT     - program_name - the name of the program (as how it is called in bash terminal) to be tested for a run

# OUTPUT    - True/False - True if the program test run does not exit in error. False if the run exits in error.

def check_program(program_name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([program_name], stdout = devnull, stderr = devnull, shell = True).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True



# Set arguments to be parsed when the script is called
print('Setting argument parser')

parser = argparse.ArgumentParser()

parser.add_argument('--analysis', 
                    help = '<Required> Bulk sample or single cell analysis', 
                    required = True, 
                    choices = ['bulk', 'single_cell'])

parser.add_argument('--mode', 
                    help = '<Required> Single-end or paired-end sequencing read.', 
                    required = True, 
                    choices = ['single', 'paired'])

parser.add_argument('--genome', 
                    help = '<Required> Your genome folder.', 
                    required = True)

parser.add_argument('--output', 
                    help = '<Required> Your desired output folder.', 
                    required = True)

parser.add_argument('--setname', 
                    help = '<Required> The prefix that will be used to label output and intermediate files.', 
                    required = True)

parser.add_argument('--chipR1', 
                    nargs = '+', 
                    help = '<Optional> Your ChIP datasets: ordered by replicate number, separated by space.')

parser.add_argument('--chipR2', 
                    nargs = '+', 
                    help = '<Optional. Paired-end mode only> Your ChIP second read datasets: ordered by replicate number, separated by space.')

parser.add_argument('--ctrlR1', 
                    nargs = '+', 
                    help = '<Optional> Your control datasets: ordered by replicate number, separated by space.')

parser.add_argument('--ctrlR2', 
                    nargs = '+', 
                    help = '<Optional. Paired-end mode only> Your control second read datasets: ordered by replicate number, separated by space.')

parser.add_argument('--sample_table', 
                    help = '<Optional> Your 4*n-sized (with n = number of replicates) table file containing the absolute paths to each of your ChIP and control sample replicates. When used, this table will disregard any assigned values to --chipR1, --chipR2, --ctrlR1, and --ctrlR2')

parser.add_argument('--custom_setting_table', 
                    help = '<Optional> Expert mode. Your 2*n-sized table file containing custom arguments for every modules.')

parser.add_argument('--fcmerge', 
                    help = '<Optional> Use to force fold change analysis based on merged replicates instead of on each replicate. Only relevant for bulk analysis', 
                    action = 'store_true')

parser.add_argument('--clustnum', 
                    help = '<Optional> Expected number of clusters to be formed by all single-cell samples. When not provided, CaRAS will attempt to calculate the optimal number of clusters by itself (auto). Only relevant for single-cell analysis. Default value is 0 (auto).', 
                    type = int,
                    default = 0)

parser.add_argument('--motif', 
                    help = '<Optional> Your predicted/known motif file, in HOMER matrix format, .motif extension')

parser.add_argument('--ref', 
                    help = '<Optional> Your sample organism genome reference build. Default is hg38 (human).', 
                    choices = ['hg19', 'hg38', 'mm9', 'mm10', 'mm39', 'dm6', 'sacCer3'], 
                    default = 'hg38')

parser.add_argument('--norm_ref', 
                    help = '<Optional> The spiked-in or carry-over DNA genome used in the Cut and Run protocol. For peak calling signal normalization. Default is no normalization (none).', 
                    choices = ['none', 'sacCer3', 'eColiK12'])

parser.add_argument('--goann', 
                    help = '<Optional> Use to annotate peaks with all relevant GO terms.', 
                    action = 'store_true')

parser.add_argument('--pathann', 
                    help = '<Optional> Use to annotate peaks with all common pathways, interactions, and occurences.', 
                    action = 'store_true')

parser.add_argument('--homer_motif', 
                    help = '<Optional> Use to perform motif enrichment analysis on selected peak set with HOMER motif analysis suite.', 
                    nargs = '+', 
                    choices = ['consensus', 'union', '1', '2', '3', '4', '5', '6'])

parser.add_argument('--meme_motif', 
                    help = '<Optional> Use to perform motif enrichment analysis on selected peak set with MEME motif analysis suite.', 
                    nargs = '+', 
                    choices = ['consensus', 'union', '1', '2', '3', '4', '5', '6'])

parser.add_argument('--thread', 
                    help = '<Optional> Maximum number of processes to use. Default is half the maximum available.', 
                    type = int, 
                    choices = range(1, (multiprocessing.cpu_count() + 1), 1),
                    metavar = "[1-{}]".format(multiprocessing.cpu_count()),
                    default = int(multiprocessing.cpu_count() / 2))

parser.add_argument('--deltemp', 
                    help = '<Optional> Use to immediately delete large intermediary sequencing read files right after the point where they are not going to be used for any further process. Example: fastq files, sam files', 
                    action = 'store_true')

parser.add_argument('--run', 
                    help = '<Optional> Use to immediately run the suite by running the master script. When not used, the generated bash master script (MASTER_script.sh) in the output folder can be run manually by user', 
                    action = 'store_true')

args = parser.parse_args()



########################################################################################################################
### PARSED FLAG ARGUMENTS VARIABLE ASSIGNMENT
########################################################################################################################

# Assign parsed arguments into variables
print('Parsing arguments')
caras_program_name = 'caras.py'

dataset_name = args.setname # Resulting files will be named based in this
read_mode = args.mode # Sequencing read mode of the fastq files. Single-end or paired-end
peak_type = 'narrow'

genome_dir = os.path.abspath(args.genome) # Absolute path to the genome folder
output_folder = os.path.abspath(args.output) # Absolute path to the folder [dataset_name] that contains all the outputs of the suite
output_dir = '{}/{}'.format(os.path.abspath(args.output), dataset_name) # Absolute path to all the outputs of the suite

# Known motif for the target IP protein is assigned here. If no motif was assigned, 
#   all motif overlaps number will be all zero later in the analysis
if args.motif:
    motif_file_full_path = os.path.abspath(args.motif) # Absolute path to the .motif file (HOMER motif matrix file)
if not args.motif:
    motif_file_full_path = None

complete_genome = ['hg19', 'hg38', 'mm9', 'mm10', 'dm6', 'sacCer3']

# Reference genome to be used depending on sample organism. For now, the suite can take in human and mouse samples
if args.ref == 'hg19':
    genome_ref = 'hg19' # For human samples
if args.ref == 'hg38':
    genome_ref = 'hg38' # For human samples
if args.ref == 'mm9':
    genome_ref = 'mm9' # For mice samples
if args.ref == 'mm10':
    genome_ref = 'mm10' # For mice samples
if args.ref == 'mm39':
    genome_ref = 'mm39' # For mice samples
if args.ref == 'dm6':
    genome_ref = 'dm6' # For fruitfly samples
if args.ref == 'sacCer3':
    genome_ref = 'sacCer3' # For yeast samples

# Reference genome to be used depending on the spiked-in or carry-over DNA genome used in the Cut and Run protocol.
if args.norm_ref:
    norm_ref = args.norm_ref
else:
    norm_ref = 'none'

force_merge = 1 if args.fcmerge else 0

analysis_mode = args.analysis

expected_cluster = args.clustnum

if args.homer_motif:
    homer_motif_peakset = [str(peakset) for peakset in args.homer_motif]

if args.meme_motif:
    meme_motif_peakset = [str(peakset) for peakset in args.meme_motif]

print('CaRAS set on {} mode'.format(analysis_mode))

home_dir = os.path.expanduser('~') # Absolute path to the user's home directory
root_dir = os.path.expanduser('/') # Absolute path to the root directory

# Assigning a general number of threads to be used by all called function in the suite
# If not user-determined, the suite will use all available processing cores available in the system where it is run
cpu_count = args.thread

# Getting the directory on which the script was called
current_dir = str(pathlib.Path(__file__).parent.absolute())



# ########################################################################################################################
# ### SOFTWARE UPDATE CHECK
# ########################################################################################################################

print('\nChecking for updates on GitHub\n')

remote_directory = 'https://raw.githubusercontent.com/JSuryatenggara/CaRAS/main/caras_installation'

caras_scripts_dir = str(pathlib.Path(__file__).parent.absolute())
caras_installation_dir = '/'.join((caras_scripts_dir.split('/'))[:-1])
local_directory = caras_installation_dir

program_update_check_list = ['caras_scripts/caras.py',
                            'caras_scripts/caras_dashboard.py',
                            'caras_scripts/caras_wizard.py',
                            'caras_scripts/caras_reads_normalizer.py',
                            'caras_scripts/caras_Genrich.py',
                            'caras_scripts/caras_upset_plotter.py',
                            'caras_scripts/caras_peak_feature_extractor.py',
                            'caras_scripts/caras_fold_change_calculator.py',
                            'caras_scripts/caras_IDR_integrator.py',
                            'caras_scripts/caras_GO_annotator.py',
                            'caras_scripts/caras_pathway_annotator.py',
                            'caras_scripts/caras_peak_caller_stats_calculator.py',
                            'caras_scripts/caras_meme_sequence_extractor.py',
                            'caras_scripts/caras_default_settings_table.tsv',
                            'caras_scripts/SEACR_1.3.sh',
                            'caras_scripts/SEACR_1.3.R',
                            'caras_env_linux.yml',
                            'caras_env_macos.yml',
                            'caras_installer.py',
                            'homer_genome_update.sh',
                            'idr.py',
                            'chromfa_splitter.py']


update_counter = 0

for program in program_update_check_list:
    
    try:
        remote_file = requests.get('{}/{}'.format(remote_directory, program))

    except:
        print('{} is no longer required in the latest implementation of CaRAS'.format(program.split('/')[-1]))
        continue
        
    try:
        local_file = open('{}/{}'.format(local_directory, program), 'r')

    except:
        print('WARNING: {} is not found in your local system'.format(program.split('/')[-1]))
        continue


    if remote_file.text == local_file.read():
        print('{} is up to date'.format(program.split('/')[-1]))

    elif remote_file.text != local_file.read():
        update_counter += 1
        print('Newer version of {} is available on our github'.format(program.split('/')[-1]))


if update_counter == 0:
    print('\nYour CaRAS is up to date\n')

if update_counter > 0:
    print('\n{} updates available on our GitHub (https://github.com/JSuryatenggara/CaRAS)'.format(update_counter))
    print('Update all at once by downloading our latest release (https://github.com/JSuryatenggara/CaRAS/releases)\n')



########################################################################################################################
### EFFECTIVE GENOME SIZE ASSIGNMENT
########################################################################################################################

# Defining effective genome sizes for all CaRAS supported genomes:
# Source: https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
# sacCer3_effective_genome_size was manually calculated based on the number of non-N characters in the genome's FASTA file

hg19_effective_genome_size      = '2861327131'
hg38_effective_genome_size      = '2937639113'
mm9_effective_genome_size       = '2558509480'
mm10_effective_genome_size      = '2647521431'
mm39_effective_genome_size      = '2649921816'
dm6_effective_genome_size       = '137057575'
sacCer3_effective_genome_size   = '12071326'

if genome_ref == 'hg19':
    effective_genome_size = hg19_effective_genome_size
if genome_ref == 'hg38':
    effective_genome_size = hg38_effective_genome_size
if genome_ref == 'mm9':
    effective_genome_size = mm9_effective_genome_size
if genome_ref == 'mm10':
    effective_genome_size = mm10_effective_genome_size
if genome_ref == 'mm39':
    effective_genome_size = mm39_effective_genome_size
if genome_ref == 'dm6':
    effective_genome_size = dm6_effective_genome_size
if genome_ref == 'sacCer3':
    effective_genome_size = sacCer3_effective_genome_size



########################################################################################################################
### SAMPLE FILES LISTING
########################################################################################################################

# List of all samples absolute path. r2 only applies for paired end reads, empty in single-end reads.
chip_r1_sample_list = []
ctrl_r1_sample_list = []
chip_r2_sample_list = []
ctrl_r2_sample_list = []


# Getting the absolute paths to the samples from direct input to the terminal (when not using --sample_table option)
if not args.sample_table:
    chip_r1_sample_list = [os.path.abspath(chip_r1_sample) for chip_r1_sample in args.chipR1] if args.chipR1 != None else []
    chip_r2_sample_list = [os.path.abspath(ctrl_r1_sample) for ctrl_r1_sample in args.chipR2] if args.chipR2 != None else []
    ctrl_r1_sample_list = [os.path.abspath(chip_r2_sample) for chip_r2_sample in args.ctrlR1] if args.ctrlR1 != None else []
    ctrl_r2_sample_list = [os.path.abspath(ctrl_r2_sample) for ctrl_r2_sample in args.ctrlR2] if args.ctrlR2 != None else []


# Getting the absolute paths to the samples from the sample tab separated table file (only if using --sample_table option)
# Will override samples directly inputted to the terminal.
if args.sample_table:
    samples_table_absolute_path = os.path.abspath(args.sample_table)
    samples_table_df = pd.read_csv(samples_table_absolute_path, delimiter='\t')

    # Fill in the list according to the provided sample inputs. Remove nan from entry in case the lists length are not all equal.
    chip_r1_sample_list = [os.path.abspath(chip_r1_sample) for chip_r1_sample in samples_table_df['chip_read_1'] if str(chip_r1_sample) != 'nan']
    ctrl_r1_sample_list = [os.path.abspath(ctrl_r1_sample) for ctrl_r1_sample in samples_table_df['ctrl_read_1'] if str(ctrl_r1_sample) != 'nan']
    chip_r2_sample_list = [os.path.abspath(chip_r2_sample) for chip_r2_sample in samples_table_df['chip_read_2'] if str(chip_r2_sample) != 'nan']
    ctrl_r2_sample_list = [os.path.abspath(ctrl_r2_sample) for ctrl_r2_sample in samples_table_df['ctrl_read_2'] if str(ctrl_r2_sample) != 'nan']

    if chip_r1_sample_list == None:
        chip_r1_sample_list = []
    if chip_r2_sample_list == None:
        chip_r2_sample_list = []
    if ctrl_r1_sample_list == None:
        ctrl_r1_sample_list = []
    if ctrl_r2_sample_list == None:
        ctrl_r2_sample_list = []

sample_table_output_dict = {'chip_read_1' : chip_r1_sample_list,
                            'chip_read_2' : chip_r2_sample_list,
                            'ctrl_read_1' : ctrl_r1_sample_list,
                            'ctrl_read_2' : ctrl_r2_sample_list}


# Check the length of resulting sample lists. If any is empty, notify the user and then exit the suite.
if read_mode == 'single':

    if not chip_r1_sample_list or len(chip_r1_sample_list) == 0:
        print('You have not assigned the ChIP sample')
        exit()

    if len(chip_r2_sample_list) > 0:
        print('Single-end read mode does not accept ChIP second reads (R2) file')
        exit()

# Check the length of resulting sample lists. If any is empty, notify the user and then exit the suite.
if read_mode == 'paired':

    if not chip_r1_sample_list or len(chip_r1_sample_list) == 0:
        print('You have not assigned the ChIP sample (read 1)')
        exit()

    # Check for r2 only if the input files are not aligned
    if not chip_r2_sample_list or len(chip_r2_sample_list) == 0:
        if not all('.bam' in chip_r1_sample for chip_r1_sample in chip_r1_sample_list):
            print('You have not assigned the ChIP sample (read 2)')
            exit()

if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    if len(chip_r1_sample_list) != len(ctrl_r1_sample_list):
        print('If controls are to be used, every ChIP sample need to be paired with a control sample')
        exit()



########################################################################################################################
### NEW SAMPLE FILENAMES ASSIGNING
########################################################################################################################

if read_mode == 'single':
    # Parsing the absolute path, basename, and extension of the ChIP original files (single-end reads):
    chip_r1_absolute_path, chip_r1_original_name, chip_r1_original_extension = file_basename_parsing(chip_r1_sample_list)

    # Renaming rule: --setname + '_chip_rep_' + replicate number
    chip_r1_name = ['{}_chip_rep{}_R1'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]
    chip_name = ['{}_chip_rep{}'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        # Parsing the absolute path, basename, and extension of the control original files (single-end reads):
        ctrl_r1_absolute_path, ctrl_r1_original_name, ctrl_r1_original_extension = file_basename_parsing(ctrl_r1_sample_list)

        # Renaming rule: --setname + '_ctrl_rep_' + replicate number
        ctrl_r1_name = ['{}_ctrl_rep{}_R1'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(ctrl_r1_original_name))]
        ctrl_name = ['{}_ctrl_rep{}'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(ctrl_r1_original_name))]



elif read_mode == 'paired':
    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        if not all('.bam' in chip_r1_sample for chip_r1_sample in chip_r1_sample_list):
            if not all('.bam' in ctrl_r1_sample for ctrl_r1_sample in ctrl_r1_sample_list):
                # Parsing the absolute path, basename, and extension of the ChIP original files (paired-end reads r1 and r2):
                chip_r1_absolute_path, chip_r1_original_name, chip_r1_original_extension = file_basename_parsing(chip_r1_sample_list)
                chip_r2_absolute_path, chip_r2_original_name, chip_r2_original_extension = file_basename_parsing(chip_r2_sample_list)

                # Renaming rule: --setname + '_chip_rep_' + replicate number for unpaired file (without r2)
                # Renaming rule: --setname + '_chip_rep_' + replicate number + _R1 or _R2 for paired files (with r2)
                chip_r1_name = ['{}_chip_rep{}_R1'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]
                chip_r2_name = ['{}_chip_rep{}_R2'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r2_original_name))]
                chip_name = ['{}_chip_rep{}'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]

                # Parsing the absolute path, basename, and extension of the control original files (paired-end reads r1 and r2):
                ctrl_r1_absolute_path, ctrl_r1_original_name, ctrl_r1_original_extension = file_basename_parsing(ctrl_r1_sample_list)
                ctrl_r2_absolute_path, ctrl_r2_original_name, ctrl_r2_original_extension = file_basename_parsing(ctrl_r2_sample_list)

                # Renaming rule: --setname + '_ctrl_rep_' + replicate number for unpaired file (without r2)
                # Renaming rule: --setname + '_ctrl_rep_' + replicate number + _R1 or _R2 for paired files (with r2)
                ctrl_r1_name = ['{}_ctrl_rep{}_R1'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(ctrl_r1_original_name))]
                ctrl_r2_name = ['{}_ctrl_rep{}_R2'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(ctrl_r2_original_name))]
                ctrl_name = ['{}_ctrl_rep{}'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(ctrl_r1_original_name))]

    else:
        if not all('.bam' in chip_r1_sample for chip_r1_sample in chip_r1_sample_list):
            # Parsing the absolute path, basename, and extension of the ChIP original files (paired-end reads r1 and r2):
            chip_r1_absolute_path, chip_r1_original_name, chip_r1_original_extension = file_basename_parsing(chip_r1_sample_list)
            chip_r2_absolute_path, chip_r2_original_name, chip_r2_original_extension = file_basename_parsing(chip_r2_sample_list)

            # Renaming rule: --setname + '_chip_rep_' + replicate number for unpaired file (without r2)
            # Renaming rule: --setname + '_chip_rep_' + replicate number + _R1 or _R2 for paired files (with r2)
            chip_r1_name = ['{}_chip_rep{}_R1'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]
            chip_r2_name = ['{}_chip_rep{}_R2'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r2_original_name))]
            chip_name = ['{}_chip_rep{}'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]

    # Check for r2 only if the input files are not aligned
    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        if all('.bam' in chip_r1_sample for chip_r1_sample in chip_r1_sample_list):
            if all('.bam' in ctrl_r1_sample for ctrl_r1_sample in ctrl_r1_sample_list):
                # Parsing the absolute path, basename, and extension of the ChIP original files (single-end reads):
                chip_r1_absolute_path, chip_r1_original_name, chip_r1_original_extension = file_basename_parsing(chip_r1_sample_list)

                # Renaming rule: --setname + '_chip_rep_' + replicate number
                chip_r1_name = ['{}_chip_rep{}_R1'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]
                chip_name = ['{}_chip_rep{}'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]


                # Parsing the absolute path, basename, and extension of the control original files (single-end reads):
                ctrl_r1_absolute_path, ctrl_r1_original_name, ctrl_r1_original_extension = file_basename_parsing(ctrl_r1_sample_list)

                # Renaming rule: --setname + '_ctrl_rep_' + replicate number
                ctrl_r1_name = ['{}_ctrl_rep{}_R1'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(ctrl_r1_original_name))]
                ctrl_name = ['{}_ctrl_rep{}'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(ctrl_r1_original_name))]

    else:
        if all('.bam' in chip_r1_sample for chip_r1_sample in chip_r1_sample_list):
            # Parsing the absolute path, basename, and extension of the ChIP original files (single-end reads):
            chip_r1_absolute_path, chip_r1_original_name, chip_r1_original_extension = file_basename_parsing(chip_r1_sample_list)

            # Renaming rule: --setname + '_chip_rep_' + replicate number
            chip_r1_name = ['{}_chip_rep{}_R1'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]
            chip_name = ['{}_chip_rep{}'.format(dataset_name, (list_counter + 1)) for list_counter in range(len(chip_r1_original_name))]



if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    if '.bam' in chip_r1_original_extension and '.bam' in ctrl_r1_original_extension:
        if len(set(chip_r1_original_extension)) == 1 and len(set(ctrl_r1_original_extension)) == 1:
            start_from_bam = True

            # Update the required samples to waive the read 2 when all sample files are in .bam format
            sample_table_output_dict = {'chip_read_1' : chip_r1_sample_list,
                                        'chip_read_2' : [],
                                        'ctrl_read_1' : ctrl_r1_sample_list,
                                        'ctrl_read_2' : []}

        else:
            print('You cannot combine aligned read files (.bam) with raw read files (e.g. fastq)')
            print('Exiting program')
            quit()

    else:
        start_from_bam = False

else:
    if '.bam' in chip_r1_original_extension:
        if len(set(chip_r1_original_extension)) == 1:
            start_from_bam = True

            # Update the required samples to waive the read 2 when all sample files are in .bam format
            sample_table_output_dict = {'chip_read_1' : chip_r1_sample_list,
                                        'chip_read_2' : [],
                                        'ctrl_read_1' : [],
                                        'ctrl_read_2' : []}

        else:
            print('You cannot combine aligned read files (.bam) with raw read files (e.g. fastq)')
            print('Exiting program')
            quit()

    else:
        start_from_bam = False



########################################################################################################################
### CUSTOM SETTINGS TABLE READING
########################################################################################################################

if args.custom_setting_table:
    # Read the custom settings table, if provided in the command
    custom_settings_table_full_path = os.path.abspath(args.custom_setting_table)

if not args.custom_setting_table:
    # If not provided, suite will read the default custom settings table, currently provided in the genome folder
    if os.path.isfile('{}/caras_default_settings_table.tsv'.format(args.genome)):
        custom_settings_table_full_path = os.path.abspath('{}/caras_default_settings_table.tsv'.format(args.genome))

    else:
        caras_path_check = subprocess.Popen('which caras.py', shell = True, stdout = subprocess.PIPE) # Get the full path to caras.py
        caras_path = caras_path_check.communicate()[0].decode('utf-8')
        caras_dir = '/'.join(caras_path.split('/')[:-1])
        print('{}/caras_default_settings_table.tsv'.format(caras_dir))
        
        if os.path.isfile('{}/caras_default_settings_table.tsv'.format(caras_dir)):
            custom_settings_table_full_path = '{}/caras_default_settings_table.tsv'.format(caras_dir)

        else:
            print('Could not find default_settings_table.tsv in default paths. Please provide the path manually using --custom_setting_table flag')
            exit()

custom_settings_table_df        = pd.read_csv(custom_settings_table_full_path, delimiter='\t')
custom_settings_table_df        = custom_settings_table_df.replace(np.nan, '', regex = True)
custom_settings_table_header    = custom_settings_table_df.columns.values.tolist()

# List of the programs in the suite, which default settings are subject to be customized by the table.
suite_program_list = ["fastqc1",
                        "clumpify",
                        "bbduk",
                        "trimmomatic",
                        "fastqc2",
                        "bwa_mem",
                        "samtools_view",
                        "plotfingerprint",
                        "reads_normalizer",
                        "fastqc3",
                        "macs2_callpeak",
                        "gem",
                        "sicer2",
                        "homer_findPeaks",
                        "genrich",
                        "seacr",
                        "homer_mergePeaks",
                        "homer_annotatePeaks",
                        "fold_change_calculator",
                        "peak_feature_extractor",
                        "homer_findMotifsGenome",
                        "meme_chip"]

# Check the formatting of the custom settings table, to ensure correct program-argument readings.
# Check if the table consists of two columns 
if len(custom_settings_table_header) != 2:
    print('Custom settings table parsing error. Custom settings table should consist of 2 columns')
    print('Exiting program')
    exit()

# Check if the table headers are 'program' and 'argument'. Check first if they are both strings to avoid TypeError.
if isinstance(custom_settings_table_header[0], str) and isinstance(custom_settings_table_header[1], str):
    if custom_settings_table_header[0].strip().lower() != 'program' or custom_settings_table_header[1].strip().lower() != 'argument':
        print("Please make sure the custom settings table header of the two columns are 'program' and 'argument' (case sensitive)")
        # Exit program if not standard table format
        exit()

else:
    print("Please make sure the custom settings table header of the two columns are 'program' and 'argument' (case sensitive)")
    # Exit program if not standard table format
    exit()

# Parse the location of 'program' and 'argument' columns
custom_settings_table_program_colnum    = custom_settings_table_header.index('program')
custom_settings_table_argument_colnum   = custom_settings_table_header.index('argument')
custom_settings_table_array             = custom_settings_table_df.values.tolist()

argument_dict = {}

for suite_program in suite_program_list:
    argument_dict[suite_program] = []

# For each entry in the 'program' column that matched with any of the program name in suite_program_list, 
#   assign the argument as the value and program as the key in dictionary argument_dict
for custom_settings_table_array_row in custom_settings_table_array:
    if custom_settings_table_array_row[custom_settings_table_program_colnum] in suite_program_list:
        current_custom_settings_table_program = custom_settings_table_array_row[custom_settings_table_program_colnum]
        current_custom_settings_table_argument = custom_settings_table_array_row[custom_settings_table_argument_colnum]
        argument_dict[current_custom_settings_table_program].append(current_custom_settings_table_argument)
        if current_custom_settings_table_argument != '' and peak_type != 'unsure':
            # Declare all inserted arguments into the vanilla program call. Inserted = defaults, customs, tweaked defaults.
            print('Adding "{}" to {} call command as optional argument(s)'.format(
                current_custom_settings_table_argument, 
                current_custom_settings_table_program))

# Finally, join all arguments value within each program key with a single space and 
#   assign the joined string into their own variable for easier calling later downstream
fastqc1_arg                 = ' ' + ' '.join(argument_dict['fastqc1']).strip(' ')
clumpify_arg                = ' ' + ' '.join(argument_dict['clumpify']).strip(' ')
bbduk_arg                   = ' ' + ' '.join(argument_dict['bbduk']).strip(' ')
trimmomatic_arg             = ' ' + ' '.join(argument_dict['trimmomatic']).strip(' ')
fastqc2_arg                 = ' ' + ' '.join(argument_dict['fastqc2']).strip(' ')
bwa_mem_arg                 = ' ' + ' '.join(argument_dict['bwa_mem']).strip(' ')
samtools_view_arg           = ' ' + ' '.join(argument_dict['samtools_view']).strip(' ')
plotfingerprint_arg         = ' ' + ' '.join(argument_dict['plotfingerprint']).strip(' ')
reads_normalizer_arg        = ' ' + ' '.join(argument_dict['reads_normalizer']).strip(' ')
fastqc3_arg                 = ' ' + ' '.join(argument_dict['fastqc3']).strip(' ')
macs2_callpeak_arg          = ' ' + ' '.join(argument_dict['macs2_callpeak']).strip(' ')
gem_arg                     = ' ' + ' '.join(argument_dict['gem']).strip(' ')
sicer2_arg                  = ' ' + ' '.join(argument_dict['sicer2']).strip(' ')
homer_findPeaks_arg         = ' ' + ' '.join(argument_dict['homer_findPeaks']).strip(' ')
genrich_arg                 = ' ' + ' '.join(argument_dict['genrich']).strip(' ')
seacr_arg                   = ' ' + ' '.join(argument_dict['seacr']).strip(' ')
homer_mergePeaks_arg        = ' ' + ' '.join(argument_dict['homer_mergePeaks']).strip(' ')
homer_annotatePeaks_arg     = ' ' + ' '.join(argument_dict['homer_annotatePeaks']).strip(' ')
fold_change_calculator_arg  = ' ' + ' '.join(argument_dict['fold_change_calculator']).strip(' ')
peak_feature_extractor_arg  = ' ' + ' '.join(argument_dict['peak_feature_extractor']).strip(' ')
homer_findMotifsGenome_arg  = ' ' + ' '.join(argument_dict['homer_findMotifsGenome']).strip(' ')
meme_chip_arg               = ' ' + ' '.join(argument_dict['meme_chip']).strip(' ')

suite_program_arg = [
    fastqc1_arg,
    clumpify_arg,
    bbduk_arg,
    trimmomatic_arg,
    fastqc2_arg,
    bwa_mem_arg,
    samtools_view_arg,
    plotfingerprint_arg,
    reads_normalizer_arg,
    fastqc3_arg,
    macs2_callpeak_arg,
    gem_arg,
    sicer2_arg,
    homer_findPeaks_arg,
    genrich_arg,
    seacr_arg,
    homer_mergePeaks_arg,
    homer_annotatePeaks_arg,
    fold_change_calculator_arg,
    peak_feature_extractor_arg,
    homer_findMotifsGenome_arg,
    meme_chip_arg]



########################################################################################################################
### COMMAND LINE STRING GENERATION
########################################################################################################################

analysis_mode_arg = ' --analysis {}'.format(analysis_mode)

read_mode_arg = ' --mode {}'.format(read_mode)

output_folder_arg = ' --output {}'.format(args.output)

dataset_name_arg = ' --setname {}'.format(dataset_name)

genome_ref_arg = ' --ref {}'.format(genome_ref)

genome_dir_arg = ' --genome {}'.format(genome_dir)

if args.norm_ref:
    norm_ref_arg = ' --norm_ref {}'.format(norm_ref)
else:
    norm_ref_arg = ''

if args.sample_table:
    sample_table_arg = ' --sample_table {}'.format(samples_table_absolute_path)
else:
    sample_table_arg = ''

if args.custom_setting_table:
    setting_table_arg = ' --custom_setting_table {}'.format(custom_settings_table_full_path)
else:
    setting_table_arg = ''

if args.motif:
    motif_file_arg = ' --motif {}'.format(motif_file_full_path)
else:
    motif_file_arg = ''

if args.fcmerge:
    fcmerge_arg = ' --fcmerge'
else:
    fcmerge_arg = ''

if args.goann:
    goann_arg = ' --goann'
else:
    goann_arg = ''

if args.pathann:
    pathann_arg = ' --pathann'
else:
    pathann_arg = ''

if args.homer_motif:
    homer_motif_arg = ' --homer_motif {}'.format(' '.join(homer_motif_peakset))
else:
    homer_motif_arg = ''

if args.meme_motif:
    meme_motif_arg = ' --meme_motif {}'.format(' '.join(meme_motif_peakset))
else:
    meme_motif_arg = ''

if args.deltemp:
    deltemp_arg = ' --deltemp'
else:
    deltemp_arg = ''

if args.thread:
    cpu_count_arg = ' --thread {}'.format(cpu_count)
else:
    cpu_count_arg = ''

if args.run:
    run_arg = ' --run'
else:
    run_arg = ''

if chip_r1_sample_list:
    chip_r1_sample_string = ' '.join(chip_r1_sample_list)
    chip_r1_arg = ' --chipR1 {}'.format(chip_r1_sample_string)
else:
    chip_r1_arg = ''

if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    ctrl_r1_sample_string = ' '.join(ctrl_r1_sample_list)
    ctrl_r1_arg = ' --ctrlR1 {}'.format(ctrl_r1_sample_string)
else:
    ctrl_r1_arg = ''

if chip_r2_sample_list:
    chip_r2_sample_string = ' '.join(chip_r2_sample_list)
    chip_r2_arg = ' --chipR2 {}'.format(chip_r2_sample_string)
else:
    chip_r2_arg = ''

if ctrl_r2_sample_list or len(ctrl_r2_sample_list) != 0:
    ctrl_r2_sample_string = ' '.join(ctrl_r2_sample_list)
    ctrl_r2_arg = ' --ctrlR2 {}'.format(ctrl_r2_sample_string)
else:
    ctrl_r2_arg = ''



# Store the (nigh-identical to the one used) CaRAS command line that will produce identical pipeline processes
if not args.sample_table: # The command line when samples are assigned manually in the command line
    command_line_string = '{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}'.format(caras_program_name,
                                                                                analysis_mode_arg,
                                                                                read_mode_arg,
                                                                                chip_r1_arg,
                                                                                ctrl_r1_arg,
                                                                                chip_r2_arg,
                                                                                ctrl_r2_arg,
                                                                                output_folder_arg,
                                                                                dataset_name_arg,
                                                                                genome_ref_arg,
                                                                                genome_dir_arg,
                                                                                norm_ref_arg,
                                                                                setting_table_arg,
                                                                                motif_file_arg,
                                                                                fcmerge_arg,
                                                                                goann_arg,
                                                                                pathann_arg,
                                                                                homer_motif_arg,
                                                                                meme_motif_arg,
                                                                                deltemp_arg,
                                                                                cpu_count_arg,
                                                                                run_arg)

if args.sample_table: # The command line when samples are assigned automatically by loading the sample table
    command_line_string = '{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}'.format(caras_program_name,
                                                                        analysis_mode_arg,
                                                                        read_mode_arg,
                                                                        output_folder_arg,
                                                                        dataset_name_arg,
                                                                        genome_ref_arg,
                                                                        genome_dir_arg,
                                                                        norm_ref_arg,
                                                                        sample_table_arg,
                                                                        setting_table_arg,
                                                                        motif_file_arg,
                                                                        fcmerge_arg,
                                                                        goann_arg,
                                                                        pathann_arg,
                                                                        homer_motif_arg,
                                                                        meme_motif_arg,
                                                                        deltemp_arg,
                                                                        cpu_count_arg,
                                                                        run_arg)

print('\nNow processing:')
print(command_line_string + '\n')



########################################################################################################################
### SUITE REQUIREMENTS CHECKPOINT
########################################################################################################################

print('Checking for suite requirements')
# 0 for error_status means the suite is not lacking anything necessary and is good to go
error_status = 0

# Running shutil.which with the program name is to check whether the program exists in the PATH
# Running the check_program function with the program name is to check whether the program 
#   can be called and executed in the shell without returning any error
# Running some function os.path.isfile with some file names are to check the existence of
#   the essential files required for the suite to work from start to end

# If any of the suite requirements is not met, it will set the error_status to 1 and 
#   the script will print some troubleshooting instructions for the user

if which('fastqc') is None:
    print('Please make sure fastqc is installed, in PATH, and marked as executable')
    error_status = 1

if check_program('fastqc --help') is False:
    print('Pre-run test of fastqc failed. Please try running fastqc individually to check for the problem')
    error_status = 1

if which('clumpify.sh') is None:
    print('Please make sure clumpify.sh is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('clumpify.sh') is False:
    print('Pre-run test of clumpify.sh failed. Please try running clumpify.sh individually to check for the problem')
    error_status = 1

if which('bbduk.sh') is None:
    print('Please make sure bbduk.sh is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('bbduk.sh') is False:
    print('Pre-run test of bbduk.sh failed. Please try running bbduk.sh individually to check for the problem')
    error_status = 1

# Try to look for trimmomatic in the home directory first
trimmomatic_exist, trimmomatic_full_path = find_program('trimmomatic.jar', home_dir)
if not trimmomatic_exist:
    # If trimmomatic does not exist in the home directory, try looking in the root directory
    trimmomatic_exist, trimmomatic_full_path = find_program('trimmomatic.jar', root_dir)
if not trimmomatic_exist:
    print('Please make sure trimmomatic is installed in your computer')
    error_status = 1

if check_program('java -jar {}'.format(trimmomatic_full_path)) is False:
    print('Pre-run test of trimmomatic failed. Please try running trimmomatic individually to check for the problem')
    error_status = 1

if which('bwa') is None:
    print('Please make sure bwa is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('bwa') is False:
    print('Pre-run test of bwa failed. Please try running bwa individually to check for the problem')
    error_status = 1

if which('samtools') is None:
    print('Please make sure samtools is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('samtools') is False:
    print('Pre-run test of samtools failed. Please try running samtools individually to check for the problem')
    error_status = 1

if which('bamCoverage') is None:
    print('Please make sure deeptools is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('bamCoverage') is False:
    print('Pre-run test of bamCoverage failed. Please try running bamCoverage individually to check for the problem')
    error_status = 1

if which('macs2') is None:
    print('Please make sure MACS2 is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('macs2') is False:
    print('Pre-run test of macs2 failed. Please try running macs2 individually to check for the problem')
    error_status = 1

# Try to look for gem in the home directory first
gem_exist, gem_full_path = find_program('gem.jar', home_dir)
if not gem_exist:
    # If gem does not exist in the home directory, try looking in the root directory
    gem_exist, gem_full_path = find_program('gem.jar', root_dir)
if not gem_exist:
    print('Please make sure gem is installed in your computer')
    error_status = 1

if check_program('java -jar {}'.format(gem_full_path)) is False:
    print('Pre-run test of gem failed. Please try running gem individually to check for the problem')
    error_status = 1

if which('makeTagDirectory') is None:
    print('Please make sure HOMER is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('makeTagDirectory') is False:
    print('Pre-run test of makeTagDirectory failed. Please try running HOMER makeTagDirectory individually to check for the problem')
    error_status = 1

if which('findPeaks') is None:
    print('Please make sure HOMER is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('findPeaks') is False:
    print('Pre-run test of findPeaks failed. Please try running HOMER findPeaks individually to check for the problem')
    error_status = 1

if which('Genrich') is None:
    print('Please make sure Genrich is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('Genrich') is False:
    print('Pre-run test of Genrich failed. Please try running Genrich individually to check for the problem')
    error_status = 1

if which('caras_Genrich.py') is None:
    print('Please make sure caras_Genrich.py is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('caras_Genrich.py') is False:
    print('Pre-run test of caras_Genrich.py failed. Please try running caras_Genrich.py individually to check for the problem')
    error_status = 1

# Try to look for SEACR in the home directory first
seacr_r_exist, seacr_r_full_path = find_program('SEACR_1.3.R', home_dir)
if not seacr_r_exist:
    # If seacr does not exist in the home directory, try looking in the root directory
    seacr_r_exist, seacr_r_full_path = find_program('SEACR_1.3.R', root_dir)

if not seacr_r_exist:
    print('Please make sure SEACR is installed in your computer')
    error_status = 1

if which('SEACR_1.3.sh') is None:
    print('Please make sure SEACR_1.3.sh is installed, in PATH, and marked as executable')
    error_status = 1 
    
if check_program('SEACR_1.3.sh') is False:
    print('Pre-run test of SEACR_1.3.sh failed. Please try running SEACR_1.3.sh individually to check for the problem')
    error_status = 1
    
if which('mergePeaks') is None:
    print('Please make sure HOMER is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('mergePeaks') is False:
    print('Pre-run test of mergePeaks failed. Please try running HOMER mergePeaks individually to check for the problem')
    error_status = 1

# # Try to look for upSet plot in the home directory first
# upset_plot_exist, upset_plot_full_path = find_program('multi_peaks_UpSet.R', home_dir)
# if not upset_plot_exist:
#     # If upSet plot does not exist in the home directory, try looking in the root directory
#     upset_plot_exist, upset_plot_full_path = find_program('multi_peaks_UpSet.R', root_dir)
# if not upset_plot_exist:
#     print('Please make sure multi_peaks_UpSet.R script exists in your computer')
#     error_status = 1

# if check_program(upset_plot_full_path) is False:
#     print('Pre-run test of UpSet plot failed. Please try running multi_peaks_UpSet.R individually to check for the problem')
#     error_status = 1

### MULTIVENN SKIPPED FOR CARAS (CANNOT TAKE MORE THAN 5 SETS WHILE CARAS HAS 6) ###
# # Try to look for multi Venn in the home directory first
# multi_venn_exist, multi_venn_full_path = find_program('multi_peaks_Venn.R', home_dir)
# if not multi_venn_exist:
#     # If multi Venn does not exist in the home directory, try looking in the root directory
#     multi_venn_exist, multi_venn_full_path = find_program('multi_peaks_Venn.R', root_dir)
# if not multi_venn_exist:
#     print('Please make sure multi Venn is installed in your computer')
#     error_status = 1

# if check_program(multi_venn_full_path) is False:
#     print('Pre-run test of multi Venn failed. Please try running multi Venn individually to check for the problem')
#     error_status = 1

if which('annotatePeaks.pl') is None:
    print('Please make sure HOMER is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('annotatePeaks.pl') is False:
    print('Pre-run test of annotatePeaks.pl failed. Please try running HOMER annotatePeaks.pl individually to check for the problem')
    error_status = 1

if which('caras_reads_normalizer.py') is None:
    print('Please make sure caras_reads_normalizer.py is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('caras_reads_normalizer.py') is False:
    print('Pre-run test of caras_reads_normalizer.py failed. Please try running caras_reads_normalizer.py individually to check for the problem')
    error_status = 1

if which('caras_upset_plotter.py') is None:
    print('Please make sure caras_upset_plotter.py is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('caras_upset_plotter.py') is False:
    print('Pre-run test of caras_upset_plotter.py failed. Please try running caras_upset_plotter.py individually to check for the problem')
    error_status = 1

if which('caras_fold_change_calculator.py') is None:
    print('Please make sure caras_fold_change_calculator.py is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('caras_fold_change_calculator.py') is False:
    print('Pre-run test of caras_fold_change_calculator.py failed. Please try running caras_fold_change_calculator.py individually to check for the problem')
    error_status = 1

if which('caras_peak_feature_extractor.py') is None:
    print('Please make sure caras_peak_feature_extractor.py is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('caras_peak_feature_extractor.py') is False:
    print('Pre-run test of caras_peak_feature_extractor.py failed. Please try running caras_peak_feature_extractor.py individually to check for the problem')
    error_status = 1

if which('caras_peak_caller_stats_calculator.py') is None:
    print('Please make sure caras_peak_caller_stats_calculator.py is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('caras_peak_caller_stats_calculator.py') is False:
    print('Pre-run test of caras_peak_caller_stats_calculator.py failed. Please try running caras_peak_caller_stats_calculator.py individually to check for the problem')
    error_status = 1

if which('caras_GO_annotator.py') is None:
    print('Please make sure caras_GO_annotator.py is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('caras_GO_annotator.py') is False:
    print('Pre-run test of caras_GO_annotator.py failed. Please try running caras_GO_annotator.py individually to check for the problem')
    error_status = 1

if which('caras_pathway_annotator.py') is None:
    print('Please make sure caras_pathway_annotator.py is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('caras_pathway_annotator.py') is False:
    print('Pre-run test of caras_pathway_annotator.py failed. Please try running caras_pathway_annotator.py individually to check for the problem')
    error_status = 1

if which('findMotifsGenome.pl') is None:
    print('Please make sure HOMER is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('findMotifsGenome.pl') is False:
    print('Pre-run test of findMotifsGenome.pl failed. Please try running HOMER findMotifsGenome.pl individually to check for the problem')
    error_status = 1

if which('caras_meme_sequence_extractor.py') is None:
    print('Please make sure caras_meme_sequence_extractor.py is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('caras_meme_sequence_extractor.py') is False:
    print('Pre-run test of caras_meme_sequence_extractor.py failed. Please try running caras_meme_sequence_extractor.py individually to check for the problem')
    error_status = 1

if which('meme-chip') is None:
    print('Please make sure MEME is installed, in PATH, and marked as executable')
    error_status = 1 

if check_program('meme-chip') is False:
    print('Pre-run test of meme-chip failed. Please try running meme-chip individually to check for the problem')
    error_status = 1

if not os.path.isfile('{}/bbmap/adapters.fa'.format(genome_dir)):
    print('Please make sure bbduk adapter sequences reference fasta file "adapter.fa" exists, inside folder named "bbmap", inside your genome directory folder')
    print('Example: [YOUR GENOME DIRECTORY]/bbmap/adapters.fa')
    error_status = 1 

if not os.path.isfile('{}/bwa/{}.fa'.format(genome_dir, genome_ref)):
    print('Please make sure bwa aligner genome reference fasta file "{}.fa" exists, inside folder named "bwa", inside your genome directory folder'.format(genome_ref))
    print('Example: [YOUR GENOME DIRECTORY]/bwa/{}.fa'.format(genome_ref))
    error_status = 1 

if not os.path.isfile('{}/GEM/Read_Distribution_default.txt'.format(genome_dir)):
    print('Please make sure GEM default read distribution reference file "Read_Distribution_default.txt" exists, inside folder named "GEM", inside your genome directory folder')
    print('Example: [YOUR GENOME DIRECTORY]/GEM/Read_Distribution_default.txt')
    error_status = 1 

if not os.path.isfile('{}/GEM/{}.chrom.sizes'.format(genome_dir, genome_ref)):
    print('Please make sure GEM {} chromosome sizes reference file "{}.chrom.sizes" exists, inside folder named "GEM", inside your genome directory folder'.format(genome_ref, genome_ref))
    print('Example: [YOUR GENOME DIRECTORY]/GEM/{}.chrom.sizes'.format(genome_ref))
    error_status = 1 

if not os.path.isdir('{}/GEM/{}_Chr_FASTA'.format(genome_dir, genome_ref)):
    print('Please make sure GEM genome reference FOLDER "{}_Chr_FASTA" exists (containing multiple individual fasta files for every chromosome), inside folder named "GEM", inside your genome directory folder'.format(genome_ref))
    print('Example: [YOUR GENOME DIRECTORY]/GEM/{}_Chr_FASTA/'.format(genome_ref))
    error_status = 1 

# error_status of 1 will cause the script to terminate at this point, before going into the main program 
if error_status == 1:
    print('Exiting suite')
    exit()
 




# Every script has #!/bin/bash on top to let the shell know which interpreter to run on the corresponding script
# Every script has set -euxo pipefail following #!/bin/bash to make the scripts run safer (for now)
# The -e option will cause a bash script to exit immediately when a command fails
# The -o pipefail option sets the exit code of a pipeline to that of the rightmost command to exit with a non-zero status, 
#   or to zero if all commands of the pipeline exit successfully
# The -u option causes the bash shell to treat unset variables as an error and exit immediately
# The -x option causes bash to print each command before executing it. For debugging. Will remove from the finished product.

# Creating output directory folder
print('Creating output directory folder')

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print('Generating bash scripts')

########################################################################################################################
### 00_raw_data_script.sh
########################################################################################################################

# Create a directory "00_raw_data" for the compressed raw sequencing reads files (in .fq.gz format)
# Create a script "00_raw_data.sh" within, that copies into this directory and renames all raw sequencing reads files
# All raw sequencing reads files: chip_r1_sample_list, ctrl_r1_sample_list, chip_r2_sample_list, 
#   ctrl_r2_sample_list (assigned directly as command argument or samples_table)
# In: original_filename.fq.gz in original directory
# Out: filename.fq.gz in suite output directory (input for 02_adapter_trimming_script.sh)

raw_data_dir = '{}/00_raw_data'.format(output_dir)
raw_data_script_name = '{}/00_raw_data_script.sh'.format(raw_data_dir)
if not os.path.exists(raw_data_dir):
    os.makedirs(raw_data_dir)
raw_data_script = open(raw_data_script_name, 'w')

raw_data_script.write('#!/bin/bash\n\n')
raw_data_script.write('set -euxo pipefail\n\n')



if start_from_bam == False:
    # Writing bash commands to copy and rename the raw sequence read files into the dataset working directory
    if read_mode == 'single':
        for list_counter in range(len(chip_r1_name)):
            if chip_r1_original_extension[list_counter] == '.fq.gz' or chip_r1_original_extension[list_counter] == '.fastq.gz':

                raw_data_script.write('cp {} {}/{}.fq.gz &\n'.format(
                    chip_r1_absolute_path[list_counter],
                    raw_data_dir,
                    chip_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_r1_name) - 1):
                    raw_data_script.write('\nwait\n\n')

        for list_counter in range(len(chip_r1_name)):
            if chip_r1_original_extension[list_counter] == '.fq' or chip_r1_original_extension[list_counter] == '.fastq':

                raw_data_script.write('cp {} {}/{}.fq &\n'.format(
                    chip_r1_absolute_path[list_counter],
                    raw_data_dir,
                    chip_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_r1_name) - 1):
                    raw_data_script.write('\nwait\n\n')

        for list_counter in range(len(chip_r1_name)):
            if chip_r1_original_extension[list_counter] == '.fq' or chip_r1_original_extension[list_counter] == '.fastq':

                raw_data_script.write('gzip -f {}/{}.fq &\n'.format(
                    raw_data_dir, 
                    chip_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_r1_name) - 1):
                    raw_data_script.write('\nwait\n\n')


        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_r1_name)):
                if ctrl_r1_original_extension[list_counter] == '.fq.gz' or ctrl_r1_original_extension[list_counter] == '.fastq.gz':

                    raw_data_script.write('cp {} {}/{}.fq.gz &\n'.format(
                        ctrl_r1_absolute_path[list_counter], 
                        raw_data_dir, 
                        ctrl_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_r1_name) - 1):
                        raw_data_script.write('\nwait\n\n')

            for list_counter in range(len(ctrl_r1_name)):
                if ctrl_r1_original_extension[list_counter] == '.fq' or ctrl_r1_original_extension[list_counter] == '.fastq':

                    raw_data_script.write('cp {} {}/{}.fq &\n'.format(
                        ctrl_r1_absolute_path[list_counter],
                        raw_data_dir,
                        ctrl_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_r1_name) - 1):
                        raw_data_script.write('\nwait\n\n')

            for list_counter in range(len(ctrl_r1_name)):
                if ctrl_r1_original_extension[list_counter] == '.fq' or ctrl_r1_original_extension[list_counter] == '.fastq':

                    raw_data_script.write('gzip -f {}/{}.fq &\n'.format(
                        raw_data_dir, 
                        ctrl_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_r1_name) - 1):
                        raw_data_script.write('\nwait\n\n')



    elif read_mode == 'paired':
        for list_counter in range(len(chip_r1_name)):
            if chip_r1_original_extension[list_counter] == '.fq.gz' or chip_r1_original_extension[list_counter] == '.fastq.gz':

                raw_data_script.write('cp {} {}/{}.fq.gz &\n'.format(
                    chip_r1_absolute_path[list_counter],
                    raw_data_dir,
                    chip_r1_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_r1_name) - 1):
                    raw_data_script.write('\nwait\n\n')

        for list_counter in range(len(chip_r1_name)):
            if chip_r1_original_extension[list_counter] == '.fq' or chip_r1_original_extension[list_counter] == '.fastq':

                raw_data_script.write('cp {} {}/{}.fq &\n'.format(
                    chip_r1_absolute_path[list_counter],
                    raw_data_dir,
                    chip_r1_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_r1_name) - 1):
                    raw_data_script.write('\nwait\n\n')

        for list_counter in range(len(chip_r1_name)):
            if chip_r1_original_extension[list_counter] == '.fq' or chip_r1_original_extension[list_counter] == '.fastq':

                raw_data_script.write('gzip -f {}/{}.fq &\n'.format(
                    raw_data_dir, 
                    chip_r1_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_r1_name) - 1):
                    raw_data_script.write('\nwait\n\n')


        for list_counter in range(len(chip_r2_name)):
            if chip_r2_original_extension[list_counter] == '.fq.gz' or chip_r2_original_extension[list_counter] == '.fastq.gz':

                raw_data_script.write('cp {} {}/{}.fq.gz &\n'.format(
                    chip_r2_absolute_path[list_counter],
                    raw_data_dir,
                    chip_r2_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_r2_name) - 1):
                    raw_data_script.write('\nwait\n\n')

        for list_counter in range(len(chip_r2_name)):
            if chip_r2_original_extension[list_counter] == '.fq' or chip_r2_original_extension[list_counter] == '.fastq':

                raw_data_script.write('cp {} {}/{}.fq &\n'.format(
                    chip_r2_absolute_path[list_counter],
                    raw_data_dir,
                    chip_r2_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_r2_name) - 1):
                    raw_data_script.write('\nwait\n\n')

        for list_counter in range(len(chip_r2_name)):
            if chip_r2_original_extension[list_counter] == '.fq' or chip_r2_original_extension[list_counter] == '.fastq':

                raw_data_script.write('gzip -f {}/{}.fq &\n'.format(
                    raw_data_dir, 
                    chip_r2_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_r2_name) - 1):
                    raw_data_script.write('\nwait\n\n')


        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_r1_name)):
                if ctrl_r1_original_extension[list_counter] == '.fq.gz' or ctrl_r1_original_extension[list_counter] == '.fastq.gz':

                    raw_data_script.write('cp {} {}/{}.fq.gz &\n'.format(
                        ctrl_r1_absolute_path[list_counter], 
                        raw_data_dir, 
                        ctrl_r1_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_r1_name) - 1):
                        raw_data_script.write('\nwait\n\n')

            for list_counter in range(len(ctrl_r1_name)):
                if ctrl_r1_original_extension[list_counter] == '.fq' or ctrl_r1_original_extension[list_counter] == '.fastq':

                    raw_data_script.write('cp {} {}/{}.fq &\n'.format(
                        ctrl_r1_absolute_path[list_counter],
                        raw_data_dir,
                        ctrl_r1_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_r1_name) - 1):
                        raw_data_script.write('\nwait\n\n')

            for list_counter in range(len(ctrl_r1_name)):
                if ctrl_r1_original_extension[list_counter] == '.fq' or ctrl_r1_original_extension[list_counter] == '.fastq':

                    raw_data_script.write('gzip -f {}/{}.fq &\n'.format(
                        raw_data_dir, 
                        ctrl_r1_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_r1_name) - 1):
                        raw_data_script.write('\nwait\n\n')


            for list_counter in range(len(ctrl_r2_name)):
                if ctrl_r2_original_extension[list_counter] == '.fq.gz' or ctrl_r2_original_extension[list_counter] == '.fastq.gz':

                    raw_data_script.write('cp {} {}/{}.fq.gz &\n'.format(
                        ctrl_r2_absolute_path[list_counter], 
                        raw_data_dir, 
                        ctrl_r2_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_r2_name) - 1):
                        raw_data_script.write('\nwait\n\n')

            for list_counter in range(len(ctrl_r2_name)):
                if ctrl_r2_original_extension[list_counter] == '.fq' or ctrl_r2_original_extension[list_counter] == '.fastq':

                    raw_data_script.write('cp {} {}/{}.fq &\n'.format(
                        ctrl_r2_absolute_path[list_counter],
                        raw_data_dir,
                        ctrl_r2_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_r2_name) - 1):
                        raw_data_script.write('\nwait\n\n')

            for list_counter in range(len(ctrl_r2_name)):
                if ctrl_r2_original_extension[list_counter] == '.fq' or ctrl_r2_original_extension[list_counter] == '.fastq':

                    raw_data_script.write('gzip -f {}/{}.fq &\n'.format(
                        raw_data_dir, 
                        ctrl_r2_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_r2_name) - 1):
                        raw_data_script.write('\nwait\n\n')

raw_data_script.close() # Closing the script '00_raw_data_script.sh'. Flushing the write buffer



########################################################################################################################
### {dataset_name}_run_info.txt
########################################################################################################################

# Writes a text file "{dataset_name}_run_info.txt" that summarizes the assignment of the files 
#   (IP sample or control, read 1 or 2; replicate number)
# Also summarizes the file name conversion for every sequencing reads (.fq.gz) to be processed
# Basically, each line tells the user of what their original files have been renamed into
# This is the list that user will be coming back to when they cannot remember which file they assigned 
#   as what and which (IP sample or control, read 1 or 2; replicate number)

run_info_file_name = '{}/{}_run_info.txt'.format(output_dir, dataset_name)
run_info_file = open(run_info_file_name, 'w')

if start_from_bam == False:

    if read_mode == 'single':
        for list_counter in range(len(chip_r1_name)):
            run_info_file.write('Chromatin IP dataset replicate {} : Original filename = {}{} --> New filename = {}.fq.gz\n\n'.format(
                (list_counter + 1), 
                chip_r1_original_name[list_counter], 
                chip_r1_original_extension[list_counter], 
                chip_name[list_counter]))

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_r1_name)):
                run_info_file.write('Control dataset replicate {} : Original filename = {}{} --> New filename = {}.fq.gz\n\n'.format(
                    (list_counter + 1), 
                    ctrl_r1_original_name[list_counter], 
                    ctrl_r1_original_extension[list_counter], 
                    ctrl_name[list_counter]))

    elif read_mode == 'paired':
        for list_counter in range(len(chip_r1_name)):
            run_info_file.write('Chromatin IP dataset replicate {}, 1st read : Original filename = {}{} --> New filename = {}.fq.gz\n\n'.format(
                (list_counter + 1), 
                chip_r1_original_name[list_counter], 
                chip_r1_original_extension[list_counter], 
                chip_r1_name[list_counter]))

        for list_counter in range(len(chip_r2_name)):
            run_info_file.write('Chromatin IP dataset replicate {}, 2nd read : Original filename = {}{} --> New filename = {}.fq.gz\n\n'.format(
                (list_counter + 1), 
                chip_r2_original_name[list_counter], 
                chip_r2_original_extension[list_counter], 
                chip_r2_name[list_counter]))

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_r1_name)):
                run_info_file.write('Control dataset replicate {}, 1st read : Original filename = {}{} --> New filename = {}.fq.gz\n\n'.format(
                    (list_counter + 1), 
                    ctrl_r1_original_name[list_counter], 
                    ctrl_r1_original_extension[list_counter], 
                    ctrl_r1_name[list_counter]))

            for list_counter in range(len(ctrl_r2_name)):
                run_info_file.write('Control dataset replicate {}, 2nd read : Original filename = {}{} --> New filename = {}.fq.gz\n\n'.format(
                    (list_counter + 1), 
                    ctrl_r2_original_name[list_counter], 
                    ctrl_r2_original_extension[list_counter], 
                    ctrl_r2_name[list_counter]))



if start_from_bam == True:

    for list_counter in range(len(chip_r1_name)):
        run_info_file.write('Chromatin IP dataset replicate {} : Original filename = {}{} --> New filename = {}.bam\n\n'.format(
            (list_counter + 1), 
            chip_r1_original_name[list_counter], 
            chip_r1_original_extension[list_counter], 
            chip_name[list_counter]))

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        for list_counter in range(len(ctrl_r1_name)):
            run_info_file.write('Control dataset replicate {} : Original filename = {}{} --> New filename = {}.bam\n\n'.format(
                (list_counter + 1), 
                ctrl_r1_original_name[list_counter], 
                ctrl_r1_original_extension[list_counter], 
                ctrl_name[list_counter]))

run_info_file.close() # Closing the file 'run_info_file.txt'. Flushing the write buffer



########################################################################################################################
### {dataset_name}_command_line.txt     {dataset_name}_sample_table.tsv     {dataset_name}_setting_table.tsv
########################################################################################################################

# Writes the input terminal single command line that was used to call the pipeline
#   in a text file: "{dataset_name}_command_line.txt in the output directory"

command_line_file_name = '{}/{}_command_line.txt'.format(output_dir, dataset_name)

command_line_file = open(command_line_file_name, 'w')

command_line_file.write((command_line_string) + '\n')

command_line_file.close() # Closing the file '{dataset_name}_command_line.txt'. Flushing the write buffer



# Writes the full path of each inputted ChIP and control samples in the pipeline run 
#   in a tab-separated value file: "{dataset_name}_sample_table.tsv in the output directory"
#   in CaRAS sample table format (see documentation)

sample_table_output_df = pd.DataFrame.from_dict(sample_table_output_dict, orient = 'index')

sample_table_output_df = sample_table_output_df.transpose()

sample_table_output_df.to_csv('{}/{}_sample_table.tsv'.format(output_dir, dataset_name), sep = '\t', index = False)



# Writes the flags and argument values for each program calls used in the pipeline run 
#   in a tab-separated value file: "{dataset_name}_setting_table.tsv in the output directory"
#   in CaRAS setting table format (see documentation)

setting_table_output_dict = {'program' : suite_program_list, 'argument' : suite_program_arg}

setting_table_output_df = pd.DataFrame.from_dict(setting_table_output_dict)

setting_table_output_df['argument'] = setting_table_output_df['argument'].apply(lambda x: x.lstrip(' '))

setting_table_output_df.to_csv('{}/{}_setting_table.tsv'.format(output_dir, dataset_name), sep = '\t', index = False)



####################################################################################################################
### 01_raw_reads_quality_control_script.sh
####################################################################################################################

# Create a directory "01_raw_reads_quality_control" for the quality control report files of the raw sequencing reads
# Create a script "01_raw_reads_quality_control_script.sh" within, 
#   that calls for sequencing reads QC analysis program: fastqc
# In: filename.fq.gz (output from 00_raw_data.sh)
# Out: filename_fastqc.html (Reports. Not to be further processed downstreams)

raw_reads_quality_control_dir = '{}/01_raw_reads_quality_control'.format(output_dir)
raw_reads_quality_control_script_name = '{}/01_raw_reads_quality_control_script.sh'.format(raw_reads_quality_control_dir)
if not os.path.exists(raw_reads_quality_control_dir):
    os.makedirs(raw_reads_quality_control_dir)
raw_reads_quality_control_script = open(raw_reads_quality_control_script_name, 'w')

raw_reads_quality_control_script.write('#!/bin/bash\n\n')
raw_reads_quality_control_script.write('set -euxo pipefail\n\n')

if start_from_bam == False:

    # Bash commands to call fastqc to run quality control analysis on all raw fq.gz
    if read_mode == 'single':
        for list_counter in range(len(chip_name)):
            raw_reads_quality_control_script.write('fastqc {} -t {} {}/{}.fq.gz -o {}\n\n'.format(
                fastqc1_arg, 
                cpu_count, 
                raw_data_dir, 
                chip_name[list_counter], 
                raw_reads_quality_control_dir))

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                raw_reads_quality_control_script.write('fastqc {} -t {} {}/{}.fq.gz -o {}\n\n'.format(
                    fastqc1_arg, 
                    cpu_count, 
                    raw_data_dir, 
                    ctrl_name[list_counter], 
                    raw_reads_quality_control_dir))
                    


    elif read_mode == 'paired':
        for list_counter in range(len(chip_r1_name)):
            raw_reads_quality_control_script.write('fastqc {} -t {} {}/{}.fq.gz -o {}\n\n'.format(
                fastqc1_arg, 
                cpu_count, 
                raw_data_dir, 
                chip_r1_name[list_counter], 
                raw_reads_quality_control_dir))

        for list_counter in range(len(chip_r2_name)):
            raw_reads_quality_control_script.write('fastqc {} -t {} {}/{}.fq.gz -o {}\n\n'.format(
                fastqc1_arg, 
                cpu_count, 
                raw_data_dir, 
                chip_r2_name[list_counter], 
                raw_reads_quality_control_dir))

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_r1_name)):
                raw_reads_quality_control_script.write('fastqc {} -t {} {}/{}.fq.gz -o {}\n\n'.format(
                    fastqc1_arg, 
                    cpu_count, 
                    raw_data_dir, 
                    ctrl_r1_name[list_counter], 
                    raw_reads_quality_control_dir))

            for list_counter in range(len(ctrl_r2_name)):
                raw_reads_quality_control_script.write('fastqc {} -t {} {}/{}.fq.gz -o {}\n\n'.format(
                    fastqc1_arg, 
                    cpu_count, 
                    raw_data_dir, 
                    ctrl_r2_name[list_counter], 
                    raw_reads_quality_control_dir))

raw_reads_quality_control_script.close() # Closing the script '01_raw_reads_quality_control_script.sh'. Flushing the write buffer



####################################################################################################################
### 02_deduplicating_script.sh
####################################################################################################################

# Create a directory "02_deduplicating" for the deduplicated sequencing reads, 
#   free of optical and tile-edge duplicate reads
# Create a script "02_deduplicating_script.sh" within, that calls for duplicates removal program: clumpify.sh
# In: filename.fq.gz (output from 00_raw_data.sh)
# Out: filename.deduped.fq.gz (input for 02_adapter_trimming_script.sh)

deduplicating_dir = '{}/02_deduplicating'.format(output_dir)
deduplicating_script_name = '{}/02_deduplicating_script.sh'.format(deduplicating_dir)
if not os.path.exists(deduplicating_dir + '/logs'):
    os.makedirs(deduplicating_dir + '/logs')
deduplicating_script = open(deduplicating_script_name, 'w')

deduplicating_script.write('#!/bin/bash\n\n')
deduplicating_script.write('set -euxo pipefail\n\n')

if start_from_bam == False:

    # Bash commands to call clumpify.sh to remove optical and tile-edge duplicate reads from the fq.gz files
    if read_mode == 'single':
        for list_counter in range(len(chip_name)):
            deduplicating_script.write('clumpify.sh in={}/{}.fq.gz out={}/{}.deduped.fq.gz {} 1>{}/logs/{}.deduplicating.out 2>{}/logs/{}.deduplicating.err\n\n'.format(
                raw_data_dir, 
                chip_name[list_counter], 
                deduplicating_dir, 
                chip_name[list_counter], 
                clumpify_arg, 
                deduplicating_dir, 
                chip_name[list_counter], 
                deduplicating_dir, 
                chip_name[list_counter]))

        if args.deltemp:
            # Immediately deletes the raw fastq after deduplication
            for list_counter in range(len(chip_name)):
                deduplicating_script.write('rm -r -f {}/{}.fq.gz &\n'.format(
                    raw_data_dir, 
                    chip_name[list_counter]))


                    
        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                deduplicating_script.write('clumpify.sh in={}/{}.fq.gz out={}/{}.deduped.fq.gz {} 1>{}/logs/{}.deduplicating.out 2>{}/logs/{}.deduplicating.err\n\n'.format(
                    raw_data_dir, 
                    ctrl_name[list_counter], 
                    deduplicating_dir, 
                    ctrl_name[list_counter], 
                    clumpify_arg, 
                    deduplicating_dir, 
                    ctrl_name[list_counter], 
                    deduplicating_dir, 
                    ctrl_name[list_counter]))

            if args.deltemp:        
                for list_counter in range(len(ctrl_name)):
                    deduplicating_script.write('rm -r -f {}/{}.fq.gz &\n'.format(
                        raw_data_dir, 
                        ctrl_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                        deduplicating_script.write('\nwait\n\n')
                        


    elif read_mode == 'paired':
        for list_counter in range(len(chip_name)):
            deduplicating_script.write('clumpify.sh in={}/{}.fq.gz in2={}/{}.fq.gz out={}/{}.deduped.fq.gz out2={}/{}.deduped.fq.gz {} 1>{}/logs/{}.deduplicating.out 2>{}/logs/{}.deduplicating.err\n\n'.format(
                raw_data_dir, 
                chip_r1_name[list_counter], 
                raw_data_dir, 
                chip_r2_name[list_counter], 
                deduplicating_dir, 
                chip_r1_name[list_counter], 
                deduplicating_dir, 
                chip_r2_name[list_counter], 
                clumpify_arg, 
                deduplicating_dir, 
                chip_name[list_counter], 
                deduplicating_dir, 
                chip_name[list_counter]))

        if args.deltemp:        
            for list_counter in range(len(chip_name)):
                deduplicating_script.write('rm -r -f {}/{}.fq.gz {}/{}.fq.gz &\n'.format(
                    raw_data_dir, 
                    chip_r1_name[list_counter], 
                    raw_data_dir, 
                    chip_r2_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                    deduplicating_script.write('\nwait\n\n')
                    
        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                deduplicating_script.write('clumpify.sh in={}/{}.fq.gz in2={}/{}.fq.gz out={}/{}.deduped.fq.gz out2={}/{}.deduped.fq.gz {} 1>{}/logs/{}.deduplicating.out 2>{}/logs/{}.deduplicating.err\n\n'.format(
                    raw_data_dir, 
                    ctrl_r1_name[list_counter], 
                    raw_data_dir, 
                    ctrl_r2_name[list_counter], 
                    deduplicating_dir, 
                    ctrl_r1_name[list_counter], 
                    deduplicating_dir, 
                    ctrl_r2_name[list_counter], 
                    clumpify_arg, 
                    deduplicating_dir, 
                    ctrl_name[list_counter], 
                    deduplicating_dir, 
                    ctrl_name[list_counter]))
                    
            if args.deltemp:        
                for list_counter in range(len(ctrl_name)):
                    deduplicating_script.write('rm -r -f {}/{}.fq.gz {}/{}.fq.gz &\n'.format(
                        raw_data_dir, 
                        ctrl_r1_name[list_counter], 
                        raw_data_dir, 
                        ctrl_r2_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                        deduplicating_script.write('\nwait\n\n')
                        
deduplicating_script.close() # Closing the script '02_deduplicating_script.sh'. Flushing the write buffer



####################################################################################################################
### 03_adapter_trimming_script.sh
####################################################################################################################

# Create a directory "03_adapter_trimming" for the adapter-trimmed sequencing reads, 
#   free of left-flanking adapter sequences 
# Create a script "03_adapter_trimming_script.sh" within, that calls for adapter sequence trimming program: bbduk.sh
# In: filename.deduped.fq.gz (output from 03_deduplicating_script.sh)
# Out: filename.adaptertrimmed.fq.gz (input for 03_quality_trimming_script.sh)

adapter_trimming_dir = '{}/03_adapter_trimming'.format(output_dir)
adapter_trimming_script_name = '{}/03_adapter_trimming_script.sh'.format(adapter_trimming_dir)
if not os.path.exists(adapter_trimming_dir + '/logs'):
    os.makedirs(adapter_trimming_dir + '/logs')
adapter_trimming_script = open(adapter_trimming_script_name, 'w')

adapter_trimming_script.write('#!/bin/bash\n\n')
adapter_trimming_script.write('set -euxo pipefail\n\n')

if start_from_bam == False:

    # Bash commands to call clumpify.sh to trim adapter sequences from the left flank of reads in the fq.gz files
    if read_mode == 'single':
        for list_counter in range(len(chip_name)):
            adapter_trimming_script.write('bbduk.sh in={}/{}.deduped.fq.gz out={}/{}.adaptertrimmed.fq.gz ref={}/bbmap/adapters.fa {} overwrite=t 1>{}/logs/{}.adapter_trimming.out 2>{}/logs/{}.adapter_trimming.err\n\n'.format(
                deduplicating_dir, 
                chip_name[list_counter], 
                adapter_trimming_dir, 
                chip_name[list_counter], 
                genome_dir, 
                bbduk_arg, 
                adapter_trimming_dir, 
                chip_name[list_counter], 
                adapter_trimming_dir, 
                chip_name[list_counter]))

        if args.deltemp:
            # Immediately deletes the deduplicated fastq after adapter trimming
            for list_counter in range(len(chip_name)):
                adapter_trimming_script.write('rm -r -f {}/{}.deduped.fq.gz &\n'.format(
                    deduplicating_dir, 
                    chip_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                    adapter_trimming_script.write('\nwait\n\n')
                    
        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                adapter_trimming_script.write('bbduk.sh in={}/{}.deduped.fq.gz out={}/{}.adaptertrimmed.fq.gz ref={}/bbmap/adapters.fa {} overwrite=t 1>{}/logs/{}.adapter_trimming.out 2>{}/logs/{}.adapter_trimming.err\n\n'.format(
                    deduplicating_dir, 
                    ctrl_name[list_counter], 
                    adapter_trimming_dir, 
                    ctrl_name[list_counter], 
                    genome_dir, 
                    bbduk_arg, 
                    adapter_trimming_dir, 
                    ctrl_name[list_counter], 
                    adapter_trimming_dir, 
                    ctrl_name[list_counter]))

            if args.deltemp:        
                for list_counter in range(len(ctrl_name)):
                    adapter_trimming_script.write('rm -r -f {}/{}.deduped.fq.gz &\n'.format(
                        deduplicating_dir, 
                        ctrl_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                        adapter_trimming_script.write('\nwait\n\n')


    if read_mode == 'paired':
        for list_counter in range(len(chip_name)):
            adapter_trimming_script.write('bbduk.sh in={}/{}.deduped.fq.gz in2={}/{}.deduped.fq.gz out={}/{}.adaptertrimmed.fq.gz out2={}/{}.adaptertrimmed.fq.gz ref={}/bbmap/adapters.fa {} overwrite=t 1>{}/logs/{}.adapter_trimming.out 2>{}/logs/{}.adapter_trimming.err\n\n'.format(
                deduplicating_dir, 
                chip_r1_name[list_counter], 
                deduplicating_dir, 
                chip_r2_name[list_counter], 
                adapter_trimming_dir, 
                chip_r1_name[list_counter], 
                adapter_trimming_dir, 
                chip_r2_name[list_counter], 
                genome_dir, 
                bbduk_arg, 
                adapter_trimming_dir, 
                chip_name[list_counter], 
                adapter_trimming_dir, 
                chip_name[list_counter]))

        if args.deltemp:        
            for list_counter in range(len(chip_name)):
                adapter_trimming_script.write('rm -r -f {}/{}.deduped.fq.gz {}/{}.deduped.fq.gz &\n'.format(
                    deduplicating_dir, 
                    chip_r1_name[list_counter], 
                    deduplicating_dir, 
                    chip_r2_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                    adapter_trimming_script.write('\nwait\n\n')

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                adapter_trimming_script.write('bbduk.sh in={}/{}.deduped.fq.gz in2={}/{}.deduped.fq.gz out={}/{}.adaptertrimmed.fq.gz out2={}/{}.adaptertrimmed.fq.gz ref={}/bbmap/adapters.fa {} overwrite=t 1>{}/logs/{}.adapter_trimming.out 2>{}/logs/{}.adapter_trimming.err\n\n'.format(
                    deduplicating_dir, 
                    ctrl_r1_name[list_counter], 
                    deduplicating_dir, 
                    ctrl_r2_name[list_counter], 
                    adapter_trimming_dir, 
                    ctrl_r1_name[list_counter], 
                    adapter_trimming_dir, 
                    ctrl_r2_name[list_counter], 
                    genome_dir, 
                    bbduk_arg, 
                    adapter_trimming_dir, 
                    ctrl_name[list_counter], 
                    adapter_trimming_dir, 
                    ctrl_name[list_counter]))

            if args.deltemp:        
                for list_counter in range(len(ctrl_name)):
                    adapter_trimming_script.write('rm -r -f {}/{}.deduped.fq.gz {}/{}.deduped.fq.gz &\n'.format(
                        deduplicating_dir, 
                        ctrl_r1_name[list_counter], 
                        deduplicating_dir, 
                        ctrl_r2_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                        adapter_trimming_script.write('\nwait\n\n')
                            
adapter_trimming_script.close() # Closing the script '03_adapter_trimming_script.sh'. Flushing the write buffer



####################################################################################################################
### 04_quality_trimming_script.sh
####################################################################################################################

# Create a directory "04_quality_trimming" for the quality-trimmed sequencing reads, 
#   free of sequence reads with low confidence bases 
# Create a script "04_quality_trimming_script.sh" within, 
#   that calls for phred-score-based quality trimming program: trimmomatic
# In: filename.adaptertrimmed.fq.gz (output from 02_adapter_trimming_script.sh)
# Out: filename.qualitytrimmed.fq.gz (input for 06_bwa_mem_aligning_script.sh)

quality_trimming_dir = '{}/04_quality_trimming'.format(output_dir)
quality_trimming_script_name = '{}/04_quality_trimming_script.sh'.format(quality_trimming_dir)
if not os.path.exists(quality_trimming_dir + '/logs'):
    os.makedirs(quality_trimming_dir + '/logs')

quality_trimming_script = open(quality_trimming_script_name, 'w')

quality_trimming_script.write('#!/bin/bash\n\n')
quality_trimming_script.write('set -euxo pipefail\n\n')

if start_from_bam == False:

    # Bash commands to call clumpify.sh to trim the bases in sequence reads from 
    #   both left and right flanks when their phred score is lower than 20
    # Also runs a moving scan to cut apart the sequence reads where the average 
    #   phred score of 20 bp lenght window is lower than 20
    # Finally, it discards all sequence reads with post-trimming lenght below 20 bp

    if read_mode == 'single':
        for list_counter in range(len(chip_name)):
            quality_trimming_script.write('java -jar {} SE -threads {} -phred33 {}/{}.adaptertrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz {} 1>{}/logs/{}.quality_trimming.out 2>{}/logs/{}.quality_trimming.err\n\n'.format(
                trimmomatic_full_path, 
                cpu_count, 
                adapter_trimming_dir, 
                chip_name[list_counter], 
                quality_trimming_dir, 
                chip_name[list_counter], 
                trimmomatic_arg, 
                quality_trimming_dir, 
                chip_name[list_counter], 
                quality_trimming_dir, 
                chip_name[list_counter]))

        if args.deltemp:        
            # Immediately deletes the adapter trimmed fastq after quality trimming
            for list_counter in range(len(chip_name)):
                quality_trimming_script.write('rm -r -f {}/{}.adaptertrimmed.fq.gz &\n'.format(
                    adapter_trimming_dir, 
                    chip_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                    quality_trimming_script.write('\nwait\n\n')
                    
        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                quality_trimming_script.write('java -jar {} SE -threads {} -phred33 {}/{}.adaptertrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz {} 1>{}/logs/{}.quality_trimming.out 2>{}/logs/{}.quality_trimming.err\n\n'.format(
                    trimmomatic_full_path, 
                    cpu_count, 
                    adapter_trimming_dir, 
                    ctrl_name[list_counter], 
                    quality_trimming_dir, 
                    ctrl_name[list_counter], 
                    trimmomatic_arg, 
                    quality_trimming_dir, 
                    ctrl_name[list_counter], 
                    quality_trimming_dir, 
                    ctrl_name[list_counter]))

            if args.deltemp:        
                for list_counter in range(len(ctrl_name)):
                    quality_trimming_script.write('rm -r -f {}/{}.adaptertrimmed.fq.gz &\n'.format(
                        adapter_trimming_dir, 
                        ctrl_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                        quality_trimming_script.write('\nwait\n\n')
                        


    # When run in paired-end reads mode, trimmomatic also separates paired reads 
    #   (labeled .quality_trimmed here) from orphan reads (labeled .unpaired here)
    if read_mode == 'paired':
        for list_counter in range(len(chip_name)):
            quality_trimming_script.write('java -jar {} PE -threads {} -phred33 {}/{}.adaptertrimmed.fq.gz {}/{}.adaptertrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz {}/{}.unpaired.fq.gz {}/{}.qualitytrimmed.fq.gz {}/{}.unpaired.fq.gz {} 1>{}/logs/{}.quality_trimming.out 2>{}/logs/{}.quality_trimming.err\n\n'.format(
                trimmomatic_full_path, 
                cpu_count, 
                adapter_trimming_dir, 
                chip_r1_name[list_counter], 
                adapter_trimming_dir, 
                chip_r2_name[list_counter], 
                quality_trimming_dir, 
                chip_r1_name[list_counter], 
                quality_trimming_dir, 
                chip_r1_name[list_counter], 
                quality_trimming_dir, 
                chip_r2_name[list_counter], 
                quality_trimming_dir, 
                chip_r2_name[list_counter], 
                trimmomatic_arg, 
                quality_trimming_dir, 
                chip_name[list_counter], 
                quality_trimming_dir, 
                chip_name[list_counter]))

        if args.deltemp:        
            for list_counter in range(len(chip_name)):
                quality_trimming_script.write('rm -r -f {}/{}.adaptertrimmed.fq.gz {}/{}.adaptertrimmed.fq.gz &\n'.format(
                    adapter_trimming_dir, 
                    chip_r1_name[list_counter], 
                    adapter_trimming_dir, 
                    chip_r2_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                    quality_trimming_script.write('\nwait\n\n')
                    
        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                quality_trimming_script.write('java -jar {} PE -threads {} -phred33 {}/{}.adaptertrimmed.fq.gz {}/{}.adaptertrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz {}/{}.unpaired.fq.gz {}/{}.qualitytrimmed.fq.gz {}/{}.unpaired.fq.gz {} 1>{}/logs/{}.quality_trimming.out 2>{}/logs/{}.quality_trimming.err\n\n'.format(
                    trimmomatic_full_path, 
                    cpu_count, 
                    adapter_trimming_dir, 
                    ctrl_r1_name[list_counter], 
                    adapter_trimming_dir, 
                    ctrl_r2_name[list_counter], 
                    quality_trimming_dir, 
                    ctrl_r1_name[list_counter], 
                    quality_trimming_dir, 
                    ctrl_r1_name[list_counter], 
                    quality_trimming_dir, 
                    ctrl_r2_name[list_counter], 
                    quality_trimming_dir, 
                    ctrl_r2_name[list_counter], 
                    trimmomatic_arg, 
                    quality_trimming_dir, 
                    ctrl_name[list_counter], 
                    quality_trimming_dir, 
                    ctrl_name[list_counter]))

            if args.deltemp:        
                for list_counter in range(len(ctrl_name)):
                    quality_trimming_script.write('rm -r -f {}/{}.adaptertrimmed.fq.gz {}/{}.adaptertrimmed.fq.gz &\n'.format(
                        adapter_trimming_dir, 
                        ctrl_r1_name[list_counter], 
                        adapter_trimming_dir, 
                        ctrl_r2_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                        quality_trimming_script.write('\nwait\n\n')

quality_trimming_script.close() # Closing the script '04_quality_trimming_script.sh'. Flushing the write buffer



####################################################################################################################
### 05_preprocessed_reads_quality_control_script.sh
####################################################################################################################

# Create a directory "05_preprocessed_reads_quality_control" for the quality control 
#   report files of the deduplicated and trimmed sequencing reads
# Create a script "05_preprocessed_reads_quality_control_script.sh" within, 
#   that calls for sequencing reads QC analysis program: fastqc
# In: filename.qualitytrimmed.fq.gz (output from 03_quality_trimming_script.sh)
# Out: filename_fastqc.html (Reports. not to be further processed downstreams)

preprocessed_reads_quality_control_dir = '{}/05_preprocessed_reads_quality_control'.format(output_dir)
preprocessed_reads_quality_control_script_name = '{}/05_preprocessed_reads_quality_control_script.sh'.format(preprocessed_reads_quality_control_dir)
if not os.path.exists(preprocessed_reads_quality_control_dir):
    os.makedirs(preprocessed_reads_quality_control_dir)
preprocessed_reads_quality_control_script = open(preprocessed_reads_quality_control_script_name, 'w')

preprocessed_reads_quality_control_script.write('#!/bin/bash\n\n')
preprocessed_reads_quality_control_script.write('set -euxo pipefail\n\n')

if start_from_bam == False:

    # Bash commands to call fastqc to run quality control analysis on all deduplicated and trimmed fq.gz
    if read_mode == 'single':
        for list_counter in range(len(chip_name)):
            preprocessed_reads_quality_control_script.write('fastqc {} -t {} {}/{}.qualitytrimmed.fq.gz -o {}\n\n'.format(
                fastqc2_arg, 
                cpu_count, 
                quality_trimming_dir, 
                chip_name[list_counter], 
                preprocessed_reads_quality_control_dir))

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                preprocessed_reads_quality_control_script.write('fastqc {} -t {} {}/{}.qualitytrimmed.fq.gz -o {}\n\n'.format(
                    fastqc2_arg, 
                    cpu_count, 
                    quality_trimming_dir, 
                    ctrl_name[list_counter], 
                    preprocessed_reads_quality_control_dir))



    elif read_mode == 'paired':
        for list_counter in range(len(chip_r1_name)):
            preprocessed_reads_quality_control_script.write('fastqc {} -t {} {}/{}.qualitytrimmed.fq.gz -o {}\n\n'.format(
                fastqc2_arg, 
                cpu_count, 
                quality_trimming_dir, 
                chip_r1_name[list_counter], 
                preprocessed_reads_quality_control_dir))
                
        for list_counter in range(len(chip_r2_name)):
            preprocessed_reads_quality_control_script.write('fastqc {} -t {} {}/{}.qualitytrimmed.fq.gz -o {}\n\n'.format(
                fastqc2_arg, 
                cpu_count, 
                quality_trimming_dir, 
                chip_r2_name[list_counter], 
                preprocessed_reads_quality_control_dir))

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_r1_name)):
                preprocessed_reads_quality_control_script.write('fastqc {} -t {} {}/{}.qualitytrimmed.fq.gz -o {}\n\n'.format(
                    fastqc2_arg, 
                    cpu_count, 
                    quality_trimming_dir, 
                    ctrl_r1_name[list_counter], 
                    preprocessed_reads_quality_control_dir))

            for list_counter in range(len(ctrl_r2_name)):
                preprocessed_reads_quality_control_script.write('fastqc {} -t {} {}/{}.qualitytrimmed.fq.gz -o {}\n\n'.format(
                    fastqc2_arg, 
                    cpu_count, 
                    quality_trimming_dir, 
                    ctrl_r2_name[list_counter], 
                    preprocessed_reads_quality_control_dir))

preprocessed_reads_quality_control_script.close() # Closing the script '05_preprocessed_reads_quality_control_script.sh'. Flushing the write buffer



####################################################################################################################
### 06_bwa_mem_aligning_script.sh
####################################################################################################################

# Create a directory "06_bwa_mem_aligning" for the aligned reads to the reference genome
# Create a script "06_bwa_mem_aligning_script.sh" within, that calls for sequence reads aligning program: bwa mem
# The mem algorithm was chosen due to its significantly better performance in aligning shorter reads
# Per version 4.0, the alignments are pipe-sorted to samtools sort, resulting in .bam output instead of .sam
#   Saves space: .bam files are typically 4-6 times smaller in size than .sam files
#   Needs index: unlike .sam files, .bam files needs to be indexed before MAPQ filtering by samtools view
#       Per version 4.0, samtools index are performed at the end of this step to generate indices for the aligned bam files
# In: filename.qualitytrimmed.fq.gz (output from 03_quality_trimming_script.sh)
# Out: filename.mapped.bam (input for 07_MAPQ_filtering_script.sh)
# Out: filename.mapped.bam.bai (indices for the aligned bam file. Input for 07_MAPQ_filtering_script.sh)

bwa_mem_aligning_dir = '{}/06_bwa_mem_aligning'.format(output_dir)
bwa_mem_aligning_script_name = '{}/06_bwa_mem_aligning_script.sh'.format(bwa_mem_aligning_dir)
if not os.path.exists(bwa_mem_aligning_dir + '/logs'):
    os.makedirs(bwa_mem_aligning_dir + '/logs')
bwa_mem_aligning_script = open(bwa_mem_aligning_script_name, 'w')

bwa_mem_aligning_script.write('#!/bin/bash\n\n')
bwa_mem_aligning_script.write('set -euxo pipefail\n\n')

if start_from_bam == False:

    # Bash commands to call bwa with mem algorithm to align all sequence reads in 
    #   the preprocessed .fz.gz files to bwa reference genome
    if read_mode == 'single':
        for list_counter in range(len(chip_name)):
            bwa_mem_aligning_script.write('bwa mem {} -t {} {}/bwa/{}.fa {}/{}.qualitytrimmed.fq.gz 2> {}/logs/{}.bwamemaligning.err | samtools sort -@ {} -o {}/{}.aligned.bam 2> /dev/null\n\n'.format(
                bwa_mem_arg, 
                cpu_count, 
                genome_dir, 
                genome_ref, 
                quality_trimming_dir, 
                chip_name[list_counter], 
                bwa_mem_aligning_dir, 
                chip_name[list_counter],
                cpu_count, 
                bwa_mem_aligning_dir, 
                chip_name[list_counter]))

        if norm_ref != 'none':
            for list_counter in range(len(chip_name)):
                bwa_mem_aligning_script.write('bwa mem {} -t {} {}/bwa/{}.fa {}/{}.qualitytrimmed.fq.gz 2> {}/logs/{}.normaligning.err | samtools sort -@ {} -o {}/{}.normaligned.bam 2> /dev/null\n\n'.format(
                    bwa_mem_arg, 
                    cpu_count, 
                    genome_dir, 
                    norm_ref, 
                    quality_trimming_dir, 
                    chip_name[list_counter], 
                    bwa_mem_aligning_dir, 
                    chip_name[list_counter],
                    cpu_count, 
                    bwa_mem_aligning_dir, 
                    chip_name[list_counter]))

        if args.deltemp:        
            # Immediately deletes the quality trimmed fastq after alignment to genome
            for list_counter in range(len(chip_name)):
                bwa_mem_aligning_script.write('rm -r -f {}/{}.qualitytrimmed.fq.gz &\n'.format(
                    quality_trimming_dir, 
                    chip_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                    bwa_mem_aligning_script.write('\nwait\n\n')

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                bwa_mem_aligning_script.write('bwa mem {} -t {} {}/bwa/{}.fa {}/{}.qualitytrimmed.fq.gz 2> {}/logs/{}.bwamemaligning.err | samtools sort -@ {} -o {}/{}.aligned.bam 2> /dev/null\n\n'.format(
                    bwa_mem_arg, 
                    cpu_count, 
                    genome_dir, 
                    genome_ref, 
                    quality_trimming_dir, 
                    ctrl_name[list_counter],
                    bwa_mem_aligning_dir, 
                    ctrl_name[list_counter],
                    cpu_count, 
                    bwa_mem_aligning_dir, 
                    ctrl_name[list_counter]))

            if norm_ref != 'none':
                for list_counter in range(len(ctrl_name)):
                    bwa_mem_aligning_script.write('bwa mem {} -t {} {}/bwa/{}.fa {}/{}.qualitytrimmed.fq.gz 2> {}/logs/{}.normaligning.err | samtools sort -@ {} -o {}/{}.normaligned.bam 2> /dev/null\n\n'.format(
                        bwa_mem_arg, 
                        cpu_count, 
                        genome_dir, 
                        norm_ref, 
                        quality_trimming_dir, 
                        ctrl_name[list_counter],
                        bwa_mem_aligning_dir, 
                        ctrl_name[list_counter],
                        cpu_count, 
                        bwa_mem_aligning_dir, 
                        ctrl_name[list_counter]))

            if args.deltemp:        
                for list_counter in range(len(ctrl_name)):
                    bwa_mem_aligning_script.write('rm -r -f {}/{}.qualitytrimmed.fq.gz &\n'.format(
                        quality_trimming_dir, 
                        ctrl_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                        bwa_mem_aligning_script.write('\nwait\n\n')

    # Only the sequence reads determined as paired by trimmomatic algorithm are processed 
    #   for bwa mem alignment. Orphans are ignored.
    elif read_mode == 'paired':
        for list_counter in range(len(chip_name)):
            bwa_mem_aligning_script.write('bwa mem {} -t {} {}/bwa/{}.fa {}/{}.qualitytrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz 2> {}/logs/{}.bwamemaligning.err | samtools sort -@ {} -o {}/{}.aligned.bam 2> /dev/null\n\n'.format(
                bwa_mem_arg, 
                math.ceil(cpu_count/2), 
                genome_dir, 
                genome_ref, 
                quality_trimming_dir, 
                chip_r1_name[list_counter], 
                quality_trimming_dir, 
                chip_r2_name[list_counter],
                bwa_mem_aligning_dir, 
                chip_name[list_counter], 
                math.ceil(cpu_count/2),  
                bwa_mem_aligning_dir, 
                chip_name[list_counter]))

        if norm_ref != 'none':
            for list_counter in range(len(chip_name)):
                bwa_mem_aligning_script.write('bwa mem {} -t {} {}/bwa/{}.fa {}/{}.qualitytrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz 2> {}/logs/{}.normaligning.err | samtools sort -@ {} -o {}/{}.normaligned.bam 2> /dev/null\n\n'.format(
                    bwa_mem_arg, 
                    math.ceil(cpu_count/2), 
                    genome_dir, 
                    norm_ref, 
                    quality_trimming_dir, 
                    chip_r1_name[list_counter], 
                    quality_trimming_dir, 
                    chip_r2_name[list_counter],
                    bwa_mem_aligning_dir, 
                    chip_name[list_counter], 
                    math.ceil(cpu_count/2),  
                    bwa_mem_aligning_dir, 
                    chip_name[list_counter]))

        if args.deltemp:        
            for list_counter in range(len(chip_name)):
                bwa_mem_aligning_script.write('rm -r -f {}/{}.qualitytrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz {}/{}.unpaired.fq.gz {}/{}.unpaired.fq.gz &\n'.format(
                    quality_trimming_dir, 
                    chip_r1_name[list_counter], 
                    quality_trimming_dir, 
                    chip_r2_name[list_counter], 
                    quality_trimming_dir, 
                    chip_r1_name[list_counter], 
                    quality_trimming_dir, 
                    chip_r2_name[list_counter]))
                
                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                    bwa_mem_aligning_script.write('\nwait\n\n')

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                bwa_mem_aligning_script.write('bwa mem {} -t {} {}/bwa/{}.fa {}/{}.qualitytrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz 2> {}/logs/{}.bwamemaligning.err | samtools sort -@ {} -o {}/{}.aligned.bam 2> /dev/null\n\n'.format(
                    bwa_mem_arg, 
                    math.ceil(cpu_count/2), 
                    genome_dir, 
                    genome_ref, 
                    quality_trimming_dir, 
                    ctrl_r1_name[list_counter], 
                    quality_trimming_dir, 
                    ctrl_r2_name[list_counter],
                    bwa_mem_aligning_dir, 
                    ctrl_name[list_counter],  
                    math.ceil(cpu_count/2), 
                    bwa_mem_aligning_dir, 
                    ctrl_name[list_counter]))

            if norm_ref != 'none':
                for list_counter in range(len(ctrl_name)):
                    bwa_mem_aligning_script.write('bwa mem {} -t {} {}/bwa/{}.fa {}/{}.qualitytrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz 2> {}/logs/{}.normaligning.err | samtools sort -@ {} -o {}/{}.normaligned.bam 2> /dev/null\n\n'.format(
                        bwa_mem_arg, 
                        math.ceil(cpu_count/2), 
                        genome_dir, 
                        norm_ref, 
                        quality_trimming_dir, 
                        ctrl_r1_name[list_counter], 
                        quality_trimming_dir, 
                        ctrl_r2_name[list_counter],
                        bwa_mem_aligning_dir, 
                        ctrl_name[list_counter],  
                        math.ceil(cpu_count/2), 
                        bwa_mem_aligning_dir, 
                        ctrl_name[list_counter]))

            if args.deltemp:        
                for list_counter in range(len(ctrl_name)):
                    bwa_mem_aligning_script.write('rm -r -f {}/{}.qualitytrimmed.fq.gz {}/{}.qualitytrimmed.fq.gz {}/{}.unpaired.fq.gz {}/{}.unpaired.fq.gz &\n'.format(
                        quality_trimming_dir, 
                        ctrl_r1_name[list_counter], 
                        quality_trimming_dir, 
                        ctrl_r2_name[list_counter], 
                        quality_trimming_dir, 
                        ctrl_r1_name[list_counter], 
                        quality_trimming_dir, 
                        ctrl_r2_name[list_counter]))

                    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                        bwa_mem_aligning_script.write('\nwait\n\n')



    # Writing bash commands to generate an index file for all individual .bam files
    for list_counter in range(len(chip_name)):
        bwa_mem_aligning_script.write('samtools view -F 4 {}/{}.aligned.bam -bS > {}/{}.mapped.bam &\n'.format(
            bwa_mem_aligning_dir, 
            chip_name[list_counter],
            bwa_mem_aligning_dir, 
            chip_name[list_counter]))

        if ((list_counter + 1) % int(cpu_count / 2)) == 0 or list_counter == (len(chip_name) - 1):
            bwa_mem_aligning_script.write('\nwait\n\n')
            
    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        for list_counter in range(len(ctrl_name)):
            bwa_mem_aligning_script.write('samtools view -F 4 {}/{}.aligned.bam -bS > {}/{}.mapped.bam &\n'.format(
                bwa_mem_aligning_dir, 
                ctrl_name[list_counter],
                bwa_mem_aligning_dir, 
                ctrl_name[list_counter]))
    
            if ((list_counter + 1) % int(cpu_count / 2)) == 0 or list_counter == (len(ctrl_name) - 1):
                bwa_mem_aligning_script.write('\nwait\n\n')

    if norm_ref != 'none':
        for list_counter in range(len(chip_name)):
            bwa_mem_aligning_script.write('samtools view -F 4 {}/{}.normaligned.bam -bS > {}/{}.normmapped.bam &\n'.format(
                bwa_mem_aligning_dir, 
                chip_name[list_counter],
                bwa_mem_aligning_dir, 
                chip_name[list_counter]))

            if ((list_counter + 1) % int(cpu_count / 2)) == 0 or list_counter == (len(chip_name) - 1):
                bwa_mem_aligning_script.write('\nwait\n\n')
                
        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                bwa_mem_aligning_script.write('samtools view -F 4 {}/{}.normaligned.bam -bS > {}/{}.normmapped.bam &\n'.format(
                    bwa_mem_aligning_dir, 
                    ctrl_name[list_counter],
                    bwa_mem_aligning_dir, 
                    ctrl_name[list_counter]))
        
                if ((list_counter + 1) % int(cpu_count / 2)) == 0 or list_counter == (len(ctrl_name) - 1):
                    bwa_mem_aligning_script.write('\nwait\n\n')



    # Writing bash commands to generate an index file for all individual .bam files
    for list_counter in range(len(chip_name)):
        bwa_mem_aligning_script.write('rm -r -f {}/{}.aligned.bam &\n'.format(
            bwa_mem_aligning_dir, 
            chip_name[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            bwa_mem_aligning_script.write('\nwait\n\n')
            
    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        for list_counter in range(len(ctrl_name)):
            bwa_mem_aligning_script.write('rm -r -f {}/{}.aligned.bam &\n'.format(
                bwa_mem_aligning_dir, 
                ctrl_name[list_counter]))
    
            if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                bwa_mem_aligning_script.write('\nwait\n\n')

    if norm_ref != 'none':
        for list_counter in range(len(chip_name)):
            bwa_mem_aligning_script.write('rm -r -f {}/{}.normaligned.bam &\n'.format(
                bwa_mem_aligning_dir, 
                chip_name[list_counter]))

            if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                bwa_mem_aligning_script.write('\nwait\n\n')
                
        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                bwa_mem_aligning_script.write('rm -r -f {}/{}.normaligned.bam &\n'.format(
                    bwa_mem_aligning_dir, 
                    ctrl_name[list_counter]))
        
                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                    bwa_mem_aligning_script.write('\nwait\n\n')



    # Writing bash commands to generate an index file for all individual .bam files
    for list_counter in range(len(chip_name)):
        bwa_mem_aligning_script.write('samtools index -@ {} {}/{}.mapped.bam &\n'.format(
            cpu_count, 
            bwa_mem_aligning_dir, 
            chip_name[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            bwa_mem_aligning_script.write('\nwait\n\n')
            
    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        for list_counter in range(len(ctrl_name)):
            bwa_mem_aligning_script.write('samtools index -@ {} {}/{}.mapped.bam &\n'.format(
                cpu_count, 
                bwa_mem_aligning_dir, 
                ctrl_name[list_counter]))
    
            if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                bwa_mem_aligning_script.write('\nwait\n\n')

    if norm_ref != 'none':
        for list_counter in range(len(chip_name)):
            bwa_mem_aligning_script.write('samtools index -@ {} {}/{}.normmapped.bam &\n'.format(
                cpu_count, 
                bwa_mem_aligning_dir, 
                chip_name[list_counter]))

            if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                bwa_mem_aligning_script.write('\nwait\n\n')
                
        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            for list_counter in range(len(ctrl_name)):
                bwa_mem_aligning_script.write('samtools index -@ {} {}/{}.normmapped.bam &\n'.format(
                    cpu_count, 
                    bwa_mem_aligning_dir, 
                    ctrl_name[list_counter]))
        
                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                    bwa_mem_aligning_script.write('\nwait\n\n')

bwa_mem_aligning_script.close() # Closing the script '06_bwa_mem_aligning_script.sh'. Flushing the write buffer



####################################################################################################################
### 07_MAPQ_filtering_script.sh
####################################################################################################################

# Create a directory "07_MAPQ_filtering" for the filtered reads, free of sequence reads with low alignment confidence
# Create a script "07_MAPQ_filtering_script.sh" within, that calls for .bam file reading program 
#   that enables MAPQ score thresholding and removes scaffold chromosomes: samtools view
# In: filename.mapped.bam (output from 06_bwa_mem_aligning_script.sh)
# Out: filename.mapqfiltered.bam (input for 08_results_script.sh)

mapq_filtering_dir = '{}/07_MAPQ_filtering'.format(output_dir)
mapq_filtering_script_name = '{}/07_MAPQ_filtering_script.sh'.format(mapq_filtering_dir)
if not os.path.exists(mapq_filtering_dir + '/logs'):
    os.makedirs(mapq_filtering_dir + '/logs')
mapq_filtering_script = open(mapq_filtering_script_name, 'w')

mapq_filtering_script.write('#!/bin/bash\n\n')
mapq_filtering_script.write('set -euxo pipefail\n\n')

if start_from_bam == False:

    # Writing bash commands to filter out all aligned reads with score less than user-determined MAPQ score from the .bam files
    for list_counter in range(len(chip_name)):
        mapq_filtering_script.write('samtools view {} -@ {} -h -b {}/{}.mapped.bam > {}/{}.mapqfiltered.bam 2> {}/logs/{}.mapqfiltering.err\n\n'.format(
            samtools_view_arg, 
            cpu_count, 
            bwa_mem_aligning_dir, 
            chip_name[list_counter], 
            mapq_filtering_dir, 
            chip_name[list_counter], 
            mapq_filtering_dir, 
            chip_name[list_counter]))

    if args.deltemp:        
        # Immediately deletes the aligned bam after MAPQ based filtering
        for list_counter in range(len(chip_name)):
            mapq_filtering_script.write('rm -r -f {}/{}.mapped.bam &\n'.format(
                bwa_mem_aligning_dir, 
                chip_name[list_counter]))

            if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                mapq_filtering_script.write('\nwait\n\n')

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        for list_counter in range(len(ctrl_name)):
            mapq_filtering_script.write('samtools view {} -@ {} -h -b {}/{}.mapped.bam > {}/{}.mapqfiltered.bam 2> {}/logs/{}.mapqfiltering.err\n\n'.format(
                samtools_view_arg, 
                cpu_count,
                bwa_mem_aligning_dir, 
                ctrl_name[list_counter], 
                mapq_filtering_dir, 
                ctrl_name[list_counter], 
                mapq_filtering_dir, 
                ctrl_name[list_counter]))

        if args.deltemp:        
            for list_counter in range(len(ctrl_name)):
                mapq_filtering_script.write('rm -r -f {}/{}.mapped.bam &\n'.format(
                    bwa_mem_aligning_dir, 
                    ctrl_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                    mapq_filtering_script.write('\nwait\n\n')

mapq_filtering_script.close() # Closing the script '07_MAPQ_filtering_script.sh'. Flushing the write buffer



####################################################################################################################
### 08_results_script.sh
####################################################################################################################

# Create a directory "08_results" for the final aligned reads, complete with their merged versions,
#   their indices, and bigwig coverage files for track viewing 
# Create a script "08_results_script.sh" within, that calls for sequence reads merging: samtools merge, 
#   sequence reads sorting: samtools sort, sequence reads indexing: samtools index, fingerprint plot generation:
#   bedtools plotFingerprint, and bigwig file generation out of the aligned sequence reads: bedtools bamCoverage
# In: filename.mapqfiltered.bam (output from 07_MAPQ_filtering_script.sh)
# Out: filename.bam (sorted input for 11_macs2_peak_calling_script.sh, 12_gem_peak_calling_script.sh, 
#   13_homer_peak_calling_script.sh, 22_peaks_processing_script.sh)
# Out: filename.bam.bai (Indices for the bam file. Required, but not to be further processed downstreams)
# Out: fingerprint_plots/filename.png (Fingerprint plots of the aligned reads. Not to be further processed downstreams)
# Out: filename.bw (BigWig files for peak visualization purposes. Not to be further processed downstreams)
# Out: filename.namesorted.bam (input for 14_genrich_peak_calling_script.sh)

results_dir = '{}/08_results'.format(output_dir)
results_script_name = '{}/08_results_script.sh'.format(results_dir)
if not os.path.exists(results_dir + '/logs'):
    os.makedirs(results_dir + '/logs')
results_script = open(results_script_name, 'w')

if not os.path.exists(results_dir + '/fingerprint_plots'):
    os.makedirs(results_dir + '/fingerprint_plots')

results_script.write('#!/bin/bash\n\n')
results_script.write('set -euxo pipefail\n\n')



if start_from_bam == False:

    # Writing bash commands to sort each of the aligned sequence read files
    for list_counter in range(len(chip_name)):
        results_script.write('samtools sort -@ {} {}/{}.mapqfiltered.bam > {}/{}.bam\n\n'.format(
            cpu_count, 
            mapq_filtering_dir, 
            chip_name[list_counter], 
            results_dir, 
            chip_name[list_counter]))

    if args.deltemp:        
        # Immediately deletes the mapq filtered bam after sorting
        for list_counter in range(len(chip_name)):
            results_script.write('rm -r -f {}/{}.mapqfiltered.bam &\n'.format(
                mapq_filtering_dir, 
                chip_name[list_counter]))

            if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                results_script.write('\nwait\n\n')
                
    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        for list_counter in range(len(ctrl_name)):
            results_script.write('samtools sort -@ {} {}/{}.mapqfiltered.bam > {}/{}.bam\n\n'.format(
                cpu_count, 
                mapq_filtering_dir, 
                ctrl_name[list_counter], 
                results_dir, 
                ctrl_name[list_counter]))

        if args.deltemp:        
            for list_counter in range(len(ctrl_name)):
                results_script.write('rm -r -f {}/{}.mapqfiltered.bam &\n'.format(
                    mapq_filtering_dir, 
                    ctrl_name[list_counter]))

                if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
                    results_script.write('\nwait\n\n')



if start_from_bam == True:

    # Writing bash commands to sort each of the aligned sequence read files
    for list_counter in range(len(chip_name)):
        # Only do this if the input files and paths are not the output files and paths of this command
        if chip_r1_absolute_path[list_counter] != '{}/{}.bam'.format(results_dir, chip_name[list_counter]):
            results_script.write('samtools sort -@ {} {} > {}/{}.bam\n\n'.format(
                cpu_count, 
                chip_r1_absolute_path[list_counter], 
                results_dir, 
                chip_name[list_counter]))

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        # Writing bash commands to sort each of the aligned sequence read files
        for list_counter in range(len(ctrl_name)):
            # Only do this if the input files and paths are not the output files and paths of this command
            if ctrl_r1_absolute_path[list_counter] != '{}/{}.bam'.format(results_dir, ctrl_name[list_counter]):
                results_script.write('samtools sort -@ {} {} > {}/{}.bam\n\n'.format(
                    cpu_count, 
                    ctrl_r1_absolute_path[list_counter], 
                    results_dir, 
                    ctrl_name[list_counter]))



# Writing bash commands to generate an index file for all individual .bam files
for list_counter in range(len(chip_name)):
    results_script.write('samtools index -@ {} {}/{}.bam &\n'.format(
        cpu_count, 
        results_dir, 
        chip_name[list_counter]))
    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
        results_script.write('\nwait\n\n')

if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    for list_counter in range(len(ctrl_name)):
        results_script.write('samtools index -@ {} {}/{}.bam &\n'.format(
            cpu_count, 
            results_dir, 
            ctrl_name[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
            results_script.write('\nwait\n\n')



bam_chip_list = ['{}/{}.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]
bam_chip_string = ' --bam ' + ' '.join(bam_chip_list)

if norm_ref != 'none':
    norm_bam_chip_list = ['{}/{}.normmapped.bam'.format(bwa_mem_aligning_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]
    norm_bam_chip_string = ' --norm_bam ' + ' '.join(norm_bam_chip_list)

else:
    norm_bam_chip_string = ' '

if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    bam_ctrl_list = ['{}/{}.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]
    bam_ctrl_string = ' ' + ' '.join(bam_ctrl_list)

    if norm_ref != 'none':
        norm_bam_ctrl_list = ['{}/{}.normmapped.bam'.format(bwa_mem_aligning_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]
        norm_bam_ctrl_string = ' ' + ' '.join(norm_bam_ctrl_list)

    else:
        norm_bam_ctrl_string = ' '

else:
    bam_ctrl_string = ' '
    norm_bam_ctrl_string = ' '

# Writing bash commands to generate bedgraph and bigwig coverage files out of all individual .bam files
results_script.write('caras_reads_normalizer.py --thread {} {} --make_bam --make_bdg --make_bw {}{}{}{} --log_dir {}/logs\n\n'.format(
    cpu_count,
    reads_normalizer_arg,
    bam_chip_string,
    bam_ctrl_string,
    norm_bam_chip_string,
    norm_bam_ctrl_string,
    results_dir))



# # Writing bash commands to sort each of the aligned sequence read files
# for list_counter in range(len(chip_name)):
#     results_script.write('samtools sort -@ {} {}/{}.normalized.bam > {}/{}.normalized.sorted.bam\n\n'.format(
#         cpu_count, 
#         mapq_filtering_dir, 
#         chip_name[list_counter], 
#         results_dir, 
#         chip_name[list_counter]))

# if args.deltemp:        
#     # Immediately deletes the mapq filtered bam after sorting
#     for list_counter in range(len(chip_name)):
#         results_script.write('rm -r -f {}/{}.normalized.bam &\n'.format(
#             mapq_filtering_dir, 
#             chip_name[list_counter]))

#         if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
#             results_script.write('\nwait\n\n')
            
# if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
#     for list_counter in range(len(ctrl_name)):
#         results_script.write('samtools sort -@ {} {}/{}.normalized.bam > {}/{}.normalized.sorted.bam\n\n'.format(
#             cpu_count, 
#             mapq_filtering_dir, 
#             ctrl_name[list_counter], 
#             results_dir, 
#             ctrl_name[list_counter]))

#     if args.deltemp:        
#         for list_counter in range(len(ctrl_name)):
#             results_script.write('rm -r -f {}/{}.normalized.bam &\n'.format(
#                 mapq_filtering_dir, 
#                 ctrl_name[list_counter]))

#             if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
#                 results_script.write('\nwait\n\n')



# Writing bash commands to generate an index file for all individual .bam files
for list_counter in range(len(chip_name)):
    results_script.write('samtools index -@ {} {}/{}.normalized.bam &\n'.format(
        cpu_count, 
        results_dir, 
        chip_name[list_counter]))
    if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
        results_script.write('\nwait\n\n')

if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    for list_counter in range(len(ctrl_name)):
        results_script.write('samtools index -@ {} {}/{}.normalized.bam &\n'.format(
            cpu_count, 
            results_dir, 
            ctrl_name[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
            results_script.write('\nwait\n\n')



plotfingerprint_chip_list = ['{}/{}.normalized.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]
plotfingerprint_chip_label_list = [chip_name[list_counter] for list_counter in range(len(chip_name))]

if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    plotfingerprint_ctrl_list = ['{}/{}.normalized.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]
    plotfingerprint_ctrl_label_list = [ctrl_name[list_counter] for list_counter in range(len(ctrl_name))]

else:
    plotfingerprint_ctrl_list = []
    plotfingerprint_ctrl_label_list = []


plotfingerprint_chip_nested_list = []
plotfingerprint_chip_label_nested_list = []

plotfingerprint_ctrl_nested_list = []
plotfingerprint_ctrl_label_nested_list = []


for counter in range(math.ceil(len(plotfingerprint_chip_list) / 10)): 
    try:
        plotfingerprint_chip_nested_list.append(plotfingerprint_chip_list[(counter * 10) : ((counter + 1) * 10)]) 

    except:
        plotfingerprint_chip_nested_list.append(plotfingerprint_chip_list[(counter + 10) : ])

    try:
        plotfingerprint_chip_label_nested_list.append(plotfingerprint_chip_label_list[(counter * 10) : ((counter + 1) * 10)]) 
    
    except:
        plotfingerprint_chip_label_nested_list.append(plotfingerprint_chip_label_list[(counter + 10) : ])


for counter in range(math.ceil(len(plotfingerprint_ctrl_list) / 10)): 
    try:
        plotfingerprint_ctrl_nested_list.append(plotfingerprint_ctrl_list[(counter * 10) : ((counter + 1) * 10)])

    except:
        plotfingerprint_ctrl_nested_list.append(plotfingerprint_ctrl_list[(counter + 10) : ])

    try:
        plotfingerprint_ctrl_label_nested_list.append(plotfingerprint_ctrl_label_list[(counter * 10) : ((counter + 1) * 10)])

    except:
        plotfingerprint_ctrl_label_nested_list.append(plotfingerprint_ctrl_label_list[(counter + 10) : ])


while len(plotfingerprint_chip_nested_list) > len(plotfingerprint_ctrl_nested_list):
    plotfingerprint_ctrl_nested_list.append([])

while len(plotfingerprint_chip_nested_list) < len(plotfingerprint_ctrl_nested_list):
    plotfingerprint_chip_nested_list.append([])


while len(plotfingerprint_chip_label_nested_list) > len(plotfingerprint_ctrl_label_nested_list):
    plotfingerprint_ctrl_label_nested_list.append([])

while len(plotfingerprint_chip_label_nested_list) < len(plotfingerprint_ctrl_label_nested_list):
    plotfingerprint_chip_label_nested_list.append([])


plotfingerprint_combined_string_list = []
plotfingerprint_combined_label_string_list = []


for nested_list_counter in range(len(plotfingerprint_chip_nested_list)):
    plotfingerprint_combined_list = plotfingerprint_chip_nested_list[nested_list_counter] + plotfingerprint_ctrl_nested_list[nested_list_counter]
    plotfingerprint_combined_string = ' '.join(plotfingerprint_combined_list)
    plotfingerprint_combined_string_list.append(plotfingerprint_combined_string)

for label_nested_list_counter in range(len(plotfingerprint_chip_label_nested_list)):
    plotfingerprint_combined_label_list = plotfingerprint_chip_label_nested_list[label_nested_list_counter] + plotfingerprint_ctrl_label_nested_list[label_nested_list_counter]
    plotfingerprint_combined_label_string = ' '.join(plotfingerprint_combined_label_list)
    plotfingerprint_combined_label_string_list.append(plotfingerprint_combined_label_string)


for nested_list_counter in range(len(plotfingerprint_chip_nested_list)):
    # Writing bash commands to generate a fingerprint plot .png file for all individual .bam files
    results_script.write('plotFingerprint {} -p {} -b {} -l {} -o {}/fingerprint_plots/{}_rep{}-{}.png &\n\n'.format(
        plotfingerprint_arg,
        cpu_count,
        plotfingerprint_combined_string_list[nested_list_counter],
        plotfingerprint_combined_label_string_list[nested_list_counter],
        results_dir,
        dataset_name,
        (nested_list_counter * 10) + 1,
        (nested_list_counter + 1) * 10))

    # Writing bash commands to generate a fingerprint plot .svg file for all individual .bam files
    results_script.write('plotFingerprint {} -p {} -b {} -l {} -o {}/fingerprint_plots/{}_rep{}-{}.svg &\n\n'.format(
        plotfingerprint_arg,
        cpu_count,
        plotfingerprint_combined_string_list[nested_list_counter],
        plotfingerprint_combined_label_string_list[nested_list_counter],
        results_dir,
        dataset_name,
        (nested_list_counter * 10) + 1,
        (nested_list_counter + 1) * 10))

    results_script.write('\nwait\n\n')



if analysis_mode == 'bulk':

    # Converting both ChIP and control .bam lists from vertical lists into space-separated 
    #   serial string (to be used as an input argument in samtools merge)
    samtools_merge_chip_list = ['{}/{}.normalized.bam'.format(results_dir, 
        chip_name[list_counter]) for list_counter in range(len(chip_name))]
    samtools_merge_chip_string = ' '.join(samtools_merge_chip_list)

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        samtools_merge_ctrl_list = ['{}/{}.normalized.bam'.format(results_dir, 
            ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]
        samtools_merge_ctrl_string = ' '.join(samtools_merge_ctrl_list)


    # Writing bash commands to merge all ChIP and all control .bam files
    # The replicate-merged .bam files are going to be used for fold change calculation downstream
    results_script.write('samtools merge -@ {} -f {}/{}_chip_merged.normalized.bam {} &\n'.format(
        cpu_count, 
        results_dir, 
        dataset_name, 
        samtools_merge_chip_string))

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        results_script.write('samtools merge -@ {} -f {}/{}_ctrl_merged.normalized.bam {} &\n'.format(
            cpu_count, 
            results_dir, 
            dataset_name, 
            samtools_merge_ctrl_string))

    results_script.write('\nwait\n\n')


    # Writing bash commands to generate an index file for all replicate-merged .bam files
    results_script.write('samtools index -@ {} {}/{}_chip_merged.normalized.bam &\n'.format(
        cpu_count, 
        results_dir, 
        dataset_name))

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        results_script.write('samtools index -@ {} {}/{}_ctrl_merged.normalized.bam &\n'.format(
            cpu_count, 
            results_dir, 
            dataset_name))

    results_script.write('\nwait\n\n')



    merged_bam_chip_string = ' {}/{}_chip_merged.normalized.bam'.format(results_dir, dataset_name)

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        merged_bam_ctrl_string = ' {}/{}_ctrl_merged.normalized.bam'.format(results_dir, dataset_name)
    else:
        merged_bam_ctrl_string = ' '

    # Writing bash commands to generate bedgraph and bigwig coverage files out of all individual .bam files
    results_script.write('caras_reads_normalizer.py --thread {} {} --make_bdg --make_bw --bam{}{} --log_dir {}/logs\n\n'.format(
        cpu_count,
        reads_normalizer_arg,
        merged_bam_chip_string,
        merged_bam_ctrl_string,
        results_dir))



    # Writing bash commands to generate a fingerprint plot .png file for all replicate-merged .bam files
    results_script.write('plotFingerprint {} -p {} -b {}/{}_chip_merged.normalized.bam -l {}_chip_merged -o {}/fingerprint_plots/{}_merged.png &\n'.format(
        plotfingerprint_arg,
        cpu_count,
        results_dir,
        dataset_name,
        dataset_name,
        results_dir,
        dataset_name))

    # Writing bash commands to generate a fingerprint plot .svg file for all replicate-merged .bam files
    results_script.write('plotFingerprint {} -p {} -b {}/{}_chip_merged.normalized.bam -l {}_chip_merged -o {}/fingerprint_plots/{}_merged.svg &\n'.format(
        plotfingerprint_arg,
        cpu_count,
        results_dir,
        dataset_name,
        dataset_name,
        results_dir,
        dataset_name))

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        # Writing bash commands to generate a fingerprint plot .png file for all replicate-merged .bam files
        results_script.write('plotFingerprint {} -p {} -b {}/{}_chip_merged.normalized.bam {}/{}_ctrl_merged.normalized.bam -l {}_chip_merged {}_ctrl_merged -o {}/fingerprint_plots/{}_merged.png &\n'.format(
            plotfingerprint_arg,
            cpu_count,
            results_dir,
            dataset_name,
            results_dir,
            dataset_name,
            dataset_name,
            dataset_name,
            results_dir,
            dataset_name))

        # Writing bash commands to generate a fingerprint plot .svg file for all replicate-merged .bam files
        results_script.write('plotFingerprint {} -p {} -b {}/{}_chip_merged.normalized.bam {}/{}_ctrl_merged.normalized.bam -l {}_chip_merged {}_ctrl_merged -o {}/fingerprint_plots/{}_merged.svg &\n'.format(
            plotfingerprint_arg,
            cpu_count,
            results_dir,
            dataset_name,
            results_dir,
            dataset_name,
            dataset_name,
            dataset_name,
            results_dir,
            dataset_name))
            
    results_script.write('\nwait\n\n')



# Writing bash commands to sort all individual .bam files according to the reads name
# The name-sorted individual .bam files are necessary as inputs for Genrich peak calling
for list_counter in range(len(chip_name)):
    results_script.write('samtools sort -@ {} -n {}/{}.normalized.bam > {}/{}.namesorted.bam\n\n'.format(
        cpu_count, 
        results_dir, 
        chip_name[list_counter], 
        results_dir, 
        chip_name[list_counter]))

if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    for list_counter in range(len(ctrl_name)):
        results_script.write('samtools sort -@ {} -n {}/{}.normalized.bam > {}/{}.namesorted.bam\n\n'.format(
            cpu_count, 
            results_dir, 
            ctrl_name[list_counter], 
            results_dir, 
            ctrl_name[list_counter]))

results_script.close() # Closing the script '08_results_script.sh'. Flushing the write buffer



####################################################################################################################
### 09_aligned_reads_quality_control_script.sh
####################################################################################################################

# Create a directory "09_aligned_reads_quality_control" for the quality control 
#   report files of the aligned and filtered sequencing reads
# Create a script "09_aligned_reads_quality_control_script.sh" within, 
#   that calls for sequencing reads QC analysis program: fastqc
# In: filename.bam (output from 08_results_script.sh)
# Out: filename_fastqc.html (Reports. not to be further processed downstreams)

aligned_reads_quality_control_dir = '{}/09_aligned_reads_quality_control'.format(output_dir)
aligned_reads_quality_control_script_name = '{}/09_aligned_reads_quality_control_script.sh'.format(aligned_reads_quality_control_dir)
if not os.path.exists(aligned_reads_quality_control_dir):
    os.makedirs(aligned_reads_quality_control_dir)
aligned_reads_quality_control_script = open(aligned_reads_quality_control_script_name, 'w')

aligned_reads_quality_control_script.write('#!/bin/bash\n\n')
aligned_reads_quality_control_script.write('set -euxo pipefail\n\n')


# Bash commands to call fastqc to run quality control analysis on all individual .bam files 
for list_counter in range(len(chip_name)):
    aligned_reads_quality_control_script.write('fastqc {} -t {} {}/{}.normalized.bam -o {}\n\n'.format(
        fastqc3_arg, 
        cpu_count, 
        results_dir, 
        chip_name[list_counter], 
        aligned_reads_quality_control_dir))

if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    for list_counter in range(len(ctrl_name)):
        aligned_reads_quality_control_script.write('fastqc {} -t {} {}/{}.normalized.bam -o {}\n\n'.format(
            fastqc3_arg, 
            cpu_count, 
            results_dir, 
            ctrl_name[list_counter], 
            aligned_reads_quality_control_dir))


if analysis_mode == 'bulk':

    # Bash commands to call fastqc to run quality control analysis on all the replicate-merged .bam files
    aligned_reads_quality_control_script.write('fastqc {} -t {} {}/{}_chip_merged.normalized.bam -o {}\n\n'.format(
        fastqc3_arg, 
        cpu_count, 
        results_dir, 
        dataset_name, 
        aligned_reads_quality_control_dir))

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        aligned_reads_quality_control_script.write('fastqc {} -t {} {}/{}_ctrl_merged.normalized.bam -o {}\n\n'.format(
            fastqc3_arg, 
            cpu_count, 
            results_dir, 
            dataset_name, 
            aligned_reads_quality_control_dir))

aligned_reads_quality_control_script.close() # Closing the script '09_aligned_reads_quality_control_script.sh'. Flushing the write buffer



####################################################################################################################
### 11_macs2_peak_calling_script.sh
####################################################################################################################

# Create a directory "11_macs2_peak_calling" for the MACS2-called peak list file
# Create a script "11_macs2_peak_calling_script.sh" within, that calls for ChIP-seq peak calling program: MACS2
# Various output will be generated. The peak list file that is going to be used 
#   in subsequent analysis is the one in narrowPeak format, .narrowPeak extension
# In: filename.bam (output from 08_results_script.sh)
# Out: setname_MACS2_peaks.narrowPeak

macs2_dir = '{}/11_macs2_peak_calling'.format(output_dir)
macs2_peak_calling_script_name = '{}/11_macs2_peak_calling_script.sh'.format(macs2_dir)
if not os.path.exists(macs2_dir + '/logs'):
    os.makedirs(macs2_dir + '/logs')
macs2_peak_calling_script = open(macs2_peak_calling_script_name, 'w')

macs2_peak_calling_script.write('#!/bin/bash\n\n')
macs2_peak_calling_script.write('set -euxo pipefail\n\n')



# The -f BAM flag and argument is used to process single-end reads. 
#   It tells MACS that the input .bam files consists of one directional reads. 
if read_mode == 'single':
    macs2_read_mode_arg = ' -f BAM'

# The -f BAMPE flag and argument is used to process paired-end reads. 
#   It tells MACS that the input .bam files consists of pairs of two directional reads. 
elif read_mode == 'paired':
    macs2_read_mode_arg = ' -f BAMPE'

# Default mode (no flag) is chosen here, which is optimized for transcription factor peak calling.
# Only when CaRAS is running in narrow peak mode (ChIP protein is a transcription factor)
if peak_type == 'narrow':
    macs2_peak_type_arg = ''

# Broad peak mode (--broad flag) is chosen here, which is optimized for histone modifier peak calling. 
# Only when CaRAS is running in broad peak mode (ChIP protein is a histone modifier)
elif peak_type == 'broad':
    macs2_peak_type_arg = ' --broad'



if analysis_mode == 'bulk':

    macs2_chip_list = ['{}/{}.normalized.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]
    macs2_chip_string = ' -t ' + ' '.join(macs2_chip_list)

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        macs2_ctrl_list = ['{}/{}.normalized.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]
        macs2_ctrl_string = ' -c ' + ' '.join(macs2_ctrl_list)

    else:
        macs2_ctrl_string = ''

    macs2_peak_calling_script.write('macs2 callpeak{}{}{}{}{} --keep-dup all -g {} --name {}_MACS2 --outdir {} 1> {}/logs/{}.MACS2.out 2> {}/logs/{}.MACS2.err\n'.format(
        macs2_callpeak_arg,
        macs2_read_mode_arg,
        macs2_peak_type_arg,
        macs2_chip_string,
        macs2_ctrl_string, 
        effective_genome_size, 
        dataset_name,
        macs2_dir, 
        macs2_dir, 
        dataset_name,
        macs2_dir, 
        dataset_name))



if analysis_mode == 'single_cell':

    macs2_chip_list = [' -t {}/{}.normalized.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        macs2_ctrl_list = [' -c {}/{}.normalized.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]

    else:
        macs2_ctrl_list = [' ' for list_counter in range(len(chip_name))]

    for list_counter in range(len(chip_name)):
        macs2_peak_calling_script.write('macs2 callpeak{}{}{}{}{} --keep-dup all -g {} --name {}_rep{}_MACS2 --outdir {} 1> {}/logs/{}_rep{}.MACS2.out 2> {}/logs/{}_rep{}.MACS2.err &\n'.format(
            macs2_callpeak_arg,
            macs2_read_mode_arg,
            macs2_peak_type_arg,
            macs2_chip_list[list_counter],
            macs2_ctrl_list[list_counter], 
            effective_genome_size, 
            dataset_name,
            (list_counter + 1),
            macs2_dir, 
            macs2_dir, 
            dataset_name,
            (list_counter + 1),
            macs2_dir, 
            dataset_name,
            (list_counter + 1)))

    macs2_peak_calling_script.write('\nwait\n\n')

macs2_peak_calling_script.close() # Closing the script '11_macs2_peak_calling_script.sh'. Flushing the write buffer



####################################################################################################################
### 12_gem_peak_calling_script.sh
####################################################################################################################

# Create a directory "12_gem_peak_calling" for the GEM-called peak list file
# Create a script "12_gem_peak_calling_script.sh" within, that calls for ChIP-seq peak calling program: GEM
# There are 4 choices of read distribution that can be used in GEM peak calling. 
#   Default read distribution, which is used for general ChIP-seq datasets, is chosen here.
# Various output will be generated. The peak list file that is going to be used 
#   in subsequent analysis is the one named GPS_events
# In: filename.normalized.bam (output from 08_results_script.sh)
# Out: setname_GEM_GPS_events.txt (input for 21_peaks_merging_script.sh)

gem_dir = '{}/12_gem_peak_calling'.format(output_dir)
gem_peak_calling_script_name = '{}/12_gem_peak_calling_script.sh'.format(gem_dir)
if not os.path.exists(gem_dir + '/logs'):
    os.makedirs(gem_dir + '/logs')

gem_peak_calling_script = open(gem_peak_calling_script_name, 'w')

gem_peak_calling_script.write('#!/bin/bash\n\n')
gem_peak_calling_script.write('set -euxo pipefail\n\n')


# The -Xmx100G flag is necessary when processing paired-end reads with multiple 
#   replicates to give GEM enough memory space to work
# GEM peak caller utilizes different background signal for peak calling, depending on the experiment type. 
# Read_Distribution_default.txt is used for general ChIP-seq experiment. 
#   It has specialized read distributions for ChIP-exo and CLIP experiment.
# GEM peak caller harnesses the knowledge of the genomic sequence, which is why reference genome is required.
# GEM peak caller is calling peaks while considering all the possible binding motifs, 
#   thus it is required to set the range of motif length to be considered (k_min to k_max)



if analysis_mode == 'bulk':

    gem_chip_list = ['--expt {}/{}.normalized.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]
    gem_chip_string = ' '.join(gem_chip_list)

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        gem_ctrl_list = ['--ctrl {}/{}.normalized.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]
        gem_ctrl_string = ' '.join(gem_ctrl_list)

    else:
        gem_ctrl_string = ''

    gem_peak_calling_script.write('java -jar {} {} --nrf --t {} --d {}/GEM/Read_Distribution_default.txt --g {}/GEM/{}.chrom.sizes --genome {}/GEM/{}_Chr_FASTA --s {} {} {} --f SAM --out {}/{}_GEM 1> {}/logs/{}.GEM.out 2> {}/logs/{}.GEM.err\n\n'.format(
        gem_full_path, 
        gem_arg, 
        math.ceil(cpu_count/2), 
        genome_dir, 
        genome_dir, 
        genome_ref, 
        genome_dir, 
        genome_ref, 
        effective_genome_size,
        gem_chip_string, 
        gem_ctrl_string, 
        gem_dir, 
        dataset_name, 
        gem_dir, 
        dataset_name, 
        gem_dir, 
        dataset_name))


    gem_peak_calling_script.write('cp -r -f {}/{}_GEM/* {}\n\n'.format(
        gem_dir, 
        dataset_name, 
        gem_dir))

    gem_peak_calling_script.write('rm -r -f {}/{}_GEM\n\n'.format(
        gem_dir, 
        dataset_name))



if analysis_mode == 'single_cell':

    gem_chip_list = [' --expt {}/{}.normalized.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        gem_ctrl_list = [' --ctrl {}/{}.normalized.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]

    else:
        gem_ctrl_list = [' ' for list_counter in range(len(chip_name))]
        

    for list_counter in range(len(chip_name)):
        gem_peak_calling_script.write('java -jar {} {} --nrf --t {} --d {}/GEM/Read_Distribution_default.txt --g {}/GEM/{}.chrom.sizes --genome {}/GEM/{}_Chr_FASTA --s {}{}{} --f SAM --out {}/{}_rep{}_GEM 1> {}/logs/{}_rep{}.GEM.out 2> {}/logs/{}_rep{}.GEM.err\n\n'.format(
            gem_full_path,
            gem_arg,
            math.ceil(cpu_count/2),
            genome_dir,
            genome_dir,
            genome_ref,
            genome_dir,
            genome_ref,
            effective_genome_size,
            gem_chip_list[list_counter],
            gem_ctrl_list[list_counter],
            gem_dir,
            dataset_name,
            (list_counter + 1),
            gem_dir,
            dataset_name,
            (list_counter + 1),
            gem_dir,
            dataset_name,
            (list_counter + 1)))


        gem_peak_calling_script.write('cp -r -f {}/{}_rep{}_GEM/* {}\n\n'.format(
            gem_dir, 
            dataset_name, 
            (list_counter + 1),
            gem_dir))

        gem_peak_calling_script.write('rm -r -f {}/{}_rep{}_GEM\n\n'.format(
            gem_dir, 
            dataset_name,
            (list_counter + 1)))


gem_peak_calling_script.close() # Closing the script '12_gem_peak_calling_script.sh'. Flushing the write buffer



####################################################################################################################
### 13_homer_peak_calling_script.sh
####################################################################################################################

# Create a directory "13_homer_peak_calling" for the HOMER-called peak list file
# Create a script "13_homer_peak_calling_script.sh" within, that calls for ChIP-seq peak calling program: HOMER
# Various output will be generated. The peak list file that is going to be used 
#   in subsequent analysis is the one in HOMER format, with user-determined filename and extension
# In: filename.normalized.bam (output from 08_results_script.sh)
# Out: setname_HOMER.peaks (input for 21_peaks_merging_script.sh)

homer_dir = '{}/13_homer_peak_calling'.format(output_dir)
homer_peak_calling_script_name = '{}/13_homer_peak_calling_script.sh'.format(homer_dir)
if not os.path.exists(homer_dir + '/logs'):
    os.makedirs(homer_dir + '/logs')
homer_peak_calling_script = open(homer_peak_calling_script_name, 'w')

homer_peak_calling_script.write('#!/bin/bash\n\n')
homer_peak_calling_script.write('set -euxo pipefail\n\n')



# No special flag needed for tags generation of single end sequencing reads
if read_mode == 'single':
    homer_read_mode_arg = ' '

# The -sspe flag used here is necessary to process paired-end reads. 
# It will basically flip the sequence of the second reads to match the first reads. 
elif read_mode == 'paired':
    homer_read_mode_arg = ' -sspe'

# -style factor flag and argument are used here which makes HOMER optimized for transcription factor peak calling
# Only when CaRAS is running in narrow peak mode (ChIP protein is a transcription factor)
if peak_type == 'narrow':
    homer_peak_type_arg = ' -style factor'

# -style histone flag and argument are used here which makes HOMER optimized for histone modifier peak calling
# Only when CaRAS is running in broad peak mode (ChIP protein is a histone modifier)
elif peak_type == 'broad':
    homer_peak_type_arg = ' -style histone'



if analysis_mode == 'bulk':

    homer_peak_calling_script.write('makeTagDirectory {}/chip_tag {}/{}_chip_merged.normalized.bam -format sam{} -keepAll\n\n'.format(
    homer_dir,
    results_dir, 
    dataset_name,
    homer_read_mode_arg))

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        homer_peak_calling_script.write('makeTagDirectory {}/ctrl_tag {}/{}_ctrl_merged.normalized.bam -format sam{} -keepAll\n\n'.format(
        homer_dir,
        results_dir, 
        dataset_name,
        homer_read_mode_arg))
        

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        homer_peak_calling_script.write('findPeaks {}/chip_tag {}{} -tbp 0 -inputtbp 0 -gsize {} -o {}/{}_HOMER.peaks -i {}/ctrl_tag 1> {}/logs/{}.HOMER.out 2> {}/logs/{}.HOMER.err\n\n'.format(
            homer_dir, 
            homer_findPeaks_arg, 
            homer_peak_type_arg,
            effective_genome_size,
            homer_dir, 
            dataset_name, 
            homer_dir, 
            homer_dir, 
            dataset_name,
            homer_dir, 
            dataset_name))

    else:
        homer_peak_calling_script.write('findPeaks {}/chip_tag {}{} -tbp 0 -inputtbp 0 -gsize {} -o {}/{}_HOMER.peaks 1> {}/logs/{}.HOMER.out 2> {}/logs/{}.HOMER.err\n\n'.format(
            homer_dir, 
            homer_findPeaks_arg, 
            homer_peak_type_arg,
            effective_genome_size,
            homer_dir, 
            dataset_name,
            homer_dir, 
            dataset_name,
            homer_dir, 
            dataset_name))



if analysis_mode == 'single_cell':

    homer_chip_list = [' {}/{}.normalized.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        homer_ctrl_list = [' {}/{}.normalized.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]

    else:
        homer_ctrl_list = [' ' for list_counter in range(len(chip_name))]


    for list_counter in range(len(chip_name)):
        homer_peak_calling_script.write('makeTagDirectory {}/{}_tag{} -format sam{} -keepAll \n\n'.format(
            homer_dir, 
            chip_name[list_counter],
            homer_chip_list[list_counter],
            homer_read_mode_arg))

        # if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
        #     homer_peak_calling_script.write('\nwait\n\n')

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        for list_counter in range(len(chip_name)):
            homer_peak_calling_script.write('makeTagDirectory {}/{}_tag{} -format sam{} -keepAll \n\n'.format(
                homer_dir, 
                ctrl_name[list_counter],
                homer_ctrl_list[list_counter],
                homer_read_mode_arg))

            # if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(ctrl_name) - 1):
            #     homer_peak_calling_script.write('\nwait\n\n')


    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        for list_counter in range(len(chip_name)):
            homer_peak_calling_script.write('findPeaks {}/{}_tag {}{} -tbp 0 -inputtbp 0 -gsize {} -o {}/{}_rep{}_HOMER.peaks -i {}/{}_tag 1> {}/logs/{}_rep{}.HOMER.out 2> {}/logs/{}_rep{}.HOMER.err &\n'.format(
                homer_dir, 
                chip_name[list_counter],
                homer_findPeaks_arg, 
                homer_peak_type_arg,
                effective_genome_size,
                homer_dir, 
                dataset_name, 
                (list_counter + 1),
                homer_dir, 
                ctrl_name[list_counter],
                homer_dir, 
                dataset_name,
                (list_counter + 1),
                homer_dir, 
                dataset_name,
                (list_counter + 1)))

            if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                homer_peak_calling_script.write('\nwait\n\n')
                
    else:
        for list_counter in range(len(chip_name)):
            homer_peak_calling_script.write('findPeaks {}/{}_tag {}{} -tbp 0 -inputtbp 0 -gsize {} -o {}/{}_rep{}_HOMER.peaks 1> {}/logs/{}_rep{}.HOMER.out 2> {}/logs/{}_rep{}.HOMER.err &\n'.format(
                homer_dir, 
                chip_name[list_counter],
                homer_findPeaks_arg, 
                homer_peak_type_arg,
                effective_genome_size,
                homer_dir, 
                dataset_name,
                (list_counter + 1),
                homer_dir, 
                dataset_name,
                (list_counter + 1),
                homer_dir, 
                dataset_name,
                (list_counter + 1)))

            if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
                homer_peak_calling_script.write('\nwait\n\n')
            

homer_peak_calling_script.write('rm -r -f {}/*_tag\n\n'.format(homer_dir))
homer_peak_calling_script.write('rm -r -f {}/*_tag\n\n'.format(homer_dir))

homer_peak_calling_script.close() # Closing the script '13_homer_peak_calling_script.sh'. Flushing the write buffer



####################################################################################################################
### 14_genrich_peak_calling_script.sh
####################################################################################################################

# Create a directory "14_genrich_peak_calling" for the Genrich-called peak list file
# Create a script "14_genrich_peak_calling_script.sh" within, that calls for ChIP-seq peak calling program: Genrich
# Default mode is chosen here, which is optimized for peak calling in general ChIP-seq datasets
# Various output will be generated. The peak list file that is going to be used 
#   in subsequent analysis is the one in narrowPeak format, with .narrowPeak extension
# In: filename.namesorted.bam (output from 08_results_script.sh)
# Out: setname_Genrich.narrowPeak (input for 21_peaks_merging_script.sh)



genrich_dir = '{}/14_genrich_peak_calling'.format(output_dir)
genrich_peak_calling_script_name = '{}/14_genrich_peak_calling_script.sh'.format(genrich_dir)
if not os.path.exists(genrich_dir + '/logs'):
    os.makedirs(genrich_dir + '/logs')
genrich_peak_calling_script = open(genrich_peak_calling_script_name, 'w')

genrich_peak_calling_script.write('#!/bin/bash\n\n')
genrich_peak_calling_script.write('set -euxo pipefail\n\n')



# Genrich is processing paired-end reads by default: it is designed to discard unpaired reads.
# The -y flag is used in single mode so that Genrich will keep them.
if read_mode == 'single':
    genrich_read_mode_arg = ' --mode single -y'

# Genrich is processing paired-end reads by default: it is designed to discard unpaired reads.
# No special flag needed for tags generation of paired end sequencing reads 
elif read_mode == 'paired':
    genrich_read_mode_arg = ' --mode paired'



if analysis_mode == 'bulk':

    genrich_chip_list = ['{}/{}.namesorted.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]
    genrich_chip_string = ' -t ' + ' '.join(genrich_chip_list)

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        genrich_ctrl_list = ['{}/{}.namesorted.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]
        genrich_ctrl_string = ' -c ' + ' '.join(genrich_ctrl_list)

    else:
        genrich_ctrl_string = ''


    genrich_peak_calling_script.write('caras_Genrich.py{}{} -o {}/{}_Genrich.narrowPeak{} {} 1> {}/logs/{}.Genrich.out 2> {}/logs/{}.Genrich.err\n'.format(
        genrich_read_mode_arg,
        genrich_chip_string, 
        genrich_dir, 
        dataset_name,
        genrich_ctrl_string, 
        genrich_arg, 
        genrich_dir, 
        dataset_name,
        genrich_dir, 
        dataset_name))



if analysis_mode == 'single_cell':

    genrich_chip_list = [' -t {}/{}.namesorted.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        genrich_ctrl_list = [' -c {}/{}.namesorted.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]

    else:
        genrich_ctrl_list = [' ' for list_counter in range(len(chip_name))]
        

    for list_counter in range(len(chip_name)):
        genrich_peak_calling_script.write('caras_Genrich.py{}{} -o {}/{}_rep{}_Genrich.narrowPeak{} {} 1> {}/logs/{}_rep{}.Genrich.out 2> {}/logs/{}_rep{}.Genrich.err\n\n'.format(
            genrich_read_mode_arg,
            genrich_chip_list[list_counter], 
            genrich_dir, 
            dataset_name,
            (list_counter + 1), 
            genrich_ctrl_list[list_counter], 
            genrich_arg, 
            genrich_dir, 
            dataset_name,
            (list_counter + 1), 
            genrich_dir, 
            dataset_name,
            (list_counter + 1)))

genrich_peak_calling_script.close() # Closing the script '14_genrich_peak_calling_script.sh'. Flushing the write buffer



####################################################################################################################
### 15_seacr_peak_calling_script.sh
####################################################################################################################

# Create a directory "15_seacr_peak_calling" for the SEACR-called peak list file
# Create a script "15_seacr_peak_calling_script.sh" within, that calls for ChIP-seq peak calling program: SEACR
# Default mode is chosen here, which is optimized for peak calling in general ChIP-seq datasets
# Various output will be generated. The peak list file that is going to be used 
#   in subsequent analysis is the one in narrowPeak format, with .narrowPeak extension
# In: filename.namesorted.bam (output from 08_results_script.sh)
# Out: setname_SEACR.narrowPeak (input for 21_peaks_merging_script.sh)



seacr_dir = '{}/15_seacr_peak_calling'.format(output_dir)
seacr_peak_calling_script_name = '{}/15_seacr_peak_calling_script.sh'.format(seacr_dir)
if not os.path.exists(seacr_dir + '/logs'):
    os.makedirs(seacr_dir + '/logs')
seacr_peak_calling_script = open(seacr_peak_calling_script_name, 'w')

seacr_peak_calling_script.write('#!/bin/bash\n\n')
seacr_peak_calling_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'single_cell':
    seacr_mode_arg = ' norm relaxed'

elif ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    seacr_mode_arg = ' norm relaxed'

else:
    seacr_mode_arg = ' norm stringent'

popen_which_R = subprocess.Popen('which R', shell = True, stdout = subprocess.PIPE)
# Get the directory of R being used as default in the system
which_R_out = popen_which_R.communicate()[0]
default_R_directory = which_R_out.decode('utf-8').strip().rstrip('/R')


if analysis_mode == 'bulk':

    seacr_chip_string = ' {}/{}_chip_merged.normalized.normalized.bdg'.format(results_dir, dataset_name)

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        seacr_ctrl_string = ' {}/{}_ctrl_merged.normalized.normalized.bdg'.format(results_dir, dataset_name)

    else:
        seacr_ctrl_string = ' 0.01'

    seacr_peak_calling_script.write('SEACR_1.3.sh{}{}{}{} {}/{} {} 1> {}/logs/{}.SEACR.out 2> {}/logs/{}.SEACR.err\n'.format(
        seacr_chip_string,
        seacr_ctrl_string,
        seacr_arg,
        seacr_mode_arg, 
        seacr_dir, 
        dataset_name,
        default_R_directory,
        seacr_dir, 
        dataset_name,
        seacr_dir, 
        dataset_name))


if analysis_mode == 'single_cell':

    seacr_chip_list = [' {}/{}.normalized.normalized.bdg'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        seacr_ctrl_list = [' {}/{}.normalized.normalized.bdg'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]

    else:
        seacr_ctrl_list = [' 0.01' for list_counter in range(len(chip_name))]
        

    for list_counter in range(len(chip_name)):
        seacr_peak_calling_script.write('SEACR_1.3.sh{}{}{} {} {}/{}_rep{} {} 1> {}/logs/{}_rep{}.SEACR.out 2> {}/logs/{}_rep{}.SEACR.err &\n'.format(
            seacr_chip_list[list_counter], 
            seacr_ctrl_list[list_counter], 
            seacr_arg,
            seacr_mode_arg, 
            seacr_dir, 
            dataset_name,
            (list_counter + 1), 
            default_R_directory,
            seacr_dir, 
            dataset_name,
            (list_counter + 1), 
            seacr_dir, 
            dataset_name,
            (list_counter + 1)))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            seacr_peak_calling_script.write('\nwait\n\n')


seacr_peak_calling_script.close() # Closing the script '15_seacr_peak_calling_script.sh'. Flushing the write buffer




####################################################################################################################
### 16_sicer2_peak_calling_script.sh
####################################################################################################################

# Write this script only when CaRAS is running in broad peak mode (ChIP protein is a histone modifier)
# For narrow peaks (transcription factor ChIP protein), this script is replaced by 16_sicer2_peak_calling_script.sh 

# Create a directory "16_sicer2_peak_calling" for the SICER2-called peak list file
# Create a script "16_sicer2_peak_calling_script.sh" within, that calls for ChIP-seq peak calling program: SICER2
# In: filename.normalized.bam (output from 08_results_script.sh)
# Out: filename-W*-G*-islands-summary (input for 21_peaks_merging_script.sh)
#       * depends on -w and -g flag argument used when executing. Default is 200 and 600, respectively.

sicer2_dir = '{}/16_sicer2_peak_calling'.format(output_dir)
sicer2_peak_calling_script_name = '{}/16_sicer2_peak_calling_script.sh'.format(sicer2_dir)
if not os.path.exists(sicer2_dir + '/logs'):
    os.makedirs(sicer2_dir + '/logs')

sicer2_peak_calling_script = open(sicer2_peak_calling_script_name, 'w')

sicer2_peak_calling_script.write('#!/bin/bash\n\n')
sicer2_peak_calling_script.write('set -euxo pipefail\n\n')


# No particular optional flag is used here. Everything is on default mode.
# Default mode equals to executing with these arguments: -w 200 -rt 1 -f 150 -egf 0.74 -fdr 0.01 -g 600 -e 1000
sicer2_peak_calling_script.write('cd {} \n\n'.format(sicer2_dir)) # Change to output directory so SICER2 will dump the outputs there



if analysis_mode == 'bulk':

    sicer2_chip_string = ' -t ' + '{}/{}_chip_merged.normalized.bam'.format(results_dir, dataset_name)

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        sicer2_ctrl_string = ' -c ' + '{}/{}_ctrl_merged.normalized.bam'.format(results_dir, dataset_name)
    
    else:
        sicer2_ctrl_string = ''
    

    sicer2_peak_calling_script.write('sicer{}{} -cpu {} -s {} {} -rt 1000000000 1> {}/logs/{}.SICER2.out 2> {}/logs/{}.SICER2.err\n\n'.format(
        sicer2_chip_string,
        sicer2_ctrl_string,
        cpu_count,
        'mm10' if genome_ref == 'mm39' else genome_ref,
        sicer2_arg,
        sicer2_dir, 
        dataset_name, 
        sicer2_dir, 
        dataset_name))



if analysis_mode == 'single_cell':

    sicer2_chip_list = [' -t {}/{}.normalized.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]

    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        sicer2_ctrl_list = [' -c {}/{}.normalized.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]

    else:
        sicer2_ctrl_list = [' ' for list_counter in range(len(chip_name))]


    for list_counter in range(len(chip_name)):
        sicer2_peak_calling_script.write('sicer{}{} -cpu {} -s {} {} -rt 1000000000 1> {}/logs/{}_rep{}.SICER2.out 2> {}/logs/{}_rep{}.SICER2.err\n\n'.format(
            sicer2_chip_list[list_counter],
            sicer2_ctrl_list[list_counter],
            cpu_count,
            'mm10' if genome_ref == 'mm39' else genome_ref,
            sicer2_arg,
            sicer2_dir, 
            dataset_name, 
            (list_counter + 1),
            sicer2_dir, 
            dataset_name,
            (list_counter + 1)))

        sicer2_peak_calling_script.write('rm -rf {}.bed\n\n'.format(sicer2_chip_list[list_counter].lstrip(' -t ').rstrip('.bam')))

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            sicer2_peak_calling_script.write('rm -rf {}.bed\n\n'.format(sicer2_ctrl_list[list_counter].lstrip(' -c ').rstrip('.bam')))


sicer2_peak_calling_script.write('cd {} \n\n'.format(current_dir))

sicer2_peak_calling_script.close() # Closing the script '16_sicer2_peak_calling_script.sh'. Flushing the write buffer



####################################################################################################################
### 21_peaks_merging_script.sh
####################################################################################################################

# Create a directory "21_peaks_merging" for the reformated peak lists from all peak callers
# Create a script "21_peaks_merging_script.sh" within, that reads the different peak lists, 
#   reformating them into "HOMER custom_setting_table" with awk, and save them in this folder
# The script also calls for HOMER peak merging program: HOMER mergePeaks, that combines all the peak lists together, 
#   merging peaks within 100 bp vicinity (default settings; customizable) into one single peak 
# In: setname_MACS2_peaks.narrowPeak (output from 11_macs2_peak_calling_script.sh)
# In: setname_GEM_GPS_events.txt (output from 12_gem_peak_calling_script.sh)
# In: setname_HOMER.peaks (output from 13_homer_peak_calling_script.sh)
# In: setname_Genrich.narrowPeak (output from 14_genrich_peak_calling_script.sh)
# Out: setname_merged_peaks* (multiple files) (input for 22_peaks_processing_script.sh)

peaks_merging_dir = '{}/21_peaks_merging'.format(output_dir)
peaks_merging_script_name = '{}/21_peaks_merging_script.sh'.format(peaks_merging_dir)
if not os.path.exists(peaks_merging_dir):
    os.makedirs(peaks_merging_dir)
peaks_merging_script = open(peaks_merging_script_name, 'w')

peaks_merging_script.write('#!/bin/bash\n\n')
peaks_merging_script.write('set -euxo pipefail\n\n')


chromosomes_df = pd.read_csv('{}/chrom.sizes/{}.chrom.sizes'.format(genome_dir, genome_ref), header = None, usecols = [0], sep = '\t')
chromosomes_df.sort_values(by = [0], inplace = True)
chromosomes_list = chromosomes_df[0].values.tolist()
chromosomes_string = '\|'.join(chromosomes_list)


# Bash commands to remove all pre-existing peak caller combination files
peaks_merging_script.write('rm -r -f {}/{}_merged_peaks*\n\n'.format(peaks_merging_dir, dataset_name))


if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    gem_sorter_column = 4
    homer_sorter_column = 11
    sicer2_sorter_column = 7
    sicer2_peakset_filename = 'islands-summary'

else:
    gem_sorter_column = 5
    homer_sorter_column = 6
    sicer2_sorter_column = 4
    sicer2_peakset_filename = 'scoreisland'



if analysis_mode == 'bulk':

    # Bash commands to reformat MACS2 peak list file (.narrowPeak) into HOMER custom_setting_table format 
    #   and save it as extension-less text file in this folder under the name MACS2 
    peaks_merging_script.write('cat {}/{}_MACS2*_peaks.narrowPeak'.format(macs2_dir, dataset_name)) # Read called peaks list by MACS2
    peaks_merging_script.write(r""" | awk '$3-$2 < 10000'""")
    peaks_merging_script.write(' | sort -k 7 -n -r') # Run reversed (-r) numerical (-n) sort by the signalValue (column 7), causing the list to start with highest scoring peaks
    peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest signalValue. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1,$2,$3}'""") # Get only the chr column ($1), start column ($2), and end column ($3)
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
    peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
    peaks_merging_script.write(' > {}/MACS2'.format(peaks_merging_dir)) # Save it under a short name (as to not generate mergePeaks output filename too long)
    peaks_merging_script.write(' &\n\n')

    # Bash commands to reformat GEM peak list file (GPS_events.txt) into HOMER custom_setting_table format 
    #   and save it as extension-less text file in this folder under the name GEM 
    peaks_merging_script.write('cat {}/{}_GEM*GPS_events.txt'.format(gem_dir, dataset_name))
    peaks_merging_script.write(' | tail -n +2') # Throw out the headers
    peaks_merging_script.write(' | sort -k {} -n -r'.format(gem_sorter_column))
    peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest fold change. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_merging_script.write(r""" | awk '{split($1,a,":"); OFS="\t"; print "chr"a[1],a[2]-150,a[2]+150}'""") # Get the 1 bp coordinate generated by GEM. Extend left and right to 50 bp.
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
    peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
    peaks_merging_script.write(' > {}/GEM'.format(peaks_merging_dir)) # Save it under a short name (as to not generate mergePeaks output filename too long)
    peaks_merging_script.write(' &\n\n')

    # Bash commands to reformat HOMER peak list file (user-determined filename) into HOMER custom_setting_table format 
    #   and save it as extension-less text file in this folder under the name HOMER 
    peaks_merging_script.write('cat {}/{}_HOMER.peaks'.format(homer_dir, dataset_name))
    peaks_merging_script.write(' | tail -n +42') # Throw out the headers
    peaks_merging_script.write(' | sort -k {} -n -r'.format(homer_sorter_column))
    peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest fold change. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $2,$3,$4}'""") # Get only the chr column ($2), start column ($3), and end column ($4)
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
    peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
    peaks_merging_script.write(' > {}/HOMER'.format(peaks_merging_dir)) # Save it under a short name (as to not generate mergePeaks output filename too long)
    peaks_merging_script.write(' &\n\n')

    # Bash commands to reformat Genrich peak list file (.narrowPeak) into HOMER custom_setting_table format 
    #   and save it as extension-less text file in this folder under the name Genrich 
    peaks_merging_script.write('cat {}/{}_Genrich.narrowPeak'.format(genrich_dir, dataset_name))
    peaks_merging_script.write(' | sort -k 7 -n -r') # Run reversed (-r) numerical (-n) sort by the signalValue (column 7), causing the list to start with highest scoring peaks
    peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest signalValue. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1,$2,$3}'""") # Get only the chr column ($1), start column ($2), and end column ($3)
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
    peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
    peaks_merging_script.write(' > {}/Genrich'.format(peaks_merging_dir)) # Save it under a short name (as to not generate mergePeaks output filename too long)
    peaks_merging_script.write(' &\n\n')

    # Bash commands to reformat SEACR peak list file (.bed) into HOMER custom_setting_table format 
    #   and save it as extension-less text file in this folder under the name SEACR 
    peaks_merging_script.write('cat {}/{}*.bed'.format(seacr_dir, dataset_name))
    peaks_merging_script.write(' | sort -k 4 -n -r')
    peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest signalValue. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1,$2,$3}'""") # Get only the chr column ($1), start column ($2), and end column ($3)
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
    peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
    peaks_merging_script.write(' > {}/SEACR'.format(peaks_merging_dir)) # Save it under a short name (as to not generate mergePeaks output filename too long)
    peaks_merging_script.write(' &\n\n')

    # Bash commands to reformat GEM peak list file (GPS_events.txt) into HOMER custom_setting_table format 
    #   and save it as extension-less text file in this folder under the name GEM 
    peaks_merging_script.write('cat {}/*W*G*{}'.format(sicer2_dir, sicer2_peakset_filename))
    peaks_merging_script.write(' | sort -k {} -n -r'.format(sicer2_sorter_column))
    peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest fold change. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1,$2,$3}'""") # Get only the chr column ($1), start column ($2), and end column ($3)
    peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
    peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
    peaks_merging_script.write(' > {}/SICER2'.format(peaks_merging_dir)) # Save it under a short name (as to not generate mergePeaks output filename too long)
    peaks_merging_script.write(' &\n\n')

    peaks_merging_script.write('wait\n\n')

    # Change into the working directory where the input files are located. 
    #   HOMER mergePeaks program will create ridiculously long output file name, otherwise. 
    peaks_merging_script.write('cd {}\n\n'.format(peaks_merging_dir))

    # Bash commands to merge all peak lists, generating multiple output peak lists 
    #   that is called by each possible combinations of peak callers
    peaks_merging_script.write('mergePeaks {} MACS2 GEM HOMER Genrich SEACR SICER2 -prefix {}_merged_peaks -venn venn.txt\n\n'.format(
        homer_mergePeaks_arg, 
        dataset_name))

    peaks_merging_script.write('caras_upset_plotter.py --venn {}/venn.txt --prefix {} 1> /dev/null 2> /dev/null\n\n'.format(peaks_merging_dir, dataset_name))



if analysis_mode == 'single_cell':
    # Commands to extract peak locations from peak caller output of MACS2 (narrow peak mode), GEM, HOMER (factor peak mode), and Genrich
    # Only when CaRAS is running in narrow peak mode (ChIP protein is a transcription factor)
    for list_counter in range(len(chip_name)):
        # Bash commands to reformat MACS2 peak list file (.narrowPeak) into HOMER custom_setting_table format 
        #   and save it as extension-less text file in this folder under the name MACS2 
        peaks_merging_script.write('cat {}/{}_rep{}_MACS2*_peaks.narrowPeak'.format(macs2_dir, dataset_name, (list_counter + 1))) # Read called peaks list by MACS2
        peaks_merging_script.write(r""" | awk '$3-$2 < 10000'""")
        peaks_merging_script.write(' | sort -k 7 -n -r') # Run reversed (-r) numerical (-n) sort by the signalValue (column 7), the list to start with highest scoring peaks
        peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest signalValue. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1,$2,$3}'""") # Get only the chr column ($1), start column ($2), and end column ($3)
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
        peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
        peaks_merging_script.write(' > {}/MACS2_{}'.format(peaks_merging_dir, (list_counter + 1))) # Save it under a short name (as to not generate mergePeaks output filename too long)
        peaks_merging_script.write(' &\n')

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            peaks_merging_script.write('\nwait\n\n')

    for list_counter in range(len(chip_name)):
        # Bash commands to reformat GEM peak list file (GPS_events.txt) into HOMER custom_setting_table format 
        #   and save it as extension-less text file in this folder under the name GEM 
        peaks_merging_script.write('cat {}/{}_rep{}_GEM*GPS_events.txt'.format(gem_dir, dataset_name, (list_counter + 1)))
        peaks_merging_script.write(' | tail -n +2') # Throw out the headers
        peaks_merging_script.write(' | sort -k {} -n -r'.format(gem_sorter_column))
        peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest fold change. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
        peaks_merging_script.write(r""" | awk '{split($1,a,":"); OFS="\t"; print "chr"a[1],a[2]-150,a[2]+150}'""") # Get the 1 bp coordinate from GEM. Extend it left and right to 50 bp.
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
        peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
        peaks_merging_script.write(' > {}/GEM_{}'.format(peaks_merging_dir, (list_counter + 1))) # Save it under a short name (as to not generate mergePeaks output filename too long)
        peaks_merging_script.write(' &\n')

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            peaks_merging_script.write('\nwait\n\n')

    for list_counter in range(len(chip_name)):
        # Bash commands to reformat HOMER peak list file (user-determined filename) into HOMER custom_setting_table format 
        #   and save it as extension-less text file in this folder under the name HOMER 
        peaks_merging_script.write('cat {}/{}_rep{}_HOMER.peaks'.format(homer_dir, dataset_name, (list_counter + 1)))
        peaks_merging_script.write(' | tail -n +42') # Throw out the headers
        peaks_merging_script.write(' | sort -k {} -n -r'.format(homer_sorter_column))
        peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest fold change. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $2,$3,$4}'""") # Get only the chr column ($2), start column ($3), and end column ($4)
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
        peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
        peaks_merging_script.write(' > {}/HOMER_{}'.format(peaks_merging_dir, (list_counter + 1))) # Save it under a short name (as to not generate mergePeaks output filename too long)
        peaks_merging_script.write(' &\n')

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            peaks_merging_script.write('\nwait\n\n')

    for list_counter in range(len(chip_name)):
        # Bash commands to reformat Genrich peak list file (.narrowPeak) into HOMER custom_setting_table format 
        #   and save it as extension-less text file in this folder under the name Genrich 
        peaks_merging_script.write('cat {}/{}_rep{}_Genrich.narrowPeak'.format(genrich_dir, dataset_name, (list_counter + 1)))
        peaks_merging_script.write(' | sort -k 7 -n -r') # Run reversed (-r) numerical (-n) sort by the signalValue (column 7), the list to start with highest scoring peaks
        peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest signalValue. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1,$2,$3}'""") # Get only the chr column ($1), start column ($2), and end column ($3)
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
        peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
        peaks_merging_script.write(' > {}/Genrich_{}'.format(peaks_merging_dir, (list_counter + 1))) # Save it under a short name (as to not generate output filename too long)
        peaks_merging_script.write(' &\n')

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            peaks_merging_script.write('\nwait\n\n')

    for list_counter in range(len(chip_name)):
        # Bash commands to reformat SEACR peak list file (.bed) into HOMER custom_setting_table format 
        #   and save it as extension-less text file in this folder under the name SEACR 
        peaks_merging_script.write('cat {}/{}_rep{}*.bed'.format(seacr_dir, dataset_name, (list_counter + 1)))
        peaks_merging_script.write(' | sort -k 4 -n -r') # Run reversed (-r) numerical (-n) sort by the signalValue (column 7), causing the list to start with highest scoring peaks
        peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest signalValue. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1,$2,$3}'""") # Get only the chr column ($1), start column ($2), and end column ($3)
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
        peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
        peaks_merging_script.write(' > {}/SEACR_{}'.format(peaks_merging_dir, (list_counter + 1))) # Save it under a short name (as to not generate mergePeaks output filename too long)
        peaks_merging_script.write(' &\n')

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            peaks_merging_script.write('\nwait\n\n')

    # Change into the working directory where the input files are located. 
    #   HOMER mergePeaks program will create ridiculously long output file name, otherwise. 
    peaks_merging_script.write('cd {}\n\n'.format(peaks_merging_dir))

    for list_counter in range(len(chip_name)):
        # Bash commands to reformat SICER2 peak list file (.bed) into HOMER custom_setting_table format 
        #   and save it as extension-less text file in this folder under the name SICER2
        peaks_merging_script.write('cat {}/*rep{}*W*G*{}'.format(sicer2_dir, (list_counter + 1), sicer2_peakset_filename))
        peaks_merging_script.write(' | sort -k {} -n -r'.format(sicer2_sorter_column))
        peaks_merging_script.write(' | head -200000') # Take the 200,000 peaks with highest fold change. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1,$2,$3}'""") # Get only the chr column ($1), start column ($2), and end column ($3)
        peaks_merging_script.write(r""" | awk '{OFS="\t";print $1":"$2"-"$3,$1,$2,$3,"+"}'""") # Make it into HOMER format (chr:start-end \t chr \t start \t end \t strand)
        peaks_merging_script.write(" | grep -w '{}'".format(chromosomes_string)) # Filter out chr_alt, chr_fix, chrN_random, chrUn, and chrM
        peaks_merging_script.write(' > {}/SICER2_{}'.format(peaks_merging_dir, (list_counter + 1))) # Save it under a short name (as to not generate output filename too long)
        peaks_merging_script.write(' &\n')

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            peaks_merging_script.write('\nwait\n\n')


    for list_counter in range(len(chip_name)):
        # Bash commands to merge all peak lists, generating multiple output peak lists 
        #   that is called by each possible combinations of peak callers
        peaks_merging_script.write('mergePeaks {} MACS2_{} GEM_{} HOMER_{} Genrich_{} SEACR_{} SICER2_{} -prefix {}_rep{}_merged_peaks -venn venn_rep{}.txt\n\n'.format(
            homer_mergePeaks_arg, 
            (list_counter + 1),
            (list_counter + 1),
            (list_counter + 1),
            (list_counter + 1),
            (list_counter + 1),
            (list_counter + 1),
            dataset_name,
            (list_counter + 1),
            (list_counter + 1)))


    for list_counter in range(len(chip_name)):
        peaks_merging_script.write('caras_upset_plotter.py --venn {}/venn.txt --prefix {}_rep{} 1> /dev/null 2> /dev/null\n\n'.format(peaks_merging_dir, dataset_name, (list_counter + 1)))


# Going back to the original directory from which the suite was called, 
#   just so not to cause unnecessary confusion for the user.
peaks_merging_script.write('cd {}\n\n'.format(current_dir))

peaks_merging_script.close() # Closing the script '21_peaks_merging_script.sh'. Flushing the write buffer



####################################################################################################################
### 22_peaks_processing_script.sh
####################################################################################################################

# Create a directory "22_peaks_processing" for the full peak list with gene annotation and recalculated peak stats
# Create a script "22_peaks_processing_script.sh" within, that concatenates all peaks 
#   from all possible peak caller combination files, appends gene annotations to each peak, 
#   and appends various useful peak stats to each peak; resulting in a complete custom_setting_table 
#   in .tsv format that can be easily viewed, sorted, and filtered in spreadsheet applications.
# The script will also generate a separate summary custom_setting_table file that shows 
#   the statistics of peaks that are called each combination of peak callers.
# In: setname_merged_peaks* (multiple files) (output from 21_peaks_merging_script.sh)
# Out: setname_all_peaks_calculated.tsv (input for 23_go_annotation_script.sh or 
#   23_pathway_annotation_script.sh or 23_go_pathway_annotation_script.sh)

peaks_processing_dir = '{}/22_peaks_processing'.format(output_dir)
peaks_processing_script_name = '{}/22_peaks_processing_script.sh'.format(peaks_processing_dir)
if analysis_mode == 'bulk' and not os.path.exists(peaks_processing_dir + '/IDR_files'):
    os.makedirs(peaks_processing_dir + '/IDR_files')
peaks_processing_script = open(peaks_processing_script_name, 'w')

peaks_processing_script.write('#!/bin/bash\n\n')
peaks_processing_script.write('set -euxo pipefail\n\n')



if analysis_mode == 'bulk':

    # Bash commands to concatenate all peaks from all possible peak caller combination files
    peaks_processing_script.write('cat {}/{}_merged_peaks* | grep -v name > {}/{}_all_peaks_concatenated.tsv\n\n'.format(
        peaks_merging_dir, 
        dataset_name, 
        peaks_processing_dir, 
        dataset_name))


    # Bash commands to call HOMER annotatePeaks.pl to append gene annotations to each peak
    peaks_processing_script.write('annotatePeaks.pl {} {}/{}_all_peaks_concatenated.tsv {}{} -nmotifs -go {}/{}_gene_ontology > {}/{}_all_peaks_annotated.tsv\n\n'.format(
        homer_annotatePeaks_arg,
        peaks_processing_dir, 
        dataset_name, 
        '{}/bwa/{}.fa -gtf {}/annotation/{}.refGene.gtf'.format(genome_dir, genome_ref, genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        ' -m {}'.format(motif_file_full_path) if motif_file_full_path != None else '', 
        peaks_processing_dir, 
        dataset_name, 
        peaks_processing_dir, 
        dataset_name))


    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:

        if force_merge == 0:

            # See comments at the top, force merge == 0 means the fold change analysis later on 
            #   will be calculated based each paired replicates of ChIP and control .bam files
            # Converting both ChIP and control .bam lists from vertical lists into space-separated 
            #   serial string (to be used as an input argument in MACS2)
            fcc_chip_list = ['{}/{}.normalized.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]
            fcc_ctrl_list = ['{}/{}.normalized.bam'.format(results_dir, ctrl_name[list_counter]) for list_counter in range(len(ctrl_name))]
            fcc_chip_string = ' '.join(fcc_chip_list)
            fcc_ctrl_string = ' '.join(fcc_ctrl_list)


            # Bash commands to call the home-made script caras_fold_change_calculator.py 
            #   to append tag counts, fold changes, motif hits, and peak caller overlaps stats to each peak
            peaks_processing_script.write('caras_fold_change_calculator.py {} --thread {} --peak {} --input_tsv {}/{}_all_peaks_annotated.tsv --output_tsv {}/{}_all_peaks_calculated.tsv --chip_bam {} --ctrl_bam {}\n\n'.format(
                fold_change_calculator_arg,
                cpu_count, 
                peak_type,
                peaks_processing_dir, 
                dataset_name, 
                peaks_processing_dir, 
                dataset_name, 
                fcc_chip_string, 
                fcc_ctrl_string))


        if force_merge == 1:

            # See comments at the top, force merge == 1 means the fold change analysis later on 
            #   will be calculated solely based on merged ChIP .bam and merged control .bam files
            # Bash commands to call the home-made script caras_fold_change_calculator.py 
            #   to append tag counts, fold changes, motif hits, and peak caller overlaps stats to each peak
            peaks_processing_script.write('caras_fold_change_calculator.py {} --thread {} --peak {} --input_tsv {}/{}_all_peaks_annotated.tsv --output_tsv {}/{}_all_peaks_calculated.tsv --chip_bam {}/{}_chip_merged.normalized.bam --ctrl_bam {}/{}_ctrl_merged.normalized.bam\n\n'.format(
                fold_change_calculator_arg, 
                cpu_count, 
                peak_type,
                peaks_processing_dir, 
                dataset_name, 
                peaks_processing_dir, 
                dataset_name, 
                results_dir, 
                dataset_name, 
                results_dir, 
                dataset_name))


    else:

        if force_merge == 0:

            # See comments at the top, force merge == 0 means the fold change analysis later on 
            #   will be calculated based each paired replicates of ChIP and control .bam files
            # Converting both ChIP and control .bam lists from vertical lists into space-separated 
            #   serial string (to be used as an input argument in MACS2)
            fcc_chip_list = ['{}/{}.normalized.bam'.format(results_dir, chip_name[list_counter]) for list_counter in range(len(chip_name))]
            fcc_chip_string = ' '.join(fcc_chip_list)


            # Bash commands to call the home-made script caras_fold_change_calculator.py 
            #   to append tag counts, fold changes, motif hits, and peak caller overlaps stats to each peak
            peaks_processing_script.write('caras_fold_change_calculator.py {} --thread {} --peak {} --input_tsv {}/{}_all_peaks_annotated.tsv --output_tsv {}/{}_all_peaks_calculated.tsv --chip_bam {}\n\n'.format(
                fold_change_calculator_arg,
                cpu_count, 
                peak_type,
                peaks_processing_dir, 
                dataset_name, 
                peaks_processing_dir, 
                dataset_name, 
                fcc_chip_string))


        if force_merge == 1:

            # See comments at the top, force merge == 1 means the fold change analysis later on 
            #   will be calculated solely based on merged ChIP .bam and merged control .bam files
            # Bash commands to call the home-made script caras_fold_change_calculator.py 
            #   to append tag counts, fold changes, motif hits, and peak caller overlaps stats to each peak
            peaks_processing_script.write('caras_fold_change_calculator.py {} --thread {} --peak {} --input_tsv {}/{}_all_peaks_annotated.tsv --output_tsv {}/{}_all_peaks_calculated.tsv --chip_bam {}/{}_chip_merged.normalized.bam\n\n'.format(
                fold_change_calculator_arg, 
                cpu_count, 
                peak_type,
                peaks_processing_dir, 
                dataset_name, 
                peaks_processing_dir, 
                dataset_name, 
                results_dir, 
                dataset_name))


    # Bash commands to reformat MACS2 peak list file (.narrowPeak) into a ranked list for IDR calculation
    peaks_processing_script.write('cat {}/{}_MACS2_peaks.narrowPeak'.format(macs2_dir, dataset_name)) # Read called peaks list by MACS2
    peaks_processing_script.write(' | sort -k 7 -n -r') # Run reversed (-r) numerical (-n) sort by the signalValue (column 7), causing the list to start with highest scoring peaks
    peaks_processing_script.write(' | head -200000') # Take the 200,000 peaks with highest signalValue. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$7}'""") # Get only the chr column ($1), start column ($2), end column ($3), and signalValue column ($7)
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$1":"$2"-"$3,"0",".",$4,"-1","-1","-1"}'""") # Make it into narrowPeak format
    peaks_processing_script.write(' > {}/IDR_files/MACS2_IDR_input.narrowPeak'.format(peaks_processing_dir)) # Save it as MACS2_IDR_input.narrowPeak
    peaks_processing_script.write(' &\n\n')

    # Bash commands to reformat GEM peak list file (GPS_events.txt) into a ranked list for IDR calculation
    peaks_processing_script.write('cat {}/{}_GEM*GPS_events.txt'.format(gem_dir, dataset_name))
    peaks_processing_script.write(' | tail -n +2') # Throw out the headers
    peaks_processing_script.write(' | sort -k {} -n -r'.format(gem_sorter_column))
    peaks_processing_script.write(' | head -200000') # Take the 200,000 peaks with highest fold change. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_processing_script.write(r""" | awk '{split($1,a,":"); OFS="\t"; print "chr"a[1],a[2]-150,a[2]+150""") # Get the 1 bp coordinate generated by GEM. Extend it L&R to 50 bp.
    peaks_processing_script.write(',${}'.format(gem_sorter_column) + "}'")
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$1":"$2"-"$3,"0",".",$4,"-1","-1","-1"}'""") # Make it into narrowPeak format
    peaks_processing_script.write(' > {}/IDR_files/GEM_IDR_input.narrowPeak'.format(peaks_processing_dir)) # Save it as GEM_IDR_input.narrowPeak
    peaks_processing_script.write(' &\n\n')

    # Bash commands to reformat HOMER peak list file (user-determined filename) into a ranked list for IDR calculation
    peaks_processing_script.write('cat {}/{}_HOMER.peaks'.format(homer_dir, dataset_name))
    peaks_processing_script.write(' | tail -n +42') # Throw out the headers
    peaks_processing_script.write(' | sort -k {} -n -r'.format(homer_sorter_column))
    peaks_processing_script.write(' | head -200000') # Take the 200,000 peaks with highest fold change. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $2,$3,$4""") # Get only the chr column ($2), start column ($3), end column ($4), and fold change column ($11)
    peaks_processing_script.write(',${}'.format(homer_sorter_column) + "}'")
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$1":"$2"-"$3,"0",".",$4,"-1","-1","-1"}'""") # Make it into narrowPeak format
    peaks_processing_script.write(' > {}/IDR_files/HOMER_IDR_input.narrowPeak'.format(peaks_processing_dir)) # Save it as HOMER_IDR_input.narrowPeak
    peaks_processing_script.write(' &\n\n')

    # Bash commands to reformat Genrich peak list file (.narrowPeak) into a ranked list for IDR calculation
    peaks_processing_script.write('cat {}/{}_Genrich.narrowPeak'.format(genrich_dir, dataset_name))
    peaks_processing_script.write(' | sort -k 7 -n -r') # Run reversed (-r) numerical (-n) sort by the signalValue (column 7), causing the list to start with highest scoring peaks
    peaks_processing_script.write(' | head -200000') # Take the 200,000 peaks with highest signalValue. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$7}'""") # Get only the chr column ($1), start column ($2), end column ($3), and signalValue column ($7)
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$1":"$2"-"$3,"0",".",$4,"-1","-1","-1"}'""") # Make it into narrowPeak format
    peaks_processing_script.write(' > {}/IDR_files/Genrich_IDR_input.narrowPeak'.format(peaks_processing_dir)) # Save it as Genrich_IDR_input.narrowPeak
    peaks_processing_script.write(' &\n\n')

    # Bash commands to reformat SEACR peak list file (.bed) into HOMER custom_setting_table format 
    #   and save it as extension-less text file in this folder under the name SEACR 
    peaks_processing_script.write('cat {}/{}*.bed'.format(seacr_dir, dataset_name))
    peaks_processing_script.write(' | sort -k 4 -n -r') # Run reversed (-r) numerical (-n) sort by the signalValue (column 4), causing the list to start with highest scoring peaks
    peaks_processing_script.write(' | head -200000') # Take the 200,000 peaks with highest signalValue. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$4}'""") # Get only the chr column ($1), start column ($2), end column ($3), and signalValue column ($4)
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$1":"$2"-"$3,"0",".",$4,"-1","-1","-1"}'""") # Make it into narrowPeak format
    peaks_processing_script.write(' > {}/IDR_files/SEACR_IDR_input.narrowPeak'.format(peaks_processing_dir)) # Save it as SEACR_IDR_input.narrowPeak
    peaks_processing_script.write(' &\n\n')

    # Bash commands to reformat SICER2 peak list file (*W*G*islands-summary) into a ranked list for IDR calculation 
    peaks_processing_script.write('cat {}/{}*W*G*{}'.format(sicer2_dir, dataset_name, sicer2_peakset_filename))
    peaks_processing_script.write(' | sort -k {} -n -r'.format(sicer2_sorter_column))
    peaks_processing_script.write(' | head -200000') # Take the 200,000 peaks with highest fold change. Anomalies happen sometimes where a peak caller calls 1,000,000+ peaks.
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3""") # Get only the chr column ($1), start column ($2), end column ($3), and fold change column ($7)
    peaks_processing_script.write(',${}'.format(sicer2_sorter_column) + "}'")
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$1":"$2"-"$3,"0",".",$4,"-1","-1","-1"}'""") # Make it into narrowPeak format
    peaks_processing_script.write(' > {}/IDR_files/SICER2_IDR_input.narrowPeak'.format(peaks_processing_dir)) # Save it as SICER2_IDR_input.narrowPeak
    peaks_processing_script.write(' &\n\n')
    
    peaks_processing_script.write('wait\n\n')


    # Bash command to run IDR module to calculate peak IDR based on union peak set versus MACS2 peak set
    peaks_processing_script.write('idr -s {}/{}_all_peaks_calculated.narrowPeak {}/IDR_files/MACS2_IDR_input.narrowPeak -o {}/IDR_files/MACS2_IDR_output.tsv --input-file-type narrowPeak --rank signal.value --use-nonoverlapping-peaks &\n\n'.format(
        peaks_processing_dir, 
        dataset_name,
        peaks_processing_dir,
        peaks_processing_dir))

    # Bash command to run IDR module to calculate peak IDR based on union peak set versus GEM peak set
    peaks_processing_script.write('idr -s {}/{}_all_peaks_calculated.narrowPeak {}/IDR_files/GEM_IDR_input.narrowPeak -o {}/IDR_files/GEM_IDR_output.tsv --input-file-type narrowPeak --rank signal.value --use-nonoverlapping-peaks &\n\n'.format(
        peaks_processing_dir, 
        dataset_name,
        peaks_processing_dir,
        peaks_processing_dir))

    # Bash command to run IDR module to calculate peak IDR based on union peak set versus HOMER peak set
    peaks_processing_script.write('idr -s {}/{}_all_peaks_calculated.narrowPeak {}/IDR_files/HOMER_IDR_input.narrowPeak -o {}/IDR_files/HOMER_IDR_output.tsv --input-file-type narrowPeak --rank signal.value --use-nonoverlapping-peaks &\n\n'.format(
        peaks_processing_dir, 
        dataset_name,
        peaks_processing_dir,
        peaks_processing_dir))

    # Bash command to run IDR module to calculate peak IDR based on union peak set versus Genrich peak set
    peaks_processing_script.write('idr -s {}/{}_all_peaks_calculated.narrowPeak {}/IDR_files/Genrich_IDR_input.narrowPeak -o {}/IDR_files/Genrich_IDR_output.tsv --input-file-type narrowPeak --rank signal.value --use-nonoverlapping-peaks &\n\n'.format(
        peaks_processing_dir, 
        dataset_name,
        peaks_processing_dir,
        peaks_processing_dir))

    # Bash command to run IDR module to calculate peak IDR based on union peak set versus SEACR peak set
    peaks_processing_script.write('idr -s {}/{}_all_peaks_calculated.narrowPeak {}/IDR_files/SEACR_IDR_input.narrowPeak -o {}/IDR_files/SEACR_IDR_output.tsv --input-file-type narrowPeak --rank signal.value --use-nonoverlapping-peaks &\n\n'.format(
        peaks_processing_dir, 
        dataset_name,
        peaks_processing_dir,
        peaks_processing_dir))

    # Bash command to run IDR module to calculate peak IDR based on union peak set versus SICER2 peak set
    peaks_processing_script.write('idr -s {}/{}_all_peaks_calculated.narrowPeak {}/IDR_files/SICER2_IDR_input.narrowPeak -o {}/IDR_files/SICER2_IDR_output.tsv --input-file-type narrowPeak --rank signal.value --use-nonoverlapping-peaks &\n\n'.format(
        peaks_processing_dir, 
        dataset_name,
        peaks_processing_dir,
        peaks_processing_dir))

    peaks_processing_script.write('wait\n\n')


    # Bash command to run the caras_IDR_integrator.py to integrate the results above into the full peak list dataset_name_all_peaks_calculated.tsv
    peaks_processing_script.write('caras_IDR_integrator.py --input_tsv {}/{}_all_peaks_calculated.tsv --idr_tsv {}/IDR_files/MACS2_IDR_output.tsv {}/IDR_files/GEM_IDR_output.tsv {}/IDR_files/HOMER_IDR_output.tsv {}/IDR_files/Genrich_IDR_output.tsv {}/IDR_files/SEACR_IDR_output.tsv {}/IDR_files/SICER2_IDR_output.tsv --output_tsv {}/{}_all_peaks_calculated.tsv\n\n'.format(
        peaks_processing_dir, 
        dataset_name,
        peaks_processing_dir,
        peaks_processing_dir,
        peaks_processing_dir,
        peaks_processing_dir,
        peaks_processing_dir,
        peaks_processing_dir,
        peaks_processing_dir,
        dataset_name))    


    # Bash commands to call the home-made script caras_peak_caller_stats_calculator.py to generate 
    #   a separate summary custom_setting_table file that shows the statistics 
    #   (tag counts, fold changes, motif hits, motif hits rate, positive peaks rate) 
    #   of peaks that are called each combination of peak callers
    peaks_processing_script.write('caras_peak_caller_stats_calculator.py {}/{}_all_peaks_calculated.tsv {}/{}_peak_caller_combinations_statistics.tsv\n\n'.format(
        peaks_processing_dir, 
        dataset_name, 
        peaks_processing_dir, 
        dataset_name))



if analysis_mode == 'single_cell':

    for list_counter in range(len(chip_name)):
        # Bash commands to concatenate all peaks from all possible peak caller combination files
        peaks_processing_script.write('cat {}/{}_rep{}_merged_peaks* | grep -v name > {}/{}_rep{}_all_peaks_concatenated.tsv &\n'.format(
            peaks_merging_dir, 
            dataset_name, 
            (list_counter + 1),
            peaks_processing_dir, 
            dataset_name,
            (list_counter + 1)))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(chip_name) - 1):
            peaks_processing_script.write('\nwait\n\n')


    for list_counter in range(len(chip_name)):
        # Bash commands to call HOMER annotatePeaks.pl to append gene annotations to each peak
        peaks_processing_script.write('annotatePeaks.pl {} {}/{}_rep{}_all_peaks_concatenated.tsv {}{} -nmotifs -go {}/{}_rep{}_gene_ontology > {}/{}_rep{}_all_peaks_annotated.tsv \n\n'.format(
            homer_annotatePeaks_arg,
            peaks_processing_dir, 
            dataset_name, 
            (list_counter + 1),
            '{}/bwa/{}.fa -gtf {}/annotation/{}.refGene.gtf'.format(genome_dir, genome_ref, genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
            ' -m {}'.format(motif_file_full_path) if motif_file_full_path != None else '', 
            peaks_processing_dir, 
            dataset_name, 
            (list_counter + 1),
            peaks_processing_dir, 
            dataset_name, 
            (list_counter + 1)))


    for list_counter in range(len(chip_name)):

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:

            if force_merge == 0:

                # See comments at the top, force merge == 0 means the fold change analysis later on 
                #   will be calculated based each paired replicates of ChIP and control .bam files
                # Converting both ChIP and control .bam lists from vertical lists into space-separated 
                #   serial string (to be used as an input argument in MACS2)
                fcc_chip_string = '{}/{}.normalized.bam'.format(results_dir, chip_name[list_counter])
                fcc_ctrl_string = '{}/{}.normalized.bam'.format(results_dir, ctrl_name[list_counter])


                # Bash commands to call the home-made script caras_fold_change_calculator.py 
                #   to append tag counts, fold changes, motif hits, and peak caller overlaps stats to each peak
                peaks_processing_script.write('caras_fold_change_calculator.py {} --thread {} --peak {} --input_tsv {}/{}_rep{}_all_peaks_annotated.tsv --output_tsv {}/{}_rep{}_all_peaks_calculated.tsv --chip_bam {} --ctrl_bam {}\n\n'.format(
                    fold_change_calculator_arg,
                    cpu_count, 
                    peak_type,
                    peaks_processing_dir, 
                    dataset_name, 
                    (list_counter + 1),
                    peaks_processing_dir, 
                    dataset_name, 
                    (list_counter + 1),
                    fcc_chip_string, 
                    fcc_ctrl_string))


            if force_merge == 1:

                # See comments at the top, force merge == 1 means the fold change analysis later on 
                #   will be calculated solely based on merged ChIP .bam and merged control .bam files
                # Bash commands to call the home-made script caras_fold_change_calculator.py 
                #   to append tag counts, fold changes, motif hits, and peak caller overlaps stats to each peak
                peaks_processing_script.write('caras_fold_change_calculator.py {} --thread {} --peak {} --input_tsv {}/{}_rep{}_all_peaks_annotated.tsv --output_tsv {}/{}_rep{}_all_peaks_calculated.tsv --chip_bam {}/{}_chip_merged.normalized.bam --ctrl_bam {}/{}_ctrl_merged.normalized.bam\n\n'.format(
                    fold_change_calculator_arg, 
                    cpu_count, 
                    peak_type,
                    peaks_processing_dir, 
                    dataset_name, 
                    (list_counter + 1),
                    peaks_processing_dir, 
                    dataset_name, 
                    (list_counter + 1),
                    results_dir, 
                    dataset_name, 
                    results_dir, 
                    dataset_name))


        else:

            if force_merge == 0:

                # See comments at the top, force merge == 0 means the fold change analysis later on 
                #   will be calculated based each paired replicates of ChIP and control .bam files
                # Converting both ChIP and control .bam lists from vertical lists into space-separated 
                #   serial string (to be used as an input argument in MACS2)
                fcc_chip_string = '{}/{}.normalized.bam'.format(results_dir, chip_name[list_counter])


                # Bash commands to call the home-made script caras_fold_change_calculator.py 
                #   to append tag counts, fold changes, motif hits, and peak caller overlaps stats to each peak
                peaks_processing_script.write('caras_fold_change_calculator.py {} --thread {} --peak {} --input_tsv {}/{}_rep{}_all_peaks_annotated.tsv --output_tsv {}/{}_rep{}_all_peaks_calculated.tsv --chip_bam {}\n\n'.format(
                    fold_change_calculator_arg,
                    cpu_count, 
                    peak_type,
                    peaks_processing_dir, 
                    dataset_name, 
                    (list_counter + 1),
                    peaks_processing_dir, 
                    dataset_name, 
                    (list_counter + 1),
                    fcc_chip_string))


            if force_merge == 1:

                # See comments at the top, force merge == 1 means the fold change analysis later on 
                #   will be calculated solely based on merged ChIP .bam and merged control .bam files
                # Bash commands to call the home-made script caras_fold_change_calculator.py 
                #   to append tag counts, fold changes, motif hits, and peak caller overlaps stats to each peak
                peaks_processing_script.write('caras_fold_change_calculator.py {} --thread {} --peak {} --input_tsv {}/{}_rep{}_all_peaks_annotated.tsv --output_tsv {}/{}_rep{}_all_peaks_calculated.tsv --chip_bam {}/{}_chip_merged.normalized.bam\n\n'.format(
                    fold_change_calculator_arg, 
                    cpu_count, 
                    peak_type,
                    peaks_processing_dir, 
                    dataset_name, 
                    (list_counter + 1),
                    peaks_processing_dir, 
                    dataset_name, 
                    (list_counter + 1),
                    results_dir, 
                    dataset_name))


    peaks_processing_script.write('cat {}/{}_rep*all_peaks_concatenated.tsv | cut -f 2,3,4,7 > {}/{}_all_samples_concatenated.tsv\n\n'.format(
        peaks_processing_dir, 
        dataset_name,
        peaks_processing_dir, 
        dataset_name))

    peaks_processing_script.write('cat {}/{}_all_samples_concatenated.tsv'.format(peaks_processing_dir, dataset_name))
    peaks_processing_script.write(' | cut -f 1-3')
    peaks_processing_script.write(r""" | awk '{OFS="\t";print $1,$2,$3,$3-$2}'""")
    peaks_processing_script.write(' | sort -k1,1 -k2,2n') 
    peaks_processing_script.write(' > {}/{}_all_samples_concatenated.bed\n\n'.format(peaks_processing_dir, dataset_name))

    peaks_processing_script.write('bedtools merge -i {}/{}_all_samples_concatenated.bed -c 1,4,4 -o count,mean,collapse > {}/{}_all_samples_merged_intervals.tsv\n\n'.format(
        peaks_processing_dir,
        dataset_name, 
        peaks_processing_dir,
        dataset_name))


    calculated_peaks_list = ['{}/{}_rep{}_all_peaks_calculated.tsv'.format(
        peaks_processing_dir, 
        dataset_name, 
        (list_counter + 1), 
        results_dir, 
        dataset_name) 
        for list_counter in range(len(chip_name))]
        
    calculated_peaks_string = ' '.join(calculated_peaks_list)
    

    peaks_processing_script.write('caras_peak_feature_extractor.py --thread {} --interval_bed {}/{}_all_samples_merged_intervals.tsv --input_tsv {}/{}_all_samples_concatenated.tsv --peaklist_tsv {} --output {} --setname {} --ref {} --repnum {} --clustnum {}{}\n\n'.format(
        cpu_count,
        peaks_processing_dir,
        dataset_name,
        peaks_processing_dir,
        dataset_name,
        calculated_peaks_string,
        peaks_processing_dir,
        dataset_name,
        genome_ref,
        len(chip_name),
        expected_cluster,
        peak_feature_extractor_arg))

        
peaks_processing_script.close() # Closing the script '22_peaks_processing_script.sh'. Flushing the write buffer



####################################################################################################################
### 23_go_annotation_script.sh --- 23_pathway_annotation_script.sh --- 23_go_pathway_annotation_script.sh
####################################################################################################################

# Create a directory "23_supplementary_annotations" for the previous full peak list, 
#   furtherly appended with gene ontology and pathways annotations
supplementary_annotations_dir = '{}/23_supplementary_annotations'.format(output_dir)
if not os.path.exists(supplementary_annotations_dir):
    os.makedirs(supplementary_annotations_dir)

# Create a script "23_go_annotation_script.sh" within, that will iterate through various HOMER 
#   gene ontology database (resulting from HOMER annotatePeaks.pl program call with -go flag toggled on) 
#   and append every relevant gene ontology terms to each peak)
# In: setname_all_peaks_calculated.tsv (input from 22_peaks_processing_script.sh)
# Out: setname_all_peaks_GO_annotated.tsv

go_annotation_script_name = '{}/23_go_annotation_script.sh'.format(supplementary_annotations_dir)

go_annotation_script = open(go_annotation_script_name, 'w')

go_annotation_script.write('#!/bin/bash\n\n')
go_annotation_script.write('set -euxo pipefail\n\n')

# Bash commands to call the home-made script caras_GO_annotator.py that is appending to each peak 
#   all gene ontology terms with the same Entrez ID as the peak 
if analysis_mode == 'bulk':
    go_annotation_script.write('caras_GO_annotator.py --thread {} --input_tsv {}/{}_all_peaks_calculated.tsv --output_tsv {}/{}_all_peaks_GO_annotated.tsv --go_folder {}/{}_gene_ontology\n\n'.format(
        cpu_count, 
        peaks_processing_dir, 
        dataset_name, 
        supplementary_annotations_dir, 
        dataset_name, 
        peaks_processing_dir, 
        dataset_name))

if analysis_mode == 'single_cell':
    for list_counter in range(len(chip_name)):
        go_annotation_script.write('caras_GO_annotator.py --thread {} --input_tsv {}/{}_rep{}_all_peaks_calculated.tsv --output_tsv {}/{}_rep{}_all_peaks_GO_annotated.tsv --go_folder {}/{}_rep{}_gene_ontology\n\n'.format(
            cpu_count, 
            peaks_processing_dir, 
            dataset_name,
            (list_counter + 1), 
            supplementary_annotations_dir, 
            dataset_name,
            (list_counter + 1), 
            peaks_processing_dir, 
            dataset_name,
            (list_counter + 1)))

go_annotation_script.close() # Closing the script '23_go_annotation_script.sh'. Flushing the write buffer



# Create a script "23_pathway_annotation_script.sh" within, that will iterate through various HOMER 
#   database (resulting from HOMER annotatePeaks.pl program call with -go flag toggled on) 
#   and append every relevant pathways, cooccurences, and interactions to each peak)
# In: setname_all_peaks_calculated.tsv (input from 22_peaks_processing_script.sh)
# Out: setname_all_peaks_pathway_annotated.tsv

pathway_annotation_script_name = '{}/23_pathway_annotation_script.sh'.format(supplementary_annotations_dir)

pathway_annotation_script = open(pathway_annotation_script_name, 'w')

pathway_annotation_script.write('#!/bin/bash\n\n')
pathway_annotation_script.write('set -euxo pipefail\n\n')

# Bash commands to call the home-made script caras_pathway_annotator.py that is appending to each peak 
#   all pathways, occurences, and interactions with the same Entrez ID as the peak
if analysis_mode == 'bulk':
    pathway_annotation_script.write('caras_pathway_annotator.py --thread {} --input_tsv {}/{}_all_peaks_calculated.tsv --output_tsv {}/{}_all_peaks_pathway_annotated.tsv --go_folder {}/{}_gene_ontology\n\n'.format(
        cpu_count, 
        peaks_processing_dir, 
        dataset_name, 
        supplementary_annotations_dir, 
        dataset_name, 
        peaks_processing_dir, 
        dataset_name))

if analysis_mode == 'single_cell':
    for list_counter in range(len(chip_name)):
        pathway_annotation_script.write('caras_pathway_annotator.py --thread {} --input_tsv {}/{}_rep{}_all_peaks_calculated.tsv --output_tsv {}/{}_rep{}_all_peaks_pathway_annotated.tsv --go_folder {}/{}_rep{}_gene_ontology\n\n'.format(
            cpu_count, 
            peaks_processing_dir, 
            dataset_name,
            (list_counter + 1), 
            supplementary_annotations_dir, 
            dataset_name,
            (list_counter + 1), 
            peaks_processing_dir, 
            dataset_name,
            (list_counter + 1)))

pathway_annotation_script.close() # Closing the script '23_pathway_annotation_script.sh'. Flushing the write buffer



# Create a script "23_go_pathway_annotation_script.sh" within, that does what both 
#   23_go_annotation_script.sh and 23_pathway_annotation_script.sh do in direct succession.
# When both GO and pathway annotations are toggled on when calling the suite, 
#   the pathway annotation script will automatically process further the GO-annotated peak list file, 
#   resulting in each peak having all annotations from both functions. 
# In: setname_all_peaks_calculated.tsv (input from 22_peaks_processing_script.sh)
# Out: setname_all_peaks_GO_pathway_annotated.tsv

go_pathway_annotation_script_name = '{}/23_go_pathway_annotation_script.sh'.format(supplementary_annotations_dir)

go_pathway_annotation_script = open(go_pathway_annotation_script_name, 'w')

go_pathway_annotation_script.write('#!/bin/bash\n\n')
go_pathway_annotation_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    # Bash commands to call the home-made script caras_GO_annotator.py that is appending to each peak all gene ontology terms with the same Entrez ID as the peak 
    go_pathway_annotation_script.write('caras_GO_annotator.py --thread {} --input_tsv {}/{}_all_peaks_calculated.tsv --output_tsv {}/{}_all_peaks_GO_annotated.tsv --go_folder {}/{}_gene_ontology\n\n'.format(
        cpu_count, 
        peaks_processing_dir, 
        dataset_name, 
        supplementary_annotations_dir, 
        dataset_name, 
        peaks_processing_dir, 
        dataset_name))

    # Bash commands to call the home-made script caras_pathway_annotator.py that is appending to each peak 
    #   all pathways, occurences, and interactions with the same Entrez ID as the peak 
    go_pathway_annotation_script.write('caras_pathway_annotator.py --thread {} --input_tsv {}/{}_all_peaks_GO_annotated.tsv --output_tsv {}/{}_all_peaks_GO_pathway_annotated.tsv --go_folder {}/{}_gene_ontology\n\n'.format(
        cpu_count, 
        supplementary_annotations_dir, 
        dataset_name, 
        supplementary_annotations_dir, 
        dataset_name, 
        peaks_processing_dir, 
        dataset_name))

if analysis_mode == 'single_cell':
    for list_counter in range(len(chip_name)):
        go_pathway_annotation_script.write('caras_GO_annotator.py --thread {} --input_tsv {}/{}_rep{}_all_peaks_calculated.tsv --output_tsv {}/{}_rep{}_all_peaks_GO_annotated.tsv --go_folder {}/{}_rep{}_gene_ontology\n\n'.format(
            cpu_count, 
            peaks_processing_dir, 
            dataset_name,
            (list_counter + 1), 
            supplementary_annotations_dir, 
            dataset_name,
            (list_counter + 1), 
            peaks_processing_dir, 
            dataset_name,
            (list_counter + 1)))
        
        go_pathway_annotation_script.write('caras_pathway_annotator.py --thread {} --input_tsv {}/{}_rep{}_all_peaks_GO_annotated.tsv --output_tsv {}/{}_rep{}_all_peaks_GO_pathway_annotated.tsv --go_folder {}/{}_rep{}_gene_ontology\n\n'.format(
            cpu_count, 
            supplementary_annotations_dir, 
            dataset_name,
            (list_counter + 1), 
            supplementary_annotations_dir, 
            dataset_name,
            (list_counter + 1), 
            peaks_processing_dir, 
            dataset_name,
            (list_counter + 1)))

go_pathway_annotation_script.close() # Closing the script '23_go_pathway_annotation_script.sh'. Flushing the write buffer



####################################################################################################################
### 24_homer_motif_enrichment_script.sh
####################################################################################################################

# Create all the following scripts and parent folders only when processing datasets with narrow peaks. None below exists in datasets with broad peaks.

# Create a directory "24_homer_motif_enrichment" for motif enrichment analysis results of the selected peak set by HOMER
homer_motif_enrichment_dir = '{}/24_homer_motif_enrichment'.format(output_dir)
if not os.path.exists(homer_motif_enrichment_dir + '/logs'):
    os.makedirs(homer_motif_enrichment_dir + '/logs')

# Create a script "24_homer_motif_enrichment_script.sh" within, that reads and filters the processed peak list,
#   selecting the 6_overlap, 1_overlap, or both peak sets, and save them in this folder as inputs for HOMER findMotifsGenome.pl
# This script also calls for motif enrichment analysis program: HOMER findMotifsGenome.pl
#   that performs motif enrichment analysis to find the most potential protein-DNA binding motifs based on:
#       6_overlap peak set if "6_overlap" is given as an argument for --homer_motif flag in CaRAS command line
#       1_overlap peak set if "1_overlap" is given as an argument for --homer_motif_flag in CaRAS command line
#       Both 6_overlap and 1_overlap peak set if "both" is given as an argument for --homer_motif_flag in CaRAS command line
# In: setname_all_peaks_calculated.tsv (input from 22_peaks_processing_script.sh)
# Out: homerResults.html and knownResults.html (all the individual enriched motifs can be accessed through these .html files)


homer_motif_enrichment_1_overlap_script_name = '{}/24_homer_motif_enrichment_1_overlap_script.sh'.format(homer_motif_enrichment_dir)
homer_motif_enrichment_1_overlap_script = open(homer_motif_enrichment_1_overlap_script_name, 'w')
homer_motif_enrichment_1_overlap_script.write('#!/bin/bash\n\n')
homer_motif_enrichment_1_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    homer_motif_enrichment_1_overlap_dir = '{}/24_homer_motif_enrichment/1_overlap_peak_set'.format(output_dir)
    homer_motif_enrichment_1_overlap_script.write('mkdir -p {}\n\n'.format(homer_motif_enrichment_1_overlap_dir))

    # Bash commands to reformat Genrich peak list file (.narrowPeak) into HOMER custom_setting_table format
    homer_motif_enrichment_1_overlap_script.write('cat {}/{}_all_peaks_calculated.tsv'.format(peaks_processing_dir, dataset_name))
    homer_motif_enrichment_1_overlap_script.write(r""" | awk '$7 >= 1'""") # Get only the peaks with peak caller overlaps ($7) with 1 or more overlaps (union) 
    homer_motif_enrichment_1_overlap_script.write(r""" | awk 'BEGIN{fs="\t"}{OFS="\t";print $1,$2,$3,$4,$5}'""") # Get only the unique peak ID ($1), chr ($2), start ($3), end ($4), and strand columns ($5) 
    homer_motif_enrichment_1_overlap_script.write(' > {}/1_overlap_peaks.tsv\n\n'.format(homer_motif_enrichment_1_overlap_dir)) # Save it, as an input peak list for HOMER findMotifsGenome.pl

    # Bash commands to call HOMER findMotifsGenome.pl to perform motif enrichment analysis based on CaRAS 1_overlap peak set. One command, one run for every one replicate.
    homer_motif_enrichment_1_overlap_script.write('findMotifsGenome.pl {}/1_overlap_peaks.tsv {} {} -p {} {} 1> {}/logs/{}.HOMERmotifenrichment_1_overlap.out 2> {}/logs/{}.HOMERmotifenrichment_1_overlap.err\n\n'.format(
        homer_motif_enrichment_1_overlap_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_1_overlap_dir, 
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))
    
if analysis_mode == 'single_cell':
    homer_motif_enrichment_1_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_1_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_1_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_1_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_1_overlap\"\n\n".format(cpu_count, homer_motif_enrichment_dir))

    homer_motif_enrichment_1_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_1_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_1_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_1_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    homer_motif_enrichment_1_overlap_script.write(r""" | awk '\$7 >= 1'""")
    homer_motif_enrichment_1_overlap_script.write(r""" | awk 'BEGIN{FS="'"\t"'"}{OFS="'"\t"'";print \$1,\$2,\$3,\$4,\$5}'""")
    homer_motif_enrichment_1_overlap_script.write(" > {}/{{}}_1_overlap/1_overlap_peaks.tsv\"\n\n".format(homer_motif_enrichment_dir))
                                                  
    homer_motif_enrichment_1_overlap_script.write("find {} -name '1_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"findMotifsGenome.pl {}/{{}}/1_overlap_peaks.tsv {} {}/{{}} -p {} {} 1> {}/logs/{}_{{}}.HOMERmotifenrichment.out 2> {}/logs/{}_{{}}.HOMERmotifenrichment.err\"\n\n".format(
        homer_motif_enrichment_dir,
        homer_motif_enrichment_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_dir,
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))
    
homer_motif_enrichment_1_overlap_script.close() # Closing the script 'homer_motif_enrichment_1_overlap_script.sh'. Flushing the write buffer


homer_motif_enrichment_2_overlap_script_name = '{}/24_homer_motif_enrichment_2_overlap_script.sh'.format(homer_motif_enrichment_dir)
homer_motif_enrichment_2_overlap_script = open(homer_motif_enrichment_2_overlap_script_name, 'w')
homer_motif_enrichment_2_overlap_script.write('#!/bin/bash\n\n')
homer_motif_enrichment_2_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    homer_motif_enrichment_2_overlap_dir = '{}/24_homer_motif_enrichment/2_overlap_peak_set'.format(output_dir)
    homer_motif_enrichment_2_overlap_script.write('mkdir -p {}\n\n'.format(homer_motif_enrichment_2_overlap_dir))

    # Bash commands to reformat Genrich peak list file (.narrowPeak) into HOMER custom_setting_table format
    homer_motif_enrichment_2_overlap_script.write('cat {}/{}_all_peaks_calculated.tsv'.format(peaks_processing_dir, dataset_name))
    homer_motif_enrichment_2_overlap_script.write(r""" | awk '$7 >= 2'""") # Get only the peaks with peak caller overlaps ($7) with 2 or more overlaps
    homer_motif_enrichment_2_overlap_script.write(r""" | awk 'BEGIN{fs="\t"}{OFS="\t";print $1,$2,$3,$4,$5}'""") # Get only the unique peak ID ($1), chr ($2), start ($3), end ($4), and strand columns ($5) 
    homer_motif_enrichment_2_overlap_script.write(' > {}/2_overlap_peaks.tsv\n\n'.format(homer_motif_enrichment_2_overlap_dir)) # Save it, as an input peak list for HOMER findMotifsGenome.pl

    # Bash commands to call HOMER findMotifsGenome.pl to perform motif enrichment analysis based on CaRAS 2_overlap peak set. One command, one run for every one replicate.
    homer_motif_enrichment_2_overlap_script.write('findMotifsGenome.pl {}/2_overlap_peaks.tsv {} {} -p {} {} 1> {}/logs/{}.HOMERmotifenrichment_2_overlap.out 2> {}/logs/{}.HOMERmotifenrichment_2_overlap.err\n\n'.format(
        homer_motif_enrichment_2_overlap_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_2_overlap_dir, 
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))

if analysis_mode == 'single_cell':
    homer_motif_enrichment_2_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_2_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_2_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_2_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_2_overlap\"\n\n".format(cpu_count, homer_motif_enrichment_dir))

    homer_motif_enrichment_2_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_2_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_2_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_2_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    homer_motif_enrichment_2_overlap_script.write(r""" | awk '\$7 >= 2'""")
    homer_motif_enrichment_2_overlap_script.write(r""" | awk 'BEGIN{FS="'"\t"'"}{OFS="'"\t"'";print \$1,\$2,\$3,\$4,\$5}'""")
    homer_motif_enrichment_2_overlap_script.write(" > {}/{{}}_2_overlap/2_overlap_peaks.tsv\"\n\n".format(homer_motif_enrichment_dir))
                                                  
    homer_motif_enrichment_2_overlap_script.write("find {} -name '2_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"findMotifsGenome.pl {}/{{}}/2_overlap_peaks.tsv {} {}/{{}} -p {} {} 1> {}/logs/{}_{{}}.HOMERmotifenrichment.out 2> {}/logs/{}_{{}}.HOMERmotifenrichment.err\"\n\n".format(
        homer_motif_enrichment_dir,
        homer_motif_enrichment_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_dir,
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))
    
homer_motif_enrichment_2_overlap_script.close() # Closing the script 'homer_motif_enrichment_2_overlap_script.sh'. Flushing the write buffer


homer_motif_enrichment_3_overlap_script_name = '{}/24_homer_motif_enrichment_3_overlap_script.sh'.format(homer_motif_enrichment_dir)
homer_motif_enrichment_3_overlap_script = open(homer_motif_enrichment_3_overlap_script_name, 'w')
homer_motif_enrichment_3_overlap_script.write('#!/bin/bash\n\n')
homer_motif_enrichment_3_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    homer_motif_enrichment_3_overlap_dir = '{}/24_homer_motif_enrichment/3_overlap_peak_set'.format(output_dir)
    homer_motif_enrichment_3_overlap_script.write('mkdir -p {}\n\n'.format(homer_motif_enrichment_3_overlap_dir))

    # Bash commands to reformat Genrich peak list file (.narrowPeak) into HOMER custom_setting_table format
    homer_motif_enrichment_3_overlap_script.write('cat {}/{}_all_peaks_calculated.tsv'.format(peaks_processing_dir, dataset_name))
    homer_motif_enrichment_3_overlap_script.write(r""" | awk '$7 >= 3'""") # Get only the peaks with peak caller overlaps ($7) with 3 or more overlaps
    homer_motif_enrichment_3_overlap_script.write(r""" | awk 'BEGIN{fs="\t"}{OFS="\t";print $1,$2,$3,$4,$5}'""") # Get only the unique peak ID ($1), chr ($2), start ($3), end ($4), and strand columns ($5) 
    homer_motif_enrichment_3_overlap_script.write(' > {}/3_overlap_peaks.tsv\n\n'.format(homer_motif_enrichment_3_overlap_dir)) # Save it, as an input peak list for HOMER findMotifsGenome.pl

    # Bash commands to call HOMER findMotifsGenome.pl to perform motif enrichment analysis based on CaRAS 3_overlap peak set. One command, one run for every one replicate.
    homer_motif_enrichment_3_overlap_script.write('findMotifsGenome.pl {}/3_overlap_peaks.tsv {} {} -p {} {} 1> {}/logs/{}.HOMERmotifenrichment_3_overlap.out 2> {}/logs/{}.HOMERmotifenrichment_3_overlap.err\n\n'.format(
        homer_motif_enrichment_3_overlap_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_3_overlap_dir, 
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))

if analysis_mode == 'single_cell':
    homer_motif_enrichment_3_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_3_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_3_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_3_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_3_overlap\"\n\n".format(cpu_count, homer_motif_enrichment_dir))

    homer_motif_enrichment_3_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_3_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_3_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_3_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    homer_motif_enrichment_3_overlap_script.write(r""" | awk '\$7 >= 3'""")
    homer_motif_enrichment_3_overlap_script.write(r""" | awk 'BEGIN{FS="'"\t"'"}{OFS="'"\t"'";print \$1,\$2,\$3,\$4,\$5}'""")
    homer_motif_enrichment_3_overlap_script.write(" > {}/{{}}_3_overlap/3_overlap_peaks.tsv\"\n\n".format(homer_motif_enrichment_dir))
                                                  
    homer_motif_enrichment_3_overlap_script.write("find {} -name '3_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"findMotifsGenome.pl {}/{{}}/3_overlap_peaks.tsv {} {}/{{}} -p {} {} 1> {}/logs/{}_{{}}.HOMERmotifenrichment.out 2> {}/logs/{}_{{}}.HOMERmotifenrichment.err\"\n\n".format(
        homer_motif_enrichment_dir,
        homer_motif_enrichment_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_dir,
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))
    
homer_motif_enrichment_3_overlap_script.close() # Closing the script 'homer_motif_enrichment_3_overlap_script.sh'. Flushing the write buffer


homer_motif_enrichment_4_overlap_script_name = '{}/24_homer_motif_enrichment_4_overlap_script.sh'.format(homer_motif_enrichment_dir)
homer_motif_enrichment_4_overlap_script = open(homer_motif_enrichment_4_overlap_script_name, 'w')
homer_motif_enrichment_4_overlap_script.write('#!/bin/bash\n\n')
homer_motif_enrichment_4_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    homer_motif_enrichment_4_overlap_dir = '{}/24_homer_motif_enrichment/4_overlap_peak_set'.format(output_dir)
    homer_motif_enrichment_4_overlap_script.write('mkdir -p {}\n\n'.format(homer_motif_enrichment_4_overlap_dir))

    # Bash commands to reformat Genrich peak list file (.narrowPeak) into HOMER custom_setting_table format
    homer_motif_enrichment_4_overlap_script.write('cat {}/{}_all_peaks_calculated.tsv'.format(peaks_processing_dir, dataset_name))
    homer_motif_enrichment_4_overlap_script.write(r""" | awk '$7 >= 4'""") # Get only the peaks with peak caller overlaps ($7) with 4 or more overlaps
    homer_motif_enrichment_4_overlap_script.write(r""" | awk 'BEGIN{fs="\t"}{OFS="\t";print $1,$2,$3,$4,$5}'""") # Get only the unique peak ID ($1), chr ($2), start ($3), end ($4), and strand columns ($5) 
    homer_motif_enrichment_4_overlap_script.write(' > {}/4_overlap_peaks.tsv\n\n'.format(homer_motif_enrichment_4_overlap_dir)) # Save it, as an input peak list for HOMER findMotifsGenome.pl

    # Bash commands to call HOMER findMotifsGenome.pl to perform motif enrichment analysis based on CaRAS 4_overlap peak set. One command, one run for every one replicate.
    homer_motif_enrichment_4_overlap_script.write('findMotifsGenome.pl {}/4_overlap_peaks.tsv {} {} -p {} {} 1> {}/logs/{}.HOMERmotifenrichment_4_overlap.out 2> {}/logs/{}.HOMERmotifenrichment_4_overlap.err\n\n'.format(
        homer_motif_enrichment_4_overlap_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_4_overlap_dir, 
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))

if analysis_mode == 'single_cell':
    homer_motif_enrichment_4_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_4_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_4_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_4_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_4_overlap\"\n\n".format(cpu_count, homer_motif_enrichment_dir))

    homer_motif_enrichment_4_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_4_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_4_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_4_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    homer_motif_enrichment_4_overlap_script.write(r""" | awk '\$7 >= 4'""")
    homer_motif_enrichment_4_overlap_script.write(r""" | awk 'BEGIN{FS="'"\t"'"}{OFS="'"\t"'";print \$1,\$2,\$3,\$4,\$5}'""")
    homer_motif_enrichment_4_overlap_script.write(" > {}/{{}}_4_overlap/4_overlap_peaks.tsv\"\n\n".format(homer_motif_enrichment_dir))
                                                  
    homer_motif_enrichment_4_overlap_script.write("find {} -name '4_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"findMotifsGenome.pl {}/{{}}/4_overlap_peaks.tsv {} {}/{{}} -p {} {} 1> {}/logs/{}_{{}}.HOMERmotifenrichment.out 2> {}/logs/{}_{{}}.HOMERmotifenrichment.err\"\n\n".format(
        homer_motif_enrichment_dir,
        homer_motif_enrichment_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_dir,
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))
    
homer_motif_enrichment_4_overlap_script.close() # Closing the script 'homer_motif_enrichment_4_overlap_script.sh'. Flushing the write buffer


homer_motif_enrichment_5_overlap_script_name = '{}/24_homer_motif_enrichment_5_overlap_script.sh'.format(homer_motif_enrichment_dir)
homer_motif_enrichment_5_overlap_script = open(homer_motif_enrichment_5_overlap_script_name, 'w')
homer_motif_enrichment_5_overlap_script.write('#!/bin/bash\n\n')
homer_motif_enrichment_5_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    homer_motif_enrichment_5_overlap_dir = '{}/24_homer_motif_enrichment/5_overlap_peak_set'.format(output_dir)
    homer_motif_enrichment_5_overlap_script.write('mkdir -p {}\n\n'.format(homer_motif_enrichment_5_overlap_dir))

    # Bash commands to reformat Genrich peak list file (.narrowPeak) into HOMER custom_setting_table format
    homer_motif_enrichment_5_overlap_script.write('cat {}/{}_all_peaks_calculated.tsv'.format(peaks_processing_dir, dataset_name))
    homer_motif_enrichment_5_overlap_script.write(r""" | awk '$7 >= 5'""") # Get only the peaks with peak caller overlaps ($7) with 5 or more overlaps
    homer_motif_enrichment_5_overlap_script.write(r""" | awk 'BEGIN{fs="\t"}{OFS="\t";print $1,$2,$3,$4,$5}'""") # Get only the unique peak ID ($1), chr ($2), start ($3), end ($4), and strand columns ($5) 
    homer_motif_enrichment_5_overlap_script.write(' > {}/5_overlap_peaks.tsv\n\n'.format(homer_motif_enrichment_5_overlap_dir)) # Save it, as an input peak list for HOMER findMotifsGenome.pl

    # Bash commands to call HOMER findMotifsGenome.pl to perform motif enrichment analysis based on CaRAS 5_overlap peak set. One command, one run for every one replicate.
    homer_motif_enrichment_5_overlap_script.write('findMotifsGenome.pl {}/5_overlap_peaks.tsv {} {} -p {} {} 1> {}/logs/{}.HOMERmotifenrichment_5_overlap.out 2> {}/logs/{}.HOMERmotifenrichment_5_overlap.err\n\n'.format(
        homer_motif_enrichment_5_overlap_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_5_overlap_dir, 
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))

if analysis_mode == 'single_cell':
    # Bash commands to call HOMER findMotifsGenome.pl to perform motif enrichment analysis based on CaRAS 5_overlap peak set. One command, one run for every one replicate.
    homer_motif_enrichment_5_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_5_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_5_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_5_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_5_overlap\"\n\n".format(cpu_count, homer_motif_enrichment_dir))

    homer_motif_enrichment_5_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_5_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_5_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_5_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    homer_motif_enrichment_5_overlap_script.write(r""" | awk '\$7 >= 5'""")
    homer_motif_enrichment_5_overlap_script.write(r""" | awk 'BEGIN{FS="'"\t"'"}{OFS="'"\t"'";print \$1,\$2,\$3,\$4,\$5}'""")
    homer_motif_enrichment_5_overlap_script.write(" > {}/{{}}_5_overlap/5_overlap_peaks.tsv\"\n\n".format(homer_motif_enrichment_dir))
                                                  
    homer_motif_enrichment_5_overlap_script.write("find {} -name '5_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"findMotifsGenome.pl {}/{{}}/5_overlap_peaks.tsv {} {}/{{}} -p {} {} 1> {}/logs/{}_{{}}.HOMERmotifenrichment.out 2> {}/logs/{}_{{}}.HOMERmotifenrichment.err\"\n\n".format(
        homer_motif_enrichment_dir,
        homer_motif_enrichment_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_dir,
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))
    
homer_motif_enrichment_5_overlap_script.close() # Closing the script 'homer_motif_enrichment_5_overlap_script.sh'. Flushing the write buffer


homer_motif_enrichment_6_overlap_script_name = '{}/24_homer_motif_enrichment_6_overlap_script.sh'.format(homer_motif_enrichment_dir)
homer_motif_enrichment_6_overlap_script = open(homer_motif_enrichment_6_overlap_script_name, 'w')
homer_motif_enrichment_6_overlap_script.write('#!/bin/bash\n\n')
homer_motif_enrichment_6_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    homer_motif_enrichment_6_overlap_dir = '{}/24_homer_motif_enrichment/6_overlap_peak_set'.format(output_dir)
    homer_motif_enrichment_6_overlap_script.write('mkdir -p {}\n\n'.format(homer_motif_enrichment_6_overlap_dir))

    # Bash commands to reformat Genrich peak list file (.narrowPeak) into HOMER custom_setting_table format
    homer_motif_enrichment_6_overlap_script.write('cat {}/{}_all_peaks_calculated.tsv'.format(peaks_processing_dir, dataset_name))
    homer_motif_enrichment_6_overlap_script.write(r""" | awk '$7 >= 6'""") # Get only the 6_overlap peaks (peaks with peak caller overlaps ($7) = 6) 
    homer_motif_enrichment_6_overlap_script.write(r""" | awk 'BEGIN{fs="\t"}{OFS="\t";print $1,$2,$3,$4,$5}'""") # Get only the unique peak ID ($1), chr ($2), start ($3), end ($4), and strand columns ($5) 
    homer_motif_enrichment_6_overlap_script.write(' > {}/6_overlap_peaks.tsv\n\n'.format(homer_motif_enrichment_6_overlap_dir)) # Save it, as an input peak list for HOMER findMotifsGenome.pl

    # Bash commands to call HOMER findMotifsGenome.pl to perform motif enrichment analysis based on CaRAS 6_overlap peak set. One command, one run for every one replicate.
    homer_motif_enrichment_6_overlap_script.write('findMotifsGenome.pl {}/6_overlap_peaks.tsv {} {} -p {} {} 1> {}/logs/{}.HOMERmotifenrichment_6_overlap.out 2> {}/logs/{}.HOMERmotifenrichment_6_overlap.err\n\n'.format(
        homer_motif_enrichment_6_overlap_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_6_overlap_dir, 
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))

if analysis_mode == 'single_cell':
    # Bash commands to call HOMER findMotifsGenome.pl to perform motif enrichment analysis based on CaRAS 6_overlap peak set. One command, one run for every one replicate.
    homer_motif_enrichment_6_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_6_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_6_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_6_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_6_overlap\"\n\n".format(cpu_count, homer_motif_enrichment_dir))

    homer_motif_enrichment_6_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    homer_motif_enrichment_6_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    homer_motif_enrichment_6_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    homer_motif_enrichment_6_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    homer_motif_enrichment_6_overlap_script.write(r""" | awk '\$7 >= 6'""")
    homer_motif_enrichment_6_overlap_script.write(r""" | awk 'BEGIN{FS="'"\t"'"}{OFS="'"\t"'";print \$1,\$2,\$3,\$4,\$5}'""")
    homer_motif_enrichment_6_overlap_script.write(" > {}/{{}}_6_overlap/6_overlap_peaks.tsv\"\n\n".format(homer_motif_enrichment_dir))
                                                  
    homer_motif_enrichment_6_overlap_script.write("find {} -name '6_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"findMotifsGenome.pl {}/{{}}/6_overlap_peaks.tsv {} {}/{{}} -p {} {} 1> {}/logs/{}_{{}}.HOMERmotifenrichment.out 2> {}/logs/{}_{{}}.HOMERmotifenrichment.err\"\n\n".format(
        homer_motif_enrichment_dir,
        homer_motif_enrichment_dir,
        '{}/bwa/{}.fa'.format(genome_dir, genome_ref) if genome_ref not in complete_genome else genome_ref,
        homer_motif_enrichment_dir,
        cpu_count,
        homer_findMotifsGenome_arg,
        homer_motif_enrichment_dir,
        dataset_name,
        homer_motif_enrichment_dir,
        dataset_name))
    
homer_motif_enrichment_6_overlap_script.close() # Closing the script 'homer_motif_enrichment_6_overlap_script.sh'. Flushing the write buffer



####################################################################################################################
### 25_meme_motif_enrichment_script.sh
####################################################################################################################

# Create all the following scripts and parent folders only when processing datasets with narrow peaks. None below exists in datasets with broad peaks.

# Create a directory "25_meme_motif_enrichment" for motif enrichment analysis results of the selected peak set by MEME
meme_motif_enrichment_dir = '{}/25_meme_motif_enrichment'.format(output_dir)
if not os.path.exists(meme_motif_enrichment_dir + '/logs'):
    os.makedirs(meme_motif_enrichment_dir + '/logs')

# Create a script "24_homer_motif_enrichment_script.sh" within, that reads and filters the processed peak list,
#   selecting the 6_overlap, 1_overlap, or both peak sets, and save them in this folder as inputs for HOMER findMotifsGenome.pl
# This script also calls for motif enrichment analysis program: HOMER findMotifsGenome.pl
#   that performs motif enrichment analysis to find the most potential protein-DNA binding motifs based on:

# Create a script "25_meme_motif_enrichment_script.sh" within, that generates target and background sequences
#   based on the peak coordinates in the processed peak list, and save them in this folder as inputs for meme-chip
# This script also calls for motif enrichment analysis program: meme-chip
#   that performs motif enrichment analysis to find the most potential protein-DNA binding motifs based on:
#       6_overlap peak set if "6_overlap" is given as an argument for --meme_motif flag in CaRAS command line
#       1_overlap peak set if "1_overlap" is given as an argument for --meme_motif_flag in CaRAS command line
#       Both 6_overlap and 1_overlap peak set if "both" is given as an argument for --meme_motif_flag in CaRAS command line
# In: setname_all_peaks_calculated.tsv (input from 22_peaks_processing_script.sh)
# Out: meme-chip.html (all the enriched motifs and modular analysis results can be accessed through this .html file)

if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
    meme_sequence_extractor_background_arg = ' --background'
else:
    meme_sequence_extractor_background_arg = ''

if analysis_mode == 'bulk':
    # Put _rep[#] behind the file or folder name when there are going to be multi-replicated results from multi-replicated meme-chip runs
    if len(chip_name) == 1:
        fasta_suffix_list = ['']
    elif len(chip_name) > 1:
        if force_merge == 1:
            fasta_suffix_list = ['']
        elif force_merge == 0:
            fasta_suffix_list = ['_rep{}'.format(list_counter + 1) for list_counter in range(len(chip_name))]

if analysis_mode == 'single_cell':
    if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
        meme_chip_background_arg = ' -neg {}/{{}}/background.fa'.format(meme_motif_enrichment_dir)
    else:
        meme_chip_background_arg = ''

meme_motif_enrichment_1_overlap_script_name = '{}/25_meme_motif_enrichment_1_overlap_script.sh'.format(meme_motif_enrichment_dir)
meme_motif_enrichment_1_overlap_script = open(meme_motif_enrichment_1_overlap_script_name, 'w')
meme_motif_enrichment_1_overlap_script.write('#!/bin/bash\n\n')
meme_motif_enrichment_1_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':        
    meme_motif_enrichment_1_overlap_dir = '{}/25_meme_motif_enrichment/1_overlap_peak_set'.format(output_dir)
    meme_motif_enrichment_1_overlap_script.write('mkdir -p {}\n\n'.format(meme_motif_enrichment_1_overlap_dir))

    meme_motif_enrichment_1_overlap_script.write('caras_meme_sequence_extractor.py --input {}/{}_all_peaks_calculated.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {} --ref {}{} --filter 1\n\n'.format(
        peaks_processing_dir,
        dataset_name,
        genome_dir,
        genome_dir,
        meme_motif_enrichment_1_overlap_dir,
        genome_ref,
        meme_sequence_extractor_background_arg))

    # Bash commands to call meme-chip to perform motif enrichment analysis based on CaRAS 1_overlap peak set. One command, one run for every one replicate.
    for list_counter in range(len(fasta_suffix_list)):

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            meme_chip_background_arg = ' -neg {}/background{}.fa'.format(meme_motif_enrichment_1_overlap_dir, fasta_suffix_list[list_counter])
        else:
            meme_chip_background_arg = ''

        meme_motif_enrichment_1_overlap_script.write('meme-chip -oc {}/motif_analysis{}{} {} {}/target{}.fa 1> {}/logs/{}.MEMEmotifenrichment_1_overlap{}.out 2> {}/logs/{}.MEMEmotifenrichment_1_overlap{}.err &\n'.format(
            meme_motif_enrichment_1_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_chip_background_arg,
            meme_chip_arg,
            meme_motif_enrichment_1_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(fasta_suffix_list) - 1):
            meme_motif_enrichment_1_overlap_script.write('\nwait\n\n')

if analysis_mode == 'single_cell':
    meme_motif_enrichment_1_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_1_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_1_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_1_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_1_overlap\"\n\n".format(cpu_count, meme_motif_enrichment_dir))

    meme_motif_enrichment_1_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_1_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_1_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_1_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    meme_motif_enrichment_1_overlap_script.write(r""" | awk '\$7 >= 6'""")
    meme_motif_enrichment_1_overlap_script.write(" > {}/{{}}_1_overlap/1_overlap_peaks.tsv\"\n\n".format(meme_motif_enrichment_dir))
                                                  
    meme_motif_enrichment_1_overlap_script.write("find {} -name '1_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"caras_meme_sequence_extractor.py --input {}/{{}}/1_overlap_peaks.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {}/{{}} --ref {}{}\"\n\n".format(
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        genome_dir,
        genome_dir,
        meme_motif_enrichment_dir,
        genome_ref,
        meme_sequence_extractor_background_arg))

    meme_motif_enrichment_1_overlap_script.write("find {} -name '1_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j{} \"meme-chip -oc {}/{{}}{}{} {}/{{}}/target.fa 1> {}/logs/{}_{{}}_1_overlap.MEMEmotifenrichment.out 2> {}/logs/{}_{{}}_1_overlap.MEMEmotifenrichment.err\"\n".format(
        meme_motif_enrichment_dir,
        cpu_count,
        meme_motif_enrichment_dir,
        meme_chip_background_arg,
        meme_chip_arg,
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        dataset_name,
        meme_motif_enrichment_dir,
        dataset_name))
    
meme_motif_enrichment_1_overlap_script.close() # Closing the script 'meme_motif_enrichment_1_overlap_script.sh'. Flushing the write buffer


meme_motif_enrichment_2_overlap_script_name = '{}/25_meme_motif_enrichment_2_overlap_script.sh'.format(meme_motif_enrichment_dir)
meme_motif_enrichment_2_overlap_script = open(meme_motif_enrichment_2_overlap_script_name, 'w')
meme_motif_enrichment_2_overlap_script.write('#!/bin/bash\n\n')
meme_motif_enrichment_2_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':        
    meme_motif_enrichment_2_overlap_dir = '{}/25_meme_motif_enrichment/2_overlap_peak_set'.format(output_dir)
    meme_motif_enrichment_2_overlap_script.write('mkdir -p {}\n\n'.format(meme_motif_enrichment_2_overlap_dir))

    meme_motif_enrichment_2_overlap_script.write('caras_meme_sequence_extractor.py --input {}/{}_all_peaks_calculated.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {} --ref {}{} --filter 2\n\n'.format(peaks_processing_dir,
                                    dataset_name,
                                    genome_dir,
                                    genome_dir,
                                    meme_motif_enrichment_2_overlap_dir,
                                    genome_ref,
                                    meme_sequence_extractor_background_arg))

    # Bash commands to call meme-chip to perform motif enrichment analysis based on CaRAS 2_overlap peak set. One command, one run for every one replicate.
    for list_counter in range(len(fasta_suffix_list)):

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            meme_chip_background_arg = ' -neg {}/background{}.fa'.format(meme_motif_enrichment_2_overlap_dir, fasta_suffix_list[list_counter])
        else:
            meme_chip_background_arg = ''

        meme_motif_enrichment_2_overlap_script.write('meme-chip -oc {}/motif_analysis{}{} {} {}/target{}.fa 1> {}/logs/{}.MEMEmotifenrichment_2_overlap{}.out 2> {}/logs/{}.MEMEmotifenrichment_2_overlap{}.err &\n'.format(
            meme_motif_enrichment_2_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_chip_background_arg,
            meme_chip_arg,
            meme_motif_enrichment_2_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(fasta_suffix_list) - 1):
            meme_motif_enrichment_2_overlap_script.write('\nwait\n\n')

if analysis_mode == 'single_cell':
    meme_motif_enrichment_2_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_2_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_2_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_2_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_2_overlap\"\n\n".format(cpu_count, meme_motif_enrichment_dir))

    meme_motif_enrichment_2_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_2_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_2_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_2_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    meme_motif_enrichment_2_overlap_script.write(r""" | awk '\$7 >= 6'""")
    meme_motif_enrichment_2_overlap_script.write(" > {}/{{}}_2_overlap/2_overlap_peaks.tsv\"\n\n".format(meme_motif_enrichment_dir))
                                                  
    meme_motif_enrichment_2_overlap_script.write("find {} -name '2_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"caras_meme_sequence_extractor.py --input {}/{{}}/2_overlap_peaks.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {}/{{}} --ref {}{}\"\n\n".format(
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        genome_dir,
        genome_dir,
        meme_motif_enrichment_dir,
        genome_ref,
        meme_sequence_extractor_background_arg))

    meme_motif_enrichment_2_overlap_script.write("find {} -name '2_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j{} \"meme-chip -oc {}/{{}}{}{} {}/{{}}/target.fa 1> {}/logs/{}_{{}}_2_overlap.MEMEmotifenrichment.out 2> {}/logs/{}_{{}}_2_overlap.MEMEmotifenrichment.err\"\n".format(
        meme_motif_enrichment_dir,
        cpu_count,
        meme_motif_enrichment_dir,
        meme_chip_background_arg,
        meme_chip_arg,
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        dataset_name,
        meme_motif_enrichment_dir,
        dataset_name))
    
meme_motif_enrichment_2_overlap_script.close() # Closing the script 'meme_motif_enrichment_2_overlap_script.sh'. Flushing the write buffer


meme_motif_enrichment_3_overlap_script_name = '{}/25_meme_motif_enrichment_3_overlap_script.sh'.format(meme_motif_enrichment_dir)
meme_motif_enrichment_3_overlap_script = open(meme_motif_enrichment_3_overlap_script_name, 'w')
meme_motif_enrichment_3_overlap_script.write('#!/bin/bash\n\n')
meme_motif_enrichment_3_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    meme_motif_enrichment_3_overlap_dir = '{}/25_meme_motif_enrichment/3_overlap_peak_set'.format(output_dir)
    meme_motif_enrichment_3_overlap_script.write('mkdir -p {}\n\n'.format(meme_motif_enrichment_3_overlap_dir))

    meme_motif_enrichment_3_overlap_script.write('caras_meme_sequence_extractor.py --input {}/{}_all_peaks_calculated.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {} --ref {}{} --filter 3\n\n'.format(peaks_processing_dir,
                                    dataset_name,
                                    genome_dir,
                                    genome_dir,
                                    meme_motif_enrichment_3_overlap_dir,
                                    genome_ref,
                                    meme_sequence_extractor_background_arg))

    # Bash commands to call meme-chip to perform motif enrichment analysis based on CaRAS 3_overlap peak set. One command, one run for every one replicate.
    for list_counter in range(len(fasta_suffix_list)):

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            meme_chip_background_arg = ' -neg {}/background{}.fa'.format(meme_motif_enrichment_3_overlap_dir, fasta_suffix_list[list_counter])
        else:
            meme_chip_background_arg = ''

        meme_motif_enrichment_3_overlap_script.write('meme-chip -oc {}/motif_analysis{}{} {} {}/target{}.fa 1> {}/logs/{}.MEMEmotifenrichment_3_overlap{}.out 2> {}/logs/{}.MEMEmotifenrichment_3_overlap{}.err &\n'.format(
            meme_motif_enrichment_3_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_chip_background_arg,
            meme_chip_arg,
            meme_motif_enrichment_3_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(fasta_suffix_list) - 1):
            meme_motif_enrichment_3_overlap_script.write('\nwait\n\n')

if analysis_mode == 'single_cell':
    meme_motif_enrichment_3_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_3_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_3_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_3_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_3_overlap\"\n\n".format(cpu_count, meme_motif_enrichment_dir))

    meme_motif_enrichment_3_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_3_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_3_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_3_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    meme_motif_enrichment_3_overlap_script.write(r""" | awk '\$7 >= 6'""")
    meme_motif_enrichment_3_overlap_script.write(" > {}/{{}}_3_overlap/3_overlap_peaks.tsv\"\n\n".format(meme_motif_enrichment_dir))
                                                  
    meme_motif_enrichment_3_overlap_script.write("find {} -name '3_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"caras_meme_sequence_extractor.py --input {}/{{}}/3_overlap_peaks.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {}/{{}} --ref {}{}\"\n\n".format(
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        genome_dir,
        genome_dir,
        meme_motif_enrichment_dir,
        genome_ref,
        meme_sequence_extractor_background_arg))

    meme_motif_enrichment_3_overlap_script.write("find {} -name '3_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j{} \"meme-chip -oc {}/{{}}{}{} {}/{{}}/target.fa 1> {}/logs/{}_{{}}_3_overlap.MEMEmotifenrichment.out 2> {}/logs/{}_{{}}_3_overlap.MEMEmotifenrichment.err\"\n".format(
        meme_motif_enrichment_dir,
        cpu_count,
        meme_motif_enrichment_dir,
        meme_chip_background_arg,
        meme_chip_arg,
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        dataset_name,
        meme_motif_enrichment_dir,
        dataset_name))
    
meme_motif_enrichment_3_overlap_script.close() # Closing the script 'meme_motif_enrichment_3_overlap_script.sh'. Flushing the write buffer


meme_motif_enrichment_4_overlap_script_name = '{}/25_meme_motif_enrichment_4_overlap_script.sh'.format(meme_motif_enrichment_dir)
meme_motif_enrichment_4_overlap_script = open(meme_motif_enrichment_4_overlap_script_name, 'w')
meme_motif_enrichment_4_overlap_script.write('#!/bin/bash\n\n')
meme_motif_enrichment_4_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    meme_motif_enrichment_4_overlap_dir = '{}/25_meme_motif_enrichment/4_overlap_peak_set'.format(output_dir)
    meme_motif_enrichment_4_overlap_script.write('mkdir -p {}\n\n'.format(meme_motif_enrichment_4_overlap_dir))

    meme_motif_enrichment_4_overlap_script.write('caras_meme_sequence_extractor.py --input {}/{}_all_peaks_calculated.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {} --ref {}{} --filter 4\n\n'.format(peaks_processing_dir,
                                    dataset_name,
                                    genome_dir,
                                    genome_dir,
                                    meme_motif_enrichment_4_overlap_dir,
                                    genome_ref,
                                    meme_sequence_extractor_background_arg))

    # Bash commands to call meme-chip to perform motif enrichment analysis based on CaRAS 4_overlap peak set. One command, one run for every one replicate.
    for list_counter in range(len(fasta_suffix_list)):

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            meme_chip_background_arg = ' -neg {}/background{}.fa'.format(meme_motif_enrichment_4_overlap_dir, fasta_suffix_list[list_counter])
        else:
            meme_chip_background_arg = ''

        meme_motif_enrichment_4_overlap_script.write('meme-chip -oc {}/motif_analysis{}{} {} {}/target{}.fa 1> {}/logs/{}.MEMEmotifenrichment_4_overlap{}.out 2> {}/logs/{}.MEMEmotifenrichment_4_overlap{}.err &\n'.format(
            meme_motif_enrichment_4_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_chip_background_arg,
            meme_chip_arg,
            meme_motif_enrichment_4_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(fasta_suffix_list) - 1):
            meme_motif_enrichment_4_overlap_script.write('\nwait\n\n')

if analysis_mode == 'single_cell':
    meme_motif_enrichment_4_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_4_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_4_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_4_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_4_overlap\"\n\n".format(cpu_count, meme_motif_enrichment_dir))

    meme_motif_enrichment_4_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_4_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_4_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_4_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    meme_motif_enrichment_4_overlap_script.write(r""" | awk '\$7 >= 6'""")
    meme_motif_enrichment_4_overlap_script.write(" > {}/{{}}_4_overlap/4_overlap_peaks.tsv\"\n\n".format(meme_motif_enrichment_dir))
                                                  
    meme_motif_enrichment_4_overlap_script.write("find {} -name '4_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"caras_meme_sequence_extractor.py --input {}/{{}}/4_overlap_peaks.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {}/{{}} --ref {}{}\"\n\n".format(
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        genome_dir,
        genome_dir,
        meme_motif_enrichment_dir,
        genome_ref,
        meme_sequence_extractor_background_arg))

    meme_motif_enrichment_4_overlap_script.write("find {} -name '4_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j{} \"meme-chip -oc {}/{{}}{}{} {}/{{}}/target.fa 1> {}/logs/{}_{{}}_4_overlap.MEMEmotifenrichment.out 2> {}/logs/{}_{{}}_4_overlap.MEMEmotifenrichment.err\"\n".format(
        meme_motif_enrichment_dir,
        cpu_count,
        meme_motif_enrichment_dir,
        meme_chip_background_arg,
        meme_chip_arg,
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        dataset_name,
        meme_motif_enrichment_dir,
        dataset_name))
    
meme_motif_enrichment_4_overlap_script.close() # Closing the script 'meme_motif_enrichment_4_overlap_script.sh'. Flushing the write buffer


meme_motif_enrichment_5_overlap_script_name = '{}/25_meme_motif_enrichment_5_overlap_script.sh'.format(meme_motif_enrichment_dir)
meme_motif_enrichment_5_overlap_script = open(meme_motif_enrichment_5_overlap_script_name, 'w')
meme_motif_enrichment_5_overlap_script.write('#!/bin/bash\n\n')
meme_motif_enrichment_5_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    meme_motif_enrichment_5_overlap_dir = '{}/25_meme_motif_enrichment/5_overlap_peak_set'.format(output_dir)
    meme_motif_enrichment_5_overlap_script.write('mkdir -p {}\n\n'.format(meme_motif_enrichment_5_overlap_dir))

    meme_motif_enrichment_5_overlap_script.write('caras_meme_sequence_extractor.py --input {}/{}_all_peaks_calculated.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {} --ref {}{} --filter 5\n\n'.format(peaks_processing_dir,
                                    dataset_name,
                                    genome_dir,
                                    genome_dir,
                                    meme_motif_enrichment_5_overlap_dir,
                                    genome_ref,
                                    meme_sequence_extractor_background_arg))

    # Bash commands to call meme-chip to perform motif enrichment analysis based on CaRAS 5_overlap peak set. One command, one run for every one replicate.
    for list_counter in range(len(fasta_suffix_list)):

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            meme_chip_background_arg = ' -neg {}/background{}.fa'.format(meme_motif_enrichment_5_overlap_dir, fasta_suffix_list[list_counter])
        else:
            meme_chip_background_arg = ''

        meme_motif_enrichment_5_overlap_script.write('meme-chip -oc {}/motif_analysis{}{} {} {}/target{}.fa 1> {}/logs/{}.MEMEmotifenrichment_5_overlap{}.out 2> {}/logs/{}.MEMEmotifenrichment_5_overlap{}.err &\n'.format(
            meme_motif_enrichment_5_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_chip_background_arg,
            meme_chip_arg,
            meme_motif_enrichment_5_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(fasta_suffix_list) - 1):
            meme_motif_enrichment_5_overlap_script.write('\nwait\n\n')

if analysis_mode == 'single_cell':
    meme_motif_enrichment_5_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_5_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_5_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_5_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_5_overlap\"\n\n".format(cpu_count, meme_motif_enrichment_dir))

    meme_motif_enrichment_5_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_5_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_5_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_5_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    meme_motif_enrichment_5_overlap_script.write(r""" | awk '\$7 >= 6'""")
    meme_motif_enrichment_5_overlap_script.write(" > {}/{{}}_5_overlap/5_overlap_peaks.tsv\"\n\n".format(meme_motif_enrichment_dir))
                                                  
    meme_motif_enrichment_5_overlap_script.write("find {} -name '5_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"caras_meme_sequence_extractor.py --input {}/{{}}/5_overlap_peaks.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {}/{{}} --ref {}{}\"\n\n".format(
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        genome_dir,
        genome_dir,
        meme_motif_enrichment_dir,
        genome_ref,
        meme_sequence_extractor_background_arg))

    meme_motif_enrichment_5_overlap_script.write("find {} -name '5_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j{} \"meme-chip -oc {}/{{}}{}{} {}/{{}}/target.fa 1> {}/logs/{}_{{}}_5_overlap.MEMEmotifenrichment.out 2> {}/logs/{}_{{}}_5_overlap.MEMEmotifenrichment.err\"\n".format(
        meme_motif_enrichment_dir,
        cpu_count,
        meme_motif_enrichment_dir,
        meme_chip_background_arg,
        meme_chip_arg,
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        dataset_name,
        meme_motif_enrichment_dir,
        dataset_name))
    
meme_motif_enrichment_5_overlap_script.close() # Closing the script 'meme_motif_enrichment_5_overlap_script.sh'. Flushing the write buffer


meme_motif_enrichment_6_overlap_script_name = '{}/25_meme_motif_enrichment_6_overlap_script.sh'.format(meme_motif_enrichment_dir)
meme_motif_enrichment_6_overlap_script = open(meme_motif_enrichment_6_overlap_script_name, 'w')
meme_motif_enrichment_6_overlap_script.write('#!/bin/bash\n\n')
meme_motif_enrichment_6_overlap_script.write('set -euxo pipefail\n\n')

if analysis_mode == 'bulk':
    meme_motif_enrichment_6_overlap_dir = '{}/25_meme_motif_enrichment/6_overlap_peak_set'.format(output_dir)
    meme_motif_enrichment_6_overlap_script.write('mkdir -p {}\n\n'.format(meme_motif_enrichment_6_overlap_dir))

    meme_motif_enrichment_6_overlap_script.write('caras_meme_sequence_extractor.py --input {}/{}_all_peaks_calculated.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {} --ref {}{} --filter 6\n\n'.format(peaks_processing_dir,
                                    dataset_name,
                                    genome_dir,
                                    genome_dir,
                                    meme_motif_enrichment_6_overlap_dir,
                                    genome_ref,
                                    meme_sequence_extractor_background_arg))

    # Bash commands to call meme-chip to perform motif enrichment analysis based on CaRAS 6_overlap peak set. One command, one run for every one replicate.
    for list_counter in range(len(fasta_suffix_list)):

        if ctrl_r1_sample_list or len(ctrl_r1_sample_list) != 0:
            meme_chip_background_arg = ' -neg {}/background{}.fa'.format(meme_motif_enrichment_6_overlap_dir, fasta_suffix_list[list_counter])
        else:
            meme_chip_background_arg = ''

        meme_motif_enrichment_6_overlap_script.write('meme-chip -oc {}/motif_analysis{}{} {} {}/target{}.fa 1> {}/logs/{}.MEMEmotifenrichment_6_overlap{}.out 2> {}/logs/{}.MEMEmotifenrichment_6_overlap{}.err &\n'.format(
            meme_motif_enrichment_6_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_chip_background_arg,
            meme_chip_arg,
            meme_motif_enrichment_6_overlap_dir,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter],
            meme_motif_enrichment_dir,
            dataset_name,
            fasta_suffix_list[list_counter]))

        if ((list_counter + 1) % cpu_count) == 0 or list_counter == (len(fasta_suffix_list) - 1):
            meme_motif_enrichment_6_overlap_script.write('\nwait\n\n')

if analysis_mode == 'single_cell':
    meme_motif_enrichment_6_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_6_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_6_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_6_overlap_script.write(" | parallel -j{} \"mkdir -p {}/{{}}_6_overlap\"\n\n".format(cpu_count, meme_motif_enrichment_dir))

    meme_motif_enrichment_6_overlap_script.write("find {} -name '*all_peaks_clustered*' -type f -exec basename {{}} \\;".format(peaks_processing_dir))
    meme_motif_enrichment_6_overlap_script.write(" | sed 's/_all_peaks_clustered\.tsv//'")
    meme_motif_enrichment_6_overlap_script.write(" | sed 's/{}_//'".format(dataset_name))
    meme_motif_enrichment_6_overlap_script.write(" | parallel -j{} \"cat {}/{}_{{}}_all_peaks_clustered.tsv".format(cpu_count, peaks_processing_dir, dataset_name))
    meme_motif_enrichment_6_overlap_script.write(r""" | awk '\$7 >= 6'""")
    meme_motif_enrichment_6_overlap_script.write(" > {}/{{}}_6_overlap/6_overlap_peaks.tsv\"\n\n".format(meme_motif_enrichment_dir))
                                                  
    meme_motif_enrichment_6_overlap_script.write("find {} -name '6_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j1 \"caras_meme_sequence_extractor.py --input {}/{{}}/6_overlap_peaks.tsv --fastadir {}/bwa --chrsizedir {}/chrom.sizes --outputdir {}/{{}} --ref {}{}\"\n\n".format(
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        genome_dir,
        genome_dir,
        meme_motif_enrichment_dir,
        genome_ref,
        meme_sequence_extractor_background_arg))

    meme_motif_enrichment_6_overlap_script.write("find {} -name '6_overlap_peaks.tsv\' -type f -exec dirname {{}} \\; | xargs -L1 basename | parallel -j{} \"meme-chip -oc {}/{{}}{}{} {}/{{}}/target.fa 1> {}/logs/{}_{{}}_6_overlap.MEMEmotifenrichment.out 2> {}/logs/{}_{{}}_6_overlap.MEMEmotifenrichment.err\"\n".format(
        meme_motif_enrichment_dir,
        cpu_count,
        meme_motif_enrichment_dir,
        meme_chip_background_arg,
        meme_chip_arg,
        meme_motif_enrichment_dir,
        meme_motif_enrichment_dir,
        dataset_name,
        meme_motif_enrichment_dir,
        dataset_name))
    
meme_motif_enrichment_6_overlap_script.close() # Closing the script 'meme_motif_enrichment_6_overlap_script.sh'. Flushing the write buffer



####################################################################################################################
### MASTER_script.sh
####################################################################################################################

# Create a master script "MASTER_script.sh" in the output directory, 
#   containing bash commands that calls every individual suite-generated bash script, sequentially
# Running the master script means running through the whole suite, processing the raw fq.gz down to the very end
# Any individual suite-generated bash script can be viewed in its own folder, 
#   have its commands modified, and executed directly in bash terminal

master_script_name = '{}/MASTER_script.sh'.format(output_dir)
master_script = open(master_script_name, 'w')

master_script.write('#!/bin/bash\n\n')
master_script.write('set -euxo pipefail\n\n')
master_script.write('ulimit -n 2000\n\n')

master_script.write('echo Creating a renamed copy of raw data in the output directory.\n')
master_script.write('{}\n\n'.format(raw_data_script_name))

master_script.write('echo Running quality check on the raw sequencing reads\n')
master_script.write('{}\n\n'.format(raw_reads_quality_control_script_name))

master_script.write('echo Removing optical and tile-edge duplicates with bbmap clumpify\n')
master_script.write('{}\n\n'.format(deduplicating_script_name))

master_script.write('echo Trimming the left sequencing adapter from each read with bbmap bbduk\n')
master_script.write('{}\n\n'.format(adapter_trimming_script_name))

master_script.write('echo Trimming low quality bases and drop low quality reads using trimmomatic\n')
master_script.write('{}\n\n'.format(quality_trimming_script_name))

master_script.write('echo Running quality check on the preprocessed sequencing reads.\n')
master_script.write('{}\n\n'.format(preprocessed_reads_quality_control_script_name))

master_script.write('echo Aligning the paired-ends sequenced reads to {} using bwa mem algorithm\n'.format(genome_ref))
master_script.write('{}\n\n'.format(bwa_mem_aligning_script_name))

master_script.write('echo Filtering the aligned reads according to MAPQ score\n')
master_script.write('{}\n\n'.format(mapq_filtering_script_name))

master_script.write('echo Sorting, indexing and generating coverages of the filtered aligned reads\n')
master_script.write('{}\n\n'.format(results_script_name))

master_script.write('echo Running quality check on the aligned sequencing reads.\n')
master_script.write('{}\n\n'.format(aligned_reads_quality_control_script_name))

master_script.write('echo Running peak calling with MACS2\n')
master_script.write('{}\n\n'.format(macs2_peak_calling_script_name))

master_script.write('echo Running peak calling with GEM\n')
master_script.write('{}\n\n'.format(gem_peak_calling_script_name))

master_script.write('echo Running peak calling with HOMER\n')
master_script.write('{}\n\n'.format(homer_peak_calling_script_name))

master_script.write('echo Running peak calling with Genrich\n')
master_script.write('{}\n\n'.format(genrich_peak_calling_script_name))

master_script.write('echo Running peak calling with SEACR\n')
master_script.write('{}\n\n'.format(seacr_peak_calling_script_name))

master_script.write('echo Running peak calling with SICER2\n')
master_script.write('{}\n\n'.format(sicer2_peak_calling_script_name))

master_script.write('echo Collecting, reformating and merging all peak lists with HOMER mergePeaks\n')
master_script.write('{}\n\n'.format(peaks_merging_script_name))

master_script.write('echo Concatenating, annotating, and calculating peak stats: peak caller overlaps, fold change, and motif hits\n')
master_script.write('{}\n\n'.format(peaks_processing_script_name))



if analysis_mode == 'bulk' and genome_ref in complete_genome:
    
    if args.goann and not args.pathann:
        master_script.write('echo Annotating every individual peak with all relevant GO terms from HOMER annotatePeaks -go feature\n')
        master_script.write('{}\n\n'.format(go_annotation_script_name))

    if args.pathann and not args.goann:
        master_script.write('echo Annotating every individual peak with all relevant pathways, cooccurences, and interactions from HOMER annotatePeaks -go feature\n')
        master_script.write('{}\n\n'.format(pathway_annotation_script_name))

    # If the --goann and --pathann switchs are both toggled on, 
    #   pathway annotation will be executed after the GO annotation, 
    #   which makes the pathway annotation columns are located to the right 
    #   of the GO annotation columns in the resulting table. 
    if args.goann and args.pathann:
        master_script.write('echo Annotating every individual peak with all relevant GO terms, pathways, cooccurences, and interactions from HOMER annotatePeaks -go feature\n')
        master_script.write('{}\n\n'.format(go_pathway_annotation_script_name))


if args.homer_motif: # Motif enrichment analyses are only performed on datasets with narrow peaks
    
    if 'union' in homer_motif_peakset or '1' in homer_motif_peakset: # If user wants analysis on the 1_overlap peak set (everything)
        master_script.write('echo Performing motif enrichment analysis on the 1_overlap peak set with HOMER findMotifsGenome.pl\n')
        master_script.write('{}\n\n'.format(homer_motif_enrichment_1_overlap_script_name))
    
    if '2' in homer_motif_peakset: # If user wants analysis on the 2_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 2_overlap peak set with HOMER findMotifsGenome.pl\n')
        master_script.write('{}\n\n'.format(homer_motif_enrichment_2_overlap_script_name))
        
    if '3' in homer_motif_peakset: # If user wants analysis on the 3_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 3_overlap peak set with HOMER findMotifsGenome.pl\n')
        master_script.write('{}\n\n'.format(homer_motif_enrichment_3_overlap_script_name))

    if '4' in homer_motif_peakset: # If user wants analysis on the 4_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 4_overlap peak set with HOMER findMotifsGenome.pl\n')
        master_script.write('{}\n\n'.format(homer_motif_enrichment_4_overlap_script_name))

    if '5' in homer_motif_peakset: # If user wants analysis on the 5_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 5_overlap peak set with HOMER findMotifsGenome.pl\n')
        master_script.write('{}\n\n'.format(homer_motif_enrichment_5_overlap_script_name))

    if 'consensus' in homer_motif_peakset or '6' in homer_motif_peakset: # If user wants analysis on the 6_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 6_overlap peak set with HOMER findMotifsGenome.pl\n')
        master_script.write('{}\n\n'.format(homer_motif_enrichment_6_overlap_script_name))
    

if args.meme_motif: # Motif enrichment analyses are only performed on datasets with narrow peaks
    
    if 'union' in meme_motif_peakset or '1' in meme_motif_peakset: # If user wants analysis on the 1_overlap peak set (everything)
        master_script.write('echo Performing motif enrichment analysis on the 1_overlap peak set with MEME meme-chip &\n')
        master_script.write('{} &\n\n'.format(meme_motif_enrichment_1_overlap_script_name))
    
    if '2' in meme_motif_peakset: # If user wants analysis on the 2_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 2_overlap peak set with MEME meme-chip &\n')
        master_script.write('{} &\n\n'.format(meme_motif_enrichment_2_overlap_script_name))
        
    if '3' in meme_motif_peakset: # If user wants analysis on the 3_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 3_overlap peak set with MEME meme-chip &\n')
        master_script.write('{} &\n\n'.format(meme_motif_enrichment_3_overlap_script_name))

    if '4' in meme_motif_peakset: # If user wants analysis on the 4_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 4_overlap peak set with MEME meme-chip &\n')
        master_script.write('{} &\n\n'.format(meme_motif_enrichment_4_overlap_script_name))

    if '5' in meme_motif_peakset: # If user wants analysis on the 5_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 5_overlap peak set with MEME meme-chip &\n')
        master_script.write('{} &\n\n'.format(meme_motif_enrichment_5_overlap_script_name))

    if 'consensus' in meme_motif_peakset or '6' in meme_motif_peakset: # If user wants analysis on the 6_overlap peak set
        master_script.write('echo Performing motif enrichment analysis on the 6_overlap peak set with MEME meme-chip &\n')
        master_script.write('{} &\n\n'.format(meme_motif_enrichment_6_overlap_script_name))

    master_script.write('wait\n\n')

master_script.close() # Closing the script 'MASTER_script.sh'. Flushing the write buffer



####################################################################################################################
### SCRIPTS MARKING AS EXECUTABLE
####################################################################################################################

print('Marking all generated bash scripts as executable')
# Bash commands to mark all the suite-generated bash scripts above as executable
subprocess.run('chmod +x ' + raw_data_script_name, shell = True)
subprocess.run('chmod +x ' + raw_reads_quality_control_script_name, shell = True)
subprocess.run('chmod +x ' + deduplicating_script_name, shell = True)
subprocess.run('chmod +x ' + adapter_trimming_script_name, shell = True)
subprocess.run('chmod +x ' + quality_trimming_script_name, shell = True)
subprocess.run('chmod +x ' + preprocessed_reads_quality_control_script_name, shell = True)
subprocess.run('chmod +x ' + bwa_mem_aligning_script_name, shell = True)
subprocess.run('chmod +x ' + mapq_filtering_script_name, shell = True)
subprocess.run('chmod +x ' + results_script_name, shell = True)
subprocess.run('chmod +x ' + aligned_reads_quality_control_script_name, shell = True)
subprocess.run('chmod +x ' + macs2_peak_calling_script_name, shell = True)
subprocess.run('chmod +x ' + gem_peak_calling_script_name, shell = True)
subprocess.run('chmod +x ' + homer_peak_calling_script_name, shell = True)
subprocess.run('chmod +x ' + genrich_peak_calling_script_name, shell = True)
subprocess.run('chmod +x ' + seacr_peak_calling_script_name, shell = True)
subprocess.run('chmod +x ' + sicer2_peak_calling_script_name, shell = True)

subprocess.run('chmod +x ' + peaks_merging_script_name, shell = True)
subprocess.run('chmod +x ' + peaks_processing_script_name, shell = True)
subprocess.run('chmod +x ' + go_annotation_script_name, shell = True)
subprocess.run('chmod +x ' + pathway_annotation_script_name, shell = True)
subprocess.run('chmod +x ' + go_pathway_annotation_script_name, shell = True)

# Only mark all motif enrichment analysis scripts in datasets with narrow peaks
# Motif enrichment analysis scripts and their parent folders are not generated in datasets with broad peaks
subprocess.run('chmod +x ' + homer_motif_enrichment_1_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + homer_motif_enrichment_2_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + homer_motif_enrichment_3_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + homer_motif_enrichment_4_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + homer_motif_enrichment_5_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + homer_motif_enrichment_6_overlap_script_name, shell = True)

subprocess.run('chmod +x ' + meme_motif_enrichment_1_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + meme_motif_enrichment_2_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + meme_motif_enrichment_3_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + meme_motif_enrichment_4_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + meme_motif_enrichment_5_overlap_script_name, shell = True)
subprocess.run('chmod +x ' + meme_motif_enrichment_6_overlap_script_name, shell = True)

subprocess.run('chmod +x ' + master_script_name, shell = True)



####################################################################################################################
### THE BIG RED BUTTON
####################################################################################################################

# Automatically executes all scripts in order if --run switch is used
if args.run:
    subprocess.run(master_script_name, shell = True)
elif not args.run:
    print('\nNOTICE ::: The suite now can be started manually by running MASTER_script.sh in {} ::: NOTICE\n'.format(output_dir))


