<p align="center">
<img src="https://raw.githubusercontent.com/JSuryatenggara/CaRAS/storage/images/CaRAS_logo.png" width="1000">
</p>
<p align="center">
	<b>C</b>UT-<b>a</b>nd-<b>R</b>UN <b>A</b>nalysis <b>S</b>uite
</p>
<p align="center">
	<b>Inherently Simple CUT&RUN Analysis</b>
</p>
<p align="center">
	Developed by Jeremiah Suryatenggara and Mahmoud A. Bassal
</p>

## Introduction
CaRAS (CUT-and-RUN Analysis Suite) (https://github.com/JSuryatenggara/CaRAS) was conceived as an all-in-one analysis pipeline specifically tailored and optimized for CUT&RUN/Tag. CaRAS is a unified interface and wrapper tailored for the new generation of ChIP-Seq: CUT&RUN and CUT&Tag analysis (for the sake of conciseness, both techniques will be referred to as simply CnR). CaRAS builds on the same programming framework as ChIP-AP (https://github.com/JSuryatenggara/ChIP-AP), sharing a number of similar methodologies and working characteristics but diverges in that it tailors specific analysis steps specifically for CnR analyses that are yet unpublished for CnR aimed workflows.


CaRAS is the first analysis pipeline and framework which <b>generates normalized bam (nBAM) files</b>, which are aligned bam files that are normalized to the spiked-in or leftover DNA, which is critical in CnR analysis. Facilitated by the normalization of aligned reads, ChIP-Seq peak callers were able to be optimized to work properly for CnR peak detection. This leads to the next key feature which is the utilization of <b>five modified ChIP-Seq peak callers in combination with a CUT&RUN dedicated peak caller (SEACR)</b>, giving a total of 6 peak callers to establish consensus with. CaRAS is also equipped with <b>single cell analysis capabilities (beta)</b>, where single cell samples are clustered based on their enriched gene ontology terms (GO) or pathways, and downstream analyses will be performed on each cluster to enable capturing the states / phases / differentiation of the cells.


Just like ChIP-AP, our <b>C</b>ut-<b>a</b>nd-<b>R</b>UN <b>A</b>nalysis <b>S</b>uite (CaRAS) utilizes multiple peak callers, which enables filtering in/out results from every peak caller used in this pipeline (MACS2, GEM, HOMER, Genrich, SEACR, and SICER2) to:

1. Selectively obtain high confidence peaks based on overlaps by different peak callers.

2. Avoid the caveats of depending on a single peak caller in datasets with different peak characteristics and signal-to-noise ratios.

3. Freely choose the “sweet spot” between the consensus peaks set (intersection of detected peaks from all peak callers) and the total peaks set (the union of all detected peaks), that answers one’s biological questions which may be supported by additional evidence from other experiments. This can be done from the output without re-processing data.

Additionally, just like ChIP-AP, CaRAS is a fully automated ChIP-seq data processing and analysis pipeline with these general features: 

1. For input, it takes unaligned sequencing reads in FASTQ format (extension: .fastq / .fq / .fastq.gz / .fq.gz) or previously aligned reads in BAM format (extension: .bam).

2. CaRAS is capable of processing and analyze multiple sample replicates.

3. CaRAS is a complete, integrated workflow which performs all analysis steps (QC, cleanup, alignment, peak-calling, pathway analysis) and outputs a final integrated peak set.

4. The output of CaRAS is a detected peaks list, annotated with information on associated gene names, IDs, functions, gene ontology, and related pathways.

<br>

## System requirements
<b>OS</b> – Linux (Ubuntu-variants 20.04+ tested)

<b>CPU</b> – (minimum) Quad-Core Intel/AMD CPU, (recommended) Octa-Core Intel/AMD CPU.

<b>RAM</b> – (minimum) 16Gb, (recommended) 32Gb+

<b>Storage (SSD/HDD)</b> – Installation alone takes ~60Gb.  Roughly, an additional 50-150Gb space is required for processing and analyzing samples.

<br>

## Installing CaRAS on local machine or laptop
CaRAS has been designed to be as simple as possible to run for end-users. We have prepared full installation guides and wiki’s located here on our github for your operating system (OS) of choice (look at the top you will see a “Wiki” tab). These are complete step-by-step guides with accompanying screenshots to walk users through every step of the installation processes. For fresh installations, we have not encountered issues when following our installation instructions accurately. If you have previous installations of configurations, some tweaking may be necessary to get everything to work right – just be aware.

For advanced users, the installation instructions are as follows:

1 - Download the caras_installation zip in our github repo and unzip in the folder of your choosing.

    wget https://github.com/JSuryatenggara/CaRAS/raw/main/chipap_installation.zip

2 - Download the pre-generated, required genome files from our dropbox and unzip in the same folder as the CaRAS scripts – i.e.

    <path>…/caras_installation/caras_scripts/

3 - Anaconda3 is required for CaRAS to install and run. If you have already installed Anaconda3 on your machine, proceed to the next step. If you need to setup Anaconda3, download the required command-line installer from Anaconda using wget. We recommend isntalling anaconda with default settings and when prompted at the end of the installation, initialize the base conda environment.

Before continuing, you must close and reopen the terminal to initialize conda.

4 - In a new terminal window, you should now see (base) written to the left of the command prompt, this tells you conda has been properly initialized. If it hasn’t been initialized, refer to the conda documentation, section "Should I add Anaconda to the macOS or Linux PATH?"

5 - Next navigate to the caras_installation folder. Check that the caras_installer.py script is executable. If not, type the command:

    chmod +x ./caras_installer.py

and then execute it. The installer will then ask you a couple of questions regarding installing in its own environment (highly recommended) and then proceed to install everything and all dependencies. Depending on your CPU and internet speeds this will take a while as it needs to download ~10-15Gb of genome data for HOMER.

<br>

## Quick start – For Command line User Only
CaRAS is capable of handling multiple sample replicates in a single bulk analysis run. It is also capable of handling an imbalanced number of sample replicates (i.e., 3 ChIP, 2 Controls). It does so by merging each corresponding sample type following alignment and these merged files will be used for downstream processing, yielding one set of analysis results.

In single-cell analysis mode, multiple samples are considered as unique samples instead of replicates of one another. These samples will later be clustered into several groups based on their enriched GO terms or pathways profile similarity. Each group will be considered as unique cell population and undergo its own downstream processing, yielding at the end of the run multiple sets of analysis results.

For peak calling, peaks are called based on ChIP signal over control.

### Example : To run bulk analysis of single-end unaligned reads with default settings (sample organism human hg38, no reads normalization):
    caras.py  --analysis bulk \\
              --mode single \\
              --chipR1 [chip fastq replicate1] [chip fastq replicate2] … \\
              --ctrlR1 [control fastq replicate1] [control fastq replicate2] … \\
              --genome [path to genome folder] \\
              --output [path to output save folder] \\
              --setname [dataset name]

### Example : To run single-cell analysis of paired-end unaligned reads with default settings (sample organism mm10, reads normalization based on reads mapped to E. coli K12):
    caras.py  --analysis single_cell \\
              --mode paired \\
              --chipR1 [chip fastq replicate1, first read] [chip fastq replicate2, first read] … \\
              --chipR2 [chip fastq replicate1, second read] [chip fastq replicate2, second read] … \\
              --ctrlR1 [control fastq replicate1, first read] [control fastq replicate2, first read] … \\
              --ctrlR2 [control fastq replicate1, second read] [control fastq replicate2, second read] … \\
              --genome [path to output save folder] \\
              --ref mm10 \\
              --normref eColiK12 \\ 
              --output [path to output save folder] \\
              --setname [dataset name]

### Example : To run bulk analysis of single/paired-end aligned reads with default settings:
    caras.py  --analysis bulk \\
              --mode single / paired \\ 
              --chipR1 [chip bam replicate1] [chip bam replicate2] … \\
              --ctrlR1 [control bam replicate1] [control bam replicate2] … \\
              --genome [path to genome folder] \\
              --output [path to output save folder] \\
              --setname [dataset name]

<br>

## Usage notes and Command Line Flags / Parameters
### Required Arguments

| <nobr><b>Argument/Flag</b></nobr> | <nobr><b>Possible Values</b></nobr> | <nobr><b>Detailed Explanation</b></nobr> |
|-|-|-|
| <nobr>--analysis</nobr> | <nobr>bulk / single_cell</nobr> | Bulk or single-cell analysis type. Bulk analysis is the typical protocol for NGS-based protein-DNA genome-wide binding study, where every sample typically comes from numerous cells and will be regarded as a biological replicate. As such, bulk analysis will generate one set of results. On the other hand, single-cell analysis is the advanced protocol applicable only for a few select cutting edge experiments, where every sample comes from one single cell and will be regarded as potentially unique until clustered based on enriched GO terms or pathways profile similarity. As such, single-cell analysis will generate multiple sets of results, depending on how many distinct cell populations (clusters) are recognized in all samples (i.e., n clusters = n sets). |
| <nobr>--mode</nobr> | <nobr>single / paired</nobr> | Single-end or paired-end sequencing analysis. If a paired-end run, files will typically be labelled ending in *_R1 and *_R2 before the file extension. If these labels aren’t present then likely, you have single-ended sequencing data and select the “single” option. |
| <nobr>--genome</nobr> | <nobr>[directory]</nobr> | Your genome folder directory. Requires full path, not relative path. This is the folder where the pre-computed genome alignment and processing files are saved. These genomes can be downloaded from ( https://www.dropbox.com/s/ioqd3hwdahh9xon/genomes.zip ) or you can compute your own (in our GitHub guides directory, look for “How to Generate Your Own CaRAS Genome Files”). |
| <nobr>--output</nobr> | <nobr>[directory]</nobr> | Your desired output folder. Requires full path, not relative path. |
| <nobr>--setname</nobr> | <nobr>[text]</nobr> | The prefix to label output and intermediate files (no space allowed). CaRAS will rename all processed data files to have this prefix. |

### Optional Arguments

| <nobr><b>Argument/Flag</b></nobr> | <nobr><b>Possible Values</b></nobr> | <nobr><b>Detailed Explanation</b></nobr> |
|-|-|-|
| <nobr>--chipR1</nobr> | <nobr>[rep1 rep2 …]</nobr> | Your ChIP datasets: ordered by replicate, separated by space. Best to include full path rather than relative paths. |
| <nobr>--chipR2</nobr> | [rep1 rep2 …] | [ Paired-end Only ] Your ChIP datasets second read: ordered by replicate, separated by space. Best to include full path rather than relative paths. |
| <nobr>--ctrlR1</nobr> | <nobr>[rep1 rep2 …]</nobr> | Your control datasets: ordered by replicate, separated by space. Best to include full path rather than relative paths. |
| <nobr>--ctrlR2</nobr> | <nobr>[rep1 rep2 …]</nobr> | [ Paired-end Only ] Your control datasets second read: ordered by replicate, separated by space. Best to include full path rather than relative paths. |
| --sample<br>_table | <nobr>[file]</nobr> | Rather than including the input sample file names and paths on the command line, one can make a sample table containing the same information and input this instead. The sample-table is a 4*n-sized table (with n = number of replicates) containing the absolute paths to each of your ChIP and control replicates (See below for more information regarding this file and its layout).
When this option is used, this table will disregard any assigned values to --chipR1, --chipR2, --ctrlR1, and --ctrlR2. |
| --custom<br>_setting<br>_table | <nobr>[file]</nobr> | [ For Advanced Users ONLY ] The settings-table allows you to fine-tune the behaviour of every program as part of CaRAS. Hence, can be disastrous if you mess up! If you are unsure even a little about what you’re doing, then please stick to default settings – this goes even for bioinformaticians.<br> This .tsv file is a 2*n-sized table (with n = number of replicates) containing custom arguments for every program as part of CaRAS. The default settings table is provided in the genome folder. You can COPY this file and make changes as necessary. To use your customs settings table, provide full path to updated .tsv file. See below for more information regarding this file and its layout. |
| <nobr>--ref</nobr> | hg19 / hg38 /mm9 / mm10 / mm39 / dm6 / sacCer3 | Your sample genome reference build. Default is hg38 (human). The genomes listed to the left are provided and pre-calculated. To add your own custom genome to CaRAS, look for “How to Generate Your Own CaRAS Genome Files” in our GitHub guides directory. |
| <nobr>--norm_ref</nobr> | <nobr>eColiK12 / sacCer3</nobr> | Your read normalizer organism genome reference build. Choose according to the organism which DNA is spiked-in or carried-over to the samples during the experiment. CaRAS will attempt to align the reads to the chosen genome and use the outcome as the basis for inter-samples read normalization. <b>Important</b>: CaRAS does not perform read normalization as the correct read normalization organism genome has to be chosen carefully by users based on the CUT&RUN assay. To add your own custom genome to CaRAS, look for “How to Generate Your Own CaRAS Genome Files” in our GitHub guides directory. |
| <nobr>--motif</nobr> | <nobr>[file]</nobr> | Your predicted/known motif file, in HOMER matrix format (.motif). If provided, once peaks are called, HOMER motif discovery will be run on the total called peak set for this motif. |
| <nobr>--fcmerge</nobr> | | This flag will force fold change analysis to be computed based on merged replicates instead of on each replicate separately. This flag will be automatically activated when the number of ChIP samples is unequal to control samples. |
| --clustnum | <nobr>[integer]</nobr> | Expected number of clusters to be formed by all single-cell samples. When not provided, CaRAS will attempt to calculate the optimal number of clusters by itself (auto). Only relevant for single-cell analysis and will be ignored in bulk analysis. Use this if you preemptively know or want your cells to be divided into a certain number of groups. CaRAS will then generate the best threshold to divide between these groups based on their mapped 2D distribution based on GO terms or pathways profile similarity. |
| <nobr>--goann</nobr> | | This flag will instruct to annotate peaks with all relevant GO terms as provided by HOMER. |
| <nobr>--pathann</nobr> | | This flag will instruct to annotate peaks with all relevant pathway and interaction enrichment terms as provided by HOMER. |
| <nobr>--deltemp</nobr> | | This flag will instruct to delete large intermediary files right after they are not going to be used for further processes (eg intermediary fq files. This option will save a significant amount of space by the end of your analysis, so recommended. |
| <nobr>--thread</nobr> | <nobr>[integer]</nobr> | Maximum number of processes to use. Default is half the maximum available on your system so as to not choke it during the run. If running on a laptop or low-thread count CPU, best to push this up to maximum number threads available -1 – but this will significantly slow your laptop if attempting other tasks while CaRAS is running. |
| <nobr>--run</nobr> | | Use to immediately run the suite by running the master script. This is the big red button. Use it at your own risk.<br>When not used, the generated master script (MASTER_script.sh) in the output folder can be run manually by user. |
| <nobr>--homer_motif</nobr> | 1 / 2 / 3 / 4 / 5 / 6 | This flag will instruct HOMER (findMotifsGenome.pl) to perform motif enrichment analysis at the end of the pipeline run. The analysis will be performed on selected peaks called by 1/2/3/4/5/6 peak callers. 1 being the same as the union peak set (all called peaks), and 6 being the same as the consensus peak set (intersection of all peak caller sets). Inputting multiple arguments (space-separated) instructs HOMER to perform multiple analysis runs on selected sets. |
| <nobr>--meme_motif</nobr> | 1 / 2 / 3 / 4 / 5 / 6 | This flag will instruct MEME suite (meme-chip) to perform motif enrichment analysis at the end of the pipeline run. The analysis will be performed on selected peaks called by 1/2/3/4/5/6 peak callers. 1 being the same as the union peak set (all called peaks), and 6 being the same as the consensus peak set (intersection of all peak caller sets). Inputting multiple arguments (space-separated) instructs MEME suite to perform multiple analysis runs on selected sets. |

<br>

## CaRAS Graphical Overview
<img src="https://raw.githubusercontent.com/JSuryatenggara/CaRAS/storage/images/caras_workflow_guide.png" width="1000">

For detailed explanations of all the steps and methodologies used throughout CaRAS refer to our documentation (https://github.com/JSuryatenggara/CaRAS/wiki/CaRAS-Guide)

<br>

## Pipeline Outputs
### Final Analysis Table (including supplementary annotations)
The table below shows the contents of [setname]_all_peaks_go_pathway_annotated.tsv. Smaller sized and less verbose variants of this table are saved in the output folder with suffixes: concatenated, annotated, and calculated (see Source in the table below)

<table>
    <thead>
        <tr>
            <th>Column #</th>
            <th>Peak Attribute</th>
            <th>Source</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td><nobr><b>Column 1 (A)</td>
            <td><nobr>Peak ID (unique peak ID)</td>
            <td rowspan=6>
                <b>Pipeline script:</b><br>22_peaks_ processing _script.sh
                <br><b>Called program:</b><br>cat (Bash)
                <br><b>Output file:</b><br>[setname]_all_peaks_concatenated.tsv
                <br><b>Output folder:</b><br>22_peaks_processing
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 2 (B)</td>
            <td><nobr>Chr (chromosome)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 3 (C)</td>
            <td><nobr>Start (peak start coordinate)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 4 (D)</td>
            <td><nobr>End (peak end coordinate)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 5 (E)</td>
            <td><nobr>Strand (on which peak is found)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 6 (F)</td>
            <td><nobr>Peak Caller Combination</td>
        </tr>
        <tr>
            <td><nobr><b>Column 7 (G)</td>
            <td><nobr>Peak Caller Overlaps</td>
            <td rowspan=8>
                <b>Pipeline script:</b><br>22_peaks_processing_script.sh
                <br><b>Called script:</b><br>fold_change_calculator.py
                <br><b>Output file:</b><br>[setname]_all_peaks_calculated.tsv
                <br><b>Output folder:</b><br>22_peaks_processing
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 8 (H)</td>
            <td><nobr>ChIP Tag Count</td>
        </tr>
        <tr>
            <td><nobr><b>Column 9 (I)</td>
            <td><nobr>Control Tag Count</td>
        </tr>
        <tr>
            <td><nobr><b>Column 10 (J)</td>
            <td><nobr>Fold Change</td>
        </tr>
        <tr>
            <td><nobr><b>Column 11 (K)</td>
            <td><nobr>Peak Center</td>
        </tr>
        <tr>
            <td><nobr><b>Column 12 (L)</td>
            <td><nobr>Number of Motifs</td>
        </tr>
        <tr>
            <td><nobr><b>Column 13 (M)</td>
            <td><nobr>negLog10_IDR</td>
        </tr>
        <tr>
            <td><nobr><b>Column 14 (N)</td>
            <td><nobr>IDR</td>
        </tr>
        <tr>
            <td><nobr><b>Column 15 (O)</td>
            <td><nobr>Annotation</td>
            <td rowspan=14>
                <b>Pipeline script:</b><br>22_peaks_processing_script.sh
                <br><b>Called program:</b><br>HOMER annotatePeaks
                <br><b>Output file:</b><br>[setname]_all_peaks_annotated.tsv
                <br><b>Output folder:</b><br>22_peaks_processing
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 16 (P)</td>
            <td><nobr>Detailed Annotation</td>
        </tr>
        <tr>
            <td><nobr><b>Column 17 (Q)</td>
            <td><nobr>Distance to TSS</td>
        </tr>
        <tr>
            <td><nobr><b>Column 18 (R)</td>
            <td><nobr>Nearest PromoterID</td>
        </tr>
        <tr>
            <td><nobr><b>Column 19 (S)</td>
            <td><nobr>Entrez ID</td>
        </tr>
        <tr>
            <td><nobr><b>Column 20 (T)</td>
            <td><nobr>Nearest Unigene</td>
        </tr>
        <tr>
            <td><nobr><b>Column 21 (U)</td>
            <td><nobr>Nearest Refseq</td>
        </tr>
        <tr>
            <td><nobr><b>Column 22 (V)</td>
            <td><nobr>Nearest Ensembl</td>
        </tr>
        <tr>
            <td><nobr><b>Column 23 (W)</td>
            <td><nobr>Gene Name</td>
        </tr>
        <tr>
            <td><nobr><b>Column 24 (X)</td>
            <td><nobr>Gene Alias</td>
        </tr>
        <tr>
            <td><nobr><b>Column 25 (Y)</td>
            <td><nobr>Gene Description</td>
        </tr>
        <tr>
            <td><nobr><b>Column 26 (Z)</td>
            <td><nobr>Gene Type</td>
        </tr>
        <tr>
            <td><nobr><b>Column 27 (AA)</td>
            <td><nobr>CpG%</td>
        </tr>
        <tr>
            <td><nobr><b>Column 28 (AB)</td>
            <td><nobr>GC%</td>
        </tr>
        <tr>
            <td><nobr><b>Column 29 (AC)</td>
            <td><nobr>Biological Process</td>
            <td rowspan=3>
                <b>Pipeline script:</b><br>23_go_annotation_script.sh
                <br><b>Called script:</b><br>GO_annotator.py
                <br><b>Output file:</b><br>[setname]_all_peaks_go_annotated.tsv
                <br><b>Output folder:</b><br>23_supplementary_annotations
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 30 (AD)</td>
            <td><nobr>Molecular Function</td>
        </tr>
        <tr>
            <td><nobr><b>Column 31 (AE)</td>
            <td><nobr>Cellular Component</td>
        </tr>
        <tr>
            <td><nobr><b>Column 32 (AF)</td>
            <td><nobr>Interaction with Common Protein</td>
            <td rowspan=8>
                <b>Pipeline script:</b><br>23_pathway_annotation_script.sh
                <br><b>Called script:</b><br>pathway_annotator.py
                <br><b>Output file:</b><br>[setname]_all_peaks_pathway_annotated.tsv
                <br><b>Output folder:</b><br>23_supplementary_annotations
            </td>
        </tr>
        <tr>
            <td><nobr><b>Column 33 (AG)</td>
            <td><nobr>Somatic Mutations (COSMIC)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 34 (AH)</td>
            <td><nobr>Pathway (KEGG)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 35 (AI)</td>
            <td><nobr>Pathway (BIOCYC)</td>
        </tr>
        <tr>
           <td><nobr><b>Column 36 (AJ)</td>
            <td><nobr>Pathway (pathwayInteractionDB)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 37 (AK)</td>
            <td><nobr>Pathway (REACTOME)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 38 (AL)</td>
            <td><nobr>Pathway (SMPDB)</td>
        </tr>
        <tr>
            <td><nobr><b>Column 39 (AM)</td>
            <td><nobr>Pathway (Wikipathways)</td>
        </tr>
    </tbody>
</table>

<br>

### Pipeline Run Info
This file summarizes the assignment of the files (IP sample or control, read 1 or 2; replicate number) and the file name conversion for every unaligned or aligned sequencing reads to be processed. Each line tells the user what the original files have been renamed into. Check this file if you suspect the order of samples were incorrectly entered (e.g., swapped chip with control).

    Chromatin IP dataset replicate 1, 1st read : Original filename = a.fastq --> New filename = setname_chip_rep1_R1.fq.gz

    Chromatin IP dataset replicate 2, 1st read : Original filename = b.fastq --> New filename = setname _chip_rep2_R1.fq.gz  

    Chromatin IP dataset replicate 1, 2nd read : Original filename = c.fastq --> New filename = setname _chip_rep1_R2.fq.gz  

    Chromatin IP dataset replicate 2, 2nd read : Original filename = d.fastq --> New filename = setname _chip_rep2_R2.fq.gz  

    Control dataset replicate 1, 1st read : Original filename = e.fastq --> New filename = setname _ctrl_rep1_R1.fq.gz  

    Control dataset replicate 2, 1st read : Original filename = f.fastq --> New filename = setname _ctrl_rep2_R1.fq.gz  

    Control dataset replicate 1, 2nd read : Original filename = g.fastq --> New filename = setname _ctrl_rep1_R2.fq.gz  

    Control dataset replicate 2, 2nd read : Original filename = h.fastq --> New filename = setname _ctrl_rep2_R2.fq.gz

<br>

### Pipeline Run Command
Contains the input command line that was used to call the pipeline in a text file: [filename]_command_line.txt in the output save folder. This is useful for documentation of the run, and for re-running of the pipeline after a run failure or some tweaking if need be.

    caras.py --mode paired --ref [genome_build] --genome [path_to_computed_genome_folders] --output [full_path_to_output_save_folder] --setname [dataset name] --sample_table [path_to_sample_table_file] --custom_setting_table [path_to_setting_table_file].tsv --motif [path_to_known_motif_file] --fcmerge --goann --pathann --deltemp --thread [#_of_threads_to_use] --run


<br>

### Sample Table
Contains the full path of each input ChIP and control sample in the pipeline run in a tab-separated value file: [filename]_sample_table.tsv in the output save folder in CaRAS sample table format. This is useful for documentation of the run, and for re-running of the pipeline after a run failure or some tweaking if need be. Below is an example of sample table file content (header included), given paired-end samples with two ChIP replicates and two control replicates:

| chip_read_1 | chip_read_2 | ctrl_read_1 | ctrl_read_2 |
|-|-|-|-| 
| ... /a.fastq | ... /c.fastq | ... /e.fastq | ... /g.fastq |
| ... /b.fastq | ... /d.fastq | ... /f.fastq | ... /h.fastq |

If your sample is single-ended, then the sample table can simply be formatted as follows:

| chip_read_1 | ctrl_read_1 |
|-|-|
| ... /a.fastq | ... /e.fastq |
| ... /b.fastq | ... /f.fastq |

<br>

### Setting Table & Default Parameters
As in its predecessor, ChIP-AP, the cornerstone of CaRAS’s functionality is the settings table. CaRAS, with the raw fq files and the settings table, is able to reproduce (near) identically any analysis that was performed (provided the same program versions are used). The ‘near identically’ statements is owing to the fact that reported alignments of multi-mappers may, in some cases, give ever so slightly different results. This ambiguity can be alleviated however by filtering out alignments with low MAPQ scores in the corresponding alignment filter step, post-alignment to ensure consistent results from every analysis run. The provision of the settings table therefore ensures reproducibility of any analysis with minimal effort and bypasses the usually sparse and significantly under-detailed methods sections of publications. Science is supposed to be reproducible, yet bioinformatics analysis are typically black-boxes which are irreproducible. This 1 file, changes that!

The structure of the settings table is simple. It is a 2-column tab-separated value file with the names of the programs on the 1st column, and the necessary flags required or changed in the 2nd column. If making your own custom table, then the 1st column below must be copied as-is and not changed. These 2 columns together, list the flags and argument values for each program used in the pipeline.

When CaRAS is run, a copy of the used settings table is saved as a table file: [filename]_ setting_table.tsv in the output save folder. If you have a custom settings table made and provided it as input, then CaRAS will make a 2nd copy of this table in the same output save folder. This decision is made as it is useful documentation of the run performed. This file is also useful for re-running of the pipeline after run failure or some tweaking if necessary. If submitting an issue request on GitHub, you must provide us your settings table used as well as all other requested information. See GitHub for details regarding this.

<u>We consider the dissemination of the information of this file as vital and essential along with results obtained. The table can be included as a supplemental table in a manuscript or can be included as a processed data file when submitting data to GEO – either way, the information of this file must be presented when publishing data.</u>

Below is an example of settings table file in its default-setting state:


| program | argument |
|-|-|
| <b>fastqc1 | -q |
| <b>clumpify | dedupe spany addcount qout=33 fixjunk |
| <b>bbduk | ktrim=r k=21 mink=8 hdist=2 hdist2=1 |
| <b>trimmomatic | LEADING:20 SLIDINGWINDOW:4:20 TRAILING:20 MINLEN:20 |
| <b>fastqc2 | -q |
| <b>bwa_mem |  |
| <b>samtools_view | -q 20 |
| <b>plotfingerprint |  |
| <b>reads_normalizer | -norm properly_paired |
| <b>fastqc3 | -q  |
| <b>macs2_callpeak |  |
| <b>gem | -Xmx10G --k_min 8 --k_max 12 |
| <b>homer_findPeaks | -region |
| <b>genrich | -v |
| <b>seacr | |
| <b>sicer2 | -w 50 -g 100 |
| <b>homer_mergePeaks |  |
| <b>peak_feature_extractor | --filter 4 --top_rank 50 --database biological_process |
| <b>homer_annotatePeaks |  |
| <b>fold_change_calculator | --normfactor user_value --chip_norm 1 --ctrl_norm 1 |
| <b>homer_findMotifsGenome | -size given -mask | 
| <b>meme_chip | -meme-nmotifs 25 |

<br>

## Aligned reads (BAM) normalization formula
For normalization purpose, it has been shown that simply adding the same low amount of spiked-in DNA fragments from a different species to each sample would suffice as a valid basis of normalization and allows accurate quantification of protein occupancy (Skene & Henikoff, 2017). In a later publication on an improved CUT&RUN method, it was shown that the bacterial DNA fragments originating from the purification of pA-MNase complex that was carried over throughout the procedures would suffice to serve the same purpose (Meers, Bryson, Henikoff, & Henikoff, 2019). To calculate the normalization factor for each sample, the reads from the sequencing run need to be realigned onto the spiked-in or the carry-over DNA’s respective genomes and have the number of successfully aligned reads counted. Once the number of aligned reads is determined, it is therefore possible to perform normalization based on the number of these aligned non-sample reads as they ought to be equal across all samples and replicates. CARAS performs this alignment with subsequent normalization with formula as described below:

<img src="https://raw.githubusercontent.com/JSuryatenggara/CaRAS/storage/images/norm_equation_guide.png" width="500">

Given multiple samples involved in a single batch of analysis, the batch multiplier is determined by the largest and the second largest number of mapped spiked-in or carry-over DNA (reference) reads. The normalization factor for each sample is then calculated by dividing the batch multiplier with the corresponding sample’s number of reads. Reads normalization will then be carried out by simply multiplying the original number of reads in each sample with its corresponding normalization factor.

Once the normalization factors are calculated per sample/replicate, it is possible to perform the normalization. To date in the literature, the method for normalizing read counts is to quantify the number of reads per base position throughout the genome in a bedgraph file, and then multiply the number of quantified reads by the normalization factor, i.e., convert aligned BAM file to bedgraph file. The issue though, is that ChIP-Seq peak callers do not call peaks from normalized begraphs, they instead require aligned BAM files. This raised the issue in that to date, there is no published method for normalizing aligned BAM files. To circumvent this issue and to correctly process CnR datasets using traditional ChIP-Seq peak callers, we developed a novel approach for generating normalized BAM (nBAM) files.

The concept of normalization is relatively straight forward, one simply needs to multiply up the number of reads present by the pre-determined normalization factor. While this is a simple operation when working with bedgraphs, when working with aligned reads (.bam) files, this requires “multiplying” each read by the appropriate normalization factor. In other words, if the normalization factor is 5, and only one read is at genomic position X, then we effectively need to have 5 reads at position X instead of 1. So, our approach multiplies each read in the BAM file by the required normalization factor.

In certain datasets with low read depth (e.g., single cell sequencing datasets), the number of mapped spiked-in or carry-over DNA often fluctuate excessively between samples, resulting in a too large of a normalization factor. Due to our direct multiplication of read instances in aligned reads (.bam) files, very large normalization factor may generate a massively oversized files, which will take too much processing resource and time while providing little to no benefit for the analysis quality. Therefore, to avoid this, we applied an upper limit of 10 to the batch multiplier.

Additionally, this approach requires caution in its implementation as in that duplicate reads will appear as “PCR duplicate artefacts” to processing tools. As such, to ensure satisfactory normalization behaviour, each peak caller must be tuned to ignore “PCR duplicates” when processing nBAMs so as to not reverse the normalization of read presence.

<br>

## ChIP vs Control fold change calculation
ChIP weighted peak center coordinate is determined by the median coordinate of all reads in the peak region. Fold change was then calculated as read depth in ChIP sample divided by non-zero read depth of control sample at the weighted peak center coordinate.

To correct for read depth discrepancies due to imbalanced number of reads between ChIP and control samples, CaRAS uses one of three available normalization factors depending on the command line flag specified:

- Based on only uniquely mapped reads in ChIP vs control bam files (filtered using samtools view -F256 and counted using samtools view -c)
- Based on only successfully mapped reads in ChIP vs control bam files (filtered using samtools view -F4 and counted using samtools view -c) 
- (Default) Based on user-determined value.

Advanced users may choose to change this by changing the argument field for fold_change_calculator in the settings table as follows:
- --normfactor uniquely_mapped 
    - Default - based on uniquely mapped reads
- --normfactor mapped 
    - change to all mapped reads
- --normfactor user_value --chip_norm [x] --ctrl_norm [y] 
    - change to user-determined normalization factor, where x / y is the user-determined ratio of the number of ChIP reads to the number of control reads

The feature for read depth normalization (based on total number of reads) before fold change value is calculated was added to compensate for read depth discrepancy between paired ChIP and control samples with same replicate number. This feature was mainly developed for ChIP-AP, where read normalization was not carried out prior, and thus the feature was needed to estimate users a more accurate fold change value of each peak. In CaRAS, read normalization was already performed across all ChIP and control samples (altogether, not separately) directly after read alignment to target organism genome. This normalization method (based on the number of reads successfully mapped to the read normalizer organism genome) is regarded as the standard, and more proper normalization method for CUT&RUN and CUT&Tag. Therefore, CaRAS by default will not perform any read further read depth normalization prior to fold change calculation, as reflected in the default settings table argument for fold_change_calculator: “--normfactor user_value --chip_norm 1 --ctrl_norm 1”. Users may still opt to perform read depth normalization based on total number of reads prior to fold change calculation by modifying the settings table.

<br>

## Irreproducibility Rate (IDR) Calculation
The irreproducible discovery rate (IDR) calculations for all peaks are integrated into a ChIP AP run. As the IDR suite calculate peak reproducibility rate between two replicates (i.e., peak sets) only, we chose as “replicates” each individual peak caller set and the union peak set. The output full peak list (union peak set) is ranked by number of detecting peak callers and fold change of sample signal over input, whereas the individual peak caller sets are ranked by peak score (or fold change if peak score is not available).

The reproducibility of all peaks in the output full peak list are calculated based on their detectability by different individual peak callers and ranked accordingly. These six -log IDR values are copied into the output full peak list, providing six -log IDR values for every peak, which are then summed and converted into a final IDR value.

<img src="https://raw.githubusercontent.com/JSuryatenggara/CaRAS/storage/images/idr_equation_guide.png" width="500">

As reproducibility is not significantly relevant for single-cell analysis with its duplicitous number of cell samples, on top of the time-consuming process were this step is to be applied to every single sample, the IDR calculation step is only applicable to bulk analysis.

<br>

## Peak Feature Extraction
With CaRAS, every peak is given a gene annotation based on the closest gene to its location in the genome. Since every gene would have one or more GO terms or pathways associated to it, enriched terms or pathways can be derived from a list of called peaks, inferring the regulatory effects of the target DNA-binding protein on every single cell sample in the experiment. Therefore, a ranked list of enriched terms or pathways can be used as a profile to compare sample to sample and group them based on how similar their enriched GO terms or pathways are. With this, one can ideally segregate multiple groups of diverse single cells based on their cell types, activation states, or differentiation phases. Each of these groups are then regarded as one entity and should yield a set of results which represents every single cell within that group. Ultimately, CaRAS will follow-up with regular downstream analyses on every group, rather than on every single cell.

To go into the details of how the samples are grouped based on extracted features, CaRAS is employing the methods of rank-biased overlap (RBO) to calculate the distances between each sample’s features (ranked list of enriched terms or pathways) to one another. From these pair-wise distances, an n*n distance matrix is generated, where n is equal to the number of analyzed samples. From here, multidimensional scaling (MDS) method is applied to the distance matrix, transforming the inter-sample distances into a 2D plane where the clusters can later be discernably observed. To determine the optimal number of clusters, CaRAS calculates eigenvalues and their associated eigen vectors from the transformed distance matrix. This is then followed by, with eigengap heuristic, the identification of maximum gap value which corresponds to the number of clusters. Unless provided a pre-determined expected number of resulting clusters, CaRAS will then use the number of clusters with the maximum eigengap in segregating the samples during spectral clustering. Lastly, the spectral clustering is carried out on based on the transformed distance matrix, resulting in grouped samples based on their extracted features, which are the enriched GO terms and pathways.

By default, as recommended in default_settings_table.tsv, CaRAS will use HOMER to perform the gene ontology and pathway analysis only on the subset of peaks which at least called by 4 peak callers. This is done to increase the accuracy of the ranked GO terms or pathways list by taking only the sufficiently confident peaks. CaRAS will also by default only take the top 50 GO terms or pathways for profiling each single cell sample and compare it to one another. As the ranked list profile similarity is done by the RBO (rank-biased overlap) method, this is done with a rationale of removing poorly enriched GO terms or pathways, which has a high proclivity towards randomness in ranked order, and thus will negatively affect the estimation accuracy of profile similarity between each single cell sample. Lastly, CaRAS will by default choose to extract the biological_process GO terms ranked list as profile to represent each single cell sample. Within our limited number of acquired and publicly available datasets, these parameters above are the combination which separate different groups of cells with the clearest segregation and greatest accuracy.

<br>

## Interpreting CaRAS Output
[Some parts adapted from ChIP-AP manual (https://github.com/JSuryatenggara/ChIP-AP)]. CaRAS does report a fair amount of stuff. If you ran it locally you have a swath of folders and you have no clue what to look for and it’s all confusing. We get that. The reality though it’s very simple to know what to look for to know your experimental run worked and we are going to walk you through that!

### Did my analysis work?
There are a couple of things to look for to answer this question. 1, the fingerprint plot and 2, the venn diagram of the merged peaks. Let's begin…

1. <b>The Fingerprint Plot</b>

    The fingerprint plot tells us how well the enrichment of your samples worked. It is generated by the function from the deeptools package and is generated after the alignment files. As such, the plots are found in the “08_results” folder and are labelled “fingerprint_xxxxx.png/svg.” The PNG files allow you to view them in any image viewer, the SVG files are for opening in Adobe Illustrator or Inkscape to make HQ publication figures later if you need.

    <img src="https://raw.githubusercontent.com/JSuryatenggara/CaRAS/storage/images/fingerprint_plot_guide.png" width="500"> 

    To interpret the fingerprint plot, (more information can be found on the deeptools documentation site), but put simply, the input control should be a diagonal line as close as possible toward the 1:1 diagonal. Your ChIP sample should have a bend/kink towards the bottom right corner. The greater the separation between the input and the chip sample, the greater the enrichment you will see in the final result (i.e., lots of peaks). If the lines are overlapping, then you will see little enrichment and your experiment didn’t work that well. If you’re sample lines are switched – then you probably switched the sample names and we recommend doing the right thing and repeating the experiment and not simply switch sample names for the sake of a publication.

    In this example, there is reasonable enrichment in our chip samples. And so we are confident we can see enrichment.

2. <b>The Venn Diagram (well Venn Text)</b>

    In the folder “21_peaks_merging” folder, you will find the “venn.txt” file. This will show you a textual venn diagram of the overlap between the called peaks across all peak callers. To know your experiment worked well, you should see a full list with combinations of all peak callers and relatively large numbers for the consensus peak sets (ie peaks called by multiple peak callers) – this is the ideal case. However, from our experience, there will almost always be 1, maybe 2 peak callers that don’t like a dataset for some reason and so you may find a peak caller performed poorly but the others performed admirably. This is still a good and valid result. 
    
    <img src="https://raw.githubusercontent.com/JSuryatenggara/CaRAS/storage/images/venn_diagram_guide.png" width="500">

    If you look at this file and only see small number of peaks and little overlap, and only 1 peak caller seems to have dominated peak calling, then likely your experiment didn’t work that great. Just because only 1 peak caller performed well though, doesn’t mean the experiment is a write-off and a failure. It can still be valid and so doing some manual validations on the top fold-change differential peaks by chip-PCR might give you an indication whether there is salvageable data or not. Also if you have other confirmatory experimental evidence then even 1 peak caller getting results is fine. This is why we implemented multiple peak callers, because there are many instances where the signal:noise just creates a mess for most peak callers but generally 1 will be the super-hero of the day in such a situation.

3. <b>What Results Files Do I Look At Exactly?</b>

    Valid question. In the folder “22_peak_processing,” open the “xxxx_all_peaks_calculated.tsv” file in excel and you’re good to go. Now to open it there is a little step to do…

    Open a new blank workbook in excel

    <img src="https://raw.githubusercontent.com/JSuryatenggara/CaRAS/storage/images/excel_screenshot_1_guide.png" width="500">
    
    In the ribbon at the top, go to “Data”, then select “From Text/CSV”

    <img src="https://raw.githubusercontent.com/JSuryatenggara/CaRAS/storage/images/excel_screenshot_2_guide.png" width="500">

    In the dialog box that opens up, find and open the peaks files “xxxx_all_peaks_calculated.tsv.” Follow all the prompts and keep pressing “Next” / “Proceed” till the end and the file opens. Opening the peak file this way circumvents an issue that Excel constantly makes which is it will interpret some gene names such as OCT1 as a date, when its not. From this file, you can view all the results and data for your analysis. Refer to Interpreting CaRAS Output for the definition of what each column means.

4. <b>How Do I View My Alignments And Data?</b>
    
    People typically want to view their results on UCSC or other genome browsers. As we don’t have a web-server to host such coverage files (and making accessible UCSC hub is a real pain and we don’t want to implement that), the onus is on you to view them locally on your machine. All laptops, whether then can run CaRAS or not can run IGV Downloads | Integrative Genomics Viewer (broadinstitute.org) and view the coverage and bam files. The coverage and bam files can be located in the “08_results” folder. 
    
    Download IGV, install it (super easy) and then load the coverage and bam files needed. Make sure you load the right genome build however! That’s critical. From section Main Pipeline Output, you can copy columns B,C,D straight into IGV and it will take you to the peak section. 


5. <b>In Short, What's Relevant?</b>  
    1. Check fingerprint plot and make sure it looks good
    2. Check venn.txt file and make sure you get good spread of peaks
    Together points 1 and 2 tell you your experiment worked!
    3. Your final peak file is in “22_peak_processing” open the “xxxx_all_peaks_calculated.tsv” – This is the file you need to upload to GEO as your processed data file for your analysis and the only file you need to work with when looking through your data.
    4. – Also as part of your submission to GEO or as a supplemental table in your manuscript, you MUST include the settings table named “default_settings_table.tsv” located in the root analysis directory. This provided with the raw fq files, which must be uploaded to GEO, will ensure complete reproducibility of the analysis performed.


6. <b>What about the single-cell analysis results? What’s the difference?</b>
    As mentioned before (see section peak feature extraction above), in single-cell analysis, the differences (distance values) between all sample’s extracted features profile are to be projected into a 2D plane and clustered based on the projected their 2D positions to one another. After clustering, the peak lists of all samples which are clustered into a same group will be concatenated into one, namely [setname]_[peakset]_[GO/pathway database]_[group#]_ all_peaks _clustered.tsv. To ascertain that the samples are satisfactorily segregated as expected, user can view the MDS-projected 2D scatter plot, namely [setname]_[peakset]_[GO/pathway database]_Term_Ranking _RBO_Distance_MDS_Spectral_Clustering.png.

    <img src="https://raw.githubusercontent.com/JSuryatenggara/CaRAS/storage/images/spectral_clustering_guide.png" width="500">

    Above is an example of MDS-projected 2D plot, showing a good segregation by spectral clustering of 2 different populations with 20 samples each. Besides the fact that spectral clustering works nicely on these samples 2D distribution, this also means that the algorithm of determining the optimal number by maximum eigengap value is accurate (because in this example, the optimal number of clusters was automatically estimated by CaRAS). In cases where the estimation is off (e.g., 3 clusters), users can explicitly state the optimal number of clusters, which is 2 in this case, via --clustnum flag in peak_feature_extractor.py and rerun the script. Beforehand, if the user already knows how many clusters to expect even without looking at the plot above first, users can explicitly state the expected number of clusters via --clustnum in the CaRAS command line.

    Upon agreeable clustering results, a group-based gene ontology enrichment analysis is performed which results are stored in folder 22_peaks_processing, reflecting the collective cellular roles of the DNA-binding protein in each distinct population. Later down the workflow, motif enrichment analyses by HOMER and MEME are also performed in a group-based manner, resulting in the enrichment of DNA-binding motifs of the target protein in each distinct population. 

<br>

## Manuals and Citations
If you use CaRAS in your analysis, please cite us and all the following programs.

| Program | Reference |
|-|-|
| CaRAS<br>v1.0 | Guide: https://github.com/JSuryatenggara/CaRAS/wiki/CaRAS-Guide<br>Github: https://github.com/JSuryatenggara/CaRAS |
| Python3<br>3.7.x / 3.8.x | We have noted in our testing that there is a change in python 3.8 on macOS in how multi-threading is handled which breaks CaRAS.  As such, for macOS installs you must ensure that ptyhon3.7.x is installed.  If using our installation guides, the provided yml files will ensure all the correct dependencies and requirements are met automatically. |
| FastQC<br>v0.12.1 | Guide: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/<br>GitHub: https://github.com/s-andrews/FastQC |
| Clumpify<br>v38.18 (BBmap) | Introduction: https://www.biostars.org/p/225338/<br>Guide: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify-guide/<br>GitHub: https://github.com/BioInfoTools/BBMap/blob/master/sh/clumpify.sh<br>Citation: https://www.osti.gov/biblio/1241166-bbmap-fast-accurate-splice-aware-aligner
| BBDuk<br>v38.18 (BBmap) | Introduction: http://seqanswers.com/forums/showthread.php?t=42776<br>Guide: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/<br>GitHub: https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk.sh<br>Citation: https://www.osti.gov/biblio/1241166-bbmap-fast-accurate-splice-aware-aligner
| Trimmomatic<br>v0.39 | Guide: http://www.usadellab.org/cms/?page=trimmomatic<br>Downloadable manual page: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf<br>GitHub: https://github.com/timflutre/trimmomatic<br>Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/ |
| bwa<br>v0.7.17 | Guide: http://bio-bwa.sourceforge.net/bwa.shtml<br>GitHub: https://github.com/lh3/bwa<br>Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2705234/ | 
| samtools view<br>v1.9 (samtools) | Guide: http://www.htslib.org/doc/samtools-view.html<br>GitHub: https://github.com/samtools/samtools<br>Citation: https://pubmed.ncbi.nlm.nih.gov/19505943/<br> | 
| deeptools<br>plotFingerprint<br>v3.5.2 (deepTools) | Guide: https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html<br>Citation: https://academic.oup.com/nar/article/44/W1/W160/2499308?login=true |
| MACS2<br>v2.2.7.1 | Guide: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html<br>Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2732366/<br>GitHub: https://github.com/macs3-project/MACS/wiki |
| GEM<br>v2.7 | Guide: https://groups.csail.mit.edu/cgs/gem/<br>GitHub: https://github.com/gifford-lab/GEM<br>Citation: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002638 |
| SICER2<br>v1.0.2 | Guide: https://zanglab.github.io/SICER2/<br>GitHub: https://github.com/bioinf/SICER2<br>Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2732366/ |
| SEACR<br>v1.3 | GitHub: https://github.com/FredHutch/SEACR<br>Citation: https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4 |
| HOMER<br>findPeaks<br>v4.11 (HOMER) | Guide: http://homer.ucsd.edu/homer/ngs/peaks.html<br>Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2898526/ |
| Genrich<br>v0.6.1 | Guide: https://informatics.fas.harvard.edu/atac-seq-guidelines.html<br>GitHub: https://github.com/jsh58/Genrich |
| HOMER<br>mergePeaks<br>v4.11 (HOMER) | Guide: http://homer.ucsd.edu/homer/ngs/mergePeaks.html<br>Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2898526/ |
| HOMER<br>annotatePeaks<br>v4.11 (HOMER) | Guide: http://homer.ucsd.edu/homer/ngs/annotation.html<br>Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2898526/ |
| IDR<br>v2.0.4.2 | GitHub: https://github.com/nboley/idr<br>Citation: https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-3/Measuring-reproducibility-of-high-throughput-experiments/10.1214/11-AOAS466.full |
| HOMER<br>findMotifsGenome<br>v4.11 (HOMER) | Guide: http://homer.ucsd.edu/homer/ngs/peakMotifs.html<br>Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2898526/
| MEME<br>meme-chip<br>V5.0.5 (MEME) | Guide: https://meme-suite.org/meme/doc/meme-chip.html<br> Citation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2703892/ | 
