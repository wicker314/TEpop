TEpop
Transposable Element Population Analysis

+++++ WICKERsoft set up +++++++++++++++++
    • all scripts depend on subroutines stored in the perl module Biomodule.pm.
    • place the Biomodule.pm and all WICKERsoft scripts into a dedicated directory which is in a common path (e.g. /usr/local/bin/) so that you can run them from anywhere in your system.
    • make sure that all scripts are executable. if necessary, make them executable with chmod 755 <script_name>
    • change the path in the beginning of all WICKERsoft scripts to that directory
	example, how the header should look: 
	#!/usr/bin/perl
	use lib '/usr/local/bin/'; 
	use Biomodule;
      
Additional software needed:

From standard Ubuntu repositories:
- Emboss
- Clustalw
- blastall / ncbi-blast+-legacy (local stand-alone blast program)
- perl-tk
- libgd-graph-perl
- libgd-perl
- libgd-text-perl

R-packages:
- R-base
- library(SNPRelate)
- library(magrittr)



Comments:
    • all WICKERsoft scripts have usage instructions to tell you what type of input files and parameters are needed. simply type the script name without arguments to see usage instructions
    • Wild cards (* or ?): WICKERsoft scripts do not recognize wild cards. So, if you are asked to provide a query for multiple files, provide a string common to all file names you want to process without wild cards.
    • some WICKERsoft scripts use the perl-tk module for graphics. Thes have the additional line 'use Tk;' in the header. The perl-tk module can be obtained from Ubuntu repositories.
    • Naming convention: We name TEs in the following way:
      3-letter-code_4-letter-species_family_source-seq-copyID_comment, 
      for example: RLG_Taes_Angela_chr1A_003245_TSD 


TEpop Protocol:
The TEpop pipeline is designed for LTR retrotransposons only
The pipeline is run for each LTR retrotransposon family separately
A reliable consensus sequence of a full-length element is needed to start with. 

Step 1- Compile LTR sequences +++++++++++++++++++++++++++++
- Isolate the LTR sequence from your LTR retrotransposon consensus sequence. Pay attention that the ends are precisely determined (LTRs most often start with TG and end in CA, TA/CA is occasionally found).

- The TEpop pipeline can also be run with a single LTR sequence, at the risk of missing subfamilies with divergent LTRs 
I you want to do this, proceed to Step 2.

- Usually, LTR retrotransposon families subdivide into multiple subfamilies, which often differ considerable in their LTR sequence. Since the second step of the pipeline relies on full-length blast hits, it is advisable to use LTRs from multiple copies, or extract additional diverse LTRs from the genome, e.g. based on a blastn search of the initial LTR.

- all your input LTR fasta file names should end in _LTR




Step 2 - TEpop_RLX_2024 : isolation of full-length copies form chromosomes +++++++++++
Sine the TEpop pipeline is designed for large genomes, it requires individual chromosome fasta sequences and blast-able databases for each individual chromosome in dedicated directories. 
By default the chromosomes are in: /home/wicker/data/dir_chromosomes/
and corresponding blast-able databases (extensions .nhr, .nin, .nsq) are in:  /home/wicker/data/databases/blast_databases/
These paths have to be edited in the TEpop_RLX_2024 script to fit your local set up.

Run script  TEpop_RLX_2024 according to instructions in the usage, provide arguments. 
1.= end of file name query for LTRs (i.e. _LTR)
2.= minimal alignment size (close to size of smallest LTR used)
3.= size range of element (e.g. 6000-8000)
4.= full TE family name (e.g. RLC_HvulMorexv3_BARE1)
5.= name query for databases (e.g. Bgt_96224_v316_chr)
6.= allowed number of bases missing at LTR ends (e.g. 5)

The script TEpop_RLX_2024 does the following:
1. runs blastn searches of all LTR input sequences against the queried chromosomes
2. identifies cases where (near-complete) LTRs are found in the same orientation and in the required distance range
3. loads respective chromosome into memory and excises full-length copies
4. checks fro presence of target site duplications (TSDs)
5. Isolates solo-LTRS (which are usually not used)

Comments:
- Isolated full-length copies are written as fasta files into the working directory. They have the file name extension _TSD, indicating full-length copies. (up to 2 mismatches in TSD are allowd)
- TE copy names contain chromosome of origin and position in kb (e.g. RLG_Taes_Angela_chr1A-003452_TSD. Including the position in kb in the name has the advantage that if the same copy is identified by multiple LTRs in the the blast searches, the fasta file is simply over-written.
- it is recommended to check the size distribution of the isolated copies. Usually, it has a clear and narrow peak. We recommend to discard the shortest and longest 5% as they may contain copies with insertions of other TEs, or cases where 2 solo-LTRs were found by chance in the same orientation and in the required distance range


Step 3 - Determination of insertion age
Insertion age of individual copies is estimated by comparing the two LTRs and, counting transition and transversions and and applying a substitution rate of 1.6E-9 per site per year.

Run the program LTR with arguments:
1.= query for full-length retros (e.g. RLG_Taes)
2.= PBS or IDs to assemble PBS or ID or N for nothing (us N here)

- The program aligns the first and second half of each copy. Usually, this precisely aligne the 2 LTRs
- The program write out one lTR sequence as well as a .pair file which contains the alignment of the two LTRs

Date insertions with date_pair, using the following arguments:
1.= query for pair files
2.= optional name for log file (e.g. RLC_Hvul_BARE1)

this produces an output file, e.g. age_distribution_RLC_Hvul_BARE1 with numbers of aligned base, transitions and transversions, and estimated insertion age in million years (6th column).

This file will be needed  in later steps.

The files _LTR and _LTR.pair will not be needed anymore and should be deleted.

Step 4 – produce a variant call file (.vcf)
Here, all copies are aligned with a consensus sequence which will then be used to produce a .vcf file. 
Option 1: Use an existing consensus sequence of a full-length TE
Option 2: randomly pick 15-30 copies, align them with clustalw and produce a ne consensus using the program consensus.

Individual copies are aligned with the consensus with the program bestflash which implements the WATER program form the EMBOSS package, using the following arguments:

1.= sequence with which all others are compared (e.g. RLG_Taes_Angela_consensus-1)
2.= string common to all others (e.g. RLG_Taes_Angela_chr)
3.= gap open penalty (e.g. 50)
4.= gap extension penalty (e.g. 0.1)

This produces an alignment of each copy with the consensus (.pair file)

The vcf is produced with the program pair_to_vcf, which produces a consensus-based vcf. It ignores inDels and uses the following arguments

1.= reference sequence (e.g. RLG_Taes_Angela_consensus-1)
2.= minimum allele frequency for major allele (e.g. 0.05)
3.= length of seq to use (e.g. only LTR for veeery large dataest, but better to use just a large number, e.g. 20000)
4.= outfile name (optional, eg. RLG_Taes_Angela, default = ref_seq.vcf)

Step 5 – TE population analysis
In this step, we produce the PCA which is the main product of the TEpop pipeline.
The vcf from the previous step is used as input for the program Tepop_R_script, which uses the following arguments:

1.= vcf file
2.= number of cores
3.= norm or rnd (for randomized algorithm, norm is recommended, rnd for speed)
4.= -s to only show script (optional, if you only want to see the script)

This produces 3 files
1. pca_table that can be used fro visualization (see below)
2. pca_variance table, summariziong how muh of the variance is explained by each PC
3. a png of the PCA

The PCA can be visualized with visual_PCA_Tepop_2024, which uses the following arguments:
1.= PCA table produced in the previous step
2.= age distribution table produced in step3
3.= centromere positions (optional)
4.= fix for fixed canvas size (only necessary of you have multiple screens)

The program allows you to do the following:
- Start visualization with the process_file button
- Visualize the PCA
- adjust PCA image size with the x/y zoom entries
- visualize insertion ages (check the show_age button)
- define and isolate subfamilies (select a section of the PCA with the mouse)
- write fasta files of selected subfamilies
- export .png images of the PCA


