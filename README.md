# GuidewithTarget (Degradome Sequencing Bioinformatics Pipeline)

<img src=".\Illustration.PNG" width="100%" />

A pipeline is built to analyze “guide:target” complementarity for Ago3-catalyzed splicing of transcripts guided by piRNAs in Drosophila Melanogaster. The analysis contains three parts, i.e., identification of Ago3 associated piRNAs from RNA immunoprecipitation (RIP) sequencing libraries, identification of cleavage products from degradome sequencing libraries, and pairing analysis between piRNAs and target transcripts. The pipeline is implemented via customized home-made scripts which involve several popular bioinformatics tools. Home-made scripts were written in R environment (v4.2.0), with ggplot2 (v3.4.0) package for data visualization. 

# Dependencies
GuidewithTarget is built based on following open source packages:
- sra-tools 2.8.0 
- r-base 4.2.0
- r-ggplot2 3.4.0
- bedtools 2.30.0
- blast 2.5.
- fastp 0.22.0
- samtools  1.7

# Data Preparation
- anti-Ago1 RIP-seq (SRR2042571)
- piRNAdb.dme.v1_7_6.fa from https://www.pirnadb.org/download/archive
- ago3 mutant (SRR1926158)
- ago3 wide type (SRR1926188)

# Workflow
- Prepare Ago3_PI.fa/Ago3_bg.fa
- Analysis of Degradome Sequencing
- Blast Ago3_PI.fa/Ago3_bg.fa to products
- Data Visualization

# Usage Tutorial
bash main_dm_sf.sh -a wtago3 -b mutago3 -c ago3 -r 2 -p ago3_pi.fa -q non_ago3_pi.fa
- a: ago3 wide type file prefix
- b: ago3 mutant file prefix
- c: protein
- r: blast reference, 2 referes to genome, 1 refers to transcriptome fasta
- p: piRNA of Interest (POI)
- q: non piRNA of Interest (NPOI)

Warning:

main_dm_sf.sh is the main bash script to combi all home-made scripts with popular bioinformatics tools, whole pipeline requires five directories addresses.
- raw_data_add: raw data address which stores raw data of degradome, fastq etc.
- data_add: working directory which stores all output from pipeline.
- genome_add: genome.fa or transcriptome.fa locates in this directory.
- piRNA_add: Ago3_IP.fa or Ago3_bg.fa locates
- script_add: put all home-made scripts in.
