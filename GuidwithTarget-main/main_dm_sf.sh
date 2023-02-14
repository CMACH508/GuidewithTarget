: << 'COMMENT'
This file main pipeline to process degradome
#NGdata store has three pool
#Degradome:
  /home/zenglin/LinZeng/degradome/DrosoDataRaw
  adhago3 vs. wtago3
  mutago3 vs. wtago3
  adhaub vs. wtaub
  mutaub vs. wtaub
  ago3_pi.fa
  aub_pi.fa
  non_ago3_pi.fa
  non_aub_pi.fa
Usage Example:
	bash script_final/main_dm.sh -a wtago3 -b mutago3 -c ago3 -r 2 -p ago3_pi.fa -q non_ago3_pi.fa
Parameter Explain:
  -a wt_file_prefix: wtsptd
  -b mut_file_prefix: mutsptd
  -c cell line: ss sptd pac dip
  -r reference: genomic or transcriptomic GENOMIC/TRANSCRIPTOMIC
  -n normalize by library size or not Y/N
Preparation:
  # download SRR2042571
  wget -O ago3_krimp.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR2042571/SRR2042571
  # get ago3_krimp smallRNA sequencing data
  fastq-dump ago3_krimp.sra --split-3 --gzip -O ./
  # unzip the file
  gunzip -c ago3_krimp.fastq.gz > ago3_krimp.fastq
  # transfer to fasta
  cat ago3_krimp.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ago3_krimp.fasta
  # build bowtie index
   bowtie-build --threads 4- ago3_krimp.fasta ago3_krimp
  # dowload dm piRNA sequence fasta file from piRNAdb 
  https://www.pirnadb.org/download/archive
  # map Drosophila melanogaster - v1.7.6 fasta file to ago3_krimp.fasta
  bowtie --sam-nohead -f -p 40 -t -x ago3_krimp ../pidroso/piRNAdb.dme.v1_7_6.fa -S ago3_piRNA.sam
  # fitler sam file
  less ago3_piRNA.sam | cut -f 1,2,3,4 | awk '($2 == "0")' > ago3_filter_cnt.sam
  # fetch sequence of piRNA of interest (POI) and background sequence
  Rscript POI_BG_fetchfasta.R ago3_filter_cnt.sam ../pidroso/piRNAdb.dme.v1_7_6.fa
  # Reminder: Move XX_BG.fa and XX_IP.fa to piRNA_add file directory  
COMMENT

func() {
    echo "Usage:"
    echo "bash main_new.sh -a wtsptd -b pisptd -c sptd -c stpd -r 2 -x Y" 
    echo "Description:"
    echo "S_DIR,the path of source."
    echo "D_DIR,the path of destination."
    exit -1
}

while getopts 'a:b:c:r:e:p:q:h' flag
do
  case "${flag}" in
	a) wt_file_prefix=${OPTARG};;
	b) mut_file_prefix=${OPTARG};;
  c) cell_line=${OPTARG};;
  e) evalue=${OPTARG};;
  r) reference=${OPTARG};;
	p) piRNA=${OPTARG};;
  q) nonpiRNA=${OPTARG};;
  h) func;;
  esac
done

echo "wt_file: $wt_file_prefix"
echo "mut_file: $mut_file_prefix"
echo "cell_line: $cell_line"
echo "reference: $reference"
echo "POI: $piRNA"
echo "NPOI: $nonpiRNA"
echo "norm: $normORnot"



: << 'COMMENT'
Appointment of Address, pay attention to back splash of each address!
COMMENT
raw_data_add="/home/zenglin/LinZeng/degradome/DrosoDataRaw/"
data_add="/home/zenglin/LinZeng/degradome/DrosoDataStore/"
genome_add="/home/zenglin/LinZeng/degradome/droso_genome/"
piRNA_add="/home/zenglin/LinZeng/degradome/pidroso/"
script_add="/home/zenglin/LinZeng/degradome/script_final/"
data_add="${data_add}${cell_line}/"




: << 'COMMENT'
Directory preparation!
Here we index sequencing file as reference, and map pi.fa to reference
We witnessed pi6.fa mapping to mut sequencing sample around 15% but 70% to wt sequencing
So it proves efficiency of knock out!
COMMENT
if [ -f $data_add ]
then
echo "Directory exists"
else
mkdir $data_add
echo "Directory does not exists"
echo "Establish this for you!"
fi




: << 'COMMENT'
Preprocess Data
Using FastP to automatially detech and splice adaptor
COMMENT
fastp -w 16 -i "${raw_data_add}${wt_file_prefix}_1.fastq.gz" \
 -o "${data_add}wt${cell_line}_1_fp.fastq" \
 -I "${raw_data_add}${wt_file_prefix}_2.fastq.gz" \
 -O "${data_add}wt${cell_line}_2_fp.fastq" \
 -R "${data_add}wt${cell_line}_deg"
mv fastp.html "${data_add}wt${cell_line}_deg.html" #wtago3_deg
mv fastp.json "${data_add}wt${cell_line}_deg.json" #wtago3_deg

fastp -w 16 -i "${raw_data_add}${mut_file_prefix}_1.fastq.gz" \
 -o "${data_add}mut${cell_line}_1_fp.fastq" \
 -I "${raw_data_add}${mut_file_prefix}_2.fastq.gz" \
 -O "${data_add}mut${cell_line}_2_fp.fastq" \
 -R "${data_add}mut${cell_line}_deg"
mv fastp.html "${data_add}mut${cell_line}_deg.html" #mutago3_deg
mv fastp.json "${data_add}mut${cell_line}_deg.json" #mutago3_deg



: << 'COMMENT'
Preprocess Data
Map Data to Reference
Here 2 options for reference, genome or transcriptome
COMMENT
if [ $reference -eq 1 ]
#if reference == 1, then map to transcriptome
then
echo "Mapping to transcriptome:"
########################################################
#HARVEST products in wide type sample
/home/zenglin/miniconda3/envs/clash/bin/bowtie-align-s  \
 -p 40 \
 -t \
 -x "${genome_add}dm6_mRNA" \
 "${data_add}wt${cell_line}_1_fp.fastq" \
 -S "${data_add}wt${cell_line}.sam" 
tail +75221  "${data_add}wt${cell_line}.sam" | \
 cut -f 1,2,3,4 | \
 awk '($2 != "4")' \
 > "${data_add}wt${cell_line}_candidate.sam"
#HARVEST products in mut type sample
/home/zenglin/miniconda3/envs/clash/bin/bowtie-align-s  \
 -p 40 \
 -t \
 -x "${genome_add}dm6_genome" \
 "${data_add}mut${cell_line}_1_fp.fastq" \
 -S "${data_add}mut${cell_line}.sam"
tail +75221 "${data_add}mut${cell_line}.sam" | \
 cut -f 1,2,3,4 | \
 awk '($2 != "4")' \
 > "${data_add}mut${cell_line}_candidate.sam"
else
echo "Mapping to genome:"
########################################################
#OR map to genome
#HARVES products in wide type 
/home/zenglin/miniconda3/envs/clash/bin/bowtie-align-s \
 -p 40 \
 -t \
 -x "${genome_add}dm6_genome" \
 "${data_add}wt${cell_line}_1_fp.fastq" \
 -S "${data_add}wt${cell_line}.sam" 
tail +1873 "${data_add}wt${cell_line}.sam" | \
 cut -f 1,2,3,4 | \
 awk '($2 != "4")' \
 > "${data_add}wt${cell_line}_candidate.sam"
#HARVESTS products in mut type
/home/zenglin/miniconda3/envs/clash/bin/bowtie-align-s \
 -p 40 \
 -t \
 -x "${genome_add}dm6_genome" \
 "${data_add}mut${cell_line}_1_fp.fastq" \
 -S "${data_add}mut${cell_line}.sam"
tail +1873 "${data_add}mut${cell_line}.sam" |\
 cut -f 1,2,3,4 |\
 awk '($2 != "4")' \
 > "${data_add}mut${cell_line}_candidate.sam"
fi
samtools flagstat "${data_add}wt${cell_line}.sam" > "${data_add}wt${cell_line}_flagstat.txt"
samtools flagstat "${data_add}mut${cell_line}.sam" >  "${data_add}mut${cell_line}_flagstat.txt"
rm "${data_add}wt${cell_line}.sam" 
rm "${data_add}mut${cell_line}.sam"
#samtools view -S -b  "${data_add}wt${cell_line}.sam"> "${data_add}wt${cell_line}.bam" 
#samtools sort  "${data_add}wt${cell_line}.bam" > "${data_add}wt${cell_line}.sorted.bam"
#samtools index "${data_add}wt${cell_line}.sorted.bam"
#samtools view -S -b  "${data_add}mut${cell_line}.sam"> "${data_add}mut${cell_line}.bam"
#samtools sort  "${data_add}mut${cell_line}.bam" > "${data_add}mut${cell_line}.sorted.bam"
#samtools index "${data_add}mut${cell_line}.sorted.bam"


: << 'COMMENT'
Prepare for Cutting edge of each splicing position!
Fetch the sequence of cutting edge 
COMMENT
# We hope to get the first position of 5'end of 3'Product
# Input: samfile 1 and samfile 2
# output: cellline_insert.txt // cellline_insert.bed
Rscript "${script_add}fetch_DegraPosition_dm.R" \
 "wt${cell_line}_candidate.sam" \
 "mut${cell_line}_candidate.sam" \
 $data_add \
 $reference
if [ $reference -eq 1 ]
then
echo "Extracting Sequence from Transcriptome File:"
bedtools getfasta -s -fi "${genome_add}mrna.fa" \
	-bed "${data_add}${cell_line}_insert.bed" \
	-name > "${data_add}${cell_line}_insert.fa"
else
echo "Extracting Sequence from Genome File:"
bedtools getfasta -s -fi "${genome_add}dm6.fa" \
        -bed "${data_add}${cell_line}_insert.bed" \
        -name > "${data_add}${cell_line}_insert.fa"
fi


# Establish library
makeblastdb -in "${data_add}${cell_line}_insert.fa" \
  -dbtype nucl \
  -max_file_sz 2GB \
  -out "${data_add}${cell_line}_insert_db"

# map piRNA of Interest to 3 Cleaved Products
blastn -query "${piRNA_add}${piRNA}" \
  -num_threads 10 \
  -db "${data_add}${cell_line}_insert_db" \
  -strand minus -task blastn \
  -outfmt 0 -ungapped > "${data_add}${cell_line}_POI_Target_Pair.txt"
Rscript "${script_add}Parse_POI_Target_Pair.R" \
  "${cell_line}_POI_Target_Pair.txt" \
  $data_add \
  $reference 


# map nonpi6 need loop
for i in {1..10}
do
mkdir "${data_add}blast_result"
echo "Sample Times: ${i}"
Rscript "${script_add}RandomSampling_NPOI.R" \
  "${piRNA_add}${nonpiRNA}" \
  "${piRNA_add}${piRNA}" \
  $i \
  $data_add \
  $cell_line
blastn -query "${data_add}${cell_line}_nonPiRNA_${i}.fa" \
  -num_threads 10 \
  -db "${data_add}${cell_line}_insert_db" \
  -strand minus \
  -task blastn \
  -outfmt 0 \
  -ungapped > "${data_add}${cell_line}_NPOI_Target_Pair_${i}.txt"
Rscript "${script_add}Parse_NPOI_Target_Pair.R" \
  "${cell_line}_NPOI_Target_Pair_${i}.txt" \
  $data_add \
  $reference \
  $i 
cd "${data_add}blast_result"
cat *.txt > "${cell_line}_Parsed_Filter_NPOI_T${i}.txt"
echo "${cell_line}_Parsed_Filter_NPOI_T${i}.txt"
mv "${cell_line}_Parsed_Filter_NPOI_T${i}.txt" $data_add
cd $data_add
rm -r "${data_add}blast_result"
done

Rscript "${script_add}FullMatch_dm.R" "${cell_line}_Parsed_Filter_POI.txt" $data_add $cell_line 20 $reference
Rscript "${script_add}FullMatch_dm.R" "${cell_line}_Parsed_Filter_POI.txt" $data_add $cell_line 25 $reference
Rscript "${script_add}FullMatch_dm.R" "${cell_line}_Parsed_Filter_POI.txt" $data_add $cell_line 30 $reference

Rscript "${script_add}EasyStatProducts.R" "${data_add}${cell_line}_Parsed_Filter_POI.txt"
