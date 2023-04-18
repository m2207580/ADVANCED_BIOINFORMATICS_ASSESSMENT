#!/bin/bash

#ASSESSMENT BY STUDENT m2207580

##INPUT FILES
#$1 should be the Fastq Read 1 file
#$2 should be the Fastq Read 2 file
#$3 should be the Annotation BED file

#Firstly, we want to make sure the files we are inputting exist in the machine, otherwise the script won't run
if [ ! -f "$1" ]
then
    echo "File $1 not found"
    exit -1
fi

if [ ! -f "$2" ]
then
    echo "File $2 not found"
    exit -1
fi

if [ ! -f "$3" ]
then
    echo "File $3 not found"
    exit -1
fi


##SETTING UP OUR WORKING SPACE
#Before we start our pipeline, we need to do a series of steps to set the directories.
mkdir /home/ubuntu/ASSESSMENT
mkdir /home/ubuntu/ASSESSMENT/data /home/ubuntu/ASSESSMENT/logs /home/ubuntu/ASSESSMENT/meta /home/ubuntu/ASSESSMENT/results
mkdir /home/ubuntu/ASSESSMENT/data/aligned_data /home/ubuntu/ASSESSMENT/data/reference /home/ubuntu/ASSESSMENT/data/trimmed_fastq /home/ubuntu/ASSESSMENT/data/untrimmed_fastq
mkdir /home/ubuntu/ASSESSMENT/results/fastqc_trimmed_reads /home/ubuntu/ASSESSMENT/results/fastqc_untrimmed_reads

##QUALITY ASSESSMENT OF UNTRIMMED FILES
fastqc -t 4 $1 $2 #The option -t provides the opportunity to have a number of files run simultaneously, so we don’t waste time in waiting for the files to be analyzed serially.
##Now we  move the outputs to the correct results folder
mv /home/ubuntu/NGS0001.R1_fastqc.html /home/ubuntu/ASSESSMENT/results/fastqc_untrimmed_reads
mv /home/ubuntu/NGS0001.R1_fastqc.zip /home/ubuntu/ASSESSMENT/results/fastqc_untrimmed_reads
mv /home/ubuntu/NGS0001.R2_fastqc.html /home/ubuntu/ASSESSMENT/results/fastqc_untrimmed_reads
mv /home/ubuntu/NGS0001.R2_fastqc.zip /home/ubuntu/ASSESSMENT/results/fastqc_untrimmed_reads


##TRIMMING OF RAW FILES
#Before trimming our raw sequence data, we must unzip the files and move them to the right folder:
zcat $1 > /home/ubuntu/ASSESSMENT/data/untrimmed_fastq/NGS0001.R1.fastq
zcat $2 > /home/ubuntu/ASSESSMENT/data/untrimmed_fastq/NGS0001.R2.fastq


#Now, we run the trimming:
trimmomatic PE -threads 4  -phred33 /home/ubuntu/ASSESSMENT/data/untrimmed_fastq/NGS0001.R1.fastq /home/ubuntu/ASSESSMENT/data/untrimmed_fastq/NGS0001.R2.fastq -baseout /home/ubuntu/ASSESSMENT/data/trimmed_fastq/trimmed_data ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50


##QUALITY ASSESSMENT OF TRIMMED READS
fastqc -t 4 /home/ubuntu/ASSESSMENT/data/trimmed_fastq/trimmed_data_1P /home/ubuntu/ASSESSMENT/data/trimmed_fastq/trimmed_data_2P 

mv /home/ubuntu/*fastqc* /home/ubuntu/ASSESSMENT/results/fastqc_trimmed_reads/ #we move our resulting files to the results folder


##ALIGNMENT WITH BWA-MEM
#Firstly, we must download the reference genome file and move it to the right folder:
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv /home/ubuntu/hg19.fa.gz /home/ubuntu/ASSESSMENT/data/reference/

#Before being able to perform the aligment, BWA needs to build a reference index that takes around 1 hour to complete. For this reason, it is recommended that indexing the reference genome be run outside the script and be considered a "project set up" step. The code for this step would be:
bwa index /home/ubuntu/data/reference/hg19.fa.gz

#Now we are ready to run the alignment:
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1:111:D1375ACXX:1:1101\tSM:NGS0001|tPL:ILLUMINA\tLB:nextera-ngs0001\tDT:20233-04-10\tPU:11V6WR1' -I 250,50 /home/ubuntu/ASSESSMENT/data/reference/hg19.fa.gz /home/ubuntu/ASSESSMENT/data/trimmed_fastq/trimmed_data_1P /home/ubuntu/ASSESSMENT/data/trimmed_fastq/trimmed_data_2P > /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001.sam # "t 4" gives us four threads, "v 1" gives us verbosity level of error, "R" tells how to read the group header line, "I 250,50" specify the mean as 250 and standard deviation as 50.


##POST-ALIGNMENT QUALITY ASSESSMENT AND FILTERING
#First, we need to convert from SAM to BAM, sort and index the resulting BAM file.
samtools view -h -b /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001.sam > /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001.bam #convert our aligned SAM file into BAM file,“h” makes the header be included and “b” makes the output be in BAM format

samtools sort /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001.bam > /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted.bam #to sort the alignment file and generate a new one

samtools index /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted.bam > /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted.bam.bai #index the alignment and generate a .bai file

#Now, we perform duplicate markings on our aligned BAM file
picard MarkDuplicates I=/home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted.bam O=/home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_marked.bam M=/home/ubuntu/ASSESSMENT/data/aligned_data/marked_dup_metrics.txt # “I” indicates the input, “O” indicates the output and “M” indicates the metrics file

#Now, we filter the duplicate-marked BAM file.
samtools view -F 1796  -q 20 -o /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_filtered.bam /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_marked.bam #the number after “F” stands for the filtering criteria to skip alignments if i. The read is unmapped ii. The alignment or this read is not primary iii. The read fails platform/vendor quality checks iv. The read is a PCR or optical duplicate. The “q 20” sets the minimum MAPQ quality score to 20

#It is important to remember to index the resulting bam file too:
samtools index /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_filtered.bam > /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_filtered.bam.bai

#We now generate standard alignment statistics:
#Flagstats: Counts the number of alignments for each of the 13 FLAG types
samtools flagstats /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_filtered.bam
#Idxstats: Reports alignment summary statistics per chromosome in the indexed file
samtools idxstats /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_filtered.bam
#Depth of coverage: This tool computes both the depth and breadth of coverage of features in the annotation file (b) on the features in the BAM file (a). “d” reports the depth at each position in each A feature.
bedtools coverage -d -a /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_filtered.bam -b $3
#Distribution of insert sizes: This tool collects metrics about the insert size distribution of a paired-end library. “I” corresponds with the input BAM file, “O” corresponds with the text output file with the stats, and “H” with the histogram pdf. “M” indicates the minimum percentage option is 50%
picard CollectInsertSizeMetrics I=/home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_filtered.bam O=/home/ubuntu/ASSESSMENT/data/aligned_data/insert_size_metrics.txt H=/home/ubuntu/ASSESSMENT/data/aligned_data/insert_size_histogram.pdf M=0.5


##VARIANT CALLING
#Firstly, we must unzip the reference sequence FASTA file:
zcat /home/ubuntu/ASSESSMENT/data/reference/hg19.fa.gz > /home/ubuntu/ASSESSMENT/data/reference/hg19.fa
#Then we index it:
samtools faidx /home/ubuntu/ASSESSMENT/data/reference/hg19.fa > /home/ubuntu/ASSESSMENT/data/reference/hg19.fa.fai
#Now, we variant call with Freebayes:
freebayes -f /home/ubuntu/ASSESSMENT/data/reference/hg19.fa /home/ubuntu/ASSESSMENT/data/aligned_data/NGS0001_sorted_filtered.bam >/home/ubuntu/ASSESSMENT/results/NGS0001.vcf 

#Now we filter the vcf file for the regions in our annotation.bed file
bedtools intersect -header -wa -a /home/ubuntu/ASSESSMENT/results/NGS0001.vcf -b $3 > /home/ubuntu/ASSESSMENT/results/NGS0001_annotation.vcf

#Now we compress the resulting VCF with:
bgzip /home/ubuntu/ASSESSMENT/results/NGS0001_annotation.vcf > /home/ubuntu/ASSESSMENT/results/NGS0001_annotation.vcf.gz
#And then index the VCF:
tabix -p vcf /home/ubuntu/ASSESSMENT/results/NGS0001_annotation.vcf.gz > /home/ubuntu/ASSESSMENT/results/NGS0001_annotation.vcf.gz.tbi #where “p” indicates the format of the input file

#For variant quality filtering:
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" /home/ubuntu/ASSESSMENT/results/NGS0001_annotation.vcf.gz > /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.vcf #where we establish QUAL > 1: removes horrible sites, QUAL / AO > 10 : additional contribution of each observation should be 10 log units (~ Q10 per read), SAF > 0 & SAR > 0 : reads on both strands, RPR > 1 & RPL > 1 : at least two reads “balanced” to each side of the site

#Afterwards, we compress the filtered VCF:
bgzip /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.vcf > /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.vcf.gz
#And we index once again:
tabix -p vcf /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.vcf.gz > /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.vcf.gz.tbi


##VARIANT ANNOTATION AND PRIORITIZATION
#First, we convert our final compressed VCF file to ANNOVAR input:
/home/ubuntu/annovar/convert2annovar.pl -format vcf4 /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.vcf.gz > /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.vcf.avinput

#Now we run the ANNOVAR table function to output a csv file:
/home/ubuntu/annovar/table_annovar.pl /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.vcf.avinput /home/ubuntu/annovar/humandb/ -buildver hg19 -out /home/ubuntu/ASSESSMENT/results/NGS0001_annovar -remove -protocol refGene,ensGene,snp142,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f,f -otherinfo -nastring . -csvout


#We can also do variant annotation with another tool called snpEff (after installing it as a set-up step with conda)
snpEff -canon hg19  /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.vcf  > /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.snpeff.vcf
#To add the dbSNP database as extra annotation:
snpSift annotate db/dbSnp/dbSnp142 /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.snpeff.vcf > /home/ubuntu/ASSESSMENT/results/NGS0001_filtered_annotated.snpeff.vcf

#Now, we can perform basic variant prioritization on the csv ANNOVAR output. Firstly, we must use head to copy the header row from our ANNOVAR csv to a new output file:
head -n 1 /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.hg19_multianno.csv > /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.hg19_multianno_final.csv
#Now we filter:
awk -F, '(($6 ~ /"exonic"/) && ($11 ~ /"exonic"/) && ($16 ~ /./) {print}' /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.hg19_multianno.csv >> /home/ubuntu/ASSESSMENT/results/NGS0001_filtered.hg19_multianno_final.csv #This command tells the machine to append the resulting rows to the created csv with the header after selecting those with “exonic” in column 6 and 11 (Func.refGene and Func.ensGene) and “.” in column 16 (dbSNP) in a new csv file, and “F” indicates the type of field separator

