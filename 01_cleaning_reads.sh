#!/bin/sh

############# COMBEE 2021-2023 bumblebee genetic and microbiome data ################ 

#### using bbmap 38.92

# batch_number is the folder for each batch 

batch_number="/scratch/patrycja/COMBEE_lapidarius/batch3"
# genome=/scratch/patrycja/reference_genomes/Bombus_pascuorum/GCF_905332965.1_iyBomPasc1.1_genomic.fna
genome=/scratch/patrycja/reference_genomes/Bombus_lapidarius/GCA_964186655.1_iyBomLapi1.1_genomic.fna
prop_aligned="counts_summary.txt"
mkdir $batch_number/Clean
mkdir "$batch_number"/Aligned
# export PATH=$PATH:/home/panas/Downloads/gatk-4.2.2.0
# export PATH=$PATH:/home/panas/Downloads/bbmap

### Cleaning loop

for i in "$batch_number"/raw/*_R1_001.fastq.gz; \
	do dname=$(dirname ${i}); name=$(basename ${i} _R1_001.fastq.gz); \
	bbduk.sh in1=${dname}/${name}_R1_001.fastq.gz in2=${dname}/${name}_R2_001.fastq.gz \
	out1=$batch_number/Clean/${name}_clean_1.fastq.gz out2=$batch_number/Clean/${name}_clean_2.fastq.gz \
	ref=/home/panas/Downloads/bbmap/resources/nextera.fa.gz ktrim=r k=17 mink=8 hdist=1 tpe tbo qtrim=rl trimq=30 ordered=t threads=112 stats=trimstats.txt ;
done


## call the genome (using the curated GCF version)

## index genome
samtools faidx $genome
bwa-mem2 index $genome


#echo 'export PATH=$PATH:/home/panas/Downloads/gatk-4.2.2.0' >> ~/.bashrc

for i in "$batch_number"/Clean/*_clean_1.fastq.gz
	do dname=$(dirname "${i}"); name=$(basename "${i}" _clean_1.fastq.gz)
	echo "processing $name"
   
	in1=${dname}/${name}_clean_1.fastq.gz
	in2=${dname}/${name}_clean_2.fastq.gz
	bam=$batch_number/Aligned/${name}_aligned.bam
	sorted_bam=$batch_number/Aligned/${name}_aligned_sorted.bam
	rmd_bam=$batch_number/Aligned/${name}_rmd.bam
		  
	bwa-mem2 mem -t 110 -R "@RG\tID:${name}\tSM:${name}\tPL:illumina\tLB:$batch_number\tPU:$batch_number" $genome $in1 $in2 | samtools view -@ 110 -bSu - > $bam
    samtools sort -@ 110 -o $sorted_bam $bam
	samtools index -@ 110 $sorted_bam    
	gatk --java-options "-Xmx80G" MarkDuplicates I=$sorted_bam O=$rmd_bam REMOVE_DUPLICATES=true M=$batch_number/Aligned/"${name}".duplicates.txt
	samtools index -@ 110 $rmd_bam

	rm -f "$bam"
	rm -f "$sorted_bam"
	rm -f "${sorted_bam}".bai
done 


## get the total and the mapped only number of reads of a BAM file 
## legend
# samtools view -c = count reads, -F flag (filter out) reads: https://broadinstitute.github.io/picard/explain-flags.html
# -F 260 = filter out unmapped and secondary aligned reads
# -F 256 = filter out only secondary aligned reads


for i in "$batch_number"/Aligned/*_rmd.bam;
do
	base_name=$(basename "$i" _rmd.bam)
	count_total=$(samtools view -c "$i")
	count_mapped=$(samtools view -c -F 260 "$i")
	echo "$base_name, $count_total, $count_mapped" >> "$prop_aligned"
done
