#!/bin/sh

## move .bam files from all batches into one folder

mkdir bam_pascuorum
mv batch2/Aligned/*_rmd.bam* bam_pascuorum/
mv batch3/Aligned/*_rmd.bam* bam_pascuorum/

## check the proportion of mapped reads in each file
output_file="counts_summary.txt"

for i in bam_pascuorum/PA1-1_rmd.bam;
do
	base_name=$(basename "$i" _rmd.bam)
	count_total=$(samtools view -c "$i")
	count_mapped=$(samtools view -c -F 260 "$i")
    prop_mapped=$( echo "scale=2; $count_mapped / $count_total * 100" | bc)
	echo "$base_name, $count_total, $count_mapped" >> "$output_file"
done

ls bam_pascuorum/*_rmd.bam > all_bams.list


# Run mpileup per chromosome using parallel 
cat pascuorum_chromosomes.txt | parallel -j 40 \
    "bcftools mpileup -Ou -f ../reference_genomes/Bombus_pascuorum/GCF_905332965.1_iyBomPasc1.1_genomic.fna \
    -b all_bams.list \
    -r {} -q 20 -Q 30 -C 50 -d 250 -a FORMAT/DP,FORMAT/AD \
    -o bcf_pascuorum/raw_mpileup_{}.bcf"
    
## and run variant calling separately for each chromosome
for i in $(cat genome_chromosomes2.txt)
do name=$(basename ${i})
    echo "processing $name"
bcftools call -mv -v -Oz --threads 80 -o all_batches/vcf/${name}_raw_snps.vcf.gz all_batches/vcf/raw_mpileup_${name}.bcf
done

## merge all vcfs from different chromosomes
ls -1 *_raw_snps.vcf.gz > vcf_list.txt
bcftools concat -Oz vcf_list.txt -o Combee_BP_raw_snp.vcf.gz
bcftools index Combee_BP_raw_snp.vcf.gz

## check summary information stats from vcf file
bcftools stats  -s - Combee_BP_raw_snp.vcf.gz > Goe_raw_Combee_BP_raw_snp.sumstats