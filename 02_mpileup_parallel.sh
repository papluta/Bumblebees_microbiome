#!/bin/sh

## move .bam files from all batches into one folder

mkdir bam_pascuorum
mv batch2/Aligned/*_rmd.bam* bam_pascuorum/
mv batch3/Aligned/*_rmd.bam* bam_pascuorum/

ls bam_pascuorum/*_rmd.bam > all_bams.list


# Run mpileup per chromosome using parallel 
cat pascuorum_chromosomes.txt | parallel -j 40 \
    "bcftools mpileup -Ou -f ../reference_genomes/Bombus_pascuorum/GCF_905332965.1_iyBomPasc1.1_genomic.fna \
    -b all_bams.list \
    -r {} -q 20 -Q 30 -C 50 -d 250 -a FORMAT/DP,FORMAT/AD \
    -o bcf_pascuorum/raw_mpileup_{}.bcf"


## Use the next script to call variants as soon as the mpileup is done
    
# ## start calling for the chromosomes that are already done
# find . -type f -mmin +10 > pascuorum_chromosomes_done.txt
# cat ../pascuorum_chromosomes_done.txt | sed "s/^..raw_mpileup_//g" | sed "s/.bcf$//g" > ../pascuorum_chromosomes_done.txt

# for i in $(cat pascuorum_chromosomes_done.txt)
# do name=$(basename ${i})
#     echo "processing $name"
# bcftools call -mv -v -Oz --threads 80 -o vcf_pascuorum/${name}_raw_snps.vcf.gz bcf_pascuorum/raw_mpileup_${name}.bcf
# done

# ## merge all vcfs from different chromosomes
# ls -1 *_raw_snps.vcf.gz > vcf_list.txt
# bcftools concat -Oz vcf_list.txt -o Goe_raw_snp.vcf.gz
# tabix -p vcf Goe_raw_snp.vcf.gz

# ## check summary information stats from vcf file
# bcftools stats  -s - Goe_raw_snp.vcf.gz > Goe_raw_snp.sumstats

