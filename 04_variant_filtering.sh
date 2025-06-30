#!/bin/sh

DIR_V=/scratch/patrycja/COMBEE_lapidarius/vcf_lapidarius
DIR_F=/scratch/patrycja/COMBEE_lapidarius/vcf_filter_lapidarius
PREF="Combee_BL" # BP for pascuorum, BL for lapidarius

DIR_V=/scratch/patrycja/COMBEE_pascuorum/vcf_pascuorum
DIR_F=/scratch/patrycja/COMBEE_pascuorum/vcf_filter_pascuorum
PREF="Combee_BP" # BP for pascuorum, BL for lapidarius


# prefilter the dataset
vcftools --gzvcf $DIR_V/${PREF}_raw_snp_sorted.vcf.gz --max-missing 0.5 --minQ 30 --minDP 3 --recode --recode-INFO-all --out $DIR_F/${PREF}_f05q30dp3
bgzip -c $DIR_F/${PREF}_f05q30dp3.recode.vcf > $DIR_F/${PREF}_f05q30dp3.vcf.gz
bcftools index $DIR_F/${PREF}_raw_snp.vcf.gz


# Calculate relatedness, min. depth and missingness per individual
vcftools --vcf $DIR_F/${PREF}_f05q30dp3.recode.vcf --missing-indv --out $DIR_F/${PREF}_f05q30dp3
vcftools --vcf $DIR_F/${PREF}_f05q30dp3.recode.vcf --relatedness2 --out $DIR_F/${PREF}_f05q30dp3
vcftools --vcf $DIR_F/${PREF}_f05q30dp3.recode.vcf --depth --out $DIR_F/${PREF}_f05q30dp3

# remove poor quality and high missingness individuals
bcftools query -l $DIR_V/${PREF}_raw_snp.vcf.gz > all_samples.txt
grep -Fxv -f lapidarius_to_exclude_quality.txt all_samples.txt > samples_to_keep.txt
bcftools view -S samples_to_keep.txt -Oz -o $DIR_V/${PREF}_quality.vcf.gz $DIR_V/${PREF}_raw_snp.vcf.gz
bcftools index $DIR_V/${PREF}_quality.vcf.gz
echo "done!"

# 22M SNPs
vcftools --gzvcf $DIR_V/${PREF}_quality.vcf.gz --max-missing 0.95 --maf 0.05 --minQ 30 --minDP 4 --min-meanDP 5 \
    --max-meanDP 50 --remove-indels --recode --recode-INFO-all --out $DIR_F/${PREF}_g095maf005q30meandp5to50mindp4

# 15 M SNPs
vcftools --gzvcf $DIR_V/${PREF}_quality.vcf.gz --max-missing 0.99 --maf 0.05 --minQ 40 --minDP 5 --min-meanDP 5 \
    --max-meanDP 50 --remove-indels --recode --recode-INFO-all --out $DIR_F/${PREF}_g099maf005q40meandp5to50mindp5

# 358 674 SNPs for BP | 60 657 SNPs for BL
vcftools --gzvcf $DIR_V/${PREF}_quality.vcf.gz --max-missing 0.99 --maf 0.05 --minQ 40 --minDP 8 --min-meanDP 8 \
    --max-meanDP 80 --remove-indels --recode --recode-INFO-all --out $DIR_F/${PREF}_g099maf005q40meandp8to80mindp8

vcftools --vcf $DIR_F/${PREF}_g099maf005q40meandp8to80mindp8.recode.vcf --relatedness2 --out $DIR_F/${PREF}_g099maf005q40meandp8to80mindp8
vcftools --vcf $DIR_F/${PREF}_g099maf005q40meandp8to80mindp8.recode.vcf --missing-indv --out $DIR_F/${PREF}_g099maf005q40meandp8to80mindp8

# exclude related individuals with the R script
bcftools query -l $DIR_V/${PREF}_quality.vcf.gz > all_samples.txt
grep -Fxv -f lapidarius_to_exclude_related.txt all_samples.txt > samples_to_keep.txt
bcftools view -S samples_to_keep.txt -Oz -o $DIR_V/${PREF}_NR_quality.vcf.gz $DIR_V/${PREF}_raw_snp.vcf.gz
bcftools index $DIR_V/${PREF}_NR_quality.vcf.gz
echo "done!"

vcftools --gzvcf $DIR_V/${PREF}_NR_quality.vcf.gz --max-missing 0.99 --maf 0.05 --minQ 40 --minDP 8 --min-meanDP 8 \
    --max-meanDP 80 --remove-indels --recode --recode-INFO-all --out $DIR_F/${PREF}_NR_g099maf005q40meandp8to80mindp8

#Get HWE by snp to exclude highly heterozygotsity sites (errors)
vcftools --vcf $DIR_F/${PREF}_NR_g099maf005q40meandp8to80mindp8.recode.vcf --hardy --out $DIR_F/Goe_hardy

#run python script on both datasets
python ../../SnpsHe.py 

#remove snps identified with high heterozygosity
grep -Fwvf highHE_het60.txt $DIR_F/${PREF}_NR_g099maf005q40meandp8to80mindp8.recode.vcf > $DIR_F/${PREF}_NR_filtered_He60.recode.vcf

## pruning vcf files
# --set-missing-var-ids = assigns chromosome and position based IDs to the variants with missing IDs
# --indep-pairwise = window size (in kb), variant count to shift the window at the end of each step, paiwise r^2 threshold. Prunes the variants so that it excludes ones in linkage disequilibrium

plink --vcf $DIR_F/${PREF}_NR_filtered_He60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --indep-pairwise 50 5 0.5 --out $DIR_F/${PREF}_NR_filtered_He60_p1
plink --vcf $DIR_F/${PREF}_NR_filtered_He60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --extract $DIR_F/${PREF}_NR_filtered_He60_p1.prune.in --make-bed --out $DIR_F/${PREF}_NR_filtered_He60_p2
plink --bfile $DIR_F/${PREF}_NR_filtered_He60_p2 --recode vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --out $DIR_F/${PREF}_NR_filtered_He60_pruned

bgzip -c $DIR_F/${PREF}_NR_filtered_He60_pruned.vcf > $DIR_F/${PREF}_NR_filtered_He60_pruned.vcf.gz
bcftools index $DIR_F/${PREF}_NR_filtered_He60_pruned.vcf.gz