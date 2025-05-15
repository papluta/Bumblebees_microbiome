#!/bin/bash

## Start calling for the chromosomes that are already done

# Directory with bcf files
BCF_DIR="./bcf_pascuorum"
VCF_DIR="./vcf_pascuorum"
int=600  # in seconds!

# Loop forever
while true; do

    for bcf_file in "$BCF_DIR"/*.bcf; do
        [ -e "$bcf_file" ] || continue  # skip if no files match

        base_name=$(basename raw_mpileup_${bcf_file}.bcf)
        output_file="$VCF_DIR/${base_name}_raw_snps.vcf.gz"

        # Skip if call file already exists
        if [ -f "$output_file" ]; then
            continue
        fi

        # size must be stable for 1 minute
        size1=$(stat --format="%s" "$bcf_file")
        sleep 60
        size2=$(stat --format="%s" "$bcf_file")

        if [ "$size1" -eq "$size2" ]; then
            echo "Calling for $bcf_file"

            bcftools call -mv -v -Oz --threads 80 -o "$output_file" "$bcf_file"

            if [ $? -eq 0 ]; then
                echo "bcftools call finished for $base_name"
            else
                echo "Error running bcftools call on $bcf_file"
            fi
        else
            echo "$bcf_file nothing new"
        fi
    done

    echo "10-minute nap"
    sleep "$int"
done
