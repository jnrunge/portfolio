#!/bin/bash

# the point is to later quickly compare x, y, and autosome coverage

# Input BAM file
bamfile=$1

# Output file
outfile=$2

# Get the list of chromosomes
chromosomes=$(samtools idxstats $bamfile | grep "^N" | cut -f 1 | sort)

# Write the header to the output file
echo -e "Chromosome\tNum_MQ_30_or_greater" > $outfile

# Loop over each chromosome
for chrom in $chromosomes; do
    # Count the number of reads with MQ >= 30
    num_MQ_30_or_greater=$(samtools view -c -q 30 -F 0x100 $bamfile $chrom)

    # Write the result to the output file
    echo -e "$chrom\t$num_MQ_30_or_greater" >> $outfile
done

gzip $outfile