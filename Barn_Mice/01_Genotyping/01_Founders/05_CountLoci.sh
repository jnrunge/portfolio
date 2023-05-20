#!/bin/bash
rm -f $1.counts.txt
# Read list of VCF files
while read -r vcf_file; do
  # Count number of loci in VCF file
  count=$(bcftools view -H "${vcf_file}" | wc -l)
  # Append count to input file
  echo $count >> $1.counts.txt
done < $1

paste $1 $1.counts.txt > $1.counts.txt.tmp

mv -f $1.counts.txt.tmp $1.counts.txt