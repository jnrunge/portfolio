bcftools query -f "%CHROM\t%POS\t%INFO/DP\n" $1 | gzip > $1.DP.txt.gz