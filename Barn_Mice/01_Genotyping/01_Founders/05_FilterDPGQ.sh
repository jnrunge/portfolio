bcftools view -i 'INFO/DP>='$1' & INFO/DP<='$2'' $5 | bcftools +setGT -- -t q -n . -e 'FMT/GQ>='$3'' | bcftools +setGT -- -t q -n . -e 'FMT/DP>='$4'' | bgzip > $5.Cov$1to$2.GQ$3.DP$4.vcf.gz

bcftools index $5.Cov$1to$2.GQ$3.DP$4.vcf.gz