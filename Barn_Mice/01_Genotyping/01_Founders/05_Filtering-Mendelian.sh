bcftools +mendelian $1 -T $2 -m a -Oz -o $1.MendelianAnnotated.vcf.gz
bcftools index $1.MendelianAnnotated.vcf.gz
bcftools query -f "%INFO/MERR\n" $1.MendelianAnnotated.vcf.gz | gzip > $1.MendelianAnnotated.vcf.gz.MERR.txt.gz