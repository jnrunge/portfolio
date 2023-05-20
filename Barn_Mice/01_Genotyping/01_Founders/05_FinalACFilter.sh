bcftools view -i "AC != AN & AC != 0" -Oz -o Founders.filtered.var.vcf.gz Founders.filtered.vcf.gz
bcftools index Founders.filtered.var.vcf.gz