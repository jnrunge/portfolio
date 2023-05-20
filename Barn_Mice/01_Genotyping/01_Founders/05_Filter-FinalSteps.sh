# select founders only
# no mendelian errors
# only biallelic and zero missing genotypes

bcftools view -s SM_SW_1,SM_SW_2,SM_SW_3,SM_SW_4,SM_SW_5,SM_SW_6,SM_SW_7,SM_SW_8,SM_SW_9,SM_SW_10,SM_SW_11,SM_SW_12 -i "INFO/MERR=0" $1 | bcftools view -i "F_MISSING=0" -m2 -M2 -Oz -o $1.F0.0MERR.0Missing.m2M2.vcf.gz 



