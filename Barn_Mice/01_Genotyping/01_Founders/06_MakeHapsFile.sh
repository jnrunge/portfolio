cd $1

file=$2 # prefix before .vcf.gz

# change the notation of "unphased" homozygous sites so that for hapsample it is treated as phased
bcftools view $file.vcf.gz | sed -e 's/0\/0/0|0/g' | sed -e 's/1\/1/1|1/g' > $file.AllPhasedOrUnphasedHet.vcf

cat $file.AllPhasedOrUnphasedHet.vcf | bcftools convert --hapsample $file.AllPhasedOrUnphasedHet_hapsample

zcat $file.AllPhasedOrUnphasedHet_hapsample.hap.gz | awk '{print $1,$3}' | sed -e 's/\s/\t/g' > $file.posfile_for_mpileup

rm -f $file.AllPhasedOrUnphasedHet.vcf