# input: chrFasta chrNumeric VCF Genmap GenmapRaw samtoolswhatshap_env ped outdir reference bamfiles 

date

gotodir=`dirname "$0"`
cd $gotodir
pwd

Rscript ../find_cM.R $2 $1 $3 $4 $5

. ../../activateEnv.sh $6

cd $8

whatshap phase --ped $7 --chromosome $1 --genmap $4 -o $3.phased.vcf -r $9 $3 ${@:10}

bgzip -f $3.phased.vcf

bcftools index $3.phased.vcf.gz

date