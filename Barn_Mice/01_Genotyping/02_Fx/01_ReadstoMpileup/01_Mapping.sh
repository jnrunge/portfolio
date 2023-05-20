# input: jupyter_env samtools_env reference individualID fastQ1 fastQ2
# reference / fastq are full paths

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

bwa-mem2 mem -t 1 -r 1.5 -E 1 -w 100 -T 0 -a $3 $5 $6 > $4.sam # same as bwa mem but faster

samtools sort -n -o $4.sam.sortn.bam $4.sam
rm -f $4.sam

samtools fixmate -m $4.sam.sortn.bam $4.sam.sortn.fixmate.bam
rm -f $4.sam.sortn.bam

samtools sort -o $4.sam.sortn.fixmate.sortc.bam $4.sam.sortn.fixmate.bam
rm -f $4.sam.sortn.fixmate.bam

samtools stats $4.sam.sortn.fixmate.sortc.bam > $4.sam.sortn.fixmate.sortc.bam.stats

. ${SCRIPT_DIR}/../../../activateEnv.sh $1
Rscript ${SCRIPT_DIR}/getReadLengthsFromBamStats.R $4.sam.sortn.fixmate.sortc.bam.stats
. ${SCRIPT_DIR}/../../../activateEnv.sh $2
rl=`cat $4.sam.sortn.fixmate.sortc.bam.stats-rl.txt`
echo Read Lengths detected as $rl

rm -f $4.bam
rm -f $4.bam.bai

samtools markdup -l $rl $4.sam.sortn.fixmate.sortc.bam $4.bam
samtools index $4.bam
rm -f $4.sam.sortn.fixmate.sortc.bam

rm -f $4.sam.sortn.fixmate.sortc.bam.stats-rl.txt $4.sam.sortn.fixmate.sortc.bam.stats $5 $6