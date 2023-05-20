# input: pe-file1 pe-file2 reference individualID count-of-number-of-time-this-individual-will-be-in-this-folder

bwa-mem2 mem -t 2 -r 1.5 -E 1 -w 100 -T 0 -a $3 $1 $2 > $4.$5.sam # same as bwa mem but faster

samtools sort -n -o $4.$5.sam.sortn.bam $4.$5.sam
rm -f $4.$5.sam

samtools fixmate -m $4.$5.sam.sortn.bam $4.$5.sam.sortn.fixmate.bam
rm -f $4.$5.sam.sortn.bam

samtools sort -o $4.$5.sam.sortn.fixmate.sortc.bam $4.$5.sam.sortn.fixmate.bam
rm -f $4.$5.sam.sortn.fixmate.bam

samtools markdup -l 150 $4.$5.sam.sortn.fixmate.sortc.bam $4.$5.bam
rm -f $4.$5.sam.sortn.fixmate.sortc.bam

samtools addreplacerg -r ID:$4.$5 -r SM:SM_$4 -o $4.$5.rg.bam $4.$5.bam
mv -f $4.$5.rg.bam $4.$5.bam

samtools index $4.$5.bam

# delete pe-files
rm -f $1 $2