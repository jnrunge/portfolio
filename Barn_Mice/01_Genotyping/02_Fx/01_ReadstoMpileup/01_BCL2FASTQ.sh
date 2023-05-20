module load bcl2fastq2/v2.20

# args: inputdir, outputdir, samplesheet, corestotal, cores_read/write

bcl2fastq --runfolder-dir $1 \
--output-dir $2 \
--sample-sheet \
$3 \
--barcode-mismatches 0 \
-p $4 -w $5 -r $5