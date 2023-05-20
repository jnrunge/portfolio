This folder's notebooks get you from the engram-stored raw reads to 

* phased founder + f1 genotypes
* some qc metrics
    * Phasing success
    * Call comparison with previous standard VCF
* ancestryHMM input file (founder haplotypes)


The way it works is that you only change information in the `config.R` file and the 2nd sheet in the `1bxyye65_Q50ma6H8Avmh4aF6xT_eDk8PSwyiMUZmkwQ` Google drive document. Since the founders are done with sequencing, the Google drive document can likely be left alone. The scripts are run consecutively according to their number and expect everything to be done in the previous document. If you only need to add another run, you need to either do everything again or you still need to have each run's individual bam file (e.g. `SW_1.1.bam`) in the bam sub-folder. These files are however deleted within the pipeline to save space, so you cannot run 03_MergeAndCall and then add more data later on. So best to just run it all again to be sure (takes about 3 days).