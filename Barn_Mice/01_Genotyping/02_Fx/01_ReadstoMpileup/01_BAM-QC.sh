#!/bin/bash

cd $1 # fx head dir

bamfile=$2 # relative path with "run" folder included
outfile=$3 # given here, but should not really be changed

# Calculate coverage
coverage=$(samtools depth -a $bamfile | awk '{sum+=$3} END {print sum/NR}')

# Calculate mean mapping quality (MQ)
mean_MQ=$(samtools view -q 0 $bamfile | awk '{sum+=$5} END {if (NR>0) print sum/NR; else print 0}')
mean_MQwo0=$(samtools view -q 1 $bamfile | awk '{sum+=$5} END {if (NR>0) print sum/NR; else print 0}')

# Calculate % reads with max MQ (assuming 30 as max MQ)
total_reads=$(samtools view -c $bamfile)
max_MQ_reads=$(samtools view -c -q 30 $bamfile)
max_MQ_percent=$(bc -l <<< "100*$max_MQ_reads/$total_reads")

# Calculate mean insert size
mean_insert_size=$(samtools view -f 0x2 $bamfile | awk '{sum+=$9} END {if (NR>0) print sum/NR; else print 0}')

# Calculate mean read length
mean_read_length=$(samtools view -q 1 $bamfile | awk '{sum+=length($10)} END {if (NR>0) print sum/NR; else print 0}')

# Calculate % duplicated
total_reads=$(samtools view -F 0x900 -c $bamfile)
duplicate_reads=$(samtools view -f 0x400 -F 0x900 -c $bamfile)
duplication_percent=$(bc -l <<< "100*$duplicate_reads/$total_reads")

# Write the stats to the output file

# make sure there is no race condition

(
  flock -x 1337
  # Your command here, e.g.:
  current_time=$(date +%s)
    echo -e "$bamfile\t$current_time\t$coverage\t$mean_MQ\t$mean_MQwo0\t$max_MQ_percent\t$mean_insert_size\t$mean_read_length\t$duplication_percent" >> $outfile
    ) 1337>$outfile.lock

