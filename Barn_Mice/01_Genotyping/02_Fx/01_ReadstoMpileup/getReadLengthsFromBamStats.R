# from my own https://github.com/jnrunge/BrusselSprouts/blob/main/scripts/getReadLengthsFromBamStats.R

args = commandArgs(trailingOnly=TRUE)
library(data.table)

bamstats=args[1]

df<-data.table::fread(cmd=paste("cat ",bamstats, " | grep ^RL | cut -f 2-",sep=""), data.table=FALSE)
df
max_rl<-as.numeric(df$V1[df$V2 == max(df$V2)])
max_rl<-max_rl[max_rl==max(max_rl,na.rm=TRUE)]
writeLines(as.character(max_rl), con = paste(bamstats,"-rl.txt",sep=""))