args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(dplyr)

for(arg in args){
    if(arg==args[length(args)]){
        next
        }
    if(endsWith(arg, "gz"))
        {
            file=fread(cmd=paste("zcat ",arg,sep=""), data.table=FALSE, header=FALSE)

        }
    if(endsWith(arg, "xz"))
        {
            file=fread(cmd=paste("xzcat ",arg,sep=""), data.table=FALSE, header=FALSE)

        }
    if(arg == args[1]){
        file_merged=file
        }else{
        file_merged=bind_rows(file_merged,file)
        }
    }

file_merged=file_merged[order(file_merged$V1,file_merged$V2),]

#while(sum(duplicated(paste(file_merged$V1,file_merged$V2)) > 0)){
#    j=which(duplicated(paste(file_merged$V1,file_merged$V2)))[1]
#    i=which(file_merged$V1 == file_merged$V1[j] & file_merged$V2 == file_merged$V2[j])[1]
#    file_merged$V4[i]=file_merged$V4[i]+file_merged$V4[j]
#    file_merged$V5[i]=paste(file_merged$V5[i],file_merged$V5[j],sep="")
#    file_merged$V6[i]=paste(file_merged$V6[i],file_merged$V6[j],sep="")
#    file_merged=file_merged[-j,]
#}

while(sum(duplicated(paste(file_merged$V1,file_merged$V2)) > 0)){
    j=which(duplicated(paste(file_merged$V1,file_merged$V2)))
    file_merged_sub=file_merged[j,] #for removing >1 duplicates per repeat of this while
    j2=which(!duplicated(paste(file_merged_sub$V1,file_merged_sub$V2)))
    j=j[j2]
    i=j-1

    if(sum(file_merged$V1[i] == file_merged$V1[j] & file_merged$V2[i] == file_merged$V2[j]) != length(i)){
        stop("Order went wrong.")
    }
    file_merged$V4[i]=file_merged$V4[i]+file_merged$V4[j]
    file_merged$V5[i]=paste(file_merged$V5[i],file_merged$V5[j],sep="")
    file_merged$V6[i]=paste(file_merged$V6[i],file_merged$V6[j],sep="")
    file_merged=file_merged[-j,]
}

fwrite(file_merged, args[length(args)], sep="\t", col.names=FALSE)
print(system(command=paste("xz -vf ",args[length(args)],sep=""), intern=TRUE))
print(system(command=paste("xz -t ",args[length(args)],".xz 2> ",args[length(args)],".xz.t",sep=""), intern=TRUE))