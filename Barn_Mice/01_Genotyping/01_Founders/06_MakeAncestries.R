library(dplyr)
library(data.table)
library(doParallel)
library(foreach)
library(stringr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

mpileups=readLines(args[1])
haps_input_file=args[2]
samples_input_file=args[3]
haps_file=args[4]
cores=as.numeric(args[5])

stringr_count_bases = function(string_column,ref_column)
    {
    #if(nchar(string) > 0 & nchar(ref) > 0)
     #   {
        

        a=str_count(string_column, fixed("a"))
        a=a+str_count(string_column, fixed("A"))

        t=str_count(string_column, fixed("t"))
        t=t+str_count(string_column, fixed("T"))

        g=str_count(string_column, fixed("g"))
        g=g+str_count(string_column, fixed("G"))

        c=str_count(string_column, fixed("c"))
        c=c+str_count(string_column, fixed("C"))
    
        ref=str_count(string_column, fixed(","))
        ref=ref+str_count(string_column, fixed("."))

        a = a + (ref * (ref_column == 'a'))
        a = a + (ref * (ref_column == 'A'))
    
        t = t + (ref * (ref_column == 't'))
        t = t + (ref * (ref_column == 'T'))
    
        g = g + (ref * (ref_column == 'g'))
        g = g + (ref * (ref_column == 'G'))
    
        c = c + (ref * (ref_column == 'c'))
        c = c + (ref * (ref_column == 'C'))


        return(cbind(t,a,g,c))
      #  }
   # else
     #   {
       # return(cbind(NA,NA,NA,NA))
  #  }
    
}

ascii_score_to_numeric = function(ascii)
    {
    return(mean(utf8ToInt(ascii)-33))
}

ascii_score_to_filter = function(ascii)
    {
    return((utf8ToInt(ascii)-33) >= 30)
}


cl <- makeCluster(cores)
registerDoParallel(cl)


print(Sys.time())
library(stringr)

# I make a nice and readable df out of the founder mpileups

all_mpileup_df=foreach(mpileup=mpileups, .combine="rbind", .packages=c("data.table","stringr")) %dopar%
    {
    #gc()
    #print(mpileup)
    mpileup_df=fread(cmd = paste0("xzcat ",mpileup), header = FALSE)
    #print(Sys.time())
    if(nrow(mpileup_df) > 0)
        {
        mpileup_df$A = 0
        mpileup_df$T = 0
        mpileup_df$G = 0
        mpileup_df$C = 0

        #mpileup_df = mpileup_df[1:100,]

        #print(unique(mpileup_df$V5))
        #print(unique(mpileup_df$V3))

        mpileup_df_2=stringr_count_bases(mpileup_df$V5, mpileup_df$V3)

        mpileup_df_2 = as.data.frame(mpileup_df_2)

        mpileup_df$A = mpileup_df_2$a
        mpileup_df$T = mpileup_df_2$t
        mpileup_df$G = mpileup_df_2$g
        mpileup_df$C = mpileup_df_2$c

        mpileup_df$Individual = strsplit(basename(mpileup),".",fixed=TRUE)[[1]][1]


        return(mpileup_df)
    }
}
all_mpileup_df$A = as.character(all_mpileup_df$A)
all_mpileup_df$T = as.character(all_mpileup_df$T)
all_mpileup_df$G = as.character(all_mpileup_df$G)
all_mpileup_df$C = as.character(all_mpileup_df$C)
print(Sys.time())


# combine with haps file

haps = fread(file = haps_input_file)

samples=fread(file = samples_input_file)

samples=samples[-1,-3]

colnames(haps)[6:ncol(haps)] = c(rbind(samples$ID_1,paste(samples$ID_2, "_",sep="")))
haps_1 = haps
haps_2 = haps

haps_1[,6:ncol(haps)] = "0"
haps_2[,6:ncol(haps)] = "0"

all_mpileup_df$pos_total=paste(all_mpileup_df$V1,all_mpileup_df$V2, sep = " ")

mpileup_rows = which(paste("SM_",all_mpileup_df$Individual, sep="") == samples$ID_1[1])

fill_mpileup_df = expand.grid(V1 = NA, V2 = NA, pos_total = unique(all_mpileup_df$pos_total), Individual = unique(all_mpileup_df$Individual))

fill_mpileup_df=tidyr::separate(fill_mpileup_df, col = pos_total, into = c("V1","V2"), sep = " ")
fill_mpileup_df=fill_mpileup_df[,c(2,3,1)]
fill_mpileup_df$V2 = as.numeric(fill_mpileup_df$V2)

fill_mpileup_df = dplyr::left_join(x=fill_mpileup_df, y=all_mpileup_df, by = c("V1","V2","Individual"))

remove(all_mpileup_df)
gc()

fill_mpileup_df$A[is.na(fill_mpileup_df$A)] = 0
fill_mpileup_df$C[is.na(fill_mpileup_df$C)] = 0
fill_mpileup_df$T[is.na(fill_mpileup_df$T)] = 0
fill_mpileup_df$G[is.na(fill_mpileup_df$G)] = 0

haps=haps[paste(haps$V1,haps$V3) %in% paste(fill_mpileup_df$V1, fill_mpileup_df$V2),]
#haps_1=haps_1[paste(haps_1$V1,haps_1$V3) %in% paste(fill_mpileup_df$V1, fill_mpileup_df$V2),]
#haps_2=haps_2[paste(haps_2$V1,haps_2$V3) %in% paste(fill_mpileup_df$V1, fill_mpileup_df$V2),]
fill_mpileup_df=fill_mpileup_df[paste(fill_mpileup_df$V1, fill_mpileup_df$V2) %in% paste(haps$V1,haps$V3),]

fill_mpileup_df=fill_mpileup_df[order(paste(fill_mpileup_df$V1, fill_mpileup_df$V2)),]
haps=haps[order(paste(haps$V1,haps$V3)),]
#haps_1=haps_1[order(paste(haps_1$V1,haps_1$V3)),]
#haps_2=haps[order(paste(haps_2$V1,haps_2$V3)),]

if(print(nrow(fill_mpileup_df)/length(unique(fill_mpileup_df$Individual))) != print(nrow(haps)))
    {
    stop("haps and fill_mpileup_df lengths don't correspond. Aborting...")
}

#print(nrow(haps_1))
#print(nrow(haps_2))

sample_cols = which(colnames(haps) == samples$ID_1[1] | colnames(haps) == paste(samples$ID_1[1],"_",sep=""))

mpileup_rows = which(paste("SM_",fill_mpileup_df$Individual, sep="") == samples$ID_1[1])



haps_all = haps[,1:5]
fillup_na=(nrow(haps))
haps_all_ext = foreach(sample=samples$ID_1, .combine = "cbind") %dopar%
    {
        #print(sample)
        sample_cols = which(colnames(haps) == sample | colnames(haps) == paste(sample,"_",sep=""))
        #print(sample_cols[1])
        #print(sample_cols[2])
        if (sample %in% paste("SM_",fill_mpileup_df$Individual, sep=""))

            {

            mpileup_rows = which(paste("SM_",fill_mpileup_df$Individual, sep="") == sample)

            hap1_count_ref = 
            (as.numeric(fill_mpileup_df$A[mpileup_rows]) * as.numeric(haps$V4 == "A") * as.numeric(haps[[colnames(haps)[sample_cols[1]]]] == "0"))+
            (as.numeric(fill_mpileup_df$C[mpileup_rows]) * as.numeric(haps$V4 == "C") * as.numeric(haps[[colnames(haps)[sample_cols[1]]]] == "0"))+
            (as.numeric(fill_mpileup_df$G[mpileup_rows]) * as.numeric(haps$V4 == "G") * as.numeric(haps[[colnames(haps)[sample_cols[1]]]] == "0"))+
            (as.numeric(fill_mpileup_df$T[mpileup_rows]) * as.numeric(haps$V4 == "T") * as.numeric(haps[[colnames(haps)[sample_cols[1]]]] == "0"))

            hap1_count_alt = 
            (as.numeric(fill_mpileup_df$A[mpileup_rows]) * as.numeric(haps$V5 == "A") * as.numeric(haps[[colnames(haps)[sample_cols[1]]]] == "1"))+
            (as.numeric(fill_mpileup_df$C[mpileup_rows]) * as.numeric(haps$V5 == "C") * as.numeric(haps[[colnames(haps)[sample_cols[1]]]] == "1"))+
            (as.numeric(fill_mpileup_df$G[mpileup_rows]) * as.numeric(haps$V5 == "G") * as.numeric(haps[[colnames(haps)[sample_cols[1]]]] == "1"))+
            (as.numeric(fill_mpileup_df$T[mpileup_rows]) * as.numeric(haps$V5 == "T") * as.numeric(haps[[colnames(haps)[sample_cols[1]]]] == "1"))

            hap2_count_ref = 
            (as.numeric(fill_mpileup_df$A[mpileup_rows]) * as.numeric(haps$V4 == "A") * as.numeric(haps[[colnames(haps)[sample_cols[2]]]] == "0"))+
            (as.numeric(fill_mpileup_df$C[mpileup_rows]) * as.numeric(haps$V4 == "C") * as.numeric(haps[[colnames(haps)[sample_cols[2]]]] == "0"))+
            (as.numeric(fill_mpileup_df$G[mpileup_rows]) * as.numeric(haps$V4 == "G") * as.numeric(haps[[colnames(haps)[sample_cols[2]]]] == "0"))+
            (as.numeric(fill_mpileup_df$T[mpileup_rows]) * as.numeric(haps$V4 == "T") * as.numeric(haps[[colnames(haps)[sample_cols[2]]]] == "0"))

            hap2_count_alt = 
            (as.numeric(fill_mpileup_df$A[mpileup_rows]) * as.numeric(haps$V5 == "A") * as.numeric(haps[[colnames(haps)[sample_cols[2]]]] == "1"))+
            (as.numeric(fill_mpileup_df$C[mpileup_rows]) * as.numeric(haps$V5 == "C") * as.numeric(haps[[colnames(haps)[sample_cols[2]]]] == "1"))+
            (as.numeric(fill_mpileup_df$G[mpileup_rows]) * as.numeric(haps$V5 == "G") * as.numeric(haps[[colnames(haps)[sample_cols[2]]]] == "1"))+
            (as.numeric(fill_mpileup_df$T[mpileup_rows]) * as.numeric(haps$V5 == "T") * as.numeric(haps[[colnames(haps)[sample_cols[2]]]] == "1"))

            # add in unphased het sites

            hap1_count_ref = hap1_count_ref + 
            (as.numeric(fill_mpileup_df$A[mpileup_rows]) * as.numeric(haps$V4 == "A") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[1]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$C[mpileup_rows]) * as.numeric(haps$V4 == "C") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[1]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$G[mpileup_rows]) * as.numeric(haps$V4 == "G") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[1]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$T[mpileup_rows]) * as.numeric(haps$V4 == "T") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[1]]]], fixed = TRUE)))

            hap1_count_alt = hap1_count_alt +
            (as.numeric(fill_mpileup_df$A[mpileup_rows]) * as.numeric(haps$V5 == "A") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[1]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$C[mpileup_rows]) * as.numeric(haps$V5 == "C") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[1]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$G[mpileup_rows]) * as.numeric(haps$V5 == "G") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[1]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$T[mpileup_rows]) * as.numeric(haps$V5 == "T") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[1]]]], fixed = TRUE)))

            hap2_count_ref = hap2_count_ref + 
            (as.numeric(fill_mpileup_df$A[mpileup_rows]) * as.numeric(haps$V4 == "A") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[2]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$C[mpileup_rows]) * as.numeric(haps$V4 == "C") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[2]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$G[mpileup_rows]) * as.numeric(haps$V4 == "G") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[2]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$T[mpileup_rows]) * as.numeric(haps$V4 == "T") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[2]]]], fixed = TRUE)))

            hap2_count_alt = hap2_count_alt +
            (as.numeric(fill_mpileup_df$A[mpileup_rows]) * as.numeric(haps$V5 == "A") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[2]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$C[mpileup_rows]) * as.numeric(haps$V5 == "C") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[2]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$G[mpileup_rows]) * as.numeric(haps$V5 == "G") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[2]]]], fixed = TRUE)))+
            (as.numeric(fill_mpileup_df$T[mpileup_rows]) * as.numeric(haps$V5 == "T") * as.numeric(grepl("*", haps[[colnames(haps)[sample_cols[2]]]], fixed = TRUE)))


            return(cbind(hap1_count_ref, hap1_count_alt, hap2_count_ref, hap2_count_alt))
        }
        else
            {
            return(cbind(rep(NA,fillup_na),rep(NA,fillup_na),rep(NA,fillup_na),rep(NA,fillup_na)))

        }

}
haps_all=cbind(haps_all, haps_all_ext)
remove(haps_all_ext)
gc()


haps_all=haps_all[,c(-2)]

count_ancs=(ncol(haps_all)-4)/4

colnames(haps_all) = c("chr","bp","ref","alt",paste(rep(1:12, each=4), "_", c(1,1,2,2), c("a","b"), sep=""))

colnames(haps_all)

haps_all$cM = 0

haps_all$bp = as.numeric(haps_all$bp)

setorder(haps_all, chr, bp)

fwrite(x = haps_all, file = haps_file, sep="\t",row.names = FALSE, col.names = TRUE, quote = FALSE, nThread=cores)


system(command = paste('zcat ',haps_file,' | sed -e "s/NC_000067\\.6/1/g" | sed -e "s/NC_000068\\.7/2/g" | sed -e "s/NC_000069\\.6/3/g" | sed -e "s/NC_000070\\.6/4/g" | sed -e "s/NC_000071\\.6/5/g" | sed -e "s/NC_000072\\.6/6/g" | sed -e "s/NC_000073\\.6/7/g" | sed -e "s/NC_000074\\.6/8/g" | sed -e "s/NC_000075\\.6/9/g" | sed -e "s/NC_000076\\.6/10/g" | sed -e "s/NC_000077\\.6/11/g" | sed -e "s/NC_000078\\.6/12/g" | sed -e "s/NC_000079\\.6/13/g" | sed -e "s/NC_000080\\.6/14/g" | sed -e "s/NC_000081\\.6/15/g" | sed -e "s/NC_000082\\.6/16/g" | sed -e "s/NC_000083\\.6/17/g" | sed -e "s/NC_000084\\.6/18/g" | sed -e "s/NC_000085\\.6/19/g" | sed -e "s/NC_000086\\.7/X/g" | sed -e "s/NC_000087\\.7/Y/g" | sed -e "s/NC_005089\\.1/MT/g" | gzip > ',haps_file,'.fixChrNames.gz',sep=""), intern = TRUE)

file.remove(haps_file)

file.rename(paste0(haps_file,'.fixChrNames.gz'), haps_file)