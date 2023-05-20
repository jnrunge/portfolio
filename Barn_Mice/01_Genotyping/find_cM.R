args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=5) {
  stop("Exactly 5 arguments must be supplied (Chr, Region, VCF, outfile, genmap)", call.=FALSE)
}

chromosome = args[1] # vcf format
region = args[2] # numeric format
vcf = args[3]
outfile = args[4]

genmap_raw = read.csv(args[5]) # ColumbiaProjects/Barn_Mice/XX_Data/avg.map.csv.gz
library(data.table)
vcf_skip = system(command = paste("zcat ",vcf," | grep -o '^##' | wc -l", sep = ""), intern = TRUE)
vcf_read = fread(vcf, skip = as.numeric(vcf_skip), sep="\t")



vcf_read = subset(vcf_read, `#CHROM` == region, select=c("POS"))
genmap_raw = subset(genmap_raw, Chr == chromosome)



find_cM = function(pos)
    {
    next_pos_2nd="unset"
    next_pos=which(abs(pos-genmap_raw$Pos) == min(abs(pos-genmap_raw$Pos)))[1]
    if(next_pos != 1 & genmap_raw$Pos[next_pos] > pos)
    {next_pos_2nd = next_pos-1}
    if(next_pos != nrow(genmap_raw) & genmap_raw$Pos[next_pos] < pos)
    {next_pos_2nd = next_pos+1}
    if(next_pos_2nd == "unset")
    {next_pos_2nd = which(genmap_raw$Pos == genmap_raw$Pos[-next_pos][which(abs(pos-genmap_raw$Pos[-next_pos]) == min(abs(pos-genmap_raw$Pos[-next_pos])))[1]])
    }
    cM_Mb = NA
    cM = NA
    
    #print(genmap_raw$Pos[next_pos])
    #print(genmap_raw$Pos[next_pos_2nd])
    
    if(next_pos < next_pos_2nd & genmap_raw$Pos[next_pos] <= pos)
        {
        cM_b = (genmap_raw$cM[next_pos_2nd]-genmap_raw$cM[next_pos])/(genmap_raw$Pos[next_pos_2nd]-genmap_raw$Pos[next_pos])
        cM=genmap_raw$cM[next_pos] + cM_b * (pos-genmap_raw$Pos[next_pos])
        cM_Mb = cM_b*1000000
        #print("1")
    }
    
    if(next_pos > next_pos_2nd & genmap_raw$Pos[next_pos_2nd] <= pos)
        {
        cM_b = (genmap_raw$cM[next_pos]-genmap_raw$cM[next_pos_2nd])/(genmap_raw$Pos[next_pos]-genmap_raw$Pos[next_pos_2nd])
        cM=genmap_raw$cM[next_pos_2nd] + cM_b * (pos-genmap_raw$Pos[next_pos_2nd])
        cM_Mb = cM_b*1000000
        #print("2")
    }
    
    if(genmap_raw$Pos[next_pos] > pos & genmap_raw$Pos[next_pos_2nd] > pos)
        
        {
        cM_Mb = 0
        cM = 0
        #print("3")
    }
    
    if(next_pos > next_pos_2nd & genmap_raw$Pos[next_pos] <= pos)
        {
        
            cM_b = (genmap_raw$cM[next_pos]-genmap_raw$cM[next_pos_2nd])/(genmap_raw$Pos[next_pos]-genmap_raw$Pos[next_pos_2nd])
            cM=genmap_raw$cM[next_pos] + cM_b * (pos-genmap_raw$Pos[next_pos])
            cM_Mb = cM_b*1000000
        #print("4")
        }
    
     if(next_pos < next_pos_2nd & genmap_raw$Pos[next_pos_2nd] <= pos)
        {
        
            cM_b = (genmap_raw$cM[next_pos_2nd]-genmap_raw$cM[next_pos])/(genmap_raw$Pos[next_pos_2nd]-genmap_raw$Pos[next_pos])
            cM=genmap_raw$cM[next_pos_2nd] + cM_b * (pos-genmap_raw$Pos[next_pos_2nd])
            cM_Mb = cM_b*1000000
         #print("5")
        }
    
    #print(cM_b)
    return(cbind(pos,cM_Mb,cM))
}

library(foreach)

library(readr)

write_tsv(as.data.frame(foreach(pos=vcf_read$POS, .combine = "rbind") %do% find_cM(pos)), file = outfile, escape = NULL)

