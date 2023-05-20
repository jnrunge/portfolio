columbia_username="REDACTED"
env_jupyter="my_conda4" # jupyter env with all the R requirements (same as notebook)
env_mapping_etc="samtools-116" # environment with samtools etc
slurm_acc="REDACTED" 

ref_file="REDACTED/GCF_000001635.26_GRCm38.p6_genomic.fna"

# founders
folder_for_sequences="REDACTED" # where things should be put -- empty before!

# Fx
folder_for_sequences_fx="REDACTED" # where things should be put -- empty before!


Barn_Mice_dir="REDACTED" # where the scripts / notebooks are
jnr_general_scripts_dir="REDACTED"  # https://github.com/jnrunge/general


# dont change below

if(!endsWith(folder_for_sequences, "/")){
    folder_for_sequences=paste0(folder_for_sequences,"/")
}

if(!endsWith(folder_for_sequences_fx, "/")){
    folder_for_sequences_fx=paste0(folder_for_sequences_fx,"/")
}

if(!endsWith(Barn_Mice_dir, "/")){
    Barn_Mice_dir=paste0(Barn_Mice_dir,"/")
}

if(!endsWith(jnr_general_scripts_dir, "/")){
    jnr_general_scripts_dir=paste0(jnr_general_scripts_dir,"/")
}

trim_dir=paste0(folder_for_sequences, "trim_galore")
bam_dir=paste0(folder_for_sequences, "bam")



library(googlesheets4)
gs4_auth(email = paste0(columbia_username,"@REDACTED"),cache="~/gcache")
files_df=read_sheet("REDACTED", sheet = 2)
files_df_fx=read_sheet("REDACTED", sheet = 3)


source(paste0(jnr_general_scripts_dir,"functions.R"))