plate_to_labid=fread("~/columbia2019/firstplate_to_csv.tsv", data.table=FALSE)
plate_to_labid$Sample=str_replace(plate_to_labid$Sample, "m", "")
getLabID = function(x, plate_to_labid)
    {
    if(length(plate_to_labid$Sample[plate_to_labid$Final_ID == x | plate_to_labid$Final_ID_2 == x | plate_to_labid$Final_ID_3 == x | plate_to_labid$Final_ID_4 == x])==0){
        return(NA)
        }
    if(suppressWarnings(as.numeric(x) %in% 1:12))
        {
        return(as.character(x))
        }else
        {
        return(as.character(plate_to_labid$Sample[plate_to_labid$Final_ID == x | plate_to_labid$Final_ID_2 == x | plate_to_labid$Final_ID_3 == x | plate_to_labid$Final_ID_4 == x]))
        }
}