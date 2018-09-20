library(tidyverse)

unichem <- function(drugList, input_source, output_source, input_source_name, 
                    output_source_name){
    mapping <- RCurl::getURL(paste0("https://www.ebi.ac.uk/unichem/rest/mapping/",
                                    input_source,"/",output_source))
    mapping <- jsonlite::fromJSON(mapping) %>% setNames(c(output_source_name, 
                                                          input_source_name))
    print(drugList[!drugList[,which(colnames(drugList) == colnames(mapping)[2])] 
                   %>% unlist %in% mapping[,2],])
    left_join(drugList, mapping, by = NULL) %>% unique
}

name2CID <- function(nm){
    library(RCurl)
    library(jsonlite)
    nm <- URLencode(nm)
    name2cid <- getURL(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", 
                            nm, "/cids/JSON", sep = ""))
    name2cid <- fromJSON(name2cid)
    return(name2cid$IdentifierList$CID)
}

convertIDs <- function(drugList){
    cids <- sapply(drugList, name2CID)
    cids <- cids[!sapply(cids, is.null)]
    cids <- reshape2::melt(cids) %>%
        rename(CID = value, drugName = L1)
    noCID <- setdiff(drugList, cids$drugName)
    print(noCID)
    return(cids)
}

cid2chembl <- function(in_fx, out_fx){
    pubchem2chembl <- read_tsv("./data/raw/src1src22.txt") %>% 
        setNames(c("ChEMBLID","CID")) %>%
        mutate(CID = as.character(CID))
    read_tsv(in_fx) %>%
        mutate(CID = as.character(CID)) %>%
        left_join(pubchem2chembl) %>%
        unique() %>%
        write_tsv(out_fx)
}

pubchem2chembl <- read_tsv("./data/raw/src1src22.txt") %>% 
    setNames(c("ChEMBLID","CID"))

snell2016 <- read_delim("./data/raw/Snell2016/drugList.csv", delim = ";")
snell2016 <- unichem(snell2016, 2, 1, "DrugBankID", "ChEMBLID")  
# 1 Gamma Hydroxybutyric Acid DB01440   
# 2 Nafarelin                 DB00666   
# 3 Colistin                  DB00803
snell2016$ChEMBLID[snell2016$drugName == "Gamma Hydroxybutyric Acid"] <- "CHEMBL1342"
snell2016$ChEMBLID[snell2016$drugName == "Nafarelin"] <- "CHEMBL1201309"
snell2016$ChEMBLID[snell2016$drugName == "Colistin"] <- "CHEMBL501505"
system('mkdir -p ./data/processed/Snell2016/')
write_tsv(snell2016, path = "./data/processed/Snell2016/drugList_CHEMBL.tsv")

snell2018 <- read_delim("./data/raw/Snell2018/drugList.csv", delim = ";")
snell2018 <- unichem(snell2018, 2, 1, "DrugBankID", "ChEMBLID")
# 1 colistin                  DB00803   
# 2 salicylate-sodium         DB01398   
# 3 ivermectin                DB00602   
# 4 gadopentetate dimeglumine DB00789 
snell2018$ChEMBLID[snell2018$drugName == "colistin"] <- "CHEMBL501505"
snell2018$ChEMBLID[snell2018$drugName == "salicylate-sodium"] <- "ChEMBL447868"
snell2018$ChEMBLID[snell2018$drugName == "ivermectin"] <- "CHEMBL341047"
snell2018$ChEMBLID[snell2018$drugName == "gadopentetate dimeglumine"] <- "CHEMBL1200431"
system('mkdir -p ./data/processed/Snell2018/')
write_tsv(snell2018, path = "./data/processed/Snell2018/drugList_CHEMBL.tsv")

mofidifar2018 <- read_delim("./data/raw/Mofidifar2018/drugList.csv", 
                            delim = ";")
mofidifar2018 <- unichem(mofidifar2018, 2, 1, "DrugBankID", "ChEMBLID")
system('mkdir -p ./data/processed/Mofidifar2018/')
write_tsv(mofidifar2018, 
          path = "./data/processed/Mofidifar2018/drugList_CHEMBL.tsv")

ziehm2017 <- read_delim("./data/raw/Ziehm2017/drugList.csv", delim = ";")
ziehm2017 <- unichem(ziehm2017, 3, 1, "PDBID", "ChEMBLID")
system('mkdir -p ./data/processed/Ziehm2017/')
write_tsv(ziehm2017, path = "./data/processed/Ziehm2017/drugList_CHEMBL.tsv")

fernandes2016 <- read_delim("./data/raw/Fernandes2016/drugList.csv", 
                            delim = ";")
fernandes2016 <- convertIDs(fernandes2016$drugName)
# [1] "FK‐228"    "CHR‐3996 " "GDC‐0068 " "MK‐2206 "
fernandes2016 <- inner_join(fernandes2016, pubchem2chembl) %>% unique()
system('mkdir -p ./data/processed/Fernandes2016/')
write_tsv(fernandes2016, 
          path = "./data/processed/Fernandes2016/drugList_CHEMBL.tsv")

fuentealba2018 <- read_delim("./data/raw/Fuentealba2018/drugList.csv", 
                             delim = ";")
fuentealba2018 <- unichem(fuentealba2018, 2, 1, "DrugBankID", "ChEMBLID")
system('mkdir -p ./data/processed/Fuentealba2018/')
write_tsv(fuentealba2018, 
          path = "./data/processed/Fuentealba2018/drugList_CHEMBL.tsv")

aliper2016 <- read_delim('./data/raw/Aliper2016/drugList2.csv', col_names = F, 
                         delim = '\t')[[1]]
aliper2016_cids <- convertIDs(aliper2016)
system('mkdir -p ./data/processed/Aliper2016/')
write_tsv(aliper2016_cids, path = './data/processed/Aliper2016/drugList_CID.tsv')

barardo2017 <- read_delim('./data/raw/Barardo2017/drugList.csv', col_names = F, 
                          delim = '\t')[[1]]
barardo2017_cids <- convertIDs(barardo2017)
system('mkdir -p ./data/processed/Barardo2017/')
write_tsv(barardo2017_cids, 
          path = './data/processed/Barardo2017/drugList_CID.tsv')

calvert2016 <- read_delim('./data/raw/Calvert2016/drugList.csv', col_names = F, 
                          delim = '\t')[[1]]
calvert2016_cids <- convertIDs(calvert2016)
system('mkdir -p ./data/processed/Calvert2016/')
write_tsv(calvert2016_cids, 
          path = './data/processed/Calvert2016/drugList_CID.tsv')

donertas2018 <- read_delim('./data/raw/Donertas2018/drugList.csv', 
                           col_names = F, delim = '\t')[[1]]
donertas2018_cids <- convertIDs(donertas2018)
# [1] "15-delta prostaglandin J2"
system('mkdir -p ./data/processed/Donertas2018/')
write_tsv(donertas2018_cids, 
          path = './data/processed/Donertas2018/drugList_CID.tsv')

yang2018 <- read_delim('./data/raw/Yang2018/drugList.csv', col_names = F, 
                       delim = '\t')[[1]]
yang2018_cids <- convertIDs(gsub('conjugated ','',yang2018))
# [1] "bis(2-ethylhexyl)"
system('mkdir -p ./data/processed/Yang2018/')
write_tsv(yang2018_cids, path = './data/processed/Yang2018/drugList_CID.tsv')

liu2016 <- read_csv('./data/raw/Liu2016/drugList.csv',col_names = T) %>%
    filter(F_ratio >= 2)
liu2016_cids <- data.frame(CID = gsub("(^0)\\1+", '', 
                                     substr(liu2016$PubChem_CID, 5, 12)), 
                          drugName = liu2016$PubChem_CID)
system('mkdir -p ./data/processed/Liu2016/')
write_tsv(liu2016_cids,path = './data/processed/Liu2016/drugList_CID.tsv')

snames <- c('Aliper2016', 'Barardo2017', 'Calvert2016', 'Donertas2018', 
           'Yang2018', 'Liu2016')
in_fs <- paste('./data/processed/', snames, '/drugList_CID.tsv', sep = '')
out_fs <- paste('./data/processed/', snames, '/drugList_CHEMBL.tsv', sep = '')

sapply(1:length(in_fs), function(i)cid2chembl(in_fs[i], out_fs[i]))
