library("protti")
library("r3dmol")
library("Biostrings")
library("bioseq")
library(GenomicAlignments)             
library(UniprotR)
library(dplyr)
library(magrittr)
library(stringr)
library(tidyr)
library(ggplot2)

ecoli <- read.table("/Users/kamilkawiecki/Documents/Bioinformatics/Accession_numbers.txt")
class("ecoli")

ecolichar<- ecoli$V1 

ecoli[ , 1]

ecoligo <- GetProteinGOInfo(ecolichar)

PlotGoInfo(ecoligo)

PlotGOAll( EGO<- ecoligo,Top = 10,directorypath = getwd(),width = 8,height = 8)

ecolipatho <- GetPathology_Biotech(ecolichar)
print(ecolipatho)

Get.diseases(ecolipatho)


ecoliprot<- fetch_uniprot(ecolichar)
ecoliids <- unique(ecolichar)
ecoliinfo <- fetch_uniprot(ecoliids,columns = c("sequence","xref_pdb"))

ptsi_pgk_pdb_ids <- ecoliinfo %>% 
  mutate(pdb_id = str_split(xref_pdb, pattern = ";")) %>% 
  unnest(pdb_id) %>% 
  filter(pdb_id != "") 

ecolipdbid <- fetch_pdb(pdb_ids = unique(ptsi_pgk_pdb_ids$pdb_id))



ecolialph<- fetch_alphafold_prediction(ecolichar,return_data_frame = TRUE)
head(ecolialph,n=14)


