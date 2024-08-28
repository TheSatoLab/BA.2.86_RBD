#!/usr/bin/env R

library(stringr)
library(readr)
library(ape)
library(dplyr)
library(data.table)
library(magrittr)

#use when bash
args = commandArgs(trailingOnly=T)

##########args##########
download_date <- args[1]
download.date <- as.Date(download_date)
out_prefix <- args[2]
metadata.name <- args[3]
nextclade.name <- args[4]
out.prefix <- args[5]

#output
tsv.merged.name <- paste(out_prefix,'.merged.final.tsv',sep="")
tsv.ID.name <- paste(out_prefix,'.accessionID.list.tsv',sep="")

###read all the needed files to merge GISAID and Nextclade metadata from the previous scripts
metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)

# i) only ‘original passage’ sequences
# iii) host labelled as ‘Human’
# iv) sequence length above 28,000 base pairs
# v) proportion of ambiguous bases below 2%.

metadata.filtered <- metadata %>%
  distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         !N.Content > 0.02 | is.na(N.Content),
         str_length(Collection.date) == 10,
         Sequence.length > 28000,
         Passage.details.history == "Original",
         Pango.lineage != "",
         Pango.lineage != "None",
         Pango.lineage != "Unassigned",
         !str_detect(Additional.location.information,"[Qq]uarantine")
  )

metadata.filtered <- metadata.filtered %>%
  mutate(Collection.date = as.Date(Collection.date),
         region = str_split(Location," / ",simplify = T)[,1],
         country = str_split(Location," / ",simplify = T)[,2],
         state = str_split(Location," / ",simplify = T)[,3])

metadata.filtered <- metadata.filtered[!duplicated(metadata.filtered$Virus.name),]

## nextlade
nextclade <- fread(nextclade.name, header=T, sep="\t", quote="", check.names=T)
nextclade <- nextclade %>% mutate(seqName=str_split(seqName, "\\|", simplify=T)[,1])

nextclade.filtered <- nextclade %>% filter(seqName %in% metadata.filtered$Virus.name)
nextclade.filtered <- nextclade.filtered[!duplicated(nextclade.filtered$seqName),]

metadata.filtered <- metadata.filtered %>% filter(Virus.name %in% nextclade.filtered$seqName)

##Add Nextclade's clade information to the filtered GISAID metadata
metadata.filtered <- metadata.filtered %>% mutate(Nextclade_clade = nextclade.filtered$clade_nextstrain)

##merge nextclade and GISAID
metadata.filtered$Nextclade_clade <- nextclade.filtered[match(metadata.filtered$Virus.name, nextclade.filtered$seqName), "clade_nextstrain"]

count.pango.df <- metadata.filtered %>% group_by(Nextclade_clade) %>% summarize(count.pango = n())
count.pango.df

###sampled 20 per clade based on Nextclade_clade
num.samples <- 20
metadata.filtered.sampled <- metadata.filtered %>%
  filter(
    !is.na(Nextclade_clade), 
    Nextclade_clade != "", 
    Nextclade_clade != "recombinant", #removing any recombinant group
    Nextclade_clade != "21M" #removing 21M (Omicron, B.1.1.529)
  ) %>%
  group_by(Nextclade_clade) %>%
  slice_sample(n = num.samples)


write_tsv(x=metadata.filtered.sampled, file=tsv.merged.name)
length(unique(metadata.filtered.sampled$Accession.ID))

count.pango.sampled.df <- metadata.filtered.sampled %>% group_by(Nextclade_clade) %>% summarize(count.pango = n()) %>% ungroup()
count.pango.sampled.df

accessionID.list <- as.data.frame(metadata.filtered.sampled)
accessionID.list <- select(accessionID.list, Accession.ID)
write_tsv(x=accessionID.list, file=tsv.ID.name,col_names = F)

#done

