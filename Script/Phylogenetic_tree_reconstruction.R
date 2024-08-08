#!/usr/bin/env R

library(data.table)
library(tidyverse)
library(ape)
library(Biostrings)
library(EnvStats) # Rosner's test 
library(ggplot2)
library(ggtree)

### Change working directory
dir <- "~/Desktop/JN.1_RBD_project/"
setwd(dir)

#################### Retrieve metadata of SARS-CoV-2 subvariants ####################

metadata <- fread("metadata_tsv/metadata.tsv", header=T, sep="\t", quote="", check.names=T)

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
                              # Pango.lineage != "", Pango.lineage != "None", Pango.lineage != "Unassigned",
                              !str_detect(Additional.location.information,"[Qq]uarantine")
                              )

metadata.filtered <- metadata.filtered %>%
                       mutate(Collection.date = as.Date(Collection.date),
                              region = str_split(Location," / ",simplify = T)[,1],
                              country = str_split(Location," / ",simplify = T)[,2],
                              state = str_split(Location," / ",simplify = T)[,3])

metadata.filtered <- metadata.filtered[!duplicated(metadata.filtered$Virus.name),]

### Get 20 sequences of JN.1 (confirmed later with NextClade that they belong to NextClade clade 23I)
metadata.filtered.JN.1 <- metadata.filtered %>% filter(Pango.lineage == "JN.1") %>% slice_sample(n = 20) 

### Get accession number of SARS-CoV-2 used in the previous BA.2.86 paper
BA.2.86_paper_acc_list <- (read.delim("230919_accession.ID.v1.txt", header=FALSE, sep="\n"))$V1
selected_acc_list <- c(BA.2.86_paper_acc_list, metadata.filtered.JN.1$Accession.ID)

### Retrieve GISAID metadata of SARS-CoV-2 genomic sequences listed in previous study and those of newly added JN.1 and filter data with the same criterion
metadata.filtered <- metadata %>% filter(Accession.ID %in% selected_acc_list)
metadata.filtered <- metadata.filtered %>%
  distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         !N.Content > 0.02 | is.na(N.Content),
         str_length(Collection.date) == 10,
         Sequence.length > 28000,
         Passage.details.history == "Original",
         # Pango.lineage != "", Pango.lineage != "None", Pango.lineage != "Unassigned",
         !str_detect(Additional.location.information,"[Qq]uarantine")
  )

# write.table(metadata.filtered$Accession.ID, "240301_EPI_set_unfiltered_list.tsv", col.names=F, row.names=F, sep="\n", quote=F)
# write.table(metadata.filtered, "metadata_tsv/240301_filtered_metadata.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#################### Remove taxa that may cause long branch attraction using Rosner's test ####################

metadata.filtered <- fread("metadata_tsv/240301_filtered_metadata.tsv", header=T, sep="\t", quote="", check.names=T)
nextclade.filtered <- fread("metadata_tsv/240301_filtered_nextclade.tsv", header=T, sep="\t", quote="", check.names=T)
nextclade.filtered <- nextclade.filtered %>% mutate(seqName=str_split(seqName, "\\|", simplify=T)[,1])

metadata.filtered$Virus.name <- str_replace_all(metadata.filtered$Virus.name, " ", "_")
metadata.filtered$Nextclade_clade <- nextclade.filtered[match(metadata.filtered$Virus.name, nextclade.filtered$seqName), "clade_nextstrain"]
metadata.filtered$Nextclade_pango <- nextclade.filtered[match(metadata.filtered$Virus.name, nextclade.filtered$seqName), "Nextclade_pango"]
metadata.filtered$Virus.name <- str_replace_all(metadata.filtered$Virus.name, "'", "_")

### Remove taxa that may cause long branch attraction using Rosner's test ##########
jn.1_tree <- read.tree("gisaid_hcov-19_2024_03_01_06.fasta.edited.aln.trimal.trimmed.treefile")

fasta_file <- readDNAStringSet("gisaid_hcov-19_2024_03_01_06.fasta.edited.aln.trimal.trimmed", "fasta")
fasta_seq_name <- names(fasta_file)
fasta_seq <- paste(fasta_file)
data_fasta <- data.frame(label=fasta_seq_name, seq=fasta_seq)
data_fasta$label <- gsub(" ", "_", data_fasta$label)
data_fasta_filtered <- data_fasta

loop_index <- TRUE
while (loop_index == TRUE) {
  ### Run Rosner's generalized extreme Studentized deviate test to remove branch length outliers
  jn.1_tree.info.df <- ggtree(jn.1_tree)$data
  jn.1_tree.info.df <- jn.1_tree.info.df[,c('label','branch.length')]
  jn.1_tree.info.df$branch.length.log <- log(jn.1_tree.info.df$branch.length+1,10)
  
  rosner_test <- rosnerTest(jn.1_tree.info.df$branch.length.log, k=10, alpha=1e-4)
  rosner_test
  
  if (TRUE %in% rosner_test$all.stats$Outlier) {
    ### Filter outlier taxa out and write output FASTA file for phylogenetic tree reconstruction
    jn.1_tree <- drop.tip(jn.1_tree, jn.1_tree.info.df[rosner_test$all.stats[which(rosner_test$all.stats$Outlier==TRUE),]$Obs.Num,]$label)
    data_fasta_filtered <- data_fasta_filtered[which(!data_fasta_filtered$label %in% jn.1_tree.info.df[rosner_test$all.stats[which(rosner_test$all.stats$Outlier==TRUE),]$Obs.Num,]$label),]
    if (all(jn.1_tree.info.df[rosner_test$all.stats[which(rosner_test$all.stats$Outlier==TRUE),]$Obs.Num,]$label %in% metadata.filtered$Virus.name)==FALSE) {loop_index <- FALSE} # Solve loop problem
  } else {
    loop_index <- FALSE
  }
}

data_fasta_filtered$label <- paste('>', data_fasta_filtered$label, sep="")
# write.table(data_fasta_filtered, "gisaid_hcov-19_2024_03_01_06.fasta.edited.aln.trimal.trimmed.outlier_filtered", col.names=F, row.names=F, sep="\n", quote=F)

#################### Visualize a phylogenetic tree ####################

### Read the tree after removing outliers
jn.1_tree <- read.tree("gisaid_hcov-19_2024_03_01_06.fasta.edited.aln.trimal.trimmed.outlier_filtered.treefile")
# write.table(str_split(jn.1_tree$tip.label, "\\|", simplify=T)[,2], "240301_EPI_set_filtered_list.tsv", col.names=F, row.names=F, sep="\n", quote=F)

jn.1_tree$tip.label <- str_split(jn.1_tree$tip.label, "\\|", simplify=T)[,1]

### Plot the final tree
metadata.filtered.plot <- metadata.filtered %>% mutate(group = case_when(
  Nextclade_clade %in% c('19A', '19B', '20A', '20B', '20C', '20D', '20E', '20F', '20G', '21B', '21C', '21D', '21E', '21G', '21H', '21F') ~ 'outgroup',
  Nextclade_clade %in% c('20I') ~ 'Alpha',
  Nextclade_clade %in% c('20H') ~ 'Beta',
  Nextclade_clade %in% c('20J') ~ 'Gamma',
  Nextclade_clade %in% c('21A', '21I', '21J') ~ 'Delta',
  Nextclade_clade %in% c('21K') ~ 'BA.1',
  Nextclade_clade %in% c('22A') ~ 'BA.4',
  Nextclade_clade %in% c('21L','22C') ~ 'BA.2',
  Nextclade_clade %in% c('22B') ~ 'BA.5',
  Nextclade_clade %in% c('22D', '23C') ~ 'BA.2.75',
  Nextclade_clade %in% c('22E') ~ 'BQ.1',
  Nextclade_clade %in% c('22F', '23E') ~ 'XBB.1',
  Nextclade_clade %in% c('23A') ~ 'XBB.1.5',
  Nextclade_clade %in% c('23B') ~ 'XBB.1.16',
  Nextclade_clade %in% c('23D', '23F', '23H') ~ 'EG.5.1',
  str_detect(Nextclade_pango, 'JN.1') ~ 'JN.1',
  Nextclade_clade %in% c('23H') ~ 'HK.3',
  Nextclade_clade %in% c('23I') ~ 'BA.2.86',
  TRUE ~ NA_character_ 
))

color.matching <- c(
  "outgroup" = "grey",
  "Alpha" = "#1C67A8",
  "Beta" = "#39A04A",
  "Gamma" = "#D7C046",
  "Delta" = "#E77A25",
  "BA.1" = "#339680",
  "BA.2" = "#02686B",
  "BA.5" = "#803345",
  "BA.2.75" = "#898731",
  "BQ.1" = "#67468C",
  "XBB.1" = "#736357",
  "XBB.1.5" = "#981D2A",
  "XBB.1.16" = "#D95A24",
  "EG.5.1" = "#DE0303",
  "BA.2.86" = "#1D9FF0",
  "BA.4" = "#8C2788",
  "JN.1" = "yellow"
)

jn.1_MRCA <- MRCA(jn.1_tree, (metadata.filtered.plot %>% filter(group=="JN.1"))$Virus.name)
ba.2.86_MRCA <- MRCA(jn.1_tree, (metadata.filtered.plot %>% filter(group %in% c("BA.2.86", "JN.1")))$Virus.name)

# pdf("JN.1_phylogenetic_tree.pdf", width=6, height=12)
jn.1_tree_plot <- ggtree(jn.1_tree, layout="equal_angle", size=0.25) + labs(group="Clade") + geom_text2(aes(subset = node %in% c(ba.2.86_MRCA, jn.1_MRCA), label=str_split(label, "/", simplify = T)[,2]), size=4, hjust=1.3, alpha=1)
jn.1_tree_plot <- jn.1_tree_plot %<+% metadata.filtered.plot + geom_tippoint(aes(color=group), size=3, alpha=1) + scale_color_manual(values = color.matching) 
jn.1_tree_plot
# dev.off()
