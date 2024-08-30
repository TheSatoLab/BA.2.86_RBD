#!/usr/bin/env R

library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=T)

##########args##########
#Change when using new input
download_date <- "2024-07-01"
download.date <- as.Date(download_date)
out_prefix <- "2024_07_01"

metadata.name <- paste("/Volumes/annin-tofu/variant_monitoring/output/",out_prefix,"/metadata_tsv_",out_prefix,"/metadata.tsv",sep = "")
mut.info.name <- paste("/Volumes/annin-tofu/variant_monitoring/output/",out_prefix,"/metadata_tsv_",out_prefix,"/metadata.mut_long.tsv",sep = "")
nextclade.name <- paste("/Volumes/annin-tofu/variant_monitoring/output/",out_prefix,"/full.nextclade.tsv",sep = "")
stan_f.name <- '/Users/petadimensionlab/Dropbox/github_ijampei/next_variant_detection_220413_share/method1/multinomial_independent.ver2.stan'
#output
out.prefix <- "output/frequency_plot"
pdf.observed.name <- paste(out.prefix,".pdf",sep="")
tsv.ID.name <- paste(out_prefix,".freq.id_list.csv",sep="")
dir <- paste("/Volumes/annin-tofu/2405_Conference/slide_figures")
setwd(dir)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########parameters##########
#general
core.num <- 4

#period to be analyzed
date.end <- download.date
date.start <- "2021-11-01"

#min numbers
limit.count.analyzed <- 50 #was 50 then 20, now 50 again

#Transmissibility
bin.size <- 1
generation_time <- 2.1

#Number of top variants to get mutation
mut_num <- 25


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
date.start
date.end


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#model


# multi_nomial_model <- cmdstan_model(stan_f.name)

##########data preprocessing & QC##########
###fasta file should be filtered using python script in advance
nextclade.table <- fread(nextclade.name,header=T,sep="\t",quote="",check.names=T)
nextclade.table <- nextclade.table %>% mutate(Virus.name = sub("\\|.*", "", seqName))

metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)
#mut.info <- fread(mut.info.name,header=T,sep="\t",quote="",check.names=T)

temp <- nextclade.table %>% left_join(metadata,by=c('seqName'='Virus.name')) # Merge the nextclade and Gisaid data
temp <- temp %>% rename("Virus.name" = "seqName") %>% select(-Pango.lineage) %>% rename("Pango.lineage"="Nextclade_pango")
metadata <- temp

metadata.filtered <- metadata %>%
  distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         str_length(Collection.date) == 10,
         Pango.lineage != "",
         Pango.lineage != "None",
         !str_detect(Additional.location.information,"[Qq]uarantine")
  )

metadata.filtered <- metadata.filtered %>%
  mutate(Collection.date = as.Date(Collection.date),
         region = str_split(Location," / ",simplify = T)[,1],
         country = str_split(Location," / ",simplify = T)[,2])

metadata.filtered.analyzed <- metadata.filtered %>% filter(Collection.date >= date.start, Collection.date <= date.end)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### using the nextclade "clade" column, we categorise all variants into 5 groups
#plot observed
metadata.filtered.analyzed.2 <- metadata.filtered.analyzed %>% select(clade_display, Collection.date)
metadata.filtered.analyzed.2 <- metadata.filtered.analyzed.2 %>% filter(!clade_display=="") %>% mutate(clade = gsub(".*\\(([^()]*)\\).*", "\\1", clade_display))

##kokokara Yurinchi
clade.interest <- c("Delta","BA.1","BA.2","BA.4","BA.5","BA.2.75","BQ.1","XBB","XBB.1.5","EG.5.1","BA.2.86","HK.3","JN.1","KP.2","KP.3","LB.1","KP.2.3","KP.3.1.1")
colorpallet <- fread("G2P-Japan_colorpallet.txt" ,header=T,sep="\t",quote="",check.names=T)
colorpallet <- as.data.frame(clade.interest) %>% left_join(colorpallet, by = c('clade.interest'='Pango.lineage'))

#metadata.filtered.analyzed.2 <- fread("output/metadata.filtered.analyzed.txt" ,header=T,sep="\t",quote="",check.names=T)
metadata.filtered.analyzed.2 <- metadata.filtered.analyzed.2 %>% 
                                mutate(col= ifelse(metadata.filtered.analyzed.2$clade %in% colorpallet$clade.interest, metadata.filtered.analyzed.2$clade, "other"))

metadata.filtered.analyzed.2$col <- factor(metadata.filtered.analyzed.2$col, levels = c("other", colorpallet$clade.interest))

colorpallet <- data.frame(clade.interest="other", col_code="#D2D1D3") %>% rbind(colorpallet)
colorpallet <- colorpallet %>%
               mutate(col_code = if_else(row_number() == 9, "#736357", col_code))

barplot(rep(1,13), col = colorpallet$col_code) #check color

col.v <- colorpallet$col_code

g <- ggplot(metadata.filtered.analyzed.2,aes(x=Collection.date,fill=col))
g <- g + geom_bar(stat = 'count', position = "fill", aes(fill=col))
g <- g + scale_x_date(date_labels = "%y-%m", date_breaks = "3 months", date_minor_breaks = "2 month")
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + scale_fill_manual(values=col.v)
g

#plot observed
pdf(pdf.observed.name,width=10,height=6)
plot(g)
dev.off()

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

write.table(metadata.filtered.analyzed.2,"output/metadata.filtered.analyzed.txt",row.names=F,sep="\t",quote=F)
