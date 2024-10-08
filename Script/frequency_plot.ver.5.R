#!/usr/bin/env R

library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)


##########args##########
#Change when using new input
download_date <- "2024-02-19"
out_prefix <- "2024_02_19"
download.date <- as.Date(download_date) 
#input
metadata.name <- "metadata.tsv"
nextclade.name <- "nextclade.cut.tsv"

#output
pdf.observed.name <- paste(out_prefix,".method1.observed.pdf",sep="")
Id_list.name <- paste(out_prefix,".id_list.csv",sep="")
table.name <- paste(out_prefix,".variant_proportion",sep="") 

dir <- ""
setwd(dir)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########parameters##########
#general
core.num <- 4
#variant.ref <- "BQ.1.1" #"BA.5" #"BA.1.1" # This is now set below after the metadata has been filtered.

#period to be analyzed
date.end <- download.date
date.start <- "2022-12-01"

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

metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)

temp <- nextclade.table %>% left_join(metadata,by=c('seqName'='Virus.name'), relationship = "many-to-many") # Merge the nextclade and Gisaid data
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
clade.name <- read_csv("clade_name_list.csv")

variant.top10.v <- c("other", "XBB.1.5", "XBB.1.16", "EG.5.1", "BA.2.86", "JN.1")
clade.name.v <- c("19A","19B","20A","20B","20C","20D","20D","20E","20F","20G","20H","20I","20J","21A","21B","21C","21D","21E","21F","21G","21H",
                  "21I","21J","21K","21L","21M","22A","22B","22C","22D","22E","22F","23A","23B","23C","23D","23E","23F","23G","23H","23I")

metadata.filtered.analyzed.2 <- metadata.filtered.analyzed %>% filter(clade %in% clade.name.v)
metadata.plot <- left_join(metadata.filtered.analyzed.2, clade.name, by = "clade")

metadata.plot <- metadata.plot %>% mutate(cladename.jn1 = ifelse(cladename == "BA.2.86" & str_detect(AA.Substitutions, "Spike_L455S"), "JN.1", cladename))

col.v <- c(brewer.pal(length(variant.top10.v), "Paired"),"gray70")
metadata.plot$cladename.jn1 <- factor(metadata.plot$cladename.jn1,  levels=c("other","XBB.1.5","XBB.1.16","EG.5.1","BA.2.86","JN.1"))

g <- ggplot(metadata.plot,aes(x=Collection.date,fill=cladename.jn1))
g <- g + geom_bar(stat = 'count', position = "fill")
g <- g + scale_x_date(date_labels = "%y-%m", date_breaks = "1 months", date_minor_breaks = "1 month")
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

##### extract ID
Id_list <- c(nextclade.table$seqName)
Id_list <- metadata.filtered.analyzed.2 %>% filter(Virus.name %in% Id_list)
Id_list <- c(Id_list$Accession.ID)

write.csv(Id_list, file=Id_list.name, row.names=FALSE)
dev.off()

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

variant_counts <- metadata.plot %>%
  group_by(Collection.date, cladename.jn1) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the total count per date
total_counts <- variant_counts %>%
  group_by(Collection.date) %>%
  summarise(total = sum(count), .groups = 'drop')

# Merge the counts with the totals and calculate proportions
variant_proportions <- variant_counts %>%
  left_join(total_counts, by = "Collection.date") %>%
  mutate(proportion = count / total) %>%
  select(Collection.date, cladename.jn1, proportion)

# Pivot the data to get the desired format
variant_proportion_wide <- variant_proportions %>%
  pivot_wider(names_from = cladename.jn1, values_from = proportion, values_fill = 0)

# View the final data frame
print(variant_proportion_wide)

write.csv(variant_proportion_wide, file=table.name)
dev.off()



