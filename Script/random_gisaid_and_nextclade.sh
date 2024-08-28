#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 18 ##12
#$ -l s_vmem=5G
#$ -q '!mjobs_rerun.q'
#$ #-o /dev/null
#$ #-e /dev/null

####script
freq_script="frequency_plot.ver.5.R"
gis_and_nc_script="BA.2.86_gisaid_and_nextclade_merged.R"

####input
nextclade="nextclade"
download_date="2024-02-19" #data download date (actually upload date). Use hyphens.
out_prefix="2024_02_19" 
working_dir="" #path for the output prefix. Date should match Gisaid file name.
gisaid_metadata="metadata.tsv" #path to the GISAID metadata file
nextclade="nextclade.tsv"
nextclade_cut="nextclade.cut.tsv"

####Command
cd ${working_dir}
mkdir -p ${out_prefix}_random

####nextclade
## tar -xf sequences_fasta_${out_prefix}.tar.xz

## ${nextclade} run -d sars-cov-2 -j 8 --output-tsv=nextclade.tsv sequences.fasta

####frequency plot

module load R/4.3.0

R --vanilla --slave --args \
  "${download_date}" \
  ${out_prefix} \
  ${gisaid_metadata} \
  ${nextclade_cut} \
  ${out_prefix}_random/${out_prefix} \
  < ${freq_script}

####build phylogenetic tree
####get id
cut -f 1-8 ${nextclade} > nextclade.cut.tsv

R --vanilla --slave --args \
  "${download_date}" \
  ${out_prefix} \
  ${gisaid_metadata} \
  ${nextclade_cut} \
  ${out_prefix}_random/${out_prefix} \
  < ${gis_and_nc_script}

