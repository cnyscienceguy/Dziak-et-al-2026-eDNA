# Novel observations of Earth’s ultimate hadal zone: Passive Acoustic records and environmental DNA samples from Challenger Deep, Mariana Trench
# Dziak et al., in prep 
# eDNA research / analysis component
# Script generated 2026-01-28
# Written by Charles Nye, PhD Candidate (Oregon State University)

#### Start  ####
# Import packages #
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  ape, data.table, ggpubr, ggrepel, ggtext, kableExtra, knitr, mgcv, pacman,
  tidyverse, vegan, viridis, webshot2
)
# This function will clear out the working environment at the end of the script #
gommage = function(){
  gctorture(TRUE)
  rm(list = ls(all.names = TRUE, envir = sys.frame(-1)),
     envir = sys.frame(-1))
  gctorture(FALSE)
}
# For those who come after #

#### Import eDNA reads ####
read.table(
  "raw_data/conv-feature-table.txt",
  sep = "\t",
  header = TRUE
) -> eDNA.reads

#### Import DNA sequences from .fasta file, append to reads ####
read.FASTA(
  "raw_data/dna-sequences.fasta", type = "AA"
) %>% as.character() %>% sapply(paste, collapse = "") %>%
  as.data.frame() %>% rename("sequence" = ".") %>%
  rownames_to_column(var = "ASV_ID") %>%
  left_join(eDNA.reads, ., by = "ASV_ID") %>%
  select(c(ASV_ID, sequence, everything())) -> eDNA.reads
# Let's have a look #
head(eDNA.reads)
colnames(eDNA.reads) # List of samples
prefilter.richness = length(unique(eDNA.reads$ASV_ID))
prefilter.totRA = sum(rowSums(eDNA.reads[,-c(1:2)]))
prefilter.meanRA = mean(rowSums(eDNA.reads[,-c(1:2)])) 

#### Filter DNA sequences for on-target amplicons ####
eDNA.reads %>% mutate(seq_length = nchar(sequence)) %>%
  filter(
    seq_length > 300 & seq_length < 330
  ) %>% select(c(1:2, seq_length, everything())) -> eDNA.reads

#### Descriptive metrics ####
prefilter.richness # ASV richness prior to filtering
length(unique(eDNA.reads$ASV_ID)) # ASV richness post-filtering
prefilter.richness - length(unique(eDNA.reads$ASV_ID)) # This is the difference

prefilter.totRA # Total read abundance prior to filtering
sum(rowSums(eDNA.reads[,-c(1:3)])) # Post-filtering
prefilter.totRA - sum(rowSums(eDNA.reads[,-c(1:3)])) # This is the difference

prefilter.meanRA # Average per-sample read abundance prior to filtering
mean(rowSums(eDNA.reads[,-c(1:3)])) # Post-filtering


#### Read in metadata ####
read.csv(
  "raw_data/metadata.csv",
  header = TRUE, check.names = FALSE
) -> metadata
# Append to reads table to create working dataframe #

#### Read in BLAST results ####
read.csv(
  "raw_data/blast.csv",
  header = TRUE, check.names = FALSE,
  row.names = 1
) -> taxonomy

#### End ####
gommage()
# Get out of this canvas! #
