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
  tidyverse, vegan, viridis, webshot2, worrms, taxize
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

# Pivot longer #
eDNA.reads %>% pivot_longer(
  cols = 4:ncol(.),
  names_to = "sample",
  values_to = "RA"
) %>% filter(RA > 0) -> eDNA.reads # Remove zeroes to reduce dataframe size #

#### Append metadata to eDNA reads ####
read.csv(
  "raw_data/metadata.csv",
  header = TRUE, check.names = FALSE
) %>% right_join(., eDNA.reads, by = "sample") -> eDNA.reads

#### Calculate triplicates ####
# Calculate adjusted read counts based on triplicates
eDNA.reads %>% select(c(
  ASV_ID, sample, RA, replicate, triplicate
)) -> df # Working dataframe
df %>% select(c(
  ASV_ID, replicate, RA, triplicate
)) %>% group_by(ASV_ID, replicate) %>%
  filter(triplicate == "yes") %>%
  select(c(ASV_ID, replicate, RA)) %>%
  summarise(RA = sum(RA)/3) %>% # Average to triplicate
  mutate(RA = floor(RA)) -> rep.df # Round down
rep.samples = c("MTeDNA_S15", "MTeDNA_S15_2", "MTeDNA_S25",
                "MTeDNA_S25_2", "MTeDNA_S31", "MTeDNA_S31_2")
rep.df %>% filter(
  replicate %in% rep.samples
) %>%
  mutate(
    dup = case_when(
      replicate == "MTeDNA_S15" ~ "MTeDNA_S15",
      replicate == "MTeDNA_S15_2" ~ "MTeDNA_S15",
      replicate == "MTeDNA_S25" ~ "MTeDNA_S25",
      replicate == "MTeDNA_S25_2" ~ "MTeDNA_S25",
      replicate == "MTeDNA_S31" ~ "MTeDNA_S31",
      replicate == "MTeDNA_S31_2" ~ "MTeDNA_S31"      
    )
  ) -> rep.df2
grp = c("ASV_ID", "dup", "RA")
rep.df2 <- rep.df2 %>%
  group_by(across(all_of(grp))) %>%
  summarise(RA = (sum(RA))/2) %>% # Treat duplicates #
  rename(replicate = dup) %>%
  mutate(RA = floor(RA))
rep.df %>% filter(!(replicate %in% rep.samples)) %>%
  rbind(., rep.df2) -> rep.df

# Produce working dataframe #
read.csv(
  "raw_data/metadata.csv",
  header = TRUE, check.names = FALSE) %>%
  select(-c(sample, triplicate)) %>%
  unique() %>%
  right_join(., rep.df, by = "replicate", relationship = "many-to-many") %>%
  rename(sample = replicate) -> df
remove(rep.df)
remove(rep.df2)
remove(eDNA.reads)

#### Read in BLAST results ####
read.csv(
  "raw_data/blast.csv",
  header = TRUE, check.names = FALSE,
  row.names = 1
) -> blast

#### Use the 'worrms' and 'taxize' packages to query the BLAST results ####
blast$species %>% unique() -> worms.match # Pull unique taxa #
length(worms.match) # Number of unique taxa
wm_records_names(name = worms.match[c(1:121)]) -> query.A
wm_records_names(name = worms.match[c(122:242)]) -> query.B
wm_records_names(name = worms.match[c(243:length(worms.match))]) -> query.C
# The WoRMS query seems to break at around ~150 species, so I split this into thirds #
# Depending on your filter parameters, you may need to adjust the length of these queries #

# Use lapply and bind_rows() to turn these lists into dataframes to add together #
bind_rows(lapply(query.A, as.data.frame)) -> query.A 
bind_rows(lapply(query.B, as.data.frame)) -> query.B
bind_rows(lapply(query.C, as.data.frame)) -> query.C 
rbind(query.A, query.B, query.C) -> worms.match # Replace other object with the matches #
remove(query.A)
remove(query.B)
remove(query.C)

# Clean up the match dataframe #
worms.match %>% select(c(kingdom:genus, scientificname)) %>% 
  rename(species = scientificname) -> worms.match

# Refine BLAST results #
blast %>% filter(species %in% worms.match$species) %>%
  left_join(., worms.match, by = "species", relationship = "many-to-many") %>%
  filter(!is.na(ASV_ID)) -> blast
remove(worms.match)

#### Append taxonomy to reads ####
df %>% right_join(., blast, by = "ASV_ID", relationship = "many-to-many") -> df

#### Filter out bad taxonomy reads ####
df %>% filter(bit_score > 250 & percent_identity > 95.0) %>%
       filter(kingdom != "Fungi") -> df

#### Create community dataframe, normalize ####
df %>% select(c(sample, species, RA)) %>%
  pivot_wider(
    names_from = species,
    values_from = RA,
    values_fn = sum
  ) %>% mutate_at(c(2:ncol(.)), ~replace_na(., 0)) %>%
  filter(!is.na(sample)) %>% column_to_rownames(var = "sample") %>%
  decostand(., method = "log") %>% decostand(., method = "total") -> comm
# Used "decostand" from the vegan package to apply community transformations

#### Make metadata table for NMDS / BCD calcs ####
read.csv(
  "raw_data/metadata.csv",
  header = TRUE, check.names = FALSE) %>%
  select(c(
  replicate, type:ncol(.))) %>% unique() %>%
  rename("sample" = "replicate") %>%
  filter(sample %in% rownames(comm)) -> nmds.meta

#### PERMANOVA ####
set.seed(12031996)
prop.adonis <- adonis2(comm ~ depth_cat + month,
                       data = nmds.meta,
                       permutations = 999,
                       method = "bray",
                       by = "terms")
prop.adonis # Model output

#### Ordination ####
set.seed(12031996)
comm.mds <- metaMDS(comm, dist = "bray", k = 5, trymax = 100)
comm.mds$stress # Low stress is good!

#### NMDS plot ####
data_scores <- as.data.frame(scores(comm.mds, "sites"))
sp_scores <- as.data.frame(scores(comm.mds, "species"))
sp_scores$species <- rownames(sp_scores)
scrs <- scores(comm.mds, display = "sites") %>%
  as.data.frame()
scrs$sample <- rownames(scrs)
scrs <- scrs %>%
  as.data.frame() %>%
  left_join(.,
            nmds.meta,
            by = "sample")

cent <- aggregate(cbind(NMDS1, NMDS2) ~ depth_cat, data = scrs, FUN = mean)
segs <- merge(scrs, setNames(cent, c("depth_cat","oNMDS1","oNMDS2")),
              by = "depth_cat", sort = FALSE)
ggplot(scrs, aes(x = NMDS1, y = NMDS2, color = depth_cat)) +
  geom_segment(data = segs, mapping = aes(xend = oNMDS1, yend = oNMDS2), alpha = 0.3) +
  geom_point(data = cent, size = 4) + 
  geom_point(size = 2, alpha = 0.75) + coord_fixed() +
  scale_color_viridis(discrete = TRUE, option = "A", direction = -1) +
  theme(
    panel.background = element_rect(fill = "#82B1DB"),
    axis.line = element_line(colour = "black")) +
  theme(legend.position = "bottom") +
  labs(title = "Challenger Deep eDNA Ordination",
       subtitle = "Bray-Curtis Dissimilarity on Log-Transformed Read Abundance",
       xlab = "NMDS 1", ylab = "NMDS 2",
       color = "Depth Categories")

#### Barplot figure, facet wrap ####
df %>% 
  group_by(depth_cat, phylum) %>%
  filter(!is.na(depth_cat)) %>%
  summarise(RA = sum((log1p(RA)))) %>%
  arrange(desc(RA)) %>%
  ggplot(aes(x = as.factor(phylum), y = RA, fill = phylum)) +
  geom_bar(stat = "identity", colour = "black", width = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "turbo", direction = -1) +
  theme(
    panel.background = element_rect(fill = "#D8EDE6"),
    axis.line = element_line(colour = "black")
  ) +
  facet_wrap(vars(depth_cat)) +
  coord_flip() + scale_x_discrete(limits = rev) +
  labs(
    title = "Challenger Deep eDNA Diversity by Depth",
    subtitle = "Phyla by depth category",
    y = "Log-Transformed Read Abundance", x = "", fill = "eDNA-Observed Phyla") +
  theme(legend.position = "none")
ggsave(
  "figures/eDNA_by_depth.png",
  units = "px",
  dpi = 300,
  width = 3000,
  height = 2000,
  bg = "white"
)

#### End ####
gommage()
dev.off()
# Get out of this canvas! #
