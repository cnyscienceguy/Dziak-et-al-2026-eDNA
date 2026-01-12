# Marianas Trench eDNA research component
# Last edited October 2025
# Dziak et al., in prep

rm(list = ls()) # Refresh everything
if (!require("pacman")) install.packages("pacman")
pacman::p_load(pacman, tidyverse, data.table, ggtext, ggpubr, viridis, vegan, knitr, devtools) #Some packages might need some finesse
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")
#Ensure all packages load

# Reads table
reads.raw = read_table(
  "inputs/conv-feature-table.txt"
) %>%
  as.data.frame()
# Blast output
blast.raw = read.csv(
  "inputs/blast.csv",
  header = TRUE
)
# Metadata
meta.raw = read.csv(
  "inputs/metadata.csv",
  header = TRUE
)
# Convert fasta into sequence list
FastaToTabular(
  "inputs/dna-sequences.fasta"
)
seqs.raw = read.csv(
  "dna_table.csv",
  header = TRUE
) %>%
  rename(ASV_ID = name)

# Clean-up input sequence input
seqs.raw$ASV_ID = sub(".", "", seqs.raw$ASV_ID)
seqs = seqs.raw[,-1]

# Create reads dataframe, use script kung-fu to split and resplice
reads <- reads.raw
rownames(reads) <- reads$ASV_ID
reads$ASV_ID <- rownames(reads)
reads <- reads %>%
  select(ASV_ID,everything())

# Transform into long format
reads <- reads %>%
  pivot_longer(
    cols = 2:ncol(reads),
    names_to = "sample",
    values_to = "read_abundance"
  ) %>%
  filter(read_abundance > 0)

# Calculate pre-filtering ASV richness
ASV.richness.pre <- length(unique(reads$ASV_ID))

# Filter ASVs by certain criteria
reads2 <- reads %>%
  select(c(ASV_ID, read_abundance))
reads2 <- reads2 %>%
  group_by(ASV_ID) %>%
  summarise(read_abundance = sum(read_abundance)) %>%
  filter(read_abundance > 1) #Remove all singletons

# List of ASVs that pass filtering
passing_ASVs <- as.list(reads2$ASV_ID)

# Filter reads dataframe to only include passing ASVs
reads <- reads %>%
  .[reads$ASV_ID %in% c(passing_ASVs), ]

# Calculate passing ASV richness
ASV.richness.post <- length(unique(reads$ASV_ID))
# Calculate percentage of passing ASVs
ASV.percent.pass <- ASV.richness.post / ASV.richness.pre

# Calculate adjusted read counts based on triplicates
# Append metadata
reads <- reads %>%
  left_join(.,
            meta.raw,
            by = "sample")
# Remove superfluous columns, for now
reads <- reads %>%
  select(c(ASV_ID, sample, read_abundance, replicate, triplicate))

# Calculate adjustments
rep.df <- reads %>%
  select(ASV_ID, replicate, read_abundance, triplicate) %>%
  group_by(ASV_ID, replicate) %>%
  filter(triplicate == "yes") %>%
  select(ASV_ID, replicate, read_abundance) %>%
  summarise(read_abundance = sum(read_abundance)/3)
# Round down
rep.df$read_abundance <- floor(rep.df$read_abundance)

rep.samples = c(
  "MTeDNA_S15",
  "MTeDNA_S15_2",
  "MTeDNA_S25",
  "MTeDNA_S25_2",
  "MTeDNA_S31",
  "MTeDNA_S31_2"
)
rep.df2 = rep.df %>% filter(
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
  )
grp = c("ASV_ID", "dup", "read_abundance")
rep.df2 <- rep.df2 %>%
  group_by(across(all_of(grp))) %>%
  summarise(
    read_abundance = (sum(read_abundance))/2
  ) %>%
  rename(
    replicate = dup
  )
rep.df2$read_abundance = floor(rep.df2$read_abundance)
rep.df <- rep.df %>%
  filter(!(replicate %in% rep.samples))
rep.df <- rbind(rep.df, rep.df2)

# Append sequence data to reads dataframe to create primary dataframe
df <- rep.df %>%
  left_join(.,
            seqs,
            by = "ASV_ID")

# Append BLAST results
blast <- blast.raw[,-1]

blast <- blast %>%
  mutate(
    tax.assignment =
      case_when(
        percent_identity < 97.9 ~ phylum,
        TRUE ~ species
      )
  ) %>%
  filter(kingdom != "Fungi" & class != "Arachnida" & class!= "Insecta")

df <- df %>%
  left_join(
    .,
    blast,
    by = "ASV_ID"
  ) %>%
  filter(bit_score > 250) %>%
  filter_at(vars(c(kingdom, phylum, class, order, family)), all_vars(!is.na(.))) %>%
  filter(kingdom != "Fungi" & class != "Arachnida" & class!= "Insecta") #Remove fungi; there are marine fungi, but they are not well-studied or necessarily of interest

meta <- meta.raw %>%
  select(-sample) %>%
  unique()

df <- df %>%
  left_join(.,
            meta,
            by = "replicate") %>% 
  rename(sample = replicate) %>%
  unique() %>%
  as.data.frame()
meta <- meta %>%
  rename(sample = replicate)

#Calculate ASV richness post-taxonomic filtering
df.ASV.richness <- length(unique(df$ASV_ID))
#Calculate species richness post-taxonomic filtering
df.sp.richness <- length(unique(df$species))
#Calculate taxa assignment richness post-taxonomic filtering
df.taxon.richness <- length(unique(df$tax.assignment))

#Set analysis taxonomic rank
taxon <- "tax.assignment"
grp <- c("sample", taxon)

#Create separate dataframes
df$read_abundance <- as.integer(df$read_abundance)
df.long = df %>% #Long format for plotting
  filter_at(vars(c(lat, long)), all_vars(!is.na(.)))

df.table <- df %>%
  select(c(sample, read_abundance, taxon)) %>%
  group_by(across(all_of(grp))) %>%
  summarise(read_abundance = sum(read_abundance))
df.wide <- df.table %>% #Wide format to create community data
  pivot_wider(
    names_from = taxon,
    values_from = read_abundance
  )
df.wide <- df.wide %>%
  mutate_at(c(2:ncol(df.wide)), ~replace_na(., 0)) %>%
  as.data.frame()

#Create community datasets
comm <- df.wide[,-1] #Read abundance
rownames(comm) <- df.wide$sample
comm.prop <- decostand(comm, method = "total") #Read abundance proportions
comm.log <- decostand(comm, method = "log") #Log-transformed read abundance
comm.log.prop <- decostand(comm.log, method = "total") #Log-transformed read abundance proportions
comm.pa <- decostand(comm, method = "pa") #Presence/absence
comm.hell <- decostand(comm, method = "hellinger") #Hellinger transformation (sqrt of proportions)
comm.ra = comm

remove(comm)
# Set analysis level #
comm = comm.log

#Create metadata for wide format
comm.meta <- df.wide %>%
  left_join(.,
            meta,
            by = "sample") %>%
  select(sample, triplicate:Comments)

#Shannon Diversity
SW.out <- diversity(comm.ra, index = "shannon")
SW.table <- data.frame(SW.out)
SW.table$sample <- rownames(SW.table)
SW.table <- SW.table %>%
  rename(Shannon_index = SW.out)

comm.meta <- comm.meta %>%
  left_join(.,
            SW.table,
            by = "sample")

#Permutational analysis of variance
set.seed(12031996)
prop.adonis <- adonis2(comm ~ depth_cat + month,
                       data = comm.meta,
                       permutations = 999,
                       method = "bray",
                       by = "terms")
prop.adonis

#Ordination/Nonmetric multi-dimensional scaling
set.seed(12031996)
comm.mds <- metaMDS(comm, dist = "bray", k = 5, trymax = 100)
comm.mds$stress

#Environmental loadings
comm.meta$depth <- as.integer(comm.meta$depth)
env <- comm.meta %>%
  select(depth, Shannon_index)
en = envfit(comm.mds, env, permutations = 999, na.rm = TRUE)
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

data_scores <- as.data.frame(scores(comm.mds, "sites"))
sp_scores <- as.data.frame(scores(comm.mds, "species"))
sp_scores$species <- rownames(sp_scores)
scrs <- scores(comm.mds, display = "sites") %>%
  as.data.frame()
scrs$sample <- rownames(scrs)
scrs <- scrs %>%
  as.data.frame() %>%
  left_join(.,
            comm.meta,
            by = "sample")

cent <- aggregate(cbind(NMDS1, NMDS2) ~ depth_cat, data = scrs, FUN = mean)
segs <- merge(scrs, setNames(cent, c("depth_cat","oNMDS1","oNMDS2")),
              by = "depth_cat", sort = FALSE)
NMDS.plot = ggplot(scrs, aes(x = NMDS1, y = NMDS2, color = depth_cat)) +
  geom_segment(data = segs, mapping = aes(xend = oNMDS1, yend = oNMDS2), alpha = 0.3) +
  geom_point(data = cent, size = 4) + 
  geom_point(size = 2, alpha = 0.3) + coord_fixed() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, linewidth = 1, alpha = 1, colour = "black") +
  #geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
  #           shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  #geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04), 
  #          label = row.names(en_coord_cat), colour = "navy", fontface = "bold") + 
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  ggtitle("eDNA NMDS Ordination") +
  xlab("NMDS 1") + ylab("NMDS 2") +
  labs(subtitle = "Bray-Curtis Dissimilarity on Log-Transformed Reads") +
  labs(color = "Depth Categories")
NMDS.plot

plot1 = df.long %>%
  group_by(depth_cat, phylum) %>%
  summarise(read_abundance = sum(((read_abundance)))) %>%
  arrange(desc(read_abundance)) %>%
  ggplot(aes(x = as.factor(depth_cat), y = read_abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = "black", position = "fill") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  #facet_grid(vars(depth_cat)) +
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  coord_flip() + scale_x_discrete(limits = rev) +
  xlab("") +
  ylab("eDNA Read Abundance (RA proportions)") +
  labs(fill = "eDNA Detected Taxa") + ggtitle("A")

plot2 = df.long %>% 
  group_by(depth_cat, phylum) %>%
  summarise(read_abundance = sum((log1p(read_abundance)))) %>%
  arrange(desc(read_abundance)) %>%
  ggplot(aes(x = as.factor(depth_cat), y = read_abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = "black", position = "fill") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  #facet_grid(vars(depth_cat)) +
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  coord_flip() + scale_x_discrete(limits = rev) +
  xlab("") +
  ylab("eDNA Log RA (LRA, proportions)") +
  labs(fill = "eDNA Detected Taxa") + ggtitle("B")

log.plot.01 = df.long %>% 
  group_by(depth_cat, phylum) %>%
  summarise(read_abundance = sum((log1p(read_abundance)))) %>%
  arrange(desc(read_abundance)) %>%
  ggplot(aes(x = as.factor(depth_cat), y = read_abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = "black", position = "fill", width = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  theme(
    panel.background = element_rect(fill = "#D8EDE6"),
    axis.line = element_line(colour = "black")
  ) +
  coord_flip() + scale_x_discrete(limits = rev) +
  labs(
    title = "Challenger Deep eDNA (Log-Transformed Reads)",
    subtitle = "A. Taxa by depth category, LRA proportions",
    y = "LRA; proportions", x = "", fill = "eDNA-Observed Phyla"
  ) +
  theme(legend.position = "none")
log.plot.02 = df.long %>% 
  group_by(depth_cat, phylum) %>%
  summarise(read_abundance = sum((log1p(read_abundance)))) %>%
  arrange(desc(read_abundance)) %>%
  ggplot(aes(x = as.factor(phylum), y = read_abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(vars(depth_cat)) +
  theme(
    panel.background = element_rect(fill = "#D8EDE6"),
    axis.line = element_line(colour = "black")
  ) +
  coord_flip() + scale_x_discrete(limits = rev) +
  labs(
    title = "",
    subtitle = "B. Taxa by depth category, LRA",
    x = "", y = "LRA", fill = "eDNA-Observed Phyla"
  ) +
  theme(legend.position = "none")
ggarrange(log.plot.01, log.plot.02,
          common.legend = TRUE,
          legend = "bottom")

df.table <- df.table[,-1] %>%
  group_by(tax.assignment) %>%
  summarise(read_abundance = sum(read_abundance)) %>%
  arrange(desc(read_abundance)) %>%
  slice(1:15)

simper.general <- summary(simper(comm, permutations = 999), ordered = TRUE)
simper.depth <- summary(with(comm.meta, simper(comm, depth_cat, permutations = 999)))

simper.depth.table = simper.depth$`01. Epipelagic (0-200m)_04. Benthos (10km)` %>%
  as.data.frame() %>%
  mutate(
    species = rownames(.)
  ) %>%
  select(
    c(
      species,
      average,
      p
    )
  ) %>%
  head(
    n = 15
  ) %>%
  as_tibble() %>%
  tibble::remove_rownames() %>%
  kable(
    "simple",
    col.names = c(
      "Species",
      "Avg. Dissim. Contribution",
      "Permuted P-Value"
    )
  )

deep.df = df %>% filter(depth_cat != "01. Epipelagic (0-200m)")
plot3 = df %>%
  group_by(tax.assignment, depth_cat) %>%
  summarise(read_abundance = sum((read_abundance))) %>%
  ggplot(
    aes(
      x = reorder(tax.assignment, -(read_abundance), sum),
      y = (read_abundance),
      fill = depth_cat
    )
  ) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_viridis(discrete = TRUE, option = "plasma") +
  coord_flip() + scale_x_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("eDNA Taxa") +
  ylab("RA") + labs(fill = "Depth Category") + ggtitle(
    "Marianas Trench eDNA Taxa Abundance by Depth"
  )

eDNA_barplots = ggarrange(
  plot1, plot2, common.legend = TRUE, legend = "right", nrow = 2
)
annotate_figure(eDNA_barplots, top = text_grob("mtCOI eDNA Phyla by Depth", 
                                               color = "black", size = 14))
plot3


# Redefine df.wide and df.long #
remove(df.wide)
remove(df.long)
remove(comm)

comm = comm.ra
comm$sample = rownames(comm.ra)
df.wide = comm %>% left_join(
  comm.meta,
  by = "sample"
) %>%
  select(sample, triplicate:Shannon_index, everything())
df.long = df.wide %>%
  pivot_longer(
    cols = 18:ncol(df.wide),
    names_to = "tax.assignment",
    values_to = "comm.index"
  )

df.long = df.long %>% 
  filter_at(vars(c(tax.assignment, comm.index,
                   depth_cat)), all_vars(!is.na(.)))

grp = c("tax.assignment", "depth_cat")
df.long %>%
  group_by(across(all_of(grp))) %>%
  summarise(comm.index = sum((comm.index))) %>%
  filter(comm.index > 0) %>%
  write.csv(
    "taxa_by_depth.csv"
  )

##### End script ####
