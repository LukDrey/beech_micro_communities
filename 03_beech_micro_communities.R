#########################################################################
## Tree size drives diversity and community structure of microbial     ##
## communities on the bark of beech (Fagus sylvatica)                  ##
#########################################################################

# This is the code accompanying the analysis for our Paper
# in Frontiers in Forests and Global Change - Temperate and Boreal Forests.


#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(here); packageVersion("here")

library(dplyr); packageVersion("dplyr")

library(decontam); packageVersion("decontam")

library(phyloseq); packageVersion("phyloseq")

library(ggplot2); packageVersion("ggplot2")

library(lulu); packageVersion("lulu")

library(Biostrings); packageVersion("Biostrings")

library(microbiome); packageVersion("microbiome")

library(fantaxtic); packageVersion("fantaxtic")

library(paletteer); packageVersion("paletteer")

library(igraph); packageVersion("igraph")

library(SpiecEasi); packageVersion("SpiecEasi")

library(rgexf); packageVersion("rgexf")

library(tibble); packageVersion("tibble")

library(vegan); packageVersion("vegan")

library(ggpubr); packageVersion("ggpubr")

##----------------------------------------------------------------
##                        Custom Functions                       -
##----------------------------------------------------------------

# This is a wrapper function round ape{moran.I()} to allow comparison of multiple,
# within plot, tree distances as described in the Material & Methods section.

moran_wrapper <- function(variable, coord_x, coord_y) {
  require(ape)
  dst <- as.matrix(dist(cbind(coord_x, coord_y)))
  dst_inv <- 1/dst
  diag(dst_inv) <- 0 
  morans_I <- ape::Moran.I(variable, dst_inv)
  if(morans_I$p.value > 0.05) {
    print('No autospatial correlation')
  } else {
    print('Autospatial correlation!')
  }
  return(morans_I)
}

# Setting the colors for the community bar graphs. 
my_cols <- paletteer_d('ggsci::default_igv')

# Setting the colors for the raincloud plots comparing alpha diversity. 
box_cols <- c( "#5050FFFF", "#CE3D32FF", "#008099FF")

##################################################################
##                          Section 2                           ##
##                Decontam Run and LULU curation                ##
##################################################################

##---------
##  Algae  
##---------

###
# Decontam
###

# Load ASV table for algae (available as supplementary data).
algae_otu <- read.csv(here('asv_table_algae.csv'), header = F, sep = ',')

# Set column names from first row. 
colnames(algae_otu) <- algae_otu[1,]
algae_otu <- algae_otu[-1,]

class(algae_otu)
algae_otu <- as.matrix(algae_otu)

# Naming the Samples. 
colnames(algae_otu) <- paste0("Sample_", colnames(algae_otu))

# set ASV ID as the rownames.
rownames(algae_otu) <- algae_otu[,1]
algae_otu <- algae_otu[,-1]

class(algae_otu) <- "numeric" 

# Load the DNA concentration necessary for decontam (available as supplementary data). 
algae_conc <- read.csv(here('algae_decontam_conc.csv'), header = T, sep =',')
rownames(algae_conc) <- algae_conc$sample

# Make the parts of the phyloseq object.
ASV_algae_decontam <- otu_table(algae_otu, taxa_are_rows = TRUE)
sampledata_algae_decontam <- sample_data(algae_conc)

# Combine with phyloseq
ps_algae_decontam <- phyloseq(ASV_algae_decontam, sampledata_algae_decontam)
ps_algae_decontam

# Check the library sizes of the Samples and Controls.
df_algae <- as.data.frame(sample_data(ps_algae_decontam)) # Put sample_data into a ggplot-friendly data.frame
df_algae$LibrarySize <- sample_sums(ps_algae_decontam)
df_algae <- df_algae[order(df_algae$LibrarySize),]
df_algae$Index <- seq(nrow(df_algae))
ggplot(data=df_algae, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) +
  geom_point()

# Decontam with combined frequency and prevalence approach. 
sample_data(ps_algae_decontam)$is.neg <- sample_data(ps_algae_decontam)$Sample_or_Control == "Control Sample"
contam_algae_combi <- isContaminant(ps_algae_decontam,
                                    method = 'combined',
                                    neg = 'is.neg',
                                    conc = 'quant_reading')
table(contam_algae_combi$contaminant)
which(contam_algae_combi$contaminant)

# Trim the identified contaminants from the phyloseq object.
# In this case there are no identified contaminants. So the non-contaminant phyloseq object,
# will be the same as before.
ps_algae_noncontam <- ps_algae_decontam
ps_algae_noncontam

###
#LULU curation
###

# Remove taxa without reads.
ps_algae_noncontam_pruned <- prune_taxa(taxa_sums(ps_algae_noncontam) > 0,
                                        ps_algae_noncontam)
ASV_table_algae <- as.data.frame(otu_table(ps_algae_noncontam_pruned))
algae_matchlist <- read.table(here('match_list_algae.txt'),
                              header = F,
                              as.is = T,
                              stringsAsFactors = F)

# Run the LULU algorithm. 

ASV_table_algae_cur <- lulu(ASV_table_algae, algae_matchlist)


ASV_table_algae_cur$curated_count
ASV_table_algae_cur$discarded_count

##---------
##  Fungi  
##---------

###
# Decontam
###

# Load ASV table for fungi (available as supplementary data).
fungi_otu <- read.csv(here('asv_table_fungi.txt'), header = F, sep = ',')

# Set column names from first row. 
colnames(fungi_otu) <- fungi_otu[1,]
fungi_otu <- fungi_otu[-1,]

class(fungi_otu)
fungi_otu <- as.matrix(fungi_otu)


# Naming the Samples. 
colnames(fungi_otu) <- paste0("Sample_", colnames(fungi_otu))

# set ASV ID as the rownames.
rownames(fungi_otu) <- fungi_otu[,1]
fungi_otu <- fungi_otu[,-1]

class(fungi_otu) <- "numeric" 

# Load the DNA concentration necessary for decontam (available as supplementary data). 
fungi_conc <- read.csv(here('fungi_decontam_conc.csv'), header = T, sep =',')
rownames(fungi_conc) <- fungi_conc$sample

# Make the parts of the phyloseq object.
ASV_fungi_decontam <- otu_table(fungi_otu, taxa_are_rows = TRUE)
sampledata_fungi_decontam <- sample_data(fungi_conc)

# Combine with phyloseq
ps_fungi_decontam <- phyloseq(ASV_fungi_decontam, sampledata_fungi_decontam)
ps_fungi_decontam

# Check the library sizes of the Samples and Controls.
df_fungi <- as.data.frame(sample_data(ps_fungi_decontam)) # Put sample_data into a ggplot-friendly data.frame
df_fungi$LibrarySize <- sample_sums(ps_fungi_decontam)
df_fungi <- df_fungi[order(df_fungi$LibrarySize),]
df_fungi$Index <- seq(nrow(df_fungi))
ggplot(data=df_fungi, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) +
  geom_point()

# Decontam with combined frequency and prevalence approach. 
sample_data(ps_fungi_decontam)$is.neg <- sample_data(ps_fungi_decontam)$Sample_or_Control == "Control Sample"
contam_fungi_combi <- isContaminant(ps_fungi_decontam,
                                    method = 'combined',
                                    neg = 'is.neg',
                                    conc = 'quant_reading')
table(contam_fungi_combi$contaminant)
which(contam_fungi_combi$contaminant)

# Trim the identified contaminants from the phyloseq object 
ps_fungi_noncontam <- prune_taxa(!contam_fungi_combi$contaminant,
                                 ps_fungi_decontam)
ps_fungi_noncontam

###
#LULU curation
###

# Remove taxa without reads.
ps_fungi_noncontam_pruned <- prune_taxa(taxa_sums(ps_fungi_noncontam) > 0,
                                        ps_fungi_noncontam)
ASV_table_fungi <- as.data.frame(otu_table(ps_fungi_noncontam_pruned))
fungi_matchlist <- read.table(here('match_list_fungi.txt'),
                              header = F,
                              as.is = T,
                              stringsAsFactors = F)

# Run the LULU algorithm. 

ASV_table_fungi_cur <- lulu(ASV_table_fungi, fungi_matchlist)


ASV_table_fungi_cur$curated_count
ASV_table_fungi_cur$discarded_count

##---------
## Bacteria  
##---------

###
# Decontam
###

# Load ASV table for bacteria (available as supplementary data).
bacteria_otu <- read.csv(here('asv_table_bacteria.txt'), header = F, sep = '\t')

# Set column names from first row. 
colnames(bacteria_otu) <- bacteria_otu[1,]
bacteria_otu <- bacteria_otu[-1,]

class(bacteria_otu)
bacteria_otu <- as.matrix(bacteria_otu)

# Naming the Samples. 
colnames(bacteria_otu) <- paste0("Sample_", colnames(bacteria_otu))

# set ASV ID as the rownames.
rownames(bacteria_otu) <- bacteria_otu[,1]
bacteria_otu <- bacteria_otu[,-1]

class(bacteria_otu) <- "numeric" 

# Load the DNA concentration necessary for decontam (available as supplementary data). 
bacteria_conc <- read.csv(here('bacteria_decontam_conc.csv'), header = T, sep =',')
rownames(bacteria_conc) <- bacteria_conc$sample

# Make the parts of the phyloseq object.
ASV_bacteria_decontam <- otu_table(bacteria_otu, taxa_are_rows = TRUE)
sampledata_bacteria_decontam <- sample_data(bacteria_conc)

# Combine with phyloseq
ps_bacteria_decontam <- phyloseq(ASV_bacteria_decontam, sampledata_bacteria_decontam)
ps_bacteria_decontam

# Check the library sizes of the Samples and Controls.
df_bacteria <- as.data.frame(sample_data(ps_bacteria_decontam)) # Put sample_data into a ggplot-friendly data.frame
df_bacteria$LibrarySize <- sample_sums(ps_bacteria_decontam)
df_bacteria <- df_bacteria[order(df_bacteria$LibrarySize),]
df_bacteria$Index <- seq(nrow(df_bacteria))
ggplot(data=df_bacteria, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) +
  geom_point()

# Decontam with combined frequency and prevalence approach. 
sample_data(ps_bacteria_decontam)$is.neg <- sample_data(ps_bacteria_decontam)$Sample_or_Control == "Control Sample"
contam_bacteria_combi <- isContaminant(ps_bacteria_decontam,
                                       method = 'combined',
                                       neg = 'is.neg',
                                       conc = 'quant_reading')
table(contam_bacteria_combi$contaminant)
which(contam_bacteria_combi$contaminant)

# Trim the identified contaminants from the phyloseq object. 
ps_bacteria_noncontam <- prune_taxa(!contam_bacteria_combi$contaminant,
                                    ps_bacteria_decontam)
ps_bacteria_noncontam

###
#LULU curation
###

# Remove taxa without reads.
ps_bacteria_noncontam_pruned <- prune_taxa(taxa_sums(ps_bacteria_noncontam) > 0,
                                           ps_bacteria_noncontam)
ASV_table_bacteria <- as.data.frame(otu_table(ps_bacteria_noncontam_pruned))
bacteria_matchlist <- read.table(here('match_list_bacteria.txt'),
                                 header = F,
                                 as.is = T,
                                 stringsAsFactors = F)

# Run the LULU algorithm. 

ASV_table_bacteria_cur <- lulu(ASV_table_bacteria, bacteria_matchlist)


ASV_table_bacteria_cur$curated_count
ASV_table_bacteria_cur$discarded_count



##################################################################
##                          Section 3                           ##
##         Data Loading, Data Cleaning and Initialization       ##
##################################################################

##----------------------------------------------------------------
##                        Metadata table                         -
##----------------------------------------------------------------

# Loading Forest Management Data (available at https://www.bexis.uni-jena.de/PublicData/PublicData.aspx?DatasetId=16466)
formi <- read.table(here('16466.txt'), header = T, sep = '\t')

# renaming the PlotID to match other columns, extracting only data for Hainich
HAI_formi <- formi %>% 
  dplyr::rename(., Plot_ID = EP_Plotid) %>%
  filter(.,grepl('HEW', Plot_ID)) %>%
  select(., c("Plot_ID", "ForMI")) %>% 
  add_row(Plot_ID = "HEW51", ForMI = NA)

# Loading tree measurements (available as supplementary data)
tree_meas <- read.csv(here('beech_micro_communities_metadata.csv'), sep = ',', header = T)

# Remove the corner points. Necessary for calculation of x & y coordinates within plots
# but do not represent trees. 
tree_meas <- filter(tree_meas, sample_num != 'corner')

# Join the tables to get one metadata table.
metadata <- inner_join(tree_meas, HAI_formi, by = 'Plot_ID')

# Create management intensity categories.
metadata$intensity_formi <- ifelse(metadata$ForMI >= 1, 'high', 'low') 

# Create tree size categories. 

metadata$size_cat <- ifelse(metadata$DBH_cm < 15, 'small',
                            ifelse(metadata$DBH_cm >= 30, 'big', 'medium'))

# Set rownames to allow turning dataframe into sample_data object in phyloseq. 
rownames(metadata) <- metadata$Sample_ID
metadata$sample_num <- NULL
metadata$Sample_ID <- NULL

##----------------------------------------------------------------
##                          ASV tables                           -
##----------------------------------------------------------------

##---------
##  Algae  
##---------

# Keep only samples that do represent real tree swabs. Cut Controls.
asv_algae <- ASV_table_algae_cur$curated_table %>% dplyr::select(num_range('Sample_', 1:96))


##---------
##  Bacteria  
##---------

# Keep only samples that do represent real tree swabs. Cut Controls.
asv_bacteria <- ASV_table_bacteria_cur$curated_table %>% dplyr::select(num_range('Sample_', 1:96))

##---------
##  Fungi  
##---------

# Keep only samples that do represent real tree swabs. Cut Controls. 
asv_fungi <- ASV_table_fungi_cur$curated_table %>% dplyr::select(num_range('Sample_', 1:96))

##---------------------------------------------------------------
##                        Taxonomy Tables                       -
##---------------------------------------------------------------

##---------
##  Algae  
##---------

# Load the algal taxonomy table as obtained from Seed2 and blast. 
tax_algae <- read.csv(here('seed_taxonomy.csv'), header = T, sep = "\t")
tax_algae$ASV_ID <- gsub("_A", "", tax_algae$ASV_ID)

# Load algal sequences. 
algae_seqs_fasta <- readDNAStringSet(here('ASVs_algae.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_algae <- names(algae_seqs_fasta)
sequence_algae <- paste(algae_seqs_fasta)
algae_rep_seqs <- data.frame(seq_name_algae, sequence_algae)

algae_rep_seqs <- algae_rep_seqs %>% dplyr::rename(ASV_ID = seq_name_algae)

# Join the taxonomy table and the representative sequences
tax_clean_algae <- left_join(tax_algae, algae_rep_seqs, by = 'ASV_ID')

# Set rownames.
row.names(tax_clean_algae) <- tax_clean_algae$ASV_ID

# Remove unwanted columns. 

tax_clean_algae$ASV_ID <- NULL
tax_clean_algae$sequence_algae <- NULL

# Remove Fungi. 
tax_clean_algae <- dplyr::filter(tax_clean_algae, kingdom != 'Fungi')

##---------
##  Bacteria  
##---------

# Load the bacterial taxonomy table. 
# (Available as supplementary data)
tax_bacteria <- readRDS(here('tax_table_bacteria.rds'))
tax_bacteria <- as.data.frame(tax_bacteria) %>% tibble::rownames_to_column('sequence')
tax_bacteria <- tax_bacteria %>%
  dplyr::rename(sequence_bacteria = sequence)

# Load bacterial sequences. 
bacteria_seqs_fasta <- readDNAStringSet(here('ASVs_bacteria.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_bacteria <- names(bacteria_seqs_fasta)
sequence_bacteria <- paste(bacteria_seqs_fasta)
bacteria_rep_seqs <- data.frame(seq_name_bacteria, sequence_bacteria)

# Join the taxonomy table and the representative sequences
tax_clean_bacteria <- left_join(tax_bacteria, bacteria_rep_seqs, by = 'sequence_bacteria')

# Remove the underscore to make the ASV ID the same as in the ASV table. 
tax_clean_bacteria$seq_name_bacteria <- gsub("_", "", tax_clean_bacteria$seq_name_bacteria)

# Set rownames.
row.names(tax_clean_bacteria) <- tax_clean_bacteria$seq_name_bacteria

# Remove any Chloroplast Sequences
bacteria_tax_fin_raw <- dplyr::filter(tax_clean_bacteria, Order != 'Chloroplast')
bacteria_tax_fin_raw <- dplyr::filter(tax_clean_bacteria, Family != 'Mitochondria')

bacteria_tax_fin_raw$seq_name_bacteria <- NULL
bacteria_tax_fin_raw$sequence_bacteria <- NULL

##---------
##  Fungi  
##---------

# Load the fungal taxonomy table.
# (Available as supplementary data)
tax_fungi <- readRDS(here('tax_table_fungi.rds'))
tax_fungi <- as.data.frame(tax_fungi) %>% tibble::rownames_to_column('sequence')
tax_fungi <- tax_fungi %>%
  dplyr::rename(sequence_fungi = sequence)

# Load the fungal reads.
fungi_seqs_fasta <- readDNAStringSet(here('ASVs_fungi.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_fungi <- names(fungi_seqs_fasta)
sequence_fungi <- paste(fungi_seqs_fasta)
fungi_rep_seqs <- data.frame(seq_name_fungi, sequence_fungi)

# Join the taxonomy table and the representative sequences
tax_clean_fungi <- left_join(tax_fungi, fungi_rep_seqs, by = 'sequence_fungi')

# Split the taxonomy into different columns of taxonomic levels.

fungi_tax_fin <- tidyr::separate(tax_clean_fungi, Kingdom, c(NA, 'Kingdom') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Phylum, c(NA, 'Phylum') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Class, c(NA, 'Class') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Order, c(NA, 'Order') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Family, c(NA, 'Family') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Genus, c(NA, 'Genus') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Species, c(NA, 'Species') , sep = '__')

# Set the rownames.
rownames(fungi_tax_fin) <- fungi_tax_fin$seq_name_fungi 
fungi_tax_fin$seq_name_fungi <- NULL 

# Remove the sequence column.
fungi_tax_fin$sequence_fungi <- NULL

# Remove sequences that could not be assigned at Phylum level.
fungi_tax_fin_assigned <- dplyr::filter(fungi_tax_fin, !is.na(Phylum))

##---------------------------------------------------------------
##               Create the Phyloseq Objects                    -
##---------------------------------------------------------------

##---------
##  Algae  
##---------

# Transform dataframe to matrix.
asvmat_algae <- as.matrix(asv_algae) 

# Transform dataframe to matrix.
taxmat_algae <- as.matrix(tax_clean_algae) 

# Create ASV table for phyloseq.
ASV_ALG <- otu_table(asvmat_algae, taxa_are_rows = T)

# Create taxonomy table for phyloseq.
TAX_ALG <- tax_table(taxmat_algae) 

# Metadata for phyloseq. 
sampledata <- sample_data(metadata) 

# Combine in phyloseq object. 
phy_algae <- phyloseq(ASV_ALG, TAX_ALG, sampledata) 
phy_algae

# Seed2 gives slightly different Taxonomy levels. We correct this here.

colnames(tax_table(phy_algae)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
rank_names(phy_algae)

##---------
##  Bacteria  
##---------

# Transform dataframe to matrix.
asvmat_bacteria <- as.matrix(asv_bacteria) 

# Transform dataframe to matrix.
taxmat_bacteria <- as.matrix(bacteria_tax_fin_raw) 

# Create OTU table for phyloseq.
ASV_BAC <- otu_table(asvmat_bacteria, taxa_are_rows = T) 

# Create taxonomy table for phyloseq.
TAX_BAC <- tax_table(taxmat_bacteria) 

# Metadata for phyloseq.
sampledata <- sample_data(metadata)  

# Combine in phyloseq object. 
phy_bacteria <- phyloseq(ASV_BAC, TAX_BAC, sampledata) 
phy_bacteria

##---------
##  Fungi  
##---------

# Transform dataframe to matrix.
asvmat_fungi <- as.matrix(asv_fungi)

# Transform dataframe to matrix.
taxmat_fungi <- as.matrix(fungi_tax_fin_assigned) 

# Create ASV table for phyloseq.
ASV_FUN <- otu_table(asvmat_fungi, taxa_are_rows = T) 

# Create taxonomy table for phyloseq.
TAX_FUN <- tax_table(taxmat_fungi) 

# Metadata for phyloseq. 
sampledata <- sample_data(metadata) 

# Combine in phyloseq object. 
phy_fungi <- phyloseq(ASV_FUN, TAX_FUN, sampledata) 
phy_fungi

#################################################################
##                          Section 3                          ##
##                   Intra-Group Diversities                   ##
#################################################################

##---------
##  Algae  
##---------

###
#Alpha Diversity 
###

# Calculate the Shannon Diversity.
shannon_alg <- estimate_richness(phy_algae, split = T, measures = 'Shannon')

# Make a new sample_data object and merge with existing phyloseq object.
new_var_alg <- sample_data(shannon_alg)
phy_algae <- merge_phyloseq(phy_algae, new_var_alg)

# Turn the sample_data object into a dataframe for ANOVA.
alg_eval_df <- data.frame(sample_data(phy_algae))

# Conduct ANOVA and verify results with Tukey HSD.
anova_alg_alpha <- aov(Shannon ~ intensity_formi + size_cat, data = alg_eval_df)
summary(anova_alg_alpha)

TukeyHSD(anova_alg_alpha)

###
#Raincloud Plots of Alpha Diversity
###

# Code was adopted from Cédric Scherer. You can find a tutorial at: 
# https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/

# Extract the sample data from the phyloseq object.
algae_sampledata <- sample_data(phy_algae)

# Comparing the Management Categories.
alg_rain_manage <- ggplot(algae_sampledata,
                          aes(x = intensity_formi,
                              y = Shannon,
                              fill = intensity_formi,
                              colour = intensity_formi)) + 
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .5, 
    ## move geom to the right
    justification = -.2, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) + 
  geom_boxplot(
    width = .14, 
    ## remove outliers
    outlier.color = NA,  
    alpha = 0.5,
    size = 0.1
  ) +
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .5, 
    size = 0.5
  )  +
  ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        axis.line = element_line(colour = "black", size = .2),
        axis.ticks = element_line(size = 0.2),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        legend.position = 'none', 
        plot.title = element_text(vjust = 0, hjust = 0.007),
        plot.subtitle = element_text(vjust = -0.5, hjust = 0.03, size = 9))+
  ylab( 'Alpha Diversity (Shannon)') +
  xlab('Management Intensity') + 
  labs(subtitle = '(A)') +
  scale_fill_manual(values = box_cols) + 
  scale_colour_manual(values = box_cols)
alg_rain_manage

# Comparing the Tree Sizes.
alg_rain_size <- ggplot(algae_sampledata,
                        aes(x = size_cat,
                            y = Shannon, 
                            fill = size_cat,
                            colour = size_cat)) + 
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .5, 
    ## move geom to the right
    justification = -.25, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) + 
  geom_boxplot(
    width = .15, 
    ## remove outliers
    outlier.color = NA,  
    alpha = 0.5, 
    size = 0.1
  ) +
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .3, 
    ## add some transparency
    alpha = .5,
    size = .5
  )  +
  ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        axis.line = element_line(colour = "black", size = .2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none', 
        plot.subtitle = element_text(vjust = -0.5, hjust = .03, size = 9)) +
  ylab( 'Alpha Diversity (Shannon)') +
  xlab('Tree Size') + 
  labs(subtitle = '(D)') +
  scale_fill_manual(values = box_cols) + 
  scale_colour_manual(values = box_cols) 
alg_rain_size

###
#Spatial Autocorrelation of alpha diversity
###

# Make a list that includes the sample data split by plots.
plotwise_alg <- sample_data(phy_algae) %>% 
  dplyr::group_by(Plot_ID) %>%
  group_split()

# Initiate empty lists to fill.
results_list_alg <- list()
plot_ID <- NULL

# Calculate Morans I for the trees within each plot. 
for (j in 1:length(plotwise_alg)) {
  results <- moran_wrapper(plotwise_alg[[j]]$Shannon,
                           plotwise_alg[[j]]$x_coord,
                           plotwise_alg[[j]]$y_coord)
  results_list_alg[[j]] <- results
  plot_ID[j] <- unique(plotwise_alg[[j]]$Plot_ID)
}

# Set the Plot ID as a name for the list.
names(results_list_alg) <- plot_ID  

# Combine the results in a dataframe.
morans_i_alg_plots <- data.frame(do.call(rbind.data.frame, results_list_alg))
morans_i_alg_plots

###
#Community Barplots
###

# Aggregate the taxa at order level. This is to allow better visualization, 
# while still retaining the most information.
phy_order_alg <- phy_algae %>%
  aggregate_taxa(level = "Order")

# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_alg_ord_top25 <- get_top_taxa(phy_order_alg, 24,
                                  relative = TRUE,
                                  other_label = "Others")

# Transform the subset dataset to compositional (relative) abundances.
phy_alg_ord_top25_plot <-  phy_alg_ord_top25 %>%
  aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
taxa_names(phy_alg_ord_top25_plot) <- tax_table(phy_alg_ord_top25_plot)[, 4]

# Turn the columns necessary for the ploting into sorted factors to set the order they are plotted in.
sample_data(phy_alg_ord_top25_plot)$intensity_formi_f <- factor(sample_data(phy_alg_ord_top25_plot)$intensity_formi, 
                                                                levels = c("low", "high")) 
sample_data(phy_alg_ord_top25_plot)$Plot_ID <- factor(sample_data(phy_alg_ord_top25_plot)$Plot_ID,
                                                      levels = c('HEW5', 'HEW11', 'HEW20', 'HEW27', 'HEW34','HEW35','HEW36','HEW37',
                                                                 'HEW8','HEW26','HEW28','HEW31','HEW32','HEW33','HEW43','HEW49')) 
sample_data(phy_alg_ord_top25_plot)$size_cat <- factor(sample_data(phy_alg_ord_top25_plot)$size_cat,
                                                       levels = c('big', 'medium', 'small')) 
sample_names(phy_alg_ord_top25_plot) <- c(paste0('T', 1:96))

# Sort the taxa names alphabetically. 
taxa_names_alg_ord <- sort(taxa_names(phy_alg_ord_top25_plot))

taxa_names_alg_ord <- c("Chaetophorales", "Chlorellales", "Klebsormidiales",
                        "Prasiolales", "Sphaeropleales", "Trebouxiales", "Unknown")

# Custom plotting to make a nice stacked barplot. 
alg_ord_plots <- phy_alg_ord_top25_plot %>%
  plot_composition( group_by =  'Plot_ID', otu.sort = taxa_names_alg_ord, sample.sort = 'size_cat') +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Order') +
  guides(fill = guide_legend(title.position = 'top')) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_text(colour = "black", size = 5, hjust = 1),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black')) + 
  xlab('Tree ID') +
  ylab('Relative Abundance')  +
  scale_y_continuous(label = scales::percent)  + 
  labs( subtitle = '(A)')

alg_ord_plots


##---------
##  Bacteria  
##-------------------------------------------------------------------------------------------------------

# Calculate the Shannon Diversity.
shannon_bac <- estimate_richness(phy_bacteria, split = T, measures = 'Shannon')

# Make a new sample_data object and merge with existing phyloseq object.
new_var_bac <- sample_data(shannon_bac)
phy_bacteria <- merge_phyloseq(phy_bacteria, new_var_bac)

# Turn the sample_data object into a dataframe for ANOVA.
bac_eval_df <- data.frame(sample_data(phy_bacteria))

# Conduct ANOVA and verify results with Tukey HSD.
anova_bac_alpha <- aov(Shannon ~ intensity_formi + size_cat, data = bac_eval_df)
summary(anova_bac_alpha)

TukeyHSD(anova_bac_alpha)

###
#Raincloud Plots of Alpha Diversity
###

# Code was adopted from Cédric Scherer. You can find a tutorial at: 
# https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/

# Compare Management Categories.

bacteria_sampledata <- sample_data(phy_bacteria)

bac_rain_manage <- ggplot(bacteria_sampledata, aes(x = intensity_formi, y = Shannon, fill = intensity_formi, colour = intensity_formi)) + 
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .7, 
    ## move geom to the right
    justification = -.15, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) + 
  geom_boxplot(
    width = .14, 
    ## remove outliers
    outlier.color = NA, ## `outlier.shape = NA` works as well 
    alpha = 0.5, 
    size = 0.1
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .5,
    size = 0.5
  )  +
  ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        axis.line = element_line(colour = "black", size = .2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none', 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9),
        plot.subtitle = element_text(vjust = -0.5, hjust = 0.03, size = 9)) +
  xlab('Management Intensity') + 
  labs(subtitle = '(B)') +
  scale_fill_manual(values = box_cols) + 
  scale_colour_manual(values = box_cols)
bac_rain_manage

# Compare Tree Sizes.

bac_rain_size <- ggplot(bacteria_sampledata, aes(x = size_cat, y = Shannon, fill = size_cat, colour = size_cat)) + 
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .5, 
    ## move geom to the right
    justification = -.25, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) + 
  geom_boxplot(
    width = .15, 
    ## remove outliers
    outlier.color = NA, ## `outlier.shape = NA` works as well 
    alpha = 0.5, 
    size = 0.1
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .3, 
    ## add some transparency
    alpha = .5,
    size = 0.5
  )  +
  ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.text = element_text(colour = "black", size = 9),
        axis.line = element_line(colour = "black", size = .2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none', 
        plot.subtitle = element_text(vjust = -0.5, hjust = .03, size = 9)) +
  xlab('Tree Size') + 
  labs(subtitle = '(E)') +
  scale_fill_manual(values = box_cols) + 
  scale_colour_manual(values = box_cols) 
bac_rain_size

###
#Spatial Autocorrelation of alpha diversity
###

# Make a list that includes the sample data split by plots.
plotwise_bac <- sample_data(phy_bacteria) %>% 
  dplyr::group_by(Plot_ID) %>%
  group_split()

# Initiate empty lists to fill.
results_list_bac <- list()
plot_ID <- NULL

# Calculate Morans I for the trees within each plot. 
for (j in 1:length(plotwise_bac)) {
  results <- moran_wrapper(plotwise_bac[[j]]$Shannon,
                           plotwise_bac[[j]]$x_coord,
                           plotwise_bac[[j]]$y_coord)
  results_list_bac[[j]] <- results
  plot_ID[j] <- unique(plotwise_bac[[j]]$Plot_ID)
}

# Set the Plot ID as a name for the list.
names(results_list_bac) <- plot_ID  

# Combine the results in a dataframe.
morans_i_bac_plots <- data.frame(do.call(rbind.data.frame, results_list_bac))
morans_i_bac_plots

###
#Community Barplots
###

# Aggregate the taxa at order level. This is to allow better visualization, 
# while still retaining the most information.
phy_ord_bac <- phy_bacteria %>%
  aggregate_taxa(level = "Order")

# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_bac_ord_top25 <- get_top_taxa(phy_ord_bac,
                                  24, relative = TRUE,
                                  other_label = "Others")

# Transform the subset dataset to compositional (relative) abundances.
phy_bac_ord_top25_plot <-  phy_bac_ord_top25 %>%
  aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
taxa_names(phy_bac_ord_top25_plot) <- tax_table(phy_bac_ord_top25_plot)[, 4]

# Turn the columns necessary for the ploting into sorted factors to set the order they are plotted in.
sample_data(phy_bac_ord_top25_plot)$intensity_formi_f <- factor(sample_data(phy_bac_ord_top25_plot)$intensity_formi, 
                                                                levels = c("low", "high")) 
sample_data(phy_bac_ord_top25_plot)$Plot_ID <- factor(sample_data(phy_bac_ord_top25_plot)$Plot_ID,
                                                      levels = c('HEW5', 'HEW11', 'HEW20', 'HEW27', 'HEW34','HEW35','HEW36','HEW37',
                                                                 'HEW8','HEW26','HEW28','HEW31','HEW32','HEW33','HEW43','HEW49')) 
sample_data(phy_bac_ord_top25_plot)$size_cat <- factor(sample_data(phy_bac_ord_top25_plot)$size_cat,
                                                       levels = c('big', 'medium', 'small')) 

sample_names(phy_bac_ord_top25_plot) <- c(paste0('T', 1:96))

# Sort the taxa names alphabetically.
taxa_names_bac <- sort(taxa_names(phy_bac_ord_top25_plot))

taxa_names_bac <- c("Abditibacteriales", "Acetobacterales", "Acidobacteriales", "Blastocatellales",
                    "Burkholderiales", "Caulobacterales", "Chitinophagales", "Chthoniobacterales",
                    "Cytophagales", "Deinococcales", "Fimbriimonadales", "Frankiales", "Isosphaerales",
                    "Micrococcales", "Micromonosporales", "Propionibacteriales", "Pseudomonadales",
                    "Pseudonocardiales", "Rhizobiales", "Solirubrobacterales", "Sphingobacteriales",
                    "Sphingomonadales", "Tepidisphaerales", "Thermomicrobiales", 'Others')

# Custom plotting to make a nice stacked barplot. 
bac_ord_plots <- phy_bac_ord_top25_plot %>%
  plot_composition( group_by =  'Plot_ID', otu.sort = taxa_names_bac, sample.sort = 'size_cat') +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Order') +
  guides(fill = guide_legend(title.position = 'top')) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_text(colour = "black", size = 5, hjust = 1),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black')) + 
  xlab('Tree ID') +
  ylab('Relative Abundance')  +
  scale_y_continuous(label = scales::percent)  + 
  labs( subtitle = '(B)')

bac_ord_plots

##---------
##  Fungi  
##----------------------------------------------------------------------------------------------

# Calculate the Shannon Diversity.
shannon_fun <- estimate_richness(phy_fungi, split = T, measures = 'Shannon')

# Make a new sample_data object and merge with existing phyloseq object.
new_var_fun <- sample_data(shannon_fun)
phy_fungi <- merge_phyloseq(phy_fungi, new_var_fun)

# Turn the sample_data object into a dataframe for ANOVA.
fun_eval_df <- data.frame(sample_data(phy_fungi))

# Conduct ANOVA and verify results with Tukey HSD.
anova_fun_alpha <- aov(Shannon ~ intensity_formi + size_cat, data = fun_eval_df)
summary(anova_fun_alpha)

TukeyHSD(anova_fun_alpha)

###
#Raincloud Plots of Alpha Diversity
###

# Code was adopted from Cédric Scherer. You can find a tutorial at: 
# https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/

# Compare Management Categories.

fungi_sampledata <- sample_data(phy_fungi)

fun_rain_manage <- ggplot(fungi_sampledata, aes(x = intensity_formi, y = Shannon, fill = intensity_formi, colour = intensity_formi)) + 
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .7, 
    ## move geom to the right
    justification = -.15, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) + 
  geom_boxplot(
    width = .14, 
    ## remove outliers
    outlier.color = NA, ## `outlier.shape = NA` works as well 
    alpha = 0.5,
    size = 0.1
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .5,
    size = .5
  )  +
  ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        axis.line = element_line(colour = "black", size = .2),
        axis.ticks = element_line(size = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9),
        legend.position = 'none', 
        plot.title = element_text(vjust = 0, hjust = 0.007),
        plot.subtitle = element_text(vjust = -0.5, hjust = 0.03, size = 9))+
  xlab('Management Intensity') + 
  labs(subtitle = '(C)') +
  scale_fill_manual(values = box_cols) + 
  scale_colour_manual(values = box_cols) + 
  scale_x_discrete(labels=c("low" = "low", "medium" = "high"))
fun_rain_manage

# Compare Tree Sizes.
fun_rain_size <- ggplot(fungi_sampledata, aes(x = size_cat, y = Shannon, fill = size_cat, colour = size_cat)) + 
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .5, 
    ## move geom to the right
    justification = -.25, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) + 
  geom_boxplot(
    width = .15, 
    ## remove outliers
    outlier.color = NA, ## `outlier.shape = NA` works as well 
    alpha = 0.5,
    size = 0.1
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .3, 
    ## add some transparency
    alpha = .5,
    size = .5
  )  +
  ## remove white space on the left
  coord_cartesian(xlim = c(1.2, NA)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.text = element_text(colour = "black", size = 9),
        axis.line = element_line(colour = "black", size = .2),
        axis.ticks = element_line(size = 0.2),
        legend.position = 'none', 
        plot.subtitle = element_text(vjust = -0.5, hjust = .03, size = 9)) +
  xlab('Tree Size') + 
  labs(subtitle = '(F)') +
  scale_fill_manual(values = box_cols) + 
  scale_colour_manual(values = box_cols) 
fun_rain_size

###
#Spatial Autocorrelation of alpha diversity
###

# Make a list that includes the sample data split by plots.
plotwise_fun <- sample_data(phy_fungi) %>% 
  dplyr::group_by(Plot_ID) %>%
  group_split()

# Initiate empty lists to fill.
results_list_fun <- list()
plot_ID <- NULL

# Calculate Morans I for the trees within each plot. 
for (j in 1:length(plotwise_fun)) {
  results <- moran_wrapper(plotwise_fun[[j]]$Shannon,
                           plotwise_fun[[j]]$x_coord,
                           plotwise_fun[[j]]$y_coord)
  results_list_fun[[j]] <- results
  plot_ID[j] <- unique(plotwise_fun[[j]]$Plot_ID)
}

# Set the Plot ID as a name for the list.
names(results_list_fun) <- plot_ID  

# Combine the results in a dataframe.
morans_i_fun_plots <- data.frame(do.call(rbind.data.frame, results_list_fun))
morans_i_fun_plots

###
#Community Barplots
###

# Aggregate the taxa at order level. This is to allow better visualization, 
# while still retaining the most information.
phy_order_fun <- phy_fungi %>%
  aggregate_taxa(level = "Order")

# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_fun_ord_top25 <- get_top_taxa(phy_order_fun, 24, relative = TRUE, other_label = "Others")

# Transform the subset dataset to compositional (relative) abundances.
phy_fun_ord_top25_plot <-  phy_fun_ord_top25 %>%
  aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
taxa_names(phy_fun_ord_top25_plot) <- tax_table(phy_fun_ord_top25_plot)[, 4]

# Turn the columns necessary for the ploting into sorted factors to set the order they are plotted in.
sample_data(phy_fun_ord_top25_plot)$Plot_ID <- factor(sample_data(phy_fun_ord_top25_plot)$Plot_ID,
                                                      levels = c('HEW5', 'HEW11', 'HEW20', 'HEW27', 'HEW34','HEW35','HEW36','HEW37',
                                                                 'HEW8','HEW26','HEW28','HEW31','HEW32','HEW33','HEW43','HEW49')) 
sample_data(phy_fun_ord_top25_plot)$size_cat <- factor(sample_data(phy_fun_ord_top25_plot)$size_cat,
                                                       levels = c('big', 'medium', 'small')) 

sample_names(phy_fun_ord_top25_plot) <- c(paste0('T', 1:96))

# Sort the taxa names alphabetically.
taxa_names_fun_ord <- sort(taxa_names(phy_fun_ord_top25_plot))

taxa_names_fun_ord <- c("Agaricales", "Arthoniales", "Caliciales", "Candelariales",                 
                        "Capnodiales", "Chaetothyriales", "Cystobasidiales", "Diaporthales",                         
                        "Dothideales", "Erythrobasidiales", "Exobasidiales", "Helotiales",                           
                        "Hypocreales", "Lecanorales", "Microbotryomycetes_ord_Incertae_sedis", "Myriangiales",                         
                        "Orbiliales", "Ostropales", "Phallales", "Pleosporales", "Polyporales",                          
                        "Sordariales", "Tremellales", "Unknown", "Others")

# Custom plotting to make a nice stacked barplot. 
fun_ord_plots <- phy_fun_ord_top25_plot %>%
  plot_composition( group_by =  'Plot_ID', otu.sort = taxa_names_fun_ord, sample.sort = 'size_cat') +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Order') +
  guides(fill = guide_legend(title.position = 'top')) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_text(colour = "black", size = 5, hjust = 1),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black')) + 
  xlab('Tree ID') +
  ylab('Relative Abundance')  +
  scale_y_continuous(label = scales::percent) + 
  labs( subtitle = '(C)')

fun_ord_plots

#################################################################
##                          Section 4                          ##
##                   Intra-Group Interactions                  ##
#################################################################

##---------
##  Algae  
##-------------------------------------------------------------------------------------------------------

# Subset the raw dataset to ASVs that contribute at least one percent to all reads of a sample.
phy_algae_trans  = transform_sample_counts(phy_algae, function(x) x / sum(x) )
phy_algae_abundfilt = filter_taxa(phy_algae_trans, function(x) sum(x) > .01, TRUE)
keep_names_alg <- taxa_names(phy_algae_abundfilt)

algae_subset <- subset(otu_table(phy_algae), rownames(otu_table(phy_algae)) %in% keep_names_alg)
phy_algae_raw_abundfilt <- merge_phyloseq(algae_subset, tax_table(phy_algae), sample_data(phy_algae))

# Save the phyloseq object to run the network inference on the server.
saveRDS(phy_algae_raw_abundfilt, 'phy_algae_raw_abundfilt.rds')

# The following code is outcommented, because it is very computationally intense. Hence, we recommend to not run this on a personal computer.
# Rather run it on a server or massive ram machine. 
#
# Load the necessary package on the servers R installation.
# library(SpiecEasi)
# 
# Import the saved phyloseq object.
# phy_algae_raw_abundfilt <- readRDS('phy_algae_raw_abundfilt.rds')
# 
# Set the arguments for the Spiec-Easi algorithm, Note that it runs on 20 cores simultaneously. See the tutorial at https://github.com/zdk123/SpiecEasi.
# pargs2 <- list(rep.num=50, seed=10010, ncores=20)
# 
# Run the main Spiec-Easi algorithm.
# se_alg_raw_abundfilt <- spiec.easi(phy_algae_raw_abundfilt, method='mb',
#                                    lambda.min.ratio=1e-1, nlambda=100,
#                                    sel.criterion='bstars', pulsar.select=TRUE,
#                                    pulsar.params=pargs2)
# 
# Save the network and continue working on the local computer.
# saveRDS(se_alg_raw_abundfilt, 'se_alg_raw_abundfilt.rds')

# Import the network and turn it into an igraph object for further processing.
se_alg_raw_abundfilt <- readRDS(here('se_alg_raw_abundfilt.rds'))
alg.mb <- adj2igraph(getRefit(se_alg_raw_abundfilt),  
                     vertex.attr=list(name=taxa_names(phy_algae_raw_abundfilt)))

# Obtain the networks stability to judge good fit. Close to 0.05 is good. 
getStability(se_alg_raw_abundfilt)

###
# Plot the network.
###

# First we need to convert it to a gephi readable format. 
alg.mb_gephi <- igraph.to.gexf(alg.mb)

write.gexf(alg.mb_gephi, output = here('algae.gexf'))

# Export the network and read it into Gephi. Then obtain the modules and do the plotting.

# Get each taxas, betweenness centrality values.
alg_between <- betweenness(alg.mb, directed = F)

# Subset it to the top five taxa to obtain hub taxa.
alg_top5_between <- sort(alg_between, decreasing = T)[1:5]

# Subset the phyloseq object, to only contain the hub taxa to ease displaying their taxonomy.
hub_taxa_alg <- subset_taxa(phy_algae, taxa_names(phy_algae) %in% names(alg_top5_between))
tax_table(hub_taxa_alg)

####
#Module composition
####

# Read in the data on the modules obtained from Gephi.
modules_alg <- read.csv(here('modules_algae.csv'), sep = ",")

# Transform ASV information to dataframe to make handling easier.
otu_tab_alg <- as.data.frame(otu_table(phy_algae_raw_abundfilt))

# Calculate the total read count of each ASV.
otu_tab_alg$total <- rowSums(otu_tab_alg)

# Keep only the total read count.
otu_tab_alg <- otu_tab_alg %>% 
  select(.,total) %>%
  rownames_to_column('Label')

# Merge module table and read count information.
otu_modules_alg <- inner_join(otu_tab_alg, modules_alg)

# Initiate an empty dataframe were we can store the total read counts per module.
otu_mod_alg <- data.frame(Label = otu_modules_alg$Label,
                              Mod0 = rep(0, 129),
                              Mod1 = rep(0, 129),
                              Mod2 = rep(0, 129),
                              Mod3 = rep(0, 129),
                              Mod4 = rep(0, 129),
                              Mod5 = rep(0, 129),
                              Mod6 = rep(0, 129),
                              Mod7 = rep(0, 129),
                              Mod8 = rep(0, 129),
                              Mod9 = rep(0, 129),
                              Mod11 = rep(0, 129))

# If the ASV belongs to a module then write its read count, if not write zero. 
otu_mod_alg$Mod0 <- if_else(otu_modules_alg$modularity_class == 0, otu_modules_alg$total, 0)
otu_mod_alg$Mod1 <- if_else(otu_modules_alg$modularity_class == 1, otu_modules_alg$total, 0)
otu_mod_alg$Mod2 <- if_else(otu_modules_alg$modularity_class == 2, otu_modules_alg$total, 0)
otu_mod_alg$Mod3 <- if_else(otu_modules_alg$modularity_class == 3, otu_modules_alg$total, 0)
otu_mod_alg$Mod4 <- if_else(otu_modules_alg$modularity_class == 4, otu_modules_alg$total, 0)
otu_mod_alg$Mod5 <- if_else(otu_modules_alg$modularity_class == 5, otu_modules_alg$total, 0)
otu_mod_alg$Mod6 <- if_else(otu_modules_alg$modularity_class == 6, otu_modules_alg$total, 0)
otu_mod_alg$Mod7 <- if_else(otu_modules_alg$modularity_class == 7, otu_modules_alg$total, 0)
otu_mod_alg$Mod8 <- if_else(otu_modules_alg$modularity_class == 8, otu_modules_alg$total, 0)
otu_mod_alg$Mod9 <- if_else(otu_modules_alg$modularity_class == 9, otu_modules_alg$total, 0)
otu_mod_alg$Mod11 <- if_else(otu_modules_alg$modularity_class == 11, otu_modules_alg$total, 0)

# set the ASV name back to the rownames.
rownames(otu_mod_alg) <- otu_mod_alg$Label
otu_mod_alg$Label <- NULL

# Create the phyloseq object to easily obtain the taxonomic information.
otu_mat_mod_alg <- as.matrix(otu_mod_alg)
tax_mat_mod_alg <- as.matrix(tax_table(phy_algae_raw_abundfilt))

OTU_mod_alg <- otu_table(otu_mat_mod_alg, taxa_are_rows = T)
TAX_mod_alg <- tax_table(tax_mat_mod_alg)

phy_modules_alg <- phyloseq(OTU_mod_alg, TAX_mod_alg)

# Find out which modules contain more than 10 ASVs.
module_asv_count_alg <- colSums(otu_mat_mod_alg != 0)
which(module_asv_count_alg > 10)
sort(module_asv_count_alg)

####
# Find out which taxa are the top ones in the modules with more than ten ASVs.
####

# Subset by Module name.
phy_algae_mod2 <- prune_samples(sample_names(phy_modules_alg) == 'Mod2', phy_modules_alg)
phy_algae_mod4 <- prune_samples(sample_names(phy_modules_alg) == 'Mod4', phy_modules_alg)
phy_algae_mod5 <- prune_samples(sample_names(phy_modules_alg) == 'Mod5', phy_modules_alg)
phy_algae_mod7 <- prune_samples(sample_names(phy_modules_alg) == 'Mod7', phy_modules_alg)
phy_algae_mod9 <- prune_samples(sample_names(phy_modules_alg) == 'Mod9', phy_modules_alg)

# Remove any taxa without reads.
phy_algae_mod2_clean <- prune_taxa(taxa_sums(phy_algae_mod2) > 0, phy_algae_mod2)
phy_algae_mod4_clean <- prune_taxa(taxa_sums(phy_algae_mod4) > 0, phy_algae_mod4)
phy_algae_mod5_clean <- prune_taxa(taxa_sums(phy_algae_mod5) > 0, phy_algae_mod5)
phy_algae_mod7_clean <- prune_taxa(taxa_sums(phy_algae_mod7) > 0, phy_algae_mod7)
phy_algae_mod9_clean <- prune_taxa(taxa_sums(phy_algae_mod9) > 0, phy_algae_mod9)


# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 2. 
top_algae_mod2 <- get_top_taxa(phy_algae_mod2_clean, 1, discard_other = T)
tax_table(top_algae_mod2)
abundances(top_algae_mod2) / sum(abundances(phy_algae_mod2_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 4. 
top_algae_mod4 <- get_top_taxa(phy_algae_mod4_clean, 1, discard_other = T)
tax_table(top_algae_mod4)
abundances(top_algae_mod4) / sum(abundances(phy_algae_mod4_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 5. 
top_algae_mod5 <- get_top_taxa(phy_algae_mod5_clean, 1, discard_other = T)
tax_table(top_algae_mod5)
abundances(top_algae_mod5) / sum(abundances(phy_algae_mod5_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 7. 
top_algae_mod7 <- get_top_taxa(phy_algae_mod7_clean, 1, discard_other = T)
tax_table(top_algae_mod7)
abundances(top_algae_mod7) / sum(abundances(phy_algae_mod7_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 9. 
top_algae_mod9 <- get_top_taxa(phy_algae_mod9_clean, 1, discard_other = T)
tax_table(top_algae_mod9)
abundances(top_algae_mod9) / sum(abundances(phy_algae_mod9_clean))

##---------
##  Bacteria  
##-------------------------------------------------------------------------------------------------------

# Subset the raw dataset to ASVs that contribute at least one percent to all reads of a sample.
phy_bacteria_trans  = transform_sample_counts(phy_bacteria, function(x) x / sum(x) )
phy_bacteria_abundfilt = filter_taxa(phy_bacteria_trans, function(x) sum(x) > .01, TRUE)
keep_names <- taxa_names(phy_bacteria_abundfilt)

my_subset <- subset(otu_table(phy_bacteria), rownames(otu_table(phy_bacteria)) %in% keep_names)
phy_bacteria_raw_abundfilt <- merge_phyloseq(my_subset, tax_table(phy_bacteria), sample_data(phy_bacteria))

# Save the phyloseq object to run the network inference on the server.
saveRDS(phy_bacteria_raw_abundfilt, 'phy_bacteria_raw_abundfilt.rds')

# The following code is outcommented, because it is very computationally intense. Hence, we recommend to not run this on a personal computer.
# Rather run it on a server or massive ram machine. 
#
# Load the necessary package on the servers R installation.
# library(SpiecEasi)
# 
# Import the saved phyloseq object.
# phy_bacteria_raw_abundfilt <- readRDS('phy_bacteria_raw_abundfilt.rds')
# 
# Set the arguments for the Spiec-Easi algorithm, Note that it runs on 20 cores simultaneously. See the tutorial at https://github.com/zdk123/SpiecEasi.
# pargs2 <- list(rep.num=50, seed=10010, ncores=20)
# 
# Run the main Spiec-Easi algorithm.
# se_bac_raw_abundfilt <- spiec.easi(phy_bacteria_raw_abundfilt, method='mb',
#                                    lambda.min.ratio=1e-1, nlambda=100,
#                                    sel.criterion='bstars', pulsar.select=TRUE,
#                                    pulsar.params=pargs2)
# 
# Save the network and continue working on the local computer.
# saveRDS(se_bac_raw_abundfilt, 'se_bac_raw_abundfilt.rds')

# Import the network and turn it into an igraph object for further processing.
se_bac_raw_abundfilt <- readRDS(here('se_bac_raw_abundfilt.rds'))
bac.mb <- adj2igraph(getRefit(se_bac_raw_abundfilt),  
                     vertex.attr=list(name=taxa_names(phy_bacteria_raw_abundfilt)))

# Obtain the networks stability to judge good fit. Close to 0.05 is good. 
getStability(se_bac_raw_abundfilt)

###
#plot the network
###

# First we need to convert it to a gephi readable format. 

bac.mb_gephi <- igraph.to.gexf(bac.mb)

write.gexf(bac.mb_gephi, output = here('bacteria.gexf'))

# Export the network and read it into Gephi. Then obtain the modules and do the plotting.

# Get each taxas, betweenness centrality values.
bac_between <- betweenness(bac.mb, directed = F)

# Subset it to the top five taxa to obtain hub taxa.
bac_top5_between <- sort(bac_between, decreasing = T)[1:5]

# Subset the phyloseq object, to only contain the hub taxa to ease displaying their taxonomy.
hub_taxa_bac <- subset_taxa(phy_bacteria, taxa_names(phy_bacteria) %in% names(bac_top5_between))
tax_table(hub_taxa_bac)

####
#Module composition
####

# Read in the data on the modules obtained from Gephi.
modules_bac <- read.csv(here('modules_bacteria.csv'))

# Transform ASV information to dataframe to make handling easier.
otu_tab_bac <- as.data.frame(otu_table(phy_bacteria_raw_abundfilt))

# Calculate the total read count of each ASV.
otu_tab_bac$total <- rowSums(otu_tab_bac)

# Keep only the total read count.
otu_tab_bac <- otu_tab_bac %>% 
  select(.,total) %>%
  rownames_to_column('Label')

# Merge module table and read count information.
otu_modules_bac <- inner_join(otu_tab_bac, modules_bac)

# Initiate an empty dataframe were we can store the total read counts per module.
otu_mod_bac <- data.frame(Label = otu_modules_bac$Label,
                          Mod0 = rep(0, 628),
                          Mod1 = rep(0, 628),
                          Mod2 = rep(0, 628),
                          Mod3 = rep(0, 628),
                          Mod4 = rep(0, 628),
                          Mod5 = rep(0, 628),
                          Mod6 = rep(0, 628))

# If the ASV belongs to a module then write its read count, if not write zero. 
otu_mod_bac$Mod0 <- if_else(otu_modules_bac$modularity_class == 0, otu_modules_bac$total, 0)
otu_mod_bac$Mod1 <- if_else(otu_modules_bac$modularity_class == 1, otu_modules_bac$total, 0)
otu_mod_bac$Mod2 <- if_else(otu_modules_bac$modularity_class == 2, otu_modules_bac$total, 0)
otu_mod_bac$Mod3 <- if_else(otu_modules_bac$modularity_class == 3, otu_modules_bac$total, 0)
otu_mod_bac$Mod4 <- if_else(otu_modules_bac$modularity_class == 4, otu_modules_bac$total, 0)
otu_mod_bac$Mod5 <- if_else(otu_modules_bac$modularity_class == 5, otu_modules_bac$total, 0)
otu_mod_bac$Mod6 <- if_else(otu_modules_bac$modularity_class == 6, otu_modules_bac$total, 0)

# set the ASV name back to the rownames.
rownames(otu_mod_bac) <- otu_mod_bac$Label
otu_mod_bac$Label <- NULL

# Create the phyloseq object to easily obtain the taxonomic information.
otu_mat_mod_bac <- as.matrix(otu_mod_bac)
tax_mat_mod_bac <- as.matrix(tax_table(phy_bacteria_raw_abundfilt))

OTU_mod_bac <- otu_table(otu_mat_mod_bac, taxa_are_rows = T)
TAX_mod_bac <- tax_table(tax_mat_mod_bac)

phy_modules_bac <- phyloseq(OTU_mod_bac, TAX_mod_bac)

# Find out which modules contain more than 10 ASVs.
module_asv_count_bac <- colSums(otu_mat_mod_bac != 0)
which(module_asv_count_bac > 10)
sort(module_asv_count_bac)

####
# Find out which taxa are the top ones in the modules with more than ten ASVs.
####

# Subset by Module name.
phy_bacteria_mod0 <- prune_samples(sample_names(phy_modules_bac) == 'Mod0', phy_modules_bac)
phy_bacteria_mod1 <- prune_samples(sample_names(phy_modules_bac) == 'Mod1', phy_modules_bac)
phy_bacteria_mod2 <- prune_samples(sample_names(phy_modules_bac) == 'Mod2', phy_modules_bac)
phy_bacteria_mod3 <- prune_samples(sample_names(phy_modules_bac) == 'Mod3', phy_modules_bac)
phy_bacteria_mod4 <- prune_samples(sample_names(phy_modules_bac) == 'Mod4', phy_modules_bac)
phy_bacteria_mod5 <- prune_samples(sample_names(phy_modules_bac) == 'Mod5', phy_modules_bac)
phy_bacteria_mod6 <- prune_samples(sample_names(phy_modules_bac) == 'Mod6', phy_modules_bac)

# Remove any taxa without reads.
phy_bacteria_mod0_clean <- prune_taxa(taxa_sums(phy_bacteria_mod0) > 0, phy_bacteria_mod0)
phy_bacteria_mod1_clean <- prune_taxa(taxa_sums(phy_bacteria_mod1) > 0, phy_bacteria_mod1)
phy_bacteria_mod2_clean <- prune_taxa(taxa_sums(phy_bacteria_mod2) > 0, phy_bacteria_mod2)
phy_bacteria_mod3_clean <- prune_taxa(taxa_sums(phy_bacteria_mod3) > 0, phy_bacteria_mod3)
phy_bacteria_mod4_clean <- prune_taxa(taxa_sums(phy_bacteria_mod4) > 0, phy_bacteria_mod4)
phy_bacteria_mod5_clean <- prune_taxa(taxa_sums(phy_bacteria_mod5) > 0, phy_bacteria_mod5)
phy_bacteria_mod6_clean <- prune_taxa(taxa_sums(phy_bacteria_mod6) > 0, phy_bacteria_mod6)

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 0. 
top_bacteria_mod0 <- get_top_taxa(phy_bacteria_mod0_clean, 1, discard_other = T)
tax_table(top_bacteria_mod0)
abundances(top_bacteria_mod0) / sum(abundances(phy_bacteria_mod0_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 1. 
top_bacteria_mod1 <- get_top_taxa(phy_bacteria_mod1_clean, 1, discard_other = T)
tax_table(top_bacteria_mod1)
abundances(top_bacteria_mod1) / sum(abundances(phy_bacteria_mod1_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 2. 
top_bacteria_mod2 <- get_top_taxa(phy_bacteria_mod2_clean, 1, discard_other = T)
tax_table(top_bacteria_mod2)
abundances(top_bacteria_mod2) / sum(abundances(phy_bacteria_mod2_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 3. 
top_bacteria_mod3 <- get_top_taxa(phy_bacteria_mod3_clean, 1, discard_other = T)
tax_table(top_bacteria_mod3)
abundances(top_bacteria_mod3) / sum(abundances(phy_bacteria_mod3_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 4. 
top_bacteria_mod4 <- get_top_taxa(phy_bacteria_mod4_clean, 1, discard_other = T)
tax_table(top_bacteria_mod4)
abundances(top_bacteria_mod4) / sum(abundances(phy_bacteria_mod4_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 5. 
top_bacteria_mod5 <- get_top_taxa(phy_bacteria_mod5_clean, 1, discard_other = T)
tax_table(top_bacteria_mod5)
abundances(top_bacteria_mod5) / sum(abundances(phy_bacteria_mod5_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 6. 
top_bacteria_mod6 <- get_top_taxa(phy_bacteria_mod6_clean, 1, discard_other = T)
tax_table(top_bacteria_mod6)
abundances(top_bacteria_mod6) / sum(abundances(phy_bacteria_mod6_clean))



##---------
##  Fungi  
##-------------------------------------------------------------------------------------------------------

# Subset the raw dataset to ASVs that contribute at least one percent to all reads of a sample.
phy_fungi_trans  = transform_sample_counts(phy_fungi, function(x) x / sum(x) )
phy_fungi_abundfilt = filter_taxa(phy_fungi_trans, function(x) sum(x) > .01, TRUE)
keep_names <- taxa_names(phy_fungi_abundfilt)

my_subset <- subset(otu_table(phy_fungi), rownames(otu_table(phy_fungi)) %in% keep_names)
phy_fungi_raw_abundfilt <- merge_phyloseq(my_subset, tax_table(phy_fungi), sample_data(phy_fungi))

# Save the phyloseq object to run the network inference on the server.
saveRDS(phy_fungi_raw_abundfilt, 'phy_fungi_raw_abundfilt.rds')

# The following code is outcommented, because it is very computationally intense. Hence, we recommend to not run this on a personal computer.
# Rather run it on a server or massive ram machine. 
#
# Load the necessary package on the servers R installation.
# library(SpiecEasi)
# 
# Import the saved phyloseq object.
# phy_fungi_raw_abundfilt <- readRDS('phy_fungi_raw_abundfilt.rds')
# 
# Set the arguments for the Spiec-Easi algorithm, Note that it runs on 20 cores simultaneously. See the tutorial at https://github.com/zdk123/SpiecEasi.
# pargs2 <- list(rep.num=50, seed=10010, ncores=20)
# 
# Run the main Spiec-Easi algorithm.
# se_fun_raw_abundfilt <- spiec.easi(phy_fungi_raw_abundfilt, method='mb',
#                                    lambda.min.ratio=1e-1, nlambda=100,
#                                    sel.criterion='bstars', pulsar.select=TRUE,
#                                    pulsar.params=pargs2)
# 
# Save the network and continue working on the local computer.
# saveRDS(se_fun_raw_abundfilt, 'se_fun_raw_abundfilt.rds')

# Import the network and turn it into an igraph object for further processing.
se_fun_raw_abundfilt <- readRDS(here('se_fun_raw_abundfilt.rds'))
fun.mb <- adj2igraph(getRefit(se_fun_raw_abundfilt),  
                     vertex.attr=list(name=taxa_names(phy_fungi_raw_abundfilt)))

# Obtain the networks stability to judge good fit. Close to 0.05 is good. 
getStability(se_fun_raw_abundfilt)

###
#plot the network
###

# First we need to convert it to a gephi readable format. 
fun.mb_gephi <- igraph.to.gexf(fun.mb)

write.gexf(fun.mb_gephi, output = here('fungi.gexf'))


# Export the network and read it into Gephi. Then obtain the modules and do the plotting.

# Get each taxas, betweenness centrality values.
fun_between <- betweenness(fun.mb, directed = F)

# Subset it to the top five taxa to obtain hub taxa.
fun_top5_between <- sort(fun_between, decreasing = T)[1:5]

# Subset the phyloseq object, to only contain the hub taxa to ease displaying their taxonomy.
hub_taxa_fun <- subset_taxa(phy_fungi, taxa_names(phy_fungi) %in% names(fun_top5_between))
tax_table(hub_taxa_fun)

####
#Module composition
####

# Read in the data on the modules obtained from Gephi.
modules_fungi <- read.csv(here('modules_fungi.csv'))

# Transform ASV information to dataframe to make handling easier.
otu_tab_fungi <- as.data.frame(otu_table(phy_fungi_raw_abundfilt))

# Calculate the total read count of each ASV.
otu_tab_fungi$total <- rowSums(otu_tab_fungi)

# Keep only the total read count.
otu_tab_fungi <- otu_tab_fungi %>% 
  select(.,total) %>%
  rownames_to_column('Label')

# Merge module table and read count information.
otu_modules_fungi <- inner_join(otu_tab_fungi, modules_fungi)

# Initiate an empty dataframe were we can store the total read counts per module.
otu_mod_fungi <- data.frame(Label = otu_modules_fungi$Label,
                          Mod0 = rep(0, 289),
                          Mod1 = rep(0, 289),
                          Mod2 = rep(0, 289),
                          Mod3 = rep(0, 289),
                          Mod4 = rep(0, 289),
                          Mod5 = rep(0, 289),
                          Mod6 = rep(0, 289),
                          Mod7 = rep(0, 289),
                          Mod8 = rep(0, 289))

# If the ASV belongs to a module then write its read count, if not write zero. 
otu_mod_fungi$Mod0 <- if_else(otu_modules_fungi$modularity_class == 0, otu_modules_fungi$total, 0)
otu_mod_fungi$Mod1 <- if_else(otu_modules_fungi$modularity_class == 1, otu_modules_fungi$total, 0)
otu_mod_fungi$Mod2 <- if_else(otu_modules_fungi$modularity_class == 2, otu_modules_fungi$total, 0)
otu_mod_fungi$Mod3 <- if_else(otu_modules_fungi$modularity_class == 3, otu_modules_fungi$total, 0)
otu_mod_fungi$Mod4 <- if_else(otu_modules_fungi$modularity_class == 4, otu_modules_fungi$total, 0)
otu_mod_fungi$Mod5 <- if_else(otu_modules_fungi$modularity_class == 5, otu_modules_fungi$total, 0)
otu_mod_fungi$Mod6 <- if_else(otu_modules_fungi$modularity_class == 6, otu_modules_fungi$total, 0)
otu_mod_fungi$Mod7 <- if_else(otu_modules_fungi$modularity_class == 7, otu_modules_fungi$total, 0)
otu_mod_fungi$Mod8 <- if_else(otu_modules_fungi$modularity_class == 8, otu_modules_fungi$total, 0)

# set the ASV name back to the rownames.
rownames(otu_mod_fungi) <- otu_mod_fungi$Label
otu_mod_fungi$Label <- NULL

# Create the phyloseq object to easily obtain the taxonomic information.
otu_mat_mod_fungi <- as.matrix(otu_mod_fungi)
tax_mat_mod_fungi <- as.matrix(tax_table(phy_fungi_raw_abundfilt))

OTU_mod_fun <- otu_table(otu_mat_mod_fungi, taxa_are_rows = T)
TAX_mod_fun <- tax_table(tax_mat_mod_fungi)

phy_modules_fun <- phyloseq(OTU_mod_fun, TAX_mod_fun)

# Find out which modules contain more than 10 ASVs.
module_asv_count_fun <- colSums(otu_mat_mod_fungi != 0)
which(module_asv_count_fun > 10)
sort(module_asv_count_fun)

####
# Find out which taxa are the top ones in the modules with more than ten ASVs.
####

# Subset by Module name.
phy_fungi_mod0 <- prune_samples(sample_names(phy_modules_fun) == 'Mod0', phy_modules_fun)
phy_fungi_mod1 <- prune_samples(sample_names(phy_modules_fun) == 'Mod1', phy_modules_fun)
phy_fungi_mod2 <- prune_samples(sample_names(phy_modules_fun) == 'Mod2', phy_modules_fun)
phy_fungi_mod3 <- prune_samples(sample_names(phy_modules_fun) == 'Mod3', phy_modules_fun)
phy_fungi_mod4 <- prune_samples(sample_names(phy_modules_fun) == 'Mod4', phy_modules_fun)
phy_fungi_mod5 <- prune_samples(sample_names(phy_modules_fun) == 'Mod5', phy_modules_fun)
phy_fungi_mod6 <- prune_samples(sample_names(phy_modules_fun) == 'Mod6', phy_modules_fun)
phy_fungi_mod7 <- prune_samples(sample_names(phy_modules_fun) == 'Mod7', phy_modules_fun)
phy_fungi_mod8 <- prune_samples(sample_names(phy_modules_fun) == 'Mod8', phy_modules_fun)

# Remove any taxa without reads.
phy_fungi_mod0_clean <- prune_taxa(taxa_sums(phy_fungi_mod0) > 0, phy_fungi_mod0)
phy_fungi_mod1_clean <- prune_taxa(taxa_sums(phy_fungi_mod1) > 0, phy_fungi_mod1)
phy_fungi_mod2_clean <- prune_taxa(taxa_sums(phy_fungi_mod2) > 0, phy_fungi_mod2)
phy_fungi_mod3_clean <- prune_taxa(taxa_sums(phy_fungi_mod3) > 0, phy_fungi_mod3)
phy_fungi_mod4_clean <- prune_taxa(taxa_sums(phy_fungi_mod4) > 0, phy_fungi_mod4)
phy_fungi_mod5_clean <- prune_taxa(taxa_sums(phy_fungi_mod5) > 0, phy_fungi_mod5)
phy_fungi_mod6_clean <- prune_taxa(taxa_sums(phy_fungi_mod6) > 0, phy_fungi_mod6)
phy_fungi_mod7_clean <- prune_taxa(taxa_sums(phy_fungi_mod7) > 0, phy_fungi_mod7)
phy_fungi_mod8_clean <- prune_taxa(taxa_sums(phy_fungi_mod8) > 0, phy_fungi_mod8)

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 0. 
top_fungi_mod0 <- get_top_taxa(phy_fungi_mod0_clean, 1, discard_other = T)
tax_table(top_fungi_mod0)
abundances(top_fungi_mod0) / sum(abundances(phy_fungi_mod0_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 1. 
top_fungi_mod1 <- get_top_taxa(phy_fungi_mod1_clean, 1, discard_other = T)
tax_table(top_fungi_mod1)
abundances(top_fungi_mod1) / sum(abundances(phy_fungi_mod1_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 2. 
top_fungi_mod2 <- get_top_taxa(phy_fungi_mod2_clean, 1, discard_other = T)
tax_table(top_fungi_mod2)
abundances(top_fungi_mod2) / sum(abundances(phy_fungi_mod2_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 3. 
top_fungi_mod3 <- get_top_taxa(phy_fungi_mod3_clean, 1, discard_other = T)
tax_table(top_fungi_mod3)
abundances(top_fungi_mod3) / sum(abundances(phy_fungi_mod3_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 4. 
top_fungi_mod4 <- get_top_taxa(phy_fungi_mod4_clean, 1, discard_other = T)
tax_table(top_fungi_mod4)
abundances(top_fungi_mod4) / sum(abundances(phy_fungi_mod4_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 5. 
top_fungi_mod5 <- get_top_taxa(phy_fungi_mod5_clean, 1, discard_other = T)
tax_table(top_fungi_mod5)
abundances(top_fungi_mod5) / sum(abundances(phy_fungi_mod5_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 6. 
top_fungi_mod6 <- get_top_taxa(phy_fungi_mod6_clean, 1, discard_other = T)
tax_table(top_fungi_mod6)
abundances(top_fungi_mod6) / sum(abundances(phy_fungi_mod6_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 7. 
top_fungi_mod7 <- get_top_taxa(phy_fungi_mod7_clean, 1, discard_other = T)
tax_table(top_fungi_mod7)
abundances(top_fungi_mod7) / sum(abundances(phy_fungi_mod7_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 8. 
top_fungi_mod8 <- get_top_taxa(phy_fungi_mod8_clean, 1, discard_other = T)
tax_table(top_fungi_mod8)
abundances(top_fungi_mod8) / sum(abundances(phy_fungi_mod8_clean))

#################################################################
##                          Section 5                          ##
##                   Inter-Group Interactions                  ##
#################################################################

# Append the ending _B to the bacterial ASV names to indicate they are bacteria.
phy_bac_full <- phy_bacteria_raw_abundfilt
taxa_names(phy_bac_full) <- paste0(taxa_names(phy_bac_full),'_B')

# Append the ending _F to the fungal ASV names to indicate they are bacteria.
phy_fun_full <- phy_fungi_raw_abundfilt
taxa_names(phy_fun_full) <- paste0(taxa_names(phy_fun_full),'_F')

# Append the ending _A to the algal ASV names to indicate they are bacteria.
phy_alg_full <- phy_algae_raw_abundfilt
taxa_names(phy_alg_full) <- paste0(taxa_names(phy_alg_full),'_A')


# Merge the organismal groups into one phyloseq object.
full_physeq <- merge_phyloseq(phy_fun_full, 
                              phy_bac_full,
                              phy_alg_full)

# Save the phyloseq object to run the network inference on the server.
saveRDS(full_physeq, 'full_physeq.rds')

# The following code is outcommented, because it is very computationally intense. Hence, we recommend to not run this on a personal computer.
# Rather run it on a server or massive ram machine. 
#
# Load the necessary package on the servers R installation.
# library(spieceasi)
#
# Import the saved phyloseq object.
# full_physeq <- readRDS('full_physeq.rds')
# 
# Set the arguments for the Spiec-Easi algorithm, Note that it runs on 20 cores simultaneously. See the tutorial at https://github.com/zdk123/SpiecEasi.
# pargs2 <- list(rep.num=50, seed=10010, ncores=20)
# 
# Run the main Spiec-Easi algorithm.
# se_cross_raw_abundfilt <- spiec.easi(full_physeq, method='mb',
#                                    lambda.min.ratio=1e-1, nlambda=100,
#                                    sel.criterion='bstars', pulsar.select=TRUE,
#                                    pulsar.params=pargs2)
# 
# Save the network and continue working on the local computer.
# saveRDS(se_cross_raw_abundfilt, 'se_cross_raw_abundfilt.rds')

# Import the network and turn it into an igraph object for further processing.
se_cross_raw_abundfilt <- readRDS(here('se_cross_raw_abundfilt.rds'))
cross.mb <- adj2igraph(getRefit(se_cross_raw_abundfilt),  
                       vertex.attr=list(name=taxa_names(full_physeq)))

# Obtain the networks stability to judge good fit. Close to 0.05 is good. 
getStability(se_cross_raw_abundfilt)

###
#plot the network
###

# First we need to convert it to a gephi readable format. 
cross.mb_gephi <- igraph.to.gexf(cross.mb)

write.gexf(cross.mb_gephi, output = here('cross.gexf'))

# Export the network and read it into Gephi. Then obtain the modules and do the plotting.

# Get each taxas betweenness centrality values.
cross_between <- betweenness(cross.mb, directed = F)

# Subset it to the top five taxa to obtain hub taxa.
cross_top5_between <- sort(cross_between, decreasing = T)[1:5]

# Subset the phyloseq object, to only contain the hub taxa to ease displaying their taxonomy.
hub_taxa_cross <- subset_taxa(full_physeq, taxa_names(full_physeq) %in% names(cross_top5_between))
tax_table(hub_taxa_cross)

####
#Module composition
####

# Read in the data on the modules obtained from Gephi.
modules <- read.csv(here('modules_cross.csv'))

# To obtain the modules relative abundances we need all the ASV counts.
# Append the endings _B, _F, _A to the raw reads to match the ASV names to the names used for the network.
phy_bac_rel <- phy_bacteria
taxa_names(phy_bac_rel) <- paste0(taxa_names(phy_bac_rel),'_B')

phy_fun_rel <- phy_fungi
taxa_names(phy_fun_rel) <- paste0(taxa_names(phy_fun_rel),'_F')

phy_alg_rel <- phy_algae
taxa_names(phy_alg_rel) <- paste0(taxa_names(phy_alg_rel),'_A')

##---------
##  Algae  
##-------------------------------------------------------------------------------------------------------

# Transform to dataframe to make handling easier.
otu_abundance_alg <- as.data.frame(otu_table(phy_alg_rel))

# Calculate the total read count of each ASV.
otu_abundance_alg$total <- rowSums(otu_abundance_alg)

# Keep only the total read count.
otu_abundance_alg <- otu_abundance_alg %>% 
  select(.,total) %>%
  rownames_to_column('Label')

# Merge module table and read count information.
otu_modules_alg <- inner_join(otu_abundance_alg, modules)

# Initiate an empty dataframe were we can store the total read counts per module.
otu_tab_mod_alg <- data.frame(Label = otu_modules_alg$Label,
                              Mod0 = rep(0,129),
                              Mod1 = rep(0,129),
                              Mod2 = rep(0,129),
                              Mod3 = rep(0,129),
                              Mod4 = rep(0,129),
                              Mod5 = rep(0,129),
                              Mod6 = rep(0,129),
                              Mod7 = rep(0,129))

# If the ASV belongs to a module then write its read count, if not write zero. 
otu_tab_mod_alg$Mod0 <- if_else(otu_modules_alg$modularity_class == 0, otu_modules_alg$total, 0)
otu_tab_mod_alg$Mod1 <- if_else(otu_modules_alg$modularity_class == 1, otu_modules_alg$total, 0)
otu_tab_mod_alg$Mod2 <- if_else(otu_modules_alg$modularity_class == 2, otu_modules_alg$total, 0)
otu_tab_mod_alg$Mod3 <- if_else(otu_modules_alg$modularity_class == 3, otu_modules_alg$total, 0)
otu_tab_mod_alg$Mod4 <- if_else(otu_modules_alg$modularity_class == 4, otu_modules_alg$total, 0)
otu_tab_mod_alg$Mod5 <- if_else(otu_modules_alg$modularity_class == 5, otu_modules_alg$total, 0)
otu_tab_mod_alg$Mod6 <- if_else(otu_modules_alg$modularity_class == 6, otu_modules_alg$total, 0)
otu_tab_mod_alg$Mod7 <- if_else(otu_modules_alg$modularity_class == 7, otu_modules_alg$total, 0)

# set the ASV name back to the rownames.
rownames(otu_tab_mod_alg) <- otu_tab_mod_alg$Label
otu_tab_mod_alg$Label <- NULL

# Sum up the columns to obtain module total read counts.
algae_mod_rel <- as.data.frame(colSums(otu_tab_mod_alg))
algae_mod_rel

##---------
##  Bacteria  
##-------------------------------------------------------------------------------------------------------

# Transform to dataframe to make handling easier.
otu_abundance_bac <- as.data.frame(otu_table(phy_bac_rel))

# Calculate the total read count of each ASV.
otu_abundance_bac$total <- rowSums(otu_abundance_bac)

# Keep only the total read count.
otu_abundance_bac <- otu_abundance_bac %>% 
  select(total) %>%
  rownames_to_column('Label')

# Merge module table and read count information.
otu_modules_bac <- inner_join(otu_abundance_bac, modules)

# Initiate an empty dataframe were we can store the total read counts per module.
otu_tab_mod_bac <- data.frame(Label = otu_modules_bac$Label,
                              Mod0 = rep(0,628),
                              Mod1 = rep(0,628),
                              Mod2 = rep(0,628),
                              Mod3 = rep(0,628),
                              Mod4 = rep(0,628),
                              Mod5 = rep(0,628),
                              Mod6 = rep(0,628),
                              Mod7 = rep(0,628))

# If the ASV belongs to a module then write its read count, if not write zero. 
otu_tab_mod_bac$Mod0 <- if_else(otu_modules_bac$modularity_class == 0, otu_modules_bac$total, 0)
otu_tab_mod_bac$Mod1 <- if_else(otu_modules_bac$modularity_class == 1, otu_modules_bac$total, 0)
otu_tab_mod_bac$Mod2 <- if_else(otu_modules_bac$modularity_class == 2, otu_modules_bac$total, 0)
otu_tab_mod_bac$Mod3 <- if_else(otu_modules_bac$modularity_class == 3, otu_modules_bac$total, 0)
otu_tab_mod_bac$Mod4 <- if_else(otu_modules_bac$modularity_class == 4, otu_modules_bac$total, 0)
otu_tab_mod_bac$Mod5 <- if_else(otu_modules_bac$modularity_class == 5, otu_modules_bac$total, 0)
otu_tab_mod_bac$Mod6 <- if_else(otu_modules_bac$modularity_class == 6, otu_modules_bac$total, 0)
otu_tab_mod_bac$Mod7 <- if_else(otu_modules_bac$modularity_class == 7, otu_modules_bac$total, 0)

# set the ASV name back to the rownames.
rownames(otu_tab_mod_bac) <- otu_tab_mod_bac$Label
otu_tab_mod_bac$Label <- NULL

# Sum up the columns to obtain module total read counts.
bacteria_mod_rel <- as.data.frame(colSums(otu_tab_mod_bac))
bacteria_mod_rel

##---------
##  Fungi  
##-------------------------------------------------------------------------------------------------------

# Transform to dataframe to make handling easier.
otu_abundance_fun <- as.data.frame(otu_table(phy_fun_rel))

# Calculate the total read count of each ASV.
otu_abundance_fun$total <- rowSums(otu_abundance_fun)

# Keep only the total read count.
otu_abundance_fun <- otu_abundance_fun %>% 
  select(total) %>%
  rownames_to_column('Label')

# Merge module table and read count information.
otu_modules_fun <- inner_join(otu_abundance_fun, modules)

# Initiate an empty dataframe were we can store the total read counts per module.
otu_tab_mod_fun <- data.frame(Label = otu_modules_fun$Label,
                              Mod0 = rep(0,289),
                              Mod1 = rep(0,289),
                              Mod2 = rep(0,289),
                              Mod3 = rep(0,289),
                              Mod4 = rep(0,289),
                              Mod5 = rep(0,289),
                              Mod6 = rep(0,289),
                              Mod7 = rep(0,289))

# If the ASV belongs to a module then write its read count, if not write zero. 
otu_tab_mod_fun$Mod0 <- if_else(otu_modules_fun$modularity_class == 0, otu_modules_fun$total, 0)
otu_tab_mod_fun$Mod1 <- if_else(otu_modules_fun$modularity_class == 1, otu_modules_fun$total, 0)
otu_tab_mod_fun$Mod2 <- if_else(otu_modules_fun$modularity_class == 2, otu_modules_fun$total, 0)
otu_tab_mod_fun$Mod3 <- if_else(otu_modules_fun$modularity_class == 3, otu_modules_fun$total, 0)
otu_tab_mod_fun$Mod4 <- if_else(otu_modules_fun$modularity_class == 4, otu_modules_fun$total, 0)
otu_tab_mod_fun$Mod5 <- if_else(otu_modules_fun$modularity_class == 5, otu_modules_fun$total, 0)
otu_tab_mod_fun$Mod6 <- if_else(otu_modules_fun$modularity_class == 6, otu_modules_fun$total, 0)
otu_tab_mod_fun$Mod7 <- if_else(otu_modules_fun$modularity_class == 7, otu_modules_fun$total, 0)

# set the ASV name back to the rownames.
rownames(otu_tab_mod_fun) <- otu_tab_mod_fun$Label
otu_tab_mod_fun$Label <- NULL

# Sum up the columns to obtain module total read counts.
fungi_mod_rel <- as.data.frame(colSums(otu_tab_mod_fun))
fungi_mod_rel

##---------
##  Combined  
##-------------------------------------------------------------------------------------------------------

# Set the ASV name back as a column.
algae_mod_rel <- rownames_to_column(algae_mod_rel, var = 'Label')
fungi_mod_rel <- rownames_to_column(fungi_mod_rel, var = 'Label')
bacteria_mod_rel <- rownames_to_column(bacteria_mod_rel, var = 'Label')

# Join the tables of the modules total read counts.
modules_relative <- inner_join(algae_mod_rel, fungi_mod_rel) %>%
  inner_join(bacteria_mod_rel)

# Calculate the relative abundances per organismal group per module. 
modules_relative_only <- modules_relative %>%
  dplyr::rename(total_alg = "colSums(otu_tab_mod_alg)",
         total_fun = "colSums(otu_tab_mod_fun)",
         total_bac = "colSums(otu_tab_mod_bac)") %>%
  mutate(relative_alg = total_alg / sum(algae_mod_rel$`colSums(otu_tab_mod_alg)`),
         relative_fun = total_fun / sum(fungi_mod_rel$`colSums(otu_tab_mod_fun)`),
         relative_bac = total_bac / sum(bacteria_mod_rel$`colSums(otu_tab_mod_bac)`)) %>%
  select(Label, relative_alg, relative_fun, relative_bac) 

# Set the module names as the rownames
rownames(modules_relative_only) <- modules_relative_only$Label
modules_relative_only$Label <- NULL

# Calculate the modules relative abundance.

sort(rowSums(modules_relative_only)/3, decreasing = T)

# Top 2 are the most relevant modules. 
# Module3 ~ 29 %
# Module2 ~ 26 %

####
# Find out which orders are the top ones in the most important modules.
####

# Put ASV table as dataframe to ease handling.
otu_full <- as.data.frame(otu_table(full_physeq))

# Calculate the total ASV read counts. 
otu_full$total <- rowSums(otu_full)

# Subset to only contain the names and total read counts. 
otu_full <- otu_full %>% 
  select(total) %>%
  rownames_to_column('Label')

# Join the modules and ASV total read counts. 
full_modules <- left_join(otu_full, modules)

# Initiate an empty dataframe were we can store the total read counts per module.
otu_tab_mod <- data.frame(Label = full_modules$Label,
                          Mod0 = rep(0,1046),
                          Mod1 = rep(0,1046),
                          Mod2 = rep(0,1046),
                          Mod3 = rep(0,1046),
                          Mod4 = rep(0,1046),
                          Mod5 = rep(0,1046),
                          Mod6 = rep(0,1046),
                          Mod7 = rep(0,1046))


# If the ASV belongs to a module then write its read count, if not write zero. 
otu_tab_mod$Mod0 <- if_else(full_modules$modularity_class == 0, full_modules$total, 0)
otu_tab_mod$Mod1 <- if_else(full_modules$modularity_class == 1, full_modules$total, 0)
otu_tab_mod$Mod2 <- if_else(full_modules$modularity_class == 2, full_modules$total, 0)
otu_tab_mod$Mod3 <- if_else(full_modules$modularity_class == 3, full_modules$total, 0)
otu_tab_mod$Mod4 <- if_else(full_modules$modularity_class == 4, full_modules$total, 0)
otu_tab_mod$Mod5 <- if_else(full_modules$modularity_class == 5, full_modules$total, 0)
otu_tab_mod$Mod6 <- if_else(full_modules$modularity_class == 6, full_modules$total, 0)
otu_tab_mod$Mod7 <- if_else(full_modules$modularity_class == 7, full_modules$total, 0)

# Set ASV names as rownames.
rownames(otu_tab_mod) <- otu_tab_mod$Label
otu_tab_mod$Label <- NULL

# Create the phyloseq object to easily obtain the taxonomic information.
otu_mat <- as.matrix(otu_tab_mod)
tax_mat <- as.matrix(tax_table(full_physeq))

OTU <- otu_table(otu_mat, taxa_are_rows = T)
TAX <- tax_table(tax_mat)

phy_modules <- phyloseq(OTU, TAX)

# Aggregate on Order level.
phy_modules_ord <- phy_modules %>%
  aggregate_taxa('Order')

# Subset by Module name.
phy_mod0_ord <- prune_samples(sample_names(phy_modules_ord) == 'Mod0', phy_modules_ord)
phy_mod1_ord <- prune_samples(sample_names(phy_modules_ord) == 'Mod1', phy_modules_ord)
phy_mod2_ord <- prune_samples(sample_names(phy_modules_ord) == 'Mod2', phy_modules_ord)
phy_mod3_ord <- prune_samples(sample_names(phy_modules_ord) == 'Mod3', phy_modules_ord)
phy_mod4_ord <- prune_samples(sample_names(phy_modules_ord) == 'Mod4', phy_modules_ord)
phy_mod5_ord <- prune_samples(sample_names(phy_modules_ord) == 'Mod5', phy_modules_ord)
phy_mod6_ord <- prune_samples(sample_names(phy_modules_ord) == 'Mod6', phy_modules_ord)
phy_mod7_ord <- prune_samples(sample_names(phy_modules_ord) == 'Mod7', phy_modules_ord)

# Remove any taxa without reads.
phy_mod0_ord_clean <- prune_taxa(taxa_sums(phy_mod0_ord) > 0, phy_mod0_ord)
phy_mod1_ord_clean <- prune_taxa(taxa_sums(phy_mod1_ord) > 0, phy_mod1_ord)
phy_mod2_ord_clean <- prune_taxa(taxa_sums(phy_mod2_ord) > 0, phy_mod2_ord)
phy_mod3_ord_clean <- prune_taxa(taxa_sums(phy_mod3_ord) > 0, phy_mod3_ord)
phy_mod4_ord_clean <- prune_taxa(taxa_sums(phy_mod4_ord) > 0, phy_mod4_ord)
phy_mod5_ord_clean <- prune_taxa(taxa_sums(phy_mod5_ord) > 0, phy_mod5_ord)
phy_mod6_ord_clean <- prune_taxa(taxa_sums(phy_mod6_ord) > 0, phy_mod6_ord)
phy_mod7_ord_clean <- prune_taxa(taxa_sums(phy_mod7_ord) > 0, phy_mod7_ord)

# Subset Module 2 by the organismal group.
phy_mod2_ord_clean_alg <- subset_taxa(phy_mod2_ord_clean, Kingdom == 'Eukaryota')
phy_mod2_ord_clean_bac <- subset_taxa(phy_mod2_ord_clean, Kingdom == 'Bacteria')
phy_mod2_ord_clean_fun <- subset_taxa(phy_mod2_ord_clean, Kingdom == 'Fungi')

# Relative abundance of algae, bacteria and fungi in the module.
sum(abundances(phy_mod2_ord_clean_alg)) / sum(abundances(phy_mod2_ord_clean))
sum(abundances(phy_mod2_ord_clean_bac)) / sum(abundances(phy_mod2_ord_clean))
sum(abundances(phy_mod2_ord_clean_fun)) / sum(abundances(phy_mod2_ord_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 2. 
top_alg_mod2 <- get_top_taxa(phy_mod2_ord_clean_alg, 1, discard_other = T)
tax_table(top_alg_mod2)
# Relative abundance within the module.
abundances(top_alg_mod2) / sum(abundances(phy_mod2_ord_clean))
# Relative abundance within algae in the module. 
abundances(top_alg_mod2) / sum(abundances(phy_mod2_ord_clean_alg))

top_bac_mod2 <- get_top_taxa(phy_mod2_ord_clean_bac, 1, discard_other = T)
tax_table(top_bac_mod2)
# Relative abundance within the module.
abundances(top_bac_mod2) / sum(abundances(phy_mod2_ord_clean))
# Relative abundance within bacteria in the module. 
abundances(top_bac_mod2) / sum(abundances(phy_mod2_ord_clean_bac))

top_fun_mod2 <- get_top_taxa(phy_mod2_ord_clean_fun, 1, discard_other = T)
tax_table(top_fun_mod2)
# Relative abundance within the module.
abundances(top_fun_mod2) / sum(abundances(phy_mod2_ord_clean))
# Relative abundance within fungi in the module. 
abundances(top_fun_mod2) / sum(abundances(phy_mod2_ord_clean_fun))

# Subset Module 3 by the organismal group.
phy_mod3_ord_clean_alg <- subset_taxa(phy_mod3_ord_clean, Kingdom == 'Eukaryota')
phy_mod3_ord_clean_bac <- subset_taxa(phy_mod3_ord_clean, Kingdom == 'Bacteria')
phy_mod3_ord_clean_fun <- subset_taxa(phy_mod3_ord_clean, Kingdom == 'Fungi')

# Relative abundance of algae, bacteria and fungi in the module.
sum(abundances(phy_mod3_ord_clean_alg)) / sum(abundances(phy_mod3_ord_clean))
sum(abundances(phy_mod3_ord_clean_bac)) / sum(abundances(phy_mod3_ord_clean))
sum(abundances(phy_mod3_ord_clean_fun)) / sum(abundances(phy_mod3_ord_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 3.
top_alg_mod3 <- get_top_taxa(phy_mod3_ord_clean_alg, 1, discard_other = T)
tax_table(top_alg_mod3)
# Relative abundance within the module.
abundances(top_alg_mod3) / sum(abundances(phy_mod3_ord_clean))
# Relative abundance within algae in the module.
abundances(top_alg_mod3) / sum(abundances(phy_mod3_ord_clean_alg))

top_bac_mod3 <- get_top_taxa(phy_mod3_ord_clean_bac, 1, discard_other = T)
tax_table(top_bac_mod3)
# Relative abundance within the module.
abundances(top_bac_mod3) / sum(abundances(phy_mod3_ord_clean))
# Relative abundance within bacteria in the module.
abundances(top_bac_mod3) / sum(abundances(phy_mod3_ord_clean_bac))

top_fun_mod3 <- get_top_taxa(phy_mod3_ord_clean_fun, 1, discard_other = T)
tax_table(top_fun_mod3)
# Relative abundance within the module.
abundances(top_fun_mod3) / sum(abundances(phy_mod3_ord_clean))
# Relative abundance within fungi in the module.
abundances(top_fun_mod3) / sum(abundances(phy_mod3_ord_clean_fun))

####
# Find out which taxa are the top ones in the modules.
####

# Subset by Module name.
phy_mod0 <- prune_samples(sample_names(phy_modules) == 'Mod0', phy_modules)
phy_mod1 <- prune_samples(sample_names(phy_modules) == 'Mod1', phy_modules)
phy_mod2 <- prune_samples(sample_names(phy_modules) == 'Mod2', phy_modules)
phy_mod3 <- prune_samples(sample_names(phy_modules) == 'Mod3', phy_modules)
phy_mod4 <- prune_samples(sample_names(phy_modules) == 'Mod4', phy_modules)
phy_mod5 <- prune_samples(sample_names(phy_modules) == 'Mod5', phy_modules)
phy_mod6 <- prune_samples(sample_names(phy_modules) == 'Mod6', phy_modules)
phy_mod7 <- prune_samples(sample_names(phy_modules) == 'Mod7', phy_modules)

# Remove any taxa without reads.
phy_mod0_clean <- prune_taxa(taxa_sums(phy_mod0) > 0, phy_mod0)
phy_mod1_clean <- prune_taxa(taxa_sums(phy_mod1) > 0, phy_mod1)
phy_mod2_clean <- prune_taxa(taxa_sums(phy_mod2) > 0, phy_mod2)
phy_mod3_clean <- prune_taxa(taxa_sums(phy_mod3) > 0, phy_mod3)
phy_mod4_clean <- prune_taxa(taxa_sums(phy_mod4) > 0, phy_mod4)
phy_mod5_clean <- prune_taxa(taxa_sums(phy_mod5) > 0, phy_mod5)
phy_mod6_clean <- prune_taxa(taxa_sums(phy_mod6) > 0, phy_mod6)
phy_mod7_clean <- prune_taxa(taxa_sums(phy_mod7) > 0, phy_mod7)

# Retrieve the top taxa, its taxonomy and relative abundance for Module 0. 
top_mod0 <- get_top_taxa(phy_mod0_clean, 1, discard_other = T)
tax_table(top_mod0)
abundances(top_mod0) / sum(abundances(phy_mod0_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for Module 1. 
top_mod1 <- get_top_taxa(phy_mod1_clean, 1, discard_other = T)
tax_table(top_mod1)
abundances(top_mod1) / sum(abundances(phy_mod1_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for Module 2. 
top_mod2 <- get_top_taxa(phy_mod2_clean, 1, discard_other = T)
tax_table(top_mod2)
abundances(top_mod2) / sum(abundances(phy_mod2_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for Module 3. 
top_mod3 <- get_top_taxa(phy_mod3_clean, 1, discard_other = T)
tax_table(top_mod3)
abundances(top_mod3) / sum(abundances(phy_mod3_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for Module 4. 
top_mod4 <- get_top_taxa(phy_mod4_clean, 1, discard_other = T)
tax_table(top_mod4)
abundances(top_mod4) / sum(abundances(phy_mod4_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for Module 5. 
top_mod5 <- get_top_taxa(phy_mod5_clean, 1, discard_other = T)
tax_table(top_mod5)
abundances(top_mod5) / sum(abundances(phy_mod5_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for Module 6. 
top_mod6 <- get_top_taxa(phy_mod6_clean, 1, discard_other = T)
tax_table(top_mod6)
abundances(top_mod6) / sum(abundances(phy_mod6_clean))

# Retrieve the top taxa, its taxonomy and relative abundance for Module 7. 
top_mod7 <- get_top_taxa(phy_mod7_clean, 1, discard_other = T)
tax_table(top_mod7)
abundances(top_mod7) / sum(abundances(phy_mod7_clean))

#################################################################
##                          Section 6                          ##
##         Drivers of changes in community composition         ##
#################################################################

##---------
##  Algae  
##-------------------------------------------------------------------------------------------------------

# CLR transformation of the raw data. 
phy_algae_clr <- microbiome::transform(phy_algae, transform = 'clr')

###
##Ordination                            -
###

# Principal Component Analysis via phyloseq. RDA is the same as PCA for CLR transformed data. 
phy_algae_ord_clr <- phyloseq::ordinate(phy_algae_clr, "RDA")

# Scree plot to judge the proportion of variance explained by each Principal Components.
phyloseq::plot_scree(phy_algae_ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

# Examine eigenvalues and % prop. variance explained by the Principal Components.
head(phy_algae_ord_clr$CA$eig)     

# Normalize the eigenvalues.
sapply(phy_algae_ord_clr$CA$eig[1:5], function(x) x / sum(phy_algae_ord_clr$CA$eig))  

# Scale axes to the proportion of variance the PC explains.
clr1 <- phy_algae_ord_clr$CA$eig[1] / sum(phy_algae_ord_clr$CA$eig)
clr2 <- phy_algae_ord_clr$CA$eig[2] / sum(phy_algae_ord_clr$CA$eig)

# Plot ordination grouped by management intensity.
alg_ord_plot_manage <- phyloseq::plot_ordination(phy_algae_clr, phy_algae_ord_clr, type="samples", color="intensity_formi") + 
  geom_point(size = .3) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = intensity_formi), linetype = 2) +
  scale_colour_manual(values = ggplot2::alpha(box_cols, 0.9),
                      name = 'Management\nIntensity',
                      labels = c("low", "high")) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.text = element_text(colour = 'black', size = 7),
        axis.title = element_text(size = 7),
        title = element_text(size = 7),
        legend.key = element_rect('transparent'),
        legend.key.size = unit(2, 'mm'),
        legend.position = "left",
        legend.background = element_rect('transparent'), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7))  + 
  labs( subtitle = '(A)')
alg_ord_plot_manage  

# Plot ordination grouped by tree size.
alg_ord_plot_size <- phyloseq::plot_ordination(phy_algae_clr, phy_algae_ord_clr, type="samples", color="size_cat") + 
  geom_point(size = 0.3) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = size_cat), linetype = 2) +
  scale_colour_manual(values = ggplot2::alpha(box_cols, 0.9),
                      name = 'Tree Size',
                      labels = c("large", "medium", "small")) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.text = element_text(colour = 'black', size = 7),
        axis.title = element_text(size = 7),
        title = element_text(size = 7),
        legend.key = element_rect('transparent'),
        legend.key.size = unit(2, 'mm'),
        legend.position = "left", 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        legend.background = element_rect('transparent'))  + 
  labs( subtitle = '(D)')
alg_ord_plot_size

###
##PERMANOVA
###

# Generate the distance matrix. Euclidean distances are the choice for CLR transformed data.
clr_dist_matrix_alg <- phyloseq::distance(phy_algae_clr, method = "euclidean") 

# Run the PERMANOVA analysis thorugh vegans adonis2() including both effects of management and tree sizes.
vegan::adonis2(clr_dist_matrix_alg ~ factor(phyloseq::sample_data(phy_algae_clr)$intensity_formi) +
                 factor(phyloseq::sample_data(phy_algae_clr)$size_cat), by = 'margin')

# Test the within group dispersion for the management intensity.
dispr_management_alg <- vegan::betadisper(clr_dist_matrix_alg, factor(phyloseq::sample_data(phy_algae_clr)$intensity_formi))
dispr_management_alg

plot(dispr_management_alg, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

boxplot(dispr_management_alg, main = "", xlab = "")

# Test for significant differences in dispersion.
permutest(dispr_management_alg)


# Test the within group dispersion for the tree sizes.
dispr_treesize_alg <- vegan::betadisper(clr_dist_matrix_alg, factor(phyloseq::sample_data(phy_algae_clr)$size_cat))
dispr_treesize_alg

plot(dispr_treesize_alg, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

boxplot(dispr_treesize_alg, main = "", xlab = "")

# Test for significant differences in dispersion.
permutest(dispr_treesize_alg)

##---------
##  Bacteria  
##-------------------------------------------------------------------------------------------------------

# CLR transformation of the raw data. 
phy_bacteria_clr <- microbiome::transform(phy_bacteria, transform = 'clr')

###
##Ordination                            -
###

# Principal Component Analysis via phyloseq. RDA is the same as PCA for CLR transformed data. 
phy_bacteria_ord_clr <- phyloseq::ordinate(phy_bacteria_clr, "RDA")

# Scree plot to judge the proportion of variance explained by each Principal Components.
phyloseq::plot_scree(phy_bacteria_ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

# Examine eigenvalues and % prop. variance explained by the Principal Components.
head(phy_bacteria_ord_clr$CA$eig)     

# Normalize the eigenvalues.
sapply(phy_bacteria_ord_clr$CA$eig[1:5], function(x) x / sum(phy_bacteria_ord_clr$CA$eig))  

# Scale axes to the proportion of variance the PC explains.
clr1 <- phy_bacteria_ord_clr$CA$eig[1] / sum(phy_bacteria_ord_clr$CA$eig)
clr2 <- phy_bacteria_ord_clr$CA$eig[2] / sum(phy_bacteria_ord_clr$CA$eig)

# Plot ordination grouped by management intensity.
bac_ord_plot_manage <- phyloseq::plot_ordination(phy_bacteria_clr, phy_bacteria_ord_clr, 
                                                 type="samples", color="intensity_formi") + 
  geom_point(size = .3) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = intensity_formi), linetype = 2) +
  scale_colour_manual(values = ggplot2::alpha(box_cols, 0.9),
                      name = 'Management Intensity') +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.text = element_text(colour = 'black', size = 7),
        axis.title = element_text(size = 7),
        title = element_text(size = 7),
        legend.key = element_rect('transparent'),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        legend.background = element_rect('transparent'), 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9))  + 
  labs( subtitle = '(B)')
bac_ord_plot_manage  

# Plot ordination grouped by tree size.
bac_ord_plot_size <- phyloseq::plot_ordination(phy_bacteria_clr, phy_bacteria_ord_clr, type="samples", color="size_cat") + 
  geom_point(size = .3) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = size_cat), linetype = 2) +
  scale_colour_manual(values = ggplot2::alpha(box_cols, 0.9),
                      name = 'Tree Size') +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.text = element_text(colour = 'black', size = 7),
        axis.title = element_text(size = 7),
        legend.key = element_rect('transparent'),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none", 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        title = element_text(size = 7),
        legend.background = element_rect('transparent'))  + 
  labs( subtitle = '(E)')
bac_ord_plot_size

###
##PERMANOVA
###

# Generate the distance matrix. Euclidean distances are the choice for CLR transformed data.
clr_dist_matrix_bac <- phyloseq::distance(phy_bacteria_clr, method = "euclidean") 

# Run the PERMANOVA analysis thorugh vegans adonis2() including both effects of management and tree sizes.
vegan::adonis2(clr_dist_matrix_bac ~ factor(phyloseq::sample_data(phy_bacteria_clr)$intensity_formi) +
                 factor(phyloseq::sample_data(phy_bacteria_clr)$size_cat),
               by = 'margin', 
               permutations = 10000)

# Test the within group dispersion for the management intensity.
dispr_management_bac <- vegan::betadisper(clr_dist_matrix_bac, factor(phyloseq::sample_data(phy_bacteria_clr)$intensity_formi))
dispr_management_bac

plot(dispr_management_bac, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")


boxplot(dispr_management_bac, main = "", xlab = "")

# Test for significant differences in dispersion.
permutest(dispr_management_bac, pairwise = T)

# Test the within group dispersion for the tree sizes.
dispr_treesize_bac <- vegan::betadisper(clr_dist_matrix_bac, factor(phyloseq::sample_data(phy_bacteria_clr)$size_cat))
dispr_treesize_bac

plot(dispr_treesize_bac, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")


boxplot(dispr_treesize_bac, main = "", xlab = "")

# Test for significant differences in dispersion.
permutest(dispr_treesize_bac)

##---------
##  Fungi  
##-------------------------------------------------------------------------------------------------------

# CLR transformation of the raw data. 
phy_fungi_clr <- microbiome::transform(phy_fungi, transform = 'clr')

###
##Ordination                            -
###

# Principal Component Analysis via phyloseq. RDA is the same as PCA for CLR transformed data.
phy_fungi_ord_clr <- phyloseq::ordinate(phy_fungi_clr, "RDA")

# Scree plot to judge the proportion of variance explained by each Principal Components.
phyloseq::plot_scree(phy_fungi_ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

# Examine eigenvalues and % prop. variance explained by the Principal Components.
head(phy_fungi_ord_clr$CA$eig)     

# Normalize the eigenvalues.
sapply(phy_fungi_ord_clr$CA$eig[1:5], function(x) x / sum(phy_fungi_ord_clr$CA$eig))  

# Scale axes to the proportion of variance the PC explains.
clr1 <- phy_fungi_ord_clr$CA$eig[1] / sum(phy_fungi_ord_clr$CA$eig)
clr2 <- phy_fungi_ord_clr$CA$eig[2] / sum(phy_fungi_ord_clr$CA$eig)

# Plot ordination grouped by management intensity.
fun_ord_plot_manage <- phyloseq::plot_ordination(phy_fungi_clr, phy_fungi_ord_clr, 
                                                 type="samples", color="intensity_formi") + 
  geom_point(size = .3) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = intensity_formi), linetype = 2) +
  scale_colour_manual(values = ggplot2::alpha(box_cols, 0.9),
                      name = 'Management Intensity') +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.text = element_text(colour = 'black', size = 7),
        axis.title = element_text(size = 7),
        title = element_text(size = 7),
        legend.key = element_rect('transparent'),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        legend.background = element_rect('transparent'), 
        legend.text = element_blank(),
        legend.title = element_text(size = 9))  + 
  labs( subtitle = '(C)')
fun_ord_plot_manage  

# Plot ordination grouped by tree size.
fun_ord_plot_size <- phyloseq::plot_ordination(phy_fungi_clr, phy_fungi_ord_clr, type="samples", color="size_cat") + 
  geom_point(size = .3) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = size_cat), linetype = 2) +
  scale_colour_manual(values = ggplot2::alpha(box_cols, 0.9),
                      name = 'Tree Size') +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.text = element_text(colour = 'black', size = 7),
        axis.title = element_text(size = 7),
        title = element_text(size = 7),
        legend.key = element_rect('transparent'),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none", 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        legend.background = element_rect('transparent'))  + 
  labs( subtitle = '(F)')
fun_ord_plot_size

###
##PERMANOVA
###

# Generate the distance matrix. Euclidean distances are the choice for CLR transformed data.
clr_dist_matrix_fun <- phyloseq::distance(phy_fungi_clr, method = "euclidean") 

# Run the PERMANOVA analysis thorugh vegans adonis2() including both effects of management and tree sizes.
vegan::adonis2(clr_dist_matrix_fun ~ factor(phyloseq::sample_data(phy_fungi_clr)$intensity_formi) +
                 factor(phyloseq::sample_data(phy_fungi_clr)$size_cat), by = 'margin',
               permutations = 999)

# Test the within group dispersion for the management intensity.
dispr_management_fun <- vegan::betadisper(clr_dist_matrix_fun,
                                          factor(phyloseq::sample_data(phy_fungi_clr)$intensity_formi))
dispr_management_fun

plot(dispr_management_fun, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")


boxplot(dispr_management_fun, main = "", xlab = "")

# Test for significant differences in dispersion.
permutest(dispr_management_fun)

# Test the within group dispersion for the tree sizes.
dispr_treesize_fun <- vegan::betadisper(clr_dist_matrix_fun, 
                                        factor(phyloseq::sample_data(phy_fungi_clr)$size_cat))
dispr_treesize_fun

plot(dispr_treesize_fun, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")


boxplot(dispr_treesize_fun, main = "", xlab = "")

# Test for significant differences in dispersion.
permutest(dispr_treesize_fun)

#################################################################
##                          Section 7                          ##
##         Differential Abundance Analysis with ALDEx2         ##
#################################################################

##---------
##  Algae  
##-------------------------------------------------------------------------------------------------------

# Extract the taxonomic information so we can later add taxonomic information to the differentially abundant taxa.
taxa_info_alg <- data.frame(tax_table(phy_algae))
taxa_info_alg <- taxa_info_alg %>% rownames_to_column(var = "OTU")

###
#Management intensity
###

# Using the standard ALDEx2 modules.
aldex2_ma_alg <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_algae)),
                               phyloseq::sample_data(phy_algae)$intensity_formi,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_ma_alg, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexmanage_alg <- aldex2_ma_alg %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU, wi.ep, wi.eBH)
sig_aldexmanage_alg <- left_join(sig_aldexmanage_alg, taxa_info_alg)

sig_aldexmanage_alg

###
#Tree size 
###

# Subsetting the data to include only two of the tree sizes, so pairwise t.test assumptions are met.

big_medium_alg <- prune_samples(sample_data(phy_algae)$size_cat %in% c("big", "medium"), phy_algae)

big_small_alg <- prune_samples(sample_data(phy_algae)$size_cat %in% c("big", "small"), phy_algae)

medium_small_alg <- prune_samples(sample_data(phy_algae)$size_cat %in% c("medium", "small"), phy_algae)

# Pair of big and medium trees. 

# Using the standard ALDEx2 modules.
aldex2_bm_alg <- ALDEx2::aldex(data.frame(phyloseq::otu_table(big_medium_alg)),
                               phyloseq::sample_data(big_medium_alg)$size_cat,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_bm_alg, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexbm_alg <- aldex2_bm_alg %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(wi.eBH) %>%
  dplyr::select(OTU,  wi.ep, wi.eBH) 
sig_aldexbm_alg <- left_join(sig_aldexbm_alg, taxa_info_alg)

sig_aldexbm_alg

# Pair of big and small trees. 

# Using the standard ALDEx2 modules.
aldex2_bs_alg <- ALDEx2::aldex(data.frame(phyloseq::otu_table(big_small_alg)),
                               phyloseq::sample_data(big_small_alg)$size_cat,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_bs_alg, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexbs_alg <- aldex2_bs_alg %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU,  wi.ep, wi.eBH)
sig_aldexbs_alg <- left_join(sig_aldexbs_alg, taxa_info_alg)

sig_aldexbs_alg

## Pair of small and medium trees. 

# Using the standard ALDEx2 modules.
aldex2_ms_alg <- ALDEx2::aldex(data.frame(phyloseq::otu_table(medium_small_alg)),
                               phyloseq::sample_data(medium_small_alg)$size_cat,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_ms_alg, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexms_alg <- aldex2_ms_alg %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU, wi.ep, wi.eBH)
sig_aldexms_alg <- left_join(sig_aldexms_alg, taxa_info_alg)

sig_aldexms_alg

##---------
##  Bacteria  
##-------------------------------------------------------------------------------------------------------

# Extract the taxonomic information so we can later add taxonomic information to the differentially abundant taxa.
taxa_info_bac <- data.frame(tax_table(phy_bacteria))

taxa_info_bac <- taxa_info_bac %>% rownames_to_column(var = "OTU")

###
#Management intensity
###

# Using the standard ALDEx2 modules. 
aldex2_ma_bac <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_bacteria)),
                               phyloseq::sample_data(phy_bacteria)$intensity_formi,
                               test="t",
                               effect = TRUE,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_ma_bac, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexmanage_bac <- aldex2_ma_bac %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig_aldexmanage_bac <- left_join(sig_aldexmanage_bac, taxa_info_bac)

sig_aldexmanage_bac

###
#Tree size 
###

# Subsetting the data to include only two of the tree sizes, so pairwise t.test assumptions are met.

big_medium_bac <- prune_samples(sample_data(phy_bacteria)$size_cat %in% c("big", "medium"), phy_bacteria)

big_small_bac <- prune_samples(sample_data(phy_bacteria)$size_cat %in% c("big", "small"), phy_bacteria)

medium_small_bac <- prune_samples(sample_data(phy_bacteria)$size_cat %in% c("medium", "small"), phy_bacteria)

# Pair of big and medium trees. 

# Using the standard ALDEx2 modules.
aldex2_bm_bac <- ALDEx2::aldex(data.frame(phyloseq::otu_table(big_medium_bac)),
                               phyloseq::sample_data(big_medium_bac)$size_cat,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_bm_bac, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexbm_bac <- aldex2_bm_bac %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU, wi.ep, wi.eBH)
sig_aldexbm_bac <- left_join(sig_aldexbm_bac, taxa_info_bac)

sig_aldexbm_bac

# Pair of big and small trees. 

# Using the standard ALDEx2 modules.
aldex2_bs_bac <- ALDEx2::aldex(data.frame(phyloseq::otu_table(big_small_bac)),
                               phyloseq::sample_data(big_small_bac)$size_cat,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_bs_bac, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexbs_bac <- aldex2_bs_bac %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU,  wi.ep, wi.eBH)
sig_aldexbs_bac <- left_join(sig_aldexbs_bac, taxa_info_bac)

sig_aldexbs_bac

## Pair of small and medium trees.

# Using the standard ALDEx2 modules.
aldex2_ms_bac <- ALDEx2::aldex(data.frame(phyloseq::otu_table(medium_small_bac)),
                               phyloseq::sample_data(medium_small_bac)$size_cat,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_ms_bac, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexms_bac <- aldex2_ms_bac %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU,  wi.ep, wi.eBH)
sig_aldexms_bac <- left_join(sig_aldexms_bac, taxa_info_bac)

sig_aldexms_bac

##---------
##  Fungi  
##-------------------------------------------------------------------------------------------------------

# Extract the taxonomic information so we can later add taxonomic information to the differentially abundant taxa.
taxa_info_fun <- data.frame(tax_table(phy_fungi))
taxa_info_fun <- taxa_info_fun %>% rownames_to_column(var = "OTU")

###
#Management intensity
###

# Using the standard ALDEx2 modules.
aldex2_ma_fun <- ALDEx2::aldex(data.frame(phyloseq::otu_table(phy_fungi)),
                               phyloseq::sample_data(phy_fungi)$intensity_formi,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_ma_fun, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexmanage_fun <- aldex2_ma_fun %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU,  wi.ep, wi.eBH)
sig_aldexmanage_fun <- left_join(sig_aldexmanage_fun, taxa_info_fun)

sig_aldexmanage_fun

###
#Tree size 
###

# Subsetting the data to include only two of the tree sizes, so pairwise t.test assumptions are met.

big_medium_fun <- prune_samples(sample_data(phy_fungi)$size_cat %in% c("big", "medium"), phy_fungi)

big_small_fun <- prune_samples(sample_data(phy_fungi)$size_cat %in% c("big", "small"), phy_fungi)

medium_small_fun <- prune_samples(sample_data(phy_fungi)$size_cat %in% c("medium", "small"), phy_fungi)

# Pair of big and medium trees. 

# Using the standard ALDEx2 modules.
aldex2_bm_fun <- ALDEx2::aldex(data.frame(phyloseq::otu_table(big_medium_fun)),
                               phyloseq::sample_data(big_medium_fun)$size_cat,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_bm_fun, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexbm_fun <- aldex2_bm_fun %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU,  wi.ep, wi.eBH)
sig_aldexbm_fun <- left_join(sig_aldexbm_fun, taxa_info_fun)

sig_aldexbm_fun

# Pair of big and small trees. 

# Using the standard ALDEx2 modules.
aldex2_bs_fun <- ALDEx2::aldex(data.frame(phyloseq::otu_table(big_small_fun)),
                               phyloseq::sample_data(big_small_fun)$size_cat,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_bs_fun, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexbs_fun <- aldex2_bs_fun %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU,  wi.ep, wi.eBH)
sig_aldexbs_fun <- left_join(sig_aldexbs_fun, taxa_info_fun)

sig_aldexbs_fun

## Pair of small and medium trees.

# Using the standard ALDEx2 modules.
aldex2_ms_fun <- ALDEx2::aldex(data.frame(phyloseq::otu_table(medium_small_fun)),
                               phyloseq::sample_data(medium_small_fun)$size_cat,
                               test="t",
                               effect = T,
                               denom="all")

# Plot to identify differentially abundant taxa. Red color potentially indicates differentially abundant taxa.
ALDEx2::aldex.plot(aldex2_ms_fun, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Make a table with the differentially abundant taxa with a cutoff p < 0.05.
sig_aldexms_fun <- aldex2_ms_fun %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange( wi.eBH) %>%
  dplyr::select(OTU,  wi.ep, wi.eBH)

sig_aldexms_fun <- left_join(sig_aldexms_fun, taxa_info)

sig_aldexms_fun

#################################################################
##                          Section 8                          ##
##           Plotting the final figures for the paper          ##
#################################################################

##----------------------------------------------------------------
##                      Community Barplots                       -
##----------------------------------------------------------------

# Final size tweaking of the facet labels.
fun_ord_plots <- fun_ord_plots + theme(strip.text.x = element_text(
  size = 6))

bac_ord_plots <- bac_ord_plots + theme(strip.text.x = element_text(
  size = 6))

alg_ord_plots <- alg_ord_plots + theme(strip.text.x = element_text(
  size = 6))

# List the plots and combine them to make the final figure.
barplots <- list(alg_ord_plots, bac_ord_plots, fun_ord_plots)

finished_barplots <- gridExtra::grid.arrange(grobs = barplots, nrow = 3, ncol = 1)

# Save high quality figure.
ggsave('community_barplots.tiff', device = 'tiff', finished_barplots, width = 180, height = 267,
       units = 'mm', dpi = 300)

##----------------------------------------------------------------
##                     Raincloud Plots                           -
##----------------------------------------------------------------

# Set the label names to large instead of big for better grammar.
fun_rain_size <- fun_rain_size + 
  scale_x_discrete(labels=c("big" = "large", "medium" = "medium", "small" = "small")) 

bac_rain_size <- bac_rain_size + 
  scale_x_discrete(labels=c("big" = "large", "medium" = "medium", "small" = "small"))

alg_rain_size <- alg_rain_size + 
  scale_x_discrete(labels=c("big" = "large", "medium" = "medium", "small" = "small"))

# List the plots and combine them to make the final figure.
rainclouds <- list(alg_rain_manage, bac_rain_manage, fun_rain_manage,
                   alg_rain_size, bac_rain_size, fun_rain_size)

rainclouds_finished <- gridExtra::grid.arrange(grobs = rainclouds, nrow = 2, ncol = 3)

# Save high quality figure.
ggsave('alpha_diversity_rainclouds.tiff', device = 'tiff', rainclouds_finished, width = 180,
       units = 'mm', dpi = 300)

##----------------------------------------------------------------
##                     Ordination Plots                          -
##----------------------------------------------------------------

# Extract the legend from the Plots and store it as a ggplot object.
size_leg <- get_legend(alg_ord_plot_size)  
size_leg <- as_ggplot(size_leg)

# Remove the legend from the Plot.
alg_ord_plot_size <- alg_ord_plot_size + 
  theme(legend.position = 'none')

# Extract the legend from the Plots and store it as a ggplot object.
manage_leg <- get_legend(alg_ord_plot_manage)  
manage_leg <- as_ggplot(manage_leg)

# Remove the legend from the Plot.
alg_ord_plot_manage <- alg_ord_plot_manage + 
  theme(legend.position = 'none')

# List the plots and combine them to make the final figure.
ordinations <- list(manage_leg, alg_ord_plot_manage, bac_ord_plot_manage, fun_ord_plot_manage, 
                    size_leg, alg_ord_plot_size, bac_ord_plot_size, fun_ord_plot_size)

ordinations_finished <- gridExtra::grid.arrange(grobs = ordinations, nrow = 2, ncol = 4,
                                                widths = c(0.8,2,2,2))

# Save high quality figure.
ggplot2::ggsave('ordinations_2.tiff', device = "tiff", ordinations_finished, width = 180,
                units = 'mm', dpi = 300)


#################################################################
##                          Section 9                          ##
##                  Export data for publication                ##
#################################################################

# Append the sequence information back to the phyloseq objects.

###
#Bacteria
###
bacteria_rep_seqs_new <- bacteria_rep_seqs
base::rownames(bacteria_rep_seqs_new) <- bacteria_rep_seqs_new$seq_name_bacteria 
bacteria_rep_seqs_new$seq_name_bacteria <- NULL
seq_tax_bac_mat <- as.matrix(bacteria_rep_seqs_new)
seq_tax_bac <- tax_table(seq_tax_bac_mat)

phy_bac_complete <- phy_bacteria
taxa_names(phy_bac_complete) <- gsub("ASV", "ASV_", taxa_names(phy_bac_complete))
phy_bac_complete <- merge_phyloseq(phy_bac_complete, seq_tax_bac)
colnames(tax_table(phy_bac_complete)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Sequence")

# Append the ending _B to the bacterial ASV names to indicate they are bacteria.
taxa_names(phy_bac_complete) <- paste0(taxa_names(phy_bac_complete),'_B')

###
#Fungi
###
fungi_rep_seqs_new <- fungi_rep_seqs
base::rownames(fungi_rep_seqs_new) <- fungi_rep_seqs_new$seq_name_fungi 
fungi_rep_seqs_new$seq_name_fungi  <- NULL
seq_tax_fun_mat <- as.matrix(fungi_rep_seqs_new)
seq_tax_fun <- tax_table(seq_tax_fun_mat)

phy_fun_complete <- phy_fungi
phy_fun_complete <- merge_phyloseq(phy_fun_complete, seq_tax_fun)
colnames(tax_table(phy_fun_complete)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus",  "Species", "Sequence")

# Append the ending _F to the fungal ASV names to indicate they are bacteria.

taxa_names(phy_fun_complete) <- paste0(taxa_names(phy_fun_complete),'_F')

###
#Algae
###
algae_rep_seqs_new <- algae_rep_seqs
base::rownames(algae_rep_seqs_new) <- algae_rep_seqs_new$seq_name_algae 
algae_rep_seqs_new$seq_name_algae   <- NULL
seq_tax_alg_mat <- as.matrix(algae_rep_seqs_new)
seq_tax_alg <- tax_table(seq_tax_alg_mat)

phy_alg_complete <- phy_algae
phy_alg_complete <- merge_phyloseq(phy_alg_complete, seq_tax_alg)
colnames(tax_table(phy_alg_complete)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Sequence")

# Append the ending _A to the algal ASV names to indicate they are bacteria.

taxa_names(phy_alg_complete) <- paste0(taxa_names(phy_alg_complete),'_A')

# Merge them in one phyloseq object. 
phy_complete <- merge_phyloseq(phy_alg_complete,
                               phy_bac_complete,
                               phy_fun_complete)


# Transform the taxonomy table to a dataframe and make rownames the first column. 
complete_taxonomy <- data.frame(tax_table(phy_complete)) %>%
  tibble::rownames_to_column(var = "ASV_ID")

write.csv(complete_taxonomy, "taxonomy_table_full.csv") 


