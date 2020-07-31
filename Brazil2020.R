##################################################################################
## This is code for calculating pairwise SNP distances, Euclidean distances and ##
## geographic distances between M. bovis genomes from Brazil.                   ##
##################################################################################

###############################
### Load required libraries ###
###############################

library(tidyverse)
library(reshape2)
library(fastbaps)
library(phytools)
library(ggmap)
library(ggrepel)
library(geosphere)
library(patchwork)
library(contoureR)
library(ade4)
library(tidygraph)
library(ggraph)
library(igraph)
library(ggtree)

#############################
### Set working directory ###
#############################

BrazilDir <- "/Users/avt/Documents/Brazil/Rodrigues/"
setwd(BrazilDir)

####################
### Read in data ###
####################

######################
## Read in metadata ##
######################

BrazilMetadata <- read_csv("Brazil_metadata_new.csv")

########################################################################
## Read in pairwise SNP distances and convert data into useful format ##
########################################################################

BrazilPWsnps <- read_csv("New_data/Brazil_all_masked_snps.csv")

BrazilPWsnpsDecon <- data.frame( t(combn(names(BrazilPWsnps),2)), dist=t(BrazilPWsnps)[lower.tri(BrazilPWsnps)] )

colnames(BrazilPWsnpsDecon) <- c("Taxon1", "Taxon2", "dist")

########################################
## Read in midpoint rooted phylogeny  ##
########################################

BrazilTree <- read.newick("New_data/Brazil_all_masked_snps.aln.treefile")

BrazilTreeRooted <- midpoint.root(BrazilTree)

############################################
## Import SNP alignment (regions removed) ##
############################################

BrazilAlignment <- import_fasta_sparse_nt("New_data/Brazil_all_masked_snps.fasta")

################
### Analyses ###
################

#########################################################
### Run fastBAPS on Brazil dataset to cluster genomes ###
### Constrain clusters to phylogeny                   ###
#########################################################

BrazilBestPartition <- best_baps_partition(BrazilAlignment, BrazilTreeRooted)

#########################################################
## Create dataframe with samples and fastBAPS clusters ##
#########################################################

BrazilFastBAPS <- data.frame(id = BrazilTreeRooted$tip.label, fastbaps = BrazilBestPartition, stringsAsFactors = FALSE)

##########################################
## Add fastBAPS clusters to metadata df ##
##########################################

BrazilMetadataClusters <- BrazilMetadata %>%
  left_join(BrazilFastBAPS, by = c("Lane_id" = "id"))

#################################################################
## Write out csv file with metadata to create iTOL data strips ##
#################################################################

BrazilMetadataClustersWrite <- BrazilMetadataClusters %>% 
  select(Lane_id, Spoligotype, Clonal_complex, Herd, fastbaps)

write_csv(BrazilMetadataClustersWrite, "New_data/Brazil_new_iTOL_metadata.csv")

########################################################################################
## Create a plot of geographic coordinates and colour by fastbaps cluster (Figure 2B) ##
########################################################################################

BrazilCoordPlot <- ggplot(data = BrazilMetadataClusters) +
  geom_point(aes(x = Longitude, y = Latitude, color = as.factor(fastbaps))) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  labs(color = "Cluster")

#########################################
## Define study lat/long for stamenmap ##
#########################################

BrazilCoords <- c(left = -53.5, bottom = -31, right = -51, top = -28)

##############################################
## Download stamenmap for above coordinates ##
##############################################

Brazilstamenmap <- get_stamenmap(BrazilCoords, zoom = 8, maptype = "terrain-background", color = "bw", force = TRUE)

#########################################################
## Extract single example of herd locations for labels ##
#########################################################

BrazilFarmCoords <- BrazilMetadataClusters %>% 
  group_by(Herd) %>% 
  filter(row_number()==1)

######################################################
## Plot map and overlay herd locations (Figure 2A ) ##
######################################################

BrazilMapPlot <- ggmap(Brazilstamenmap) + 
  geom_point(data = BrazilMetadataClusters, aes(x = Longitude, y = Latitude, color = Herd), alpha =0.5, fill = NA) +
  geom_label_repel(data = BrazilFarmCoords, aes(x = Longitude, y = Latitude, label = Herd)) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(colour = 'Herd') + 
  theme(legend.position = "none")

##############################################################
## Calculate and plot convex hull for each fastbaps cluster ##
##############################################################

subset <- NULL

hull <- NULL

polygon <- NULL

for(cluster in unique(BrazilMetadataClusters$fastbaps)){
  
  ######################################################
  # Get the subset of the data for the current cluster #
  ######################################################
  
  subset[[cluster]] <- BrazilMetadataClusters[BrazilMetadataClusters$fastbaps == cluster, ]
  
  ##############################################
  # Calculate the convex hull for each cluster #
  ##############################################
  
  hull[[cluster]] <- getConvexHull(x=subset[[cluster]]$Longitude, y=subset[[cluster]]$Latitude)
  
  ##################################################
  # Calculate polygon coordinates for each cluster #
  ##################################################
  
  polygon[[cluster]] <- cbind(subset[[cluster]][hull[[cluster]], "Longitude"], subset[[cluster]][hull[[cluster]], "Latitude"])
  
}

#################################################
## Collapse polygon list into single dataframe ##
#################################################

PolygonAll <- bind_rows(polygon, .id = "Cluster")

##############################################################
## Overlay polygons on plot of coordinates for each isolate ##
##############################################################

BrazilCoordPlotPolygon <- ggplot() +
  geom_point(data = BrazilMetadataClusters, aes(x = Longitude, y = Latitude, color = as.factor(fastbaps))) +
  geom_polygon(data = PolygonAll, aes(x = Longitude, y = Latitude, group = Cluster, fill = Cluster), alpha = 0.3) +
  theme_bw() +
  labs(color = "fastBAPS", fill = "fastBAPS")

##########################################
### Geographical localization analyses ###
##########################################

################################################################
## Create a copy of deconvoluted pairwise SNP distance matrix ##
################################################################

BrazilSNPsGeoDists <- BrazilPWsnpsDecon

#########################################################################
## Create a new dataframe keeping metadata required for later analyses ##
#########################################################################

BrazilMetadataClustersGeoDists <- BrazilMetadataClusters %>% 
  select(Lane_id, Isolate_ID, Herd, Animal, Spoligotype, Clonal_complex, Longitude, Latitude,  fastbaps)

###########################################################
## Add metadata for each taxon in pairwise SNP dataframe ##
###########################################################

BrazilSNPsGeoDists <- BrazilSNPsGeoDists %>% 
  left_join(BrazilMetadataClustersGeoDists, by = c("Taxon1" = "Lane_id")) %>% 
  left_join(BrazilMetadataClustersGeoDists, by = c("Taxon2" = "Lane_id"))

######################################################
## Calculate pw geographical distance for each pair ##
######################################################

BrazilSNPsGeoDists <- BrazilSNPsGeoDists %>% 
  mutate(geo_dist = distHaversine(cbind(Longitude.x, Latitude.x), cbind(Longitude.y, Latitude.y)))

#################################################################
## Create new dataframe for geographical localization analyses ##
#################################################################

BrazilSNPsGeoDistsFiltered <- BrazilSNPsGeoDists %>%
  select(Taxon1, Taxon2, Isolate_ID.x, Isolate_ID.y, Animal.x,  Animal.y,Herd.x, Herd.y, fastbaps.x, fastbaps.y, dist, geo_dist) %>% 
  mutate(herd_match = ifelse(Herd.x == Herd.y, paste("Within herd"), paste("Between herd"))) %>% 
  mutate(baps_match = ifelse(fastbaps.x == fastbaps.y, paste("Same fastbaps"), paste("Different fastbaps"))) %>% 
  mutate(animal_match = ifelse(Animal.x == Animal.y, paste("Within host"), paste("Between host")))

################################################################################
## Plot histogram of pairwise SNP distances, colour by herd_match (Figure 3A) ##
################################################################################

BrazilSNPHistogram <- ggplot(data = BrazilSNPsGeoDistsFiltered, aes(x = dist, fill = herd_match)) + 
  geom_histogram(bins = 200, alpha = 0.8, position = "identity") + 
  theme_bw()  + 
  xlab("Pairwise SNP distance") + 
  ylab("Frequency") + 
  labs(fill = "Within/between herd") + 
  theme(legend.position = "none")

#####################################################################
## Create new dataframe to keep only within herd pairs of isolates ##
#####################################################################

BrazilSNPsGeoDistsFilteredSameFarm <- BrazilSNPsGeoDistsFiltered %>% 
  filter(herd_match == "Within herd")

################################################################
## Create box plot of pairwise distances per herd (Figure 3B) ##
################################################################

BrazilSNPsFarmBox <- ggplot(data = BrazilSNPsGeoDistsFilteredSameFarm, 
                            aes(x = Herd.x, y = dist, fill = Herd.x)) + 
  geom_jitter(color="black", size=0.4, alpha=0.6) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() + 
  xlab("Herd") + 
  ylab("Pairwise SNP distance") + 
  theme(legend.position = "none") 

####################################################################################
## Create new dataframe to keep only pairs of isolates from the same BAPS cluster ##
####################################################################################

BrazilSNPsGeoDistsFilteredSameBAPS <- BrazilSNPsGeoDistsFiltered %>% 
  filter(baps_match == "Same fastbaps")

##################################################
## Conduct Mantel test on each fastBAPS cluster ##
##################################################

BAPSsubset <- NULL

BAPSSNPdists <- NULL

BAPSGeoDists <- NULL

BAPSMantel <- NULL

for(cluster in unique(BrazilSNPsGeoDistsFilteredSameBAPS$fastbaps.x)){
  
  ######################################################
  # Get the subset of the data for the current cluster #
  ######################################################
  
  BAPSsubset[[cluster]] <- BrazilSNPsGeoDistsFilteredSameBAPS[BrazilSNPsGeoDistsFilteredSameBAPS$fastbaps.x == cluster, ]
  
  #########################################################
  # Create a matrix of the SNP distances for each cluster #
  #########################################################
  
  BAPSSNPdists[[cluster]] <- dist(BAPSsubset[[cluster]]$dist)
  
  ################################################################
  # Create a matrix of the geographic distances for each cluster #
  ################################################################
  
  BAPSGeoDists[[cluster]] <- dist(BAPSsubset[[cluster]]$geo_dist)
  
  #######################################
  # Conduct Mantel test on each cluster #
  #######################################
  
  BAPSMantel[[cluster]] <- mantel.rtest(BAPSSNPdists[[cluster]], BAPSGeoDists[[cluster]], nrepet = 999)
}

############################################################
## Plot pairwise SNP distance against geographic distance ##
## and facet wrap by BAPS cluster (Figure 3C)             ##
############################################################

BrazilDistancePlot <- ggplot(data = BrazilSNPsGeoDistsFilteredSameBAPS) + 
  geom_point(aes(x = dist, y = geo_dist/1000)) +
  xlab("Pairwise SNP distance") +
  ylab("Geographic distance (kilometers)") +
  theme_bw() +
  facet_wrap(~fastbaps.x, ncol = 2) + 
  theme(legend.position = "none")

###################################
## Examine within host diversity ##
###################################

BrazilSNPsGeoDistsFilteredSameHost <- BrazilSNPsGeoDistsFiltered %>% 
  filter(animal_match == "Within host") %>% 
  select(Taxon1, Taxon2, Animal.x, Herd.x, dist, animal_match)

BrazilSNPsGeoDistsFilteredSameHostSummary <- BrazilSNPsGeoDistsFilteredSameHost %>% 
  group_by(animal.x) %>% 
  summarize(min_snp_dist = min(dist)) %>% 
  summarize(max_snp_dist = max(dist)) 

###############################################################
### Create network of all samples using SNP threshold of 15 ###
###############################################################

################################################
## Create new df containing nodes for network ##
################################################

BrazilMetadataClusters_15_snps <- BrazilMetadataClusters %>%
  mutate(id = row_number()) %>% 
  select(id, Lane_id, Herd, Clonal_complex, fastbaps, Longitude, Latitude)

##############################################################
## Extract pairwise SNP comparisons < 16 (edges of network) ##
##############################################################

BrazilSNPsGeoDistsFiltered_edges_15_snp <- BrazilSNPsGeoDistsFiltered %>%
  filter(dist < 16) %>% 
  select(Taxon1, Taxon2, dist)

####################################
## Create network edges dataframe ##
####################################

BrazilSNPsGeoDistsFiltered_edges_15_snp <- BrazilSNPsGeoDistsFiltered_edges_15_snp %>%
  left_join(BrazilMetadataClusters_15_snps[,c(1,2)], by = c("Taxon1" = "Lane_id")) %>% 
  dplyr::rename(from.id = id)

BrazilSNPsGeoDistsFiltered_edges_15_snp <- BrazilSNPsGeoDistsFiltered_edges_15_snp %>% 
  left_join(BrazilMetadataClusters_15_snps[,c(1,2)], by = c("Taxon2" = "Lane_id")) %>% 
  select(from.id, to.id = id, dist)

###########################
## Create network object ##
###########################

Brazil_15_snp_routes <- tbl_graph(nodes = BrazilMetadataClusters_15_snps, edges = BrazilSNPsGeoDistsFiltered_edges_15_snp, directed = TRUE)

#############################################################
## Create new nodes dataframe that includes the network id ##
#############################################################

BrazilMetadataClusters_15_snpsNetworks <- BrazilMetadataClusters_15_snps %>% 
  mutate(network_id = components(Brazil_15_snp_routes)$membership)

##########################################################
## Plot the network coloured by host and clonal complex ##
## as shape (Figure 3D)                                 ##
##########################################################

Brazil_15_snp_network_map_host <- ggraph(Brazil_15_snp_routes, layout = "nicely") + 
  geom_edge_link(edge_colour = 'black', edge_alpha = 0.8, edge_width = 0.4, 
                 angle_calc = "along", label_dodge = unit(2.5,'mm'), label_size = 3) + 
  geom_node_point(aes(colour = Herd, shape = Clonal_complex), size = 4) + 
  theme_graph() + 
  labs(colour = "Herd", shape = "Clonal complex")

#############################
## Calculate network stats ##
#############################

BrazilMetadataClustersNetworks <- BrazilMetadataClusters_15_snpsNetworks %>% 
  group_by(network_id) %>% 
  filter(n() > 1) %>%
  ungroup() %>% 
  select(Lane_id, Herd, network_id)

BrazilSNPsGeoDistsFiltered_15_snp <- BrazilSNPsGeoDistsFiltered %>%
  filter(dist < 16) %>% 
  select(Taxon1, Taxon2, dist)

BrazilMetadataClustersNetworksDists <- BrazilSNPsGeoDistsFiltered_15_snp %>% 
  left_join(BrazilMetadataClustersNetworks, by = c("Taxon1" = "Lane_id")) %>% 
  left_join(BrazilMetadataClustersNetworks, by = c("Taxon2" = "Lane_id")) %>% 
  mutate(herd_match = ifelse(Herd.x == Herd.y, paste("Within herd"), paste("Between herd"))) %>% 
  filter(herd_match == "Within herd")

BrazilMetadataClustersNetworksDistsSummary <- BrazilMetadataClustersNetworksDists %>% 
  group_by(network_id.x) %>% 
  summarize(min_snp_dist = min(dist))

##############################################
## Zoom in on network clusters 8, 18 and 20 ##
##############################################

############################################################################
## Create new nodes dataframe filtering out network clusters 8, 18 and 20 ##
############################################################################

BrazilMetadataClusters_15_snpsNetworksZoom <- BrazilMetadataClusters_15_snpsNetworks %>% 
  filter(network_id == 8 |  network_id == 18 | network_id == 20)%>%
  mutate(id = row_number())

##############################################################
## Extract pairwise SNP comparisons < 16 (edges of network) ##
##############################################################

BrazilSNPsGeoDistsFiltered_edges_15_snpZoom <- BrazilSNPsGeoDistsFiltered %>%
  filter(dist < 16) %>% 
  select(Taxon1, Taxon2, dist) 

########################################
## Create new network edges dataframe ##
########################################

BrazilSNPsGeoDistsFiltered_edges_15_snpZoom <- BrazilSNPsGeoDistsFiltered_edges_15_snpZoom %>%
  left_join(BrazilMetadataClusters_15_snpsNetworksZoom[,c(1,2)], by = c("Taxon1" = "Lane_id")) %>% 
  dplyr::rename(from.id = id)%>% 
  filter(!is.na(from.id))

BrazilSNPsGeoDistsFiltered_edges_15_snpZoom <- BrazilSNPsGeoDistsFiltered_edges_15_snpZoom %>% 
  left_join(BrazilMetadataClusters_15_snpsNetworksZoom[,c(1,2)], by = c("Taxon2" = "Lane_id")) %>% 
  select(from.id, to.id = id, dist)

###############################
## Create new network object ##
###############################

Brazil_15_snp_routesZoom <- tbl_graph(nodes = BrazilMetadataClusters_15_snpsNetworksZoom, edges = BrazilSNPsGeoDistsFiltered_edges_15_snpZoom, directed = TRUE)

#####################################################################
## Plot network clusters 8, 18 and 20 coloured by host (Figure 3E) ##
#####################################################################

Brazil_15_snp_network_map_hostZoom <- ggraph(Brazil_15_snp_routesZoom, layout = "nicely") + 
  geom_edge_link(aes(label = dist), edge_colour = 'black', edge_alpha = 0.8, edge_width = 0.4, 
                 angle_calc = "along", label_dodge = unit(2.5,'mm'), label_size = 3) + 
  geom_node_point(aes(colour = Herd), size = 4) + 
  theme_graph() + 
  labs(colour = "Herd")

#############################
### Eu2 specific analyses ###
#############################

##############################################
## Read in TreeTime time-scaled newick tree ##
##############################################

EU2treetime <- read.tree("New_data/nextstrain_all_Eu2_290720_timetree.nwk")

###########################################
## Read in metadata for all Eu2 isolates ##
###########################################

EU2Metadata <- read_csv("New_data/EU2_metadata.csv")

###########################################################
## Create new dataframe ready for adding metadata strips ##
## using ggtree                                          ##
###########################################################

EU2MetadataDF <- as.data.frame(EU2Metadata[,3:5])

rownames(EU2MetadataDF) = EU2Metadata$isolate_id

##########################################################
## Relabel tips change sequence file ids to isolate ids ##
##########################################################

EU2treetimeTipLabels <- as.data.frame(EU2treetime$tip.label)

EU2treetimeTipLabelsRelabelled <- EU2treetimeTipLabels %>% 
  left_join(EU2Metadata[,1:2], by = c("EU2treetime$tip.label" = "lane_id"))

EU2treetime$tip.label <- EU2treetimeTipLabelsRelabelled$isolate_id

#################################################################
## Plot TreeTime tree adding lines for time period 1200 - 2020 ##
#################################################################

EU2TTggtree <- ggtree(EU2treetime, right=TRUE, mrsd="2019-01-01") +
  geom_tiplab(size = 1, align = TRUE, offset = 5) +
  theme_tree2() +
  scale_x_continuous(breaks=seq(1200, 2020, 50)) +
  theme(panel.grid.major = element_line(color="black", size=.2),
        panel.grid.minor   = element_line(color="grey", size=.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  xlim_tree(2020)

########################################################
## Add metadata strips to above tree plot (Figure 4B) ##
########################################################

EU2TTggtreeMeta <- gheatmap(EU2TTggtree, EU2MetadataDF, colnames_position = "top", font.size=3, offset=25, width=0.15, colnames_angle = 90, hjust = 0)