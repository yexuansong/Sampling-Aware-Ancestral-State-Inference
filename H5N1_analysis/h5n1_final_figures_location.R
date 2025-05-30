##### Comparing results from ape/BiSSE/our implementations

# load a bunch of packages that might come in handy 
library("ape")
library("deSolve")
library("strap")
library("stringr")
library("ips")
require("treeio")
library("ggplot2")
require("ggtree")
library("reshape2")
library("ggstance")
library("diversitree")
library("ggimage")
library("ggpubr")
library("dplyr")
library("patchwork")
library("tidyverse")
library("janitor")
library("stringr")
library("skimr")
library("phytools")
library(dplyr)
library(rstatix)
library(jsonlite)

source(file.path("h5n1_phy_data_cleaning.R"))

####################
# unmodified phylogeny
unmod <- ggtree(tree)

# modified phylogeny
tree_modified <- di2multi(tree, tol = 1e-08)

# now since we are trying to compare the same tree, maybe it's a good idea
# to also modify the tree as well (by deleting human tip, only one), and
# deleting the root (ace estimates complex transition rates)

drop_tip <- c(1,104)
tree_modified <- drop.tip(tree_modified,drop_tip)
tree_info <- tree_info[-c(1,104),]

tips_to_drop <- grep("human", tree_modified$tip.label, value = TRUE, ignore.case = TRUE)
tree_modified <- drop.tip(tree_modified, tips_to_drop)
tree_info <- tree_info %>% 
  filter(species != "human")

# plot the phylogeny and color the tips based on species
phy_species <- ggtree(tree_modified)
phy_species <- phy_species %<+% tree_info + geom_tippoint(aes(color=species))

# plot the phylogeny and color the tips based on locations
phy_locations <- ggtree(tree_modified)
phy_locations <- phy_locations %<+% tree_info + geom_tippoint(aes(color=location))

# plot the phylogeny and color all the nodes based on beast phy
phy_beast <- ggtree(tree_modified)
phy_beast <- phy_beast %<+% tree_withnode + geom_point(aes(color=nodelabels))

#ggsave("phy_locations.pdf", phy_locations, width = 16, height = 9, units = "in", dpi = 300)

####################
# Analysis (Ancestral State Reconstruction ASR)
####################

####################
# ASR using ML
# adding tip labels (species only)
ace_phy <- multi2di(tree_modified)
ace_phy$state <- tree_info$location

ans_ace<-ace(ace_phy$state, ace_phy,type = "discrete",method = "ML", model="ER")

# ace estimate the transition rate
qij_rate <- ans_ace$rates

# ace infers the internal nodes state/species probabilities
ace_node_lik <- as.data.frame(ans_ace$lik.anc)

# now plot the inferred phylogeny 

# insert the likelihood pie chart to the tree
ace_node_lik$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
ace_pie <- nodepie(ace_node_lik,cols=1:13)

# plot the phylogeny and color the tips based on species
lo <- c( 
  "California", "Indiana", "Kansas","Maryland","Michigan","Minnesota","Montana",
  "New_Mexico","North_Carolina","Ohio","Oklahoma","Texas","Wisconsin")

#ace_species <- ggtree(ace_phy)
# adding more features to the phylogeny
ace_species <- ggtree(ace_phy,mrsd = "2024-04-06",options(ignore.negative.edge=TRUE))
ace_species <- ace_species %<+% tree_info + geom_tippoint(aes(color=location))+
  theme_tree2()+ggtitle("Maximum Likelihood (ace) Location") +
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(x = "Time")


p1 <- inset(ace_species, ace_pie,width = 0.025,height = 0.025,hjust=0.005)
#ggsave("ace_h5n1_location.pdf", p1, width = 16, height = 9, units = "in", dpi = 300)
# 


############################
# now plot the alluvial plot
# Extract node states

node_states <- apply(ans_ace$lik.anc, 1, which.max)
node_states <- factor(node_states, levels = 1:14, labels = c("Arizona", 
                                                             "California", "Indiana", "Kansas","Maryland","Michigan","Minnesota","Montana",
                                                             "New_Mexico","North_Carolina","Ohio","Oklahoma","Texas","Wisconsin"))

all_states <- c(ace_phy$state, node_states)

# Create a data frame for plotting
plot_data <- data.frame(
  node = c(1:length(ace_phy$state), 1:ace_phy$Nnode + length(ace_phy$state)),
  state = all_states,
  is_tip = c(rep(TRUE, length(ace_phy$state)), rep(FALSE, ace_phy$Nnode))
)

# Add parent information
plot_data$parent <- sapply(plot_data$node, function(n) {
  if (n == length(ace_phy$state) + 1) return(NA)  # root node
  parent <- ace_phy$edge[ace_phy$edge[,2] == n, 1]
  if (length(parent) == 0) return(NA)  # tip
  return(parent)
})

# Create a data frame for alluvial plot
alluvial_data <- plot_data[!plot_data$is_tip,]
alluvial_data <- alluvial_data[!is.na(alluvial_data$parent),]
alluvial_data$parent_state <- plot_data$state[match(alluvial_data$parent, plot_data$node)]

# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data, cols = c(state, parent_state), 
                              names_to = "generation", values_to = "location")

# Create the alluvial plot
alluvial_ace <- ggplot(alluvial_long,
                       aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "Generation", y = "Count", fill = "Locations") +
  ggtitle("Movement between location - ace") +
  scale_fill_discrete(labels = c(lo[1],lo[10],lo[11],lo[12],lo[13],lo[2],lo[3],lo[4],lo[5],lo[6],lo[7],lo[9]))

#ggsave("alluvial_ace_location.pdf", alluvial_ace, width = 16, height = 9, units = "in", dpi = 300)

####################
# ASR using Stochastic Character Mapping (Phytool)

# fit the model
simmap_phy <- multi2di(tree_modified)
simmap_phy$tip.label = ace_phy$state

lab <- simmap_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab

# ####################
# ASR using our method that accounts for sampling bias

###### Our method

#set up the parameters
ans_ace$rates[ans_ace$rates==0] <- 0.1
qij_rate <- ans_ace$rates


replace_matrix_with_vector <- function(matrix, vector) {
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      matrix[i,j] <- vector[matrix[i,j]]
    }
  }
  return(matrix)
}


#assuming texas 2 times more 

qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)
mm<-qij_matrix

pars <- c(21.1,  21.1, 21.1, 21.1, 21.1, 21.1,21.1,21.1,21.1,21.1,21.1,21.1,21.1,    
          6.7, 6.7, 6.7, 6.7, 6.7, 6.7,6.7,6.7,6.7,6.7,6.7,6.7,6.7, 
          5,  5,  5,  5,  5,  5, 5,5,5,5,5,5,5) 

vect <- rep(0.1,13*13-13)
pars <- c(pars,vect)


our_phy <- multi2di(tree_modified)
lab <- simmap_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab

our_phy$tip.state <- numeric_vector

asr_our <- saasi(pars,our_phy,mm)


colnames(asr_our) = c(1,10,11,12,13,2,3,4,5,6,7,8,9)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:13,color=c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
  "#FFD92F", "#E5C494", "#B3B3B3"
))

our_species <- ggtree(ace_phy,mrsd = "2024-04-06",options(ignore.negative.edge=TRUE))
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=location))+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
    "#FFD92F", "#E5C494", "#B3B3B3"
  ),labels = c("New_Mexico" = "New Mexico", "North_Carolina" = "North Carolina"))+
  theme_tree2()+ggtitle("A: Equal sampling") +labs(color = "Location")+
  geom_hilight(node=120, fill="purple",alpha=0.05) +
  geom_vline(xintercept = 2023.97, color = "red", linetype = "dashed", size = 0.5)+
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15), 
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))


p3 <- inset(our_species, our_pie,width = 0.07,height = 0.07,hjust=0.005)

#ggsave("our_h5n1_location_texas_more.pdf", p3,  width = 210, height = 297, units = "mm")

############################
# now plot the alluvial plot
# Extract node states

asr_our <- saasi(pars,our_phy,mm)

node_states <- apply(asr_our, 1, which.max)
node_states <- factor(node_states, levels = 1:13, labels = c(
  "California", "Indiana", "Kansas","Maryland","Michigan","Minnesota","Montana",
  "New_Mexico","North_Carolina","Ohio","Oklahoma","Texas","Wisconsin"))

all_states <- c(ace_phy$state, node_states)
# 
# Create a data frame for plotting
plot_data <- data.frame(
  node = c(1:length(ace_phy$state), 1:ace_phy$Nnode + length(ace_phy$state)),
  descendant = all_states,
  is_tip = c(rep(TRUE, length(ace_phy$state)), rep(FALSE, ace_phy$Nnode))
)
# 
# Add parent information
plot_data$parent <- sapply(plot_data$node, function(n) {
  if (n == length(ace_phy$state) + 1) return(NA)  # root node
  parent <- ace_phy$edge[ace_phy$edge[,2] == n, 1]
  if (length(parent) == 0) return(NA)  # tip
  return(parent)
})


node_to_node <- plot_data[!plot_data$is_tip & !is.na(plot_data$parent),]
node_to_node$ancestor <- plot_data$descendant[match(node_to_node$parent, plot_data$node)]

# 2. Now get node to tip transitions (new addition)
node_to_tip <- plot_data[plot_data$is_tip & !is.na(plot_data$parent),]
node_to_tip$ancestor <- plot_data$descendant[match(node_to_tip$parent, plot_data$node)]

# 3. Combine both types of transitions
alluvial_data <- rbind(node_to_node, node_to_tip)

# 4. Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data, cols = c(descendant, ancestor),
                              names_to = "generation", values_to = "location")

unique_location <- unique(alluvial_long$location)
lo <- c(
  "California", "Indiana", "Kansas", "Maryland", "Michigan", "Minnesota", "Montana",
  "New_Mexico", "North_Carolina", "Ohio", "Oklahoma", "Texas", "Wisconsin")

alluvial_long$location[alluvial_long$location == lo[1]] <- 1
alluvial_long$location[alluvial_long$location == lo[2]] <- 2
alluvial_long$location[alluvial_long$location == lo[3]] <- 3
alluvial_long$location[alluvial_long$location == lo[4]] <- 4
alluvial_long$location[alluvial_long$location == lo[5]] <- 5
alluvial_long$location[alluvial_long$location == lo[6]] <- 6
alluvial_long$location[alluvial_long$location == lo[7]] <- 7
alluvial_long$location[alluvial_long$location == lo[8]] <- 8
alluvial_long$location[alluvial_long$location == lo[9]] <- 9
alluvial_long$location[alluvial_long$location == lo[10]] <- 10
alluvial_long$location[alluvial_long$location == lo[11]] <- 11
alluvial_long$location[alluvial_long$location == lo[12]] <- 12
alluvial_long$location[alluvial_long$location == lo[13]] <- 13


# Create the alluvial plot
# Create the alluvial plot without the alpha aesthetic
alluvial_our3 <- ggplot(alluvial_long,
                       aes(x = generation, stratum = location, alluvium = node, 
                           fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Location") +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  ggtitle("C") +
  scale_fill_manual(values = c(
    "#E41A1C", "#66C2A5", "#FC8D62", "#8DA0CB", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#FFFF33", "#F781BF", "#999999","#1B9E77" ,"#D95F02"),
    labels=c(lo[1], lo[10], lo[11], lo[12],lo[13], lo[2], lo[3], lo[4], lo[5], lo[6],lo[7], lo[8], lo[9])
  ) +
  theme(text = element_text(size = 15, family = "serif"), plot.title = element_text(size=15),legend.position = "none")
# 


#######################
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)
mm<-qij_matrix

pars <- c(21.1,  21.1, 21.1, 21.1, 21.1, 21.1,21.1,21.1,21.1,21.1,21.1,21.1,21.1,    
          6.7, 6.7, 6.7, 6.7, 6.7, 6.7,6.7,6.7,6.7,6.7,6.7,6.7,6.7, 
          2,  2,  2,  2,  2,  2, 2,2,2,2,2,10,2) 
vect <- rep(0.1,13*13-13)
pars <- c(pars,vect)
our_phy <- multi2di(tree_modified)
lab <- simmap_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab
our_phy$tip.state <- numeric_vector
#mm<-mm/1.4
asr_our <- saasi(pars,our_phy,mm)

# internal_node_ids <- root_node:nnode

colnames(asr_our) = c(1,10,11,12,13,2,3,4,5,6,7,8,9)

#colnames(asr_our) = c(1:13)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:13,color=c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
  "#FFD92F", "#E5C494", "#B3B3B3"
))

our_species <- ggtree(ace_phy,mrsd = "2024-04-06",options(ignore.negative.edge=TRUE))
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=location))+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
    "#FFD92F", "#E5C494", "#B3B3B3"
  ),labels = c("New_Mexico" = "New Mexico", "North_Carolina" = "North Carolina"))+
  theme_tree2()+ggtitle("B: Texas at 5x sampling") +labs(color = "Location")+
  geom_hilight(node=120, fill="purple",alpha=0.05) +
  geom_vline(xintercept = 2023.97, color = "red", linetype = "dashed", size = 0.5)+
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15), 
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))


p4 <- inset(our_species, our_pie,width = 0.07,height = 0.07,hjust=0.005)

#ggsave("our_h5n1_location_equal_sampling.pdf", p4, width = 210, height = 297, units = "mm")

############################
# now plot the alluvial plot
# Extract node states

asr_our <- saasi(pars,our_phy,mm)

node_states <- apply(asr_our, 1, which.max)
node_states <- factor(node_states, levels = 1:13, labels = c(
  "California", "Indiana", "Kansas", "Maryland", "Michigan", "Minnesota", "Montana",
  "New_Mexico", "North_Carolina", "Ohio", "Oklahoma", "Texas", "Wisconsin"))
all_states <- c(ace_phy$state, node_states)

# Create a data frame for plotting
plot_data <- data.frame(
  node = c(1:length(ace_phy$state), 1:ace_phy$Nnode + length(ace_phy$state)),
  descendant = all_states,
  is_tip = c(rep(TRUE, length(ace_phy$state)), rep(FALSE, ace_phy$Nnode))
)

# Add parent information
plot_data$parent <- sapply(plot_data$node, function(n) {
  if (n == length(ace_phy$state) + 1) return(NA)  # root node
  parent <- ace_phy$edge[ace_phy$edge[,2] == n, 1]
  if (length(parent) == 0) return(NA)  # tip
  return(parent)
})

node_to_node <- plot_data[!plot_data$is_tip & !is.na(plot_data$parent),]
node_to_node$ancestor <- plot_data$descendant[match(node_to_node$parent, plot_data$node)]

# 2. Now get node to tip transitions (new addition)
node_to_tip <- plot_data[plot_data$is_tip & !is.na(plot_data$parent),]
node_to_tip$ancestor <- plot_data$descendant[match(node_to_tip$parent, plot_data$node)]

# 3. Combine both types of transitions
alluvial_data <- rbind(node_to_node, node_to_tip)

# 4. Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data, cols = c(descendant, ancestor),
                              names_to = "generation", values_to = "location")

unique_location <- unique(alluvial_long$location)
lo <- c(
  "California", "Indiana", "Kansas", "Maryland", "Michigan", "Minnesota", "Montana",
  "New_Mexico", "North_Carolina", "Ohio", "Oklahoma", "Texas", "Wisconsin")

alluvial_long$location[alluvial_long$location == lo[1]] <- 1
alluvial_long$location[alluvial_long$location == lo[2]] <- 2
alluvial_long$location[alluvial_long$location == lo[3]] <- 3
alluvial_long$location[alluvial_long$location == lo[4]] <- 4
alluvial_long$location[alluvial_long$location == lo[5]] <- 5
alluvial_long$location[alluvial_long$location == lo[6]] <- 6
alluvial_long$location[alluvial_long$location == lo[7]] <- 7
alluvial_long$location[alluvial_long$location == lo[8]] <- 8
alluvial_long$location[alluvial_long$location == lo[9]] <- 9
alluvial_long$location[alluvial_long$location == lo[10]] <- 10
alluvial_long$location[alluvial_long$location == lo[11]] <- 11
alluvial_long$location[alluvial_long$location == lo[12]] <- 12
alluvial_long$location[alluvial_long$location == lo[13]] <- 13


# Create the alluvial plot
# Create the alluvial plot without the alpha aesthetic
alluvial_our4 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, 
                            fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Location") +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  ggtitle("D") +
  scale_fill_manual(values = c(
    "#E41A1C", "#66C2A5", "#FC8D62", "#8DA0CB", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#FFFF33", "#F781BF", "#999999","#1B9E77" ,"#D95F02"),
    labels=c(lo[1], lo[10], lo[11], lo[12],lo[13], lo[2], lo[3], lo[4], lo[5], lo[6],lo[7], lo[8], lo[9])
  ) +
  theme(text = element_text(size = 15, family = "serif"), plot.title = element_text(size= 15),legend.position = "none")


ggarrange(p3,p4,alluvial_our3,alluvial_our4,heights = c(1.75,1), nrow=2,ncol = 2, common.legend = TRUE,legend = "bottom")


ggsave(file="fig4_updated.pdf", width = 210, height = 297, units = "mm")
