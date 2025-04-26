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
library(ggalluvial)
library(rsvg)
library(rphylopic)
library(cowplot)

source(file.path("h5n1_phy_data_cleaning.R"))


####################
# unmodified phylogeny
unmod <- ggtree(tree)

# modified phylogeny

tree_modified <- di2multi(tree, tol = 1e-08)

drop_tip <- c(1,104)
tree_modified <- drop.tip(tree_modified,drop_tip)
tree_info <- tree_info[-c(1,104),]


phy_species <- ggtree(tree_modified)
phy_species <- phy_species %<+% tree_info + geom_tippoint(aes(color=species))

# drop the human tips, combine mammal and wild mammal
tips_to_drop <- grep("human", tree_modified$tip.label, value = TRUE, ignore.case = TRUE)
tree_modified <- drop.tip(tree_modified, tips_to_drop)

tree_info <- tree_info %>% 
  filter(species != "human")
tree_info$species <- tree_info$species %>%
  gsub("domestic-mammal", "mammal",.)

phy_species <- ggtree(tree_modified)
phy_species <- phy_species %<+% tree_info + geom_tippoint(aes(color=species))

####################
# Analysis (Ancestral State Reconstruction ASR)
####################

####################
# ASR using ML
# adding tip labels (species only)


ace_phy <- multi2di(tree_modified)
ace_phy$state <- tree_info$species

# ace_phy$edge.length <- ace_phy$edge.length*10

ans_ace<-ace(ace_phy$state, ace_phy,type = "discrete",method = "ML", model="ER")

# ace estimate the transition rate
qij_rate <- ans_ace$rates
ans_ace$index.matrix


replace_matrix_with_vector <- function(matrix, vector) {
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      matrix[i,j] <- vector[matrix[i,j]]
    }
  }
  return(matrix)
}

ace_node_lik <- as.data.frame(ans_ace$lik.anc)

# now plot the inferred phylogeny 

# insert the likelihood pie chart to the tree
ace_node_lik$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
ace_pie <- nodepie(ace_node_lik,cols=1:4)


ace_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-01")
ace_species <- ace_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("A: ace") +labs(color = "Species")+
  #geom_hilight(node=120, fill="purple",alpha=0.05) +
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15),
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))


p1 <- inset(ace_species, ace_pie,width = 0.07,height = 0.07,hjust=0.003)


# ####################
# SAASI - equal sampling , WB adjusted by a factor of 2

simmap_phy <- multi2di(tree_modified)
simmap_phy$tip.label = ace_phy$state

lab <- simmap_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab

our_phy <- multi2di(tree_modified)
lab <- simmap_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab

our_phy$tip.state <- numeric_vector

phy <- our_phy


qij_rate <- 2
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)
qij_matrix

pars <- c(21.1,21.1,21.1,21.1,
          6.7,6.7,6.7,6.7,
          10,10,10,10,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1)
qij_matrix2 <- replace_matrix_with_vector(ans_ace$index.matrix,ans_ace$rates)

mm<-qij_matrix
mm[,4] = mm[,4]/2
mm[4,] = mm[4,]*2
asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:4)

our_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-01")
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("A: Equal sampling \nAdjusting WB transition rates \nby a factor of 2") +labs(color = "Species")+
  #geom_hilight(node=120, fill="purple",alpha=0.05) +
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15),
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))

p2 <- inset(our_species, our_pie,width = 0.07,height = 0.07,hjust=0.005)

asr_our <- saasi(pars, our_phy, mm)
colnames(asr_our) = c(1, 2, 3, 4)
asr_our <- apply(asr_our, 1, which.max)
all_states <- c(ace_phy$state, asr_our)
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
# For internal-to-internal node transitions
alluvial_data_internal <- plot_data[!plot_data$is_tip, ]
alluvial_data_internal <- alluvial_data_internal[!is.na(alluvial_data_internal$parent), ]
alluvial_data_internal$parent_state <- plot_data$state[match(alluvial_data_internal$parent, plot_data$node)]
# For internal node-to-tip transitions
alluvial_data_tips <- plot_data[plot_data$is_tip, ]
alluvial_data_tips <- alluvial_data_tips[!is.na(alluvial_data_tips$parent), ]
alluvial_data_tips$parent_state <- plot_data$state[match(alluvial_data_tips$parent, plot_data$node)]
# Combine both sets
alluvial_data_combined <- rbind(alluvial_data_internal, alluvial_data_tips)
# Add a transition type column
alluvial_data_combined$transition_type <- ifelse(alluvial_data_combined$is_tip, "To Tip", "To Internal")
# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data_combined, 
                              cols = c(state, parent_state), 
                              names_to = "generation", 
                              values_to = "location")
alluvial_long$location[alluvial_long$location == "cattle"] <- 1
alluvial_long$location[alluvial_long$location == "mammal"] <- 2
alluvial_long$location[alluvial_long$location == "poultry"] <- 3
alluvial_long$location[alluvial_long$location == "wild-bird"] <- 4
# Create the alluvial plot
alluvial_our2 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Species") +
  ggtitle("D") +
  scale_fill_discrete(labels = c("cattle", "mammal", "poultry", "wild bird")) +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  theme(text = element_text(size = 15, family = "serif",), plot.title = element_text(size=15),
        legend.position = "none")

# Equal sampling, equal transition
our_phy <- multi2di(tree_modified)
lab <- simmap_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab

our_phy$tip.state <- numeric_vector

phy <- our_phy


qij_rate <- 2
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)
qij_matrix

pars <- c(21.1,21.1,21.1,21.1,
          6.7,6.7,6.7,6.7,
          10,10,10,10,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1)
qij_matrix2 <- replace_matrix_with_vector(ans_ace$index.matrix,ans_ace$rates)

mm<-qij_matrix
asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:4)

our_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-01")
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("A: Equal sampling \nEqual transition rates") +labs(color = "Species")+
  #geom_hilight(node=120, fill="purple",alpha=0.05) +
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15),
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))

p3 <- inset(our_species, our_pie,width = 0.07,height = 0.07,hjust=0.005)

asr_our <- saasi(pars, our_phy, mm)
colnames(asr_our) = c(1, 2, 3, 4)
asr_our <- apply(asr_our, 1, which.max)
all_states <- c(ace_phy$state, asr_our)
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
# For internal-to-internal node transitions
alluvial_data_internal <- plot_data[!plot_data$is_tip, ]
alluvial_data_internal <- alluvial_data_internal[!is.na(alluvial_data_internal$parent), ]
alluvial_data_internal$parent_state <- plot_data$state[match(alluvial_data_internal$parent, plot_data$node)]
# For internal node-to-tip transitions
alluvial_data_tips <- plot_data[plot_data$is_tip, ]
alluvial_data_tips <- alluvial_data_tips[!is.na(alluvial_data_tips$parent), ]
alluvial_data_tips$parent_state <- plot_data$state[match(alluvial_data_tips$parent, plot_data$node)]
# Combine both sets
alluvial_data_combined <- rbind(alluvial_data_internal, alluvial_data_tips)
# Add a transition type column
alluvial_data_combined$transition_type <- ifelse(alluvial_data_combined$is_tip, "To Tip", "To Internal")
# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data_combined, 
                              cols = c(state, parent_state), 
                              names_to = "generation", 
                              values_to = "location")
alluvial_long$location[alluvial_long$location == "cattle"] <- 1
alluvial_long$location[alluvial_long$location == "mammal"] <- 2
alluvial_long$location[alluvial_long$location == "poultry"] <- 3
alluvial_long$location[alluvial_long$location == "wild-bird"] <- 4
# Create the alluvial plot
alluvial_our3 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Species") +
  ggtitle("D") +
  scale_fill_discrete(labels = c("cattle", "mammal", "poultry", "wild bird")) +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  theme(text = element_text(size = 15, family = "serif",), plot.title = element_text(size=15),
        legend.position = "none")




# WB at 1/10 sampling, WB adjusted by a factor of 2

qij_rate <- 2
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)

pars <- c(21.1,21.1,21.1,21.1,
          6.7,6.7,6.7,6.7,
          10,10,10,1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1)
qij_matrix2 <- replace_matrix_with_vector(ans_ace$index.matrix,ans_ace$rates)

mm<-qij_matrix
mm[,4] = mm[,4]/2
mm[4,] = mm[4,]*2
asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:4)

our_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-01")
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("B: WB at 1/10 sampling \nAdjusting WB transition rates \nby a factor of 2") +labs(color = "Species")+
  #geom_text2(aes(subset = (node == 127), label = "Key"), hjust = 1.5,vjust=2)+
  geom_vline(xintercept = 2024.055, color = "red", linetype = "dashed", size = 0.5)+
  geom_hilight(node=120, fill="purple",alpha=0.05) +
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15),
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))
our_species
p4 <- inset(our_species, our_pie,width = 0.07,height = 0.07,hjust=0.005)
p4
asr_our <- saasi(pars, our_phy, mm)
colnames(asr_our) = c(1, 2, 3, 4)
asr_our <- apply(asr_our, 1, which.max)
all_states <- c(ace_phy$state, asr_our)
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
# For internal-to-internal node transitions
alluvial_data_internal <- plot_data[!plot_data$is_tip, ]
alluvial_data_internal <- alluvial_data_internal[!is.na(alluvial_data_internal$parent), ]
alluvial_data_internal$parent_state <- plot_data$state[match(alluvial_data_internal$parent, plot_data$node)]
# For internal node-to-tip transitions
alluvial_data_tips <- plot_data[plot_data$is_tip, ]
alluvial_data_tips <- alluvial_data_tips[!is.na(alluvial_data_tips$parent), ]
alluvial_data_tips$parent_state <- plot_data$state[match(alluvial_data_tips$parent, plot_data$node)]
# Combine both sets
alluvial_data_combined <- rbind(alluvial_data_internal, alluvial_data_tips)
# Add a transition type column
alluvial_data_combined$transition_type <- ifelse(alluvial_data_combined$is_tip, "To Tip", "To Internal")
# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data_combined, 
                              cols = c(state, parent_state), 
                              names_to = "generation", 
                              values_to = "location")
alluvial_long$location[alluvial_long$location == "cattle"] <- 1
alluvial_long$location[alluvial_long$location == "mammal"] <- 2
alluvial_long$location[alluvial_long$location == "poultry"] <- 3
alluvial_long$location[alluvial_long$location == "wild-bird"] <- 4
# Create the alluvial plot
alluvial_our4 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Species") +
  ggtitle("E") +
  scale_fill_discrete(labels = c("cattle", "mammal", "poultry", "wild bird")) +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  theme(text = element_text(size = 15, family = "serif",), plot.title = element_text(size= 15),
        legend.position = "none")




# WB at 1/10 sampling, equal transition rates
qij_rate <- 2
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)

pars <- c(21.1,21.1,21.1,21.1,
          6.7,6.7,6.7,6.7,
          10,10,10,1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1)
qij_matrix2 <- replace_matrix_with_vector(ans_ace$index.matrix,ans_ace$rates)

mm<-qij_matrix
asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:4)

our_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-01")
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("B: WB at 1/10 sampling \nEqual transition rates") +labs(color = "Species")+
  # geom_text2(aes(subset = (node == 127), label = "Key"), hjust = 1.5,vjust=2)+
  #geom_hilight(node=120, fill="purple",alpha=0.05) +
  geom_vline(xintercept = 2024.055, color = "red", linetype = "dashed", size = 0.5)+
  geom_hilight(node=120, fill="purple",alpha=0.05) +
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15),
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))
our_species
p5 <- inset(our_species, our_pie,width = 0.07,height = 0.07,hjust=0.005)
p5
asr_our <- saasi(pars, our_phy, mm)
colnames(asr_our) = c(1, 2, 3, 4)
asr_our <- apply(asr_our, 1, which.max)
all_states <- c(ace_phy$state, asr_our)
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
# For internal-to-internal node transitions
alluvial_data_internal <- plot_data[!plot_data$is_tip, ]
alluvial_data_internal <- alluvial_data_internal[!is.na(alluvial_data_internal$parent), ]
alluvial_data_internal$parent_state <- plot_data$state[match(alluvial_data_internal$parent, plot_data$node)]
# For internal node-to-tip transitions
alluvial_data_tips <- plot_data[plot_data$is_tip, ]
alluvial_data_tips <- alluvial_data_tips[!is.na(alluvial_data_tips$parent), ]
alluvial_data_tips$parent_state <- plot_data$state[match(alluvial_data_tips$parent, plot_data$node)]
# Combine both sets
alluvial_data_combined <- rbind(alluvial_data_internal, alluvial_data_tips)
# Add a transition type column
alluvial_data_combined$transition_type <- ifelse(alluvial_data_combined$is_tip, "To Tip", "To Internal")
# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data_combined, 
                              cols = c(state, parent_state), 
                              names_to = "generation", 
                              values_to = "location")
alluvial_long$location[alluvial_long$location == "cattle"] <- 1
alluvial_long$location[alluvial_long$location == "mammal"] <- 2
alluvial_long$location[alluvial_long$location == "poultry"] <- 3
alluvial_long$location[alluvial_long$location == "wild-bird"] <- 4
# Create the alluvial plot
alluvial_our5 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Species") +
  ggtitle("E") +
  scale_fill_discrete(labels = c("cattle", "mammal", "poultry", "wild bird")) +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  theme(text = element_text(size = 15, family = "serif",), plot.title = element_text(size= 15),
        legend.position = "none")



# WB at 1/100 sampling, WB adjusted by a factor of 2

qij_rate <- 2
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)


pars <- c(21.1,21.1,21.1,21.1,
          6.7,6.7,6.7,6.7,
          10,10,10,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1)
qij_matrix2 <- replace_matrix_with_vector(ans_ace$index.matrix,ans_ace$rates)

mm<-qij_matrix
mm[,4] = mm[,4]/2
mm[4,] = mm[4,]*2
asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:4)

our_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-01")
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("C: WB at 1/100 sampling \nAdjusting WB transition rates \nby a factor of 2") +labs(color = "Species")+
  #geom_hilight(node=120, fill="purple",alpha=0.05) +
  # geom_text2(aes(subset = (node == 127), label = "Key"), hjust = 1.5)+
  geom_vline(xintercept = 2024.105, color = "red", linetype = "dashed", size = 0.5)+
  geom_hilight(node=120, fill="purple",alpha=0.05) +
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15),
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))
our_species
p6 <- inset(our_species, our_pie,width = 0.07,height = 0.07,hjust=0.005)
p6
asr_our <- saasi(pars, our_phy, mm)
colnames(asr_our) = c(1, 2, 3, 4)
asr_our <- apply(asr_our, 1, which.max)
all_states <- c(ace_phy$state, asr_our)
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
# For internal-to-internal node transitions
alluvial_data_internal <- plot_data[!plot_data$is_tip, ]
alluvial_data_internal <- alluvial_data_internal[!is.na(alluvial_data_internal$parent), ]
alluvial_data_internal$parent_state <- plot_data$state[match(alluvial_data_internal$parent, plot_data$node)]
# For internal node-to-tip transitions
alluvial_data_tips <- plot_data[plot_data$is_tip, ]
alluvial_data_tips <- alluvial_data_tips[!is.na(alluvial_data_tips$parent), ]
alluvial_data_tips$parent_state <- plot_data$state[match(alluvial_data_tips$parent, plot_data$node)]
# Combine both sets
alluvial_data_combined <- rbind(alluvial_data_internal, alluvial_data_tips)
# Add a transition type column
alluvial_data_combined$transition_type <- ifelse(alluvial_data_combined$is_tip, "To Tip", "To Internal")
# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data_combined, 
                              cols = c(state, parent_state), 
                              names_to = "generation", 
                              values_to = "location")
alluvial_long$location[alluvial_long$location == "cattle"] <- 1
alluvial_long$location[alluvial_long$location == "mammal"] <- 2
alluvial_long$location[alluvial_long$location == "poultry"] <- 3
alluvial_long$location[alluvial_long$location == "wild-bird"] <- 4
# Create the alluvial plot
alluvial_our6 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Species") +
  ggtitle("F") +
  scale_fill_discrete(labels = c("cattle", "mammal", "poultry", "wild bird")) +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  theme(text = element_text(size = 15, family = "serif",), plot.title = element_text(size= 15),
        legend.position = "none")




# WB at 1/100 sampling, equal transition rates

qij_rate <- 2
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)
pars <- c(21.1,21.1,21.1,21.1,
          6.7,6.7,6.7,6.7,
          10,10,10,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1)
qij_matrix2 <- replace_matrix_with_vector(ans_ace$index.matrix,ans_ace$rates)

mm<-qij_matrix
asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:4)

our_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-01")
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48,  2023.61,  2023.74, 2023.87,
               2024.0,  2024.13,  2024.26), 
    labels = c("04/2023","06/2023","08/2023","10/2023","12/2023",
               "02/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("C: WB at 1/100 sampling \nEqual transition rates") +labs(color = "Species")+
  #geom_hilight(node=120, fill="purple",alpha=0.05) +
  # geom_text2(aes(subset = (node == 127), label = "Key"), hjust = 1.5)+
  geom_vline(xintercept = 2024.105, color = "red", linetype = "dashed", size = 0.5)+
  geom_hilight(node=120, fill="purple",alpha=0.05) +
  theme(text = element_text(size = 15,family = "serif",hjust = 1),plot.title = element_text(size=15),
        axis.text.x = element_text(angle = 60),legend.title = element_text(size = 15, face = "bold"))
our_species
p7 <- inset(our_species, our_pie,width = 0.07,height = 0.07,hjust=0.005)
p7
asr_our <- saasi(pars, our_phy, mm)
colnames(asr_our) = c(1, 2, 3, 4)
asr_our <- apply(asr_our, 1, which.max)
all_states <- c(ace_phy$state, asr_our)
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
# For internal-to-internal node transitions
alluvial_data_internal <- plot_data[!plot_data$is_tip, ]
alluvial_data_internal <- alluvial_data_internal[!is.na(alluvial_data_internal$parent), ]
alluvial_data_internal$parent_state <- plot_data$state[match(alluvial_data_internal$parent, plot_data$node)]
# For internal node-to-tip transitions
alluvial_data_tips <- plot_data[plot_data$is_tip, ]
alluvial_data_tips <- alluvial_data_tips[!is.na(alluvial_data_tips$parent), ]
alluvial_data_tips$parent_state <- plot_data$state[match(alluvial_data_tips$parent, plot_data$node)]
# Combine both sets
alluvial_data_combined <- rbind(alluvial_data_internal, alluvial_data_tips)
# Add a transition type column
alluvial_data_combined$transition_type <- ifelse(alluvial_data_combined$is_tip, "To Tip", "To Internal")
# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data_combined, 
                              cols = c(state, parent_state), 
                              names_to = "generation", 
                              values_to = "location")
alluvial_long$location[alluvial_long$location == "cattle"] <- 1
alluvial_long$location[alluvial_long$location == "mammal"] <- 2
alluvial_long$location[alluvial_long$location == "poultry"] <- 3
alluvial_long$location[alluvial_long$location == "wild-bird"] <- 4
# Create the alluvial plot
alluvial_our7 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Species") +
  ggtitle("F") +
  scale_fill_discrete(labels = c("cattle", "mammal", "poultry", "wild bird")) +
  scale_x_discrete(labels = c("Ancestor", "Descendent")) +
  theme(text = element_text(size = 15, family = "serif",), plot.title = element_text(size= 15),
        legend.position = "none")


ggarrange(p3,p5,p7,alluvial_our3,alluvial_our5,alluvial_our7,nrow=2,ncol=3,heights = c(1.75,1),common.legend = TRUE,legend = "bottom")
ggsave(file="fig3_updated.pdf", width = 210, height = 297, units = "mm")

ggarrange(p2,p4,p6,alluvial_our2,alluvial_our4,alluvial_our6,nrow=2,ncol=3,heights = c(1.75,1),common.legend = TRUE,legend = "bottom")
ggsave(file="figs4_updated.pdf", width = 210, height = 297, units = "mm")

cc <- as.tibble(our_phy)  






