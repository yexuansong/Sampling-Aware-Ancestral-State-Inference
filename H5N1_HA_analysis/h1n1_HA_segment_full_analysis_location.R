##### Comparing results from ape/BiSSE/our implementations
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


# save the figure 
# ggsave("phy_locations.pdf", phy_locations, width = 16, height = 9, units = "in", dpi = 300)

####################
# Analysis (Ancestral State Reconstruction ASR)
####################

####################
# ASR using ML
# adding tip labels (species only)
ace_phy <- multi2di(tree_modified)
ace_phy$state <- tree_info$location

# ace using symmetric transistion rate model
ans_ace<-ace(ace_phy$state, ace_phy,type = "discrete",method = "ML", model="SYM")

# ace estimate the transition rate
qij_rate <- ans_ace$rates

# ace infers the internal nodes state/species probabilities
ace_node_lik <- as.data.frame(ans_ace$lik.anc)

# now plot the inferred phylogeny 

# insert the likelihood pie chart to the tree
ace_node_lik$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
ace_pie <- nodepie(ace_node_lik,cols=1:13)

# plot the phylogeny and color the tips based on locations
lo <- c( 
        "California", "Indiana", "Kansas","Maryland","Michigan","Minnesota","Montana",
        "New_Mexico","North_Carolina","Ohio","Oklahoma","Texas","Wisconsin")

# adding more features to the phylogeny
ace_species <- ggtree(ace_phy,mrsd = "2024-04-06",options(ignore.negative.edge=TRUE))
ace_species <- ace_species %<+% tree_info + geom_tippoint(aes(color=location))+
  theme_tree2()+ggtitle("Maximum Likelihood (ace) Location") +
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(x = "Time")


p1 <- inset(ace_species, ace_pie,width = 0.025,height = 0.025,hjust=0.005)

# save the result
# ggsave("ace_h5n1_location.pdf", p1, width = 16, height = 9, units = "in", dpi = 300)



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

# save the result
# ggsave("alluvial_ace_location.pdf", alluvial_ace, width = 16, height = 9, units = "in", dpi = 300)

####################
# ASR using Stochastic Character Mapping (Phytool)

# fit the model
simmap_phy <- multi2di(tree_modified)
simmap_phy$tip.label = ace_phy$state

lab <- simmap_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab
smap_model <- fitMk(simmap_phy,numeric_vector,model="SYM",pi="fitzjohn")

# do stochastic character mapping
phy_smap<-simmap(smap_model,nsim=100)

summary <- summary(phy_smap)

# extract the state probabilities

simmap_lik <- as.data.frame(summary$ace)
simmap_lik <- head(simmap_lik,ace_phy$Nnode)

simmap_lik$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
simmap_pie <- nodepie(simmap_lik,cols=1:14)

simmap_species <- ggtree(ace_phy,mrsd = "2024-04-06",options(ignore.negative.edge=TRUE))
simmap_species <- simmap_species %<+% tree_info + geom_tippoint(aes(color=location))+
  theme_tree2()+ggtitle("Stochastic Character Mapping") +
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(x = "Time")

p2 <- inset(simmap_species, simmap_pie,width = 0.025,height = 0.025,hjust=0.005)

# save the result
# ggsave("simmap_species_location.pdf", p2, width = 16, height = 9, units = "in", dpi = 300)

############################
# now plot the alluvial plot
# Extract node states

simmap_lik <- as.data.frame(summary$ace)
simmap_lik <- head(simmap_lik,ace_phy$Nnode)

node_states <- apply(simmap_lik, 1, which.max)
node_states <- factor(node_states, levels = 1:13, labels = c("Arizona", 
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
alluvial_simmap <- ggplot(alluvial_long,
                          aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "Generation", y = "Count", fill = "location") +
  ggtitle("Movement between Species - Simmap") +
  scale_fill_discrete(labels = c(lo[1],lo[10],lo[11],lo[12],lo[14],lo[2],lo[3],lo[4],lo[5],lo[7],lo[8],lo[9]))

# save the result
# ggsave("alluvial_simmap_location.pdf", alluvial_simmap, width = 16, height = 9, units = "in", dpi = 300)



# ####################
# ASR using SAASI 

phy <- our_phy

# estimate speciation and extinction rates using Stadler et al. 2012
branching_time <- extract_tree_times(phy)
X <- as.vector(branching_time$internal_nodes)
X <- X$time
Y <- as.vector(branching_time$tips)
Y <- Y$time


nLL <- function(lambda, mu) {
  # Return negative log likelihood
  # We use negative because mle() minimizes
  -(probability_density(lambda, mu,psi=10, m=phy$Nnode, X=X, Y=Y))
}

fit <- mle(nLL,
           start = list(lambda = 1, mu = 0.1),
           method = "L-BFGS-B",
           lower = c(0.1, 0.05),
           upper = c(100, 100))


# this gives us the speciation rate 21.1 & extinction rate 6.7

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
    breaks = c(2023.48, 2023.545, 2023.61, 2023.675, 2023.74, 2023.805, 2023.87, 2023.935,
               2024.0, 2024.065, 2024.13, 2024.195, 2024.26), 
    labels = c("04/2023","05/2023","06/2023","07/2023","08/2023","09/2023","10/2023","11/2023","12/2023",
               "01/2024","02/2024","03/2024","04/2024"),
  )+  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
    "#FFD92F", "#E5C494", "#B3B3B3"
  ))+
  theme_tree2()+ggtitle("") +labs(color = "Location")+
  geom_vline(xintercept = 2023.97, color = "red", linetype = "dashed", size = 0.5)+
  theme(text = element_text(size = 12.0,family = "serif",hjust = 1),plot.title = element_text(size=12), 
        axis.text.x = element_text(angle = 45),legend.title = element_text(size = 12, face = "bold"))


p3 <- inset(our_species, our_pie,width = 0.04,height = 0.05,hjust=0.005)

# save the figure
# ggsave("our_h5n1_location_texas_more.pdf", p3,  width = 210, height = 297, units = "mm")

############################
# now plot the alluvial plot
# Extract node states

asr_our <- saasi(pars,our_phy,mm)

node_states <- apply(asr_our, 1, which.max)
node_states <- factor(node_states, levels = 1:13, 
                      labels = c("California", "Indiana", "Kansas","Maryland",
                                 "Michigan","Minnesota","Montana","New_Mexico",
                                 "North_Carolina","Ohio","Oklahoma","Texas",
                                 "Wisconsin"))

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

# Create a data frame for alluvial plot
alluvial_data <- plot_data[!plot_data$is_tip,]
alluvial_data <- alluvial_data[!is.na(alluvial_data$parent),]
alluvial_data$ancestor <- plot_data$descendant[match(alluvial_data$parent, plot_data$node)]

# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data, cols = c(descendant, ancestor), 
                              names_to = "generation", values_to = "location")

unique_location <- unique(alluvial_long$location)
lo <- c(
        "California", "Indiana", "Kansas","Maryland","Michigan","Minnesota","Montana",
        "New_Mexico","North_Carolina","Ohio","Oklahoma","Texas","Wisconsin")


# Create the alluvial plot
alluvial_our <- ggplot(alluvial_long,
                       aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Location") +
  scale_x_discrete(labels = c("Ancestor", "Descendent"))+
  ggtitle("A")+
  scale_fill_manual(values = c(
    "#E41A1C", "#66C2A5", "#FC8D62", "#8DA0CB", "#377EB8","#4DAF4A","#984EA3",
    "#FF7F00", "#FFFF33","#F781BF","#999999"),  
    labels=c(lo[1],lo[10],lo[11],lo[12],lo[2],lo[3],lo[4],lo[5],lo[6],lo[8],lo[9]))+
theme(text = element_text(size = 12,family = "serif"),plot.title = element_text(size=12))


#######################

# now using SAASI assuming texas oversampled 5 times

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
mm<-mm/1.4
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
    breaks = c(2023.48, 2023.545, 2023.61, 2023.675, 2023.74, 2023.805, 2023.87, 2023.935,
               2024.0, 2024.065, 2024.13, 2024.195, 2024.26), 
    labels = c("04/2023","05/2023","06/2023","07/2023","08/2023","09/2023","10/2023","11/2023","12/2023",
               "01/2024","02/2024","03/2024","04/2024"),
  )+  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", 
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", 
    "#FFD92F", "#E5C494", "#B3B3B3"
  ))+
  theme_tree2()+ggtitle("") +labs(color = "Location")+
  geom_vline(xintercept = 2023.97, color = "red", linetype = "dashed", size = 0.5)+
  theme(text = element_text(size = 12.0,family = "serif",hjust = 1),plot.title = element_text(size=12), 
        axis.text.x = element_text(angle = 45),legend.title = element_text(size = 12, face = "bold"))


p4 <- inset(our_species, our_pie,width = 0.04,height = 0.05,hjust=0.005)

# save the result 
# ggsave("our_h5n1_location_equal_sampling.pdf", p4, width = 210, height = 297, units = "mm")

############################
# now plot the alluvial plot
# Extract node states

asr_our <- saasi(pars,our_phy,mm)

node_states <- apply(asr_our, 1, which.max)
node_states <- factor(node_states, levels = 1:13, labels = c(
  "California", "Indiana", "Kansas","Maryland","Michigan","Minnesota","Montana",
  "New_Mexico","North_Carolina","Ohio","Oklahoma","Texas","Wisconsin"))

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

# Create a data frame for alluvial plot
alluvial_data <- plot_data[!plot_data$is_tip,]
alluvial_data <- alluvial_data[!is.na(alluvial_data$parent),]
alluvial_data$ancestor <- plot_data$descendant[match(alluvial_data$parent, plot_data$node)]

# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data, cols = c(descendant, ancestor), 
                              names_to = "generation", values_to = "location")

unique_location <- unique(alluvial_long$location)
lo <- c(
  "California", "Indiana", "Kansas","Maryland","Michigan","Minnesota","Montana",
  "New_Mexico","North_Carolina","Ohio","Oklahoma","Texas","Wisconsin")

# Create the alluvial plot
alluvial_our2 <- ggplot(alluvial_long,
                       aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Location") +
  scale_x_discrete(labels = c("Ancestor", "Descendent"))+
  ggtitle("B")+
  scale_fill_manual(values = c(
    "#E41A1C", "#66C2A5", "#FC8D62", "#8DA0CB", "#377EB8","#4DAF4A","#984EA3",
    "#FF7F00", "#FFFF33","#A65628","#F781BF","#999999"),
    labels=c(lo[1],lo[10],lo[11],lo[12],lo[2],lo[3],lo[4],lo[5],lo[6],lo[7],lo[8],lo[9]))+
 theme(text = element_text(size = 12,family = "serif"),plot.title = element_text(size=12))

alluvial_our3 <- alluvial_our + theme(legend.position = "none")
alluvial_our4 <- alluvial_our2 + theme(legend.position = "none")
pall <- ggarrange(alluvial_our3,alluvial_our2,widths = c(1, 1, 1.5))

# comparing the two alluvial plots, and save them
# ggsave("h5n1_location_all.pdf", pall,width = 210, height = 297, units = "mm")




