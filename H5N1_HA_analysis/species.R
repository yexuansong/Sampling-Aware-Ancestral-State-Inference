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
  gsub("domestic-mammal", "mammal", .)

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
ans_ace<-ace(ace_phy$state, ace_phy,type = "discrete",method = "ML", model="SYM")

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

# plot the phylogeny and color the tips based on species

# adding more features to the phylogeny
ace_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-01")
ace_species <- ace_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48, 2023.545, 2023.61, 2023.675, 2023.74, 2023.805, 2023.87, 2023.935,
               2024.0, 2024.065, 2024.13, 2024.195, 2024.26), 
    labels = c("04/2023","05/2023","06/2023","07/2023","08/2023","09/2023","10/2023","11/2023","12/2023",
               "01/2024","02/2024","03/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("") +labs(color = "Species")+
  theme(text = element_text(size = 12.0,family = "serif",hjust = 1),plot.title = element_text(size=12),
        axis.text.x = element_text(angle = 45),legend.title = element_text(size = 12, face = "bold"))

p1 <- inset(ace_species, ace_pie,width = 0.04,height = 0.05,hjust=0.003)
####################
# ASR using Stochastic Character Mapping (Phytool)

# fit the model
simmap_phy <- multi2di(tree_modified)
simmap_phy$tip.label = ace_phy$state

lab <- simmap_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab
smap_model <- fitMk(simmap_phy,numeric_vector,model="ARD",pi="fitzjohn")

# do stochastic character mapping
phy_smap<-simmap(smap_model,nsim=100)

summary <- summary(phy_smap)

# extract the state probabilities

simmap_lik <- as.data.frame(summary$ace)
simmap_lik <- head(simmap_lik,ace_phy$Nnode)

simmap_lik$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
simmap_pie <- nodepie(simmap_lik,cols=1:4)

simmap_species <- ggtree(ace_phy,mrsd = "2024-04-06",options(ignore.negative.edge=TRUE))
simmap_species <- simmap_species %<+% tree_info + geom_tippoint(aes(color=species))+
  theme_tree2()+ggtitle("Stochastic Character Mapping") +
  theme(text = element_text(size = 10.0,family = "serif"),plot.title = element_text(size=12))+labs(x = "Time")


p2 <- inset(simmap_species, simmap_pie,width = 0.025,height = 0.025,hjust=0.005)
p2
# ####################
# ASR using our method that accounts for sampling bias


our_phy <- multi2di(tree_modified)
lab <- simmap_phy$tip.label
numeric_vector <- as.numeric(factor(lab))
names(numeric_vector) <- lab

our_phy$tip.state <- numeric_vector

phy <- our_phy
# phy$edge.length <- phy$edge.length*10
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


#estimate transition rates using ace
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)

pars <- c(21.1,21.1,21.1,21.1,
          6.7,6.7,6.7,6.7,
          10,10,10,10,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1)

mm<-qij_matrix

asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:4)

our_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-01")
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48, 2023.545, 2023.61, 2023.675, 2023.74, 2023.805, 2023.87, 2023.935,
               2024.0, 2024.065, 2024.13, 2024.195, 2024.26), 
    labels = c("04/2023","05/2023","06/2023","07/2023","08/2023","09/2023","10/2023","11/2023","12/2023",
               "01/2024","02/2024","03/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("") +labs(color = "Species")+
  theme(text = element_text(size = 12.0,family = "serif",hjust = 1),plot.title = element_text(size=12),
       axis.text.x = element_text(angle = 45),legend.title = element_text(size = 12, face = "bold"))

p5 <- inset(our_species, our_pie,width = 0.04,height = 0.05,hjust=0.005)

asr_our <- our_phisse_back_forward_ailene(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)
asr_our <- apply(asr_our, 1, which.max)
asr_our <- factor(asr_our, levels = 1:4, labels = c("cattle", "mammal","poultry","wild bird"))

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

# Create a data frame for alluvial plot
alluvial_data <- plot_data[!plot_data$is_tip,]
alluvial_data <- alluvial_data[!is.na(alluvial_data$parent),]
alluvial_data$parent_state <- plot_data$state[match(alluvial_data$parent, plot_data$node)]

# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data, cols = c(state, parent_state), 
                              names_to = "generation", values_to = "location")

# Create the alluvial plot
alluvial_our5 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Species") +
  ggtitle("A")+
  scale_fill_discrete(labels = c("cattle", "mammal","poultry","wild bird"))+
  scale_x_discrete(labels = c("Ancestor", "Descendent"))+
  theme(text = element_text(size = 12,family = "serif"),plot.title = element_text(size=12),
        legend.title = element_text(size = 12, face = "bold"))



#######
# now assume wild bird samples 10 times less than other species

pars <- c(21.1,21.1,21.1,21.1,
          6.7,6.7,6.7,6.7,
          10,10,10,1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1)

mm<-qij_matrix
mm<-mm/4
mm[mm==0] <- 1

mm[4,]<- mm[4,]/2
mm[,4]<- mm[,4]*2

asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:4)

our_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-06")
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48, 2023.545, 2023.61, 2023.675, 2023.74, 2023.805, 2023.87, 2023.935,
               2024.0, 2024.065, 2024.13, 2024.195, 2024.26), 
    labels = c("04/2023","05/2023","06/2023","07/2023","08/2023","09/2023","10/2023","11/2023","12/2023",
               "01/2024","02/2024","03/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+ labs(color = "Species")+
  theme_tree2()+ggtitle("") + geom_vline(xintercept = 2024.025, color = "red", linetype = "dashed", size = 0.5)+
  theme(text = element_text(size = 12.0,family = "serif",hjust = 1),plot.title = element_text(size=12), 
        axis.text.x = element_text(angle = 45),legend.title = element_text(size = 12, face = "bold"))



p6 <- inset(our_species, our_pie,width = 0.04,height = 0.05,hjust=0.005)
p6
#######
# now assuming wild bird samples 100 times less

pars <- c(21.1,21.1,21.1,21.1,
          6.7,6.7,6.7,6.7,
          10,10,10,0.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1,
          .1,.1,.1)

mm<-qij_matrix
mm<-mm/4
mm[mm==0] <- 1

mm[4,]<- mm[4,]/5
mm[,4]<- mm[,4]*5

asr_our <- saasi(pars,our_phy,mm)
colnames(asr_our) = c(1,2,3,4)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:4)

our_species <- ggtree(ace_phy,options(ignore.negative.edge=TRUE),mrsd = "2024-04-06")
our_species <- our_species %<+% tree_info + geom_tippoint(aes(color=species),size=2)+
  scale_x_continuous(
    breaks = c(2023.48, 2023.545, 2023.61, 2023.675, 2023.74, 2023.805, 2023.87, 2023.935,
               2024.0, 2024.065, 2024.13, 2024.195, 2024.26), 
    labels = c("04/2023","05/2023","06/2023","07/2023","08/2023","09/2023","10/2023","11/2023","12/2023",
               "01/2024","02/2024","03/2024","04/2024"),
  )+scale_color_discrete(labels = c("wild-bird" = "wild bird"))+
  theme_tree2()+ggtitle("") +labs(color = "Species")+
  geom_vline(xintercept = 2024.025, color = "red", linetype = "dashed", size = 0.5)+
  geom_vline(xintercept = 2024.08, color = "red", linetype = "dashed", size = 0.5)+
  theme(text = element_text(size = 12.0,family = "serif",hjust = 1),plot.title = element_text(size=12), 
        axis.text.x = element_text(angle = 45),legend.title = element_text(size = 12, face = "bold"))

p7 <- inset(our_species, our_pie,width = 0.04,height = 0.05,hjust=0.005)
p7


asr_our <- saasi(pars,our_phy,mm)

colnames(asr_our) = c(1,2,3,4)
asr_our <- apply(asr_our, 1, which.max)
asr_our <- factor(asr_our, levels = 1:4, labels = c("cattle", "mammal","poultry","wild bird"))

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

# Create a data frame for alluvial plot
alluvial_data <- plot_data[!plot_data$is_tip,]
alluvial_data <- alluvial_data[!is.na(alluvial_data$parent),]
alluvial_data$parent_state <- plot_data$state[match(alluvial_data$parent, plot_data$node)]

# Reshape data for ggalluvial
alluvial_long <- pivot_longer(alluvial_data, cols = c(state, parent_state), 
                              names_to = "generation", values_to = "location")

# Create the alluvial plot
alluvial_our7 <- ggplot(alluvial_long,
                        aes(x = generation, stratum = location, alluvium = node, fill = location)) +
  geom_alluvium() +
  geom_stratum() +
  theme_minimal() +
  labs(x = "", y = "", fill = "Species") +
  ggtitle("B")+
  scale_fill_discrete(labels = c("cattle", "mammal","poultry","wild bird"))+
  scale_x_discrete(labels = c("Ancestor", "Descendent"))+
  theme(text = element_text(size = 12,family = "serif"),plot.title = element_text(size=12))


p10 <- ggarrange(alluvial_our5,alluvial_our7,widths = c(1, 1, 1.5),common.legend = TRUE,legend = "right")


