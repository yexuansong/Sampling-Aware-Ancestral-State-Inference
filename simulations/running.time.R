library(microbenchmark)
library(patchwork)

pars <- c(10,  10,
          .45, .45,
          .01,  .5,
          .1,
          .1)

k = 2

qij_matrix <- function(k) {
  mat <- matrix(0.05, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2)

#q_matrix[1,2]=0.3
#############################
# simulation
set.seed(1) 


phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
# Ensure the tree has at least 30 tips
while(length(phy$tip.label) < 10000) {
  phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
}
phy <- prune(phy)
h <- history.from.sim.discrete(phy, 1:k)
#plot(h, phy, cex=.7)
plot(h, phy)

true_phy_info <- as_tibble(phy)
phy_data <- c(factor(h$tip.state),factor(h$node.state))
true_phy_info$State <- phy_data
true_phy <- ggtree(phy)
true_phy <- true_phy  %<+% true_phy_info + geom_point(aes(color=State)) +
  ggtitle("True Phylogeny") +
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))
#true_phy
node_depths <- node.depth.edgelength(phy)
tmrca <- max(node_depths)
tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
new_phy <- drop.tip(phy, tips_to_drop)
phy_our <- new_phy
true_phy_info_new <- as_tibble(new_phy)
phy_data <- c(factor(h$tip.state),factor(h$node.state))
true_phy_info_new <- true_phy_info_new %>% mutate(State = phy_data[label])
new_phy$tip.state <- new_phy$tip.state[setdiff(names(new_phy$tip.state), tips_to_drop)]
phy <- new_phy
true_phy_new <- ggtree(phy)
true_phy_new <- true_phy_new  %<+% true_phy_info_new + geom_point(aes(color=State)) +
  ggtitle("True Phylogeny") +
  theme(text = element_text(size = 12,family = "serif",face="bold"),plot.title = element_text(size=12))
true_phy_new


phy$Nnode


# asr_our <- our_phisse_back_forward_ailene(pars,phy,q_matrix)
# 
# asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
# our_pie <- nodepie(asr_our,cols=1:k)



ace_phy <- phy
ace_phy$node.label <- NULL
# Note: Do not have this problem if use earlier version `ape`
# Error in names(obj$ace) <- phy$node.label : 
#   attempt to set an attribute on NULL

ace_phy$tip.state <- ace_phy$tip.state[setdiff(names(ace_phy$tip.state), tips_to_drop)]
asr_ace<-ace(ace_phy$tip.state, ace_phy,type = "discrete", model="SYM")

ace_node_lik <- as.data.frame(asr_ace$lik.anc)
ace_node_lik$node <- 1:new_phy$Nnode + Ntip(new_phy)
#ace_node_lik$node <- 1:phy$Nnode + Ntip(phy)

ace_pie <- nodepie(ace_node_lik,cols=1:k)

p1 <- ggtree(ace_phy)
p1 <- p1 %<+% true_phy_info_new + geom_tippoint(aes(color=State))+
  ggtitle("ASR-ace") +
  theme(text = element_text(size = 12,family = "serif",face="bold"),plot.title = element_text(size=12))
p1 <- inset(p1, ace_pie,width = 0.05,height = 0.05)
p1

# running time tree ~ 100 tips
microbenchmark(
  our_phisse_back_forward_ailene(pars, phy, q_matrix),
  times = 10  # Run multiple times for accuracy
)

asr_our <- our_phisse_back_forward_ailene(pars,phy,q_matrix)
asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:k)

p3 <- ggtree(ace_phy)
p3 <- p3 %<+% true_phy_info_new + geom_tippoint(aes(color=State),size=2)+
  ggtitle("SAAI - Unequal sampling") +
  theme(text = element_text(size = 12,family = "serif",face="bold"),plot.title = element_text(size=12))
p3 <- inset(p3, our_pie,width = 0.1,height = 0.1)
p3





# now need larger tree, 500?
pars <- c(10,  10,
          .45, .45,
          .05,  .5,
          .1,
          .1)

k = 2

qij_matrix <- function(k) {
  mat <- matrix(0.05, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2)

set.seed(1) 

phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
while(length(phy$tip.label) < 10000) {
  phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
}
phy <- prune(phy)
h <- history.from.sim.discrete(phy, 1:k)
node_depths <- node.depth.edgelength(phy)
tmrca <- max(node_depths)
tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
new_phy <- drop.tip(phy, tips_to_drop)
phy_our <- new_phy
true_phy_info_new <- as_tibble(new_phy)
phy_data <- c(factor(h$tip.state),factor(h$node.state))
true_phy_info_new <- true_phy_info_new %>% mutate(State = phy_data[label])
new_phy$tip.state <- new_phy$tip.state[setdiff(names(new_phy$tip.state), tips_to_drop)]
phy <- new_phy


phy$Nnode


# running time tree ~ 200 tips
microbenchmark(
  our_phisse_back_forward_ailene(pars, phy, q_matrix),
  times = 10  # Run multiple times for accuracy
)



# now need larger tree, 1000
pars <- c(10,  10,
          .45, .45,
          .5,  1,
          .1,
          .1)

k = 2

qij_matrix <- function(k) {
  mat <- matrix(0.05, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2)

set.seed(1) 

phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
while(length(phy$tip.label) < 10000) {
  phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
}
phy <- prune(phy)
h <- history.from.sim.discrete(phy, 1:k)
node_depths <- node.depth.edgelength(phy)
tmrca <- max(node_depths)
tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
new_phy <- drop.tip(phy, tips_to_drop)
phy_our <- new_phy
true_phy_info_new <- as_tibble(new_phy)
phy_data <- c(factor(h$tip.state),factor(h$node.state))
true_phy_info_new <- true_phy_info_new %>% mutate(State = phy_data[label])
new_phy$tip.state <- new_phy$tip.state[setdiff(names(new_phy$tip.state), tips_to_drop)]
phy <- new_phy


phy$Nnode


# running time tree ~ 200 tips
microbenchmark(
  our_phisse_back_forward_ailene(pars, phy, q_matrix),
  times = 10  # Run multiple times for accuracy
)





# now need larger tree, 2000
pars <- c(10,  10,
          .45, .45,
          1,  1,
          .1,
          .1)

k = 2

qij_matrix <- function(k) {
  mat <- matrix(0.05, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2)

set.seed(1) 

phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
while(length(phy$tip.label) < 10000) {
  phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
}
phy <- prune(phy)
h <- history.from.sim.discrete(phy, 1:k)
node_depths <- node.depth.edgelength(phy)
tmrca <- max(node_depths)
tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
new_phy <- drop.tip(phy, tips_to_drop)
phy_our <- new_phy
true_phy_info_new <- as_tibble(new_phy)
phy_data <- c(factor(h$tip.state),factor(h$node.state))
true_phy_info_new <- true_phy_info_new %>% mutate(State = phy_data[label])
new_phy$tip.state <- new_phy$tip.state[setdiff(names(new_phy$tip.state), tips_to_drop)]
phy <- new_phy


phy$Nnode


# running time tree ~ 200 tips
microbenchmark(
  our_phisse_back_forward_ailene(pars, phy, q_matrix),
  times = 10  # Run multiple times for accuracy
)




# now need larger tree, 5000?
pars <- c(10,  10,
          .45, .45,
          2,  2,
          .1,
          .1)

k = 2

qij_matrix <- function(k) {
  mat <- matrix(0.05, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2)

set.seed(1) 

phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
while(length(phy$tip.label) < 10000) {
  phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
}
phy <- prune(phy)
h <- history.from.sim.discrete(phy, 1:k)
node_depths <- node.depth.edgelength(phy)
tmrca <- max(node_depths)
tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
new_phy <- drop.tip(phy, tips_to_drop)
phy_our <- new_phy
true_phy_info_new <- as_tibble(new_phy)
phy_data <- c(factor(h$tip.state),factor(h$node.state))
true_phy_info_new <- true_phy_info_new %>% mutate(State = phy_data[label])
new_phy$tip.state <- new_phy$tip.state[setdiff(names(new_phy$tip.state), tips_to_drop)]
phy <- new_phy


phy$Nnode


# running time tree ~ 200 tips
microbenchmark(
  our_phisse_back_forward_ailene(pars, phy, q_matrix),
  times = 10  # Run multiple times for accuracy
)





# now need larger tree, 5000?
pars <- c(10,  10,
          .45, .45,
          8,  8,
          .1,
          .1)

k = 2

qij_matrix <- function(k) {
  mat <- matrix(0.05, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2)

set.seed(1) 

phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
while(length(phy$tip.label) < 10000) {
  phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 20000, 20000)
}
phy <- prune(phy)
h <- history.from.sim.discrete(phy, 1:k)
node_depths <- node.depth.edgelength(phy)
tmrca <- max(node_depths)
tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
new_phy <- drop.tip(phy, tips_to_drop)
phy_our <- new_phy
true_phy_info_new <- as_tibble(new_phy)
phy_data <- c(factor(h$tip.state),factor(h$node.state))
true_phy_info_new <- true_phy_info_new %>% mutate(State = phy_data[label])
new_phy$tip.state <- new_phy$tip.state[setdiff(names(new_phy$tip.state), tips_to_drop)]
phy <- new_phy


phy$Nnode


# running time tree ~ 100,000 tips
microbenchmark(
  our_phisse_back_forward_ailene(pars, phy, q_matrix),
  times = 1  # Run multiple times for accuracy
)

# shows the mean running time
run.time <- data.frame(nodes = c(100,200,1000), time = c(2,4,30))

# plot using the ggplot

ggplot(run.time, aes(x = nodes, y = time)) +
  geom_point(size = 3, color = "#69b3a2") +  # Soft teal green for points
  geom_smooth(method = "lm", color = "#FFB6C1", se = FALSE) +  # Soft pink for line
  labs(title = "",
       x = "",
       y = "") +
  theme(text = element_text(size = 30.0,family = "serif"),plot.title = element_text(size=30))+
  theme_minimal()

run.time <- data.frame(nodes = c(100,200,1000,2000,5000,20000,100000), time = c(2,4,30,60,135,700,3500))

# plot using the ggplot

ggplot(run.time, aes(x = nodes, y = time)) +
  geom_point(size = 3, color = "#69b3a2") +  # Soft teal green for points
  geom_smooth(method = "lm", color = "#FFB6C1", se = FALSE) +  # Soft pink for line
  labs(title = "",
       x = "Number of Nodes",
       y = "Time (seconds)") +
  theme(text = element_text(size = 17.0,family = "serif"),plot.title = element_text(size=12))+
  theme_minimal()





# Install required packages if not already installed
# install.packages("ggplot2")
# install.packages("patchwork")

library(ggplot2)
library(patchwork)

# First dataset (smaller)
run.time1 <- data.frame(nodes = c(100, 200, 1000), 
                        time = c(2, 4, 30))

# Second dataset (larger)
run.time2 <- data.frame(nodes = c(100, 200, 1000, 2000, 5000, 20000, 100000), 
                        time = c(2, 4, 30, 60, 135, 700, 3500))

# Create the first plot with larger font size
p1 <- ggplot(run.time1, aes(x = nodes, y = time)) + 
  geom_point(size = 3, color = "#69b3a2") + 
  geom_smooth(method = "lm", color = "#FFB6C1", se = FALSE) + 
  labs(title = "", x = "", y = "") + 
  theme_minimal() +
  theme(text = element_text(size = 15, family = "serif"),
        plot.title = element_text(size = 15))

# Create the second plot
p2 <- ggplot(run.time2, aes(x = nodes, y = time)) + 
  geom_point(size = 3, color = "#69b3a2") + 
  geom_smooth(method = "lm", color = "#FFB6C1", se = FALSE) + 
  labs(title = "", x = "Number of Nodes", y = "Time (seconds)") + 
  theme_minimal() +
  theme(text = element_text(size = 15.0, family = "serif"),
        plot.title = element_text(size = 15))

# A better approach using viewport
library(grid)

# Create the combined plot
combined_plot <- function() {
  # Set up the main plot (p2) to fill the entire device
  print(p2)
  
  # Create a viewport for p1 in the top-left corner
  vp <- viewport(x = 0.25, y = 0.75, width = 0.4, height = 0.4, just = c("center", "center"))
  
  # Push the viewport and print p1
  pushViewport(vp)
  print(p1, newpage = FALSE)
  popViewport()
}

# Display the combined plot
combined_plot()
ggsave(file="fig6_updated.pdf", width = 210, height = 297, units = "mm")

