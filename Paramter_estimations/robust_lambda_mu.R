replace_matrix_with_vector <- function(matrix, vector) {
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      matrix[i,j] <- vector[matrix[i,j]]
    }
  }
  return(matrix)
}

tree_simulation_w_error <- function(pars,q_matrix){
  phy <- sim.bds.tree(pars, x0=1, include.extinct=FALSE, 1000, 1000)
  while(length(phy$tip.label) < 100) {
    phy <- sim.bds.tree(pars, x0=1, include.extinct=FALSE, 1000, 1000)
  }
  phy <- prune(phy)
  print('make the tree')
  h2 <- history.from.sim.discrete(phy, 1:2)
  
  node_depths <- node.depth.edgelength(phy)
  tmrca <- max(node_depths)
  tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.001]
  new_phy <- drop.tip(phy, tips_to_drop)
  
  plot(new_phy)
  phy_our <- new_phy
  
  true_phy_info_new <- as_tibble(new_phy)
  phy_data <- c(factor(h2$tip.state),factor(h2$node.state))
  true_phy_info_new <- true_phy_info_new %>% mutate(State = phy_data[label])
  new_phy$tip.state <- new_phy$tip.state[setdiff(names(new_phy$tip.state), tips_to_drop)]
  
  phy <- new_phy
  
  branching_time <- extract_tree_times(phy)
  X <- as.vector(branching_time$internal_nodes)
  X <- X$time
  #print(X)
  Y <- as.vector(branching_time$tips)
  Y <- Y$time
  
  nLL <- function(lambda, mu) {
    # Return negative log likelihood
    # We use negative because mle() minimizes
    -(probability_density(lambda, mu,psi=0.5, m=phy$Nnode, X=X, Y=Y))
  }
  
  fit <- mle(nLL,
             start = list(lambda = 2, mu = 0.5),
             method = "L-BFGS-B",
             lower = c(0.1, 0.05),
             upper = c(10, 10))
  print(coef(fit))
  # View results
  lm_rates <- unname(c(coef(fit)[1],coef(fit)[2]))
  return(list(
    lm_rates=lm_rates
  ))
}

pars <- c(5,5,      
          1, 1,
          .5, .5, 
          .1,
          .1)       
# generate 10 trees
results <- list()

qij_matrix <- function(k,val) {
  mat <- matrix(val, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2,0.05)
for (i in 1:500){
  print(i)
  set.seed(i)
  results[[i]] <- tree_simulation_w_error(pars,q_matrix)
}

####################
# plots
####################
lm_rate <- c()
for(i in 1:length(results)){
  lm_rate=c(lm_rate,results[[i]]$lm_rates)
}


separate_entries <- function(vector) {
  # Check if length is even
  if (length(vector) %% 2 != 0) {
    stop("Input vector must have an even number of elements")
  }
  
  # Get odd entries (1st, 3rd, 5th, etc.)
  odd_entries <- vector[seq(1, length(vector), by = 2)]
  
  # Get even entries (2nd, 4th, 6th, etc.)
  even_entries <- vector[seq(2, length(vector), by = 2)]
  
  # Return as a list
  return(list(x1 = odd_entries, x2 = even_entries))
}



data1 <- separate_entries(lm_rate)
data1raw <- data1
data1$x1 <- c(5,data1$x1)
data1$x2 <- c(1,data1$x2)
data1$special <- data1$x1 == 5 & data1$x2 == 1

data1 <- data.frame(x1 = data1$x1,x2=data1$x2,special=data1$special)
mean(data1$x1)
mean(data2$x2)


plot1 <- ggplot() +
  # Add scatter points
  geom_point(data = data1, aes(y = x1, x = seq_along(x1)), color = "black") +
  # Add horizontal line at 0.05
  geom_hline(yintercept = 5, color = "red") +
  geom_hline(yintercept = mean(data1$x1), color = "brown") +
  # Customize the theme
  theme_minimal() +
  # Labels
  labs(y = "Values", x = "Index")

plot2 <- ggplot() +
  # Add scatter points
  geom_point(data = data1, aes(y = x2, x = seq_along(x2)), color = "black") +
  # Add horizontal line at 0.1
  geom_hline(yintercept = 1, color = "red") +
  geom_hline(yintercept = mean(data1$x2), color = "brown") +
  # Customize the theme
  theme_minimal() +
  # Labels
  labs(y = "Values", x = "Index")



# First plot for x1
plot1 <- ggplot() +
  # Add scatter points
  geom_point(data = data1, aes(y = x1, x = seq_along(x1)), color = "black") +
  # Add horizontal lines with labels for legend
  geom_hline(aes(yintercept = 5, color = "True rate")) +
  geom_hline(aes(yintercept = mean(data1$x1), color = "Mean")) +
  # Customize colors and legend
  scale_color_manual(values = c("True rate" = "red", "Mean" = "blue")) +
  # Customize the theme
  theme_minimal() +
  theme(text = element_text(size = 12,family = "serif",face="bold"),plot.title = element_text(size=12))+
  # Labels
  labs(y = "lambda Estimates", 
       x="",
       color = "") # This changes the legend title

# Second plot for x2
plot2 <- ggplot() +
  # Add scatter points
  geom_point(data = data1, aes(y = x2, x = seq_along(x2)), color = "black") +
  # Add horizontal lines with labels for legend
  geom_hline(aes(yintercept = 1, color = "True rate")) +
  geom_hline(aes(yintercept = mean(data1$x2), color = "Mean")) +
  # Customize colors and legend
  scale_color_manual(values = c("True rate" = "red", "Mean" = "blue")) +
  # Customize the theme
  theme_minimal() +
  theme(text = element_text(size = 12,family = "serif",face="bold"),plot.title = element_text(size=12))+
  # Labels
  labs(y = "mu Estimates", 
       x="",
       color = "")

# Display plots
print(plot1)
print(plot2)

lm_estimates <- ggarrange(plot1,plot2,nrow = 1, common.legend = TRUE)





