# An implementation from Stadler et al 2012 & Stadler 2010
library(stats4)
# inputs: node time (yi) and tip time (xi), maybe psi
# outputs: ML estimates for speciation and extinction rate

# Here we calculate the probability density of the tree conditioned on 
# 1 extand sampled individuals.

# See eq.9 , Theorem 3.11 for detail.

# Note that in our model we do not consider sampled individual with 
# sampled descendants, therefore k = 0.

#c1
c1 <- function(lambda,mu,psi){
  return(abs(sqrt((lambda-mu-psi)^2+4*lambda*psi)))
}
#c2
c2<- function(lambda,mu,psi){
  return(-(lambda-mu-psi)/(c1(lambda,mu,psi)))
}

p0 <- function(lambda,mu,psi,t){
  d1 <- lambda + mu + psi
  d2 <- 2 * lambda
  d3 <- exp(-c1(lambda,mu,psi)*t)*(1-c2(lambda,mu,psi)) - (1+c2(lambda,mu,psi))
  d4 <- exp(-c1(lambda,mu,psi)*t)*(1-c2(lambda,mu,psi)) + (1+c2(lambda,mu,psi))
  d5 <- (c1(lambda,mu,psi)*d3)/d4
  return((d1+d5)/d2)
}

q <- function(lambda,mu,psi,t){
  d1 <- 2 * (1-(c2(lambda,mu,psi))^2)
  d2 <- exp(-c1(lambda,mu,psi)*t)*(1-c2(lambda,mu,psi))^2
  d3 <- exp(c1(lambda,mu,psi)*t)*(1+c2(lambda,mu,psi))^2
  return(d1+d2+d3)
}

#test
# t=1
# p0(lambda,mu,psi,rho,t)
# p1(lambda,mu,psi,rho,t)


# probability density 

probability_density <- function(lambda,mu,psi,m,X,Y){
  # Add checks for invalid inputs
  sortX <- sort(X)
  sortY <- sort(Y)
  d1 <- m*log(lambda)
  d2 <- sum(-log(q(lambda,mu,psi,sortX)))
  d3 <- sum(log(psi)+log(q(lambda,mu,psi,sortY)))
  # d1 <- lambda^m
  # d2 <- prod(1/q(lambda,mu,psi,sortX))
  # d3 <- prod(psi*q(lambda,mu,psi,sortY))
  # 
  # result <- d1*d2*d3
  result <- d1+d2+d3
  # Return small positive number if result is invalid
  #if(!is.finite(result) || result <= 0) return(1e-10)
  
  return(result)
}

lambda = 1 
mu = 0.5732560
psi = 0.14
probability_density(lambda=5,mu,psi,phy$Nnode,X,Y)



##### 
#test with simulated trees

pars <- c(2,  2,
          0.05, 0.05,
          .1,  .1,
          .05,
          .05)

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

# phy <- tree.bisse_test(pars, x0=1, include.extinct=FALSE, max.taxa = 1000, max.t = 1000)
# length(phy$tip.label)
# h <- history.from.sim.discrete(phy, 1:2)
# plot(h, phy)
phy <- tree.bisse_test(pars, x0=1, include.extinct=FALSE, 1000, 1000)
# Ensure the tree has at least 30 tips
phy
# while(length(phy$tip.label) < 5) {
#   phy <- tree.bisse_test(pars, x0=1, include.extinct=FALSE, 10, 10)
# }
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
true_phy
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
  ggtitle("True Phylogeny") + theme_tree2()+ geom_tiplab()+geom_nodelab()+
  theme(text = element_text(size = 12,family = "serif",face="bold"),plot.title = element_text(size=12))
true_phy_new

# extract branching time

time <- branching.times(phy)




extract_tree_times <- function(tree) {
  # Get total tree height (max path length from root to tip)
  tree_height <- max(node.depth.edgelength(tree))
  
  # Get all node depths (distance from root)
  all_nodes_depth <- node.depth.edgelength(tree)
  
  # Convert depths to times (time from present)
  all_times <- tree_height - all_nodes_depth
  
  # Separate internal nodes and tips
  n_tips <- length(tree$tip.label)
  
  # Extract tip times
  tip_times <- all_times[1:n_tips]
  names(tip_times) <- tree$tip.label
  
  # Extract internal node times
  internal_times <- all_times[(n_tips + 1):length(all_times)]
  names(internal_times) <- (n_tips + 1):length(all_times)
  
  # Create data frames
  internal_df <- data.frame(
    node = names(internal_times),
    time = as.numeric(internal_times),
    row.names = NULL
  )
  
  tip_df <- data.frame(
    tip = names(tip_times),
    time = as.numeric(tip_times),
    row.names = NULL
  )
  
  return(list(
    internal_nodes = internal_df,
    tips = tip_df
  ))
}


branching_time <- extract_tree_times(phy)
X <- as.vector(branching_time$internal_nodes)
X <- X$time
Y <- as.vector(branching_time$tips)
Y <- Y$time

negative_log_likelihood <- function(params, m, X, Y) {
  lambda <- params[1]
  mu <- params[2]
  psi <- params[3]
  
  # Ensure parameters are positive
  if (lambda <= 0 || mu <= 0 || psi <= 0) return(Inf)
  
  # Compute the probability density
  likelihood <- probability_density(lambda, mu, psi, m, X, Y)
  likelihood
  # if (!is.finite(likelihood) || likelihood <= 0) {
  #   likelihood <- 1e-10  # Replace with a small positive number
  # }
  # Return negative log-likelihood
  return(-likelihood)
}


# Initial parameter guesses
# init_params <- c(lambda = 1, mu = 0.05, psi = 0.1)

# # Perform MLE using optim
# mle_result <- mle(
#   par = init_params, 
#   fn = negative_log_likelihood, 
#   m = m, 
#   X = X, 
#   Y = Y, 
#   method = "L-BFGS-B", 
#   lower = c(1e-5, 1e-5, 1e-5),  # Ensure positivity
#   upper = c(2,1,1)
# )
# mle_result


#Find MLE with multiple starting points


# Define the negative log likelihood function
nLL <- function(lambda, mu,psi) {
  # Return negative log likelihood
  # We use negative because mle() minimizes
  -(probability_density(lambda, mu,psi, m=phy$Nnode, X=X, Y=Y))
}

# Find MLE
fit <- mle(nLL, 
           start = list(lambda = 1, mu = 1,psi=0.2),
           method = "L-BFGS-B",
           lower = c(0.0001, 0.0001,0.0001),
           upper = c(5, 5,5))

# View results
coef(fit)  # Get parameter estimates

# nLL <- function(params) {
#   lambda <- params[1]
#   mu <- params[2]
#   psi <- params[3]
#   if(lambda <= 0 || mu <= 0) return(1e10)
#   ll <- try({
#     -(probability_density(lambda, mu, psi, m=phy$Nnode, X=X, Y=Y))
#   }, silent=TRUE)
#   if(inherits(ll, "try-error") || !is.finite(ll)) return(1e10)
#   return(ll)
# }
# 
# fits <- list()
# starts <- list(
#   c(0.9, 0.05,0.1),
#   c(1.5, 0.1,0.2),
#   c(2.0, 0.2,0.4)
# )


# for(i in seq_along(starts)) {
#   fits[[i]] <- try({
#     mle(nLL,
#         start = list(lambda = starts[[i]][1], mu = starts[[i]][2],psi=starts[[i]][3]),
#         method = "L-BFGS-B",
#         lower = c(0.5, 0.01,0.1),
#         upper = c(10, 1,0.2))
#   }, silent=TRUE)
# }
# 
# # Use the best fit that converged
# valid_fits <- fits[!sapply(fits, inherits, "try-error")]
# if(length(valid_fits) > 0) {
#   best_fit <- valid_fits[[which.min(sapply(valid_fits, function(x) AIC(x)))]]
# }
# best_fit
