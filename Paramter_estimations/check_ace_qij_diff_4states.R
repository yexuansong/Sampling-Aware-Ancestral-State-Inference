replace_matrix_with_vector <- function(matrix, vector) {
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      matrix[i,j] <- vector[matrix[i,j]]
    }
  }
  return(matrix)
}


tree_simulation_w_error <- function(pars,q_matrix){
  # first generate the tree
  phy <- sim.bds.tree(pars, x0=1, include.extinct=FALSE, 15000, 15000)
  # Ensure the tree has at least 200 tips
  while(length(phy$tip.label) < 200) {
    phy <- sim.bds.tree(pars, x0=1, include.extinct=FALSE, 15000, 15000)
  }
  phy <- prune(phy)
  print('make the tree')
  h2 <- history.from.sim.discrete(phy, 1:4)
  
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
  
  
  tip_states <- h2$tip.state
  tip_states <- tip_states[!names(tip_states) %in% tips_to_drop]
  #####################
  # perform ace
  #####################
  phy <- as.phylo(phy)
  phy$node.label <- NULL
  # Perform ASR
  asr_ace <- ace(tip_states, phy, type="discrete", model="ARD")
  asr_rates <- asr_ace$rates
  print('finish ace')
  return(list(
    asr_rates=asr_rates
  ))
}

pars <- c(2,  2, 2, 2,
          .045, .045,.045,.045,
          .5,  .5, .5,.5,
          .1,.1,.1,.1,
          .1,.1,.1,.1,
          .1,.1,.1,.1)


results <- list()

qij_matrix <- function(k,val) {
  mat <- matrix(val, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2,0.05)
# generate 500 trees
for (i in 1:500){
  print(i)
  set.seed(i)
  results[[i]] <- tree_simulation_w_error(pars,q_matrix)
}

####################
# accuracy plots
####################
ace_rate <- c()
for(i in 1:length(results)){
  ace_rate=c(ace_rate,results[[i]]$asr_rates)
}


separate_entries <- function(vector) {

  # Get odd entries (1st, 3rd, 5th, etc.)
  entry1 <- vector[seq(1, length(vector), by = 12)]
  
  # Get even entries (2nd, 4th, 6th, etc.)
  entry2 <- vector[seq(2, length(vector), by = 12)]
  
  # Get odd entries (1st, 3rd, 5th, etc.)
  entry3 <- vector[seq(3, length(vector), by = 12)]
  
  # Get even entries (2nd, 4th, 6th, etc.)
  entry4 <- vector[seq(4, length(vector), by = 12)]
  
  # Get odd entries (1st, 3rd, 5th, etc.)
  entry5 <- vector[seq(5, length(vector), by = 12)]
  
  # Get even entries (2nd, 4th, 6th, etc.)
  entry6 <- vector[seq(6, length(vector), by = 12)]
  
  # Get odd entries (1st, 3rd, 5th, etc.)
  entry7 <- vector[seq(7, length(vector), by = 12)]
  
  # Get even entries (2nd, 4th, 6th, etc.)
  entry8 <- vector[seq(8, length(vector), by = 12)]
  
  # Get odd entries (1st, 3rd, 5th, etc.)
  entry9 <- vector[seq(9, length(vector), by = 12)]
  
  # Get even entries (2nd, 4th, 6th, etc.)
  entry10 <- vector[seq(10, length(vector), by = 12)]
  
  # Get odd entries (1st, 3rd, 5th, etc.)
  entry11 <- vector[seq(11, length(vector), by = 12)]
  
  # Get even entries (2nd, 4th, 6th, etc.)
  entry12 <- vector[seq(12, length(vector), by = 12)]
  
  # Return as a list
  return(list(x1 = entry1, x2 = entry2, x3 = entry3, x4=entry4,
              x5=entry5,x6=entry6,x7=entry7,x8=entry8,x9=entry9,
              x10=entry10,x11=entry11,x12=entry12))
}



data1 <- separate_entries(ace_rate)
data1raw <- data1

###########

pars <- c(2,  2, 2, 2,
          .045, .045,.045,.045,
          .25,  .5, .5,.5,
          .1,.1,.1,.1,
          .1,.1,.1,.1,
          .1,.1,.1,.1)


# generate 10 trees
results <- list()

qij_matrix <- function(k,val) {
  mat <- matrix(val, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2,0.05)
#q_matrix[2,1] = 0.2
for (i in 1:500){
  print(i)
  set.seed(i)
  results[[i]] <- tree_simulation_w_error(pars,q_matrix)
}

ace_rate <- c()
for(i in 1:length(results)){
  ace_rate=c(ace_rate,results[[i]]$asr_rates)
}


data2 <- separate_entries(ace_rate)
data2raw <- data2


###########
  
pars <- c(2,  2, 2, 2,
          .045, .045,.045,.045,
          .1,  .5, .5,.5,
          .1,.1,.1,.1,
          .1,.1,.1,.1,
          .1,.1,.1,.1)

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

ace_rate <- c()
for(i in 1:length(results)){
  ace_rate=c(ace_rate,results[[i]]$asr_rates)
}
data3 <- separate_entries(ace_rate)
data3raw <- data3

###########
pars <- c(2,  2, 2, 2,
          .045, .045,.045,.045,
          .05,  .5, .5,.5,
          .1,.1,.1,.1,
          .1,.1,.1,.1,
          .1,.1,.1,.1)

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

ace_rate <- c()
for(i in 1:length(results)){
  ace_rate=c(ace_rate,results[[i]]$asr_rates)
}
data4 <- separate_entries(ace_rate)
data4raw <- data4


###########
pars <- c(2,  2, 2, 2,
          .045, .045,.045,.045,
          .01,  .5, .5,.5,
          .1,.1,.1,.1,
          .1,.1,.1,.1,
          .1,.1,.1,.1)

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

ace_rate <- c()
for(i in 1:length(results)){
  ace_rate=c(ace_rate,results[[i]]$asr_rates)
}


data5 <- separate_entries(ace_rate)
data5raw <- data5

###########
pars <- c(2,  2, 2, 2,
          .045, .045,.045,.045,
          .005,  .5, .5,.5,
          .1,.1,.1,.1,
          .1,.1,.1,.1,
          .1,.1,.1,.1)

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

ace_rate <- c()
for(i in 1:length(results)){
  ace_rate=c(ace_rate,results[[i]]$asr_rates)
}


data7 <- separate_entries(ace_rate)
data7raw <- data7


replace_matrix_with_vector <- function(matrix, vector) {
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      matrix[i,j] <- vector[matrix[i,j]]
    }
  }
  return(matrix)
}


qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,qij_rate)
results[[1]]$asr_rates
qij_matrix <- replace_matrix_with_vector(ans_ace$index.matrix,results[[1]]$asr_rates)

save(data1raw, data2raw,data3raw,data4raw,data5raw,data6raw,data7raw, file="qij_rawdata_updated_4states.RData")
