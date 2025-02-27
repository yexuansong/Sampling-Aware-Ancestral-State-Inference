replace_matrix_with_vector <- function(matrix, vector) {
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      matrix[i,j] <- vector[matrix[i,j]]
    }
  }
  return(matrix)
}


tree_simulation_w_error <- function(pars,q_matrix){
  phy <- sim.bds.tree(pars, x0=1, include.extinct=FALSE, 15000, 15000)
  # Ensure the tree has at least 1000 tips
  while(length(phy$tip.label) < 1000) {
    phy <- sim.bds.tree(pars, x0=1, include.extinct=FALSE, 15000, 15000)
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
  
  
  tip_states <- h2$tip.state
  tip_states <- tip_states[!names(tip_states) %in% tips_to_drop]
  #####################
  # perform ace
  #####################
  phy <- as.phylo(phy)
  phy$node.label <- NULL
  asr_ace <- ace(tip_states, phy, type="discrete", model="ARD")
  asr_rates <- asr_ace$rates
  print('finish ace')
  return(list(
    asr_rates=asr_rates
  ))
}

pars <- c(1,  1,      
          .045, .045,
          .5,  .5, 
          .1,
          .05)       
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



data1 <- separate_entries(ace_rate)
data1raw <- data1
data1$x1 <- c(0.05,data1$x1)
data1$x2 <- c(0.1,data1$x2)
data1$special <- data1$x1 == 0.05 & data1$x2 == 0.1

data1 <- data.frame(x1 = data1$x1,x2=data1$x2,special=data1$special)

robust_qij_equal<-ggplot(data1, aes(x = x1, y = x2)) +
  geom_point(data1 = subset(data1, !special), color = "black", size = 2) +
  geom_point(data1 = subset(data1, special), color = "red", shape = 8, size = 4) +
  ggtitle("A") + 
  theme(text = element_text(size = 12,family = "serif",face="bold"),plot.title = element_text(size=12))+
  theme_light()+
  labs(x = "q12", y = "q21") + ylim(0,1)+xlim(0,0.2)

robust_qij_equal

###########

pars <- c(1,  1,      
          .045, .045,
          .25,  .5, 
          .1,
          .05)       

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


data2 <- separate_entries(ace_rate)
data2raw <- data2

data2$x1 <- c(0.05,data2$x1)
data2$x2 <- c(0.1,data2$x2)
data2$special <- data2$x1 == 0.05 & data2$x2 == 0.1

data2 <- data.frame(x1 = data2$x1,x2=data2$x2,special=data2$special)

robust_qij_2x<-ggplot(data2, aes(x = x1, y = x2)) +
  geom_point(data2 = subset(data2, !special), color = "black", size = 2) +
  geom_point(data2 = subset(data2, special), color = "red", shape = 8, size = 4) +
  ggtitle("B") + 
  theme(text = element_text(size = 12,family = "serif",face="bold"),plot.title = element_text(size=12))+
  theme_light()+
  labs(x = "q12", y = "q21") + ylim(0,1)+xlim(0,0.2)

robust_qij_2x

###########

pars <- c(1,  1,      
          .045, .045,
          .1,  .5, 
          .1,
          .05)       


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

data3$x1 <- c(0.05,data3$x1)
data3$x2 <- c(0.1,data3$x2)
data3$special <- data3$x1 == 0.05 & data3$x2 == 0.1

data3 <- data.frame(x1 = data3$x1,x2=data3$x2,special=data3$special)

robust_qij_5x<-ggplot(data3, aes(x = x1, y = x2)) +
  geom_point(data3 = subset(data3, !special), color = "black", size = 2) +
  geom_point(data3 = subset(data3, special), color = "red", shape = 8, size = 4) +
  ggtitle("C") + 
  theme(text = element_text(size = 12,family = "serif",face="bold"),plot.title = element_text(size=12))+
  theme_light()+
  labs(x = "q12", y = "q21") + ylim(0,1)+xlim(0,0.2)

robust_qij_5x




###########

pars <- c(1,  1,      
          .045, .045,
          .05,  .5, 
          .1,
          .05)       
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

data4$x1 <- c(0.05,data4$x1)
data4$x2 <- c(0.1,data4$x2)
data4$special <- data4$x1 == 0.05 & data4$x2 == 0.1

data4 <- data.frame(x1 = data4$x1,x2=data4$x2,special=data4$special)



robust_qij_10x<-ggplot(data4, aes(x = x1, y = x2)) +
  geom_point(data4 = subset(data4, !special), color = "black", size = 2) +
  geom_point(data4 = subset(data4, special), color = "red", shape = 8, size = 4) +
  ggtitle("D") + 
  theme(text = element_text(size = 12,family = "serif"),plot.title = element_text(size=12))+
  theme_light()+
  labs(x = "q12", y = "q21") + ylim(0,1)+xlim(0,0.2)

robust_qij_10x
###########

pars <- c(1,  1,      
          .045, .045,
          .01,  .5, 
          .1,
          .05)       


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

data5$x1 <- c(0.05,data5$x1)
data5$x2 <- c(0.1,data5$x2)
data5$special <- data5$x1 == 0.05 & data5$x2 == 0.1

data5 <- data.frame(x1 = data5$x1,x2=data5$x2,special=data5$special)

robust_qij_50x<-ggplot(data5, aes(x = x1, y = x2)) +
  geom_point(data5 = subset(data5, !special), color = "black", size = 2) +
  geom_point(data5 = subset(data5, special), color = "red", shape = 8, size = 4) +
  ggtitle("E") + 
  theme(text = element_text(size = 12,family = "serif"),plot.title = element_text(size=12))+
  theme_light()+
  labs(x = "q12", y = "q21") + ylim(0,1)+xlim(0,0.2)

robust_qij_50x


###########

pars <- c(1,  1,      
          .045, .045,
          .005,  .5, 
          .1,
          .05)       
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


data6 <- separate_entries(ace_rate)
data6raw <- data6

data6$x1 <- c(0.05,data6$x1)
data6$x2 <- c(0.1,data6$x2)
data6$special <- data6$x1 == 0.05 & data6$x2 == 0.1

data6 <- data.frame(x1 = data6$x1,x2=data6$x2,special=data6$special)



robust_qij_100x<-ggplot(data6, aes(x = x1, y = x2)) +
  geom_point(data6 = subset(data6, !special), color = "black", size = 2) +
  geom_point(data6 = subset(data6, special), color = "red", shape = 8, size = 4) +
  ggtitle("F") + 
  theme(text = element_text(size = 12,family = "serif"),plot.title = element_text(size=12))+
  theme_light()+
  labs(x = "q12", y = "q21") + ylim(0,1)+xlim(0,0.2)

robust_qij_100x

new_plot_robust_comb <- ggarrange(robust_qij_equal,robust_qij_2x,robust_qij_5x,robust_qij_10x,robust_qij_50x,robust_qij_100x,
                                  nrow=2, ncol = 3)



pars <- c(3,  3,      
            .045, .045,
            .01,  1, 
            .1,
            .05)       
  
  
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
  
  data7$x1 <- c(0.05,data7$x1)
  data7$x2 <- c(0.1,data7$x2)
  data7$special <- data7$x1 == 0.05 & data7$x2 == 0.1
  
  data7 <- data.frame(x1 = data7$x1,x2=data7$x2,special=data7$special)

  
save(data1raw, data2raw,data3raw,data4raw,data5raw,data6raw,data7raw, file="qij_rawdata_updated.RData")



