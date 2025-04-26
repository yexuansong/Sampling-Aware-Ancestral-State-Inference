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
  #set.seed(NULL)  # Reset seed for each simulation
  phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 1000, 1000)
  # Ensure the tree has at least 30 tips
  while(length(phy$tip.label) < 50) {
    phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 1000, 1000)
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
  # Extract tip states
  #tip_states <- phy$tip.state
  # Perform ASR
  asr_ace <- ace(tip_states, phy, type="discrete", model="SYM")
  asr_rates <- asr_ace$rates
  
  qij_matrix <- replace_matrix_with_vector(asr_ace$index.matrix,asr_ace$rates)
  
  print('finish ace')
  #####################
  # perform simmap
  #####################
  # fit the model
  tip_states <- h2$tip.state
  tip_states <- tip_states[!names(tip_states) %in% tips_to_drop]
  smap_model <- fitMk(phy,tip_states,model="ER",pi="fitzjohn")
  phy_smap<-simmap(smap_model)
  
  sm<-summary(phy_smap)
  simmap_lik <- as.data.frame(sm$ace)
  print('finish simmap')
  
  #####################
  # perform our
  #####################
  asr_our <- saasi(pars,phy,q_matrix)
  
  print('finish our correct')
  
  
  modified_phy <- ace_phy
  modified_phy$orig <- NULL
  modified_phy$hist <- NULL
  modified_phy$tip.state <- NULL
  modified_phy$node.state <- NULL
  modified_phy$edge.state <- NULL
  netg <- ltt(modified_phy)
  lmn <- lm(log(netg$ltt)~netg$times)
  lmn$coefficients[2]
  #assume the relation of l and m is 10, 10m = l, l-m =0.12
  mu <- lmn$coefficients[2]/9
  lambda <- 10*lmn$coefficients[2]/9
  mu<-unname(mu)
  lambda<- unname(lambda)
  pars[1:2] <-lambda
  pars[3:4] <- mu
  
  asr_our_est <- saasi(pars,phy,qij_matrix)
  
  print('finish our estimate')
  
  #####################
  # perform accuracy check
  #####################  
  #true_tip <- c(factor(h2$tip.state),factor(h2$node.state))
  
  
  acc_asr <- accuracy_helper(asr_ace$lik.anc,true_phy_info_new$State,tip_states)
  acc_sim <- accuracy_helper(head(simmap_lik,n=length(h2$node.state)- length(tips_to_drop)),true_phy_info_new$State,tip_states)
  acc_our <- accuracy_helper(asr_our[,1:2],true_phy_info_new$State,tip_states)
  acc_our_est <- accuracy_helper(asr_our_est[,1:2],true_phy_info_new$State,tip_states)
  
  #return the accuracies
  
  return(list(
    asr=acc_asr,
    sim=acc_sim,
    our=acc_our,
    our_est=acc_our_est,
    asr_rates=asr_rates
  ))
}

pars <- c(1,  1,      
          .045, .045,
          .01,  .5, 
          .05,
          .05)       



# generate 10 trees
results <- list()

qij_matrix <- function(k,val) {
  mat <- matrix(val, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2,0.05)
#q_matrix[2,1] = 0.2
for (i in 1:100){
  print(i)
  set.seed(i)
  results[[i]] <- tree_simulation_w_error(pars,q_matrix)
}


####################
# accuracy plots
####################

# ace accuracy
out_acc_ace <- c()
out_ka_ace <- c()
out_prob_ace <- c()
ace_rate <- c()
for(i in 1:length(results)){
  out_acc_ace=c(out_acc_ace,results[[i]]$asr$acc)
  out_ka_ace=c(out_ka_ace, results[[i]]$asr$ka)
  out_prob_ace=c(out_prob_ace, results[[i]]$asr$prob_acc)
  ace_rate=c(ace_rate,results[[i]]$asr_rates)
}

# simap accuracy
out_acc_sim <- c()
out_ka_sim <- c()
out_prob_sim <- c()
for(i in 1:length(results)){
  out_acc_sim=c(out_acc_sim,results[[i]]$sim$acc)
  out_ka_sim=c(out_ka_sim, results[[i]]$sim$ka)
  out_prob_sim=c(out_prob_sim, results[[i]]$sim$prob_acc)
}


# our accuracy
out_acc_our <- c()
out_ka_our <- c()
out_prob_our <- c()
for(i in 1:length(results)){
  out_acc_our=c(out_acc_our,results[[i]]$our$acc)
  out_ka_our=c(out_ka_our, results[[i]]$our$ka)
  out_prob_our=c(out_prob_our, results[[i]]$our$prob_acc)
}


# our accuracy est
out_acc_our_est <- c()
out_ka_our_est <- c()
out_prob_our_est <- c()
for(i in 1:length(results)){
  out_acc_our_est=c(out_acc_our_est,results[[i]]$our_est$acc)
  out_ka_our_est=c(out_ka_our_est, results[[i]]$our_est$ka)
  out_prob_our_est=c(out_prob_our_est, results[[i]]$our_est$prob_acc)
}

#overall acc comparison 
plot_data_acc <- data.frame(
  Ace = out_acc_ace,
  Simmap = out_acc_sim,
  SAASI = out_acc_our,
  SAASI_est = out_acc_our_est
)

acc_data <- pivot_longer(plot_data_acc, 
                         cols = c(Ace,Simmap,SAASI,SAASI_est),
                         names_to = "Metric",
                         values_to = "Accuracy")
acc_data$Metric <- factor(acc_data$Metric, levels = c("Ace", "Simmap", "SAASI","SAASI_est"))


p1 <- ggplot(acc_data, aes(x = Metric, y = Accuracy)) +
  geom_violin(width = 1.1, fill = NA, color = "black") +  # Remove fill, keep border
  geom_jitter(height = 0, width = 0.1) +
  geom_boxplot(width = 0.2, fill = "white", color = "black") +  # Keep boxplot visible
  labs(x = "", y = "Accuracy", fill = "Methods") +
  ggtitle("A: Absolute accuracy") +
  scale_x_discrete(labels = c("Ace" = "ace", 
                              "Simmap" = "simmap", 
                              "SAASI" = "SAASI", 
                              "SAASI_est" = "SAASI - estimated paramters")) +
  theme_minimal() +
  ylim(0.25, 1.0) +
  theme(text = element_text(size = 15, family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 15))  # Rotate x-axis labels

p1


#probabilistic acc comparison
plot_data_prob <- data.frame(
  Ace = out_prob_ace,
  Simmap = out_prob_sim,
  SAASI = out_prob_our,
  SAASI_est = out_prob_our_est
)

acc_prob <- pivot_longer(plot_data_prob, 
                         cols = c(Ace,Simmap,SAASI,SAASI_est),
                         names_to = "Metric",
                         values_to = "Accuracy")
acc_prob$Metric <- factor(acc_prob$Metric, levels = c("Ace", "Simmap", "SAASI","SAASI_est"))

# Create the plot
p3 <- ggplot(acc_prob, aes(x = Metric, y = Accuracy)) +
  geom_violin(width = 1.1, fill = NA, color = "black") +  # Remove fill, keep border
  geom_jitter(height = 0, width = 0.1) +
  geom_boxplot(width = 0.2, fill = "white", color = "black") +  # Keep boxplot white
  ggtitle("B: Probability accuracy") +
  labs(x = "", y = "", fill = "Methods") +
  scale_x_discrete(labels = c("Ace" = "ace", 
                              "Simmap" = "simmap", 
                              "SAASI" = "SAASI", 
                              "SAASI_est" = "SAASI - estimated parameters")) +
  theme_minimal() + 
  ylim(0.25, 1.0) +
  theme(text = element_text(size = 15, family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 15))  # Rotate labels correctly
p3


accuracy_comparison <- ggarrange(p1,p3,widths = c(5,5),nrow = 1)
ggarrange(p1,p3,nrow = 1)

ggsave(file="fig2_updated.pdf", width = 210, height = 297, units = "mm")


