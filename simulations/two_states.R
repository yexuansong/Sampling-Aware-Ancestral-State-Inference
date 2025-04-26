source(file.path("simulation.R"))
source(file.path("test_bisse.R"))
source(file.path("accuracy_helper.R"))

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
library("janitor")
library("stringr")
library("phytools")
library(dplyr)
library(rstatix)
library(jsonlite)
library(rootSolve)
##################################### Parameter specification
pars <- c(1,  1,      
          .045, .045,
          .01,  .5, 
          .05,
          .05)   

k = 2

qij_matrix <- function(k) {
  mat <- matrix(0.05, nrow = k, ncol = k)  
  diag(mat) <- NA  
  return(mat)
}
q_matrix = qij_matrix(2)

#############################
# simulation
set.seed(4) 

phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 1000, 1000)
# Ensure the tree has at least 30 tips
while(length(phy$tip.label) < 20) {
  phy <- tree.sim(pars, x0=1, include.extinct=FALSE, 1000, 1000)
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
true_phy_new <- true_phy_new  %<+% true_phy_info_new + geom_point(aes(color=State),size=2) +
  ggtitle("A: True Phylogeny") +
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
true_phy_new


##################################### ASR using ace
# delete the column phy$node.label
#ace_phy <- new_phy
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
p1 <- p1 %<+% true_phy_info_new + geom_tippoint(aes(color=State),size=2)+
  ggtitle("C: ace") +
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
p1 <- inset(p1, ace_pie,width = 0.07,height = 0.07,hjust=0.005)
p1


asr_our <- saasi(pars,phy,q_matrix)


##################################### ASR Stochastic Character Mapping in Phytool
# fit the model
tip_states <- h$tip.state
tip_states <- tip_states[!names(tip_states) %in% tips_to_drop]
smap_model <- fitMk(phy,tip_states,model="SYM",pi="fitzjohn")

# do stochastic character mapping
phy_smap<-simmap(smap_model)

# print a summary of the stochastic mapping
sm<-summary(phy_smap)

simmap_lik <- as.data.frame(sm$ace)
simmap_lik <- head(simmap_lik,ace_phy$Nnode)

simmap_lik$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
simmap_pie <- nodepie(simmap_lik,cols=1:k)

p2 <- ggtree(ace_phy)
p2 <- p2 %<+% true_phy_info_new + geom_tippoint(aes(color=State))+
  ggtitle("ASR-Simmap") +
  theme(text = element_text(size = 15,family = "serif",face = "bold"),plot.title = element_text(size=15))
p2 <- inset(p2, ace_pie,width = 0.07,height = 0.07,hjust=0.005)
p2
######

asr_our <- saasi(pars,phy,q_matrix)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:k)

p3 <- ggtree(ace_phy)
p3 <- p3 %<+% true_phy_info_new + geom_tippoint(aes(color=State),size=2)+
  ggtitle("B: SAASI - true sampling") +
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
p3 <- inset(p3, our_pie,width = 0.07,height = 0.07,hjust=0.005)
p3


pars <- c(1,  1,      
          .045, .045,
          .5,  .5, 
          .05,
          .05)   

asr_our <- saasi(pars,phy,q_matrix)

asr_our$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
our_pie <- nodepie(asr_our,cols=1:k)

p4 <- ggtree(ace_phy)
p4 <- p4 %<+% true_phy_info_new + geom_tippoint(aes(color=State),size=2)+
  ggtitle("B: SAASI - equal sampling") +
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
p4 <- inset(p4, our_pie,width = 0.07,height = 0.07,hjust=0.005)
p4




ggarrange(true_phy_new,p4,p1, nrow=1,common.legend = TRUE,legend = "bottom")
ggsave(file="fig1_updated.pdf", width = 210, height = 297, units = "mm")

ggarrange(true_phy_new,p3,p1, nrow=1,common.legend = TRUE,legend = "bottom")
ggsave(file="figs-whatever_updated.pdf", width = 210, height = 297, units = "mm")




acc_our <- accuracy_helper(asr_our[,1:2],true_phy_info_new$State,tip_states)


pars <- c(1,  1,
          .045, .045,
          .05,  .5,
          .05,
          .05)

asr_our1 <- saasi(pars,phy,q_matrix)
acc_our1 <- accuracy_helper(asr_our1[,1:2],true_phy_info_new$State,tip_states)


pars <- c(1,  1,
          .045, .045,
          .1,  .5,
          .05,
          .05)

asr_our18 <- saasi(pars,phy,q_matrix)
acc_our18 <- accuracy_helper(asr_our18[,1:2],true_phy_info_new$State,tip_states)


pars <- c(1,  1,
          .045, .045,
          .15,  .5,
          .05,
          .05)

asr_our2 <- saasi(pars,phy,q_matrix)
acc_our2 <- accuracy_helper(asr_our2[,1:2],true_phy_info_new$State,tip_states)


pars <- c(1,  1,
          .045, .045,
          .2,  .5,
          .05,
          .05)

asr_our3 <- saasi(pars,phy,q_matrix)
acc_our3 <- accuracy_helper(asr_our3[,1:2],true_phy_info_new$State,tip_states)


pars <- c(1,  1,
          .045, .045,
          .25,  .5,
          .05,
          .05)

asr_our4 <- saasi(pars,phy,q_matrix)
acc_our4 <- accuracy_helper(asr_our4[,1:2],true_phy_info_new$State,tip_states)


pars <- c(1,  1,
          .045, .045,
          .3,  .5,
          .05,
          .05)

asr_our5 <- saasi(pars,phy,q_matrix)
acc_our5 <- accuracy_helper(asr_our5[,1:2],true_phy_info_new$State,tip_states)



pars <- c(1,  1,
          .045, .045,
          .35,  .5,
          .05,
          .05)

asr_our6 <- saasi(pars,phy,q_matrix)
acc_our6 <- accuracy_helper(asr_our6[,1:2],true_phy_info_new$State,tip_states)

pars <- c(1,  1,
          .045, .045,
          .4,  .5,
          .05,
          .05)

asr_our7 <- saasi(pars,phy,q_matrix)
acc_our7 <- accuracy_helper(asr_our7[,1:2],true_phy_info_new$State,tip_states)

pars <- c(1,  1,
          .045, .045,
          .45,  .5,
          .05,
          .05)

asr_our8 <- saasi(pars,phy,q_matrix)
acc_our8 <- accuracy_helper(asr_our8[,1:2],true_phy_info_new$State,tip_states)

pars <- c(1,  1,
          .045, .045,
          .5,  .5,
          .05,
          .05)

asr_our9 <- saasi(pars,phy,q_matrix)
acc_our9 <- accuracy_helper(asr_our9[,1:2],true_phy_info_new$State,tip_states)


pars <- c(1,  1,
          .045, .045,
          .55,  .5,
          .05,
          .05)

asr_our10 <- saasi(pars,phy,q_matrix)
acc_our10 <- accuracy_helper(asr_our10[,1:2],true_phy_info_new$State,tip_states)


pars <- c(1,  1,
          .045, .045,
          .6,  .5,
          .05,
          .05)

asr_our11 <- saasi(pars,phy,q_matrix)
acc_our11 <- accuracy_helper(asr_our11[,1:2],true_phy_info_new$State,tip_states)

pars <- c(1,  1,
          .045, .045,
          .65,  .5,
          .05,
          .05)

asr_our12 <- saasi(pars,phy,q_matrix)
acc_our12 <- accuracy_helper(asr_our12[,1:2],true_phy_info_new$State,tip_states)

pars <- c(1,  1,
          .045, .045,
          .7,  .5,
          .05,
          .05)

asr_our13 <- saasi(pars,phy,q_matrix)
acc_our13 <- accuracy_helper(asr_our13[,1:2],true_phy_info_new$State,tip_states)

pars <- c(1,  1,
          .045, .045,
          .75,  .5,
          .05,
          .05)

asr_our14 <- saasi(pars,phy,q_matrix)
acc_our14 <- accuracy_helper(asr_our14[,1:2],true_phy_info_new$State,tip_states)

pars <- c(1,  1,
          .045, .045,
          .8,  .5,
          .05,
          .05)

asr_our15 <- saasi(pars,phy,q_matrix)
acc_our15 <- accuracy_helper(asr_our15[,1:2],true_phy_info_new$State,tip_states)



pars <- c(1,  1,
          .045, .045,
          1,  .5,
          .05,
          .05)

asr_our16 <- saasi(pars,phy,q_matrix)
acc_our16 <- accuracy_helper(asr_our16[,1:2],true_phy_info_new$State,tip_states)


df <- data.frame(x = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,2),
                 y = c(acc_our1$acc,acc_our18$acc,
                       acc_our2$acc,acc_our3$acc,acc_our4$acc,
                       acc_our5$acc,acc_our6$acc,acc_our7$acc,
                       acc_our8$acc,acc_our9$acc,acc_our10$acc,
                       acc_our11$acc,acc_our12$acc,acc_our13$acc,
                       acc_our14$acc,acc_our15$acc,acc_our16$acc))

df$interesting <- "guess"
df$interesting[1] <- "true"
p_prop<-ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color=interesting),size =6) +labs(x = "Sampling Ratio", y = "Accuracy")+
  theme_light()+ ggtitle("A: Absolute accuracy")+
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15),
        legend.title=element_blank(),legend.position = "none")+
  ylim(0.0, 1.0)
p_prop
#ggsave("proportion_scen.pdf", p_prop,width = 210, height = 297, units = "mm")

df <- data.frame(x = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,2),
                 y = c(acc_our1$prob_acc,acc_our18$prob_acc,
                       acc_our2$prob_acc,acc_our3$prob_acc,acc_our4$prob_acc,
                       acc_our5$prob_acc,acc_our6$prob_acc,acc_our7$prob_acc,
                       acc_our8$prob_acc,acc_our9$prob_acc,acc_our10$prob_acc,
                       acc_our11$prob_acc,acc_our12$prob_acc,acc_our13$prob_acc,
                       acc_our14$prob_acc,acc_our15$prob_acc,acc_our16$prob_acc))

df$interesting <- "guess"
df$interesting[1] <- "true"
p_prop2<-ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color=interesting),size =6) +labs(x = "Sampling Ratio", y = "Accuracy")+
  theme_light()+ ggtitle("B: Probability accuracy")+
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15),
        legend.title=element_blank(),legend.position = "none")+
  ylim(0.0, 1.0)
p_prop2

ggarrange(p_prop,p_prop2,nrow=1)

ggsave(file="figs5_updated.pdf", width = 210, height = 297, units = "mm")





