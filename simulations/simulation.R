make.tree.bisse_test <- function(pars, max.taxa=Inf, max.t=Inf, x0,
                                 single.lineage=TRUE) {
  
  k <- (sqrt(4+4*length(pars))-2)/2
  #Warning message
  if ( !isTRUE(all.equal(k, as.integer(k))) )
    stop("Invalid parameter length: must be k(k+2) long")
  #check.pars.musse(pars, k)
  #if ( x0 < 1 || x0 > k )
  #stop("x0 must be an integer in [1,k]")
  
  # pars <- cbind(matrix(pars[seq_len(2*k)], k, 2), 
  #               matrix(pars[(2*k+1):(k*(k+1))], k, k-1, TRUE))
  
  #case with Sampling rate
  pars <- cbind(matrix(pars[seq_len(3*k)], k, 3),
                matrix(pars[(3*k+1):(k*(k+2))], k, k-1, TRUE))
  #print(pars)
  
  # to is a matrix representation of state changes according to the transition
  
  to <- matrix(unlist(lapply(1:k, function(i) (1:k)[-i])),
               k, k-1, TRUE)
  
  r.i <- rowSums(pars) #total rate for each state i 
  #print(pars)
  #print(r.i)
  
  #initialize all the parameters
  #remains the same
  extinct <- FALSE
  split   <- FALSE
  #Adding Sampling
  sam <- FALSE
  
  parent <- 0
  n.i <- rep(0, k) #return 0 0 0
  #print(k)
  len <- 0
  t <- 0
  hist <- list()
  
  #let the first starting states as x0
  #remains the same
  if ( single.lineage ) {
    states <- x0
    n.taxa <- lineages <- n.i[x0] <- 1 #n.taxa = 1, lineages = 1, n.i = 0 1 0 
    start <- 0
  } else {
    ##states <- rep(x0, 2)
    ##n.taxa <- lineages <- n.i[x0] <- 2
    stop("Nope.")
  }
  
  #while loop getting the time
  while ( n.taxa <= max.taxa && n.taxa > 0 ) {
    r.n <- r.i * n.i #r.n first iteration 0,0.295,0
    #print(r.n)
    r.tot <- sum(r.n) #total rate 
    #print(r.tot)
    dt <- rexp(1, r.tot) #using the total rate to get the time of next event 
    #print('dt')
    #print(dt)
    t <- t + dt #preceding the time
    
    if ( t > max.t ) { #if greater than max time, brake the while loop and return
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }
    len[lineages] <- len[lineages] + dt
    #print('here')
    #print(lineages)
    #print(len[lineages])
    
    ## Pick a state, and then a lineage in that state, to modify, then
    ## an event type: 1: speciation, 2: extinction, >2: char change.
    
    # with sampling TODO: ADDING A NEW EVENT TYPE: SAMPLING
    # type: 1: speciation, 2: extinction, 3: sampling, >3: char change
    # sampling assumes the sample is tested and hence become a taxa
    # in the phylogeny, this allows the phylogeny be non-ultrametric
    #print('passing here')
    state <- sample(k, 1, FALSE, r.n/r.tot) #from k states, pick one, without replacement, with the given probability
    #print('state')
    #print(state)
    # print('r.n/r.tot')
    # print(r.n/r.tot)
    #print('passing as well')
    j <- sample(n.i[state], 1)
    #print('passing same here')
    # print('j')
    # print(j)
    # print('n.i')
    # print(n.i)
    # print(n.i[state])
    # print('lineages')
    # print(lineages)
    # print('states')
    # print(states)
    # print('states[lineages]')
    # print(states[lineages])
    # print('lineages[states[lineages]==state]')
    # print(lineages[states[lineages]==state])
    # print(lineages[states[lineages]==state][j])
    lineage <- lineages[states[lineages] == state][j]
    # print('lineage')
    # print(lineage)
    
    #TODO: check again if needs to change, don't think so..
    #The above step gets the lineage that can be altered
    #print('maybe here')
    #print(pars[state,])
    # type <- sample(k+1, 1, FALSE, pars[state,])
    #print('pars')
    #print(pars[state,])
    type <- sample(k+2, 1, FALSE, pars[state,])
    #print('?')
    #print('type')
    #print(type)
    
    if ( type == 1 ) { # Speciation:
      if ( n.taxa == max.taxa )
        break
      new.i <- length(extinct) + 1:2 #adding two new lineages
      split[lineage] <- TRUE #if the lineage speciate, convert them to TRUE
      extinct[new.i] <- split[new.i] <- FALSE #the lineage cannot extinct, convert them to FALSE
      #add sampling
      sam[new.i] <- FALSE
      states[new.i] <- state #the new speciated state should have the same state
      parent[new.i] <- lineage #discover their parent node
      start[new.i] <- t #start time for the speciated new lineages
      len[new.i] <- 0 #length 0
      n.i[state] <- n.i[state] + 1 #update number of lineages of state i
      n.taxa <- n.taxa + 1 #update number of taxa
      #lineages <- which(!split & !extinct) #update lineages, should not be split nor extinct
      #TODO should change here if consider sampling
      lineages <- which(!split & !extinct & !sam)
      #print('lineage after speciation')
      #print(lineage)
    } else if ( type == 2 ) { # Extinction:
      extinct[lineage] <- TRUE # no speciation event, update the lineage to extinct
      #lineages <- which(!split & !extinct) 
      #TODO add sampling
      lineages <- which(!split & !extinct & !sam)
      n.i[state] <- n.i[state] - 1
      n.taxa <- n.taxa - 1
      # TODO ADDING A SAMPLING EVENT
      #print('lineage after extinction')
      #print(lineage)
    } else if (type == 3 ) { # Sampling:
      sam[lineage] <- TRUE
      #extinct[lineage] <- TRUE
      n.taxa <- n.taxa-1
      n.i[state] <- n.i[state] - 1
      lineages <- which(!split & !extinct & !sam)
      #print('lineage after sampling')
      #print(lineage)
    }
    else {
      #states[lineage] <- state.new <- to[state,type - 2]
      states[lineage] <- state.new <- to[state,type - 3]
      n.i[c(state.new, state)] <- n.i[c(state.new, state)] + c(1,-1)
      hist[[length(hist)+1]] <- c(lineage, t, state, state.new)
      #print('lineage after transition')
      #print(lineage)
    }
  }
  # info <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
  #                    start=start, state=states, extinct=extinct,
  #                    split=split)
  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
                     start=start, state=states, extinct=extinct,
                     split=split, sam=sam)
  
  hist <- as.data.frame(do.call(rbind, hist))
  if ( nrow(hist) == 0 )
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0
  
  attr(info, "t") <- t
  attr(info, "hist") <- hist
  info
}


me.to.ape.bisse <- function(x, root.state) {
  if ( nrow(x) == 0 )
    return(NULL)
  Nnode <- sum(!x$split) - 1
  n.tips <- sum(!x$split)
  
  x$idx2 <- NA
  x$idx2[!x$split] <- 1:n.tips
  x$idx2[ x$split] <- order(x$idx[x$split]) + n.tips + 1
  
  i <- match(x$parent, x$idx)
  x$parent2 <- x$idx2[i]
  x$parent2[is.na(x$parent2)] <- n.tips + 1
  
  tip.label <- ifelse(subset(x, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)
  
  x$name <- NA
  x$name[!x$split] <- tip.label
  ## More useful, but I don't want to clobber anything...
  x$name2 <- c(tip.label, node.label)[x$idx2]
  
  tip.state <- x$state[match(1:n.tips, x$idx2)]
  names(tip.state) <- tip.label
  
  node.state <- x$state[match(1:Nnode + n.tips, x$idx2)]
  names(node.state) <- node.label
  node.state["nd1"] <- root.state
  
  hist <- attr(x, "hist")
  if ( !is.null(hist) ) {
    hist$idx2 <- x$idx2[match(hist$idx, x$idx)]
    hist$name2 <- x$name2[match(hist$idx, x$idx)]
    if ( nrow(hist) > 0 )
      hist <- hist[order(hist$idx2),]
  }
  
  phy <- reorder(structure(list(edge=cbind(x$parent2, x$idx2),
                                Nnode=Nnode,
                                tip.label=tip.label,
                                tip.state=tip.state,
                                node.label=node.label,
                                node.state=node.state,
                                edge.length=x$len,
                                orig=x,
                                hist=hist),
                           class="phylo"))
  phy$edge.state <- x$state[match(phy$edge[,2], x$idx2)]
  phy
}



tree.sim <- function(pars, max.taxa=Inf, max.t=Inf,
                            include.extinct=FALSE, x0=NA) {
  # k <- (sqrt(1 + 4*length(pars))-1)/2
  # if ( !isTRUE(all.equal(k, as.integer(k))) )
  #   stop("Invalid parameter length: must be k(k+1) long")
  # if ( length(x0) != 1 || is.na(x0) || x0 < 1 || x0 > k )
  #   stop("Invalid root state")
  k <- NULL
  while(is.null(k)){
    info <- make.tree.bisse_test(pars, max.taxa, max.t, x0)
    k <- me.to.ape.bisse(info[-1,], info$state[1])
  }
  #print(info)
  ## This is a bit simple-minded, as extinction that erodes the root
  ## node will affect this once pruned.  However, it is always correct
  ## for trees without extinct taxa pruned.
  ## attr(info, "t") <- attr(info, "t") - info$len[1]
  phy <- me.to.ape.bisse(info[-1,], info$state[1])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}


