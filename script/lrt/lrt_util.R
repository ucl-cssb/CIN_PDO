library(ggtree)
library(ape)
library(treeio)
library(phytools)
library(tidyverse)



# estimate lambda based on a set of subtrees with the same birth rate
estimate_lambda <- function(subtrees, use_root = T){
  tlen = 0
  nbirth = 0
  for(i in 1:length(subtrees)){
    #i = 1
    phy = subtrees[[i]]
    #print(i)
    #print(phy)
    X <- sum(phy$edge.length)
    tlen = tlen + X
    nb.node <- phy$Nnode - 1

    tname = names(subtrees)[i]
    # only applied on old trees and simulated trees
    if(!use_root && !startsWith(tname, "early")){
      edges = as.data.frame(phy$edge)
      edges$id = seq(1:nrow(edges))
      root = phy$Nnode + 2
      ids = edges %>% filter(V1 == root)
      rlen = sum(phy$edge.length[ids$id])
      tlen = tlen - rlen
      nb.node = nb.node - 2
    }

    nbirth = nbirth + nb.node
  }
  lambda = nbirth / tlen

  return(lambda)
}


# based on yule from ape, estimate lambda based on whole tree
get_lnl_h0_yule <- function(phy){
  #phy = t1
  X <- sum(phy$edge.length)
  nb.node <- phy$Nnode - 1
  lambda <- nb.node/X

  loglik <- -lambda * X + lfactorial(phy$Nnode) + nb.node * log(lambda)

  return (list(loglik = loglik, lambda = lambda))
}


# lambda estimated from known normal part of the tree
get_lnl_h0 <- function(phy, lambda){
  #phy = t1
  X <- sum(phy$edge.length)
  nb.node <- phy$Nnode - 1
  #lambda <- nb.node/X

  loglik <- -lambda * X + lfactorial(phy$Nnode) + nb.node * log(lambda)

  return (list(loglik = loglik, lambda = lambda))
}


# get branch lengths and #birth events before error subtree
get_elen_nbirth_norm <- function(phy, nodes_norm, is_lbl = F, exclude_enodes = T){
  blens_norm = c()
  edges = as.data.frame(phy$edge)
  edges$id = seq(1:nrow(edges))

  # Find all branches starting from normal nodes
  if(is_lbl){
    # convert nodes ID from label in the edge matrix
    nodes_sel = c()
    root = phy$Nnode + 2
    for(n in nodes_norm){
      #print(n)
      n_id = which(phy$node.label==n)
      #print(n_id)
      n_norm = n_id + root - 1
      nodes_sel = c(nodes_sel, n_norm)
    }
  }else{
    nodes_sel = nodes_norm
  }
  edges %>% filter(V1 %in% nodes_sel) -> norm_edges
  #print(norm_edges)
  # Find the IDs of these normal edges
  norm_edges %>% select(id) %>% unlist() -> ids
  # Sum up all these branches
  edges_norm = phy$edge.length[ids]
  #print(edges_norm)
  tlen = sum(edges_norm)

  # Find the end nodes of these normal edges
  ends_norm = norm_edges$V2
  # Find the number of birth events by finding nodes with outgoing edges
  # terminal edges -- There are no outgoing edges
  if(exclude_enodes){
    # Find the edges whose start nodes are in ends_norm, where terminal edges will be excluded
    # each normal edge is uniquely represented by its end node
    edges %>% filter(V1 %in% ends_norm) %>% select(V1) %>% unique() %>% unlist() %>% length() -> nnode
  }else{
    nnode = nrow(norm_edges)
  }
  return(list(tlen = tlen, nbirth = nnode, blens_norm = edges_norm))
}


# get partial log likelihood for branches before error
get_plnl_norm <- function(t1, lambda, nodes_norm, is_lbl = F){
  res = get_elen_nbirth_norm(t1, nodes_norm, is_lbl)
  tlen = res$tlen
  nnode = res$nbirth

  # not multiply by lambda for terminal edges
  lnlr = -lambda * tlen + nnode * log(lambda)

  return(lnlr)
}


# get partial log likelihood for subtrees without lfactorial
get_plnl <- function(phy, lambda){
  X <- sum(phy$edge.length)
  nb.node <- phy$Nnode - 1
  # print(nb.node)
  # print(X)
  s1 = -lambda * X
  s2 = nb.node * log(lambda)
  # print(s1)
  # print(s2)
  loglik <- s1 + s2

  #return (list(loglik = loglik, s1 = s1, s2 = s2))
  return(loglik)
}



# nodes_err: Node IDs with mis-segregation errors
# MLE of lambda1 in subtree with errors is #birth events / sum of branch lengths
get_lnl_subtree <- function(t1, nodes_err, lambda1, use_mle = F){
  subtrees = list()
  for(i in 1:length(nodes_err)){
    n1 = nodes_err[i]
    subtree1 = extract.clade(t1, n1)
    subtrees[[i]] = subtree1
  }
  if(use_mle){
    lambda1 = estimate_lambda(subtrees)
    #message("MLE of birth rate under H1: ", lambda1)
  }

  lambda = lambda1
  lnl_subtree = 0
  for(i in 1:length(subtrees)){
    subtree1 = subtrees[[i]]
    lnl1 = get_plnl(subtree1, lambda)
    lnl_subtree = lnl_subtree + lnl1
  }

  return(list(lnl_subtree = lnl_subtree, lambda1 = lambda1))
}



# based on yule from ape, assuming same lambda for all trees
# estimate lambda first as many trees are combined
get_lnl_h0_all_yule <- function(trees_all, use_root = T){
  loglik = 0
  tlen = 0
  nbirth = 0
  
  for(i in 1:length(trees_all)){
    #i = 1
    phy = trees_all[[i]]
    # print(i)
    # print(phy)
    X <- sum(phy$edge.length)
    tlen = tlen + X
    nb.node <- phy$Nnode - 1
    
    tname = names(trees_all)[i]
    # only applied on old trees
    if(!use_root && !startsWith(tname, "early")){
      edges = as.data.frame(phy$edge)
      edges$id = seq(1:nrow(edges))
      root = phy$Nnode + 2
      ids = edges %>% filter(V1 == root)
      rlen = sum(phy$edge.length[ids$id])
      tlen = tlen - rlen
      nb.node = nb.node - 2
    }
    
    nbirth = nbirth + nb.node
  }
  
  lambda = nbirth / tlen
  
  for(i in 1:length(trees_all)){
    #i = 1
    phy = trees_all[[i]]
    X <- sum(phy$edge.length)
    nb.node <- phy$Nnode - 1
    
    tname = names(trees_all)[i]
    if(!use_root && !startsWith(tname, "early")){
      edges = as.data.frame(phy$edge)
      edges$id = seq(1:nrow(edges))
      root = phy$Nnode + 2
      ids = edges %>% filter(V1 == root)
      rlen = sum(phy$edge.length[ids$id])
      X = X - rlen
      nb.node = nb.node - 2
    }
    
    lnl = -lambda * X + lfactorial(phy$Nnode) + nb.node * log(lambda)
    #print(lnl)
    loglik = loglik + lnl
  }
  
  return (list(loglik = loglik, lambda = lambda))
}

