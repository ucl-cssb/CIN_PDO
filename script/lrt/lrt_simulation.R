library(ggpubr)
library(snpStats)

# assuming working directory is "script"

source("./lrt/lrt_util.R")

dir = getwd()  # output will be stored in dir


# used in test_sim_batch_fixt, for box plots of branch lengths
test = "wilcox.test"
pstat = stat_compare_means(method = test, label.x.npc = 0.2)
#theme =  theme_gray(base_size = 12) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "top")
pdot = geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)


################# LRT on simulated tree ##################
# LRT of yule tree model with different birth rate
test_lr_sim <- function(res_h0, res_h1){
  lnl_h0_mle = res_h0$loglik
  lnl_h1_mle = res_h1$loglik

  #near(lnl_h0_mle, lnl_h1_mle)

  lr = -2 * (lnl_h0_mle - lnl_h1_mle)
  #print(lr)

  # test significance with chi-square distribution
  #print(lr > cutoff1)
  preject = 1 - pchisq(lr, 1)
  #print(preject)

  # if(preject < 0.05){
  #   message("H0 rejected")
  # }else{
  #   message("H0 not rejected")
  # }

  df = data.frame(lambda0=res_h0$lambda, loglik0=res_h0$loglik, lambda0_1=res_h1$lambda0, lambda1 = res_h1$lambda, loglik1 = res_h1$loglik, test_stat = lr, pvalue = preject)

  return(df)
}


# to find what differences expect to detect given the sample size
# not use root branches as they are just artificial to make Yule tree
get_lnl_h1_simtree <- function(t1, lambda, nodes_norm2 = c(22), nodes_err = c(31), use_mle = T, use_root = F){
  lnl_t1 = 0

  #print(t1)
  root = t1$Nnode + 2
  lnlr = 0
  if(use_root){
    nodes_norm = c(root)
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
  }

  lnl_t1 = lnl_t1 + lnlr + lfactorial(t1$Nnode)

  # normal subtree
  for(j in 1:length(nodes_norm2)){
    n1 = nodes_norm2[j]
    subtree1 = extract.clade(t1, n1)
    #plot_tree(subtree1, "", F)
    lnl1 = get_plnl(subtree1, lambda)
    lnl_t1 = lnl_t1 + lnl1
  }

  # error subtree
  res_subtree = get_lnl_subtree(t1, nodes_err, lambda, use_mle)
  lnl_subtree = res_subtree$lnl_subtree
  lnl_t1 = lnl_t1 + lnl_subtree

  if(use_mle){
    lambda1 = res_subtree$lambda1
    return(list(loglik = lnl_t1, lambda = lambda1))
  }else{
    return(lnl_t1)
  }
}


test_t1_sim <- function(lambda = 0.5, t = 5, ntest = 1, use_root = T){
  res_sim = data.frame()
  for(i in 1:ntest){
    t1 = rlineage(lambda, 0, t)
    #print(str(t1))
    #plot_tree(t1, "", T)

    edges = as.data.frame(t1$edge)
    edges$id = seq(1:nrow(edges))
    root = t1$Nnode + 2
    ids = edges %>% filter(V1 == root)

    # assume one subtree is normal and the other is mis-seg
    nodes_err = ids[1,]$V2
    nodes_norm2 = ids[2,]$V2
    # either one such node is the tip
    if(!nodes_norm2 %in% edges$V1 || !nodes_err %in% edges$V1){
      #return(data.frame())
      message("single normal or error node!")
      next
    }

    # ensure there are at least 3 nodes in the error tree to estimate birth rate
    terr = extract.clade(t1, nodes_err[1])
    tnorm = extract.clade(t1, nodes_norm2[1])
    tsize1 = terr$Nnode
    tsize2 = tnorm$Nnode

    if(tsize1 < 2 || tsize2 < 2){
      message("small subtree!")
      next
    }

    # LRT on this tree
    lambda0 = estimate_lambda(c(t1), use_root)
    res_h0 = get_lnl_h0(t1, lambda0)
    #res_h0
    lambda0e = estimate_lambda(c(tnorm), use_root)
    res_h1 = get_lnl_h1_simtree(t1, lambda0e, nodes_norm2, nodes_err, T, use_root)
    res_h1$lambda0 = lambda0e
    df = test_lr_sim(res_h0, res_h1)
    #df
    res_sim = rbind(res_sim, df)
    #return(df)
  }

  return(res_sim)
}



sim_tree_with_error_fixt <- function(b0, b1, t){
  # top subtree (equal branch lengths)
  simtree20 = pbtree(b0, 0, n = 2)
  # normal subtree
  simtree21 = rlineage(b0, 0, t)
  # error subtree
  simtree22 = rlineage(b1, 0, t)

  simtree2_1 = bind.tree(simtree20, simtree21, 1)
  simtree2 = bind.tree(simtree2_1, simtree22, 1)
  #yule(simtree2)

  ntip = simtree21$Nnode + simtree22$Nnode + 2
  root = ntip + 1
  nnorm = simtree21$Nnode + 1
  nodes_err = c(nnorm+root)
  nodes_norm2 = c(ntip+2)
  # nodes_norm2
  # nodes_err
  #
  # p = plot_tree(simtree20, "top tree", T)
  # print(p)
  # p = plot_tree(simtree21, "normal tree", T)
  # print(p)
  # p = plot_tree(simtree22, "error tree", T)
  # print(p)
  # p = plot_tree(simtree2_1, "half joint tree", T)
  # print(p)
  # p = plot_tree(simtree2, "joint tree", T)
  # print(p)

  return(list(simtree2 = simtree2, simtree21 = simtree21, simtree22 = simtree22, nodes_norm2 = nodes_norm2, nodes_err = nodes_err))
}


find_internal_edges <- function(t1){
  edges = as.data.frame(t1$edge)
  edges$id = seq(1:nrow(edges))
  ntip = t1$Nnode + 1
  edges %>% filter(V2 > ntip) %>% select(id) %>% unlist() -> ie_ids

  blens = t1$edge.length[ie_ids]

  return(blens)
}

# root branches are created for connection of two different subtrees with fixed lengths, not so informative
test_sim_batch_fixt <- function(b0, b1s, times, use_root = F){
  res_sim = data.frame()
  blen_all = data.frame()
  blen_plots = list()

  i = 0
  for(t in times){
    # b0 = 0.5
    #t = 3
    message("end time: ", t)
    for(b1 in b1s){
      #b1 = 0.3
      #message("birth rate under H1: ", b1)
      i = i + 1  # used to index plot

      res = sim_tree_with_error_fixt(b0, b1, t)
      simtree21 = res$simtree21
      simtree22 = res$simtree22
      simtree2 = res$simtree2
      nodes_norm2 = res$nodes_norm2
      nodes_err = res$nodes_err

      # check the branch length distributions
      ntip_norm = simtree21$Nnode + 1
      ntip_err = simtree22$Nnode + 1

      # exclude terminal branches
      blens_norm = find_internal_edges(simtree21)
      blens_err = find_internal_edges(simtree22)
      if(length(blens_norm) > 0){
        df_blens1 = data.frame(blen = blens_norm, type = "N")
      }else{
        df_blens1 = data.frame()
      }
      if(length(blens_err) > 0){
        df_blens2 = data.frame(blen = blens_err, type = "E")
      }else{
        df_blens2 = data.frame()
      }

      df_blens = rbind(df_blens1, df_blens2)
      blen_all = rbind(blen_all, df_blens)

      nnorm = 0
      nerr = 0
      if(length(blens_norm) > 0 && length(blens_err) > 0){
        df_blens %>% group_by(type) %>% tally() -> bcount
        nerr = bcount[bcount$type=="E",]$n
        nnorm = bcount[bcount$type=="N",]$n

        pblen = ggplot(blen_all, aes(x=type, y=blen, color = type)) + geom_boxplot() + ylab("cell cycle duration")  + xlab("") + theme_pubr() + pstat + scale_color_discrete(name="Division type", breaks = c("E", "N"), labels = c(paste0("After ancestral missegregation (", nerr, ")"), paste0("Normal (", nnorm, ")")))
        blen_plots[[i]] = pblen
      }

      # lambda0 = estimate_lambda(c(simtree2))
      # #lambda0
      # res_h0 = get_lnl_h0(simtree2, lambda0)
      # need to exclude top 2 branches
      tlist0 = list()
      tlist0[["t0"]] = simtree2
      res_h0 = get_lnl_h0_all_yule(tlist0, use_root)
      # res_h0

      #get_lnl_h1_simtree(simtree2, lambda0, nodes_norm2, nodes_err, F)
      # need to exclude cases when there is no data to estimate lambda0
      lambda0 = estimate_lambda(c(simtree21), T)
      # estimate lambda1 based on error subtrees extracted by nodes_err
      res_h1 = get_lnl_h1_simtree(simtree2, lambda0, nodes_norm2, nodes_err, T, use_root)
      res_h1$lambda0 = lambda0
      df = test_lr_sim(res_h0, res_h1)
      df$time = t
      df$birth_rate0 = b0
      df$birth_rate1 = b1
      df$id = i
      df$s = df$lambda1 / df$lambda0_1  - 1
      df$nnorm = nnorm
      df$nerr = nerr
      df$ntip_norm = ntip_norm
      df$ntip_err = ntip_err
      #res_h1
      res_sim = rbind(res_sim, df)

    }
  }

  return(list(res_sim=res_sim, blen_all = blen_all, blen_plots = blen_plots))
}

# Get log likelihood under H1
# estimate lambda0 based on error subtrees
get_lnl_h1_simtree_all <- function(strees, strees_norm, strees_err, lambda0, lambda1, use_root = F){
  lnl_decom = 0

  for(j in 1:length(strees)){
    # constant part for each tree
    lnl_t1 = 0
    t1 = strees[[j]]
    #print(t1)
    root = t1$Nnode + 2
    lnlr = 0
    if(use_root){
      nodes_norm = c(root)
      lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    }
    lnl_t1 = lnl_t1 + lnlr + lfactorial(t1$Nnode)

    lnl_decom = lnl_decom + lnl_t1
  }

  # normal subtree
  for(j in 1:length(strees_norm)){
    subtree1 = strees_norm[[j]]
    #plot_tree(subtree1)
    lnl1 = get_plnl(subtree1, lambda0)
    #print(lnl1)
    lnl_decom = lnl_decom + lnl1
  }

  # error subtree
  for(j in 1:length(strees_err)){
    subtree1 = strees_err[[j]]
    #plot_tree(subtree1)
    lnl1 = get_plnl(subtree1, lambda1)
    #print(lnl1)
    lnl_decom = lnl_decom + lnl1
  }
  #print(lnl_decom)

  return (list(loglik = lnl_decom, lambda = lambda1))
}


# check estimations based on n combined trees (n = length(b1s))
# root branches are created for connection of two different subtrees, not so informative
test_sim_merged_fixt <- function(b0, b1s, times, use_root = F){
  res_sim = data.frame()
  blen_plots = list()
  strees_all = list()  # store all the trees to check branch length
  strees_norm_all = list()
  strees_err_all = list()

  j = 0 # for t
  for(t in times){
    # b0 = 0.5
    #t = 3
    j = j + 1

    blen_all = data.frame()
    strees_norm = list()
    strees_err = list()
    strees = list()
    ntip_norm_all = 0
    ntip_err_all = 0
    i = 0 # for i

    message("run ", j)
    for(b1 in b1s){
      i = i + 1

      #b1 = 0.3
      #message("birth rate under H1: ", b1)

      res = sim_tree_with_error_fixt(b0, b1, t)
      simtree21 = res$simtree21
      simtree22 = res$simtree22
      simtree2 = res$simtree2
      nodes_norm2 = res$nodes_norm2
      nodes_err = res$nodes_err

      tname = paste0("t", i)
      strees_norm[[tname]] = simtree21
      strees_err[[tname]] = simtree22
      strees[[tname]] = simtree2

      # find ids of terminal edges
      # plot_tree(simtree21)
      # plot_tree(simtree22)
      ntip_norm = simtree21$Nnode + 1
      ntip_norm_all = ntip_norm_all + ntip_norm
      ntip_err = simtree22$Nnode + 1
      ntip_err_all = ntip_err_all + ntip_err

      blens_norm = find_internal_edges(simtree21)
      blens_err = find_internal_edges(simtree22)
      if(length(blens_norm) > 0){
        df_blens1 = data.frame(blen = blens_norm, type = "N")
      }else{
        df_blens1 = data.frame()
      }
      if(length(blens_err) > 0){
        df_blens2 = data.frame(blen = blens_err, type = "E")
      }else{
        df_blens2 = data.frame()
      }

      df_blens = rbind(df_blens1, df_blens2)
      blen_all = rbind(blen_all, df_blens)

    }

    # check the joint branch length distributions
    blen_all %>% group_by(type) %>% tally() -> bcount
    bcount %>% filter(type == "E") %>% nrow() -> ne
    bcount %>% filter(type == "N") %>% nrow() -> nn

    if(ne > 0 && nn > 0){
      nerr = bcount[bcount$type=="E",]$n
      nnorm = bcount[bcount$type=="N",]$n
      # exclude terminal branches
      pblen = ggplot(blen_all, aes(x=type, y=blen, color = type)) + geom_boxplot() + ylab("cell cycle duration")  + xlab("") + theme + pstat + scale_color_discrete(name="Division type", breaks = c("E", "N"), labels = c(paste0("After ancestral missegregation (", nerr, ")"), paste0("Normal (", nnorm, ")")))
      #plot(pblen)
      blen_plots[[j]] = pblen
    }

    #length(strees)
    # a set of trees whose roots are all excluded for estimations
    estimate_lambda(strees, use_root)
    res_h0 = get_lnl_h0_all_yule(strees, use_root)
    #res_h0
    # should have the same value as res_h0
    #get_lnl_h1_simtree_all(strees, strees_norm, strees_err, res_h0$lambda, res_h0$lambda, use_root)


    #get_lnl_h1_simtree(simtree2, lambda0, nodes_norm2, nodes_err, F)
    # need to exclude cases when there is no data to estimate lambda0
    lambda0 = estimate_lambda(strees_norm, T)
    lambda1 = estimate_lambda(strees_err, T)
    res_h1 = get_lnl_h1_simtree_all(strees, strees_norm, strees_err, lambda0, lambda1, use_root)
    res_h1$lambda0 = lambda0
    df = test_lr_sim(res_h0, res_h1)
    df$time = t
    df$birth_rate0 = b0
    df$birth_rate1 = b1
    df$id = j
    df$s = df$lambda1 / df$lambda0_1  - 1
    df$nnorm = nnorm
    df$nerr = nerr
    df$ntip_norm_all = ntip_norm_all
    df$ntip_err_all = ntip_err_all
    #res_h1
    res_sim = rbind(res_sim, df)

    strees_all[[j]] = strees
    strees_norm_all[[j]] = strees_norm
    strees_err_all[[j]] = strees_err
  }

  return(list(res_sim=res_sim, blen_plots = blen_plots, strees_all = strees_all, strees_norm_all = strees_norm_all, strees_err_all = strees_err_all))
}



############### check null distribution (Suppl Fig 4) #################
# check H0 for single trees
# when t < 10, observed values seem smaller than expected value
ntest = 1000
lambda = 0.5
t = 10
use_root = T
# simulate a yule tree given lambda and t
res_sim = test_t1_sim(lambda, t, ntest)
res_sim %>% filter(pvalue < 0.05) %>% nrow()
summary(res_sim)
nres = nrow(res_sim)
df = 1

# have to save manually
# may slightly different from Suppl Fig 7 in the doc due to random seed
qq.chisq(res_sim$test_stat, 1, ylab=paste0("Observed value (", nres, " trees with t = ", t, ")"), main=paste("Expected distribution: chi-squared (", df," df)"), xlab="Expected value of test statistics", sub = "")



############## get all powers for suppl doc (Suppl Table 3) ################
# Assuming H1, estimate power as fraction of simulated trees correctly inferred s < 0

# Do this for single tree, four trees combined and 60 trees combined
get_power_single <- function(b0, b1, t, nrun){
  b1s=rep(b1, nrun)
  times=c(t)
  res_all = data.frame()
  # not use root
  res = test_sim_batch_fixt(b0, b1s, times, F)
  res_sim = res$res_sim
  res_all = rbind(res_all, res_sim)
  # nrow(res_sim)
  # summary(res_sim)
  res_sim = res_sim %>% drop_na()
  #nrow(res_sim)
  #summary(res_sim)

  res_sim %>% filter(pvalue < 0.05) %>% arrange(pvalue) %>% select(id) %>% unlist() -> sig_ids
  res_sim %>% filter(pvalue >= 0.05) %>% select(id) %>% unlist() -> nsig_ids
  res_sim %>% filter(s > 0) %>% select(id) %>% unlist() -> pos_ids
  res_sim %>% filter(pvalue < 0.05) %>% filter(s < 0)  %>% arrange(pvalue) %>% select(id) %>% unlist() -> neg_sig_ids

  power = length(neg_sig_ids) / nrow(res_sim)
  #power

  df = data.frame(b0 = b0, b1 = b1, t = t, nrun = nrun, num_merged_trees = 1, power = power, nsucc = nrow(res_sim),  nsig = length(sig_ids), nnsig = length(nsig_ids) , nsigs_neg = length(neg_sig_ids), npos = length(pos_ids))

  return(list(df = df, res_all = res_all))
}


get_power_merged <- function(nm, b0, b1, t, nrun){
  b1s=rep(b1, nm)
  times=rep(t, nrun)
  res_all = data.frame()

  res = test_sim_merged_fixt(b0, b1s, times, F)
  res_sim = res$res_sim
  res_all = rbind(res_all, res_sim)
  res_sim = res_sim %>% drop_na()
  #nrow(res_sim)
  #summary(res_sim)

  res_sim %>% filter(pvalue < 0.05)  -> sig_ids
  res_sim %>% filter(pvalue >= 0.05) -> nsig_ids
  res_sim %>% filter(s > 0) -> pos_ids
  res_sim %>% filter(pvalue < 0.05) %>% filter(s < 0) -> neg_sig_ids
  power = nrow(neg_sig_ids) / nrow(res_sim)
  #power

  df = data.frame(b0 = b0, b1 = b1, t = t, nrun = nrun, num_merged_trees = nm, power = power, nsucc = nrow(res_sim),  nsig = nrow(sig_ids), nnsig = nrow(nsig_ids) , nsigs_neg = nrow(neg_sig_ids), npos = nrow(pos_ids))

  return(list(df = df, res_all = res_all))
}


# simulate a yule tree with mixed birth rates
b0 = 0.65
b1 = 0.4
nrun = 100
rep = 10

power_all = data.frame()

# single tree
power_single = data.frame()
res_single = data.frame()
for(t in c(3, 5)){
  for(i in 1:rep){
    res1 = get_power_single(b0, b1, t, nrun)
    p1 = res1$df
    p1$rep = i
    r1 = res1$res_all
    r1$rep = i
    power_single = rbind(power_single, p1)
    res_single = rbind(res_single, r1)
  }
}
power_single %>% group_by(num_merged_trees, t) %>% summarise(mean_power = mean(power))


# merged trees of different size
power_merged = data.frame()
res_merged = data.frame()
for(t in c(3, 5)){
  for(nm in c(4, 60)){
    for(i in 1:rep){
      res1 = get_power_merged(nm, b0, b1, t, nrun)
      p1 = res1$df
      p1$rep = i
      r1 = res1$res_all
      r1$rep = i
      power_merged = rbind(power_merged, p1)
      res_merged = rbind(res_merged, r1)
    }
  }
}
power_merged %>% group_by(num_merged_trees, t) %>% summarise(mean_power = mean(power))


power_all = rbind(power_single, power_merged)
power_all %>% group_by(num_merged_trees, t) %>% summarise(mean_power = mean(power)) -> sum_pow
sum_pow

fout = file.path(dir, paste0("sum_power_b0", b0, "_b1", b1, "_nrun", nrun, "_rep", rep,".tsv"))
write_tsv(sum_pow, fout)



################## check Effect of apoptosis on estimates (last paragraph of Suppl doc) #######################
# the real death rate should be very low
# count #trees with death events given different death rate (under H0)
get_ndeath <- function(dtimes, nrun, b0, t){
  df_nd = data.frame()
  for(d in dtimes){
    ndeath = 0
    for(r in 1:nrun){
      t1 = rlineage(b0, d, t)
      # find #deaths by tip depths
      edepths = node.depth.edgelength(t1)
      tdepth = max(node.depth.edgelength(t1))
      # nodes with the same max depth, should equal to #leaves
      ned = length(which(edepths==tdepth))
      nleave = t1$Nnode + 1
      extra = ned - nleave
      if(extra != 0){
        message("#leaves have same depth: ", extra, " at death rate ", d)
        #pt = plot_tree(t1, "", T)
        #print(pt)
        ndeath = ndeath + 1
      }
    }
    df = data.frame(drate = d, ndeath = ndeath, brate = b0, time = t, nrun = nrun)
    df_nd = rbind(df_nd, df)
  }

  return(df_nd)
}


# results may slightly different from those in Suppl doc due to random seed
dtimes = c(0.05, 0.04, 0.03, 0.02, 0.01)
nrep = 100
b0 = 0.5


# for new trees with shorter time
nrun = 60
t = 3
df_nd_all = data.frame()
for(i in 1:nrep){
  df_nd = get_ndeath(dtimes, nrun, b0, t)
  df_nd$rep = i
  df_nd_all = rbind(df_nd_all, df_nd)
}
df_nd_all

df_nd_all %>% group_by(drate) %>% summarise(mean_ndeath = mean(ndeath)) -> sum_ndeath
sum_ndeath

fout = file.path(dir, paste0("sum_sim_ndeath_t", t, "_nrun", nrun, "_b0", b, "_nrep", nrep, ".tsv"))
write_tsv(sum_ndeath, fout)


# for old trees with longer time
nrun = 4
t = 5
df_nd_all_old = data.frame()
for(i in 1:nrep){
  df_nd = get_ndeath(dtimes, nrun, b0, t)
  df_nd$rep = i
  df_nd_all_old = rbind(df_nd_all_old, df_nd)
}
df_nd_all_old
df_nd_all_old %>% group_by(drate) %>% summarise(mean_ndeath = mean(ndeath)) -> sum_ndeath_old
sum_ndeath_old

fout = file.path(dir, paste0("sum_sim_ndeath_t", t, "_nrun", nrun, "_b0", b, "_nrep", nrep, ".tsv"))
write_tsv(sum_ndeath_old, fout)


####################### check exponential distribution of BLs (Suppl Fig. 7)  ######################
# simulated trees under H0 fit lambda (estimated vs real)
# see how strengths of the selection affects the empirical distribution
# simulated  individual trees under H1 fit lambda (estimated vs real)
get_blen_distr <- function(b0, b1, t){
  res = sim_tree_with_error_fixt(b0, b1, t)
  simtree21 = res$simtree21
  simtree22 = res$simtree22
  simtree2 = res$simtree2
  nodes_norm2 = res$nodes_norm2
  nodes_err = res$nodes_err
  use_root = F

  # need to exclude root
  res_h0 = get_lnl_h0_yule(simtree2)
  res_h0

  lambda0 = estimate_lambda(c(simtree21), T)
  res_h1 = get_lnl_h1_simtree(simtree2, lambda0, nodes_norm2, nodes_err, T, use_root)
  res_h1$lambda0 = lambda0
  res_h1

  main = paste0("Normal subtree, simulation until time ", t)
  # plot simulated branch length distribution with estimated or real rate
  elens_n = data.frame(blens = simtree21$edge.length)
  pn = ggplot(elens_n, aes(blens)) + geom_histogram() + geom_vline(xintercept = b0, color = "red", linetype = "dashed") + geom_vline(xintercept = res_h1$lambda0, color = "blue", linetype = "dashed") + ggtitle(main)
  pn

  main = paste0("Error subtree, simulation until time ", t)
  # plot simulated branch length distribution with estimated or real rate
  elens_e = data.frame(blens = simtree22$edge.length)
  pe = ggplot(elens_e, aes(blens)) + geom_histogram() + geom_vline(xintercept = b1, color = "red", linetype = "dashed") + geom_vline(xintercept = res_h1$lambda, color = "blue", linetype = "dashed") + ggtitle(main)
  pe

  pm = ggarrange(pn, pe, nrow = 1)

  return(pm)
}


# Write out the cdf for the exponential and calculate the probability between the bins, given the value of lambda.
get_exp_per_bin <- function(br, lambda){
  df_exp = data.frame()
  for(i in 2:length(br)){
    #print(i)
    fb = pexp(br[i], lambda) - pexp(br[i-1], lambda)
    #print(fb)
    bsize = br[i] - br[i-1]
    loc = br[i-1] + (bsize)/2
    #print(loc)
    #  / bsize
    df = data.frame(bin = loc, val = fb)
    df_exp = rbind(df_exp, df)
  }

  return(df_exp)
}


# for joint trees
# under H1, simulate 60 trees and check BL distribution
get_blen_distr_joint <- function(b0, b1, t, ntree, nbins = 10, ylim = 0.5, xlim = 4, add_vline = F){
  bsize = xlim / nbins
  br = seq(0, xlim, bsize)

  b1s = rep(b1, ntree)

  blen_all = data.frame()
  strees_norm = list()
  strees_err = list()
  strees = list()
  i = 0 # for i

  for(b1 in b1s){
    i = i + 1

    #b1 = 0.3
    #message("birth rate under H1: ", b1)

    res = sim_tree_with_error_fixt(b0, b1, t)
    simtree21 = res$simtree21
    simtree22 = res$simtree22
    simtree2 = res$simtree2
    nodes_norm2 = res$nodes_norm2
    nodes_err = res$nodes_err

    tname = paste0("t", i)
    strees_norm[[tname]] = simtree21
    strees_err[[tname]] = simtree22
    strees[[tname]] = simtree2

    blens_norm = find_internal_edges(simtree21)
    blens_err = find_internal_edges(simtree22)
    if(length(blens_norm) > 0){
      df_blens1 = data.frame(blen = blens_norm, type = "N")
    }else{
      df_blens1 = data.frame()
    }
    if(length(blens_err) > 0){
      df_blens2 = data.frame(blen = blens_err, type = "E")
    }else{
      df_blens2 = data.frame()
    }

    df_blens = rbind(df_blens1, df_blens2)
    blen_all = rbind(blen_all, df_blens)
  }

  #length(strees)
  #estimate_lambda(strees, use_root)
  res_h0 = get_lnl_h0_all_yule(strees, use_root)
  res_h0

  #get_lnl_h1_simtree(simtree2, lambda0, nodes_norm2, nodes_err, F)
  # need to exclude cases when there is no data to estimate lambda0
  lambda0 = estimate_lambda(strees_norm, T)
  lambda1 = estimate_lambda(strees_err, T)
  res_h1 = get_lnl_h1_simtree_all(strees, strees_norm, strees_err, lambda0, lambda1, use_root)
  res_h1$lambda0 = lambda0

  # check the joint branch length distributions
  nblen = 0
  for(st in strees_norm){
    #print(st)
    nblen = nblen + nrow(st$edge)
  }
  #main = paste0("Normal subtree, ", nblen, " branches", "\nsimulation until time ", t, ", ", ntree, " trees", "\nestimated cell cycle duration ", sprintf("%.1f", 1 / res_h0$lambda), " (day)")
  main = paste0(ntree, " trees", ", ",  nblen, " branches", "\nsimulation until time ", t, "\nestimated cell cycle duration ", sprintf("%.1f", 1 / res_h0$lambda), " (day)")
  # plot simulated branch length distribution with estimated or real rate
  elens_n = blen_all %>% filter(type == "N")
  pn = ggplot(elens_n, aes(blen, y=stat(count) / sum(count))) + geom_histogram(breaks = br, fill = "lightblue", color = "black") + ggtitle(main) + theme_pubr() + ylab("frequency") + scale_y_continuous(limits = c(0, ylim)) + xlab("branch length before mis-segregation")
  if(add_vline){
    pn = pn + geom_vline(xintercept = b0, color = "red", linetype = "dashed") + geom_vline(xintercept = res_h1$lambda0, color = "blue", linetype = "dashed")
  }

  df_exp = get_exp_per_bin(br, 2 * b0)
  pn = pn + geom_point(data = df_exp, mapping = aes(x = bin, y=val)) + geom_line(data = df_exp, mapping = aes(x = bin, y=val), color = "red", linetype="dashed")
  df_exp_est = get_exp_per_bin(br, 2 * res_h1$lambda0)
  pn = pn + geom_point(data = df_exp_est, mapping = aes(x = bin, y=val)) + geom_line(data = df_exp_est, mapping = aes(x = bin, y=val), color = "blue", linetype="dashed")

  nblen = 0
  for(st in strees_err){
    #print(st)
    nblen = nblen + nrow(st$edge)
  }
  main = paste0(ntree, " trees", ", ",  nblen, " branches", "\nsimulation until time ", t, "\nestimated cell cycle duration ", sprintf("%.1f", 1 / res_h1$lambda), " (day)")
  # plot simulated branch length distribution with estimated or real rate
  elens_e = blen_all %>% filter(type == "E")
  pe = ggplot(elens_e, aes(blen, y=stat(count) / sum(count))) + geom_histogram(breaks = br, fill = "lightblue", color = "black")  + ggtitle(main) + theme_pubr() + ylab("frequency") + scale_y_continuous(limits = c(0, ylim)) + xlab("branch length after mis-segregation")
  if(add_vline){
    pe = pe + geom_vline(xintercept = b1, color = "red", linetype = "dashed") + geom_vline(xintercept = res_h1$lambda, color = "blue", linetype = "dashed")
  }
  df_exp = get_exp_per_bin(br, 2 * b1)
  pe = pe + geom_point(data = df_exp, mapping = aes(x = bin, y=val)) + geom_line(data = df_exp, mapping = aes(x = bin, y=val), color = "red", linetype="dashed")
  df_exp_est = get_exp_per_bin(br, 2 * res_h1$lambda)
  pe = pe + geom_point(data = df_exp_est, mapping = aes(x = bin, y=val)) + geom_line(data = df_exp_est, mapping = aes(x = bin, y=val), color = "blue", linetype="dashed")

  pm = ggarrange(pn, pe, nrow = 1)

  return(pm)
}


# b0 = 0.65
# b1 = 0.4
# t = 3
# pm1 = get_blen_distr(b0, b1, t)
# t = 5
# pm2 = get_blen_distr(b0, b1, t)
# t = 10
# pm3 = get_blen_distr(b0, b1, t)
# ggarrange(pm1, pm2, pm3, ncol = 1)


ftype=".png"

b0 = 0.65
b1 = 0.4
ylim = 1

ntree = 60
nbins = 10

t = 3
xlim = 3
pm2 = get_blen_distr_joint(b0, b1, t, ntree, nbins, ylim, xlim)
pm2


ntree = 4
t = 5
xlim = 5
pm2t4 = get_blen_distr_joint(b0, b1, t, ntree, nbins, ylim, xlim)
pm2t4

# may different from that in paper, due to random seed
ggarrange(pm2t4, pm2, ncol = 1)

fout = file.path(dir, paste0("sim_bl_hist_b0", b0, "_b1", b1, "_combined", ftype))
ggsave(fout, width = 10, height = 6)
