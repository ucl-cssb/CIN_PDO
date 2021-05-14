# Test the imaged cell division tree to see if there is evidence of (negative) selection

library(ggpubr)

# assuming working directory is "script"

source("./lrt/lrt_util.R")


#cutoff1 = qgamma(0.95, 1/2, 1/2)
cutoff1 = qchisq(0.95, 1)
cutoff1

# Read all trees
dir="../data/trees"
dir1 = file.path(dir, "4trees")
dir2 = file.path(dir, "60trees")

# use normal phylogenetic tree (leave: extant cell, internal node: birth event (cell division))
# timing at first 2 branches are artificial
# timing at the 3rd to 6th branches are underestimated since the imaging starts from two cells
infiles1 = list.files(dir1, pattern = "[0-9].nwk$")
# only consider old datasets
infiles1
infiles2 = list.files(dir1, pattern = ".*div.*[1-9].nwk$")
infiles2
infiles_t4 = setdiff(infiles1, infiles2)
infiles_t4


infiles_t60 = list.files(dir2, pattern = "[0-9].nwk$")
# only consider old datasets
infiles_t60

infiles = c()
for(i in 1:length(infiles_t4)){
  fname = file.path(dir1, infiles_t4[i])
  infiles = c(infiles, fname)
}
for(i in 1:length(infiles_t60)){
  fname = file.path(dir2, infiles_t60[i])
  infiles = c(infiles, fname)
}


####### Read all trees and fit Yule model #########
# all the trees are not  significantly unbalanced
get_trees_all <- function(infiles, replace_prefix = T, convert_bl = T){
  yule_all = data.frame()
  trees_all = list()
  trees_with_death = list()
  trees_all_orig = list()
  trees_all_yule = list()  # exclude tree with death events

  for(i in 1:length(infiles)){
    #i = 3
    print(infiles[i])
    fname = infiles[i]
    bname = basename(fname)

    fields = str_split(bname, "\\.")
    dset = fields[[1]][1]
    if(replace_prefix){
      title = str_replace(dset, "dataset_", "")
    }else{
      title = dset
    }

    t1 = treeio::read.newick(fname)
    trees_all_orig[[title]] = t1

    edepths = node.depth.edgelength(t1)
    tdepth = max(node.depth.edgelength(t1))
    # nodes with the same max depth, should equal to #leaves
    ned = length(which(edepths==tdepth))
    nleave = t1$Nnode + 1
    extra = ned - nleave
    if(extra != 0 && abs(extra) < 5){
      message("Not all leaves have same depth at tree ", title, " with ", extra, " events")
      print(edepths)
      print(tdepth)
      trees_with_death[[title]] = t1
    }

    if(convert_bl){
      if(startsWith(dset, "dataset")){
        t1$edge.length = t1$edge.length * 3 / 60 / 24
        ttime = tdepth * 3 / 60 / 24
      }else{
        t1$edge.length = t1$edge.length * 4.5 / 60 / 24
        ttime = tdepth * 4.5 / 60 / 24
      }
    }else{
      ttime = tdepth
    }

    # add one tree with artificial 1 to avoid 0
    if(extra == 0 || abs(extra) > 5){
      trees_all_yule[[title]] = t1
    }

    trees_all[[title]] = t1
    ntip = t1$Nnode + 1

    # not use root edge since it tends to be shorter than real values
    ryule = yule(t1)
    #print(ryule)

    yule_all = rbind(yule_all, data.frame(lambda=ryule$lambda, se=ryule$se, loglik=ryule$loglik, ntip = ntip, dataset=dset, depth = tdepth, time = ttime, id=i, ndeath = abs(extra)))
  }

  yule_all$nedge = 2 * yule_all$ntip - 2

  return(list(trees_all = trees_all, yule_all = yule_all, trees_all_orig = trees_all_orig, trees_with_death = trees_with_death, trees_all_yule = trees_all_yule))
}

res = get_trees_all(infiles)

trees_all = res$trees_all
trees_all_orig = res$trees_all_orig

yule_all = res$yule_all
yule_all %>% arrange(lambda)
yule_all %>% filter(ntip > 2) -> yule_l2
trees_with_death = res$trees_with_death
trees_all_yule = res$trees_all_yule
length(trees_with_death)
length(trees_all_yule)
names(trees_with_death)

yule_all %>% dplyr::group_by(ndeath) %>% tally()
# Not all leaves have same depth at tree early_P6_OE_275 with -7 events (two branched add 1)
# Not all leaves have same depth at tree early_P6_OE_93 with -3 events

# ggplot(yule_l2, aes(lambda)) + geom_histogram(aes(y=..density..), colour="black", fill="white", bins=30)+ geom_density(alpha=.5, fill="lightblue") + theme_pubr() + xlab("birth rate")
# fout = file.path(dir, "sum_yule_lambda.pdf")
# ggsave(fout, width = 7, height = 5)


####### Plot all trees #########
plot_tree <- function(t1, title = "", use_nid = F){
  p = ggtree(t1) + geom_tiplab()

  if(use_nid){
    p <- p + geom_text2(aes(subset = !isTip, label = node), hjust = -0.3)
  }else{
    p <- p + geom_text2(aes(subset = !isTip, label = label), hjust = -0.3)
  }

  edge = data.frame(t1$edge, edge_num = 1:nrow(t1$edge), edge_len = t1$edge.length)
  colnames(edge) = c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.3) + ggtitle(title)

  return(p)
}

get_all_tplots <- function(infiles){
  tlist = list()
  # all the trees are not  significantly unbalanced
  for(i in 1:length(infiles)){
    #i = 19
    message("tree ", i)
    print(infiles[i])
    fname = infiles[i]
    bname = basename(fname)

    fields = str_split(bname, "\\.")
    dset = fields[[1]][1]
    if(startsWith(dset, "dataset_")){
      title = str_replace(dset, "dataset_", "")
    }else{
      title = dset
    }

    t1 = treeio::read.newick(fname)
    #t1$edge.length = t1$edge.length * 3 / 60 / 24

    # plot the tree
    p = plot_tree(t1, title = title)
    print(p)
    tlist[[i]]= p
  }
  return(tlist)
}



tlist = get_all_tplots(infiles)

ntree = length(infiles)
ntree
size=14
ngroup = ceiling(ntree / size)
ngroup

# Plot trees by group and then combine all pdfs together
for(i in 1:ngroup){
  s = (i-1)*size + 1
  e = s + size - 1
  if(e > ntree){
    e = ntree
  }
  print(s)
  print(e)
  ggarrange(plotlist = tlist[s:e], nrow = size/2, ncol = 2)
  fout = file.path(dir, paste0("tree_all_phylo_lbl_new", i, ".pdf"))
  ggsave(fout, width = 30, height = 49)
}




############ compute likelihood of real trees ##################

# test individual tree
test_t1 <- function(t1, tname, fun_est, fun_nlnl, use_mle, estimate_lambda0 = F){
  # Likelihood under HO

  if(estimate_lambda0){
    lambda0 = fun_est(t1)
    res_h0 = get_lnl_h0(t1, lambda0)
  }else{
    res_h0 = get_lnl_h0_yule(t1)
    lambda0 = res_h0$lambda
  }

  lnl = res_h0$loglik
  # print(lnl)
  # print(lambda0)

  # when lambda == lambda0, lnl_decom = lnl
  res = fun_nlnl(lambda0, t1, lambda0, use_mle)
  if(use_mle){
    lambda1 = res$lambda1
    lnl_decom = -res$nlnl_decom
  }else{
    lnl_decom = -res
    lambda1 = lambda0
  }

  print(lnl_decom)
  if(!use_mle && !near(lnl_decom, lnl)){
    message("Error in computing likelihood!")
  }

  if(use_mle){
    lnl_h1_mle = lnl_decom
  }else{
    # when lambda != lambda0, lnl is a function of lambda1, use nlm to get mle
    mle_h1 = nlm(fun_nlnl, lambda0, t1 = t1, lambda0 = lambda0, use_mle = use_mle, print.level = 0, iterlim = 100)
    #print(mle_h1)
    lambda1 = mle_h1$estimate
    message("NLM of birth rate under H1: ", mle_h1$estimate)
    lnl_h1_mle = -mle_h1$minimum
  }

  # compute LR
  lnl_h0_mle = lnl
  lr = -2 * (lnl_h0_mle - lnl_h1_mle)
  print(lr)

  # test significance with chi-square distribution
  print(lr > cutoff1)
  #1 - pgamma(lr, 1/2, 1/2)
  preject = 1 - pchisq(lr, 1)
  print(preject)

  if(preject < 0.05){
    message("H0 rejected in ", tname)
  }else{
    message("H0 not rejected in ", tname)
  }

  df = data.frame(dataset=tname, lambda0=res_h0$lambda, loglik0=res_h0$loglik, lambda1 = lambda1, loglik1 = lnl_h1_mle, test_stat = lr, pvalue = preject, use_mle = use_mle)

  return(df)
}


# Compute the likelihood of all trees under H0 (correct by comparing with yule())
get_lnl_h0_all <- function(trees_all, lambda, use_root = T){
  loglik = 0

  for(i in 1:length(trees_all)){
    #i = 5
    #lambda = yule_all[i, ]$lambda
    phy = trees_all[[i]]
    tname = names(trees_all)[i]
    X <- sum(phy$edge.length)
    nb.node <- phy$Nnode - 1

    if(!use_root && !startsWith(tname, "ear")){
      # exclude first two branches and birth events, only for four old datasets
      root = phy$Nnode + 2
      edges = as.data.frame(phy$edge)
      edges$id = seq(1:nrow(edges))
      edges %>% filter(V1 == root) -> root_edges
      root_edges %>% select(id) %>% unlist() -> ids
      edges_root = phy$edge.length[ids]
      #print(edges_root)
      tlen = sum(edges_root)
      X = X - tlen
      nb.node = nb.node - 2
    }

    lnl = -lambda * X + lfactorial(phy$Nnode) + nb.node * log(lambda)
    #print(lnl)
    loglik = loglik + lnl
    # loglik
    # yule_all[i,]
  }

  return (list(loglik = loglik, lambda = lambda))
}


# Test computation of lambda0 and lambda1 based on first four trees (excluding root branches)
test_lambda_on_four_tree <- function(trees_all_orig){
  # estimation of lambda0, based on normal branches before any error
  blens1 = c(672, 639)
  blens2 = c(576, 728)
  blens3 = c(501, 488, 507, 500, 580, 684, 456, 336, 462, 491, 383, 361, 395, 412, 367, 512, 216, 216, 243, 243, 226, 226, 293, 293, 148, 148)
  blens4 = c(665, 638, 825, 1274, 539, 853, 1196, 1597)
  blens = c(blens1, blens2, blens3, blens4) # used for manual check the normal branches considered
  length(blens)
  sblen = sum(blens) * 3 / 60 / 24
  sblen = sum(blens)
  nenodes = 0 + 0 + 10 + 2
  nbl = (length(blens) - nenodes)
  est0 =  nbl / sblen
  print(sblen)
  print(nbl)
  print(est0)

  # estimation of lambda1, based on branches after any error
  # when using lambda0 as lambda1, the loglik is the same
  # find all branches in the four trees
  blen_all = c(trees_all_orig[[1]]$edge.length, trees_all_orig[[2]]$edge.length, trees_all_orig[[3]]$edge.length, trees_all_orig[[4]]$edge.length)
  blens_root = c(287, 611, 266, 323, 197, 159, 250, 280)
  # cannot use setdiff to find error branches as some normal and error branches may have same values
  blens_err0 = sum(blen_all) - sum(blens_root) - sum(blens)
  sblen_err = blens_err0 * 3 / 60 / 24
  sprintf("%.5f", sblen_err)
  #sblen = sum(blens)
  nenodes_err = 0
  nedge_err = 0
  for(i in 1:4){
    nenodes_err = nenodes_err + trees_all_orig[[i]]$Nnode + 1
    nedge_err = nedge_err + length(trees_all_orig[[i]]$edge.length)
  }
  nenodes_err = nenodes_err - nenodes
  nedge_err = nedge_err - length(blens_root) - length(blens)
  nbl_err = nedge_err - nenodes_err
  est1 =  nbl_err / sblen_err
  print(sblen_err)
  print(nbl_err)
  print(est1)
}


# get tlen and nbirth of branches before error from each new tree in a batch manner
# all root branches are included when reading _full.txt for early_ trees
get_tlen_nbirth <- function(dir, trees_all, suffix = "_full.txt", use_root = T){
  tlen = 0
  nbirth = 0
  blens_norm = c()

  for(i in 1:length(trees_all)){
    #i = yule_all %>% filter(dataset=="early_P5_OE_98") %>% select(id) %>% unlist()
    #i = 54
    tname = names(trees_all)[i]
    if(suffix == "_full.txt" && !startsWith(tname, "early")){
      next
    }
   # print(tname)

    t1 = trees_all[[tname]]
    root = t1$Nnode + 2
    # read full data to determine which nodes are normal
    fname = paste0(tname, suffix)
    f = file.path(dir, fname)

    dfull = read_delim(f, col_names = T, delim=" ", col_types = list(col_character(), col_integer(), col_integer(), col_integer()))
    dfull %>% group_by(type) %>% tally() -> n_type
    #print(n_type)
    n_norm = n_type[n_type$type=="N",]
    #print(n_norm)

    if(nrow(n_norm) > 0){
      nodes_norm = dfull %>% filter(type=="N") %>% select(parent) %>% unique() %>% unlist()
      if(!use_root){
        nodes_norm = setdiff(nodes_norm, root)
      }
      # all the nodes here are manual labels
      res = get_elen_nbirth_norm(t1, nodes_norm, T)
      blens_norm = c(blens_norm, res$blens_norm)
      #print(res)
      tlen = tlen + res$tlen
      nbirth = nbirth + res$nbirth
    }else{
     # message("no normal nodes for ", tname)
    }
  }

  return(list(tlen = tlen, nbirth = nbirth, blens_norm = blens_norm))
}


# estimate normal birth rate in the 8 early trees with >=10 leaves
estimate_lambda0_early_large <- function(trees_all){
  tlen = 0
  nbirth = 0

  if("early_P5_OE_98" %in% names(trees_all)){
    t1 = trees_all[["early_P5_OE_98"]]
    plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root)
    res = get_elen_nbirth_norm(t1, nodes_norm)
    print("early_P5_OE_98")
    print(res)
    tlen = tlen + res$tlen
    nbirth = nbirth + res$nbirth
  }

  if("early_P5_OE_210" %in% names(trees_all)){
    t1 = trees_all[["early_P5_OE_210"]]
    #print(t1)
    root = t1$Nnode + 2
    nodes_norm = c(root, 13, 14)
    res = get_elen_nbirth_norm(t1, nodes_norm)
    print("early_P5_OE_210")
    print(res)
    tlen = tlen + res$tlen
    nbirth = nbirth + res$nbirth
  }

  if("early_P7_OE_268" %in% names(trees_all)){
    t1 = trees_all[["early_P7_OE_268"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root)
    res = get_elen_nbirth_norm(t1, nodes_norm)
    print("early_P7_OE_268")
    print(res)
    tlen = tlen + res$tlen
    nbirth = nbirth + res$nbirth
  }

  if("early_P7_OE_84" %in% names(trees_all)){
    tlen1 = 0
    nbirth1 = 0
    t1 = trees_all[["early_P7_OE_84"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root, 18)
    res = get_elen_nbirth_norm(t1, nodes_norm)
    tlen1 = tlen1 + res$tlen
    nbirth1 = nbirth1 + res$nbirth
    # normal subtree
    nodes_norm2 = c("17", "14")
    for(j in 1:length(nodes_norm2)){
      n1 = nodes_norm2[j]
      phy = extract.clade(t1, n1)
      #plot(subtree)
      X <- sum(phy$edge.length)
      tlen1 = tlen1 + X
      nb.node <- phy$Nnode - 1
      #print(nb.node )
      nbirth1 = nbirth1 + nb.node
    }
    print("early_P7_OE_84")
    print(tlen1)
    print(nbirth1)
    tlen = tlen + tlen1
    nbirth = nbirth + nbirth1
  }

  if("early_P7_OE_195" %in% names(trees_all)){
    tlen1 = 0
    nbirth1 = 0
    t1 = trees_all[["early_P7_OE_195"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root, 15, 20, 24, 16)
    res = get_elen_nbirth_norm(t1, nodes_norm)
    tlen1 = tlen1 + res$tlen
    nbirth1 = nbirth1 + res$nbirth
    # normal subtree
    nodes_norm2 = c("25")
    for(j in 1:length(nodes_norm2)){
      n1 = nodes_norm2[j]
      phy = extract.clade(t1, n1)
      #plot(subtree)
      X <- sum(phy$edge.length)
      tlen1 = tlen1 + X
      nb.node <- phy$Nnode - 1
      nbirth1 = nbirth1 + nb.node
    }
    print("early_P7_OE_195")
    print(tlen1)
    print(nbirth1)
    tlen = tlen + tlen1
    nbirth = nbirth + nbirth1
  }

  if("early_P6_OE_8" %in% names(trees_all)){
    tlen1 = 0
    nbirth1 = 0
    t1 = trees_all[["early_P6_OE_8"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root, 12, 14)
    res = get_elen_nbirth_norm(t1, nodes_norm)
    tlen1 = tlen1 + res$tlen
    nbirth1 = nbirth1 + res$nbirth
    # normal subtree
    nodes_norm2 = c("12", "19", "18")
    for(j in 1:length(nodes_norm2)){
      n1 = nodes_norm2[j]
      phy = extract.clade(t1, n1)
      #plot(subtree)
      X <- sum(phy$edge.length)
      tlen1 = tlen1 + X
      nb.node <- phy$Nnode - 1
      nbirth1 = nbirth1 + nb.node
    }
    print("early_P6_OE_8")
    print(tlen1)
    print(nbirth1)
    tlen = tlen + tlen1
    nbirth = nbirth + nbirth1
  }

  return(list(tlen = tlen, nbirth = nbirth))
}


# compare computation of lambda0 by manual and systematic approach
test_lambda0_early_large <- function(dir, yule_all, trees_all){
  # test computation on new large trees
  yule_all %>% filter(ntip >= 10) %>% filter(!startsWith(dataset, "dataset")) %>% select(dataset) %>% unlist -> ntrees_large
  length(trees_all[ntrees_large])
  res_large = get_tlen_nbirth(dir, trees_all[ntrees_large])
  res_large

  # res_large == res2 when computing correctly
  res2 = estimate_lambda0_early_large(trees_all)
  res2
}


# get likelihood for new trees, root branches are included
get_lnl_h1_early <- function(trees_all, dir, lambda){
  lnl_decom = 0
  subtrees = list()
  i = 1

  for(j in 1:length(trees_all)){
    #i = yule_all %>% filter(dataset=="early_P5_OE_98") %>% select(id) %>% unlist()
    #i = 54
    tname = names(trees_all)[j]
    if(!startsWith(tname, "early")){
      next
    }
    #print(tname)
    lnl_t1 = 0
    t1 = trees_all[[tname]]
    #plot_tree(t1, "", T)
    # read full data to determine which nodes are normal
    fname = paste0(tname, "_full.txt")
    f = file.path(dir, fname)

    dfull = read_delim(f, col_names = T, delim=" ", col_types = list(col_character(), col_integer(), col_integer(), col_integer()))
    dfull %>% group_by(type) %>% tally() -> n_type
    #print(n_type)
    n_norm = n_type[n_type$type=="N",]
    n_err = n_type[n_type$type=="E",]
    #print(n_norm)

    if(nrow(n_norm) > 0){
      nodes_norm = dfull %>% filter(type=="N") %>% select(parent) %>% unique() %>% unlist()
      # all the nodes here are manual labels
      lnlr = get_plnl_norm(t1, lambda, nodes_norm, T)

      if(nrow(n_err) > 0){
        # no good way to systematically extract error subtree
        # nodes_err = dfull %>% filter(type=="E") %>% select(parent) %>% unique() %>% unlist()
        # # convert nodes ID from label in the edge matrix
        # nodes_sel = c()
        # root = t1$Nnode + 2
        # for(n in nodes_err){
        #   print(n)
        #   n_id = which(t1$node.label==n)
        #   print(n_id)
        #   n_err = n_id + root - 1
        #   nodes_sel = c(nodes_sel, n_err)
        # }
        # nodes_sel
        # node.depth(t1)[nodes_sel]
        #
        if(tname == "early_P5_OE_98"){
          nodes_err = c("1", "2")
        }else if(tname == "early_P5_OE_210"){
          nodes_err = c("13", "19", "21")
        }else if(tname == "early_P6_OE_8"){
          nodes_err = c("17")
        }else if(tname == "early_P7_OE_84"){
          nodes_err = c("15")
        }else if(tname == "early_P7_OE_195"){
          nodes_err = c("17", "19", "21", "23")
        }else if(tname == "early_P7_OE_268"){
          nodes_err = c("17", "12")
        }else if(tname == "early_P5_OE_45"){
          nodes_err = c("8", "9")
        }else if(tname == "early_P5_OE_74"){
          nodes_err = c("9", "10")
        }else if(tname == "early_P5_OE_278"){
          nodes_err = c("12", "13", "15")
        }else if(tname == "early_P5_OE_294"){
          nodes_err = c("9", "10")
        }else if(tname == "early_P6_OE_111"){
          nodes_err = c("10", "11")
        }else if(tname == "early_P6_OE_128"){
          nodes_err = c("12")
        }else if(tname == "early_P6_OE_142"){
          nodes_err = c("10", "11")
        }else if(tname == "early_P6_OE_174"){
          nodes_err = c("9", "13")
        }else if(tname == "early_P6_OE_249"){
          nodes_err = c("12")
        }else if(tname == "early_P6_OE_275"){
          nodes_err = c("13", "15", "16")
        }else if(tname == "early_P6_OE_296"){
          nodes_err = c("10", "11", "12")
        }else if(tname == "early_P6_OE_324"){
          nodes_err = c("10", "11")
        }else if(tname == "early_P7_OE_33"){
          nodes_err = c("14", "15", "16")
        }else if(tname == "early_P7_OE_51"){
          nodes_err = c("9", "13")
        }else if(tname == "early_P7_OE_116"){
          nodes_err = c("11", "12")
        }else if(tname == "early_P7_OE_151"){
          nodes_err = c("9")
        }else if(tname == "early_P7_OE_171"){
          nodes_err = c("11", "16")
        }else if(tname == "early_P7_OE_217"){
          nodes_err = c("13")
        }else if(tname == "early_P7_OE_231"){
          nodes_err = c("7")
        }else if(tname == "early_P8_OE_10"){
          nodes_err = c("11")
        }else if(tname == "early_P8_OE_74"){
          nodes_err = c("10", "15")
        }else if(tname == "early_P8_OE_128"){
          nodes_err = c("10", "11")
        }else if(tname == "early_P8_OE_145"){
          nodes_err = c("10", "14", "15")
        }else{
          message("Error!")
        }

        for(n1 in nodes_err){
          subtree1 = extract.clade(t1, n1)
          subtrees[[i]] = subtree1
          i = i + 1
        }
      }
    }else{
      #message("no normal nodes for ", tname)
      lnlr = 0
      # The whole tree is a error subtree
      subtrees[[i]] = t1
      i = i + 1
    }

    lnl_t1 = lnl_t1 + lnlr + lfactorial(t1$Nnode)
    #print(lnl_t1)
    lnl_decom = lnl_decom + lnl_t1
  } # for loop

  return(list(lnl_decom = lnl_decom, subtrees = subtrees))
}


# compute likelihood for normal part and extract error subtrees for new trees with at least 10 leaves
get_lnl_h1_early_large <- function(trees_all, lambda){
  lnl_decom = 0
  subtrees = list()
  i = 1

  if("early_P5_OE_98" %in% names(trees_all)){
    lnl_t5 = 0
    t1 = trees_all[["early_P5_OE_98"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root)
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    lnl_t5 = lnl_t5 + lnlr + lfactorial(t1$Nnode)
    nodes_err = c("1", "2")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      subtrees[[i]] = subtree1
      i = i + 1
    }
    print(lnl_t5)
    lnl_decom = lnl_decom + lnl_t5
  }

  if("early_P5_OE_234" %in% names(trees_all)){
    lnl_t6 = 0
    t1 = trees_all[["early_P5_OE_234"]]
    #plot_tree(t1, "", T)
    lnl_t6 = lnl_t6 + lfactorial(t1$Nnode)
    nodes_err = c("25")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      subtrees[[i]] = subtree1
      i = i + 1
    }
    print(lnl_t6)
    lnl_decom = lnl_decom + lnl_t6
  }

  if("early_P5_OE_155" %in% names(trees_all)){
    lnl_t7 = 0
    t1 = trees_all[["early_P5_OE_155"]]
    #plot_tree(t1, "", T)
    lnl_t7 = lnl_t7 + lfactorial(t1$Nnode)
    nodes_err = c("14")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      subtrees[[i]] = subtree1
      i = i + 1
    }
    print(lnl_t7)
    lnl_decom = lnl_decom + lnl_t7
  }

  if("early_P5_OE_210" %in% names(trees_all)){
    lnl_t8 = 0
    t1 = trees_all[["early_P5_OE_210"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root, 13, 14)
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    lnl_t8 = lnl_t8 + lnlr + lfactorial(t1$Nnode)
    nodes_err = c("13", "19", "21")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      subtrees[[i]] = subtree1
      i = i + 1
    }
    print(lnl_t8)
    lnl_decom = lnl_decom + lnl_t8
  }

  if("early_P7_OE_84" %in% names(trees_all)){
    lnl_t9 = 0
    t1 = trees_all[["early_P7_OE_84"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root, 18)
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    lnl_t9 = lnl_t9 + lnlr + lfactorial(t1$Nnode)
    # normal subtree
    nodes_norm2 = c("14", "17")
    for(j in 1:length(nodes_norm2)){
      n1 = nodes_norm2[j]
      subtree1 = extract.clade(t1, n1)
      #plot(subtree)
      lnl1 = get_plnl(subtree1, lambda)
      lnl_t9 = lnl_t9 + lnl1
    }
    nodes_err = c("15")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      subtrees[[i]] = subtree1
      i = i + 1
    }
    print(lnl_t9)
    lnl_decom = lnl_decom + lnl_t9
  }

  if("early_P7_OE_195" %in% names(trees_all)){
    lnl_t10 = 0
    t1 = trees_all[["early_P7_OE_195"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root, 15, 20, 24, 16)
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    lnl_t10 = lnl_t10 + lnlr + lfactorial(t1$Nnode)
    # normal subtree
    nodes_norm2 = c("25")
    for(j in 1:length(nodes_norm2)){
      n1 = nodes_norm2[j]
      subtree1 = extract.clade(t1, n1)
      #plot(subtree)
      lnl1 = get_plnl(subtree1, lambda)
      lnl_t10 = lnl_t10 + lnl1
    }
    nodes_err = c("17", "19", "21", "23")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      subtrees[[i]] = subtree1
      i = i + 1
    }
    print(lnl_t10)
    lnl_decom = lnl_decom + lnl_t10
  }

  if("early_P7_OE_268" %in% names(trees_all)){
    lnl_t11 = 0
    t1 = trees_all[["early_P7_OE_268"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root)
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    lnl_t11 = lnl_t11 + lnlr + lfactorial(t1$Nnode)
    nodes_err = c("17", "12")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      subtrees[[i]] = subtree1
      i = i + 1
    }
    print(lnl_t11)
    lnl_decom = lnl_decom + lnl_t11
  }

  if("early_P6_OE_8" %in% names(trees_all)){
    lnl_t12 = 0
    t1 = trees_all[["early_P6_OE_8"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    nodes_norm = c(root, 12, 14)
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    lnl_t12 = lnl_t12 + lnlr + lfactorial(t1$Nnode)
    # normal subtree
    nodes_norm2 = c("12", "18", "19")
    for(j in 1:length(nodes_norm2)){
      n1 = nodes_norm2[j]
      subtree1 = extract.clade(t1, n1)
      #plot(subtree)
      lnl1 = get_plnl(subtree1, lambda)
      lnl_t12 = lnl_t12 + lnl1
    }
    nodes_err = c("17")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      subtrees[[i]] = subtree1
      i = i + 1
    }
    print(lnl_t12)
    lnl_decom = lnl_decom + lnl_t12
  }

  if("early_P6_OE_93" %in% names(trees_all)){
    lnl_t13 = 0
    t1 = trees_all[["early_P6_OE_93"]]
    #plot_tree(t1, "", T)
    root = t1$Nnode + 2
    lnl_t13 = lnl_t13 + lfactorial(t1$Nnode)
    nodes_err = c("12")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      subtrees[[i]] = subtree1
      i = i + 1
    }
    print(lnl_t13)
    lnl_decom = lnl_decom + lnl_t13
  }

  return(list(lnl_decom = lnl_decom, subtrees = subtrees))
}


# test computatation of error subtrees and partial normal likelihood
test_lnl_h1_early_large <- function(){
  # test computation on new large trees
  yule_all %>% filter(ntip >= 10) %>% filter(!startsWith(dataset, "dataset")) %>% select(dataset) %>% unlist -> ntrees_large
  length(trees_all[ntrees_large])
  res_large = get_lnl_h1_early_large(trees_all[ntrees_large], lambda0)
  res_large$lnl_decom
  length(res_large$subtrees)

  # res_large == res2 when computing correctly
  res2 = get_lnl_h1_early(trees_all[ntrees_large], dir, lambda0)
  res2$lnl_decom
  length(res2$subtrees)
}


# estimate normal birth rates based on all normal divisions before any error
# use_root only takes effect on four old trees
estimate_lambda0_all <- function(dir, trees_all, use_root = F){
  tlen = 0
  nbirth = 0
  blens_norm = c()

  if("20181127" %in% names(trees_all)){
    t1 = trees_all[["20181127"]]
    #print(t1)
    root = t1$Nnode + 2
    nodes_norm = c(26, 34, 27, 35, 38)
    if(use_root){
      nodes_norm = c(root, nodes_norm)
    }
    res = get_elen_nbirth_norm(t1, nodes_norm)
    blens_norm = c(blens_norm, res$blens_norm)
    tlen = tlen + res$tlen
    nbirth = nbirth + res$nbirth

    # normal subtree
    nodes_norm2 = c("22", "2121")
    for(j in 1:length(nodes_norm2)){
      n1 = nodes_norm2[j]
      phy = extract.clade(t1, n1)
      #plot(subtree)
      #print(phy$edge.length)
      blens_norm = c(blens_norm, phy$edge.length)
      X <- sum(phy$edge.length)
      tlen = tlen + X
      nb.node <- phy$Nnode - 1
      nbirth = nbirth + nb.node
    }
  }

  if("20190409_1" %in% names(trees_all)){
    t1 = trees_all[["20190409_1"]]
    #print(t1)
    root = t1$Nnode + 2
    nodes_norm = c(38)
    if(use_root){
      nodes_norm = c(root, nodes_norm)
    }
    res = get_elen_nbirth_norm(t1, nodes_norm)
    blens_norm = c(blens_norm, res$blens_norm)
    tlen = tlen + res$tlen
    nbirth = nbirth + res$nbirth
  }

  if("20190409_2" %in% names(trees_all)){
    t1 = trees_all[["20190409_2"]]
    #print(t1)
    root = t1$Nnode + 2
    nodes_norm = c(15)
    if(use_root){
      nodes_norm = c(root, nodes_norm)
    }
    res = get_elen_nbirth_norm(t1, nodes_norm)
    blens_norm = c(blens_norm, res$blens_norm)
    tlen = tlen + res$tlen
    nbirth = nbirth + res$nbirth
  }

  if("20200304" %in% names(trees_all)){
    t1 = trees_all[["20200304"]]
    #print(t1)
    root = t1$Nnode + 2
    nodes_norm = c(21, 22, 31, 36)
    if(use_root){
      nodes_norm = c(root, nodes_norm)
    }
    res = get_elen_nbirth_norm(t1, nodes_norm)
    blens_norm = c(blens_norm, res$blens_norm)
    tlen = tlen + res$tlen
    nbirth = nbirth + res$nbirth
  }

  # use root branches for all new trees
  res_large = get_tlen_nbirth(dir, trees_all)
  blens_norm = c(blens_norm, res_large$blens_norm)
  tlen = tlen + res_large$tlen
  nbirth = nbirth + res_large$nbirth

  #print(tlen)
  #print(nbirth)
  lambda = nbirth / tlen

  return(list(lambda = lambda, blens_norm = blens_norm, nbirth = nbirth))
}


# based on yule from ape, correct when comparing with yule()
get_lnl_h1_all <- function(dir, trees_all, lambda0, use_root, use_mle, use_error = F, est_lambda = T, lambda1 = 0){
  lnl_decom = 0
  lambda = lambda0
  blens_err = c()

  res_early = get_lnl_h1_early(trees_all, dir, lambda0)
  #res_early
  lnl_decom = lnl_decom + res_early$lnl_decom
  i = length(res_early$subtrees) + 1
  subtrees = res_early$subtrees

  if("20181127" %in% names(trees_all)){
    lnl_t1 = 0
    t1 = trees_all[["20181127"]]
    #print(t1)
    root = t1$Nnode + 2
    nodes_norm = c(26, 34, 27, 35, 38)
    if(use_root){
      nodes_norm = c(root, nodes_norm)
    }
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    lnl_t1 = lnl_t1 + lnlr + lfactorial(t1$Nnode)
    # normal subtree
    nodes_norm2 = c("22", "2121")
    for(j in 1:length(nodes_norm2)){
      n1 = nodes_norm2[j]
      subtree1 = extract.clade(t1, n1)
      #plot(subtree)
      lnl1 = get_plnl(subtree1, lambda)
      lnl_t1 = lnl_t1 + lnl1
    }
    nodes_err = c("12", "112", "111", "211", "2122")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      # lnl1 = get_plnl(subtree1, lambda)
      # lnl_t1 = lnl_t1 + lnl1
      subtrees[[i]] = subtree1
      i = i + 1
    }
    #print(lnl_t1)
    lnl_decom = lnl_decom + lnl_t1
  }

  if("20190409_1" %in% names(trees_all)){
    lnl_t2 = 0
    t1 = trees_all[["20190409_1"]]
    #print(t1)
    root = t1$Nnode + 2
    nodes_norm = c(38)
    if(use_root){
      nodes_norm = c(root, nodes_norm)
    }
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    #print(lnlr)
    lnl_t2 = lnl_t2 + lnlr + lfactorial(t1$Nnode)
    nodes_err = c("21", "22", "1")
    for(n1 in nodes_err){
      #print(n1)
      subtree1 = extract.clade(t1, n1)
      # lnl1 = get_plnl(subtree1, lambda)
      # # print(lnl1)
      # lnl_t2 = lnl_t2 + lnl1
      subtrees[[i]] = subtree1
      i = i + 1
    }
    #print(lnl_t2)
    lnl_decom = lnl_decom + lnl_t2
  }

  if("20190409_2" %in% names(trees_all)){
    lnl_t3 = 0
    t1 = trees_all[["20190409_2"]]
    #print(t1)
    root = t1$Nnode + 2
    nodes_norm = c(15)
    if(use_root){
      nodes_norm = c(root, nodes_norm)
    }
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    lnl_t3 = lnl_t3 + lnlr + lfactorial(t1$Nnode)
    nodes_err = c("11", "12", "2")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      # lnl1 = get_plnl(subtree1, lambda)
      # lnl_t3 = lnl_t3 + lnl1
      subtrees[[i]] = subtree1
      i = i + 1
    }
    #print(lnl_t3)
    lnl_decom = lnl_decom + lnl_t3
  }

  if("20200304" %in% names(trees_all)){
    lnl_t4 = 0
    t1 = trees_all[["20200304"]]
    #print(t1)
    root = t1$Nnode + 2
    nodes_norm = c(21, 22, 31, 36)
    if(use_root){
      nodes_norm = c(root, nodes_norm)
    }
    lnlr = get_plnl_norm(t1, lambda, nodes_norm)
    lnl_t4 = lnl_t4 + lnlr + lfactorial(t1$Nnode)
    nodes_err = c("12", "112", "21", "221")
    for(n1 in nodes_err){
      subtree1 = extract.clade(t1, n1)
      # lnl1 = get_plnl(subtree1, lambda)
      # lnl_t4 = lnl_t4 + lnl1
      subtrees[[i]] = subtree1
      i = i + 1
    }
    #print(lnl_t4)
    lnl_decom = lnl_decom + lnl_t4
  }

  # estimate lambda1
  nbirth = 0
  if(est_lambda){
    if(use_mle){
      tlen = 0
      for(j in 1:length(subtrees)){
        phy = subtrees[[j]]
        blens_err = c(blens_err, phy$edge.length)
        X <- sum(phy$edge.length)
        tlen = tlen + X
        nb.node <- phy$Nnode - 1
        nbirth = nbirth + nb.node
      }
      lambda1 = nbirth / tlen
    }else{
      lambda1 = lambda0
    }
  }

  message("birth rate under H1: ", lambda1)

  # # adjust lambda1 based on error model
  # if(use_error){
  #   err = 0.5
  #   # lambda1 = 2 * lambda1 - lambda0
  #   lambda1 = err * lambda1  + (1 - err) * lambda0
  # }
  #

  for(j in 1:length(subtrees)){
    subtree1 = subtrees[[j]]
    #plot_tree(subtree1)
    lnl1 = get_plnl(subtree1, lambda1)
    #print(lnl1)
    lnl_decom = lnl_decom + lnl1
  }
  #print(lnl_decom)

  # check branch length distributions


  return (list(loglik = lnl_decom, lambda = lambda1, blens_err = blens_err, nbirth = nbirth))
}


# treat all datasets as a whole
# joint log likelihood, sum of likelihood for all trees, assuming same lambda0 and lambda1
# set use_mle = F, to check accuracy of computations
test_tree_all <- function(dir, trees_all, use_root, use_error = F, use_mle = T, dset = "joint", est_lambda = T, lambda0 = 0, lambda1 = 0){
  tall = trees_all
  df_brate = data.frame()  # for birth rate estimated from E nor N branches
  # lambda0 is estimated on the whole tree under H0
  res_h0 = get_lnl_h0_all_yule(tall, use_root)
  #res_h0 = get_lnl_h0_all(tall, lambda0, use_root)
  # res_h0
  # estimate_lambda(tall, use_root)

  # lambda0 is estimated on the normal part under H1
  if(est_lambda){
    res_lambda0 = estimate_lambda0_all(dir, tall, use_root)
    lambda0 = res_lambda0$lambda
    blens_norm = res_lambda0$blens_norm
    # hist(blens_norm)
    # sum(blens_norm)
    # length(blens_norm)
    # res_lambda0$nbirth
    brate = res_lambda0$nbirth / sum(blens_norm)
    df_norm = data.frame(sum_blens = sum(blens_norm), nbirth = res_lambda0$nbirth, birth_rate = brate, nblens = length(blens_norm), type = "N")
  }
  res_h1 = get_lnl_h1_all(dir, tall, lambda0, use_root, use_mle, use_error, est_lambda, lambda1)
  res_h1$lambda0 = lambda0
  blens_err = res_h1$blens_err
  # hist(blens_err)
  # sum(blens_err)
  # length(blens_err)
  # res_h1$nbirth
  brate = res_h1$nbirth / sum(blens_err)
  df_err = data.frame(sum_blens = sum(blens_err), nbirth = res_h1$nbirth, birth_rate = brate, nblens = length(blens_err), type = "E")
  df_brate = rbind(df_norm, df_err)

  if(use_mle){
    lr = -2 * (res_h0$loglik - res_h1$loglik)
    #print(lr)

    # test significance with chi-square distribution
    #print(lr > cutoff1)
    preject = 1 - pchisq(lr, 1)
    #print(preject)

    if(preject < 0.05){
      message("H0 rejected in joint trees")
    }else{
      message("H0 not rejected in joint trees")
    }

    df = data.frame(dataset=dset, lambda0=res_h0$lambda, loglik0=res_h0$loglik, lambda0_1=res_h1$lambda0, lambda1 = res_h1$lambda, loglik1 = res_h1$loglik, test_stat = lr, pvalue = preject)
    if(use_error){
      df$lambda_err = 2 * res_h1$lambda - res_h0$lambda
    }else{
      df$lambda_err = res_h1$lambda
    }
    # compute selection coefficient under H1
    df$s = df$lambda_err / df$lambda0_1  - 1

    if(use_root){
      df$remark = "include root branches"
      df_brate$remark = "include root branches"
    }else{
      df$remark = "exclude root branches"
      df_brate$remark = "exclude root branches"
    }
  }else{
    #near(lnl_h0_mle, lnl_h1_mle)
    df = data.frame(dataset="joint", lambda0=res_h0$lambda, loglik0=res_h0$loglik, lambda0_1=res_h1$lambda0, lambda1 = res_h1$lambda, loglik1 = res_h1$loglik, remark = "correct when loglik0=loglik1")
  }

  return(list(df = df, df_brate = df_brate))
}



############## joint test of rate change on image trees ############
get_res_joint <- function(dir, tall, use_mle, dset, use_root = F, est_lambda = T, lambda0 = 0, lambda1 = 0){
  # not use root branches
  res_joint = test_tree_all(dir, tall, F, F, use_mle, dset, est_lambda, lambda0, lambda1)
  df_joint = res_joint$df
  df_joint$lambda_err=NULL
  df_brate = res_joint$df_brate

  if(use_root){
    res1 = test_tree_all(dir, tall, T, F, use_mle, dset, est_lambda, lambda0, lambda1)
    df1 = res1$df
    df1$lambda_err=NULL
    df_joint = rbind(df_joint, df1)
    df_brate = rbind(df_brate, res1$df_brate)
    df_joint$lambda_err=NULL
  }

  return(list(df_joint = df_joint, df_brate = df_brate))
}


get_res_set <- function(tall, use_mle, use_root){
  tnames = names(tall)
  res_set = data.frame()
  df_brate = data.frame()
  sets=paste0("P", seq(5,8))
  #sets
  for(s in sets){
    #tnames_p5 = yule_all %>% filter(grepl(s, dataset)) %>% select(dataset) %>% unlist()
    tnames_p5 = tnames[grep(s, tnames)]
    ntree = length(tnames_p5)
    #print(length(tnames_p5))
    trees_p5 = tall[tnames_p5]
    #print(length(trees_p5))

    res_p5 = get_res_joint(trees_p5, use_mle, paste0(s, " (", ntree, " trees)"), use_root)
    #print(res_p5)
    res_set = rbind(res_set, res_p5$df_joint)
    df_brate = rbind(df_brate, res_p5$df_brate)
  }

  return(list(res_set = res_set, df_brate = df_brate))
}


# use root branches
use_mle = T
ntrees = length(trees_all)
ntrees
use_root = T
tall = trees_all

res_lrt = data.frame()
# only need to exclude root branches for the 4 old trees
res_t4 = get_res_joint(dir1, trees_all[1:4], use_mle, "4 trees", use_root)
res_t60 = get_res_joint(dir2, trees_all[5:64], use_mle, "60 trees", use_root)
res_lrt = rbind(res_lrt, res_t4$df_joint)
res_lrt = rbind(res_lrt, res_t60$df_joint)
res_lrt

res_lrt %>% filter((dataset == "4 trees" & remark == "exclude root branches") | (dataset == "60 trees" & remark == "include root branches")) -> res_lrt_sel
res_lrt_sel %>% mutate(across(where(is.numeric), round, 3)) -> res_lrt_self

# write results to a file ((Suppl Table 2)
fout = file.path(dir, "sum_lrt_new.tsv")
write_tsv(res_lrt_sel, fout)


# on all pure pure trees (for last paragraph in Suppl doc)
names(trees_all_yule)
res_yule = get_res_joint(dir2, trees_all_yule, use_mle, paste0(length(trees_all_yule), " trees"), use_root)
res_yule




######################## check branch length distribution in real data (Suppl Fig. 5) #########################
# For new datasets (60 trees)
get_new_blen <- function(dir, files, use_root = F){
  dblen = data.frame()
  for(i in 1:length(files)){
    #i = 2
    fname1 = file.path(dir, files[i])
    print(fname1)

    fields = str_split(files[i], "\\.")
    dset = fields[[1]][1]
    dname = dset

    dblen1 = read_delim(fname1, col_names = T, delim=" ", col_types = list(col_character(), col_integer(), col_integer(), col_integer()))
    dblen1$dataset = dname

    # Find terminal edges (end nodes does not appear in parent)
    ends = dblen1$child
    starts = unique(dblen1$parent)
    tedges_e = setdiff(ends, starts)
    dblen1 %>% mutate(lineage = if_else(child %in% tedges_e, "T", "M")) -> dblen1

    if(use_root){
      root = max(tedges_e) + 1
      dblen1 %>% mutate(lineage = if_else(parent == root, "F", lineage)) -> dblen1
    }

    # print(dblen1)
    dblen = rbind(dblen, dblen1)
  }
  return(dblen)
}



get_blen_hist_real <- function(dblens, nbins = 30, use_freq = F, ylim = 1.0, xlim = 5, add_lbl = F,  add_curve = F, add_fit = F, use_density = F, lambdas = c(), nrow = 1, color = "blue"){
  #dblens = dblen0_full_nr
  bsize = xlim / nbins
  br = seq(0, xlim, bsize)
  dblens %>% filter(type == "N") -> dblen_N
  dblens %>% filter(type == "E") -> dblen_E
  # hist(dblen_E$time)
  print(sprintf("%.5f", sum(dblen_N$time)))
  print(sprintf("%.5f", sum(dblen_E$time)))
  prec = 3

  ntree = length(unique(dblens$dataset))
  nblen = nrow(dblen_N)
  main = paste0(ntree, " trees, ", nblen, " branches")
  #main2 = paste0(ntree, " trees, ", nblen, " branches", "\nestimated birth rate ", sprintf("%.2f", lambdas[1]))
  main2 = paste0(ntree, " trees, ", nblen, " branches", "\nestimated cell cycle duration ", sprintf("%.1f", 1 / lambdas[1]), " (day)")
  if(use_freq){
    pnr = ggplot(dblen_N, aes(time, y = stat(count) / sum(count))) + geom_histogram(breaks = br, color="black", fill="lightblue") + xlab("branch length before mis-segregation") + scale_x_continuous(limits = c(0, xlim)) + theme_pubr() + ggtitle(main) + ylab("frequency") + scale_y_continuous(limits = c(0, ylim))
    # show freq values for each bar
    if(add_lbl){
      pnr = pnr + geom_text(aes(label = round(stat(count) / sum(count), prec)), stat = 'bin', vjust = -0.5, breaks = br)
    }
  }else if(use_density){
    pnr = ggplot(dblen_N, aes(time)) + geom_histogram(aes(y = stat(count) / (sum(count) * bsize)), breaks = br, color="black", fill="lightblue") + xlab("branch length before mis-segregation") + scale_x_continuous(limits = c(0, xlim)) + theme_pubr() + ggtitle(main) + ylab("density") + scale_y_continuous(limits = c(0, ylim))
    # pnr = pnr + geom_text(aes(label = round(stat(density), prec)), stat = 'bin', vjust = -0.5, breaks = br)
  }else{
    pnr = ggplot(dblen_N, aes(time)) + geom_histogram(breaks = br, color="black", fill="lightblue") + xlab("branch length before mis-segregation") + scale_x_continuous(limits = c(0, xlim)) + theme_pubr() + ggtitle(main)
  }
  if(add_curve){
    message("add exponential distribution with rate ", 2 * lambdas[1])
    pnr = pnr + stat_function(fun=dexp, geom = "line", size=1, col=color, linetype = "dashed", args = (rate= 2 * lambdas[1])) + ggtitle(main2)
  }
  if(add_fit){
    df_exp = get_exp_per_bin(br, 2 * lambdas[1])
    # ggplot(df_exp, aes(x = bin, y=val)) + geom_point() + geom_line(color = "red", linetype="dashed") + theme_pubr()
    pnr = pnr + geom_point(data = df_exp, mapping = aes(x = bin, y=val)) + geom_line(data = df_exp, mapping = aes(x = bin, y=val), color = color, linetype="dashed") + ggtitle(main2)
    #+ geom_line(color = "red", linetype="dashed")
  }


  nblen = nrow(dblen_E)
  main = paste0(ntree, " trees, ", nblen, " branches")
  # main2 = paste0(ntree, " trees, ", nblen, " branches", "\nestimated birth rate ", sprintf("%.2f", lambdas[2]))
  main2 = paste0(ntree, " trees, ", nblen, " branches", "\nestimated cell cycle duration ", sprintf("%.1f", 1 / lambdas[2]), " (day)")
  if(use_freq){
    per = ggplot(dblen_E, aes(time, y = stat(count) / sum(count))) + geom_histogram(breaks = br, color="black", fill="lightblue") + xlab("branch length after mis-segregation")  + scale_x_continuous(limits = c(0, xlim)) + theme_pubr() + ggtitle(main) + ylab("frequency") + scale_y_continuous(limits = c(0, ylim))
    if(add_lbl){
      per = per + geom_text(aes(label = round(stat(count) / sum(count), prec)), stat = 'bin', vjust = -0.5, breaks = br)
    }
  }else if(use_density){
    per = ggplot(dblen_E, aes(time)) + geom_histogram(aes(y = stat(count) / (sum(count) * bsize)), breaks = br, color="black", fill="lightblue") + xlab("branch length before mis-segregation") + scale_x_continuous(limits = c(0, xlim)) + theme_pubr() + ggtitle(main) + ylab("density") + scale_y_continuous(limits = c(0, ylim))
  }else{
    per = ggplot(dblen_E, aes(time)) + geom_histogram(breaks = br, color="black", fill="lightblue") + xlab("branch length after mis-segregation")  + scale_x_continuous(limits = c(0, xlim)) + theme_pubr() + ggtitle(main)
  }
  if(add_curve){
    message("add exponential distribution with rate ", 2 * lambdas[2])
    per = per +stat_function(fun=dexp,geom = "line",size=1, col=color, linetype = "dashed", args = (rate= 2 * lambdas[2])) + ggtitle(main2)
  }
  if(add_fit){
    df_exp = get_exp_per_bin(br, 2 * lambdas[2])
    # ggplot(df_exp, aes(x = bin, y=val)) + geom_point() + geom_line(color = "red", linetype="dashed") + theme_pubr()
    per = per + geom_point(data = df_exp, mapping = aes(x = bin, y=val)) + geom_line(data = df_exp, mapping = aes(x = bin, y=val), color = color, linetype="dashed") + ggtitle(main2)
    #+ geom_line(color = "red", linetype="dashed")
  }

  pm = ggarrange(pnr, per, nrow = nrow)

  return(pm)
}


###### read 4 trees
pattern = ".*_full_elist.txt"
files = list.files(dir1, pattern = pattern)
files
length(files)


dblen0_full = get_new_blen(dir1, files, T)
dblen0_full$time = dblen0_full$length * 3 / 60 / 24

dblen0_full_nr = dblen0_full %>% filter(lineage!="F")

###### read 60 trees
pattern = "early.*_full.txt"
files = list.files(dir2, pattern = pattern)
files
length(files)

dblen = get_new_blen(dir2, files)
dblen$time = dblen$length * 4.5 / 60 / 24
dblen$dataset = str_replace(dblen$dataset, "_full", "")

dblen$ntip = dblen
dblen %>% group_by(dataset) %>% tally() -> dblen_n
dblen_n$ntip = dblen_n$n / 2 + 1
dblen_n
length(unique(dblen$dataset))
dblen %>% group_by(dataset, lineage) %>% tally() %>% filter(lineage!="T") %>% arrange(n)



nbins = 10
use_freq = T
ylim = 1.0
add_lbl = F
add_curve = F
add_fit = T
use_density = F
nrow = 1
lambdas60 = c(0.6846149, 0.3986991)
xlim = 3
pt60o = get_blen_hist_real(dblen, nbins, use_freq, ylim, xlim, add_lbl, add_curve, add_fit, use_density, lambdas60, nrow)
pt60o
lambdas4 = c(0.6274825, 0.3559951)
xlim = 5
pt4_nro = get_blen_hist_real(dblen0_full_nr, nbins, use_freq, ylim, xlim, add_lbl, add_curve, add_fit, use_density, lambdas4, nrow)
pt4_nro
ggarrange(pt4_nro, pt60o, nrow = 2)

ftype = ".png"
fout = file.path(dir, paste0("sum_t4nr_t60_blen_hist_final_nbin", nbins, ftype))
ggsave(fout, width = 9, height = 6)
