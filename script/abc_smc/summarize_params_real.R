# This script is used to plot the estimates of parameters for real data

library(tidyverse)
library(ggpubr)
library(ggsci)

# assuming working directory is "script"

bdir = "../data/abcsmc/"  # base directory of the results
dir = file.path(bdir, "neutral")
dir_sel =  file.path(bdir, "selection")
ftype=".png"

suffix = ""
weight = "2.0"
epsilon = 0.2
suffix = paste0("_logs3_weight", weight)
suffix = paste0(suffix, "_logs4_mut0")


################## Functions to read ABC results #######################
# Read ABC SMC results based on inference type
get_res <- function(dir, dataset, ncell, epsilon, type, brate = NA, suffix=""){
  if(is.na(brate)){
    fname <- file.path(dir, paste0("smc-", type, "-", dataset, "_n", ncell, "_epsilon", epsilon, suffix))
  }else{
    fname <- file.path(dir, paste0("smc-", type, "-", dataset, "_n", ncell, "_epsilon", epsilon, "_brate", brate, suffix))
  }

  res = read.table(fname, header = T)
  if(type == "joint"){
    names(res) = c("chr_rate", "arm_rate", "division_rate", "dist", "weight")
  }else{
    names(res) = c("chr_rate", "arm_rate", "dist", "weight")
    res$division_rate = -1
  }
  res$dataset = gsub("dataset_", "", dataset)
  res$type = type
  return(res)
}


# one more parameter under selection model
get_res_withsel <- function(dir, dataset, ncell, epsilon, type, brate = NA, suffix=""){
  fname <- file.path(dir, paste0("smc-", type, "-", dataset, "_n", ncell, "_epsilon", epsilon, suffix))
  res = read.table(fname, header = T)

  names(res) = c("chr_rate", "arm_rate", "division_rate", "fitness", "dist", "weight")

  res$dataset = gsub("dataset_", "", dataset)
  res$type = type

  return(res)
}


# epsilon is specified for joint inference, which is fixed to be 0.015 when doing inference with CNA only
get_res_diploid <- function(dir, suffix="", epsilon=0.1, incl_cna_only = T, withsel = F){
  resall = data.frame()
  rateall = data.frame()

  # Get joint inference
  type="joint"
  if(withsel){   # under selection model
    type = "joint-withsel"
  }
  brate = NA

  dataset = "dataset_20190409_1"
  ncell = 25
  res = get_res(dir, dataset, ncell, epsilon, type, brate, suffix)
  resall = rbind(resall, res)
  #rates = data.frame(division_rate=median(res$division_rate), chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rates = data.frame(division_rate=weighted.mean(res$division_rate, res$weight), chr_rate=weighted.mean(res$chr_rate, res$weight), arm_rate=weighted.mean(res$arm_rate, res$weight), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)
  #names(rateall) = c("division_rate", "chr_rate", "arm_rate", "dataset")

  # dataset = "dataset_20190409_2"
  # ncell = 13
  # res = get_res(dir, dataset, ncell, epsilon, type, brate, suffix)
  # resall = rbind(resall, res)
  # #rates = data.frame(division_rate=median(res$division_rate), chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  # rates = data.frame(division_rate=weighted.mean(res$division_rate, res$weight), chr_rate=weighted.mean(res$chr_rate, res$weight), arm_rate=weighted.mean(res$arm_rate, res$weight), dataset=dataset, ncell = ncell)
  # rateall = rbind(rateall, rates)

  dataset = "dataset_20181127"
  ncell = 24
  res = get_res(dir, dataset, ncell, epsilon, type, brate, suffix)
  resall = rbind(resall, res)
  #rates = data.frame(division_rate=median(res$division_rate), chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rates = data.frame(division_rate=weighted.mean(res$division_rate, res$weight), chr_rate=weighted.mean(res$chr_rate, res$weight), arm_rate=weighted.mean(res$arm_rate, res$weight), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)

  dataset = "dataset_20200304"
  ncell = 19
  res = get_res(dir, dataset, ncell, epsilon, type, brate, suffix)
  resall = rbind(resall, res)
  #rates = data.frame(division_rate=median(res$division_rate), chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rates = data.frame(division_rate=weighted.mean(res$division_rate, res$weight), chr_rate=weighted.mean(res$chr_rate, res$weight), arm_rate=weighted.mean(res$arm_rate, res$weight), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)


  # Get CNA inference
  if(incl_cna_only){
    brate = 0.5
    #epsilon=0.01
    epsilon=0.015
    type="cna"

    dataset = "dataset_20190331"
    ncell = 67
    nseq = 62
    res = get_res(dir, dataset, ncell, epsilon, type, brate)
    resall = rbind(resall, res)
    #rates = data.frame(division_rate=brate, chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
    rates = data.frame(division_rate=brate, chr_rate=weighted.mean(res$chr_rate, res$weight), arm_rate=weighted.mean(res$arm_rate, res$weight), dataset=dataset, ncell = ncell)
    rateall = rbind(rateall, rates)

    dataset = "dataset_20190310"
    ncell = 72
    nseq = 64
    res = get_res(dir, dataset, ncell, epsilon, type, brate)
    resall = rbind(resall, res)
    #rates = data.frame(division_rate=brate, chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
    rates = data.frame(division_rate=brate, chr_rate=weighted.mean(res$chr_rate, res$weight), arm_rate=weighted.mean(res$arm_rate, res$weight), dataset=dataset, ncell = ncell)
    rateall = rbind(rateall, rates)

    #epsilon=0.015
    dataset = "dataset_20200204"
    ncell = 37
    nseq = 33
    res = get_res(dir, dataset, ncell, epsilon, type, brate)
    resall = rbind(resall, res)
    #rates = data.frame(division_rate=brate, chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
    rates = data.frame(division_rate=brate, chr_rate=weighted.mean(res$chr_rate, res$weight), arm_rate=weighted.mean(res$arm_rate, res$weight), dataset=dataset, ncell = ncell)
    rateall = rbind(rateall, rates)
  }

  return(list(resall = resall, rateall = rateall))
}


# read results on inference with CNA data only
get_res_all_cna <- function(dir, suffix=""){
  resall = data.frame()
  rateall = data.frame()

  # Get CNA inference
  brate = 0.5
  #epsilon=0.01
  epsilon=0.015
  type="cna"

  dataset = "dataset_20190331"
  ncell = 67
  nseq = 62
  ncna = 5
  res = get_res(dir, dataset, ncell, epsilon, type, brate)
  resall = rbind(resall, res)
  rates = data.frame(division_rate=brate, chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)

  dataset = "dataset_20190310"
  ncell = 72
  nseq = 64
  ncna = 4
  res = get_res(dir, dataset, ncell, epsilon, type, brate)
  resall = rbind(resall, res)
  rates = data.frame(division_rate=brate, chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)

  #epsilon=0.015
  dataset = "dataset_20200204"
  ncell = 37
  nseq = 33
  ncna = 6
  res = get_res(dir, dataset, ncell, epsilon, type, brate)
  resall = rbind(resall, res)
  rates = data.frame(division_rate=brate, chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)

  return(list(resall = resall, rateall = rateall))

}


# Get results under selection
get_res_all_withsel <- function(dir, suffix="", epsilon=0.1){
  resall = data.frame()
  rateall = data.frame()

  # Get joint inference
  type="joint-withsel"
  brate = NA

  dataset = "dataset_20190409_1"
  ncell = 25
  res = get_res_withsel(dir, dataset, ncell, epsilon, type, brate, suffix)
  resall = rbind(resall, res)
  rates = data.frame(division_rate=median(res$division_rate), chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)
  #names(rateall) = c("division_rate", "chr_rate", "arm_rate", "dataset")

  dataset = "dataset_20190409_2"
  ncell = 13
  epsilon1 = 0.1
  suffix1 = paste0("_logs3_weight", weight)
  res = get_res_withsel(dir, dataset, ncell, epsilon1, type, brate, suffix1)
  resall = rbind(resall, res)
  rates = data.frame(division_rate=median(res$division_rate), chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)

  dataset = "dataset_20181127"
  ncell = 24
  res = get_res_withsel(dir, dataset, ncell, epsilon, type, brate, suffix)
  resall = rbind(resall, res)
  rates = data.frame(division_rate=median(res$division_rate), chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)

  dataset = "dataset_20200304"
  ncell = 19
  res = get_res_withsel(dir, dataset, ncell, epsilon, type, brate, suffix)
  resall = rbind(resall, res)
  rates = data.frame(division_rate=median(res$division_rate), chr_rate=median(res$chr_rate), arm_rate=median(res$arm_rate), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)

  return(list(resall = resall, rateall = rateall))
}


get_res_all_neutral <- function(dir, suffix, epsilon){
  res_dip = get_res_diploid(dir, suffix, epsilon)
  resall = res_dip$resall
  rateall = res_dip$rateall

  # exclude results for dataset_20190409_2
  resall %>% select(dataset) %>% unique()
  resall = resall %>% filter(dataset != "20190409_2")
  resall %>% select(dataset) %>% unique()

  rateall %>% select(dataset) %>% unique()
  rateall = rateall %>% filter(dataset != "dataset_20190409_2")
  rateall %>% select(dataset) %>% unique()

  # read mixed input, _logs4_mut0 for 3 datasets, _logs3 for one dataset
  dataset = "dataset_20190409_2"
  ncell = 13
  type = "joint"
  brate = NA
  epsilon = 0.1
  suffix = paste0("_logs3_weight", weight)
  res = get_res(dir, dataset, ncell, epsilon, type, brate, suffix)
  resall = rbind(resall, res)
  rates = data.frame(division_rate=weighted.mean(res$division_rate, res$weight), chr_rate=weighted.mean(res$chr_rate, res$weight), arm_rate=weighted.mean(res$arm_rate, res$weight), dataset=dataset, ncell = ncell)
  rateall = rbind(rateall, rates)

  resall$mut_rate = resall$chr_rate + resall$arm_rate
  rateall$dataset=gsub("dataset_", "", rateall$dataset)
  rateall$ratio = rateall$arm_rate/rateall$chr_rate
  rateall$mut_rate = rateall$chr_rate + rateall$arm_rate

  dlevels = c("20181127", "20190409_1", "20190409_2", "20200304", "20190310", "20190331", "20200204")
  resall$dataset <- factor(resall$dataset, levels=dlevels)

  return(list(resall = resall, rateall = rateall))
}


################## Functions to plot ABC results #######################

# Plot mutation rate jointly or in separate
# weights in ABC SMC should always be used
plot_rate_all <- function(rateall, resall, dir, pal, suffix, ftype = ".pdf"){
  rateall %>% arrange(desc(mut_rate)) %>% select(dataset) -> d_mut
  resall$mut_rate = resall$chr_rate + resall$arm_rate
  suffix = paste0(suffix, "_weighted")

  p = ggplot(resall, aes(x=dataset, y=mut_rate, fill=dataset, weight = weight)) +
      geom_violin(scale = "width") +
      geom_boxplot(width=0.1, fill="white") + ylab("CNA rate (per division)") + theme_pubr()  + scale_fill_manual(name = "Dataset: ", labels = c("PDTO-9 #1 ", "PDTO-9 #4 ", "PDTO-9 #5 ", "PDTO-9 #7 ", "PDTO-9 #2 ", "PDTO-9 #3 ", "PDTO-9 #6 "), values = pal) + guides(fill = guide_legend(nrow = 1)) + theme_pubr() + scale_x_discrete(limits=d_mut$dataset) + coord_flip() + xlab("")  + theme(axis.text.y = element_blank()) + scale_y_continuous(n.breaks = 10) + geom_hline(yintercept = weighted.mean(resall$mut_rate, resall$weight), color = "darkgrey", linetype = "dashed")
  fout = file.path(dir, paste0("sum_mu_boxplot", suffix, ".png"))
  ggsave(fout, p, width = 10, height = 6)
  fout = file.path(dir, paste0("sum_mu_boxplot", suffix, ".pdf"))
  ggsave(fout, p, width = 10, height = 6)

  rateall %>% arrange(desc(chr_rate)) %>% select(dataset) -> d_chr
  lab1 = "chr-level CNA rate (per division)"
  rateall %>% arrange(desc(arm_rate)) %>% select(dataset) -> d_arm
  #lab2 = "sub-chromosomal CNA rate (per division)"
  lab2 = "arm-level CNA rate (per division)"

  p1 = ggplot(resall, aes(x=dataset, y=chr_rate, fill=dataset, weight = weight)) +
      geom_violin(scale = "width")+
      geom_boxplot(width=0.1, fill="white") + ylab(lab1) + theme_pubr()  + scale_fill_manual(name = "Dataset: ", labels = c("PDTO-9 #1 ", "PDTO-9 #4 ", "PDTO-9 #5 ", "PDTO-9 #7 ", "PDTO-9 #2 ", "PDTO-9 #3 ", "PDTO-9 #6 "), values = pal) + guides(fill = guide_legend(nrow = 1)) + theme_pubr() + scale_x_discrete(limits=d_chr$dataset) + coord_flip() + xlab("")  + theme(axis.text.y = element_blank()) + geom_hline(yintercept = weighted.mean(resall$chr_rate, resall$weight), color = "darkgrey", linetype = "dashed") + scale_y_continuous(n.breaks = 6)
  p2 = ggplot(resall, aes(x=dataset, y=arm_rate, fill=dataset, weight = weight)) +
      geom_violin(scale = "width")+
      geom_boxplot(width=0.1, fill="white") + ylab(lab2) + theme_pubr() + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) + theme_pubr() + theme(axis.text.y = element_blank()) + scale_x_discrete(limits=d_arm$dataset) + coord_flip() + xlab("") + geom_hline(yintercept = weighted.mean(resall$arm_rate, resall$weight), color = "darkgrey", linetype = "dashed") + scale_y_continuous(n.breaks = 6)
  ggarrange(p1, p2,
            ncol = 2, nrow = 1, common.legend = T)

  fout = file.path(dir, paste0("sum_mu_sep_boxplot", suffix, ftype))
  ggsave(fout, width = 10, height = 6)
}


# Plot each parameter separately (chr-level CNA rate, arm-level CNA rate, cell division rate)
plot_by_group <- function(rateall, resall, dir, pal, suffix, ftype = ".pdf"){
  # group by rates
  # Change color by groups
  rateall %>% arrange(desc(chr_rate)) %>% select(dataset) -> d_chr
  lab1 = "chr-level CNA rate (per division)"
  #lab1 = "chromosomal CNA rate (per division)"
  p1 = ggplot(resall, aes(x=dataset, y=chr_rate, fill=dataset, weight = weight)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab(lab1) + theme_pubr()  + scale_fill_manual(name = "Dataset: ", labels = c("PDTO-9 #1 ", "PDTO-9 #4 ", "PDTO-9 #5 ", "PDTO-9 #7 ", "PDTO-9 #2 ", "PDTO-9 #3 ", "PDTO-9 #6 "), values = pal) + guides(fill = guide_legend(nrow = 1)) + theme_pubr() + scale_x_discrete(limits=d_chr$dataset) + coord_flip() + xlab("")  + theme(axis.text.y = element_blank()) + geom_hline(yintercept = weighted.mean(resall$chr_rate, resall$weight), color = "darkgrey", linetype = "dashed") + scale_y_continuous(n.breaks = 6, limits = c(0, 0.51))

  rateall %>% arrange(desc(arm_rate)) %>% select(dataset) -> d_arm
  #lab2 = "sub-chromosomal CNA rate (per division)"
  lab2 = "arm-level CNA rate (per division)"
  p2 = ggplot(resall, aes(x=dataset, y=arm_rate, fill=dataset, weight = weight)) +
    geom_violin(scale = "width")+
    geom_boxplot(width=0.1, fill="white") + ylab(lab2) + theme_pubr() + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) + theme_pubr() + theme(axis.text.y = element_blank()) + scale_x_discrete(limits=d_arm$dataset) + coord_flip() + xlab("") + geom_hline(yintercept = weighted.mean(resall$arm_rate, resall$weight), color = "darkgrey", linetype = "dashed") + scale_y_continuous(n.breaks = 6)
  p2

  resall %>% filter(type=="joint") -> resjoint
  # summary(resjoint$division_rate)
  # var(resjoint$division_rate)
  p3 = ggplot(resjoint, aes(x=dataset, y=division_rate, fill=dataset, weight = weight)) +
    geom_violin(scale = "width") +
    geom_boxplot(width=0.1, fill="white") + ylab("cell division rate (per day)") + theme_pubr() + theme(axis.text.x = element_blank()) + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) + theme_pubr() + theme(axis.text.y = element_blank()) + coord_flip() + xlab("") + geom_hline(yintercept = weighted.mean(resjoint$division_rate, resjoint$weight), color = "darkgrey", linetype = "dashed") + scale_y_continuous(n.breaks = 6)
  p3
  ggarrange(p1, p2, p3,
            ncol = 3, nrow = 1, common.legend = T)

  suffix = paste0(suffix, "_weighted")

  fout = file.path(dir, paste0("sum_params_boxplot", suffix, ftype))
  ggsave(fout, width = 10, height = 6)

}


# Plot joint CNA rate and cell division rate
plot_mu_brate <- function(rateall, resall, dir, pal, suffix, ftype){
  # group by rates
  # Change color by groups
  lab1 = "CNA rate (per division)"
  #lab1 = "chromosomal CNA rate (per division)"
  rateall %>% arrange(desc(mut_rate)) %>% select(dataset) -> d_mut
  resall$mut_rate = resall$chr_rate + resall$arm_rate

  resall %>% filter(type=="joint") -> resjoint
  # summary(resjoint$division_rate)
  # var(resjoint$division_rate)

  p1 = ggplot(resall, aes(x=dataset, y=mut_rate, fill=dataset, weight = weight)) +
      geom_violin(scale = "width")+
      geom_boxplot(width=0.1, fill="white") + ylab("CNA rate (per division)") + theme_pubr()  + scale_fill_manual(name = "Dataset: ", labels = c("PDTO-9 #1 ", "PDTO-9 #4 ", "PDTO-9 #5 ", "PDTO-9 #7 ", "PDTO-9 #2 ", "PDTO-9 #3 ", "PDTO-9 #6 "), values = pal) + guides(fill = guide_legend(nrow = 1)) + theme_pubr() + scale_x_discrete(limits=d_mut$dataset) + coord_flip() + xlab("") + theme(axis.text.y = element_blank()) + scale_y_continuous(n.breaks = 10) + geom_hline(yintercept = weighted.mean(resall$mut_rate, resall$weight), color = "darkgrey", linetype = "dashed")
  p3 = ggplot(resjoint, aes(x=dataset, y=division_rate, fill=dataset, weight = weight)) +
      geom_violin(scale = "width")+
      geom_boxplot(width=0.1, fill="white") + ylab("cell division rate (per day)") + theme_pubr() + theme(axis.text.x = element_blank())  + scale_fill_manual(values = pal) + guides(fill = guide_legend(nrow = 1)) + theme_pubr() + theme(axis.text.y = element_blank()) + coord_flip() + xlab("") + geom_hline(yintercept = weighted.mean(resjoint$division_rate, resjoint$weight), color = "darkgrey", linetype = "dashed") + scale_y_continuous(n.breaks = 6)


  ggarrange(p1, p3,
            ncol = 2, nrow = 1, common.legend = T)

  suffix = paste0(suffix, "_weighted")
  fout = file.path(dir, paste0("sum_params2_boxplot", suffix, ftype))
  ggsave(fout, width = 10, height = 6)
}


################## Functions to check ABC results #######################
get_epsilon <- function(dir, files, prefix, fprefix = "smc-joint-", ftype = ""){
  plist = list()
  res_epsilons_all = data.frame()
  nsim_all = data.frame()
  for(i in 1:length(files)){
    #i=1
    f = files[i]
    fname = file.path(dir, f)

    fields = str_replace(f, fprefix, "")
    pname = strsplit(fields,"_")[[1]][2]
    if(grepl(pname, "20190409")){
      pname = paste(pname, strsplit(fields,"_")[[1]][3], sep="_")
    }
    #pname = strsplit(pname,"-")[[1]][2]
    print(pname)

    fcon = file(fname, "r")
    res = readLines(con = fcon, n = 4)
    close(fcon)
    nsim = strsplit(res[2], ":")[[1]][2]
    tol = res[4]
    fields = strsplit(tol, ":")[[1]][2]
    fields = str_replace(fields, "\\[", "")
    fields = str_replace(fields, "\\]", "")
    tvals = strsplit(fields, ",")[[1]]

    res_epsilons = data.frame()
    for(j in 1:length(tvals)){
      df = data.frame(index = j, epsilon = as.numeric(tvals[j]))
      res_epsilons = rbind(res_epsilons, df)
    }
    res_epsilons$patient = pname
    res_epsilons_all = rbind(res_epsilons_all, res_epsilons)

    nsim_all = rbind(nsim_all, data.frame(patient = pname, nsim = nsim))

    if(ftype != ""){
      main = paste0(pname, ",\n", nsim, " simulations")
      mindex = max(res_epsilons$index)
      plt = ggplot(res_epsilons, aes(x = (index), y = epsilon)) + geom_point() + geom_line() + theme_pubr() + ggtitle(main) + xlab("population index") + scale_x_continuous(breaks = seq(1, max(mindex), 1))
      #print(plt)
      plist[[i]] = plt
    }
  }

  if(ftype != ""){
    pm = ggarrange(plotlist = plist, ncol = 4, nrow = 1)
    #print(pm)
    fout = file.path(dir, paste0("epsilon", prefix, ftype))
    ggsave(fout, width = 12, height = 5)
  }

  return(list(res_epsilons_all = res_epsilons_all, nsim_all = nsim_all))
}




################## parse just plot datasets with only CNAs (Suppl Fig. 2) ##################
res_cna = get_res_all_cna(dir)
resall = res_cna$resall
rateall = res_cna$rateall
rateall$dataset=gsub("dataset_", "", rateall$dataset)

# compared with manual counting of a more restrict set of reciprocal events
# Mapping to fig 1d
# #3
# dataset = "dataset_20190331"
# ncell = 67
# nseq = 62
# ncna = 5
#
# #2
# dataset = "dataset_20190310"
# ncell = 72
# nseq = 64
# ncna = 4
#
# #6
# dataset = "dataset_20200204"
# ncell = 37
# nseq = 33
# ncna = 6
# suffix=""
# refvals = c(4/71, 5/66, 6/36)
suffix="_seqcells" # (only consider cells sequenced)
refvals = c(4/63, 5/61, 6/32)

rateall$mut_rate = rateall$chr_rate + rateall$arm_rate
rateall %>% arrange(desc(mut_rate)) %>% select(dataset) -> d_mut
resall$mut_rate = resall$chr_rate + resall$arm_rate
pal = c('#ff7f00','#ffff33','#a65628')
dlevels = c("20190310", "20190331", "20200204")
resall$dataset <- factor(resall$dataset, levels=dlevels)

ggplot(resall, aes(x=dataset, y=mut_rate, fill=dataset, weight = weight)) +
  geom_violin(scale = "width")+
  geom_boxplot(width=0.1, fill="white") + ylab("CNA rate (per division)") + theme_pubr()  + scale_fill_manual(name = "Dataset: ", labels = c("PDTO-9 #2 ", "PDTO-9 #3 ", "PDTO-9 #6 "), values = pal) + guides(fill = guide_legend(nrow = 1)) + theme_pubr() + scale_x_discrete(limits=d_mut$dataset) + coord_flip() + xlab("") + theme(axis.text.y = element_blank()) + scale_y_continuous(n.breaks = 10) + geom_hline(yintercept = refvals[1], color = pal[1], linetype = "dashed") + geom_hline(yintercept = refvals[2], color = pal[2], linetype = "dashed") + geom_hline(yintercept = refvals[3], color = pal[3], linetype = "dashed")

fout = file.path(dir, paste0("sum_cna_boxplot", suffix, "_weighted", ftype))
ggsave(fout, width = 10, height = 6)



################ read results from all diploid datasets  ##################

res = get_res_all_neutral(dir, suffix, epsilon)

resall = res$resall
rateall = res$rateall  # weighted mean rates for each dataset


summary(resall)
unique(resall$dataset)

summary(rateall)
write.table(format(rateall, digits=2), file = file.path(dir, "weighted_mean_rate.txt"), row.names = F, col.names = T, quote = F)

# Compute weighted mean of mutation rates
weighted.mean(resall$mut_rate, resall$weight)
var(resall$mut_rate)
#sd(resall$mut_rate)



############# PLot parameter posteriors (extended Fig 14b and Suppl Fig 1) ##################
pal = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')

fsuffix = paste0(suffix, "_epsilon", epsilon)
#fsuffix = "_final_mixed"

# plot for joint mu used in suppl doc (Suppl Fig. 1)
plot_rate_all(rateall, resall, dir, pal, fsuffix, ftype)

#plot_mu_brate(rateall, resall, dir, pal, fsuffix, ftype)

# used in extend Fig (14b)
plot_by_group(rateall, resall, dir, pal, fsuffix, ftype)
#plot_by_group(rateall, resall, dir, pal, fsuffix, ".tiff")


############## Compare inferences between models (validation)  ##############
res_sel = get_res_all_withsel(dir_sel, suffix, epsilon)
resall_sel = res_sel$resall
rateall_sel = res_sel$rateall
resall_sel$model = "selection"
unique(resall_sel$dataset)
summary(resall_sel)

resall$model = "neutral"
resall$fitness = 0
resall %>% filter(type=="joint") %>% select(names(resall_sel)) -> res_joint
head(resall_sel)
head(res_joint)
res_cmb = rbind(resall_sel, res_joint)

res_cmb$model = factor(res_cmb$model, levels=c("selection", "neutral"))
# + stat_compare_means(label = "p.signif")
p1=ggplot(res_cmb, aes(x=model, y=chr_rate, fill=model)) + geom_boxplot(width=0.6) + facet_grid(. ~ dataset, scales = "free") + theme_pubr() + theme(legend.position = "top", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("chromosomal CNA rate (per division)")+ scale_fill_npg()

p2=ggplot(res_cmb, aes(x=model, y=arm_rate, fill=model)) + geom_boxplot(width=0.6) + facet_grid(. ~ dataset, scales = "free") + theme_pubr() + theme(legend.position = "top", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("sub-chromosomal CNA rate (per division)") + scale_fill_npg()

p3=ggplot(res_cmb, aes(x=model, y=division_rate, fill=model)) + geom_boxplot(width=0.6) + facet_grid(. ~ dataset, scales = "free") + theme_pubr() + theme(legend.position = "top", axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + ylab("cell division rate (per day)") + scale_fill_npg()

ggarrange(p1, p2, p3, ncol = 1, common.legend = T)

fout=file.path(bdir, paste0("cmp_param_smc", ftype))
ggsave(fout, width = 8, height = 13)


############ Plot epsilon (validation) ####################
type = "neutral"
dir = file.path(bdir, type)
fprefix = "smc-joint-"
prefix=""
files=list.files(dir, pattern = paste0(fprefix, ".*[0-9]", prefix,"$"))
files
res=get_epsilon(dir, files, prefix, fprefix, ftype)
res_epsilons_neutral = res$res_epsilons_all
res_epsilons_neutral$type = type
nsim_neutral = res$nsim_all
nsim_neutral$type = type


type = "selection"
dir = file.path(bdir, type)
files=list.files(dir, pattern = paste0(fprefix, ".*[0-9]", prefix,"$"))
files
res=get_epsilon(dir, files, prefix, fprefix, ftype)
res_epsilons_sel = res$res_epsilons_all
res_epsilons_sel$type = "selection"
nsim_sel = res$nsim_all
nsim_sel$type = "selection"

# Plot merged epsilons, TODO: change dataset name if needed
res_epsilons = rbind(res_epsilons_sel, res_epsilons_neutral)
mindex = max(res_epsilons$index)
ggplot(res_epsilons, aes(x = (index), y = epsilon, color = type)) + geom_point() + geom_line() + theme_pubr() + xlab("population index") + scale_x_continuous(breaks = seq(1, max(mindex), 1)) + facet_grid(.~patient) + scale_color_npg() + labs(color = "model")
fout = file.path(bdir, paste0("epsilon_merged_", prefix, ftype))
ggsave(fout, width = 12, height = 6)

# write the number of total samples
nsim_all = rbind(nsim_sel, nsim_neutral)
nsim_all %>% spread(key = type, value = nsim) -> nsim_wide
fmidfix=prefix
fout = file.path(bdir, paste0("nsim_all", fmidfix,".tsv"))
write_tsv(nsim_wide, fout)
