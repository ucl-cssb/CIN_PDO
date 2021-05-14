library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)

# Check DIC on real organoid data

bdir = "../data/dic/"  # base directory of the results

get_res_neutral <- function(bdir, suffix){
  files = list.files(bdir, pattern = paste0("dic_neutral.*", suffix))
  print(files)
  res_all_neutral = data.frame()
  for (i in 1:length(files)){
    fname = file.path(bdir, files[i])
    res = read.table(fname)
    names(res)=c("run", "dataset", "ncell", "epsilon", "dic_neutral")
    res_all_neutral = rbind(res_all_neutral, res)
  }
  dim(res_all_neutral)
  head(res_all_neutral)

  return(res_all_neutral)
}

get_res_withsel <- function(bdir, suffix){
  files = list.files(bdir, pattern = paste0("dic_withsel_.*", suffix))
  print(files)
  res_all_withsel = data.frame()
  for (i in 1:length(files)){
    #i=1
    fname = file.path(bdir, files[i])
    res = read.table(fname)
    names(res)=c("run", "dataset", "ncell", "epsilon", "dic_withsel")
    res_all_withsel = rbind(res_all_withsel, res)
  }
  dim(res_all_withsel)
  head(res_all_withsel)

  return(res_all_withsel)
}


plot_dic_both <- function(res_all, bdir, suffix, ftype = ".pdf"){
  res_all %>% gather(dic_neutral, dic_withsel, key="model", value="DIC") -> res_gathered
  res_gathered$dataset = gsub("dataset_", "", res_gathered$dataset)
  res_gathered %>% mutate(dataset = if_else(dataset=="20200403", "20200304", dataset)) -> res_gathered
  res_gathered$model = gsub("dic_neutral", "neutral model", res_gathered$model)
  res_gathered$model = gsub("dic_withsel", "selection model", res_gathered$model)
  res_gathered$model = factor(res_gathered$model, levels=c("selection model", "neutral model"))

  ggplot(res_gathered, aes(x=dataset, y=DIC, fill=model))  +
    #geom_violin(scale = "width") +
    geom_boxplot(width = 0.5) +
    ylab("DIC") +
    theme_pubr()  + scale_fill_npg() + guides(fill = guide_legend(nrow = 1)) + theme_pubr()  + xlab("") + scale_y_continuous(n.breaks = 6) + scale_x_discrete(labels = c("PDTO-9 #1\n(Ext. Dat. Fig. 11)", "PDTO-9 #4\n(Fig. 3a)", "PDTO-9 #5\n(Fig. 2c)", "PDTO-9 #7\n(Fig. 4a)")) + labs(fill = "")

  prefix = paste0("sc_dic_boxplot", suffix)

  fout=file.path(bdir, paste0(prefix, ftype))
  ggsave(fout, width = 8, height = 5)

  # fout=file.path(bdir, paste0(prefix, ".png"))
  # ggsave(fout, width = 8, height = 5)
}




weight = "2.0"
epsilon = 0.2
suffix = paste0("_logs3_weight", weight)
suffix = paste0(suffix, "_mut0")
suffix = paste0(suffix, "_espilon", epsilon)

res_all_neutral = get_res_neutral(bdir, suffix)
res_all_withsel = get_res_withsel(bdir, suffix)
res_all = merge(res_all_neutral, res_all_withsel, by=c("dataset", "ncell", "run"))


# remove results for 20190409_2
res_all = res_all %>% filter(dataset!="dataset_20190409_2")
res_all %>% group_by(dataset) %>% tally()

epsilon = 0.1
suffix = paste0("_logs3_weight2_")

res_all_neutral = get_res_neutral(bdir, suffix)
res_all_withsel = get_res_withsel(bdir, suffix)
res_all2 = merge(res_all_neutral, res_all_withsel, by=c("dataset", "ncell", "run"))
res_20190409_2 = res_all2 %>% filter(dataset=="dataset_20190409_2")

res_all = rbind(res_all, res_20190409_2)

suffix="_final_mixed"
ftype=".png"
plot_dic_both(res_all, bdir, suffix, ftype)
