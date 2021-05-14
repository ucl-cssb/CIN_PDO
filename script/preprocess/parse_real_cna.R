library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)


# This script is used to reformat the real PTDO CNA data in a file suitable for ABC model fitting


################### Util function #####################
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-several-columns-on-a-subset-of-rows
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


################ Read input #########################
# read band data to determine arm boundary (hg19)
get_cyto_band <- function(fcyto){
  cyto=read.table(fcyto, header=FALSE)
  names(cyto) <- c("chrom","start","end","cytoband", "stain")

  # Find the boundary between two arms
  arm_boudaries=data.frame()
  chr_ends=data.frame()
  for (i in 1:(nrow(cyto)-1)) {
    chr1 = cyto[i,]$chrom
    arm1 = substr(cyto[i,]$cytoband, 1, 1)
    chr2 = cyto[i+1,]$chrom
    arm2 = substr(cyto[i+1,]$cytoband, 1, 1)

    if(chr1 == chr2 && arm1 != arm2){
      item=data.frame(chrom=chr1, boundary = cyto[i+1,]$start)
      arm_boudaries = rbind(arm_boudaries, item)
      print(paste("chrom", chr1, "arm boundary", cyto[i+1,]$start))
    }
    if(chr1 != chr2){
      item=data.frame(chrom=chr1, cend = cyto[i,]$end)
      chr_ends=rbind(chr_ends, item)
    }
    if(i==nrow(cyto)-1){
      item=data.frame(chrom=chr2, cend = cyto[i+1,]$end)
      chr_ends=rbind(chr_ends, item)
    }
  }

  merge(chr_ends, arm_boudaries, by="chrom") -> chr_end_arm

  return(chr_end_arm)
}


# Read real CNAs
get_real_data <- function(dir, fname, format="cn"){
  # try different formats
  # a vector of copy numbers
  if(format=="cn"){
    in_file = file.path(dir, fname)
    d_real <- read.table(in_file, header=FALSE)
    names(d_real) <- c("chrom", "start", "end", "ratio", "cn", "type", "sample", "size")
  }
  else{
    # a vector of frequencies
    in_file = file.path(dir, fname)
    d_real <- read.table(in_file, header=FALSE)
    names(d_real) <- c("chrom", "start", "end", "sample", "type")
  }

  return(d_real)
}


############### write output files ########################
write_clonal_events <- function(d_clonal, ploidy = 2){
  # Only output unique clonal events for reference
  d_clonal %>% select(chrom, arm, cn) %>% filter(cn!=ploidy) %>% unique() -> dc_uniq

  fout = file.path(dir,gsub(".txt", "_clonal.txt", fname))
  write.table(dc_uniq, fout, row.names = F, quote = F, col.names = F)
}


write_subclonal_events <- function(d_subclonal, fsummary, cutoff_frac = 0.5){
  sink(fsummary, append = T)

  d_subclonal %>% group_by(chrom, arm, cn) %>% count() -> dm_subn

  dm_subn %>% filter(arm=="whole") %>% nrow() -> nchr
  cat(paste0("\nnumber of subclonal chr-level events: ", nchr, "\n"))

  dm_subn %>% filter(arm!="whole") %>% nrow() -> narm
  cat(paste0("\nnumber of subclonal arm-level events: ", narm, "\n"))

  dm_subn %>% filter(n==1) %>% nrow() -> n1
  cat(paste0("\nnumber of CNA events appearing in one cell: ", n1, "\n"))

  dm_subn %>% filter(n>1) %>% arrange(n) %>% as.data.frame() -> tbl_cna
  if(nrow(tbl_cna) > 0){
    cat("\nCNA events appearing in more than one cell: \n")
    print(tbl_cna)
  }

  dm_subn %>% filter(arm==0) %>% filter(n>1) %>% as.data.frame() -> tbl_chr
  if(nrow(tbl_chr) > 0){
    cat("\nchr-level events appearing in more than one cell: \n")
    print(tbl_chr)
  }

  dm_subn %>% filter(arm!=0) %>% filter(n>1) %>% as.data.frame() -> tbl_arm
  if(nrow(tbl_arm) > 0){
    cat("\narm-level events appearing in more than one cell: \n")
    print(tbl_arm)
  }

  # d_subclonal %>% group_by(chrom) %>% count()
  # d_subclonal %>% group_by(sample) %>% count()
  # d_subclonal %>% select(chrom, arm, type, cn) %>% unique() %>% nrow()

  d_subclonal %>% filter(prop < cutoff_frac) %>% select(chrom, arm, start, end, boundary, cend, cn, prop) %>% unique() -> dmc_small
  dmc_small %>% nrow() -> num_small
  d_subclonal %>% filter(prop > 1)  %>% unique() %>% nrow() -> num_big
  cat(paste("\nThere are", num_small, "unique subclonal subchr-level CNAs which is of <20% of the arm\n"))
  print(dmc_small)
  cat(paste("\nThere are", num_big, "unique subclonal subchr-level CNAs which span two arms\n"))

  d_subclonal %>% select(sample, chrom, arm, cn) %>% unique() -> dsc_uniq
  fout = file.path(dir,gsub(".txt", "_subclonal.txt", fname))
  write.table(dsc_uniq, fout, row.names = F, quote = F, col.names = F)

  sink()
}


#################### Processsing input #############################
# Find cells with de novo MP-like events
get_mp_samples <- function(d_subclonal, fsummary, cutoff_mp = 10){
  sink(fsummary, append = T)
  # Exclude CNAs likely to be caused by MP. There may be complicated CNAs. For convenience, all CNVs appearing in the two cells are excluded.
  # d_subclonal %>% select(sample) %>% unique() -> sample_all
  d_subclonal %>% group_by(sample) %>% tally() %>% filter(n > cutoff_mp) %>% select(sample) %>% as.data.frame() -> sample_mp
  if(nrow(sample_mp) > 0){
    msg = paste0("\nThere are ", nrow(sample_mp), " cells with MP-like CNAs\n")
    cat(msg)
    if(nrow(sample_mp) > 0){
      print(sample_mp)
    }

    d_subclonal %>% filter(sample %in% sample_mp$sample) %>% select(chrom, arm, cn) %>% unique() -> cn_mp_like
    #print(cn_mp_like)
    cn_mp_like %>% filter(arm!=0) %>% nrow() -> n_arm
    cn_mp_like %>% filter(arm==0) %>% nrow() -> n_chr
    msg = paste0("\nThere are ", n_chr, " unique chr-level CNAs in the cells with MP-like CNAs\n")
    cat(msg)
    msg = paste0("\nThere are ", n_arm, " unique subchr-level CNAs in the cells with MP-like CNAs\n")
    cat(msg)
  }

  sink()

  return(sample_mp$sample)
}


# Extract subclonal mutation (not correctly recover known clonal events due to heterogeneity of the data)
separate_clonal_events <- function(d_real, chr_end_arm, cutoff = 0.5){
  samples = as.character(unique(d_real$sample))
  nsample <- length(samples)
  fsummary=file.path(dir,gsub(".txt", "_summary.txt", fname))
  sink(fsummary)

  cat(paste("There are", nsample, "cells in the sample\n"))

  # assign arm types with cytoband data
  dm = merge(d_real, chr_end_arm, by="chrom")

  dm$chrom=substr(dm$chrom, 4, length(d_real$chrom))
  dm %>% mutate(chrom = replace(chrom, chrom=="X","23")) -> dm

  # Exclude unrealible callings in PDO9T
  dm %>% filter(!(chrom==17|chrom==19)) -> dm

  dm$arm=""
  dm %>% mutate(arm = case_when(start>boundary & end - start  < 0.8 * cend ~ "q",
                                end - start  >= 0.8 * cend ~ "whole",
                                start<boundary & end - start  < 0.8 * cend ~ "p"
  )) -> dm
  dm %>% filter(arm=="") %>% nrow() -> num_no_arm
  if(num_no_arm > 0){
    stop("error in arm assigning!")
  }

  dm$prop = 0
  dm %>% mutate(prop = ifelse(arm=="p", (end - start ) / boundary, prop)) %>% mutate(prop = ifelse(arm=="q", (end - start ) / (cend-boundary), prop)) %>% mutate(prop = ifelse(arm=="whole", (end - start ) / (cend), prop)) -> dm

  # count the number of events to separate clonal events
  dm %>% group_by(chrom, arm, type, cn) %>% count() -> mcount
  mcount %>% mutate(clonality=ifelse(n >= cutoff * nsample, "clonal", "subclonal")) -> mcount
  dmc = merge(dm, mcount, by=c("chrom", "arm", "type", "cn"))

  # Extract clonal CNAs in the dataset
  dmc %>% filter(clonality=="clonal") %>% select(chrom, arm, type, cn) %>% unique() %>% mutate(arm = case_when(arm == "whole" ~ 0, arm == "p" ~ 1, arm == "q" ~ 2)) %>% group_by(chrom, arm, type) %>% summarise(mean_cn = Mode(cn))  -> d_clonal
  write_clonal_events(d_clonal)

  # Find de novo CNAs
  d_subclonal = dmc %>% filter(clonality=="subclonal") %>% select(sample, chrom, arm, start, end, prop, type, cn) %>% mutate(arm = case_when(arm == "whole" ~ 0, arm == "p" ~ 1, arm == "q" ~ 2))
  write_subclonal_events(d_subclonal, fsummary)

  return((list(subclonal=d_subclonal, clonal=d_clonal, nsample = nsample)))
}

# Assign arms and normal CNs to real CNs
get_cn_annot <- function(d_real, chr_end_arm, ploidy = 2){
  d_real_ext = merge(d_real, chr_end_arm, by = "chrom", all = T)

  nsample = length(unique(d_real$sample))
  d_real_ext$type = as.character(d_real_ext$type)
 # d_real_ext$sample = as.character(d_real_ext$sample)

  # Add records for samples without CNA regions (for visualization)
  #d_real_ext %>% filter(is.na(sample))
  d_real_ext %>% mutate(type=ifelse(is.na(type), "normal", type)) %>% mutate(start=ifelse(is.na(start), 1, start), end=ifelse(is.na(end), cend, end), ratio=ifelse(is.na(ratio), 1, ratio), cn=ifelse(is.na(cn), ploidy, cn), size=ifelse(is.na(size), cend, size), freq = ifelse(is.na(sample), nsample, 1)) %>% uncount(freq) -> d_real_ext_na
  nchr_normal = length(unique(d_real_ext$chrom)) - length(unique(d_real$chrom))
  d_real_ext_na %>% mutate_cond(is.na(sample), sample=rep(unique(d_real$sample), nchr_normal)) -> d_real_ext_sample
  #print("Check assignment of records to normal chromsomes")
  norm_chrs = setdiff(d_real_ext$chrom, d_real$chrom)
  for(chr in norm_chrs){
    # print(chr)
    d_real_ext_sample %>% filter(chrom == chr) %>% distinct(sample) %>% nrow() -> nr
    if(nr != nsample){
      stop(paste0("error in assignment of records to normal chromsome ", chr))
    }
  }

  # Add records for remaining part of the chromsomome to plot the whole length
  d_real_all = d_real_ext_sample
  d_real_all %>% group_by(chrom) %>% mutate(mend = max(end)) %>% select(chrom, mend, cend) %>% filter(mend < cend) %>% unique() -> chr_icmpl
  if(nrow(chr_icmpl) > 0){
    d_real_ext_cmpl = d_real_ext_sample
    for(i in 1:nrow(chr_icmpl)){
      cmpl = d_real_ext[i,]
      cmpl %>% mutate(chrom = chr_icmpl[i,]$chrom, start = chr_icmpl[i,]$mend, end = chr_icmpl[i,]$cend, ratio = ploidy, cn = ploidy, type = "icmpl_normal", sample = NA) %>% mutate(size = end - start ) -> cmpl
      d_real_ext_cmpl = rbind(d_real_ext_cmpl, cmpl)
    }
    #d_real_ext %>% filter(sample=="TBD") %>% View()
    d_real_ext_cmpl %>% mutate(freq = ifelse(is.na(sample), nsample, 1)) %>% uncount(freq) -> d_real_ext_icmpl
    d_real_ext_icmpl %>% mutate_cond(is.na(sample), sample=rep(unique(d_real$sample), nrow(chr_icmpl))) -> d_real_all
  }

  # Assign arm type
  d_real_all %>% mutate(arm = case_when(start>boundary & end - start  < 0.8 * cend ~ "q",
                                               end - start  >= 0.8 * cend ~ "whole",
                                               start<boundary & end - start  < 0.8 * cend ~ "p"
  )) -> d_real_all
  d_real_all %>% filter(arm=="") %>% nrow() -> num_no_arm
  if(num_no_arm > 0){
    stop("error in arm assigning!")
  }
  # Compute the fraction of CNAs
  d_real_all %>% mutate(prop = ifelse(arm=="p", (end - start ) / boundary, 0), qmid = 0.5 * (cend + boundary)) %>% mutate(prop = ifelse(arm=="q", (end - start ) / (cend-boundary), prop)) %>% mutate(prop = ifelse(arm=="whole", (end - start ) / (cend), prop)) -> d_real_all
  # Reassign arm if proportion of p-arm-level CNA surpass the middle point of q-arm
  d_real_all %>% mutate(arm = ifelse(arm == "p" & end > qmid, "q", arm)) %>% mutate(prop = ifelse(arm=="q", (end - start ) / (cend-boundary), prop)) -> d_real_all

  #d_real_all %>% mutate(n.probes=0, chrom=str_replace(chrom, "chr", ""), cn = cn-floor(mean(cn)))  %>% select(sample = sample, chrom = chrom, arm, start.pos = start, end.pos = end, n.probes, mean=cn)  -> d_annot
  d_real_all$sample=gsub("X4.", "", d_real_all$sample)
  d_real_all %>% mutate(n.probes=0, chrom=str_replace(chrom, "chr", ""), rcn = cn-ploidy, sample = as.factor(str_replace(sample, "\\..*", "")))  %>% select(sample, chrom = chrom, arm, prop, start.pos = start, end.pos = end, cend, boundary, mean=rcn, cn = cn, prop)  -> d_annot

  return(d_annot)
}


# Only count de novo (subclonal) reciprocal events, which are better reflections of errors iccured in organoid growth
get_reciprocal_cn <- function(d_annot, fsummary, is_complete = F, chr2excl=c("17", "19"), cutoff_clone=0.5, cutoff_mp = 10, cutoff_frac = 0.5){
  # 1. Exclude chroms with unrealible CN calls
  d_annot %>% filter(! chrom %in% chr2excl) -> d_annot_excl

  # 2. Exclude clonal events (which appear in cutoff * nsample samples)
  nsample = length(unique(d_annot$sample))
  d_annot_excl %>% group_by(chrom, arm, cn) %>% summarize(ncna = n()) -> mcount
  mcount %>% mutate(clonality=ifelse(ncna >= cutoff_clone * nsample, "clonal", "subclonal")) -> mcount
  #mcount %>% mutate(freq = ifelse(arm=="whole", 3, 1)) %>% uncount(freq) %>% group_by(chrom, cn, ncna, clonality) %>% mutate_cond(arm=="whole", arm=c("whole", "p", "q")) -> mcount_ext
  d_annot_excl_addc = merge(d_annot_excl, mcount, by=c("chrom", "arm", "cn"))
  d_annot_excl_addc %>% filter(clonality=="subclonal") -> d_subclonal
  d_annot_excl_addc %>% filter(clonality=="clonal") -> d_clonal
  # d_clonal %>% select(chrom, arm, cn) %>% unique()
  # d_subclonal %>% select(chrom, arm, cn) %>% unique()
  # d_annot_excl_sum %>% select(chrom, arm, cn) %>% unique()

  # 3. Exclude MP-like samples
  mp_samples = get_mp_samples(d_subclonal, fsummary, cutoff_mp)
  d_subclonal %>% filter(! sample %in% mp_samples) -> d_subclonal_nomp

  # 4. Exclude short singleton regions
  # 4.1 Combine adjacent or close interval on one arm and take the average if possilbe
  d_subclonal_nomp %>% group_by(sample, chrom, arm) %>% summarise(cns = paste(cn, collapse=","), starts = paste(start.pos, collapse=","), ends = paste(end.pos, collapse=","), start = min(start.pos), end = min(end.pos), nevent = n(), cn = as.integer(mean(cn)), boundary = min(boundary), cend = min(cend)) -> d_annot_excl_sum
  d_annot_excl_sum %>% group_by(sample, chrom) %>% tally() %>% filter(n>3) %>% nrow() -> nerr
  if(nerr > 0){
    stop("Error in parsing CNAs, more than three types of events on one chr!")
  }
  # 4.2 Recompute the fraction of CNAs
  d_annot_excl_sum %>% mutate(prop = ifelse(arm=="p", (end - start ) / boundary, 0), qmid = 0.5 * (cend + boundary)) %>% mutate(prop = ifelse(arm=="q", (end - start ) / (cend-boundary), prop)) %>% mutate(prop = ifelse(arm=="whole", (end - start ) / (cend), prop)) -> d_annot_excl_sum
  # Reassign arm if proportion of p-arm-level CNA surpass the middle point of q-arm
  d_annot_excl_sum %>% mutate(arm = ifelse(arm == "p" & end > qmid, "q", arm)) %>% mutate(prop = ifelse(arm=="q", (end - start ) / (cend-boundary), prop)) -> d_annot_excl_sum

  # Count the number of unique events on each location
  d_annot_excl_sum %>% group_by(chrom, arm) %>% summarize(nloc = n()) -> loc_count
  d_subclonal_addn = merge(d_annot_excl_sum, loc_count, by=c("chrom", "arm"))
  # Exclude singleton regions that are too small
  # d_subclonal_addn %>% filter((prop < cutoff_frac & nloc < 2))
  d_subclonal_addn %>% filter(!(arm!="whole" & prop < cutoff_frac & nloc < 2)) -> d_subclonal_reciprocal

  # Exclude all singleton regions if all cells are sequenced (not so useful with only one such dataset)
  if(is_complete){
    #d_subclonal_reciprocal %>% filter((nloc < 2))
    d_subclonal_reciprocal %>% filter(nloc > 1) -> d_subclonal_reciprocal
  }

  return(list(subclonal = d_subclonal_reciprocal, clonal = d_clonal))
}



###################### Plot CNA heatmap ##########################
plot_cn_heatmap_real <- function(d_seg, ploidy, fout, main, figtype){
  d_seg$pos = (d_seg$end.pos + d_seg$start.pos)/2
  d_seg$width = d_seg$end.pos - d_seg$start.pos
  d_seg$cn[d_seg$cn>5] = 5
  d_seg$cn = as.factor(d_seg$cn)
  d_seg$chrom = factor(d_seg$chrom, levels=paste("",c(c(1:22),'X','Y'),sep=""))

  colors = c("#6283A9","#bdd7e7","#f0f0f0","#FCAE91", "#B9574E", "#76000D")
  cns = c("0", "1", "2", "3", "4", "5")
  bcol = colors[ploidy + 1]

  plot<- ggplot(data=d_seg, aes(x=pos,y=sample,fill=cn)) +
    geom_tile(aes(width=width)) +
    facet_grid(sample~chrom, scales="free", space = "free", switch='y') +
    scale_color_manual(values=colors, limits=cns) +
    scale_fill_manual(values=colors, limits=cns) +
    theme(legend.position = "none",
          # strip.text.x = element_blank(),
          strip.text.y = element_blank(),
          #strip.text.y =  element_text(angle=180, size=2),
          strip.background = element_blank(),
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.spacing.y=unit(0, "lines"),
          panel.spacing.x=unit(0, "lines"),
          panel.border = element_rect(color = "#d8d8d8", fill = NA, size = .2),
          panel.background=element_rect(fill = bcol)
    ) +
    scale_y_discrete(breaks=unique(d_seg$sample[!is.na(d_seg$sample)]), expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))
  #plot

  # if(figtype == ".png"){
  #   png(fout, width = 700, height = 400)
  # }
  # if(figtype == ".tiff"){
  #   tiff(fout, units="in", width=7, height=4, res=300, compression = 'lzw')
  # }
  plot = plot + ggtitle(main)
  print(plot)
  ggsave(fout, width = 12, height = 8)
  # dev.off()
}


######### Compute summary statistics on CNAs ##############
# Extract unique subclonal CNAs
get_uniq_cn_vec <- function(d_clonal, d_subclonal, ploidy, fsummary){
  sink(fsummary, append = T)

  d_subclonal %>% select(chrom, arm, cn) %>% unique() -> dmsc_uniq

  dmsc_uniq %>% filter(arm!="whole") %>% nrow() -> n_arm
  dmsc_uniq %>% filter(arm=="whole") %>% nrow() -> n_chr
  msg = paste0("\nThere are ", n_chr, " unique chr-level CNAs in the selected cells\n")
  cat(msg)
  msg = paste0("\nThere are ", n_arm, " unique subchr-level CNAs in the selected cells\n")
  cat(msg)

  # compute relative copy number changes
  d_clonal %>% select(chrom, arm, cn) %>% unique() -> d_clonal_unique
  # Add arm-level CNs for those at whole chromosomes
  d_clonal_unique %>% mutate(freq = ifelse(arm=="whole", 3, 1)) %>% uncount(freq)  %>% group_by(chrom) %>% mutate_cond(arm=="whole", arm=c("whole", "p", "q")) -> d_clonal_ext

  d_clonal_ext %>% group_by(chrom, arm) %>% summarise(clonal_cn = as.integer(mean(cn))) -> d_clonal_ext_sum
  d_clonal_ext_sum %>% group_by(chrom) %>% tally() %>% filter(n>3) %>% nrow() -> nerr
  if(nerr > 0){
    stop("Error in parsing CNAs, more than three types of events on one chr!")
  }

  dmsc_cmb = merge(dmsc_uniq, d_clonal_ext_sum, all.x=T) %>% mutate(clonal_cn=replace_na(clonal_cn, ploidy))
  dmsc_cmb$rcn = (dmsc_cmb$cn - dmsc_cmb$clonal_cn)
  dmsc_cmb = dmsc_cmb %>% filter(rcn != 0)

  fout = file.path(dir,gsub(".txt", "_uvec.txt", fname))
  write.table(dmsc_cmb, fout, row.names = F, quote = F, col.names = T)

  # add up absolute relative CNs
  dmsc_cmb %>% mutate(arm=ifelse(arm=="whole", 0, 1)) %>% group_by(chrom, arm) %>% summarise(tcn=sum(abs(rcn))) -> dmsc_tcn
  dmsc_tcn %>% group_by(arm) %>% summarise(ttcn=sum(tcn)) %>% arrange(arm) -> dmsc_ttcn
  dmsc_ttcn %>% filter(arm==0) %>% nrow() -> n_chr_cna
  dmsc_ttcn %>% filter(arm==1) %>% nrow() -> n_arm_cna
  rcn = dmsc_ttcn$ttcn / nsample
  if(n_arm_cna==0){
    rcn=c(rcn, 0.0)
  }
  if(n_chr_cna==0){
    rcn=c(0.0, rcn)
  }

  fout = file.path(dir, gsub(".txt", "_rvec.txt", fname))
  write(rcn, file=fout, ncolumns = 1)

  sink()
}


############### Parse real dataset ################
# assuming work directory is CIN_PDO
# setwd("/Users/ucbtlux/Gdrive/git/cnv_analysis/CIN_PDO")
bdir = getwd()
bdir

format="cn"

cutoff_clone = 0.5
cutoff_mp = 10
chr2excl = c("8", "17", "19", "20")
cutoff_frac = 0.8  # Used to exclude singleton arm-level events


fcyto = file.path(bdir, "data/cytoBand_hg19.txt")
chr_end_arm = get_cyto_band(fcyto)


############### Parse single cell seq data ################
dir = file.path(bdir, "data/trees/4trees")
# list all files in the dir
infiles1 = list.files(dir, pattern = ".*[0-9]+.txt")
infiles1
infiles2 = list.files(dir, pattern = "sstat")
infiles2
infiles = setdiff(infiles1, infiles2)
infiles

################## Plot original callings of CNs ##################
plot_cn_all <- function(dir, infiles, format, figtype){
  for(i in 1:length(infiles)){
    print(i)
    fname = infiles[i]
    print(fname)

    d_real = get_real_data(dir, fname, format)

    ploidy = floor(mean(d_real$cn))
    d_annot = get_cn_annot(d_real, chr_end_arm, ploidy)
    #d_real %>% filter(chrom == "chr4") %>% select(start, cn) %>% unique()
    #d_annot %>% filter(chrom == 4) %>% select(arm, cn) %>% unique()
    fout=str_replace(file.path(dir, fname), ".txt", figtype)
    plot_cn_heatmap_real(d_annot, ploidy, fout, "", figtype)
  }
}

figtype=".png"
plot_cn_all(dir, infiles, format, figtype)


##################  Postprocess for ABC inference ##################
options(digits = 10)
# https://stackoverflow.com/questions/62140483/how-to-interpret-dplyr-message-summarise-regrouping-output-by-x-override
options(dplyr.summarise.inform=F)

# not work as a function, probably due to use of sink()
for(i in 1:length(infiles)){
  #i = 1
  print(i)
  fname = infiles[i]
  print(fname)

  is_complete = F
  if(fname=="dataset_20190409_1.txt"){
    is_complete = T
  }

  d_real = get_real_data(dir, fname, format)
  nsample = length(unique(d_real$sample))
  ploidy = floor(mean(d_real$cn))

  d_annot = get_cn_annot(d_real, chr_end_arm, ploidy)

  fsummary=file.path(dir,gsub(".txt", "_summary.txt", fname))
  sink(fsummary)
  cat(paste("There are", nsample, "cells in the sample\n"))
  cat(paste("\nThe ploidy is", ploidy, "\n"))

  res = get_reciprocal_cn(d_annot, fsummary, is_complete, chr2excl, cutoff_clone, cutoff_mp, cutoff_frac)
  d_subclonal = res$subclonal


  d_clonal = res$clonal
  write_clonal_events(d_clonal, ploidy)
  write_subclonal_events(d_subclonal, fsummary, cutoff_frac)

  sink()

  # summary statistics for ABC
  get_uniq_cn_vec(d_clonal, d_subclonal, ploidy, fsummary)
}
