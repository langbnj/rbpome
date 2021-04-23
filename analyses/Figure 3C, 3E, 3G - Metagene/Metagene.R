rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

source("~/Documents/R/mysqltc.R")
library(tidyverse)
library(reshape2)
library(glue)
library(broom)
library(scales)
library(magrittr)

# Correlation plot packages
library(corrplot)
library(Hmisc)
library(ggcorrplot)

# Plotting packages
library(ggbeeswarm)
library(ggrepel)



setwd("~/Documents/Projects/RBPome/Metagene")

# Real data
clip_region <- "clip_region"
clip_cobinding <- "clip_cobinding"

# # Random peaks
# clip_region <- "clip_region_random"
# clip_cobinding <- "clip_cobinding_random"

# Sequential random peaks (at every base of the main transcript)
# clip_region <- "clip_region_random_sequential"
# clip_cobinding <- "clip_cobinding_random_sequential"

style <- "meta"
# bincount <- 8
bincount <- 1
# type <- "eclip_encode"
# type <- "eclip_encode_12"
type <- "t_lrt_ihw_nocorr_auto_w100_s5"
species <- "human"

# Minimum number of co-bound peaks for evaluating an RBP pair (previously no minimum)
# min_cobound_peaks <- 100
min_cobound_peaks <- 1
# min_cobound_peaks <- 0

# When there used to be two levels cobound confidence levels:
# # Cobinding confidence level: either use the loose cutoff (≤117 nt, cobound=1) or the strict cutoff (≤54 nt, cobound=2)
# # cobinding_hc_level <- 1 # ≤117 nt (loose, cobound=1)
# cobinding_hc_level <- 2 # ≤54 nt  (strict, cobound=2)

# Now there's one:
cobinding_hc_level <- 1 # ≤54 nt (strict, cobound=1)

# Require cobinding to be within the same celltype
# cobinding_merge_celltypes <- 1 # Previously, no
cobinding_merge_celltypes <- 0 # Now, yes

# Resamples
# num_resamples <- 3 # Fastest
# num_resamples <- 100 # Faster
num_resamples <- 1000 # Fast
# num_resamples <- 10000 # More accurate



# # Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course)
# # SELECT symbol1, symbol2 FROM rbpome WHERE scoretype='SUM_IS' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2;
# q <- Query("SELECT symbol1, symbol2 FROM rbpome WHERE scoretype='SUM_IS' AND eclip1=1 AND eclip2=1 AND hc=1 AND homodimer=0 GROUP BY symbol1, symbol2 ORDER BY symbol1, symbol2")
# rbp_pairs <- mapply(c, q$symbol1, q$symbol2, SIMPLIFY=F)  # Type conversion workaround
# # str(rbp_pairs)

# # Mildly filtered for cobinding (resampling p-value significant for at least one orientation)
# # Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course) (52: the big Cytoscape network, Figure 5a)
# q <- Query("SELECT protein_a_unique_pairs AS symbol1, protein_b_unique_pairs AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND (resampling_p_value_a_vs_b!='N/A' OR resampling_p_value_b_vs_a!='N/A') AND (resampling_p_value_a_vs_b<0.05 OR resampling_p_value_b_vs_a<0.05) GROUP BY protein_a_unique_pairs, protein_b_unique_pairs ORDER BY protein_a_unique_pairs, protein_b_unique_pairs")
# rbp_pairs <- mapply(c, q$symbol1, q$symbol2, SIMPLIFY=F)  # Type conversion workaround
# str(rbp_pairs)

# Strictly filtered for cobinding (resampling p-value and resampling Wilcoxon p-value significant in one orientation at least)
# Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course) (21: the "strict" Cytoscape network, Figure 5a)
q <- Query("SELECT protein_a AS symbol1, protein_b AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND (resampling_p_value_a_vs_b IS NOT NULL AND resampling_p_value_b_vs_a IS NOT NULL AND resampling_wilcoxon_p_value_a_vs_b IS NOT NULL AND resampling_wilcoxon_p_value_b_vs_a IS NOT NULL) AND ((resampling_p_value_a_vs_b<0.05 AND resampling_wilcoxon_p_value_a_vs_b<0.05) OR (resampling_p_value_b_vs_a<0.05 AND resampling_wilcoxon_p_value_b_vs_a<0.05)) GROUP BY protein_a, protein_b ORDER BY protein_a, protein_b")
# #DEBUG LIMIT 3
# q <- Query("SELECT protein_a_unique_pairs AS symbol1, protein_b_unique_pairs AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND (resampling_p_value_a_vs_b!='N/A' AND resampling_p_value_b_vs_a!='N/A' AND resampling_wilcoxon_p_value_a_vs_b!='N/A' AND resampling_wilcoxon_p_value_b_vs_a!='N/A') AND ((resampling_p_value_a_vs_b<0.05 AND resampling_wilcoxon_p_value_a_vs_b<0.05) OR (resampling_p_value_b_vs_a<0.05 AND resampling_wilcoxon_p_value_b_vs_a<0.05)) GROUP BY protein_a_unique_pairs, protein_b_unique_pairs ORDER BY protein_a_unique_pairs, protein_b_unique_pairs LIMIT 3")
# #END DEBUG
rbp_pairs <- mapply(c, q$symbol1, q$symbol2, SIMPLIFY=F)  # Type conversion workaround
str(rbp_pairs)

# # Ultra-strictly filtered for cobinding (resampling p-value significant in both orientations, and resampling Wilcoxon p-values significant in both orientations as well)
# # Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course) (21: the "strict" Cytoscape network, Figure 5a)
# q <- Query("SELECT protein_a_unique_pairs AS symbol1, protein_b_unique_pairs AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND (resampling_p_value_a_vs_b!='N/A' AND resampling_p_value_b_vs_a!='N/A' AND resampling_wilcoxon_p_value_a_vs_b!='N/A' AND resampling_wilcoxon_p_value_b_vs_a!='N/A') AND (resampling_p_value_a_vs_b<0.05 AND resampling_p_value_b_vs_a<0.05 AND resampling_wilcoxon_p_value_a_vs_b<0.05 AND resampling_wilcoxon_p_value_b_vs_a<0.05) GROUP BY protein_a_unique_pairs, protein_b_unique_pairs ORDER BY protein_a_unique_pairs, protein_b_unique_pairs")
# rbp_pairs <- mapply(c, q$symbol1, q$symbol2, SIMPLIFY=F)  # Type conversion workaround
# str(rbp_pairs)




# #DEBUG
# # Calculate only a subset of pairs to speed things up:
# # Selected pairs
# rbp_pairs <- list(
#   # Positive controls
#   # SELECT DISTINCT symbol1, symbol2 FROM rbpome WHERE bg_alldirect=1 AND homodimer=0 AND symbol1 IN ('AARS', 'AATF', 'ABCF1', 'AGGF1', 'AKAP1', 'AKAP8L', 'APOBEC3C', 'AQR', 'BCCIP', 'BCLAF1', 'BUD13', 'CDC40', 'CPEB4', 'CPSF6', 'CSTF2', 'CSTF2T', 'DDX21', 'DDX24', 'DDX3X', 'DDX42', 'DDX51', 'DDX52', 'DDX55', 'DDX59', 'DDX6', 'DGCR8', 'DHX30', 'DKC1', 'DROSHA', 'EFTUD2', 'EIF3D', 'EIF3G', 'EIF3H', 'EIF4G2', 'EWSR1', 'EXOSC5', 'FAM120A', 'FASTKD2', 'FKBP4', 'FMR1', 'FTO', 'FUBP3', 'FUS', 'FXR1', 'FXR2', 'G3BP1', 'GEMIN5', 'GNL3', 'GPKOW', 'GRSF1', 'GRWD1', 'GTF2F1', 'HLTF', 'HNRNPA1', 'HNRNPC', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRNPUL1', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'ILF3', 'KHDRBS1', 'KHSRP', 'LARP4', 'LARP7', 'LIN28B', 'LSM11', 'MATR3', 'METAP2', 'MTPAP', 'NCBP2', 'NIP7', 'NIPBL', 'NKRF', 'NOL12', 'NOLC1', 'NONO', 'NPM1', 'NSUN2', 'PABPC4', 'PABPN1', 'PCBP1', 'PCBP2', 'PHF6', 'POLR2G', 'PPIG', 'PPIL4', 'PRPF4', 'PRPF8', 'PTBP1', 'PUM1', 'PUM2', 'PUS1', 'QKI', 'RBFOX2', 'RBM15', 'RBM22', 'RBM5', 'RPS11', 'RPS3', 'SAFB', 'SAFB2', 'SBDS', 'SDAD1', 'SERBP1', 'SF3A3', 'SF3B1', 'SF3B4', 'SFPQ', 'SLBP', 'SLTM', 'SMNDC1', 'SND1', 'SRSF1', 'SRSF7', 'SRSF9', 'SSB', 'STAU2', 'SUB1', 'SUGP2', 'SUPV3L1', 'TAF15', 'TARDBP', 'TBRG4', 'TIA1', 'TIAL1', 'TRA2A', 'TROVE2', 'U2AF1', 'U2AF2', 'UCHL5', 'UPF1', 'UTP18', 'UTP3', 'WDR3', 'WDR43', 'WRN', 'XPO5', 'XRCC6', 'XRN2', 'YBX3', 'YWHAG', 'ZC3H11A', 'ZC3H8', 'ZNF622', 'ZNF800', 'ZRANB2') AND symbol2 IN ('AARS', 'AATF', 'ABCF1', 'AGGF1', 'AKAP1', 'AKAP8L', 'APOBEC3C', 'AQR', 'BCCIP', 'BCLAF1', 'BUD13', 'CDC40', 'CPEB4', 'CPSF6', 'CSTF2', 'CSTF2T', 'DDX21', 'DDX24', 'DDX3X', 'DDX42', 'DDX51', 'DDX52', 'DDX55', 'DDX59', 'DDX6', 'DGCR8', 'DHX30', 'DKC1', 'DROSHA', 'EFTUD2', 'EIF3D', 'EIF3G', 'EIF3H', 'EIF4G2', 'EWSR1', 'EXOSC5', 'FAM120A', 'FASTKD2', 'FKBP4', 'FMR1', 'FTO', 'FUBP3', 'FUS', 'FXR1', 'FXR2', 'G3BP1', 'GEMIN5', 'GNL3', 'GPKOW', 'GRSF1', 'GRWD1', 'GTF2F1', 'HLTF', 'HNRNPA1', 'HNRNPC', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRNPUL1', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'ILF3', 'KHDRBS1', 'KHSRP', 'LARP4', 'LARP7', 'LIN28B', 'LSM11', 'MATR3', 'METAP2', 'MTPAP', 'NCBP2', 'NIP7', 'NIPBL', 'NKRF', 'NOL12', 'NOLC1', 'NONO', 'NPM1', 'NSUN2', 'PABPC4', 'PABPN1', 'PCBP1', 'PCBP2', 'PHF6', 'POLR2G', 'PPIG', 'PPIL4', 'PRPF4', 'PRPF8', 'PTBP1', 'PUM1', 'PUM2', 'PUS1', 'QKI', 'RBFOX2', 'RBM15', 'RBM22', 'RBM5', 'RPS11', 'RPS3', 'SAFB', 'SAFB2', 'SBDS', 'SDAD1', 'SERBP1', 'SF3A3', 'SF3B1', 'SF3B4', 'SFPQ', 'SLBP', 'SLTM', 'SMNDC1', 'SND1', 'SRSF1', 'SRSF7', 'SRSF9', 'SSB', 'STAU2', 'SUB1', 'SUGP2', 'SUPV3L1', 'TAF15', 'TARDBP', 'TBRG4', 'TIA1', 'TIAL1', 'TRA2A', 'TROVE2', 'U2AF1', 'U2AF2', 'UCHL5', 'UPF1', 'UTP18', 'UTP3', 'WDR3', 'WDR43', 'WRN', 'XPO5', 'XRCC6', 'XRN2', 'YBX3', 'YWHAG', 'ZC3H11A', 'ZC3H8', 'ZNF622', 'ZNF800', 'ZRANB2') ORDER BY symbol1, symbol2;
#   # c("EWSR1", "PCBP1"),
#   # c("EWSR1", "SF3B4"),
#   # c("FMR1", "FXR2"),
#   # c("FXR1", "FXR2"),
#   c("LARP4", "FUBP3"),
#   # c("HNRNPK", "QKI"),
#   # c("HNRNPK", "U2AF2"),
#   # c("NONO", "SFPQ"),
#   # c("PCBP1", "PTBP1"),
#   # c("PCBP1", "QKI"),
#   # c("PCBP2", "QKI"),
#   # c("PTBP1", "QKI"),
#   c("U2AF1", "U2AF2")
#   # # Negative controls
#   # c("FMR1", "U2AF1"),
#   # c("FASTKD2", "SLBP"),
#   # c("U2AF1", "FMR1"),
#   # # Screen hits
#   # c("CSTF2T", "EWSR1"),
#   # c("NONO", "EWSR1")
# )
# #END DEBUG
# #DEBUG
# rbp_pairs <- list(
#   c("LARP4", "FUBP3"),
#   c("U2AF1", "U2AF2")
# )
# #END DEBUG
# #DEBUG
# rbp_pairs <- list(
#   c("PCBP1", "APOBEC3C"),
#   c("PCBP1", "RBFOX2"),
#   c("PCBP1", "HNRNPK"),
#   c("PCBP1", "PTBP1"),
#   c("PCBP1", "IGF2BP2")
# )
# #END DEBUG
# #DEBUG
# rbp_pairs <- list(
#   # c("LARP4", "FUBP3"),
#   c("U2AF1", "U2AF2"),
#   c("PCBP1", "APOBEC3C"),
#   c("PCBP1", "RBFOX2"),
#   c("PCBP1", "HNRNPK"),
#   c("PCBP1", "PTBP1"),
#   c("PCBP1", "IGF2BP2")
# )
# #END DEBUG
# #DEBUG
# rbp_pairs <- list(
#   c("PTBP1", "IGF2BP1"),
#   c("PTBP1", "IGF2BP2"),
#   c("PTBP1", "PCBP1"),
#   c("PTBP1", "QKI"),
#   c("PTBP1", "TARDBP"),
#   c("PTBP1", "ZC3H8")
# )
# #END DEBUG
#DEBUG
rbp_pairs <- list(
  c("SLBP", "SLBP")
)
#END DEBUG










# Initialise dataframes for final plots
final_heatmap <- data.frame("rbp1" = character(), "rbp2" = character(), "pair" = character(), "region" = factor(), "bin" = factor(), "value" = numeric())
final_heatmap_2 <- data.frame("rbp1" = character(), "rbp2" = character(), "pair" = character(), "region" = factor(), "bin" = factor(), "value" = numeric())
# final_heatmap_3 <- tibble()
final_heatmap_4 <- data.frame("rbp1" = character(), "rbp2" = character(), "pair" = character(), "region" = factor(), "bin" = factor(), "value" = numeric())
final_heatmap_5 <- data.frame("rbp" = character(), "region" = factor(), "bin" = factor(), "value" = numeric())
# alldensities <- list()
alldensities <- tibble(rbp1=character(), rbp2=character(), densities=list())

# Initialise dataframe for quantifying the differences between the complex and the individual RBPs (complex vs. max(rbp1, rbp2))
pair_distance <- data.frame("rbp1" = character(), "rbp2" = character(), "pair" = character(), "region" = factor(), "value" = numeric())


# #DEBUG
# # rbp1 <- "FMR1"
# # rbp2 <- "FXR2"
# # rbp_pair <- c("FMR1", "FXR2")
# # rbp_pairs <- list(c("FMR1", "FXR2"))
# rbp1 <- "EWSR1"
# rbp2 <- "NONO"
# rbp_pair <- c("EWSR1", "NONO")
# rbp_pairs <- list(c("EWSR1", "NONO"))
# rbp1 <- "PCBP1"
# rbp2 <- "IGF2BP2"
# rbp_pair <- c("PCBP1", "IGF2BP2")
# rbp_pairs <- list(c("PCBP1", "IGF2BP2"))
# rbp1 <- "RBM15"
# rbp2 <- "SRSF9"
rbp1 <- "FMR1"
rbp2 <- "FXR2"
# #END DEBUG








# Cycle through individual RBPs in rbp_pairs and get their data (regardless of close or distant binding compared to other RBPs)
rbps <- rbp_pairs %>% tibble() %>% rename(rbp=".") %>% unnest(cols=rbp) %>% unique() %>% arrange(rbp) %>% pull(rbp)
str(rbps)
rbp <- rbps[1]
for (rbp in rbps) {
  print(glue(" >> {rbp}"))

  q <- Query(paste0("SELECT CONCAT_WS('|', chr, start, stop, strand) AS peak, ensgv, region, pos FROM ",clip_region," WHERE style='",style,"' AND type='",type,"' AND species='",species,"' AND symbol='",rbp,"'"))
  
  # Assign "title" (the RBP name)
  q$title <- factor(rbp)

  # Reorder factors
  q$region <- factor(q$region, levels=c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"), ordered=T)
  
  # Rename regions to something presentable
  q$region <- plyr::mapvalues(q$region, from = c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"), to = c("5' UTR", "first\nexon", "first\nintron", "exon", "intron", "last\nexon", "last\nintron", "3' UTR"))
  
  
  # Make "pos" a continuous value across the entire metagene ([0..1])
  q$globalpos <- NA
  for (i in 1:length(levels(q$region))) {
    # print(paste0(i,": ",levels(q$region)[i]))
    # q[q$region==levels(q$region)[i],]$globalpos <- q[q$region==levels(q$region)[i],]$pos + (i-1)
    # Scale original pos to between 0.1 .. 0.9 to leave a blank area between regions:
    q[q$region==levels(q$region)[i],]$globalpos <- ((q[q$region==levels(q$region)[i],]$pos * 8/10) + 0.1) + (i-1)
  }
  
  summary(q$globalpos)
  q$globalpos <- q$globalpos / length(levels(q$region))
  summary(q$globalpos)
  
  
  # Assign histogram bins
  q$bin <- factor(cut(q$pos, breaks = 0:bincount/bincount, include.lowest = T, labels=1:bincount/bincount), ordered=T)
  levels(q$bin)
  str(q$bin)
  
  
  # Cycle through bins
  q2 <- data.frame("title" = factor(levels = levels(q$title)), "region" = factor(levels = levels(q$region)), "bin" = factor(levels = levels(q$bin)), "value" = numeric())
  # l <- list()
  # Get counts in a separate dataframe
  for (my_title in levels(q$title)) {
    for (my_region in levels(q$region)) {
      for (my_bin in levels(q$bin)) {
        # print(paste(my_title, my_region, my_bin))
        # str(data.frame(title = my_title, region = my_region, bin = my_bin))
        q2 <- rbind(q2, data.frame(title = my_title, region = my_region, bin = my_bin, value = nrow(q[q$title == my_title & q$region == my_region & q$bin == my_bin,])))
        # +1
        # q2 <- rbind(q2, data.frame(title = my_title, region = my_region, bin = my_bin, value = nrow(q[q$title == my_title & q$region == my_region & q$bin == my_bin,]) + 1))
        # q[q$title == my_title & q$region == my_region & q$bin == my_bin,]$cyan <-
      }
    }
  }
  q2
  # Test plot: Looks exactly the same as v2 (the histograms)
  ggplot(q2, aes(x=bin, y=value)) + geom_bar(width=1, stat="identity") + facet_grid(rows=vars(title), cols=vars(region), scales="free_y") + ylab("Peak count") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  q3 <- q2
  # Scale to between 0 and 1
  for (my_title in levels(q$title)) {
    for (my_region in levels(q$region)) {
      for (my_bin in levels(q$bin)) {
        # q2[q2$title == my_title & q2$region == my_region & q2$bin == my_bin,]$value
        # max(q2[q2$title == my_title,]$value)
        q3[q3$title == my_title & q3$region == my_region & q3$bin == my_bin,]$value <- q2[q2$title == my_title & q2$region == my_region & q2$bin == my_bin,]$value / max(q2[q2$title == my_title,]$value)
      }
    }
  }
  # ggplot(q3, aes(x=bin, y=value)) + geom_bar(width=1, stat="identity") + facet_grid(rows=vars(title), cols=vars(region), scales="free_y") + ylab("Peak count") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
  # Final
  # q4 <- data.frame("region" = factor(levels = levels(q$region)), "bin" = factor(levels = levels(q$bin)), "rbp1" = numeric(), "rbp2" = numeric(), "value" = numeric())
  for (my_region in levels(q$region)) {
    for (my_bin in levels(q$bin)) {
      
      # Utility values
      tmp_value_rbp <- q3[q3$region == my_region & q3$bin == my_bin,]$value

      # q4 <- rbind(q4, data.frame(region = my_region, bin = my_bin, rbp1 = rbp, rbp2 = rbp, value = tmp_value_rbp))

      # Another heatmap variant: Simply keep the individual and pair values
      final_heatmap_5 <- rbind(final_heatmap_5, data.frame(rbp=rbp, region=my_region, bin=my_bin, value=tmp_value_rbp))
    }
  }

  # q4
  # ggplot(q4, aes(x=bin, y=value)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(strip.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle(paste0(rbp1," & ",rbp2))
  # ggsave(paste0("output-plot-metagene-v6-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")
  
  
  
  
  
  
  
  # Get density data for individual RBPs (independent of complexes)
  
  all_ensgvs <- unique(q$ensgv)
  tib <- tibble()
  set.seed(1)
  for (i in 1:num_resamples) {
    
    # # Subsample to 50% (with replacement)
    # these_ensgvs <- sample(all_ensgvs, floor(length(all_ensgvs) / 2), replace=T)
    # Subsample 100% (with replacement)
    these_ensgvs <- sample(all_ensgvs, length(all_ensgvs), replace=T)
    
    this_tib <- subset(q, ensgv %in% these_ensgvs) %>% group_by(title, .drop=FALSE) %>% select(title, globalpos)
    this_tib$title <- as.character(this_tib$title)
    
    this_tib$rbp1 <- rbp1
    this_tib$rbp2 <- rbp2
    # this_tib$resample <- i
    this_tib
    
    # Append to tib
    if (nrow(tib) > 0) {
      tib <- rbind(tib, tibble(bs=i, data=list(this_tib)))
    } else {
      tib <- tibble(bs=i, data=list(this_tib))
    }
  }
  tib
  tib[1,2] %>% unnest(cols=data)
  # summary(tib)
  
  
  # Compute within-sample density
  densities.within <-
    tib %>%
    unnest(cols=data) %>% 
    # group_by(bs, region) %>%
    # group_by(bs, title, region) %>%
    group_by(bs, title) %>%
    # do(tidy(density(.$pos, 
    do(tidy(density(.$globalpos, 
                    from = 0,
                    to = 1,
                    # n = 128,
                    # bw=1/80
                    # bw=1/160
                    bw=1/320
    )))
  densities.within
  
  # Summarize densities into quantiles
  densities.qtiles <-
    densities.within %>%
    # rename(pos = x, dens = y) %>%
    rename(globalpos = x, dens = y) %>%
    ungroup() %>%
    # group_by(title, region, pos) %>%
    # group_by(title, region, globalpos) %>%
    group_by(title, globalpos) %>%
    # group_by(region, pos) %>% 
    summarise(q05 = quantile(dens, 0.025),
              q50 = quantile(dens, 0.5),
              q95 = quantile(dens, 0.975)) 
  densities.qtiles
  

  # Visualize and compare
  
  # Set title
  densities.qtiles$title <- paste0(rbp,"\n",length(unique(q$peak))," sites\n",length(unique(q$ensgv))," genes")
  # # Convert title to character
  # densities.qtiles$title <- as.character(densities.qtiles$title)
  
  # Store densities.qtiles for later
  alldensities %<>% add_row(rbp1=rbp, rbp2=NULL, densities=list(densities.qtiles))
  
  
  
  
  
  
}
final_heatmap_5
alldensities
alldensities %>% unnest(cols=densities)













# Cycle through RBP pairs (main part)
for (rbp_pair in rbp_pairs) {
  rbp1 <- rbp_pair[1]
  rbp2 <- rbp_pair[2]
  
  cat(paste0(" >> ",rbp1,"|",rbp2))
  
  # Get data
  
  # # Categorising by ENSTVs (not by peak proximity)
  # rbp1_only <- Query(paste0("SELECT r.* FROM ",clip_region," r LEFT OUTER JOIN (SELECT DISTINCT enstv FROM ",clip_region," WHERE symbol='",rbp2,"') t2 ON r.enstv=t2.enstv WHERE r.symbol='",rbp1,"' AND t2.enstv IS NULL"))
  # rbp2_only <- Query(paste0("SELECT r.* FROM ",clip_region," r LEFT OUTER JOIN (SELECT DISTINCT enstv FROM ",clip_region," WHERE symbol='",rbp1,"') t2 ON r.enstv=t2.enstv WHERE r.symbol='",rbp2,"' AND t2.enstv IS NULL"))
  # rbp1_rbp2 <- Query(paste0("SELECT r.* FROM ",clip_region," r, (SELECT DISTINCT enstv FROM ",clip_region," WHERE symbol='",rbp1,"') t1, (SELECT DISTINCT enstv FROM ",clip_region," WHERE symbol='",rbp2,"') t2 WHERE r.enstv=t1.enstv AND r.enstv=t2.enstv AND r.symbol='",rbp1,"'"))
  # str(rbp1_only)
  # str(rbp2_only)
  # str(rbp1_rbp2)
  
  # # Categorising by peak proximity
  # rbp1_only <- Query(paste0("SELECT r.acc, r.chr, r.ensg, r.ensgv, r.enst, r.enstv, r.pos, r.region, r.species, r.start, r.stop, r.strand, r.style, r.symbol, r.type, c.celltype, c.cobound, c.mindist, c.other_acc, c.other_symbol, c.rep FROM ",clip_region," r LEFT OUTER JOIN ",clip_cobinding," c ON r.symbol=c.symbol AND r.chr=c.chr AND r.start=c.start AND r.stop=c.stop AND r.strand=c.strand WHERE r.symbol='",rbp1,"' AND c.other_symbol='",rbp2,"' AND (cobound=0 OR cobound IS NULL)"))
  # rbp2_only <- Query(paste0("SELECT r.acc, r.chr, r.ensg, r.ensgv, r.enst, r.enstv, r.pos, r.region, r.species, r.start, r.stop, r.strand, r.style, r.symbol, r.type, c.celltype, c.cobound, c.mindist, c.other_acc, c.other_symbol, c.rep FROM ",clip_region," r LEFT OUTER JOIN ",clip_cobinding," c ON r.symbol=c.symbol AND r.chr=c.chr AND r.start=c.start AND r.stop=c.stop AND r.strand=c.strand WHERE r.symbol='",rbp2,"' AND c.other_symbol='",rbp1,"' AND (cobound=0 OR cobound IS NULL)"))
  # rbp1_rbp2 <- Query(paste0("SELECT r.acc, r.chr, r.ensg, r.ensgv, r.enst, r.enstv, r.pos, r.region, r.species, r.start, r.stop, r.strand, r.style, r.symbol, r.type, c.celltype, c.cobound, c.mindist, c.other_acc, c.other_symbol, c.rep FROM ",clip_region," r, ",clip_cobinding," c WHERE r.symbol='",rbp1,"' AND c.other_symbol='",rbp2,"' AND r.symbol=c.symbol AND r.chr=c.chr AND r.start=c.start AND r.stop=c.stop AND r.strand=c.strand AND cobound=",cobinding_hc_level))
  # str(rbp1_only)
  # str(rbp2_only)
  # str(rbp1_rbp2)
  
  # Categorising by peak proximity (faster queries)
  rbp1_peaks <- Query(paste0("SELECT CONCAT_WS('|', chr, start, stop, strand) AS peak, ensgv, region, pos FROM ",clip_region," WHERE style='",style,"' AND type='",type,"' AND species='",species,"' AND symbol='",rbp1,"'"))
  rbp2_peaks <- Query(paste0("SELECT CONCAT_WS('|', chr, start, stop, strand) AS peak, ensgv, region, pos FROM ",clip_region," WHERE style='",style,"' AND type='",type,"' AND species='",species,"' AND symbol='",rbp2,"'"))
  
  if (cobinding_merge_celltypes == 1) {
    rbp1_rbp2_peaks <- Query(paste0("SELECT DISTINCT CONCAT_WS('|', chr, start, stop, strand) AS peak FROM ",clip_cobinding," WHERE celltype IS NULL AND type='",type,"' AND species='",species,"' AND symbol='",rbp1,"' AND other_symbol='",rbp2,"' AND cobound=",cobinding_hc_level))
    rbp2_rbp1_peaks <- Query(paste0("SELECT DISTINCT CONCAT_WS('|', chr, start, stop, strand) AS peak FROM ",clip_cobinding," WHERE celltype IS NULL AND type='",type,"' AND species='",species,"' AND symbol='",rbp2,"' AND other_symbol='",rbp1,"' AND cobound=",cobinding_hc_level))
  } else {
    rbp1_rbp2_peaks <- Query(paste0("SELECT DISTINCT CONCAT_WS('|', chr, start, stop, strand) AS peak FROM ",clip_cobinding," WHERE celltype IS NOT NULL AND type='",type,"' AND species='",species,"' AND symbol='",rbp1,"' AND other_symbol='",rbp2,"' AND cobound=",cobinding_hc_level))
    rbp2_rbp1_peaks <- Query(paste0("SELECT DISTINCT CONCAT_WS('|', chr, start, stop, strand) AS peak FROM ",clip_cobinding," WHERE celltype IS NOT NULL AND type='",type,"' AND species='",species,"' AND symbol='",rbp2,"' AND other_symbol='",rbp1,"' AND cobound=",cobinding_hc_level))
  }
  
  str(rbp1_peaks)
  str(rbp2_peaks)
  str(rbp1_rbp2_peaks)
  str(rbp2_rbp1_peaks)
  
  # Get region information from the rbp1_peaks and rbp2_peaks dataframes
  str(subset(rbp1_peaks, peak %in% rbp1_rbp2_peaks$peak))
  str(subset(rbp2_peaks, peak %in% rbp2_rbp1_peaks$peak))
  rbp1_rbp2 <- rbind(subset(rbp1_peaks, peak %in% rbp1_rbp2_peaks$peak), subset(rbp2_peaks, peak %in% rbp2_rbp1_peaks$peak))
  # Filter out cobound peaks from rbp1_peaks and rbp2_peaks
  rbp1_only <- subset(rbp1_peaks, !(peak %in% rbp1_rbp2_peaks$peak))
  rbp2_only <- subset(rbp2_peaks, !(peak %in% rbp2_rbp1_peaks$peak))
  str(rbp1_only)
  str(rbp2_only)
  str(rbp1_rbp2)
  
  # Set row titles
  rbp1_only_title <- paste0(rbp1," only\n",length(unique(rbp1_only$peak))," sites\n",length(unique(rbp1_only$ensgv))," genes")
  rbp2_only_title <- paste0(rbp2," only\n",length(unique(rbp2_only$peak))," sites\n",length(unique(rbp2_only$ensgv))," genes")
  # rbp1_rbp2_title <- paste0(rbp1," and ",rbp2,"\n",nrow(rbp1_rbp2)," sites\n",length(unique(rbp1_rbp2$ensgv))," genes")
  # rbp1_rbp2_title <- paste0(rbp1," & ",rbp2,"\n",nrow(rbp1_rbp2)," sites\n",length(unique(rbp1_rbp2$ensgv))," genes")
  rbp1_rbp2_title <- paste0(rbp1," & ",rbp2,"\n",length(unique(rbp1_rbp2$peak))," sites\n",length(unique(rbp1_rbp2$ensgv))," genes")
  
  if (nrow(rbp1_only) > 0) { rbp1_only$title <- rbp1_only_title }
  if (nrow(rbp2_only) > 0) { rbp2_only$title <- rbp2_only_title }
  if (nrow(rbp1_rbp2) > 0) { rbp1_rbp2$title <- rbp1_rbp2_title }
  
  # Only continue if there are enough co-bound peaks for this pair
  if (nrow(rbp1_rbp2) >= min_cobound_peaks) {
    
    q <- rbind(rbp1_only, rbp2_only, rbp1_rbp2)
    
    # Reorder factors
    q$region <- factor(q$region, levels=c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"), ordered=T)
    
    # Rename regions to something presentable
    q$region <- plyr::mapvalues(q$region, from = c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"), to = c("5' UTR", "first\nexon", "first\nintron", "exon", "intron", "last\nexon", "last\nintron", "3' UTR"))
    
    # Reorder (complex first)
    q$title <- factor(q$title, levels=c(rbp1_rbp2_title, rbp1_only_title, rbp2_only_title))
    
    # # Plot v0 (3 histogram rows)
    ggplot(q, aes(x=pos)) + geom_histogram(binwidth=0.125) + coord_cartesian(xlim=c(0, 1)) + facet_grid(rows=vars(title), cols=vars(region), scales="free_y") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("Peak count")
    # # p <- ggplot(q, aes(x=pos)) + geom_histogram(binwidth=0.125) + coord_cartesian(xlim=c(0, 1)) + facet_grid(rows=vars(title), cols=vars(region), scales="free_y") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("Peak count")
    # # get_hist <- function(p) {
    # #   d <- ggplot_build(p)$data[[1]]
    # #   data.frame(x = d$x, xmin = d$xmin, xmax = d$xmax, y = d$y)
    # # }
    # # get_hist(p)
    ggsave(paste0("output-plot-metagene-v0-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    
    # v10: Density plot (geom_density)
    
    # Make "pos" a continuous value across the entire metagene ([0..1])
    q$globalpos <- NA
    for (i in 1:length(levels(q$region))) {
      # print(paste0(i,": ",levels(q$region)[i]))
      # q[q$region==levels(q$region)[i],]$globalpos <- q[q$region==levels(q$region)[i],]$pos + (i-1)
      # Scale original pos to between 0.1 .. 0.9 to leave a blank area between regions:
      q[q$region==levels(q$region)[i],]$globalpos <- ((q[q$region==levels(q$region)[i],]$pos * 8/10) + 0.1) + (i-1)
    }
    
    summary(q$globalpos)
    q$globalpos <- q$globalpos / length(levels(q$region))
    summary(q$globalpos)
    

    # Density segmented into regions
    ggplot(q, aes(x=pos)) + geom_density() + coord_cartesian(xlim=c(0, 1)) + facet_grid(cols=vars(region), rows=vars(title), scales="free_y") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("Peak count")
    ggplot(q, aes(x=pos)) + geom_density() + coord_cartesian(xlim=c(0, 1)) + facet_grid(cols=vars(region), rows=vars(title)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("Peak count")
    # Reorder (complex last)
    q$title <- factor(q$title, levels=c(rbp1_only_title, rbp2_only_title, rbp1_rbp2_title))
    ggplot(q, aes(x=pos, colour=title, fill=title)) + 
      geom_density() + 
      coord_cartesian(xlim=c(0, 1)) + 
      facet_grid(cols=vars(region)) + 
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      ylab("Peak count")
    ggplot(q, aes(x=pos, colour=title, fill=title)) + 
      geom_density(alpha=0.3) + 
      coord_cartesian(xlim=c(0, 1)) + 
      facet_grid(cols=vars(region)) + 
      scale_colour_manual(values=c("#8ea106", "#00a1d9", "#47d9bf")) + 
      scale_fill_manual(values=c("#ffffff", "#ffffff", "#47d9bfff")) + 
      theme_minimal() + 
      theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + 
      theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      xlab("mRNA region") + 
      ylab("Probability density") + 
      guides(colour=F, fill=F) + 
      ggtitle(paste0(rbp1," & ",rbp2)) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    # Density using "globalpos" (one continuous x-axis)
    # ggplot(q, aes(x=globalpos)) + geom_density() + coord_cartesian(xlim=c(0, 1)) + facet_grid(rows=vars(title), scales="free_y") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("Peak count")
    # ggplot(q, aes(x=globalpos, colour=title, fill=title)) + geom_density(alpha=0.3) + coord_cartesian(xlim=c(0, 1)) + scale_colour_manual(aesthetics = c("colour", "fill"), values=c("#47d9bf", "#8ea106", "#00a1d9")) + facet_grid(rows=vars(title), scales="free_y") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("Probability density")
    # ggsave(paste0("output-plot-metagene-v10-density-separate-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    # ggplot(q, aes(x=globalpos, colour=title, fill=title)) + geom_density(alpha=0.3) + coord_cartesian(xlim=c(0, 1)) + scale_colour_manual(aesthetics = c("colour", "fill"), values=c("#47d9bf", "#8ea106", "#00a1d9")) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab("Probability density")
    # # ggsave(paste0("output-plot-metagene-v10-density-merged-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    # Reorder (complex last)
    q$title <- factor(q$title, levels=c(rbp1_only_title, rbp2_only_title, rbp1_rbp2_title))
    # ggplot(q, aes(x=globalpos, colour=title, fill=title)) + geom_density(alpha=0.3) + coord_cartesian(xlim=c(0, 1)) + scale_colour_manual(aesthetics = c("colour", "fill"), values=c("#8ea106", "#00a1d9", "#47d9bf")) + facet_grid(rows=vars(title), scales="free_y") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + guides(colour=F, fill=F) + ylab("Probability density")
    # ggsave(paste0("output-plot-metagene-v10-density-separate-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")
    # ggplot(q, aes(x=globalpos, colour=title, fill=title)) + geom_density(alpha=0.3) + coord_cartesian(xlim=c(0, 1)) + scale_colour_manual(values=c("#8ea106", "#00a1d9", "#47d9bf")) + scale_fill_manual(values=c("#ffffff", "#ffffff", "#47d9bfff")) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + ylab("Probability density") + guides(colour=F, fill=F)
    # ggplot(q, aes(x=globalpos, colour=title, fill=title)) + geom_density(alpha=0.3) + geom_vline(xintercept=c(1:7/8), colour="#ffffff", size=2) + scale_x_continuous(breaks=1:8/8-1/16, labels=levels(q$region)) + coord_cartesian(xlim=c(0, 1)) + scale_colour_manual(values=c("#8ea106", "#00a1d9", "#47d9bf")) + scale_fill_manual(values=c("#ffffff", "#ffffff", "#47d9bfff")) + theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + xlab("mRNA region") + ylab("Probability density") + guides(colour=F, fill=F)
    # ggplot(q, aes(x=globalpos, colour=title, fill=title)) + geom_density(alpha=0.3) + geom_vline(xintercept=c(1:7/8), linetype="dotted") + scale_x_continuous(breaks=1:8/8-1/16, labels=levels(q$region)) + coord_cartesian(xlim=c(0, 1)) + scale_colour_manual(values=c("#8ea106", "#00a1d9", "#47d9bf")) + scale_fill_manual(values=c("#ffffff", "#ffffff", "#47d9bfff")) + theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("mRNA region") + ylab("Probability density") + guides(colour=F, fill=F)
    # ggplot(q, aes(x=globalpos, colour=title, fill=title)) + geom_density(alpha=0.3) + geom_vline(xintercept=c(1:7/8), colour="#ffffff", size=1.5) + scale_x_continuous(breaks=1:8/8-1/16, labels=levels(q$region)) + coord_cartesian(xlim=c(0, 1)) + scale_colour_manual(values=c("#8ea106", "#00a1d9", "#47d9bf")) + scale_fill_manual(values=c("#ffffff", "#ffffff", "#47d9bfff")) + theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("mRNA region") + ylab("Probability density") + guides(colour=F, fill=F)
    # ggplot(q, aes(x=globalpos, colour=title, fill=title)) + geom_density(alpha=0.3) + geom_vline(xintercept=c(1:7/8), linetype="dashed") + scale_x_continuous(breaks=1:8/8-1/16, labels=levels(q$region)) + coord_cartesian(xlim=c(0, 1)) + scale_colour_manual(values=c("#8ea106", "#00a1d9", "#47d9bf")) + scale_fill_manual(values=c("#ffffff", "#ffffff", "#47d9bfff")) + theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("mRNA region") + ylab("Probability density") + guides(colour=F, fill=F)
    summary(q$globalpos)
    ggplot(q, aes(x=globalpos, colour=title, fill=title)) + 
      geom_density(alpha=0.3, bw=1/160) + 
      geom_vline(xintercept=c(1:7/8), linetype="dashed", colour="#cccccc") + 
      scale_x_continuous(breaks=1:8/8-1/16, labels=levels(q$region)) + 
      coord_cartesian(xlim=c(0, 1)) + 
      scale_colour_manual(values=c("#8ea106", "#00a1d9", "#47d9bf")) + 
      scale_fill_manual(values=c("#ffffff", "#ffffff", "#47d9bfff")) + 
      theme_minimal() + 
      theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + 
      theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      xlab("mRNA region") + 
      ylab("Probability density") + 
      guides(colour=F, fill=F) + 
      ggtitle(paste0(rbp1," & ",rbp2)) + 
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0("output-plot-metagene-v10-density-combined-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")
    
    
    
    
    
    # Resampling - Originally for final_heatmap_3 / v11 below (Bar chart with error bars)
    # Reorder (complex middle)
    q$title <- factor(q$title, levels=c(rbp1_only_title, rbp1_rbp2_title, rbp2_only_title))
    # geom_bar aes(width) doesn't seem to exist anymore
    # q$width <- 0.5
    # q[q$title==rbp1_rbp2_title,]$width <- 1
    head(q)
    summary(q)
    length(unique(q$ensgv))
    
    # Subsample genes
    # ?summarise
    summary(q)
    head(q)
    str(q)
    q
    q %>% group_by(title) %>% tally()
    q %>% group_by(title, region) %>% tally()
    q %>% group_by(region, title, ensgv) %>% tally()
    
    q %>% group_by(title, region, .drop=FALSE) %>% summarise(n = n()) %>% mutate(freq = n / sum(n)) 
    # q %>% group_by(title, region, .drop=FALSE) %>% tally() %>% mutate(freq = n / sum(n)) # identical
    
    # qt <- q %>% group_by(region, title, ensgv) %>% tally(wt=region)
    
    # qt$count <- qt$n
    # qt$n <- NULL
        
    # qt %>% summarise(
    #   mean = mean(count),
    #   sd = sd(count),
    #   n = n(),
    #   se = sd / n
    # ) 
    
    # sample
    # mean(replicate(100, sd(sample(qt, replace=T))/sqrt(length(qt))))
    # qt
    
    # q %>% group_by(title, region, .drop=FALSE) %>% summarise(n = n()) %>% mutate(freq = n / sum(n)) 
    all_ensgvs <- unique(q$ensgv)
    tib <- tibble()
    set.seed(1)
    for (i in 1:num_resamples) {
  
      # # Subsample to 50% (with replacement)
      # these_ensgvs <- sample(all_ensgvs, floor(length(all_ensgvs) / 2), replace=T)
      # Subsample 100% (with replacement)
      these_ensgvs <- sample(all_ensgvs, length(all_ensgvs), replace=T)
      
      this_tib <- subset(q, ensgv %in% these_ensgvs) %>% group_by(title, region, .drop=FALSE) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
      
      this_tib$title <- as.character(this_tib$title)
      this_tib$rbp1 <- rbp1
      this_tib$rbp2 <- rbp2
      this_tib$resample <- i
      
      # Assign simpler titles
      this_tib[this_tib$title == rbp1_only_title,]$title <- "rbp1_only"
      this_tib[this_tib$title == rbp2_only_title,]$title <- "rbp2_only"
      this_tib[this_tib$title == rbp1_rbp2_title,]$title <- "rbp1_rbp2"
      # this_tib
      
      # Append to tib
      if (nrow(tib) > 0) {
        tib <- rbind(tib, this_tib)
      } else {
        tib <- this_tib
      }

      # Append to final_heatmap_3
      # if (nrow(final_heatmap_3) > 0) {
      #   final_heatmap_3 <- rbind(final_heatmap_3, this_tib)
      # } else {
      #   final_heatmap_3 <- this_tib
      # }
    }
    # final_heatmap_3
    # summary(final_heatmap_3)
    tib
    summary(tib)
    
    tib$pair <- paste0(tib$rbp1,"|",tib$rbp2)
    tib %>% group_by(pair, region, title, .drop=F) %>% summarise(mean = mean(freq), n = n(), sd = sd(freq), se = sd / n)
    tib2 <- tib %>% group_by(pair, region, title, .drop=F) %>% summarise(mean = mean(freq), n = n(), sd = sd(freq), se = sd / n)
    tib2
    ggplot(tib2, aes(x=region, y=mean, fill=title)) + geom_col(position="dodge") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position=position_dodge(width=0.9), width=0.25) + scale_fill_manual(values=c("#8ea106", "#47d9bf", "#00a1d9")) + theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("mRNA region") + ylab("Probability density") + guides(colour=F, fill=F) + ggtitle(paste0(rbp1," & ",rbp2)) + theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0("output-plot-metagene-v11-bars-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5*2, height=91.5, units="mm")
    # Write data
    tib2$region <- plyr::mapvalues(tib2$region, from = c("5' UTR", "first\nexon", "first\nintron", "exon", "intron", "last\nexon", "last\nintron", "3' UTR"), to = c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"))
    write.table(tib2, paste0("output-plot-metagene-v11-txt-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-1bin.txt"), sep="\t", quote=F, row.names=F)
    
    
    
    
    # v12: Resampling-based confidence intervals on a density plot
    # https://stackoverflow.com/questions/45908120/add-shaded-standard-error-curves-to-geom-density-in-ggplot2
    
    # Generate samples
    
    # tmpi <- iris %>% group_by(Species) %>% sample_frac(size = 1, replace = T)
    # str(tmpi)
    # length(unique(tmpi$Species))
    
    # tibdf <- subset(q, ensgv %in% these_ensgvs)
    # str(tibdf)
    
    # # Assign simpler titles
    # tibdf$title <- as.character(tibdf$title)
    # tibdf[tibdf$title == rbp1_only_title,]$title <- "rbp1_only"
    # tibdf[tibdf$title == rbp2_only_title,]$title <- "rbp2_only"
    # tibdf[tibdf$title == rbp1_rbp2_title,]$title <- "rbp1_rbp2"
    # tibdf
    
    # # DEBUG: Filter for the complex only
    # tibdf <- subset(tibdf, title=="rbp1_rbp2")
    
    all_ensgvs <- unique(q$ensgv)
    tib <- tibble()
    set.seed(1)
    for (i in 1:num_resamples) {
      
      # # Subsample to 50% (with replacement)
      # these_ensgvs <- sample(all_ensgvs, floor(length(all_ensgvs) / 2), replace=T)
      # Subsample 100% (with replacement)
      these_ensgvs <- sample(all_ensgvs, length(all_ensgvs), replace=T)
      
      # this_tib <- subset(tibdf, ensgv %in% these_ensgvs) %>% group_by(title, region, .drop=FALSE) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
      # this_tib <- subset(tibdf, ensgv %in% these_ensgvs) %>% group_by(title, region, .drop=FALSE) %>% select(region, title, pos)
      this_tib <- subset(q, ensgv %in% these_ensgvs) %>% group_by(title, .drop=FALSE) %>% select(title, globalpos)
      this_tib$title <- as.character(this_tib$title)
      
      # #DEBUG
      # this_tib <- subset(tibdf, ensgv %in% these_ensgvs) %>% group_by(region, .drop=FALSE) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
      # #END DEBUG
      
      this_tib$rbp1 <- rbp1
      this_tib$rbp2 <- rbp2
      # this_tib$resample <- i
      this_tib
      
      # Append to tib
      if (nrow(tib) > 0) {
        tib <- rbind(tib, tibble(bs=i, data=list(this_tib)))
      } else {
        tib <- tibble(bs=i, data=list(this_tib))
      }
    }
    tib
    tib[1,2] %>% unnest(cols=data)
    # summary(tib)
    
    
    # Compute within-sample density
    densities.within <-
    tib %>%
      unnest(cols=data) %>% 
      # group_by(bs, region) %>%
      # group_by(bs, title, region) %>%
      group_by(bs, title) %>%
      # do(tidy(density(.$pos, 
      do(tidy(density(.$globalpos, 
                      from = 0,
                      to = 1,
                      # n = 128,
                      # bw=1/80
                      # bw=1/160
                      bw=1/320
      )))
    densities.within
    
    # Summarize densities into quantiles
    densities.qtiles <-
      densities.within %>%
      # rename(pos = x, dens = y) %>%
      rename(globalpos = x, dens = y) %>%
      ungroup() %>%
      # group_by(title, region, pos) %>%
      # group_by(title, region, globalpos) %>%
      group_by(title, globalpos) %>%
      # group_by(region, pos) %>% 
      summarise(q05 = quantile(dens, 0.025),
                q50 = quantile(dens, 0.5),
                q95 = quantile(dens, 0.975)) 
    densities.qtiles

    # Visualize and compare
    
    # qtib <- q %>% group_by(title, region, .drop=FALSE) %>% summarise(n = n()) %>% mutate(pos = n / sum(n))
    # qtib <- q %>% group_by(title, region, .drop=FALSE) %>% summarise(n = n()) %>% mutate(globalpos = n / sum(n))
    # qtib <- tibdf %>% group_by(title, region, .drop=FALSE)
    # qtib <- tibdf %>% group_by(title, .drop=FALSE)
    # # Assign simpler titles
    # qtib$title <- as.character(qtib$title)
    # qtib[qtib$title == rbp1_only_title,]$title <- "rbp1_only"
    # qtib[qtib$title == rbp2_only_title,]$title <- "rbp2_only"
    # qtib[qtib$title == rbp1_rbp2_title,]$title <- "rbp1_rbp2"

    # # Reorder (complex last)
    # densities.qtiles$title <- factor(densities.qtiles$title, levels=c(rbp1_only_title, rbp2_only_title, rbp1_rbp2_title), ordered=T)
    # # qtib$title <- factor(qtib$title, levels=c(rbp1_only_title, rbp1_rbp2_title, rbp2_only_title), ordered=T)
    # # q$title <- factor(q$title, levels=c(rbp1_only_title, rbp2_only_title, rbp1_rbp2_title))

    # Reorder (complex middle)
    densities.qtiles$title <- factor(densities.qtiles$title, levels=c(rbp1_only_title, rbp1_rbp2_title, rbp2_only_title), ordered=T)
    
    # #DEBUG
    # qtib <- q %>% group_by(region, .drop=FALSE) %>% summarise(n = n()) %>% mutate(pos = n / sum(n))
    # #END DEBUG
    
    # ggplot(densities.qtiles, aes(x=pos, y=q50)) +
    # ggplot(densities.qtiles, aes(x=pos, y=q50, colour=title, fill=title)) +
    p <- ggplot(densities.qtiles, aes(x=globalpos, y=q50, colour=title, fill=title)) +
      # facet_grid(cols=vars(region)) +
      # facet_grid(rows=vars(title)) +
      geom_ribbon(aes(ymin = q05, ymax = q95), alpha=0.3, colour=NA) +
      # stat_density(data = qtib,
      #              # aes(pos, ..density.., color = "raw density"),
      #              aes(globalpos, ..density.., color = "raw density"),
      #              size = 2, geom = "line") +
      geom_line() +
      # scale_color_manual(values = c("red", "black")) +
      coord_cartesian(xlim=c(0, 1)) + 
      scale_colour_manual(NULL, aesthetics=c("colour", "fill"), values=c("#8ea106", "#47d9bf", "#00a1d9")) + 
      geom_vline(xintercept=c(1:7/8), linetype="dashed", colour="#cccccc") + 
      scale_x_continuous(breaks=1:8/8-1/16, labels=levels(q$region)) + 
      labs(y = "density") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + 
      theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab("mRNA region") + 
      ylab("Probability density") + 
      # guides(colour=F, fill=F) + 
      # ggtitle(paste0(rbp1," & ",rbp2)) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position="top")
      # theme(legend.position="right")
    # Make legend smaller
    p <- p + guides(color=guide_legend(override.aes=list(size=0.5)))
    # p <- p + guides(fill=guide_legend(override.aes=list(size=0.5)))
    p <- p + theme(legend.title = element_text(size = 7), legend.text = element_text(size = 7))
    # Save wide version
    ggsave(paste0("output-plot-metagene-v12-wide-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5*2, height=91.5/2, units="mm")
    # Save square version
    # p <- p + theme(legend.position="right")
    ggsave(paste0("output-plot-metagene-v12-square-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5, height=91.5*3/4, units="mm")
    # Split
    p + facet_grid(rows=vars(title), scales="free_y") + guides(colour=F, fill=F) + theme(strip.text.y = element_text(hjust = 0))
    split_height <- 91.5 * 3/4 + (3 - 2) * 26.9434 # actual plot area
    ggsave(paste0("output-plot-metagene-v12-split-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5*2, height=split_height, units="mm")
    

    
    # Convert title to character
    densities.qtiles$title <- as.character(densities.qtiles$title)
    
    # Store densities.qtiles for later
    # alldensities[[paste0(rbp1,"|",rbp2)]] <- densities.qtiles
    alldensities %<>% add_row(rbp1=rbp1, rbp2=rbp2, densities=list(densities.qtiles))
    
    
    
    
    # Prepare cyan/magenta heatmap plot (condensing the 3 histogram rows into one heatmap row)
    
    # Assign histogram bins
    # q$bin <- factor(cut(q$pos, breaks = 0:8/8, include.lowest = T, labels=paste0("<",1:8/8)), ordered=T)
    # q$bin <- factor(cut(q$pos, breaks = 0:8/8, include.lowest = T, labels=1:8/8), ordered=T)
    q$bin <- factor(cut(q$pos, breaks = 0:bincount/bincount, include.lowest = T, labels=1:bincount/bincount), ordered=T)
    # # Same look as the ggplot2 histogram (9 bins, actually) (it's always 0.125 centered on something, beginning with 0 and ending at 1. This meant the bins didn't sample the same amount of space, actually, so I'm not using this)
    # q$bin <- factor(cut(q$pos, breaks = (0:8/8)-1/16, include.lowest = T), ordered=T)
    levels(q$bin)
    str(q$bin)
    
    # Plot v1 (to test out the histogram bins)
    ggplot(q, aes(x=bin)) + geom_bar(width=1) + facet_grid(rows=vars(title), cols=vars(region), scales="free_y") + ylab("Peak count") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v1-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=183, height=91.5, units="mm")
    
    
    
    
    
    
    all_ensgvs <- unique(q$ensgv)
    tib <- tibble()
    set.seed(1)
    for (i in 1:num_resamples) {
      
      # # Subsample to 50% (with replacement)
      # these_ensgvs <- sample(all_ensgvs, floor(length(all_ensgvs) / 2), replace=T)
      # Subsample 100% (with replacement)
      these_ensgvs <- sample(all_ensgvs, length(all_ensgvs), replace=T)
      
      this_tib <- subset(q, ensgv %in% these_ensgvs) %>% group_by(title, region, bin, .drop=FALSE) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
      
      this_tib$title <- as.character(this_tib$title)
      this_tib$rbp1 <- rbp1
      this_tib$rbp2 <- rbp2
      this_tib$resample <- i
      
      # Assign simpler titles
      this_tib[this_tib$title == rbp1_only_title,]$title <- "rbp1_only"
      this_tib[this_tib$title == rbp2_only_title,]$title <- "rbp2_only"
      this_tib[this_tib$title == rbp1_rbp2_title,]$title <- "rbp1_rbp2"
      # this_tib
      
      # Append to tib
      if (nrow(tib) > 0) {
        tib <- rbind(tib, this_tib)
      } else {
        tib <- this_tib
      }
      
      # Append to final_heatmap_3
      # if (nrow(final_heatmap_3) > 0) {
      #   final_heatmap_3 <- rbind(final_heatmap_3, this_tib)
      # } else {
      #   final_heatmap_3 <- this_tib
      # }
    }
    # final_heatmap_3
    # summary(final_heatmap_3)
    tib
    summary(tib)
    
    # As boxplots:
    tib
    ggplot(tib, aes(x=bin, y=freq, fill=title, colour=title)) + geom_boxplot(position="identity", width=0.9, alpha=0.3, outlier.shape=NA) + facet_grid(cols=vars(region)) + scale_fill_manual(aesthetics=c("fill", "colour"), values=c("#8ea106", "#47d9bf", "#00a1d9")) + theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("mRNA region") + ylab("Probability density") + guides(colour=F, fill=F) + ggtitle(paste0(rbp1," & ",rbp2)) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v11-box-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5*2, height=91.5, units="mm")
    ggplot(tib, aes(x=bin, y=freq, fill=title, colour=title)) + geom_boxplot(position="dodge", width=0.9, alpha=0.3, outlier.shape=NA) + facet_grid(cols=vars(region)) + scale_fill_manual(aesthetics=c("fill", "colour"), values=c("#8ea106", "#47d9bf", "#00a1d9")) + theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("mRNA region") + ylab("Probability density") + guides(colour=F, fill=F) + ggtitle(paste0(rbp1," & ",rbp2)) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v11-boxdodged-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5*2, height=91.5, units="mm")
    # Write data
    tib$region <- plyr::mapvalues(tib$region, from = c("5' UTR", "first\nexon", "first\nintron", "exon", "intron", "last\nexon", "last\nintron", "3' UTR"), to = c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"))
    write.table(tib, paste0("output-plot-metagene-v11-txt-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-8bins.txt"), sep="\t", quote=F, row.names=F)
    
    
    
    
    
    
    
    # Cycle through bins
    q2 <- data.frame("title" = factor(levels = levels(q$title)), "region" = factor(levels = levels(q$region)), "bin" = factor(levels = levels(q$bin)), "value" = numeric())
    # l <- list()
    # Get counts in a separate dataframe
    for (my_title in levels(q$title)) {
      for (my_region in levels(q$region)) {
        for (my_bin in levels(q$bin)) {
          # print(paste(my_title, my_region, my_bin))
          # str(data.frame(title = my_title, region = my_region, bin = my_bin))
          q2 <- rbind(q2, data.frame(title = my_title, region = my_region, bin = my_bin, value = nrow(q[q$title == my_title & q$region == my_region & q$bin == my_bin,])))
          # +1
          # q2 <- rbind(q2, data.frame(title = my_title, region = my_region, bin = my_bin, value = nrow(q[q$title == my_title & q$region == my_region & q$bin == my_bin,]) + 1))
          # q[q$title == my_title & q$region == my_region & q$bin == my_bin,]$cyan <-
        }
      }
    }
    q2
    # Test plot: Looks exactly the same as v2 (the histograms)
    ggplot(q2, aes(x=bin, y=value)) + geom_bar(width=1, stat="identity") + facet_grid(rows=vars(title), cols=vars(region), scales="free_y") + ylab("Peak count") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    q3 <- q2
    # Scale to between 0 and 1
    for (my_title in levels(q$title)) {
      for (my_region in levels(q$region)) {
        for (my_bin in levels(q$bin)) {
          # q2[q2$title == my_title & q2$region == my_region & q2$bin == my_bin,]$value
          # max(q2[q2$title == my_title,]$value)
          q3[q3$title == my_title & q3$region == my_region & q3$bin == my_bin,]$value <- q2[q2$title == my_title & q2$region == my_region & q2$bin == my_bin,]$value / max(q2[q2$title == my_title,]$value)
        }
      }
    }
    ggplot(q3, aes(x=bin, y=value)) + geom_bar(width=1, stat="identity") + facet_grid(rows=vars(title), cols=vars(region), scales="free_y") + ylab("Peak count") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    # Write table for Sebastian
    head(q3)
    q3tmp <- q3
    q3tmp$title <- as.character(q3tmp$title)
    q3tmp[q3tmp$title == rbp1_rbp2_title,]$title <- "pair"
    q3tmp[q3tmp$title == rbp1_only_title,]$title <- paste0(rbp1, "_only")
    q3tmp[q3tmp$title == rbp2_only_title,]$title <- paste0(rbp2, "_only")
    q3tmp
    q3tmp$region <- plyr::mapvalues(q3tmp$region, from = c("5' UTR", "first\nexon", "first\nintron", "exon", "intron", "last\nexon", "last\nintron", "3' UTR"), to = c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"))
    write.table(q3tmp, paste0("output-plot-metagene-v1-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".tsv"), sep="\t", quote=F, row.names=F, col.names=F)
    
    
    # Get cyan (RBP1) and magenta (RBP2) contribution values
    head(q3)
    
    # 1 - absolute distance
    head(q3)
    q4 <- data.frame("region" = factor(levels = levels(q$region)), "bin" = factor(levels = levels(q$bin)), "rbp1" = numeric(), "rbp2" = numeric(), "value" = numeric())
    for (my_region in levels(q$region)) {
      for (my_bin in levels(q$bin)) {
        tmp_rbp1 <- 1 - abs(q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value - q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value)
        tmp_rbp2 <- 1 - abs(q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value - q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value)
        tmp_value <- NA
        q4 <- rbind(q4, data.frame(region = my_region, bin = my_bin, rbp1 = tmp_rbp1, rbp2 = tmp_rbp2, value = tmp_value))
      }
    }
    ggplot(q4, aes(x=bin, y=rbp1)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v2-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-1.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    ggplot(q4, aes(x=bin, y=rbp2)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v2-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-2.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    
    # absolute distance
    head(q3)
    q4 <- data.frame("region" = factor(levels = levels(q$region)), "bin" = factor(levels = levels(q$bin)), "rbp1" = numeric(), "rbp2" = numeric(), "value" = numeric())
    for (my_region in levels(q$region)) {
      for (my_bin in levels(q$bin)) {
        tmp_rbp1 <- abs(q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value - q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value)
        tmp_rbp2 <- abs(q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value - q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value)
        tmp_value <- NA
        q4 <- rbind(q4, data.frame(region = my_region, bin = my_bin, rbp1 = tmp_rbp1, rbp2 = tmp_rbp2, value = tmp_value))
      }
    }
    ggplot(q4, aes(x=bin, y=rbp1)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v3-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-1.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    ggplot(q4, aes(x=bin, y=rbp2)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v3-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-2.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    
    # Divide Both by RBPn
    head(q3)
    q4 <- data.frame("region" = factor(levels = levels(q$region)), "bin" = factor(levels = levels(q$bin)), "rbp1" = numeric(), "rbp2" = numeric(), "value" = numeric())
    for (my_region in levels(q$region)) {
      for (my_bin in levels(q$bin)) {
        tmp_rbp1 <- q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value / q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value
        tmp_rbp2 <- q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value / q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value
        tmp_value <- NA
        q4 <- rbind(q4, data.frame(region = my_region, bin = my_bin, rbp1 = tmp_rbp1, rbp2 = tmp_rbp2, value = tmp_value))
      }
    }
    ggplot(q4, aes(x=bin, y=rbp1)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v4-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-1.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    ggplot(q4, aes(x=bin, y=rbp2)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v4-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-2.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    
    # # Divide RBPn by Both
    head(q3)
    q4 <- data.frame("region" = factor(levels = levels(q$region)), "bin" = factor(levels = levels(q$bin)), "rbp1" = numeric(), "rbp2" = numeric(), "value" = numeric())
    for (my_region in levels(q$region)) {
      for (my_bin in levels(q$bin)) {
        tmp_rbp1 <- q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value / q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value
        tmp_rbp2 <- q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value / q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value
        tmp_value <- NA
        q4 <- rbind(q4, data.frame(region = my_region, bin = my_bin, rbp1 = tmp_rbp1, rbp2 = tmp_rbp2, value = tmp_value))
      }
    }
    ggplot(q4, aes(x=bin, y=rbp1)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v5-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-1.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    ggplot(q4, aes(x=bin, y=rbp2)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0("output-plot-metagene-v5-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,"-2.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
    
    # Final
    q4 <- data.frame("region" = factor(levels = levels(q$region)), "bin" = factor(levels = levels(q$bin)), "rbp1" = numeric(), "rbp2" = numeric(), "value" = numeric())
    for (my_region in levels(q$region)) {
      for (my_bin in levels(q$bin)) {
        
        # distance
        tmp_rbp1 <- q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value - q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value
        tmp_rbp2 <- q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value - q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value
        
        ### MAIN FORMULA HERE
        # # distance vs. max(rbp1, rbp2)
        # tmp_value <- q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value
        # - max(q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value,
        #       q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value)
        # tmp_value <- q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value - max(q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value, q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value)
        # q4 <- rbind(q4, data.frame(region = my_region, bin = my_bin, rbp1 = tmp_rbp1, rbp2 = tmp_rbp2, value = tmp_value))
        # # max(abs(distance vs. rbp1), abs(distance vs. rbp2))
        # # tmp_value <- max(
        # #   (q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value
        # #    - q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value),
        # #   (q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value
        # #    - q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value))
        # tmp_value <- max(abs(q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value - q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value), 
        #                  abs(q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value - q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value))
        
        # Signed distance from the highest or lowest individual RBP, whichever is greater in absolute terms
        
        # Utility values
        tmp_value_pair <- q3[q3$title == rbp1_rbp2_title & q3$region == my_region & q3$bin == my_bin,]$value
        tmp_value_rbp1 <- q3[q3$title == rbp1_only_title & q3$region == my_region & q3$bin == my_bin,]$value
        tmp_value_rbp2 <- q3[q3$title == rbp2_only_title & q3$region == my_region & q3$bin == my_bin,]$value

        # Default to 0 (if the pair is in between rbp1 and rbp2)
        tmp_value <- 0

        # If the pair is above both individual RBPs...
        if ((tmp_value_pair > tmp_value_rbp1) & (tmp_value_pair > tmp_value_rbp2)) {
          # Use the distance from the higher individual RBP
          tmp_value <- min((tmp_value_pair - tmp_value_rbp1), (tmp_value_pair - tmp_value_rbp2))
        }
        # If the pair is below both individual RBPs...
        if ((tmp_value_pair < tmp_value_rbp1) & (tmp_value_pair < tmp_value_rbp2)) {
          # Use the distance from the lower individual RBP
          tmp_value <- max((tmp_value_pair - tmp_value_rbp1), (tmp_value_pair - tmp_value_rbp2))
        }
        ### END MAIN FORMULA
        
        # Get sample size (n)
        tmp_n <- q2[q2$title == rbp1_rbp2_title & q2$region == my_region & q2$bin == my_bin,]$value

        q4 <- rbind(q4, data.frame(region = my_region, bin = my_bin, rbp1 = tmp_rbp1, rbp2 = tmp_rbp2, value = tmp_value))
        
        final_heatmap <- rbind(final_heatmap, data.frame(rbp1=rbp1, rbp2=rbp2, pair=paste0(rbp1,"|",rbp2), region=my_region, bin=my_bin, value=tmp_value, raw_value=tmp_value_pair, n=tmp_n))
        
        # Second heatmap variant: Instead of using max(rbp1, rbp2), add two rows (one for rbp1, one for rbp2, compared to the complex)
        final_heatmap_2 <- rbind(final_heatmap_2, data.frame(rbp1=rbp1, rbp2=rbp2, pair=paste0(rbp1,"|",rbp2,": ", rbp1), region=my_region, bin=my_bin, value=tmp_rbp1))
        final_heatmap_2 <- rbind(final_heatmap_2, data.frame(rbp1=rbp1, rbp2=rbp2, pair=paste0(rbp1,"|",rbp2,": ", rbp2), region=my_region, bin=my_bin, value=tmp_rbp2))

        # Another heatmap variant: Simply keep the individual and pair values
        final_heatmap_4 <- rbind(final_heatmap_4, data.frame(rbp1=rbp1, rbp2=rbp2, pair=rbp1_only_title, region=my_region, bin=my_bin, value=tmp_value_rbp1))
        final_heatmap_4 <- rbind(final_heatmap_4, data.frame(rbp1=rbp1, rbp2=rbp2, pair=rbp2_only_title, region=my_region, bin=my_bin, value=tmp_value_rbp2))
        final_heatmap_4 <- rbind(final_heatmap_4, data.frame(rbp1=rbp1, rbp2=rbp2, pair=rbp1_rbp2_title, region=my_region, bin=my_bin, value=tmp_value_pair))
      }


      # Pair distance values (averaged across bins per region)
      tmp_distance <- mean(final_heatmap[final_heatmap$rbp1==rbp1 & final_heatmap$rbp2==rbp2 & final_heatmap$region==my_region,]$value)
      pair_distance <- rbind(pair_distance, data.frame(rbp1=rbp1, rbp2=rbp2, pair=paste0(rbp1,"|",rbp2), region=my_region, value=tmp_distance))
    }
    # Pair distance values (averaged across all regions)
    tmp_distance <- mean(final_heatmap[final_heatmap$rbp1==rbp1 & final_heatmap$rbp2==rbp2,]$value)
    pair_distance <- rbind(pair_distance, data.frame(rbp1=rbp1, rbp2=rbp2, pair=paste0(rbp1,"|",rbp2), region="all", value=tmp_distance))
    
    q4
    ggplot(q4, aes(x=bin, y=value)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Enrichment factor") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(strip.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + theme(plot.title = element_text(hjust=0.5)) + ggtitle(paste0(rbp1," & ",rbp2))
    ggsave(paste0("output-plot-metagene-v6-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")
    
    # # Assign colour
    # # Cyan: Fraction of Both vs. RBP1
    # # Magenta: Fraction of Both vs. RBP2
    # q4$col <- paste0("#",)
    # ggplot(q4, aes(x=bin, y=rbp1)) + geom_bar(width=1, stat="identity") + facet_grid(cols=vars(region), scales="free_y") + ylab("Peak count") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  }
}



# Plot final heatmap
final_heatmap_backup <- final_heatmap
# Beautify
final_heatmap$pair <- fct_rev(factor(paste0(final_heatmap$rbp1," & ",final_heatmap$rbp2), ordered=T))
# final_heatmap$region <- plyr::mapvalues(final_heatmap$region, from = c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"), to = c("5' UTR", "First\nexon", "First\nintron", "Exon", "Intron", "Last\nexon", "Last\nintron", "3' UTR"))



# # Write to file
# final_heatmap_tmp <- final_heatmap
# levels(final_heatmap_tmp$region)
# final_heatmap_tmp$region <- plyr::mapvalues(final_heatmap$region, from = c("5' UTR", "first\nexon", "first\nintron", "exon", "intron", "last\nexon", "last\nintron", "3' UTR"), to = c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"))
# write.table(final_heatmap_tmp, paste0("output-table-final_heatmap.tsv"), sep="\t", quote=F, row.names=F)





# Sort the heatmap (sorting it first by column 1, then 2, etc.)
# sortmap <- final_heatmap
# sortmap$raw_value <- NULL
# asortmap <- acast(sortmap, pair~region+bin)
# asortmap
# ord.mat = function(M, decr = F, cols = NULL){
#   if(is.null(cols))
#     cols = 1: ncol(M)
#   out = do.call( "order", as.data.frame(M[,cols]))
#   if (decr)
#     out = rev(out)
#   return(M[out,])
# }
# asortmap_sorted <- ord.mat(asortmap, decr=T)
# asortmap_sorted
# str(asortmap_sorted)
# head(asortmap_sorted)
# head(rownames(asortmap_sorted))
# str(rownames(asortmap_sorted))
# # asortmap_sorted_df <- as.data.frame(asortmap_sorted)
# # asortmap_sorted_df
# # str(asortmap_sorted_df)
# target_order <- rownames(asortmap_sorted)
# # head(final_heatmap[match(target_order, final_heatmap$pair),])
# # str(final_heatmap[match(target_order, final_heatmap$pair),])
# # # Reorder final_heatmap by the "clustering" (sorting)
# # # final_heatmap_tmp <- final_heatmap[match(target_order, final_heatmap$pair),]
# final_heatmap <- final_heatmap_backup
# final_heatmap$pair <- as.character(final_heatmap$pair)
# str(final_heatmap)
# str(target_order)
# final_heatmap <- left_join(data.frame(pair=target_order, stringsAsFactors=F), final_heatmap, by="pair")
# # final_heatmap$pair <- factor(final_heatmap$pair, ordered=T)
# final_heatmap$pair <- factor(final_heatmap$pair, levels=target_order)
# str(final_heatmap)


# Hierarchically cluster the heatmap
sortmap <- final_heatmap
sortmap$raw_value <- NULL
asortmap <- acast(sortmap, pair~region+bin, value.var="value")
asortmap
set.seed(1)
asortmap_hclust <- hclust(dist(asortmap), method="complete")
plot(asortmap_hclust, hang=-1)
pdf(paste0("output-plot-metagene-v7.0-",clip_cobinding,"-",type,"-",style,"-all.pdf"), width=5, height=4)
plot(asortmap_hclust, hang=-1)
dev.off()
target_order_hclust <- asortmap_hclust$labels[c(asortmap_hclust$order)]
target_order_hclust
target_order_hclust <- factor(target_order_hclust, levels=target_order_hclust)
target_order_hclust
# Apply target order
final_heatmap$pair <- as.character(final_heatmap$pair)
target_order_hclust <- as.character(target_order_hclust)
target_order_hclust
final_heatmap <- left_join(data.frame(pair=target_order_hclust, stringsAsFactors=F), final_heatmap, by="pair")
final_heatmap$pair <- factor(final_heatmap$pair, levels=target_order_hclust)
str(final_heatmap)


# Reverse order
final_heatmap$pair <- factor(final_heatmap$pair, levels=rev(levels(final_heatmap$pair)))
# ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + scale_fill_gradient2() + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
# ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + coord_cartesian(xlim=c(0, 1)) + scale_fill_gradient2() + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank()) + theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
# ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + coord_cartesian(xlim=c(0, 1)) + scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + scale_fill_gradient2() + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
# ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + coord_cartesian(xlim=c(0, 1)) + scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + scale_fill_gradient2() + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA))
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + scale_fill_gradient2("Enrichment") + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v7.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + scale_fill_gradient2("Enrichment") + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v7.2-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
# 7.3.1
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_gradient2("Enrichment") + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0.5))
ggsave(paste0("output-plot-metagene-v7.3.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*1.5, units="mm")
# 7.3.2
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + geom_text(aes(label=n)) + scale_fill_gradient2("Enrichment") + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0.5))
ggsave(paste0("output-plot-metagene-v7.3.2-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*1.5, units="mm")
# 7.3.3 (colour-capped)
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_gradient2("Enrichment", limits=c(-0.2, 0.2), oob=scales::squish) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0.5))
ggsave(paste0("output-plot-metagene-v7.3.3-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*1.5, units="mm")
# 7.3.4 (manual multi-colour gradient)
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() +
  # scale_fill_gradientn("Enrichment of\nthe complex vs.\nindividual RBPs", breaks=c(-0.2, 0, 0.2, 0.4), colours=c("blue", "white", "orange", "red")) +
  # scale_fill_gradientn("Enrichment of\nthe complex vs.\nindividual RBPs", breaks=c(min(final_heatmap$value), -0.2, 0, 0.2, 0.4, max(final_heatmap$value)), colours=c("blue", "blue", "white", "orange", "red", "red")) +
  # scale_fill_gradientn("Enrichment of\nthe complex vs.\nindividual RBPs", breaks=c(-0.2, 0, 0.2, 0.4), colours=c("blue", "white", "orange", "red"), limits=c(-0.2, 0.4), oob=scales::squish) +
  # scale_fill_gradientn("Enrichment of\nthe complex vs.\nindividual RBPs", values=scales::rescale(c(min(final_heatmap$value), -0.2, 0, 0.2, 0.4, max(final_heatmap$value))), colours=c("blue", "blue", "white", "orange", "red", "red")) +
  scale_fill_gradientn("Enrichment of\nthe complex vs.\nindividual RBPs", values=scales::rescale(c(min(final_heatmap$value), -0.2, 0, 0.2, 0.4, max(final_heatmap$value))), colours=c("red", "red", "white", "blue", "darkblue", "darkblue")) +
  facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0.5))
ggsave(paste0("output-plot-metagene-v7.3.4-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*1.5, units="mm")
# 7.3.5 (manual multi-colour gradient based on brewer.pal(n=3, name="RdBu") and n=5 and n=7)
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() +
  scale_fill_gradientn("Enrichment of\nthe complex vs.\nindividual RBPs", values=scales::rescale(c(min(final_heatmap$value), -0.2, 0, 0.2, 0.4, max(final_heatmap$value))),
                       # colours=brewer.pal(n=7, name="RdBu")[2:7]) +
                       # colours=c("#B2182B", "#EF8A62", "#FFFFFF", "#67A9CF", "#2166AC", "#053061")) +
                       # colours=c("#EF8A62", "#EF8A62", "#FFFFFF", "#67A9CF", "#0571B0", "#0571B0")) +
                       colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
                        # colours=c("#CA0020", "#CA0020", "#FFFFFF", "#0571B0", "#053061", "#053061")) +
                        # colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
                        # colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
                        # colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0.5)) +
  theme(legend.title = element_text(size = 8.8), legend.text = element_text(size = 8.8)) +
  # guides(fill = guide_legend(override.aes = list(size = 0.5)))
  theme(legend.key.size = unit(5, "mm"))
ggsave(paste0("output-plot-metagene-v7.3.5-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*1.5, units="mm")
# 7.4
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v7.4.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", option="A", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v7.4.2-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", option="B", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v7.4.3-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", option="C", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v7.4.4-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
# # Viridis ranges from https://www.thinkingondata.com/wp-content/uploads/2018/06/sample-palette.png:
# ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_gradient2("Enrichment", low="#440154", mid="#ffffff", high="#fde725", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
# ggsave(paste0("output-plot-metagene-v7.4.5-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
# ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_gradient2("Enrichment", low="#482677", mid="#ffffff", high="#b8de29", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
# ggsave(paste0("output-plot-metagene-v7.4.6-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
# ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_gradient2("Enrichment", low="#39568c", mid="#ffffff", high="#55c667", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
# ggsave(paste0("output-plot-metagene-v7.4.7-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
# ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_distiller("Enrichment", palette="RdBu", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_distiller("Enrichment", palette="RdBu", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v7.4.8-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")



# Raw complex value version: Recluster

# Hierarchically cluster the heatmap
sortmap <- final_heatmap
sortmap$value <- NULL
asortmap <- acast(sortmap, pair~region+bin, value.var="raw_value")
asortmap
set.seed(1)
asortmap_hclust <- hclust(dist(asortmap), method="complete")
plot(asortmap_hclust, hang=-1)
pdf(paste0("output-plot-metagene-v7.5.0-",clip_cobinding,"-",type,"-",style,"-all.pdf"), width=5, height=4)
plot(asortmap_hclust, hang=-1)
dev.off()
target_order_hclust <- asortmap_hclust$labels[c(asortmap_hclust$order)]
target_order_hclust
target_order_hclust <- factor(target_order_hclust, levels=target_order_hclust)
target_order_hclust
# Apply target order
final_heatmap$pair <- as.character(final_heatmap$pair)
target_order_hclust <- as.character(target_order_hclust)
target_order_hclust
final_heatmap <- left_join(data.frame(pair=target_order_hclust, stringsAsFactors=F), final_heatmap, by="pair")
final_heatmap$pair <- factor(final_heatmap$pair, levels=target_order_hclust)
str(final_heatmap)

# 7.5
ggplot(final_heatmap, aes(x=as.numeric(as.character(bin)), y=pair, fill=raw_value)) + geom_tile() + scale_fill_gradient2("Binding") + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v7.5.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*2, units="mm")






# Plot second final heatmap variant: Instead of using max(rbp1, rbp2), plot two rows (one for rbp1, one for rbp2, compared to the complex)
# Beautify
levels(final_heatmap_2$pair) <- sub("^(\\w+)\\|(\\w+): (\\w+)$", "\\1::\\2 vs. \\3 only", levels(final_heatmap_2$pair), perl=T)
# final_heatmap_2$region <- plyr::mapvalues(final_heatmap_2$region, from = c("5utr", "exon1", "intron1", "exon", "intron", "exonL", "intronL", "3utr"), to = c("5' UTR", "first\nexon", "first\nintron", "exon", "intron", "last\nexon", "last\nintron", "3' UTR"))
final_heatmap_2$rbps <- paste0(final_heatmap_2$rbp1,"::",final_heatmap_2$rbp2)
# Reorder
final_heatmap_2$pair <- factor(final_heatmap_2$pair, levels=rev(levels(final_heatmap_2$pair)))

# Plot
ggplot(final_heatmap_2, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_gradient2() + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v8-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5, units="mm")
ggplot(final_heatmap_2, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_gradient2("Enrichment", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v8.0-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5, units="mm")
ggplot(final_heatmap_2, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v8.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5, units="mm")
ggplot(final_heatmap_2, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", option="C", limits=c(-1,1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
ggsave(paste0("output-plot-metagene-v8.4-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5, units="mm")
# v8.5
ggplot(final_heatmap_2, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_gradient2("Enrichment", limits=c(-1,1)) + facet_grid(rows=vars(rbps), cols=vars(region), scales="free_y") + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0), strip.text.y = element_blank())
ggsave(paste0("output-plot-metagene-v8.5-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183*2, height=91.5*8, units="mm")
ggplot(final_heatmap_2, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", option="D", limits=c(-1,1)) + facet_grid(rows=vars(rbps), cols=vars(region), scales="free_y") + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0), strip.text.y = element_blank())
ggsave(paste0("output-plot-metagene-v8.5.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*4, units="mm")
ggplot(final_heatmap_2, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", option="C", limits=c(-1,1)) + facet_grid(rows=vars(rbps), cols=vars(region), scales="free_y") + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0), strip.text.y = element_blank())
ggsave(paste0("output-plot-metagene-v8.5.4-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*4, units="mm")
# # v8.6
# levels(final_heatmap_2$pair) <- sub("^(\\w+)\\|(\\w+): (\\w+)$", "\\3 only vs. \\1::\\2", levels(final_heatmap_2$pair_backup), perl=T)
# final_heatmap_2$value <- -final_heatmap_2$value_backup
# ggplot(final_heatmap_2, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", option="D", limits=c(-1,1)) + facet_grid(rows=vars(rbps), cols=vars(region), scales="free_y") + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0), strip.text.y = element_blank())
# ggsave(paste0("output-plot-metagene-v8.6.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5, units="mm")
# ggplot(final_heatmap_2, aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() + scale_fill_viridis_c("Enrichment", option="C", limits=c(-1,1)) + facet_grid(rows=vars(rbps), cols=vars(region), scales="free_y") + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0), strip.text.y = element_blank())
# ggsave(paste0("output-plot-metagene-v8.6.4-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5, units="mm")









# Write pair distances to file
pair_distance
write.table(pair_distance, paste0("output-table-pair-distances.tsv"), sep="\t", quote=F, row.names=F)






# # Plot pair distances as box plots
head(pair_distance)
pair_distance$rbps <- paste0(pair_distance$rbp1,"::",pair_distance$rbp2)
head(pair_distance)
summary(pair_distance)
# ggplot(pair_distance, aes(x=1, y=value, fill=value)) + geom_boxplot(notch=T) + coord_cartesian(ylim=c(-1, 1)) + facet_grid(cols=vars(region)) + xlab("") + ylab("Pair distance scores: complex vs. individual RBPs") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0))
# ggplot(pair_distance, aes(x=region, y=value, fill=value)) + geom_boxplot(notch=T) + coord_cartesian(ylim=c(-1, 1)) + xlab("") + ylab("Pair distance scores: complex vs. individual RBPs") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(pair_distance, aes(x=region, y=value, fill=value)) + geom_boxplot(notch=T) + coord_cartesian(ylim=c(-max(abs(pair_distance$value)), max(abs(pair_distance$value)))) + xlab("") + ylab("Pair distances: complex vs. individual RBPs") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(paste0("output-plot-metagene-v9-square-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")
ggplot(pair_distance, aes(x=region, y=value, fill=value)) + geom_boxplot(notch=T) + coord_cartesian(ylim=c(-max(abs(pair_distance$value)), max(abs(pair_distance$value)))) + xlab("") + ylab("Pair distances: complex vs. individual RBPs") + ggtitle("") + theme_minimal()
ggsave(paste0("output-plot-metagene-v9-wide-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=91.5*2, height=91.5, units="mm")



# # final_heatmap_3 / v11 (Bar charts with error bars)
# # library(tidyr)
# 
# final_heatmap_3$pair <- paste0(final_heatmap_3$rbp1,"|",final_heatmap_3$rbp2)
# final_heatmap_3
# 
# final_heatmap_3 %>% group_by(pair, region, title, .drop=F) %>% summarise(mean = mean(freq), n = n(), sd = sd(freq), se = sd / n)
# final_heatmap_3_1 <- final_heatmap_3 %>% group_by(pair, region, title, .drop=F) %>% summarise(mean = mean(freq), n = n(), sd = sd(freq), se = sd / n)
# # final_heatmap_3_1 <- final_heatmap_3 %>% group_by(region, title, .drop=F) %>% summarise(mean = mean(freq), n = n(), sd = sd(freq), se = sd / n)
# 
# # ggplot(q, aes(x=region, colour=title, fill=title)) + geom_bar(aes(width=width), position="dodge") + scale_fill_manual(aesthetics = c("colour", "fill"), values=c("#8ea106", "#47d9bf", "#00a1d9")) + theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("mRNA region") + ylab("Probability density") + guides(colour=F, fill=F) + ggtitle(paste0(rbp1," & ",rbp2)) + theme(plot.title = element_text(hjust = 0.5))
# 
# ggplot(final_heatmap_3_1, aes(x=region, y=mean, fill=title)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) + scale_fill_manual(values=c("#8ea106", "#47d9bf", "#00a1d9")) + theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("mRNA region") + ylab("Probability density") + guides(colour=F, fill=F) + ggtitle(paste0(rbp1," & ",rbp2)) + theme(plot.title = element_text(hjust = 0.5))
# ggsave(paste0("output-plot-metagene-v11-bars-",clip_cobinding,"-",type,"-",style,"-",rbp1,"-",rbp2,".pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")





# v13: Plot only the pairs from the v12 RBP1-Pair-RBP2 for some RBPs of interest (hubs in the network)
print(alldensities, n=100)
# rbp_of_interest <- "PCBP1"
rbp_of_interest <- "APOBEC3C"
rbp_occurrences <- table(c(alldensities$rbp1, alldensities$rbp2))
names(rbp_occurrences[rbp_occurrences>=2])
# for (rbp_of_interest in c("PCBP1", "RBFOX2", "NONO", "CSTF2T", ""))
for (rbp_of_interest in names(rbp_occurrences[rbp_occurrences>=2])) {
  # alldensities %>% filter(rbp1==rbp_of_interest | rbp2==rbp_of_interest)
  # alldensities %>% filter(rbp1==rbp_of_interest | rbp2==rbp_of_interest) %>% unnest(cols=densities) %>% print(n=100)
  # alldensities %>% filter(rbp1==rbp_of_interest | rbp2==rbp_of_interest) %>% unnest(cols=densities) %>% filter(grepl(" & ", title, perl=T))
  # alldensities %>% filter(rbp1==rbp_of_interest | rbp2==rbp_of_interest) %>% unnest(cols=densities) %>% filter(grepl(" & ", title, perl=T)) %>%
  #   ggplot(aes(x=globalpos, y=q50, colour=title)) + geom_line()
  
  
  
  # Extract relevant rows & unnest
  tmpalldensities <- alldensities %>% filter(rbp1==rbp_of_interest | rbp2==rbp_of_interest) %>% unnest(cols=densities) %>% filter(grepl(" & ", title, perl=T))

  
  
  # Place the rbp_of_interest first in the title
  # rbp1_rbp2_title <- paste0(rbp1," & ",rbp2,"\n",nrow(rbp1_rbp2)," sites\n",length(unique(rbp1_rbp2$ensgv))," genes")
  tmpre <- paste0('^(\\S+) & ',rbp_of_interest,'\n(\\d+) sites\n(\\d+) genes$')
  tmpre
  # Not expecting all to match now (I'm only matching the ones in the incorrect order, where the rbp_of_interest comes second)
  # if (unique(str_detect(tmpalldensities$title, tmpre)) != T) {
  #   stop("Error: Couldn't match some titles in tmpalldensities")
  # }
  tmpalldensities$title <- str_replace(tmpalldensities$title, tmpre, paste0(rbp_of_interest,' & \\1\n\\2 sites\n\\3 genes'))
  # str_match(tmpalldensities$title[2500], tmpre)
  str(tmpalldensities$title)
  unique(tmpalldensities$title)
  
  
  
  # Add facet_grid strip text spacing so the "split" plots will all have the same width  
  # spacing_length <- max(nchar(tmpalldensities$title))
  # >> doesn't work because of the linebreaks, and it's a variable width font anyway
  # spacing <- "                                               "
  spacing <- paste0(rep(" ", 75), collapse="")
  
  
  
  
  p <- tmpalldensities %>%
    # Apply spacing ("exttitle" will only be used in the "split" version of the plot)
    mutate(exttitle=paste0(spacing,"\n",title,"\n",spacing)) %>%
    # Plot
    ggplot(aes(x=globalpos, y=q50, colour=title, fill=title)) +
    # ggplot(aes(x=globalpos, y=q50, colour=paste0("                             \n",title,"                             \n"), fill=paste0("                             \n",title,"                             \n"))) +
    # facet_grid(cols=vars(region)) +
    # facet_grid(rows=vars(title)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha=0.3, colour=NA) +
    # stat_density(data = qtib,
    #              # aes(pos, ..density.., color = "raw density"),
    #              aes(globalpos, ..density.., color = "raw density"),
    #              size = 2, geom = "line") +
    geom_line() +
    # scale_color_manual(values = c("red", "black")) +
    coord_cartesian(xlim=c(0, 1)) + 
    scale_colour_viridis_d(NULL) +
    scale_fill_viridis_d(NULL) +
    # scale_colour_manual(NULL, aesthetics=c("colour", "fill"), values=c("#8ea106", "#00a1d9", "#47d9bf", "#04518c")) +
    # scale_colour_manual(NULL, aesthetics=c("colour", "fill"), values=c("#47d9bf", "#00a1d9", "#04518c", "#003056")) + 
    geom_vline(xintercept=c(1:7/8), linetype="dashed", colour="#cccccc") + 
    scale_x_continuous(breaks=1:8/8-1/16, labels=levels(q$region)) + 
    labs(y = "density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + 
    theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("mRNA region") + 
    ylab("Probability density") + 
    # guides(colour=F, fill=F) + 
    # ggtitle(paste0(rbp1," & ",rbp2)) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="top")
  # theme(legend.position="right")
  p
  str(p)
  ggsave(paste0("output-plot-metagene-v13-combined-",clip_cobinding,"-",type,"-",style,"-",rbp_of_interest,".pdf"), device=cairo_pdf, width=91.5*1.5, height=91.5*3/4, units="mm")
  ggsave(paste0("output-plot-metagene-v13-combined-narrow-",clip_cobinding,"-",type,"-",style,"-",rbp_of_interest,".pdf"), device=cairo_pdf, width=91.5, height=91.5*3/4, units="mm")
  
  # Split
  # split_height <- 91.5 * 3/4 + (rbp_occurrences[rbp_of_interest] - 2) * 91.5 * 3/4 / 2
  split_height <- 91.5 * 3/4 + (rbp_occurrences[rbp_of_interest] - 2) * 26.9434 # actual plot area
  p + facet_grid(rows=vars(exttitle), scales="free_y") + guides(colour=F, fill=F) + theme(strip.text.y = element_text(hjust = 0))
  # p + facet_grid(rows=vars(title), scales="free_y") + theme(strip.background = element_blank(), strip.text.y = element_blank())
  ggsave(paste0("output-plot-metagene-v13-split-",clip_cobinding,"-",type,"-",style,"-",rbp_of_interest,".pdf"), device=cairo_pdf, width=91.5*2, height=split_height, units="mm")
}







# Larger clique of interest (NONO-SFPQ-EWSR1-CSTF2T):
rbps_of_interest <- c("NONO", "SFPQ", "EWSR1", "CSTF2T")
if (alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% nrow > 0) {
  alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest)
  alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% unnest(cols=densities) %>% print(n=100)
  alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% unnest(cols=densities) %>% filter(grepl(" & ", title, perl=T))
  alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% unnest(cols=densities) %>% filter(grepl(" & ", title, perl=T)) %>%
    ggplot(aes(x=globalpos, y=q50, colour=title)) + geom_line()
  
  
  p <- alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% unnest(cols=densities) %>% filter(grepl(" & ", title, perl=T)) %>%
    ggplot(aes(x=globalpos, y=q50, colour=title, fill=title)) +
    # facet_grid(cols=vars(region)) +
    # facet_grid(rows=vars(title)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha=0.3, colour=NA) +
    # stat_density(data = qtib,
    #              # aes(pos, ..density.., color = "raw density"),
    #              aes(globalpos, ..density.., color = "raw density"),
    #              size = 2, geom = "line") +
    geom_line() +
    # scale_color_manual(values = c("red", "black")) +
    coord_cartesian(xlim=c(0, 1)) + 
    scale_colour_viridis_d(NULL) +
    scale_fill_viridis_d(NULL) +
    # scale_colour_manual(NULL, aesthetics=c("colour", "fill"), values=c("#8ea106", "#00a1d9", "#47d9bf", "#04518c")) +
    # scale_colour_manual(NULL, aesthetics=c("colour", "fill"), values=c("#47d9bf", "#00a1d9", "#04518c", "#003056")) + 
    geom_vline(xintercept=c(1:7/8), linetype="dashed", colour="#cccccc") + 
    scale_x_continuous(breaks=1:8/8-1/16, labels=levels(q$region)) + 
    labs(y = "density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + 
    theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("mRNA region") + 
    ylab("Probability density") + 
    # guides(colour=F, fill=F) + 
    # ggtitle(paste0(rbp1," & ",rbp2)) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="top")
  # theme(legend.position="right")
  p
  ggsave(paste0("output-plot-metagene-v13-combined-",clip_cobinding,"-",type,"-",style,"-",paste0(rbps_of_interest, collapse="-"),".pdf"), device=cairo_pdf, width=91.5*1.5, height=91.5, units="mm")
  
  # Split
  # split_height <- 91.5 * 3/4 + (rbp_occurrences[rbp_of_interest] - 2) * 91.5 * 3/4 / 2
  # split_height <- 91.5 * 3/4 + (nrow(alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest)) - 2) * 26.9434 # actual plot area
  split_height <- 91.5 * 1.5
  p + facet_grid(rows=vars(title), scales="free_y") + guides(colour=F, fill=F)
  ggsave(paste0("output-plot-metagene-v13-split-",clip_cobinding,"-",type,"-",style,"-",paste0(rbps_of_interest, collapse="-"),".pdf"), device=cairo_pdf, width=91.5*1.5, height=split_height, units="mm")
}











# Larger clique of interest (PTBP1-centric: PTBP1-IGF2BP1-IGF2BP2-PCBP1-QKI-TARDBP-ZC3H8):
rbps_of_interest <- c("PTBP1", "IGF2BP1", "IGF2BP2", "PCBP1", "QKI", "TARDBP", "ZC3H8")
if (alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% nrow > 0) {
  alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest)
  alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% unnest(cols=densities) %>% print(n=100)
  alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% unnest(cols=densities) %>% filter(grepl(" & ", title, perl=T))
  alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% unnest(cols=densities) %>% filter(grepl(" & ", title, perl=T)) %>%
    ggplot(aes(x=globalpos, y=q50, colour=title)) + geom_line()
  
  
  p <- alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest) %>% unnest(cols=densities) %>% filter(grepl(" & ", title, perl=T)) %>%
    ggplot(aes(x=globalpos, y=q50, colour=title, fill=title)) +
    # facet_grid(cols=vars(region)) +
    # facet_grid(rows=vars(title)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha=0.3, colour=NA) +
    # stat_density(data = qtib,
    #              # aes(pos, ..density.., color = "raw density"),
    #              aes(globalpos, ..density.., color = "raw density"),
    #              size = 2, geom = "line") +
    geom_line() +
    # scale_color_manual(values = c("red", "black")) +
    coord_cartesian(xlim=c(0, 1)) + 
    scale_colour_viridis_d(NULL) +
    scale_fill_viridis_d(NULL) +
    # scale_colour_manual(NULL, aesthetics=c("colour", "fill"), values=c("#8ea106", "#00a1d9", "#47d9bf", "#04518c")) +
    # scale_colour_manual(NULL, aesthetics=c("colour", "fill"), values=c("#47d9bf", "#00a1d9", "#04518c", "#003056")) + 
    geom_vline(xintercept=c(1:7/8), linetype="dashed", colour="#cccccc") + 
    scale_x_continuous(breaks=1:8/8-1/16, labels=levels(q$region)) + 
    labs(y = "density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5)) + 
    theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("mRNA region") + 
    ylab("Probability density") + 
    # guides(colour=F, fill=F) + 
    # ggtitle(paste0(rbp1," & ",rbp2)) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="top")
  # theme(legend.position="right")
  p
  ggsave(paste0("output-plot-metagene-v13-combined-",clip_cobinding,"-",type,"-",style,"-",paste0(rbps_of_interest, collapse="-"),".pdf"), device=cairo_pdf, width=91.5*1.5, height=91.5, units="mm")
  
  # Split
  # split_height <- 91.5 * 3/4 + (rbp_occurrences[rbp_of_interest] - 2) * 91.5 * 3/4 / 2
  # split_height <- 91.5 * 3/4 + (nrow(alldensities %>% filter(rbp1 %in% rbps_of_interest | rbp2 %in% rbps_of_interest)) - 2) * 26.9434 # actual plot area
  split_height <- 91.5 * 1.5
  p + facet_grid(rows=vars(title), scales="free_y") + guides(colour=F, fill=F)
  ggsave(paste0("output-plot-metagene-v13-split-",clip_cobinding,"-",type,"-",style,"-",paste0(rbps_of_interest, collapse="-"),".pdf"), device=cairo_pdf, width=91.5*1.5, height=split_height, units="mm")
}











# v14: Plot a 7.3.5-style heatmap of the pairs only
# alldensities %>% 
#   unnest(cols=densities) %>% 
#   filter(grepl(" & ", title, perl=T)) %>%
#   ggplot(aes(x=globalpos, y=title)) +
#   geom_tile()

tmp_final_heatmap_4 <- tibble(final_heatmap_4)
tmp_final_heatmap_4 %<>%
  filter(grepl(" & ", pair, perl=T)) %>%
  mutate(pair=str_replace(pair, '^(\\S+) & (\\S+)\n\\d+ sites\\n\\d+ genes$', '\\1 & \\2'))

# Hierarchically cluster the heatmap
sortmap <- tmp_final_heatmap_4
sortmap$raw_value <- NULL
asortmap <- acast(sortmap, pair~region+bin, value.var="value")
asortmap
set.seed(1)
asortmap_hclust <- hclust(dist(asortmap), method="complete")
plot(asortmap_hclust, hang=-1)
pdf(paste0("output-plot-metagene-v14.0-",clip_cobinding,"-",type,"-",style,"-all.pdf"), width=5, height=4)
plot(asortmap_hclust, hang=-1)
dev.off()
target_order_hclust <- asortmap_hclust$labels[c(asortmap_hclust$order)]
target_order_hclust
target_order_hclust <- factor(target_order_hclust, levels=target_order_hclust)
target_order_hclust
# Apply target order
tmp_final_heatmap_4$pair <- as.character(tmp_final_heatmap_4$pair)
target_order_hclust <- as.character(target_order_hclust)
target_order_hclust
tmp_final_heatmap_4 <- left_join(data.frame(pair=target_order_hclust, stringsAsFactors=F), tmp_final_heatmap_4, by="pair")
tmp_final_heatmap_4$pair <- factor(tmp_final_heatmap_4$pair, levels=target_order_hclust)
str(tmp_final_heatmap_4)

# Plot
tmp_final_heatmap_4 %>%
  ggplot(aes(x=as.numeric(as.character(bin)), y=pair, fill=value)) + geom_tile() +
  # scale_fill_gradientn("Enrichment of\nthe complex vs.\nindividual RBPs", values=scales::rescale(c(min(final_heatmap$value), -0.2, 0, 0.2, 0.4, max(final_heatmap$value))),
  scale_fill_gradientn("Probability\ndensity",
                       # colours=brewer.pal(n=7, name="RdBu")[2:7]) +
                       # colours=c("#B2182B", "#EF8A62", "#FFFFFF", "#67A9CF", "#2166AC", "#053061")) +
                       # colours=c("#EF8A62", "#EF8A62", "#FFFFFF", "#67A9CF", "#0571B0", "#0571B0")) +
                       colours=c("#FFFFFF", "#4393C3", "#053061", "#053061")) +
  # colours=c("#CA0020", "#CA0020", "#FFFFFF", "#0571B0", "#053061", "#053061")) +
  # colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
  # colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
  # colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
  facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0.5)) +
  theme(legend.title = element_text(size = 8.8), legend.text = element_text(size = 8.8)) +
  # guides(fill = guide_legend(override.aes = list(size = 0.5)))
  theme(legend.key.size = unit(5, "mm")) +
  ggsave(paste0("output-plot-metagene-v14.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*1.5, units="mm")






# v15: Plot a 7.3.5-style heatmap of the individual proteins
# alldensities %>% 
#   unnest(cols=densities) %>% 
#   filter(grepl(" & ", title, perl=T)) %>%
#   ggplot(aes(x=globalpos, y=title)) +
#   geom_tile()

final_heatmap_5
tmp_final_heatmap_5 <- tibble(final_heatmap_5)

# Hierarchically cluster the heatmap
sortmap <- tmp_final_heatmap_5
sortmap$raw_value <- NULL
asortmap <- acast(sortmap, rbp~region+bin, value.var="value")
asortmap
set.seed(1)
asortmap_hclust <- hclust(dist(asortmap), method="complete")
plot(asortmap_hclust, hang=-1)
pdf(paste0("output-plot-metagene-v15.0-",clip_cobinding,"-",type,"-",style,"-all.pdf"), width=5, height=4)
plot(asortmap_hclust, hang=-1)
dev.off()
target_order_hclust <- asortmap_hclust$labels[c(asortmap_hclust$order)]
target_order_hclust
target_order_hclust <- factor(target_order_hclust, levels=target_order_hclust)
target_order_hclust
# Apply target order
tmp_final_heatmap_5$rbp <- as.character(tmp_final_heatmap_5$rbp)
target_order_hclust <- as.character(target_order_hclust)
target_order_hclust
tmp_final_heatmap_5 <- left_join(data.frame(rbp=target_order_hclust, stringsAsFactors=F), tmp_final_heatmap_5, by="rbp")
tmp_final_heatmap_5$rbp <- factor(tmp_final_heatmap_5$rbp, levels=target_order_hclust)
str(tmp_final_heatmap_5)

# Plot
tmp_final_heatmap_5 %>%
  ggplot(aes(x=as.numeric(as.character(bin)), y=rbp, fill=value)) + geom_tile() +
  # scale_fill_gradientn("Enrichment of\nthe complex vs.\nindividual RBPs", values=scales::rescale(c(min(final_heatmap$value), -0.2, 0, 0.2, 0.4, max(final_heatmap$value))),
  scale_fill_gradientn("Probability\ndensity",
                       # colours=brewer.pal(n=7, name="RdBu")[2:7]) +
                       # colours=c("#B2182B", "#EF8A62", "#FFFFFF", "#67A9CF", "#2166AC", "#053061")) +
                       # colours=c("#EF8A62", "#EF8A62", "#FFFFFF", "#67A9CF", "#0571B0", "#0571B0")) +
                       colours=c("#FFFFFF", "#4393C3", "#053061", "#053061")) +
  # colours=c("#CA0020", "#CA0020", "#FFFFFF", "#0571B0", "#053061", "#053061")) +
  # colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
  # colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
  # colours=c("#D6604D", "#D6604D", "#FFFFFF", "#4393C3", "#053061", "#053061")) +
  facet_grid(cols=vars(region)) + xlab("") + ylab("") + ggtitle("") + theme_minimal() + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(vjust = 0.5)) +
  theme(legend.title = element_text(size = 8.8), legend.text = element_text(size = 8.8)) +
  # guides(fill = guide_legend(override.aes = list(size = 0.5)))
  theme(legend.key.size = unit(5, "mm")) +
  ggsave(paste0("output-plot-metagene-v15.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=91.5*1.5, units="mm")





# v16: Correlation matrix
# Generate density correlations between individual proteins (all their sites) and complexes
alldensities
alldensities %>% unnest(cols=densities) %>% filter(grepl(" & ", title) | is.na(rbp2))
alldensities %>% unnest(cols=densities) %>% filter(grepl(" & ", title) | is.na(rbp2)) %>% pull(title) %>% unique()
alldensities %>% unnest(cols=densities) %>% filter(grepl(" & ", title) | is.na(rbp2)) %>% select(title, globalpos, q50)

m <- alldensities %>%
  unnest(cols=densities) %>%
  # Remove site and gene numbers
  mutate(title = str_replace(title, "\n\\d+ sites\n\\d+ genes$", "")) %>% 
  filter(grepl(" & ", title) | is.na(rbp2)) %>%
  select(title, globalpos, q50) %>%
  pivot_wider(names_from=title, values_from=c(q50)) %>%
  select(-globalpos) %>%
  as.matrix() %>%
  rcorr(type = "pearson")
str(m)

# corrplot: Looks cool, but need to use pdf():
corrplot(m$r, method="color", type="upper", order="hclust", diag=F, tl.cex=0.4, tl.col="black", tl.offset=0.3)
if (length(rbps) < 30) {
  # <30: Plot small
  pdf(paste0("output-plot-metagene-v16.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), width=183/2 / 25.4, height=183/2 / 25.4)
} else {
  # >=30: Plot big
  pdf(paste0("output-plot-metagene-v16.1-",clip_cobinding,"-",type,"-",style,"-all.pdf"), width=183 / 25.4, height=183 / 25.4)
}
corrplot(m$r, method="color", type="upper", order="hclust", diag=F, tl.cex=0.4, tl.col="black", tl.offset=0.3)
dev.off()


# corrplot full square:
corrplot(m$r, method="color", order="hclust", tl.cex=0.4, tl.col="black", tl.offset=0.3)
if (length(rbps) < 30) {
  # <30: Plot small
  pdf(paste0("output-plot-metagene-v16.2-",clip_cobinding,"-",type,"-",style,"-all.pdf"), width=183/2 / 25.4, height=183/2 / 25.4)
} else {
  # >=30: Plot big
  pdf(paste0("output-plot-metagene-v16.2-",clip_cobinding,"-",type,"-",style,"-all.pdf"), width=183 / 25.4, height=183 / 25.4)
}
corrplot(m$r, method="color", order="hclust", tl.cex=0.4, tl.col="black", tl.offset=0.3)
dev.off()


# Alternative using ggcorrplot:
# p.mat=m$P, 
ggcorrplot(m$r, hc.order=T, type="upper", outline.col="white", show.diag=F) + theme_minimal() +
  theme(axis.text.x=element_text(size=5, angle=90, hjust=1, vjust=0)) + 
  theme(axis.text.y=element_text(size=5)) +
  # scale_fill_gradient2(low="red", mid="white", high="blue", limits=c(-1,1))
  # scale_fill_gradient2(low="#67001f", mid="#ffffff", high="#053061", limits=c(-1,1))
  scale_fill_gradient2(low="#67001f", mid="#ffffff", high="#053061", limits=c(-1,1))
if (length(rbps) < 30) {
  # <30: Plot small
  ggsave(paste0("output-plot-metagene-v16.3-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183/2, height=183/2, units="mm")
} else {
  # >=30: Plot big
  ggsave(paste0("output-plot-metagene-v16.3-",clip_cobinding,"-",type,"-",style,"-all.pdf"), device=cairo_pdf, width=183, height=183, units="mm")
}


# colours=c("#B2182B", "#EF8A62", "#FFFFFF", "#67A9CF", "#2166AC", "#053061")) +
# colours=c("#EF8A62", "#EF8A62", "#FFFFFF", "#67A9CF", "#0571B0", "#0571B0")) +
# colours=c("#FFFFFF", "#4393C3", "#053061", "#053061")) +
  
# ?pivot_wider

as.matrix(alldensities %>% unnest(cols=densities) %>% filter(grepl(" & ", title) | is.na(rbp2)) %>% select(title, globalpos))
alldensities %>% unnest(cols=densities) %>% filter(grepl(" & ", title) | is.na(rbp2)) %>% select(title, q50)







# v17: Correlation scatter plots for RBP pairs

# mytitles <- 
#   alldensities %>% 
#   unnest(cols=densities) %>% 
#   # filter(grepl(" & ", title) | is.na(rbp2)) %>%
#   filter(is.na(rbp2)) %>%
#   # Remove site and gene numbers
#   mutate(shorttitle = str_replace(title, "\n\\d+ sites\n\\d+ genes$", "")) %>% 
#   pull(shorttitle) %>% 
#   unique
# tibble(mytitles)

# setA_df %>%
#   select(-gene_symbol) %>%
#   cor(setB_df %>% select(-gene_symbol)) %>%
#   melt()  
# 
# title1 <- "FXR1"
# title2 <- "FXR2"
# title1 <- "PCBP1"
# title2 <- "APOBEC3C"
title1 <- "U2AF1"
title2 <- "U2AF2"
# title1 <- "FUBP3"
# title2 <- "LARP4"
# title1 <- "PCBP1"
# title2 <- "PTBP1"
title1 <- "U2AF2"
title2 <- "U2AF1"
# #DEBUG
# mytitles <- c()
# #END DEBUG

mytitles <- rbps
#DEBUG
# PCBP1-HNRNPK vs PCBP1-IGF2BP2
# mytitles <- c("PCBP1 & HNRNPK", "PCBP1 & IGF2BP2")
#END DEBUG
# #DEBUG
# # PCBP1-RBFOX2 vs  PCBP1-APOBEC3C
# mytitles <- c("PCBP1 & RBFOX2", "PCBP1 & APOBEC3C")
# #END DEBUG

for (title1 in mytitles) {
  for (title2 in mytitles) {
    
    if (title2 <= title1) {
      next
    }
    if (title1 == title2) {
      next
    }
    
    print(glue("{title1} vs. {title2}"))
    complextitle <- glue("{title1} & {title2}")

    # Get correlation values
    # alldensities
    # alldensities %>% unnest(cols=densities)
    m <- alldensities %>%
      unnest(cols=densities) %>% 
      # Get complexes or individual RBPs (but not "RBP1 only")
      filter(grepl(" & ", title) | is.na(rbp2)) %>% 
      # Remove site and gene numbers
      mutate(title = str_replace(title, "\n\\d+ sites\n\\d+ genes$", "")) %>%
      # filter(title==title1 | title==title2) %>%
      filter(title==title1 | title==title2 | title==glue("{title1} & {title2}")) %>%
      select(title, globalpos, q50) %>%
      pivot_wider(names_from=title, values_from=c(q50)) %>%
      select(-globalpos) %>%
      as.matrix() %>%
      rcorr(type = "pearson") %>%
      # rcorr(type = "spearman") %>%
      tidy
    m
    m$estimate
    m$p.value
    
    # Plot RBP1 vs RBP2
    tmpd <- alldensities %>%
      unnest(cols=densities) %>% 
      # Get complexes or individual RBPs (but not "RBP1 only")
      filter(grepl(" & ", title) | is.na(rbp2)) %>% 
      # Remove site and gene numbers
      mutate(title = str_replace(title, "\n\\d+ sites\n\\d+ genes$", "")) %>%
      filter(title==title1 | title==title2) %>%
      select(title, globalpos, q50) %>%
      # Make column names generic (i.e. rbp1 instead of U2AF1) for ggplot
      mutate(title = str_replace(title, glue("^{title1}$"), "rbp1")) %>%
      mutate(title = str_replace(title, glue("^{title2}$"), "rbp2")) %>%
      pivot_wider(names_from=title, values_from=c(q50))

      ggplot(tmpd, aes(x=rbp1, y=rbp2)) +
      geom_point(shape=20) +
      stat_smooth(method="lm", colour="grey", formula=y~x) +
      # stat_smooth(method="lm", formula=y~x) +
      # geom_abline(slope=1, intercept=0, colour="grey") +
      # scale_x_log10() +
      # scale_y_log10() +
      scale_fill_brewer(palette="RdBu") +
      # coord_fixed(xlim=c(0, max(max(tmpd$rbp1, na.rm=T), max(tmpd$rbp2, na.rm=T), na.rm=T)), ylim=c(0, max(max(tmpd$rbp1, na.rm=T), max(tmpd$rbp2, na.rm=T), na.rm=T))) +
      theme_minimal() +
      theme(plot.subtitle = element_text(hjust = 1)) +
      labs(x=title1, y=title2, subtitle = glue("r={round(m %>% filter(column1==title1 & column2==title2) %>% pull(estimate), 2)}  p={sprintf('%.0g', m %>% filter(column1==title1 & column2==title2) %>% pull(p.value))}")) +
      
      ggsave(paste0("output-plot-metagene-v17.1-",clip_cobinding,"-",type,"-",style,"-",title1,"-",title2,".pdf"), device=cairo_pdf, width=91.5/2, height=91.5/2, units="mm")
    
    
    # Plot RBP1 vs complex
    tmpd <- alldensities %>%
      unnest(cols=densities) %>%
      # Get complexes or individual RBPs (but not "RBP1 only")
      filter(grepl(" & ", title) | is.na(rbp2)) %>%
      # Remove site and gene numbers
      mutate(title = str_replace(title, "\n\\d+ sites\n\\d+ genes$", "")) %>%
      filter(title==title1 | title==glue("{title1} & {title2}")) %>%
      select(title, globalpos, q50) %>%
      # Make column names generic (i.e. rbp1 instead of U2AF1) for ggplot
      mutate(title = str_replace(title, glue("^{title1}$"), "rbp1")) %>%
      mutate(title = str_replace(title, glue("^{title2}$"), "rbp2")) %>%
      mutate(title = str_replace(title, glue("^{title1} & {title2}$"), "mycomplex")) %>%
      pivot_wider(names_from=title, values_from=c(q50))
    
    # If the complex exists in this orientation, make a plot
    if (ncol(tmpd) == 3) {
      ggplot(tmpd, aes(x=rbp1, y=mycomplex)) +
        geom_point(shape=20) +
        stat_smooth(method="lm", colour="grey", formula=y~x) +
        # stat_smooth(method="lm", formula=y~x) +
        # geom_abline(slope=1, intercept=0, colour="grey") +
        # scale_x_log10() +
        # scale_y_log10() +
        scale_fill_brewer(palette="RdBu") +
        # coord_fixed(xlim=c(0, max(max(tmpd$rbp1, na.rm=T), max(tmpd$mycomplex, na.rm=T), na.rm=T)), ylim=c(0, max(max(tmpd$rbp1, na.rm=T), max(tmpd$mycomplex, na.rm=T), na.rm=T))) +
        theme_minimal() +
        theme(plot.subtitle = element_text(hjust = 1)) +
        labs(x=title1, y=complextitle, subtitle = glue("r={round(m %>% filter(column1==title1 & column2==complextitle) %>% pull(estimate), 2)}  p={sprintf('%.0g', m %>% filter(column1==title1 & column2==complextitle) %>% pull(p.value))}")) +
        
        ggsave(paste0("output-plot-metagene-v17.2-",clip_cobinding,"-",type,"-",style,"-",title1,"-",glue("{title1} & {title2}"),".pdf"), device=cairo_pdf, width=91.5/2, height=91.5/2, units="mm")
    }



    # Plot RBP2 vs complex
    tmpd <- alldensities %>%
      unnest(cols=densities) %>%
      # Get complexes or individual RBPs (but not "RBP1 only")
      filter(grepl(" & ", title) | is.na(rbp2)) %>%
      # Remove site and gene numbers
      mutate(title = str_replace(title, "\n\\d+ sites\n\\d+ genes$", "")) %>%
      filter(title==title2 | title==glue("{title1} & {title2}")) %>%
      select(title, globalpos, q50) %>%
      # Make column names generic (i.e. rbp1 instead of U2AF1) for ggplot
      mutate(title = str_replace(title, glue("^{title1}$"), "rbp1")) %>%
      mutate(title = str_replace(title, glue("^{title2}$"), "rbp2")) %>%
      mutate(title = str_replace(title, glue("^{title1} & {title2}$"), "mycomplex")) %>%
      pivot_wider(names_from=title, values_from=c(q50))

    # If the complex exists in this orientation, make a plot
    if (ncol(tmpd) == 3) {
      ggplot(tmpd, aes(x=rbp2, y=mycomplex)) +
      geom_point(shape=20) +
      stat_smooth(method="lm", colour="grey", formula=y~x) +
      # stat_smooth(method="lm", formula=y~x) +
      # geom_abline(slope=1, intercept=0, colour="grey") +
      # scale_x_log10() +
      # scale_y_log10() +
      scale_fill_brewer(palette="RdBu") +
      # coord_fixed(xlim=c(0, max(max(tmpd$rbp2, na.rm=T), max(tmpd$mycomplex, na.rm=T), na.rm=T)), ylim=c(0, max(max(tmpd$rbp2, na.rm=T), max(tmpd$mycomplex, na.rm=T), na.rm=T))) +
      theme_minimal() +
      theme(plot.subtitle = element_text(hjust = 1)) +
      labs(x=title2, y=complextitle, subtitle = glue("r={round(m %>% filter(column1==title2 & column2==complextitle) %>% pull(estimate), 2)}  p={sprintf('%.0g', m %>% filter(column1==title2 & column2==complextitle) %>% pull(p.value))}")) +
        
      ggsave(paste0("output-plot-metagene-v17.3-",clip_cobinding,"-",type,"-",style,"-",title2,"-",glue("{title1} & {title2}"),".pdf"), device=cairo_pdf, width=91.5/2, height=91.5/2, units="mm")
    }
    
  }
}

write_rds(alldensities, paste0("output-plot-metagene-v18-",clip_cobinding,"-",type,"-",style,"-alldensities.rds"))
alldensities_backup <- alldensities
alldensities <- readRDS(paste0("output-plot-metagene-v18-",clip_cobinding,"-",type,"-",style,"-alldensities.rds"))



# v19: Correlation r values
# Parameter: number of distinct pairs an RBP needs to be in to be included here (i.e. its degree in the Figure 3A network)
# min_pairs <- 3
# min_pairs <- 2
# min_pairs <- 1

# Generate density correlations between individual proteins (all their sites) and complexes
m <- alldensities %>%
  unnest(cols=densities) %>%
  # Remove site and gene numbers
  mutate(title = str_replace(title, "\n\\d+ sites\n\\d+ genes$", "")) %>% 
  filter(grepl(" & ", title) | is.na(rbp2)) %>%
  select(title, globalpos, q50) %>%
  pivot_wider(names_from=title, values_from=c(q50)) %>%
  select(-globalpos) %>%
  as.matrix() %>%
  rcorr(type = "pearson")
r <- melt(m$r)
colnames(r) <- c("rbp1", "rbp2", "r")
r <- as_tibble(r)


# Get RBPs of interest (with ≥ min_pairs interactions in the Figure 3A network, i.e. the strictly filtered one (for significance of co-binding))
rbps_of_interest <- tibble(rbp_pairs) %>% 
  unnest(cols=c(rbp_pairs)) %>% 
  rename(rbp = rbp_pairs) %>% 
  group_by(rbp) %>% 
  tally %>%
  # filter(n >= min_pairs) %>%
  arrange(desc(n))
# %>% pull(rbp)
rbps_of_interest
# rbp_of_interest <- "PCBP1"
tib <- tibble(rbp_of_interest=character(), other_rbp=character(), n=numeric(), complextitle=character(), r_vs_other_rbp=numeric(), r_vs_complex=numeric())
for (rbp_of_interest in rbps_of_interest$rbp) {
  
  n <- rbps_of_interest %>% filter(rbp==rbp_of_interest) %>% pull(n)
  
  print(glue(" >> {rbp_of_interest}"))

  # m$r["U2AF1", "U2AF2"]
  # m$r["U2AF2", "U2AF1"]
  # m$r["PCBP1",]
  # r %>% filter(rbp1 == "PCBP1")
  # r %>% filter(rbp2 == "PCBP1")
  
  rbp_pair <- rbp_pairs[[rbp_of_interest]]
  for (rbp_pair in rbp_pairs) {

    # str(rbp_pair)
    my_rbp1 <- rbp_pair[1]
    my_rbp2 <- rbp_pair[2]
    my_complextitle <- glue("{my_rbp1} & {my_rbp2}")
    # my_complextitle

    if (my_rbp1 == rbp_of_interest) {
      other_rbp <- my_rbp2
    } else if (my_rbp2 == rbp_of_interest) {
      other_rbp <- my_rbp1
    }
    
    if ((my_rbp1 == rbp_of_interest) | (my_rbp2 == rbp_of_interest)) {
      if (r %>% filter((rbp1 == rbp_of_interest) & (rbp2 == my_complextitle)) %>% nrow > 0) {
        print(glue("   >> {my_complextitle}"))
        # r %>% filter((rbp1 == my_rbp1) & (rbp2 == my_rbp2))
        
        # # # Get correlations between the individual proteins (all their sites):
        # tib %<>% add_row(rbp_of_interest = rbp_of_interest, other_rbp = other_rbp, n = n, complextitle = my_complextitle, )
        
        # Get correlations between an individual protein (all its sites) and its complexes (only when bound in proximity):
        tib %<>% add_row(rbp_of_interest = rbp_of_interest, other_rbp = other_rbp, n = n, complextitle = my_complextitle, r_vs_other_rbp = r %>% filter((rbp1 == rbp_of_interest) & (rbp2 == other_rbp)) %>% pull(r), r_vs_complex = r %>% filter((rbp1 == rbp_of_interest) & (rbp2 == my_complextitle)) %>% pull(r))
      }
    }
  }
}
tib
tib %>% filter((rbp_of_interest == "PCBP1") & (complextitle == "PCBP1 & RBFOX2"))
tib %>% filter((rbp_of_interest == "RBFOX2") & (complextitle == "PCBP1 & RBFOX2"))

# Use correlations between an individual protein (all its sites) and its complexes (only when bound in proximity):
tib$r <- tib$r_vs_complex

# v19.1, ≥3 horizontal:
tib %>% 
  filter(n >= 3) %>%
  # ggplot(aes(x=fct_reorder(rbp_of_interest, desc(n)), y=r, colour=rbp_of_interest)) + 
  ggplot(aes(x=fct_reorder(rbp_of_interest, r, .fun=mean), y=r, colour=other_rbp)) + 
  geom_beeswarm(groupOnX = T) +
  # geom_boxplot() +
  # geom_violin() +
  # geom_point(position="jitter") +
  geom_text_repel(aes(label=other_rbp), direction="x", nudge_x=0.3, size=2) +
  scale_y_reverse(limits=c(1,0), breaks=pretty_breaks(5)) +
  scale_colour_viridis_d() +
  guides(color = FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(NULL) +
  # xlab("RBP of interest") +
  ylab("pre-mRNA region preference similarity") +
  coord_flip() +
  ggsave(paste0("output-plot-metagene-v19.1-",clip_cobinding,"-",type,"-",style,".pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")

# v19.2, ≥3 vertical:
tib %>% 
  filter(n >= 3) %>%
  # ggplot(aes(x=fct_reorder(rbp_of_interest, desc(n)), y=r, colour=rbp_of_interest)) + 
  ggplot(aes(x=fct_reorder(rbp_of_interest, desc(r), .fun=mean), y=r, colour=other_rbp)) + 
  geom_beeswarm(groupOnX = T) +
  # geom_boxplot() +
  # geom_violin() +
  # geom_point(position="jitter") +
  geom_text_repel(aes(label=other_rbp), direction="both", size=2) +
  # scale_y_reverse() +
  scale_y_continuous(limits=c(0,1), breaks=pretty_breaks(5)) +
  scale_colour_viridis_d() +
  guides(color = FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  xlab(NULL) +
  # xlab("RBP of interest") +
  ylab("pre-mRNA region preference similarity") +
  # coord_flip() +
  ggsave(paste0("output-plot-metagene-v19.2-",clip_cobinding,"-",type,"-",style,".pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")

# v19.3, ≥2 horizontal:
tib %>% 
  filter(n >= 2) %>%
  # ggplot(aes(x=fct_reorder(rbp_of_interest, desc(n)), y=r, colour=rbp_of_interest)) + 
  ggplot(aes(x=fct_reorder(rbp_of_interest, r, .fun=mean), y=r, colour=other_rbp)) + 
  geom_beeswarm(groupOnX = T) +
  # geom_boxplot() +
  # geom_violin() +
  # geom_point(position="jitter") +
  geom_text_repel(aes(label=other_rbp), direction="x", nudge_x=0.4, size=2, force=10) +
  scale_y_reverse(limits=c(1,0), breaks=pretty_breaks(5)) +
  scale_colour_viridis_d() +
  guides(color = FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(NULL) +
  # xlab("RBP of interest") +
  ylab("pre-mRNA region preference similarity") +
  coord_flip() +
  ggsave(paste0("output-plot-metagene-v19.3-",clip_cobinding,"-",type,"-",style,".pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")

# v19.4, ≥2 vertical:
tib %>% 
  filter(n >= 2) %>%
  # ggplot(aes(x=fct_reorder(rbp_of_interest, desc(n)), y=r, colour=rbp_of_interest)) + 
  ggplot(aes(x=fct_reorder(rbp_of_interest, desc(r), .fun=mean), y=r, colour=other_rbp)) + 
  geom_beeswarm(groupOnX = T) +
  # geom_boxplot() +
  # geom_violin() +
  # geom_point(position="jitter") +
  geom_text_repel(aes(label=other_rbp), direction="both", size=2) +
  # scale_y_reverse() +
  scale_y_continuous(limits=c(0,1.1), breaks=pretty_breaks(5)) +
  scale_colour_viridis_d() +
  guides(color = FALSE) +
  theme_minimal() + 
  # theme(axis.text.x = element_text(angle = 30, hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  xlab(NULL) +
  # xlab("RBP of interest") +
  ylab("pre-mRNA region preference similarity") +
  # coord_flip() +
  ggsave(paste0("output-plot-metagene-v19.4-",clip_cobinding,"-",type,"-",style,".pdf"), device=cairo_pdf, width=91.5, height=91.5*1.1, units="mm")
    
# v19.5, ≥3 vertical smaller:
tib %>% 
  filter(n >= 3) %>%
  # ggplot(aes(x=fct_reorder(rbp_of_interest, desc(n)), y=r, colour=rbp_of_interest)) + 
  ggplot(aes(x=fct_reorder(rbp_of_interest, desc(r), .fun=mean), y=r, colour=other_rbp)) + 
  geom_beeswarm(groupOnX = T) +
  # geom_boxplot() +
  # geom_violin() +
  # geom_point(position="jitter") +
  geom_text_repel(aes(label=other_rbp), direction="both", size=2) +
  # scale_y_reverse() +
  scale_y_continuous(limits=c(0,1.1), breaks=pretty_breaks(5)) +
  scale_colour_viridis_d() +
  guides(color = FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  xlab(NULL) +
  # xlab("RBP of interest") +
  ylab("Region profile similarity") +
  # coord_flip() +
  # theme(axis.title.y=element_text(size=9)) +
  ggsave(paste0("output-plot-metagene-v19.5-",clip_cobinding,"-",type,"-",style,".pdf"), device=cairo_pdf, width=91.5, height=91.5/2, units="mm")
 
# v19.6, ≥2 vertical vertical labels:
tib %>% 
  filter(n >= 2) %>%
  # ggplot(aes(x=fct_reorder(rbp_of_interest, desc(n)), y=r, colour=rbp_of_interest)) + 
  ggplot(aes(x=fct_reorder(rbp_of_interest, desc(r), .fun=mean), y=r, colour=other_rbp)) + 
  geom_beeswarm(groupOnX = T) +
  # geom_boxplot() +
  # geom_violin() +
  # geom_point(position="jitter") +
  geom_text_repel(aes(label=other_rbp), direction="both", size=2) +
  # scale_y_reverse() +
  scale_y_continuous(limits=c(0,1.1), breaks=pretty_breaks(5)) +
  scale_colour_viridis_d() +
  guides(color = FALSE) +
  theme_minimal() + 
  # theme(axis.text.x = element_text(angle = 30, hjust = 0.5)) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  xlab(NULL) +
  # xlab("RBP of interest") +
  ylab("pre-mRNA region preference similarity") +
  # coord_flip() +
  ggsave(paste0("output-plot-metagene-v19.6-",clip_cobinding,"-",type,"-",style,".pdf"), device=cairo_pdf, width=91.5, height=91.5, units="mm")

