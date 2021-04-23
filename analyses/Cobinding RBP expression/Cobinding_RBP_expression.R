rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

source("~/Documents/R/mysqltc.R")
library(tidyverse)
library(reshape2)
library(glue)
library(broom)
library(scales)
library(magrittr)

# # Correlation plot packages
# library(corrplot)
# library(Hmisc)
# library(ggcorrplot)
# 
# # Plotting packages
# library(ggbeeswarm)
# library(ggrepel)

# "SELECT a.enstv, t.symbolv, t.chr, t.start, t.stop, t.strand FROM encode_abundance a, gencode_gff3_transcript t WHERE a.ensgv='$ensgv' AND a.enstv=t.enstv AND t.species='$species' AND t.biotype='protein_coding' GROUP BY a.enstv ORDER BY AVG(a.tpm) DESC LIMIT 1"

# Strictly filtered for cobinding (resampling p-value and resampling Wilcoxon p-value significant in one orientation at least)
# Real RBP pairs (SUM_IS screen hits, and with eCLIP data of course) (21: the "strict" Cytoscape network, Figure 5a)
qtmp <- Query("SELECT protein_a AS symbol1, protein_b AS symbol2 FROM rbpome_final WHERE encode_eclip_data_a=1 AND encode_eclip_data_b=1 AND (resampling_p_value_a_vs_b IS NOT NULL AND resampling_p_value_b_vs_a IS NOT NULL AND resampling_wilcoxon_p_value_a_vs_b IS NOT NULL AND resampling_wilcoxon_p_value_b_vs_a IS NOT NULL) AND ((resampling_p_value_a_vs_b<0.05 AND resampling_wilcoxon_p_value_a_vs_b<0.05) OR (resampling_p_value_b_vs_a<0.05 AND resampling_wilcoxon_p_value_b_vs_a<0.05)) GROUP BY protein_a, protein_b ORDER BY protein_a, protein_b")
# qtmp %>% select(symbol1, symbol2) %>% unique %>% transmute(pair = glue("{symbol1} {symbol2}")) -> qtmppairs
qtmp %>% select(symbol1, symbol2) %>% unique -> qtmppairs
print(qtmppairs, n=50)

# New expression query:
# q <- Query("SELECT r.protein_a AS symbol1, r.protein_b AS symbol2, a1.celltype AS celltype1, MAX(a1.tpm) AS maxtpm1, a2.celltype AS celltype2, MAX(a2.tpm) AS maxtpm2 FROM rbpome_final r, encode_abundance a1, encode_abundance a2, gencode_gff3_gene g1, gencode_gff3_gene g2 WHERE g1.species='HUMAN' AND g2.species='HUMAN' AND r.protein_a=g1.symbol AND r.protein_b=g2.symbol AND g1.ensgv=a1.ensgv AND g2.ensgv=a2.ensgv AND r.encode_eclip_data_a=1 AND r.encode_eclip_data_b=1 AND (r.resampling_p_value_a_vs_b IS NOT NULL AND r.resampling_p_value_b_vs_a IS NOT NULL AND r.resampling_wilcoxon_p_value_a_vs_b IS NOT NULL AND r.resampling_wilcoxon_p_value_b_vs_a IS NOT NULL) AND ((r.resampling_p_value_a_vs_b<0.05 AND r.resampling_wilcoxon_p_value_a_vs_b<0.05) OR (r.resampling_p_value_b_vs_a<0.05 AND r.resampling_wilcoxon_p_value_b_vs_a<0.05)) GROUP BY r.protein_a, r.protein_b, a1.celltype, a2.celltype ORDER BY r.protein_a, r.protein_b")
# q %<>% filter(celltype1==celltype2) %>% mutate(celltype = celltype1) %>% select(-celltype1, -celltype2)
q <- Query("SELECT r.protein_a AS symbol1, r.protein_b AS symbol2, g1.ensgv AS ensgv1, g2.ensgv AS ensgv2, IF(a1.celltype IS NOT NULL, a1.celltype, a2.celltype) AS celltype, MAX(a1.tpm) AS maxtpm1, MAX(a2.tpm) AS maxtpm2 FROM rbpome_final r
LEFT OUTER JOIN gencode_gff3_gene g1 ON r.protein_a=g1.symbol AND g1.species='HUMAN'
LEFT OUTER JOIN gencode_gff3_gene g2 ON r.protein_b=g2.symbol AND g2.species='HUMAN'
LEFT OUTER JOIN encode_abundance a1 ON g1.ensgv=a1.ensgv
LEFT OUTER JOIN encode_abundance a2 ON g2.ensgv=a2.ensgv AND (a2.celltype=a1.celltype OR a1.celltype IS NULL OR a2.celltype IS NULL)
WHERE r.encode_eclip_data_a=1 AND r.encode_eclip_data_b=1 AND (r.resampling_p_value_a_vs_b IS NOT NULL AND r.resampling_p_value_b_vs_a IS NOT NULL AND r.resampling_wilcoxon_p_value_a_vs_b IS NOT NULL AND r.resampling_wilcoxon_p_value_b_vs_a IS NOT NULL) AND ((r.resampling_p_value_a_vs_b<0.05 AND r.resampling_wilcoxon_p_value_a_vs_b<0.05) OR (r.resampling_p_value_b_vs_a<0.05 AND r.resampling_wilcoxon_p_value_b_vs_a<0.05)) GROUP BY r.protein_a, r.protein_b, g1.ensgv, g2.ensgv, a1.celltype, a2.celltype ORDER BY r.protein_a, r.protein_b")
q
q %>% filter((!is.na(maxtpm1)) & (!is.na(maxtpm2))) %>% select(symbol1, symbol2) %>% unique -> qpairs
qpairs

# Cobinding pairs that don't seem to be expressed in the same cell type:
setdiff(qtmppairs, qpairs)

# Get eCLIP cell types
qe <- Query("SELECT DISTINCT symbol, celltype FROM clip_raw_gene WHERE type='eclip_encode_12'")
# qe <- Query("SELECT symbol, GROUP_CONCAT(DISTINCT celltype ORDER BY celltype) AS eclip_celltypes FROM clip_raw_gene WHERE type='eclip_encode_12' GROUP BY symbol ORDER BY symbol")
qe

# Check whether the RBPs with no expression in a given celltype do have eCLIP data for it, though
q
qe
q %<>% left_join(qe %>% mutate(eclip1="Yes"), by=c("symbol1" = "symbol", "celltype")) %>% left_join(qe %>% mutate(eclip2="Yes"), by=c("symbol2" = "symbol", "celltype"))
write_tsv(q, "Table S8.tsv")
q
q %>% filter((is.na(maxtpm1) & is.na(eclip1)) | (is.na(maxtpm2) & is.na(eclip2)))
# >> No: Wherever we have a missing expression value in a given cell type, there is eCLIP data for that RBP, proving that it was in fact expressed. The expression data is unreliable.



# How many pairs have direct BioGRID evidence?
qtmppairs
# qb <- Query("SELECT DISTINCT gene1 AS symbol1, gene2 AS symbol2 FROM biogrid WHERE system='Co-crystal Structure'")
# qb <- Query("SELECT DISTINCT gene1 AS symbol1, gene2 AS symbol2 FROM biogrid WHERE system='Co-crystal Structure' OR system='Reconstituted Complex'")
qb <- Query("SELECT DISTINCT gene1 AS symbol1, gene2 AS symbol2, GROUP_CONCAT(DISTINCT system ORDER BY system) AS biogrid FROM biogrid WHERE direct=1 GROUP BY gene1, gene2")
qb
qbflip <- 
qb %<>% bind_rows(qb %>% transmute(newsymbol1=symbol2, newsymbol2=symbol1) %>% transmute(symbol1=newsymbol1, symbol2=newsymbol2)) %>% unique
qb
qtmppairs %>% left_join(qb) %>% filter(!is.na(biogrid))
# >> Only 3 pairs ("reconstituted complex").
