rm(list=setdiff(ls(), c("superreservedcon", "superreserveddrv")))
options(nwarnings = 10000)

source("~/Documents/R/mysqltc.R")
library(tidyverse)
library(broom)
library(magrittr)
library(scales)
library(glue)

setwd("~/Documents/Projects/RBPome/Comparison with other screens")

q <- tibble(Query("SELECT * FROM rbpome_final WHERE uniqueabsumis_max>=7.1"))
q

# q2
# tidy(table(q2))
# q2 %>%


# BioGRID/HIPPIE/HuRI
confirmed_name <- "Confirmed"
novel_name <- "Novel"
q %>% 
  rename(biogrid_alldirect = biogrid_all_direct_evidence_y2h_reconstituted_complex_structure) %>% 
  select(biogrid_alldirect, biogrid_all, hippie, huri) %>% 
  pivot_longer(cols=c(biogrid_alldirect, biogrid_all, hippie, huri)) %>%
  # group_by(name) %>%
  # count(name, value) %>%
  mutate(name = recode(name, biogrid_alldirect="BioGRID direct", biogrid_all="BioGRID", hippie="HIPPIE", huri="HuRI")) %>%
  # mutate(name=factor(name), value=factor(value)) %>%
  mutate(value = recode(value, '0'=novel_name, '1'=confirmed_name)) %>%
  # ggplot(aes(x=fct_rev(name), fill=fct_rev(value))) + 
  ggplot(aes(x=fct_rev(fct_relevel(name, "BioGRID", "HIPPIE")), fill=fct_rev(value))) + 
  geom_bar() + 
  scale_fill_manual(NULL, values = setNames(c("#1e3d59", "#ff7f00"), c(confirmed_name, novel_name))) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme_minimal() + 
  labs(subtitle="Expansion on existing resources") +
  xlab(NULL) + 
  ylab("Number of RBP-RBP interactions") + 
  coord_flip() + 
  # guides(fill=F) +
  # guides(fill = guide_legend(order = 1)) +
  guides(fill = guide_legend(reverse=T)) +
  theme(legend.position="top", legend.justification="left") +
ggsave("output-comparison-interactions.pdf", width=91.5, height=91.5*2/3, units="mm")

# Get number that's completely new
# -- known
# SELECT * FROM rbpome_final WHERE uniqueabsumis_max>=7.1 AND (biogrid_all_direct_evidence_y2h_reconstituted_complex_structure=1 OR biogrid_all=1 OR hippie=1 OR huri=1) LIMIT 1000000;
# SELECT DISTINCT protein_a_unique_pairs, protein_b_unique_pairs FROM rbpome_final WHERE uniqueabsumis_max>=7.1 AND (biogrid_all_direct_evidence_y2h_reconstituted_complex_structure=1 OR biogrid_all=1 OR hippie=1 OR huri=1) LIMIT 1000000;
# >> 422 unique pairs known (the DISTINCT is unnecessary)
# -- new
# SELECT * FROM rbpome_final WHERE uniqueabsumis_max>=7.1 AND biogrid_all_direct_evidence_y2h_reconstituted_complex_structure=0 AND biogrid_all=0 AND hippie=0 AND huri=0 LIMIT 1000000;
# SELECT DISTINCT protein_a_unique_pairs, protein_b_unique_pairs FROM rbpome_final WHERE uniqueabsumis_max>=7.1 AND biogrid_all_direct_evidence_y2h_reconstituted_complex_structure=0 AND biogrid_all=0 AND hippie=0 AND huri=0 LIMIT 1000000;
# >> 1994 unique pairs new (the DISTINCT is unnecessary)

  


# eCLIP
# q %>% 
#   mutate(eclip = paste0(encode_eclip_data_a,'/',encode_eclip_data_b)) %>%
#   rename(biogrid_alldirect = biogrid_all_direct_evidence_y2h_reconstituted_complex_structure) %>% 
#   select(biogrid_alldirect, biogrid_all, hippie, huri, eclip) %>% 
#   mutate(eclip = recode_factor(eclip, "0/0"="neither", "1/0"="bait only", "0/1"="prey only", "1/1"="both"))
# rbpome_name <- "RBPome"
rbpome_name <- "Present study"
q %>% 
  mutate(rbpome = 1, eclip = paste0(encode_eclip_data_a,'/',encode_eclip_data_b)) %>%
  rename(biogrid_alldirect = biogrid_all_direct_evidence_y2h_reconstituted_complex_structure) %>% 
  select(rbpome, biogrid_alldirect, biogrid_all, hippie, huri, eclip) %>% 
  mutate(eclip = recode_factor(eclip, "0/0"="neither", "1/0"="bait only", "0/1"="prey only", "1/1"="both")) %>%
  pivot_longer(cols=c(rbpome, biogrid_all, biogrid_alldirect, hippie, huri)) %>%
  mutate(name = recode(name, biogrid_alldirect="BioGRID direct", biogrid_all="BioGRID", hippie="HIPPIE", huri="HuRI", rbpome=rbpome_name)) %>%
  group_by(name, value, eclip) %>%
  filter(value==1) %>%
  tally() %>%
  filter(eclip=="both") %>%
ggplot(aes(x=fct_reorder(name, n), y=n, fill=(name==rbpome_name))) + 
  geom_col() +
  scale_fill_manual(NULL, values = c("FALSE"="#1e3d59", "TRUE"="#ff7f00")) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme_minimal() + 
  labs(subtitle="Protein-RNA interaction data (eCLIP)") +
  xlab(NULL) + 
  ylab("Interacting RBP pairs with eCLIP data") +
  # ylab("Interactions with protein-RNA interaction data") +
  # ylab(NULL) + 
  coord_flip() +
  guides(fill=F) +
  # guides(fill = guide_legend(order = 1)) +
  # guides(fill = guide_legend(reverse=T)) +
  theme(legend.position="top") +
ggsave("output-comparison-eclip.pdf", width=91.5, height=91.5*2/3, units="mm")

ours <- q %>% 
  mutate(rbpome = 1, eclip = paste0(encode_eclip_data_a,'/',encode_eclip_data_b)) %>%
  rename(symbol1 = protein_a_unique_pairs, symbol2 = protein_b_unique_pairs, biogrid_alldirect = biogrid_all_direct_evidence_y2h_reconstituted_complex_structure) %>% 
  mutate(pair = glue("{symbol1}|{symbol2}")) %>%
  select(pair, rbpome, biogrid_alldirect, biogrid_all, hippie, huri, eclip) %>% 
  mutate(eclip = recode_factor(eclip, "0/0"="neither", "1/0"="bait only", "0/1"="prey only", "1/1"="both")) %>%
  pivot_longer(cols=c(rbpome, biogrid_all, biogrid_alldirect, hippie, huri)) %>%
  # mutate(name = recode(name, biogrid_alldirect="BioGRID direct", biogrid_all="BioGRID", hippie="HIPPIE", huri="HuRI", rbpome=rbpome_name)) %>%
  mutate(name = recode(name, biogrid_alldirect="Other", biogrid_all="Other", hippie="Other", huri="Other", rbpome=rbpome_name)) %>%
  group_by(name, value, eclip) %>%
  filter(value==1) %>%
  filter(eclip=="both") %>%
  unique %>%
  ungroup %>%
  select(pair, name) %>%
  filter(name==rbpome_name) %>%
  select(pair)
others <- q %>% 
  mutate(rbpome = 1, eclip = paste0(encode_eclip_data_a,'/',encode_eclip_data_b)) %>%
  rename(symbol1 = protein_a_unique_pairs, symbol2 = protein_b_unique_pairs, biogrid_alldirect = biogrid_all_direct_evidence_y2h_reconstituted_complex_structure) %>% 
  mutate(pair = glue("{symbol1}|{symbol2}")) %>%
  select(pair, rbpome, biogrid_alldirect, biogrid_all, hippie, huri, eclip) %>% 
  mutate(eclip = recode_factor(eclip, "0/0"="neither", "1/0"="bait only", "0/1"="prey only", "1/1"="both")) %>%
  pivot_longer(cols=c(rbpome, biogrid_all, biogrid_alldirect, hippie, huri)) %>%
  # mutate(name = recode(name, biogrid_alldirect="BioGRID direct", biogrid_all="BioGRID", hippie="HIPPIE", huri="HuRI", rbpome=rbpome_name)) %>%
  mutate(name = recode(name, biogrid_alldirect="Other", biogrid_all="Other", hippie="Other", huri="Other", rbpome=rbpome_name)) %>%
  group_by(name, value, eclip) %>%
  filter(value==1) %>%
  filter(eclip=="both") %>%
  unique %>%
  ungroup %>%
  select(pair, name) %>%
  filter(name=="Other") %>%
  select(pair)

intersect(ours, others)
ours
others

nrow(ours)
nrow(others)

nrow(ours) / nrow(others)
# Our screen more than doubles the number of eCLIP-eCLIP interactions, from 28 to 71 (2.5x)!




# Get number of proteins screened for BioGRID/HIPPIE/HuRI

# rbpome
n <- q %>%
  select(protein_a_unique_pairs, protein_b_unique_pairs) %>% 
  pivot_longer(cols=c(protein_a_unique_pairs, protein_b_unique_pairs)) %>%
  mutate(name="rbpome") %>%
  select(name, value) %>%
  unique()
n

# # biogrid_all
# n %<>% add_row(tibble(Query("SELECT DISTINCT 'biogrid_all' AS 'name', gene1 AS 'value' FROM biogrid")))
# n %<>% add_row(tibble(Query("SELECT DISTINCT 'biogrid_all' AS 'name', gene2 AS 'value' FROM biogrid")))
# n %<>% unique()
# n
# 
# # biogrid_alldirect
# n %<>% add_row(tibble(Query("SELECT DISTINCT 'biogrid_all' AS 'name', gene1 AS 'value' FROM biogrid WHERE direct=1")))
# n %<>% add_row(tibble(Query("SELECT DISTINCT 'biogrid_all' AS 'name', gene2 AS 'value' FROM biogrid WHERE direct=1")))
# n %<>% unique()

# huri
# # p.bravenew=1 OR p.gerstberger=1 OR p.sonar=1 OR p.signature=1
# n %<>% add_row(tibble(Query("SELECT DISTINCT 'huri' AS 'name', h.symbol1 AS 'value' FROM huri h, catrapid_proteins p WHERE p.species='human' AND h.symbol1=p.symbol AND (p.bravenew=1 OR p.gerstberger=1 OR p.sonar=1 OR p.signature=1)")))
# n %<>% add_row(tibble(Query("SELECT DISTINCT 'huri' AS 'name', h.symbol2 AS 'value' FROM huri h, catrapid_proteins p WHERE p.species='human' AND h.symbol2=p.symbol AND (p.bravenew=1 OR p.gerstberger=1 OR p.sonar=1 OR p.signature=1)")))
# p.bravenew=1 OR p.gerstberger=1
n %<>% add_row(tibble(Query("SELECT DISTINCT 'huri' AS 'name', h.symbol1 AS 'value' FROM huri h, catrapid_proteins p WHERE p.species='human' AND h.symbol1=p.symbol AND (p.bravenew=1 OR p.gerstberger=1)")))
n %<>% add_row(tibble(Query("SELECT DISTINCT 'huri' AS 'name', h.symbol2 AS 'value' FROM huri h, catrapid_proteins p WHERE p.species='human' AND h.symbol2=p.symbol AND (p.bravenew=1 OR p.gerstberger=1)")))
n %<>% unique()
n

n %>% ggplot(aes(x=fct_rev(name), fill=name)) +
  geom_bar() +
  scale_fill_manual(NULL, values = c("#1e3d59", "#ff7f00")) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme_minimal() + 
  labs(subtitle="Interaction sampling") +
  xlab(NULL) + 
  ylab("Number of RBPs sampled as bait or prey") +
  # ylab("Interactions with protein-RNA interaction data") +
  # ylab(NULL) + 
  coord_flip() +
  guides(fill=F) +
  # guides(fill = guide_legend(order = 1)) +
  # guides(fill = guide_legend(reverse=T)) +
  theme(legend.position="top") +
  ggsave("output-comparison-input.pdf", width=91.5, height=91.5/2, units="mm")




