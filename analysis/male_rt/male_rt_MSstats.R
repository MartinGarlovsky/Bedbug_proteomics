library(tidyverse)
library(MSstats)
library(ComplexHeatmap)
library(UpSetR)
library(ggVennDiagram)
library(eulerr)
library(cowplot)

library(conflicted)
select <- dplyr::select
filter <- dplyr::filter
count <- dplyr::count
as.factor <- base::as.factor

# met brewer palette
Hiroshige.6 <- MetBrewer::met.brewer('Hiroshige')[c(1, 2, 4, 7, 9, 10)]
Hiroshige.4 <- MetBrewer::met.brewer('Hiroshige')[c(1, 4, 7, 9)]

# Load data ####
# load combined signal peptide results from Phobius/sigP for the entire bedbug proteome
signal_peps <- read_rds("output/cimex_signal_peps_combined.rds")

# DAVID uniprot ID to entrez_geneID
uniprot2entrez <- read.delim("../Bedbug_RNAseq/DAVID2entrez_geneID.txt") %>%
  mutate(gene_id = paste0("LOC", To)) %>%
  dplyr::rename(accession = From)

# flybase gene ids to uniprot ids
flybase2uniprot <- read_rds("output/flybase2uniprot.rds")

# load mass spec data
ms_data <- readxl::read_excel("data/male_rt/LFQ_v2.0.xlsx") %>% janitor::clean_names()

# # total PSMs
sum(ms_data$number_ps_ms)

# the dataset contains endosymbionts and some contaminants
# pull species
organism_protein <- gsub(" OX.*", "", x = gsub(".*OS=", "", ms_data$description))
#rle(organism_protein)
data.frame(protein = ms_data$accession,
           organism_protein) %>%
  dplyr::count(organism_protein)

# endosymbiont proteins in the dataset
endosymbiont_prots <- ms_data[grepl("endosymbiont", ms_data$description), ] %>%
  mutate(organism = str_replace(description, pattern = ".*OS=([^-]*) OX.*", replacement = "\\1"))

# remove endosymbiont proteins
no_endos <- ms_data %>% filter(!accession %in% endosymbiont_prots$accession)

# only include proteins with 2 or more unique peptides
# exclude
bbd <- no_endos %>%
  dplyr::select(accession, number_unique_peptides,
                abundance_f1_sample_ags:abundance_f17_sample_testis) %>%
  filter(number_unique_peptides >= 2)

# number/percent proteins with 2+ unique peptides
n_distinct(bbd$accession)/n_distinct(no_endos$accession) * 100

# # PSMs
sum(no_endos$number_ps_ms[no_endos$accession %in% bbd$accession])

# clean headers
colnames(bbd)[-c(1:2)] <- gsub(".*_", "", x = colnames(bbd)[-c(1:2)])
colnames(bbd)[-c(1:2)] <- paste(colnames(bbd)[-c(1:2)],
                                unlist(lapply(rle(colnames(bbd))$lengths[-c(1:2)],
                                              function(x) seq(from = 1, to = x, by = 1))), sep = "_")
colnames(bbd)[5:7] <- c("ags_3.1", "ags_3.2", "ags_4")
# replace ags name with sfc (seminal fluid containers)
colnames(bbd) <- gsub("ags", "sfc", x = colnames(bbd))

# Seminal fluid containers sample 3 was run twice on the mass spectrometer.
bbd %>% select(accession, contains("sfc")) %>%
  pivot_longer(cols = 2:6) %>%
  ggplot(aes(x = log10(value + 1))) +
  geom_histogram() +
  facet_wrap(~ name)

# technical replicates show high correlation in abundance
plot(log10(bbd$sfc_3.1), log10(bbd$sfc_3.2))
cor.test(log10(bbd$sfc_3.1), log10(bbd$sfc_3.2), method = "p")

# technical replicates show >85% overlap in proteins identified
upset(fromList(list(sfc_3.1 = bbd$accession[is.na(bbd$sfc_3.1) == FALSE],
                    sfc_3.2 = bbd$accession[is.na(bbd$sfc_3.2) == FALSE])))

# > Create MSstats input ####
bbd_reduced <- bbd %>%
  rowwise() %>%
  # take the average of the 2 reps for sfc_3
  mutate(sfc_3 = mean(c(sfc_3.1, sfc_3.2))) %>%
  #mutate(sfc_3 = sfc_3.2) %>%
  select(-c(sfc_3.1, sfc_3.2)) %>%
  relocate(sfc_3, .before = sfc_4) %>%
  select(-2) %>% dplyr::rename(ProteinName = accession) %>%
  pivot_longer(cols = 2:17) %>%
  group_by(name) %>% mutate(Run = cur_group_id()) %>% ungroup() %>%
  mutate(PeptideSequence = ProteinName,
         PrecursorCharge = NA,
         FragmentIon = NA,
         ProductCharge = NA,
         IsotopeLabelType = "L",
         Condition = str_sub(name, 1, end = -3),
         BioReplicate = name,
         Intensity = value) %>%
  select(-name, -value) %>% relocate(Run, .after = BioReplicate) %>% mutate(across(2:9, as.factor)) %>%
  data.frame()

# correlation between replicates within each sample type
cord <- bbd %>%
  rowwise() %>%
  # take the average of the 2 reps for ags_3
  mutate(sfc_3 = mean(c(sfc_3.1, sfc_3.2))) %>%
  select(-c(sfc_3.1, sfc_3.2)) %>%
  relocate(sfc_3, .before = sfc_4) %>%
  select(-2)

# all replicates show >0.95 correlation
cor(cord[, -1], use = "complete.obs", method = "p") %>% reshape2::melt() %>%
  mutate(t1 = gsub("_.", "", x = Var1),
         t2 = gsub("_.", "", x = Var2)) %>%
  filter(t1 == t2, Var1 != Var2) %>% summary() #group_by(t1) %>% summarise(mn = mean(value))

# > MSstats ####
if(!file.exists("output/male_rt_processed.rds")){

  male_rt_processed <- dataProcess(raw = bbd_reduced,
                                   normalization = 'quantile',
                                   summaryMethod = 'TMP',
                                   censoredInt = "NA",
                                   MBimpute = TRUE,
                                   maxQuantileforCensored = 0.999,
                                   use_log_file = FALSE)

  saveRDS(male_rt_processed, "output/male_rt_processed.rds")
} else {
  male_rt_processed <- readRDS('output/male_rt_processed.rds')
}

# some proteins are dropped at this processing stage by MSstats default filtering
n_distinct(male_rt_processed$ProteinLevelData$Protein)/n_distinct(bbd$accession) * 100

# Normalised data
norm_data <- quantification(male_rt_processed)
colnames(norm_data)[-1] <- gsub("^[A-Za-z1-9]+_", "", x = colnames(norm_data)[-1])

# number of replicates/protein
protein_counts <- norm_data %>%
  pivot_longer(cols = 2:17) %>%
  separate(name, into = c("tissue", "repl")) %>%
  drop_na(value) %>%
  group_by(tissue, Protein) %>% dplyr::count() %>%
  pivot_wider(names_from = tissue, values_from = n) %>%
  rowwise() %>%
  mutate(sumVar = sum(c_across(mesa:testis), na.rm = TRUE))

# find proteins only identified in one tissue/sperm and in at least 2 replicates
sfc_only <- protein_counts %>% filter(sumVar == sfc, sumVar >= 2)
mes_only <- protein_counts %>% filter(sumVar == mesa, sumVar >= 2)
spe_only <- protein_counts %>% filter(sumVar == sperm, sumVar >= 2)
tes_only <- protein_counts %>% filter(sumVar == testis, sumVar >= 2)

# PCA
# extract the 500 most variable proteins for PCA
most_var <- data.frame(Protein = norm_data[, 1], vars = apply(norm_data[, -1], 1, var)) %>%
  slice_max(vars, n = 500)

most_var_prot <- norm_data %>% filter(Protein %in% most_var$Protein) %>% dplyr::select(-Protein)

pca <- prcomp(t(as.matrix(most_var_prot)), center = TRUE, scale. = TRUE)
#summary(pca) # PCs 1 and 2 explain >80% of variance in the data
#pca$rotation

PCA_dat <- as.data.frame(pca$x)[, 1:3] %>%
  rownames_to_column() %>%
  separate(rowname, into = c('tissue', "bio"), sep = '_', remove = FALSE)

# Plot for figure
PCA_dat %>%
  ggplot(aes(x = PC1, y = PC2, colour = tissue)) +
  geom_point(size = 3, alpha = .7) +
  labs(x = paste0('PC1 (', (100*round(summary(pca)$importance[2, 1], 3)), '%)'),
       y = paste0('PC2 (', (100*round(summary(pca)$importance[2, 2], 3)), '%)')) +
  # lims(x = c(-35, 45), y = c(-45, 45)) +
  scale_colour_manual(values = rev(Hiroshige.4), name = "Tissue:") +
  #scale_shape_manual(values = c(15:18), name = "Replicate:") +
  theme_bw() +
  theme(legend.position = "",
        legend.text = element_text(size = 12, hjust = 0),
        legend.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  #ggsave('plots/male_rt/PCA_12.pdf', height = 2, width = 2.4, dpi = 600, useDingbats = FALSE) +
  NULL

# pairwise plots of PCs 1, 2, 3
rbind(as.matrix(PCA_dat[, c(4, 5)]),
      as.matrix(PCA_dat[, c(4, 6)]),
      as.matrix(PCA_dat[, c(5, 6)])) %>%
  bind_cols(tissue = rep(PCA_dat$tissue, 3),
            bio = rep(PCA_dat$bio, 3),
            pc = rep(c(paste0('PC1 (',100*round(summary(pca)$importance[2, 1], 3), '%) vs. PC2 (',
                              100*round(summary(pca)$importance[2, 2], 3), '%)'),
                       paste0('PC1 (',100*round(summary(pca)$importance[2, 1], 3), '%) vs. PC3 (',
                              100*round(summary(pca)$importance[2, 3], 3), '%)'),
                       paste0('PC2 (',100*round(summary(pca)$importance[2, 2], 3), '%) vs. PC3 (',
                              100*round(summary(pca)$importance[2, 3], 3), '%)')),
                     each = 16)) %>%
  ggplot(aes(x = PC1, y = PC2, colour = tissue, alpha = .5)) +
  geom_point(size = 8, alpha = .7) +
  scale_colour_manual(values = rev(Hiroshige.4), name = "Tissue:") +
  scale_shape_manual(values = c(15:18)) +
  facet_wrap(~pc) +
  theme_bw() +
  theme(legend.position = "",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, hjust = 0),
        legend.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_blank(),
        strip.text = element_text(size = 15)) +
  NULL


### Make contrasts ###
sampleInfo <- data.frame(Sample = levels(male_rt_processed$ProteinLevelData$SUBJECT)) %>%
  separate(Sample, into = c('tissue', 'REP'), remove = FALSE)

model1 <- model.matrix(~ 0 + tissue, data = sampleInfo)
colnames(model1) <- gsub('tissue', '', colnames(model1))
rownames(model1) <- sampleInfo$REP

# set pairwise comparisons between each group
sfc_vs_mes <- matrix(c(1, -1, 0, 0), nrow = 1)
sfc_vs_spe <- matrix(c(1, 0, -1, 0), nrow = 1)
sfc_vs_tes <- matrix(c(1, 0, 0, -1), nrow = 1)
mes_vs_spe <- matrix(c(0, 1, -1, 0), nrow = 1)
mes_vs_tes <- matrix(c(0, 1, 0, -1), nrow = 1)
spe_vs_tes <- matrix(c(0, 0, 1, -1), nrow = 1)

all_comps <- rbind(sfc_vs_mes, sfc_vs_spe, sfc_vs_tes,
                   mes_vs_spe, mes_vs_tes, spe_vs_tes)

rownames(all_comps) <- c("sfc-mesa", "sfc-sperm", "sfc-testis",
                         "mesa-sperm", "mesa-testis", "sperm-testis")
colnames(all_comps) <- c("sfc", "mesa", "sperm", "testis")

# do tests
if(!file.exists("output/male_rt_comparisons.rds")){

  male_rt_comparisons <- groupComparison(contrast.matrix = all_comps, data = male_rt_processed)

  saveRDS(male_rt_comparisons, "output/male_rt_comparisons.rds")
} else {
  male_rt_comparisons <- readRDS('output/male_rt_comparisons.rds')
}

# combine results
da_results <- male_rt_comparisons$ComparisonResult %>%
  mutate(DA = if_else(abs(log2FC) > 1 & adj.pvalue < 0.05, "SD", "NS")) %>%
  drop_na(DA)

da_results %>%
  ggplot(aes(x = log2FC, y = -log10(pvalue), colour = DA)) +
  geom_point() +
  facet_wrap(~ Label) +
  theme_bw() +
  theme(legend.position = "") +
  NULL
#ggsave('plots/male_rt/tissue_bias_volcanos.pdf', height = 4, width = 10, dpi = 600, useDingbats = FALSE)


# identify proteins most abundant in each tissue.
bias_cats <- da_results %>% filter(adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  mutate(category =
           case_when(`sfc-mesa` > 1 & `sfc-sperm` > 1 & `sfc-testis` > 1 ~ "sfc",
                     `sfc-mesa` < -1 & `mesa-sperm` > 1 & `mesa-testis` > 1 ~ "mes",
                     `sfc-sperm` < -1 & `mesa-sperm` < -1 & `sperm-testis` > 1 ~ "spe",
                     `sfc-testis` < -1 & `mesa-testis` < -1 & `sperm-testis` < -1 ~ "tes",
                     # proteins higher in sfcs and/or mesadenes
                     `sfc-sperm` > 1 & `sfc-testis` > 1 & `mesa-sperm` > 1 & `mesa-testis` > 1 ~ "sem",
                     # proteins higher in testis and/or sperm
                     `sfc-sperm` < -1 & `sfc-testis` < -1 & `mesa-sperm` < -1 & `mesa-testis` < -1 ~ "ger"
           ),
         only_in = case_when(is.na(category) == FALSE & Protein %in% c(sfc_only$Protein, mes_only$Protein, spe_only$Protein, tes_only$Protein) ~ "only",
                             is.na(category) == FALSE ~ "diff"),
         category = fct_relevel(category, 'tes', 'spe', 'ger', 'sem', 'sfc', "mes"),
         bias = case_when(category == "sem" ~ "soma",
                          category == "mes" ~ "soma",
                          category == "sfc" ~ "soma",
                          category == "spe" ~ "germ",
                          category == "tes" ~ "germ",
                          category == "ger" ~ "germ"))

bias_cats %>% count(bias)
bias_cats %>% count(category)
#write_csv(bias_cats, "output/male_rt/male_bias_cats.csv")

# number of proteins biased/unique to each tissue
bias_cats %>%
  count(category, only_in) %>% drop_na() %>%
  ggplot(aes(x = category, y = n, fill = only_in)) +
  geom_col() +
  #geom_col(position = "fill") +
  theme_bw() +
  theme() +
  NULL
#ggsave('plots/male_rt/tissue_bias_cols.pdf', height = 4, width = 5, dpi = 600, useDingbats = FALSE)

# percentage of proteins only identified in that tissue
bias_cats %>% count(category, only_in) %>% drop_na() %>%
  pivot_wider(names_from = only_in, values_from = n) %>%
  mutate(total = (diff + only),
         perc = (only/total) * 100)

# make lists of tissue-biased proteins
sfc.list <- da_results %>% filter(str_detect(Label, "sfc"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`sfc-mesa` > 1, `sfc-sperm` > 1, `sfc-testis` > 1) %>%
  mutate(only_in = if_else(Protein %in% sfc_only$Protein, "only", "diff"))

mes.list <- da_results %>% filter(str_detect(Label, "mesa"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`sfc-mesa` < -1, `mesa-sperm` > 1, `mesa-testis` > 1) %>%
  mutate(only_in = if_else(Protein %in% mes_only$Protein, "only", "diff"))

spe.list <- da_results %>% filter(str_detect(Label, "sperm"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`sfc-sperm` < -1, `mesa-sperm` < -1, `sperm-testis` > 1) %>%
  mutate(only_in = if_else(Protein %in% spe_only$Protein, "only", "diff"))

tes.list <- da_results %>% filter(str_detect(Label, "testis"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`sfc-testis` < -1, `mesa-testis` < -1, `sperm-testis` < -1) %>%
  mutate(only_in = if_else(Protein %in% tes_only$Protein, "only", "diff"))

# proteins more abundant in sfc/MES compared to sperm/testes
sem.list <- da_results %>% filter(str_detect(Label, "sfc|mes"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`sfc-sperm` > 1, `sfc-testis` > 1, `mesa-sperm` > 1, `mesa-testis` > 1) %>%
  mutate(only_in = "diff") %>%
  filter(!Protein %in% c(sfc.list$Protein, mes.list$Protein))

# proteins more abundant in sperm/testis compared to sfc/mes
ger.list <- da_results %>% filter(str_detect(Label, "sperm|testis"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`sfc-sperm` < -1, `sfc-testis` < -1, `mesa-sperm` < -1, `mesa-testis` < -1) %>%
  mutate(only_in = "diff") %>%
  filter(!Protein %in% c(spe.list$Protein, tes.list$Protein))


# sanity check sfc biased
norm_data %>%
  filter(Protein %in% sfc.list$Protein) %>%
  sample_n(20) %>%
  pivot_longer(cols = 2:17) %>%
  separate(name, into = c("tissue", "repl")) %>%
  ggplot(aes(x = tissue, y = value, colour = tissue)) +
  geom_jitter() +
  facet_wrap(~Protein, ncol = 10, scales = "free_y")


# count the number of proteins specific for each tissue.
tissue_biased_proteins <- bind_rows(sfc.list %>% select(Protein, only_in) %>% mutate(bias = "sfc"),
                                    mes.list %>% select(Protein, only_in) %>% mutate(bias = "mes"),
                                    spe.list %>% select(Protein, only_in) %>% mutate(bias = "spe"),
                                    tes.list %>% select(Protein, only_in) %>% mutate(bias = "tes"),

                                    sem.list %>% select(Protein, only_in) %>% mutate(bias = "sem"),
                                    ger.list %>% select(Protein, only_in) %>% mutate(bias = "ger")) %>%
  mutate(bias2 = case_when(bias == "sem" ~ "soma",
                           bias == "mes" ~ "soma",
                           bias == "sfc" ~ "soma",
                           bias == "spe" ~ "germ",
                           bias == "tes" ~ "germ",
                           bias == "ger" ~ "germ"))

#write_csv(tissue_biased_proteins, "output/topGO_lists/male_rt/tissue_biased_proteins.csv")

# number biased towards (mesadenes + sfcs) or (sperm + testes)
tissue_biased_proteins %>% dplyr::count(bias2)
# number biased towards each tissue
tissue_biased_proteins %>% dplyr::count(bias)


# overlaps
upset(fromList(split(tissue_biased_proteins$Protein[tissue_biased_proteins$bias2 == "germ"],
                     tissue_biased_proteins$bias[tissue_biased_proteins$bias2 == "germ"], drop = FALSE)),
      keep.order = TRUE)

#pdf('plots/germline_overlap.pdf', height = 4, width = 4)
plot(euler(c('spe' = 169 + 318, "tes" = 465 + 70,
             'spe&tes' = 204)),
     quantities = TRUE,
     fills = list(fill = Hiroshige.6[c(3, 1, 2)], alpha = .8))
#dev.off()

# overlaps
upset(fromList(split(tissue_biased_proteins$Protein[tissue_biased_proteins$bias2 == "soma"],
                     tissue_biased_proteins$bias[tissue_biased_proteins$bias2 == "soma"], drop = FALSE)),
      keep.order = TRUE)

#pdf('plots/sfp_overlap.pdf', height = 4, width = 4)
plot(euler(c("sfc" = 139 + 100, "mes" = 659 + 49,
             "sfc&mes" = 75)),
     quantities = TRUE,
     fills = list(fill = Hiroshige.6[4:6], alpha = .5))
#dev.off()


# >> make data frame ####
bias.frame <- norm_data %>% as.data.frame() %>%
  left_join(tissue_biased_proteins, by = "Protein") %>%
  mutate(sigp = if_else(Protein %in% signal_peps, "sigp", "not"),
         bias = fct_relevel(bias, 'tes', 'spe', 'ger', 'sem', 'sfc', "mes"))

# long format
bias.frame.long <- bias.frame %>%
  pivot_longer(cols = 2:17) %>% drop_na(value) %>%
  separate(name, into = c("tissue", "repl"), remove = FALSE)

# normalised means per tissue
wide_mns <- quantification(male_rt_processed, type = "group") %>%
  left_join(tissue_biased_proteins, by = c("Protein"))

# count number of proteins identified in each tissue
wide_mns %>%
  pivot_longer(cols = 2:5) %>%
  drop_na(value) %>%
  count(name)


# > Signal peptides ####
# total number of proteins w/ signal peptides - to calculated 'expected'
TotalGeneNumber <- bias.frame %>% dplyr::count(bias2, name = "N")

# proportion of all proteins containing a signal peptide
prop.total <- bias.frame %>% dplyr::count(sigp, name = "N") %>%
  mutate(pr.total = N/sum(N))

# Number of proteins in each tissue with signal peptide - 'observed'
gene.no <- bias.frame %>% dplyr::count(bias2, sigp, name = "obs_genes") %>% filter(sigp == "sigp")

# Calculate observed and expected no. genes in each comparison to do X^2 test
tiss_dist <- TotalGeneNumber %>% inner_join(gene.no) %>% inner_join(prop.total, by = "sigp") %>%
  # Calculate X^2 statistics
  mutate(exp_genes = N.x * pr.total, # expected no. genes
         obs_exp = obs_genes/exp_genes, # observed / expected no. genes
         X2 = (obs_genes-exp_genes)^2/exp_genes, # calculate X^2 stat
         pval = 1 - (pchisq(X2, df = 1)) # get pvalue
  )

# FDR corrected pval
tiss_dist$FDR <- p.adjust(tiss_dist$pval, method = 'fdr')
tiss_dist <- tiss_dist %>%
  mutate(sigLabel = case_when(FDR > 0.05 ~ "",
                              FDR <= 0.05 & FDR > 0.01 ~ "*",
                              FDR <= 0.01 & FDR > 0.001 ~ "**",
                              FDR <= 0.001 ~ "***"),
         across(where(is.double), ~round(.x, 3)))

signal_pep_cols <- tiss_dist %>% drop_na(bias2) %>%
  ggplot(aes(x = bias2, y = obs_exp, fill = bias2)) +
  geom_col() +
  geom_hline(yintercept = 1, lty = 2) +
  scale_fill_manual(values = Hiroshige.6[c(1, 6)]) +
  scale_x_discrete(labels = c("Sperm + \nTestes", "Seminal fluid\nproducing")) +
  labs(x = "Tissue", y = "observed/expected\n no. proteins") +
  theme_bw() +
  theme(legend.position = "",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5)) +
  geom_text(aes(label = sigLabel),
            size = 5, colour = "black") +
  geom_text(aes(y = -0.05, label = paste0(obs_genes, '/', N.x)),
            size = 3, colour = "black") +
  #ggsave('plots/male_rt/signal_pep_cols.pdf', height = 2, width = 2, dpi = 600, useDingbats = FALSE)
  NULL
signal_pep_cols


# > Define secretome/sperm proteome ####

## subset proteins from SFCs/Mesadenes with or without signal peptide sequence
cimex_secretome <- tissue_biased_proteins %>%
  mutate(sigp = if_else(Protein %in% signal_peps, "sigp", "not")) %>%
  filter(sigp == "sigp", bias %in% c("sfc", "mes", "sem"))
#write_csv(cimex_secretome, "output/male_rt/cimex_secretome.csv")

cimex_ags <- tissue_biased_proteins %>%
  filter(bias %in% c("sfc", "mes", "sem"), !Protein %in% cimex_secretome$Protein)
#write_csv(cimex_ags, "output/male_rt/cimex_ags.csv")

## proteins most abundant in sperm
sperm_proteome <- tissue_biased_proteins %>%
  filter(!bias %in% c("tes", "sfc", "mes", "sem"))
#write_csv(sperm_proteome, "output/male_rt/sperm_proteome.csv")

## remaining sperm/testes biased proteins
germ_proteome <- tissue_biased_proteins %>%
  filter(bias %in% c("tes"), !Protein %in% sperm_proteome$Protein)
#write_csv(germ_proteome, "output/male_rt/germ_proteome.csv")

# combine lists of tissue-biased proteins/sfps
male_cats <- list(sfp = cimex_secretome$Protein,
                  ags = cimex_ags$Protein,
                  spe = sperm_proteome$Protein,
                  tes = germ_proteome$Protein)
#write_rds(male_cats, "output/male_rt/male_cats.rds")

male_cats_long <- male_cats %>%
  reshape2::melt() %>% dplyr::rename(accession = value, category = L1)


### heatmap ###

hm_data <- pheatmap:::scale_rows(bias.frame %>%
                                   dplyr::select(2:17) %>%
                                   mutate(across(everything(), ~replace_na(.x, 0)))) %>%
  mutate(Protein = bias.frame$Protein) %>%
  relocate(Protein, .before = mesa_1) %>%
  left_join(male_cats_long, c("Protein" = "accession"))

# top annotations
vmha <- HeatmapAnnotation(
  tissue = str_sub(colnames(hm_data[, 2:17]), 1, 3),
  col = list(tissue = setNames(rev(Hiroshige.4),
                               unique(str_sub(colnames(hm_data[, 2:17]), 1, 3))))
)

# right annotations
right_lab <- rowAnnotation(
  category = hm_data$category,
  col = list(category = setNames(Hiroshige.4,
                                 #rainbow(n_distinct(na.omit(bias.frame$bias))),
                                 unique(na.omit(hm_data$category)))
  ),
  na_col = NA,
  title = NULL,
  #                          show_legend = FALSE,
  show_annotation_name = FALSE)

### plot
#pdf('plots/male_rt/bias_scaled_hm.pdf', height = 6, width = 5)
Heatmap(as.matrix(hm_data[, 2:17]),
        na_col = NA,
        heatmap_legend_param = list(title = "log2(intensities)",
                                    direction = "horizontal"),
        #left_annotation = left_lab,
        right_annotation = right_lab,
        top_annotation = vmha,
        show_row_names = FALSE,
        show_column_names = FALSE,
        #show_row_dend = FALSE,
        clustering_method_rows = "complete",
        column_split = 4,
        column_gap = unit(0, "mm"),
        row_title = NULL,
        column_title = NULL,
        use_raster = FALSE)
#dev.off()


# >> expression of Sfps in other tissues ####
wide_mns %>% filter(Protein %in% cimex_secretome$Protein) %>%
  pivot_longer(cols = 2:5) %>%
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_boxplot() +
  scale_fill_manual(values = rev(Hiroshige.4)) +
  #scale_x_discrete(labels = c("Sperm + \nTestes", "Seminal fluid\nproducing")) +
  labs(x = "Tissue", y = "log2(abundance)") +
  theme_bw() +
  theme(legend.position = "") +
  NULL

# find proteins in the top 20% ranked abundance in sperm/testes that belong to the putative seminal fluid proteome
wide_mns %>%
  select(Protein, sperm, bias) %>%
  arrange(desc(sperm)) %>%
  mutate(rnk = percent_rank(sperm) * 100) %>%
  filter(Protein %in% cimex_secretome$Protein, rnk > 80)

wide_mns %>%
  select(Protein, testis, bias) %>%
  arrange(desc(testis)) %>%
  mutate(rnk = percent_rank(testis) * 100) %>%
  filter(Protein %in% cimex_secretome$Protein, rnk > 80)

# histograms showing abundance of putative sfps in each tissue/sperm
wide_mns %>% filter(Protein %in% cimex_secretome$Protein) %>%
  pivot_longer(cols = 2:5) %>%
  ggplot(aes(x = value, fill = name)) +
  geom_histogram() +
  scale_fill_manual(values = rev(Hiroshige.4)) +
  # dotted line shows average abundance of all proteins in that tissue
  geom_vline(data = apply(wide_mns[, 2:5], 2, median, na.rm = TRUE) %>% plyr::ldply(.id = "name"), aes(xintercept = V1), lty = 2) +
  facet_wrap(~ name) +
  theme_bw() +
  theme(legend.position = "") +
  NULL


# > Evolutionary dynamics ####

# >> chromosomal locations ####
chm_locs <- read_csv("output/chrom_locs.csv")

level_order <- c("chm_1", "chm_2", "chm_3", "chm_4", "chm_5", "chm_6",
                 "chm_7", "chm_8", "chm_9", "chm_10", "chm_11", "chm_12",
                 "chm_13", "chm_X1", "chm_X2")

# count total number and proportion of proteins on each chromosome
TotalProtNumber <- chm_locs %>%
  filter(accession %in% da_results$Protein) %>%
  count(chm, name = "N") %>%
  mutate(pr.total = N/sum(N))

# number of genes in data on each chromosome (N observed)
prot.no <- male_cats_long %>%
  left_join(chm_locs, by = "accession") %>%
  distinct(accession, .keep_all = TRUE)%>%
  #drop_na(category) %>%
  group_by(chm, category) %>%
  summarise(obs_genes = n())

# count number of tissue biased genes to calculate expected
obs_N_male <- male_cats_long %>% count(category)

# Calculate observed and expected no. genes in each comparison on each chromosome to do X^2 test
chm_dist_male <- prot.no %>%
  left_join(obs_N_male) %>% inner_join(TotalProtNumber) %>%
  # Calculate X^2 statistics
  mutate(exp_genes = round(n * pr.total), # expected no. genes
         obs_exp = obs_genes/exp_genes, # observed / expected no. genes
         X2 = (obs_genes-exp_genes)^2/exp_genes, # calculate X^2 stat
         pval = 1 - (pchisq(X2, df = 1)) # get pvalue
  )

# FDR corrected pval
chm_dist_male$FDR <- p.adjust(chm_dist_male$pval, method = 'fdr')
chm_dist_male <- chm_dist_male %>%
  mutate(sigLabel = case_when(FDR >= 0.05 ~ "",
                              FDR < 0.05 & FDR > 0.01 ~ "*",
                              FDR <= 0.01 & FDR > 0.001 ~ "**",
                              FDR <= 0.001 ~ "***"),
         across(where(is.double), ~round(.x, 3)))

#plot
chm_plot <- bind_rows(chm_dist_male %>% filter(str_detect(chm, "X")),
                      chm_dist_male %>% filter(!str_detect(chm, "X")) %>%
                        group_by(category) %>%
                        summarise(obs_exp = mean(obs_exp)) %>% mutate(chm = "A")) %>%
  ggplot(aes(x = chm, y = obs_exp, fill = category)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, lty = 2) +
  scale_fill_manual(values = rev(Hiroshige.4),
                    labels = c("Mesadenes", "Sfps", "Sperm", "Testes")
  ) +
  scale_x_discrete(labels = c("Autos.", "X1", "X2")) +
  labs(x = "Chromosome", y = "observed/expected\n no. proteins") +
  theme_bw() +
  theme(#legend.position = "",
    legend.title = element_blank()) +
  geom_text(aes(label = sigLabel), position = position_dodge(width = .9),
            size = 5, colour = "black") +
  # geom_text(aes(y = -0.05, label = paste0(obs_genes, '/', N)),
  #           size = 5, colour = "black") +
  #facet_wrap(~chm) +
  #ggsave('plots/male_rt/chm_dist_X-auto.pdf', height = 2.5, width = 4, dpi = 600, useDingbats = FALSE)
  NULL
chm_plot

# plot all chromosomes
chm_dist_male %>% mutate(chr = gsub("chm_", "", x = chm)) %>%
  ggplot(aes(x = factor(chr, level = c(1:13, "X1", "X2")), y = obs_exp, fill = category)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, lty = 2) +
  scale_fill_manual(values = rev(Hiroshige.4),
                    labels = c("Mesadenes", "Sfps", "Sperm", "Testes")
  ) +
  #scale_x_discrete(labels = c("Autos.", "X1", "X2")) +
  labs(x = "Chromosome", y = "observed/expected\n no. proteins") +
  theme_bw() +
  theme(#legend.position = "",
    legend.title = element_blank()) +
  geom_text(aes(label = sigLabel), position = position_dodge(width = .9),
            size = 5, colour = "black") +
  # geom_text(aes(y = -0.05, label = paste0(obs_genes, '/', N)),
  #           size = 5, colour = "black") +
  #facet_wrap(~chm) +
  #ggsave('plots/male_rt/chm_dist_all.pdf', height = 2.5, width = 13, dpi = 600, useDingbats = FALSE)
  NULL


# >> dN/dS ####
one_ratio <- read.csv('../Bedbug_genomics/run_dnds/codeml_model0_full_results.csv') %>%
  # filter dS < 2...
  filter(dS < 2) %>%
  # filter S x dS <= 1
  mutate(S.dS = S * dS) %>%
  filter(S.dS >= 1) %>%
  left_join(uniprot2entrez, by = c("GeneID" = "gene_id")) %>%
  # bin evo rates - each bin contains same number of genes - higher bin = higher rate
  mutate(bin = ntile(omega, 10)) %>%
  left_join(male_cats_long)

# plot
one_ratio %>% filter(category %in% c("tes", "spe", "sfp", "ags")) %>%
  ggplot(aes(x = factor(category, level = c("tes", "spe", "sfp", "ags")), y = omega, colour = category)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  geom_jitter(alpha = .05) +
  scale_x_discrete(labels = c("Testes", "Sperm", "Sfps", "Mesadenes")) +
  scale_colour_manual(values = Hiroshige.4[c(4, 3, 2, 1)]) +
  labs(x = 'Gene set', y = 'Pairwise dN/dS') +
  theme_bw() +
  theme(legend.position = "") +
  NULL

# perform Mann-Whitney U tests comparing each category to rest of genes
witney_tests <- do.call(bind_rows, lapply(c('ags', 'sfp', 'spe', 'tes'), function(x) {
  one_ratio %>%
    mutate(category = if_else(grepl(x, category), x, 'other')) %>%
    do(fit = broom::tidy(wilcox.test(omega ~ category, data = .))) %>%
    unnest(fit) %>%
    mutate(category = x)
}))

witney_tests$FDR <- p.adjust(witney_tests$p.value, method = 'fdr')
witney_tests <- witney_tests %>%
  mutate(pval = ifelse(FDR < 0.001, '< 0.001', round(p.value, 3)),
         sigLabel = case_when(FDR < 0.001 ~ "***",
                              FDR < 0.01 ~ "**",
                              FDR < 0.05 ~ "*",
                              TRUE ~ ''))

witney_tests

# average genome wide dN/dS
one_ratio %>%
  summarise(N = n(),
            mn = mean(omega, na.rm = TRUE),
            se = sd(omega, na.rm = TRUE)/sqrt(N))

evo_plot_mns <- one_ratio %>%
  #drop_na(category) %>%
  group_by(category) %>%
  summarise(N = n(),
            mn = mean(omega, na.rm = TRUE),
            se = sd(omega, na.rm = TRUE)/sqrt(N)) %>%

  left_join(witney_tests) %>%

  mutate(category = fct_relevel(category, 'tes', 'spe', "sfp", 'ags')) %>%
  filter(category != "other") %>%

  #plot
  ggplot(aes(x = category, y = mn, colour = category)) +
  # add genome average
  geom_rect(aes(ymin= 0.1673432 - 0.00136304, ymax = 0.1673432 + 0.00136304, xmin = -Inf, xmax = Inf),
            fill='grey', colour = NA, alpha=0.1) +
  geom_hline(yintercept = 0.1673432, lty = 2, alpha = 0.5) +

  # points and errorbars
  geom_errorbar(aes(ymin = mn - se, ymax = mn + se), width = .25, lwd = 1.2) +
  geom_point(size = 4, pch = 21, stroke = 2, fill = "white") +
  #geom_hline(yintercept = mean(one_ratio$omega), lty = 2) +
  scale_x_discrete(labels = c("Testes", "Sperm", "Sfps", "Mesadenes")) +
  scale_colour_manual(values = Hiroshige.4) +
  labs(x = 'Gene set', y = 'Pairwise dN/dS') +
  theme_bw() +
  theme(legend.position = "") +
  geom_text(aes(y = 0.1, label = N),
            position = position_dodge(width = .9),
            size = 3, colour = "black") +
  geom_text(aes(y = mn + se + 0.0025, label = sigLabel),
            size = 5, colour = "black") +
  #ggsave('plots/male_rt/dnds.plot.pdf', height = 3.5, width = 3.5, dpi = 600, useDingbats = FALSE)
  NULL
evo_plot_mns

# pairwise tests
broom::tidy(pairwise.wilcox.test(x = one_ratio$omega, g = one_ratio$category, p.adjust.method = "fdr")) %>%
  mutate(p.value = ifelse(p.value < 0.001, '< 0.001', round(p.value, 3)))


# combined plot for figure 1
plot_grid(
  plot_grid(signal_pep_cols + theme(axis.text.x = element_text(size = 10)),
            NULL,
            ncol = 1, nrow = 2, rel_heights = c(1.6, 1)),
  chm_plot +
    theme(legend.position = c(0.2, 0.83),
          legend.background = element_blank()),
  evo_plot_mns,
  nrow = 1,
  ncol = 3, rel_widths = c(0.65, 1, 1)
)
#ggsave('plots/male_rt/sig_chm_rates.pdf', height = 4, width = 10, dpi = 600, useDingbats = FALSE)



# > Orthology comparisons ----------------------------------------------------------
hogs <- read_tsv('../proteomes/OrthoFinder/Results_27-Jun-2025/Phylogenetic_Hierarchical_Orthogroups/N0.tsv')
hogs_long <- hogs %>%
  pivot_longer(cols = 4:12) %>%
  separate_rows(value, sep = ", ") %>%
  drop_na(value) %>%
  dplyr::rename(Protein = value, species = name)

# get HOGs for cimex proteins
cimex_long <- hogs_long %>% filter(species == "Cimex")

### Fly stuff ###
# load Dmel info from flybase
flybase2uniprot <- read_tsv("data/FlyBase/FlyBaseAllGenes_2_Uniprot_2025_05_07.tsv") %>%
  mutate(gene_name = gsub("Dmel\\\\", "", x = gsub(" .*", "", x = `Gene Names`))) %>%
  dplyr::rename(FBgn = From, Accession = Entry) %>%
  left_join(hogs_long %>%
              filter(species == "Drosophila") %>%
              dplyr::select(HOG, Protein),
            by = c("Accession" = "Protein")) %>%
  # rename odd names with their CG##### instead
  mutate(gene_name = case_when(grepl("anon", x = gene_name) ~ str_extract(`Gene Names`, "\\S+$") %>% str_remove("^Dmel_"),
                               grepl("\\.|:", x = gene_name) ~ str_extract(`Gene Names`, "\\S+$") %>% str_remove("^Dmel_"),
                               TRUE ~ gene_name))

## Dmel sperm proteome (DmSP3) downloaded from Garlovsky et al. 2022 Mol. Cell. Prot.
DmSP3 <- read.csv('../../../KB_10MSE/output/DmSP_which_proteome.csv') %>%
  dplyr::select(FBgn, SYMBOL, Chrm = LOCATION_ARM, Sfp, mean.perc, grand.mean) %>%
  left_join(flybase2uniprot %>% dplyr::select(FBgn, Accession, gene_name, HOG)) %>%
  distinct(Accession, .keep_all = TRUE) %>% drop_na(HOG)

# dmel slaps - S-Lap 6 a.k.a. loopin
dmel_slaps <- DmSP3 %>% filter(str_detect(SYMBOL, "S-Lap|loopin"))

# Dmel SFPs downloaded from Wigby et al. 2020 Phil. Trans. B. supplementary material
wigbySFP <- read.csv('../../../DmSPIII_Garlovsky_Karr/DmSP3/data/dmel_SFPs_wigby_etal2020.csv') %>%
  # only high confidence Sfps
  filter(category == 'highconf') %>%
  left_join(flybase2uniprot %>% dplyr::select(FBgn, Accession, gene_name, HOG)) %>%
  drop_na(HOG)

### ### ###

### Mosquito stuff ###
## Aedes egypti sperm and seminal fluid proteome from Degner et al. 2019 Mol. Cell. Prot.
Ae_sp <- readxl::read_excel("../Degner_etal_2019/140496_1_supp_242861_p8z3rl.xlsx") %>%
  left_join(read.delim("../Degner_etal_2019/Degner_RefSeq2UniProt_2023_08_10.tsv"),
            by = c("RefSeq ID" = "From")) %>%
  janitor::clean_names() %>%
  left_join(hogs_long %>% filter(species == "Aedes") %>% select(-species), by = c("entry" = "Protein")) %>%
  dplyr::select(vector_base_id:ref_seq_id, Accession = entry, HOG,
                classification_in_current_study, orthologs_in_drosophila_melanogaster) %>%
  drop_na(HOG)

aedes_sperm <- Ae_sp %>% filter(str_detect(classification_in_current_study, "Sperm"))
aedes_sfps <- Ae_sp %>% filter(str_detect(classification_in_current_study, "SFP"))

### ### ###

# get all the ortholog keepings the first distinct FBgn
cimex2mel <- hogs_long %>%
  filter(species %in% c("Cimex", "Drosophila")) %>%
  left_join(flybase2uniprot %>% select(FBgn, Accession, gene_name, HOG),
            by = c("HOG")) %>%
  filter(species == "Cimex") %>%
  drop_na(FBgn) %>%
  distinct(FBgn, .keep_all = TRUE) %>%
  select(-c(OG, `Gene Tree Parent Clade`, species))


# > S-Laps ####
# number of S-Lap orthologs in each species
hogs_long %>%
  filter(HOG %in% na.omit(dmel_slaps$HOG)) %>%
  count(species)

# Cimex SLAPs based on HOGs
cimex_slaps <- hogs_long %>%
  filter(HOG %in% na.omit(dmel_slaps$HOG), species == "Cimex")

# count orthologs belonging to each SLAP group
cimex_slaps %>% dplyr::count(HOG)

# look at SLAPs in original data - how many unique peptides each?
ms_data %>% filter(accession %in% cimex_slaps$Protein) %>%
  select(accession, number_unique_peptides, coverage_percent, mw_k_da)

# BLAST results from NCBI giving top ranked orthologs
ncbi_slaps <- read_csv("output/male_rt/SLAP_stuff/NCBI_BLASTp_SLAPs.csv") %>%
  inner_join(cimex_slaps %>% dplyr::select(-species))

# reciprocal best BLAST results
RBB <- read_tsv("output/male_rt/blast_results/reciprocal_best_hits.tsv") %>%
  left_join(dmel_slaps %>% dplyr::select(FBgn, Accession, gene_name), by = c("Drosophila_ID" = "Accession"))


### heatmap ###
slap_hm_data <- hm_data %>%
  filter(Protein %in% cimex_slaps$Protein) %>%
  left_join(RBB, by = c("Protein" = "Cimex_ID")) %>%
  mutate(Protein2 = if_else(is.na(gene_name) == FALSE, paste0(gene_name, " (", Protein, ")"), Protein))
rownames(slap_hm_data) <- slap_hm_data$Protein

# make tree for S-Laps
library(ggtree)
# load tree of MAFFT alignment including Dmel and LAP3 then manually pruned to only show Cimex
slap_trim <- read.tree("output/male_rt/SLAP_stuff/slap_trim.nwk")
slap_tree <- ggtree(slap_trim, branch.length = "none") %<+% slap_hm_data +
  geom_tiplab(aes(label = Protein2))

# z-transformed abundance of Cimex SLAPs
gheatmap(slap_tree,
         slap_hm_data %>% select(17:2),
         offset = 2, width = 1, color = "black",
         colnames_angle=45, hjust=1,
         colnames = TRUE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  NULL
#ggsave('plots/male_rt/SLAP_PLOTS/SLAP_heatmap.pdf', height = 3.5, width = 8, dpi = 600, useDingbats = FALSE)


# >> SLAP catalytic activity ####

# load MAFFT bootstrapped tree
best_tree_bs <- read_rds("output/male_rt/SLAP_stuff/mafft_output/best_tree_bs.rds")

# labels for all SLAP/grannysmith orthologs
slap_lab <- data.frame(Protein = best_tree_bs$tip.label) %>%
  left_join(hogs_long %>% select(-OG, -`Gene Tree Parent Clade`), by = c("Protein")) %>%
  left_join(dmel_slaps %>% dplyr::select(Accession, gene_name), by = c("Protein" = "Accession")) %>%
  mutate(#gene_name = if_else(is.na(gene_name) == TRUE, Protein, gene_name),
    gene_name = case_when(Protein == "Q9V3D8" ~ "granny-smith",
                          Protein == "A0A0C4FEI8" ~ "granny-smith F",
                          is.na(gene_name) == FALSE ~ gene_name,
                          is.na(gene_name) == TRUE ~ Protein,
                          TRUE ~ Protein),
    #species = if_else(is.na(species) == TRUE, "taurus", species),
    species = case_when(species == "Chemipterus" ~ "hemipterus",
                        is.na(species) == TRUE ~ "taurus",
                        TRUE ~ species),
    Order = case_when(species %in% c("Drosophila", "Aedes") ~ "Diptera",
                      species == "Callosobruchus" ~ "Coleoptera",
                      species == "Gryllus" ~ "Orthoptera",
                      species == "taurus" ~ "Artiodactyla",
                      TRUE ~ "Hemiptera"),
    Genus = case_when(species == "Cimex"| species == "hemipterus" ~ "Cimex",
                      species == "taurus" ~ "Bos",
                      TRUE ~ species)) %>%
  left_join(RBB %>% select(Protein = Cimex_ID, Protein2 = gene_name), by = c("Protein")) %>%
  mutate(Protein2 = case_when(species == "Drosophila" ~ gene_name,
                              is.na(Protein2) == FALSE ~ paste0(Protein2, " (", Protein, ")"),
                              TRUE ~ Protein),
         species = factor(species, level = c("taurus", "Aedes", "Drosophila",
                                             "Callosobruchus",
                                             "hemipterus", "Cimex", "Rhodnius", "Triatoma", "Acyrthosiphon",
                                             "Gryllus")))

# One letter to Three letter aa table
residue_123 <- data.frame(one_letter = c("A", "R", "N", "D",
                                         "C", "Q", "E", "G",
                                         "H", "I", "L", "K",
                                         "M", "F", "P", "S",
                                         "T", "W", "Y", "V",
                                         "-", "X"),
                          three_letter = c("Ala", "Arg", "Asn", "Asp",
                                           "Cys", "Gln", "Glu", "Gly",
                                           "His", "Ile", "Leu", "Lys",
                                           "Met", "Phe", "Pro", "Ser",
                                           "Thr", "Trp", "Tyr", "Val",
                                           "---", "Unk"),
                          stringsAsFactors = FALSE)

# load results of motif based analysis
motifs <- read_tsv("SLAP_activity/region_constrained_analysis/region_constrained_details.tsv") %>%
  left_join(residue_123, by = c("residue" = "one_letter"))

# matrix of amino acids at each locus
residue_matrix <- motifs %>%
  select(protein_id, site_name, three_letter) %>%
  pivot_wider(names_from = site_name,
              values_from = three_letter) %>%
  relocate(cat_1, .after = metal_5)
rownames(residue_matrix) <- residue_matrix$protein_id
# make matrix
resmatrix <- as.matrix(residue_matrix[, -1])
rownames(resmatrix) <- residue_matrix$protein_id
colnames(resmatrix) <- c("Lys 327", "Asp 332", "Asp 350", "Lys 409", "Glu 411", "Lys 339", "Arg 413")

# matrix of scores
score_matrix <- motifs %>%
  select(protein_id, site_name, score) %>%
  pivot_wider(names_from = site_name,
              values_from = score) %>%
  relocate(cat_1, .after = metal_5)
rownames(score_matrix) <- score_matrix$protein_id

# factorised scores
score_matrix <- motifs %>%
  select(protein_id, site_name, tier) %>%
  pivot_wider(names_from = site_name,
              values_from = tier) %>%
  relocate(cat_1, .after = metal_5)
rownames(score_matrix) <- score_matrix$protein_id
# make matrix
motifmatrix <- as.matrix(score_matrix[, -1])
rownames(motifmatrix) <- score_matrix$protein_id
colnames(motifmatrix) <- c("Lys 327", "Asp 332", "Asp 350", "Lys 409", "Glu 411", "Lys 339", "Arg 413")

# Compare Drosophila S-Lap results to those from Dorus et al. 2011 PLOS Genetics. There are some changes to the categorisation, likely arising from updates to the Drosophila aa sequences in the last 15 years
dmel_tree <- read.tree(text = "((((((Q7K2S9,A1Z9G3),Q7K5K9),(Q961W5,(Q95R35,Q500X4)))),((Q9VSM6,Q9VSM7))),(Q9V3D8));")

# make tree and rotate to match Dorus et al.
dmel_tree_plot <- ggtree(dmel_tree, branch.length = "none") %<+% slap_lab +
  #geom_text(aes(label=node), hjust=-.3) +
  geom_tiplab(aes(label = gene_name))
dmel_tree_plot <- ggtree::rotate(dmel_tree_plot, 19)
dmel_tree_plot <- ggtree::rotate(dmel_tree_plot, 15)
dmel_tree_plot <- ggtree::rotate(dmel_tree_plot, 14)
dmel_tree_plot <- ggtree::rotate(dmel_tree_plot, 13)

dmel_res <- gheatmap(dmel_tree_plot, resmatrix, offset=0.5, width=1,
         colnames = TRUE, legend_title="substitution", color = "black",
         colnames_angle = 45, hjust=1) +
  scale_fill_viridis_d(option = "H") +
  NULL
dmel_res
#ggsave("plots/male_rt/SLAP_PLOTS/DMEL_SLAPs_RES.pdf", width = 10, height = 12)

gheatmap(dmel_res, motifmatrix, offset=7, width=1,
         colnames = TRUE, legend_title="substitution", color = "black",
         colnames_angle = 45, hjust=1) +
  scale_fill_viridis_d(option = "H") +
  NULL
#ggsave("plots/male_rt/SLAP_PLOTS/DMEL_SLAPs_COMP.pdf", width = 20, height = 8)

gheatmap(dmel_tree_plot, motifmatrix, offset=0.5, width=1,
         colnames = TRUE, legend_title="substitution", color = "black",
         colnames_angle = 45, hjust=1) +
  scale_fill_manual(values = c("lightpink",
                               "red", "white")) +
  NULL


# tree of all S-Lap/grsm orthologs
all_slap_orths_tree <- ggtree(best_tree_bs, branch.length = "none") %<+% slap_lab +
  #geom_text(aes(label=node), hjust=-.3) +
  geom_cladelabel(node = 169, label = "granny-smith", align = TRUE, offset = 11, barsize = 1.5) +
  geom_cladelabel(node = 155, label = "S-Lap cluster 2", align = TRUE, offset = 11, barsize = 1.5) +
  geom_cladelabel(node = 122, label = "S-Lap cluster 1", align = TRUE, offset = 11, barsize = 1.5) +
  geom_tiplab(aes(label = Protein2))

# rotate grannysmith to top
all_slap_orths_tree <- ggtree::rotate(all_slap_orths_tree, 101)
all_slap_orths_tree <- ggtree::rotate(all_slap_orths_tree, 87)

# add motif matrix to tree
all_slap_orths_map <- gheatmap(all_slap_orths_tree, motifmatrix, offset = 4.6, width = 0.4,
         colnames = TRUE, legend_title="substitution", color = "black",
         colnames_angle = 45, hjust=1) +
  scale_fill_manual(values = c("lightpink",
                               "red", "white")) +
  NULL

# >> Alphafold results ####
# load zinc ligand results from alphafold
zinc_binding <- read_tsv("SLAP_activity/SLAP_orths_zn_zn/key_file.tsv") %>%
  mutate(PROT = gsub(".json", "", x = gsub(".*_", "", x = filename))) %>%
  select(PROT, protein_id) %>%
  inner_join(read_csv("SLAP_activity/af3_comparisons/zinc_summary_chimerax.csv") %>%
               mutate(PROT = gsub("_.*", "", x = gsub("af_protein_ligand_", "", x = protein))),
             by = "PROT") %>%
  select(-PROT)

# make matrix
zn <- zinc_binding %>%
  select(protein_id, zn_index, call) %>%
  pivot_wider(names_from = zn_index, values_from = call)
zn_max <- as.matrix(zn[, -1])
rownames(zn_max) <- zn$protein_id

# add alphafold results to the plot
gheatmap(all_slap_orths_map, zn_max, offset=9.5, width=0.1,
         colnames = TRUE, legend_title="binding", color = "black",
         colnames_angle = 45, hjust=1) +
  #scale_fill_viridis_d() +
  scale_fill_manual(values = c("lightpink", "red",
                               #"red", "white",
                               viridis::viridis(n = 4)[2:3],
                               "white",
                               #"grey",
                               viridis::viridis(n = 4)[4])) +
  NULL
#ggsave("plots/male_rt/SLAP_PLOTS/ml_tree_SLAPs.pdf", width = 10, height = 12)

# slightly truncate the long branch to granny-smith for plotting
best_tree_bs2 <- best_tree_bs
best_tree_bs2$edge.length[168] <- 0.6

ggtree(treeio::drop.tip(best_tree_bs2, "P00727"), layout = "equal_angle") %<+% slap_lab +
  #geom_tiplab(aes(label = Protein2)) +
  geom_hilight(node=93, fill = "#D81B60", type="encircle") +
  geom_hilight(node=97, fill = "#1E88E5", type="encircle") +
  geom_hilight(node=158, fill = "#FFC107", type="encircle") +
  geom_tippoint(aes(colour = Order, shape = species, fill = Order), stroke = 1.5, size = 5) +
  #scale_colour_viridis_d(option = "H") +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  scale_shape_manual(guide = "none",
                     values = c(16, # Aedes
                                21, # Drosphila
                                16, # beetle
                                16, # hempiterus
                                21, # Cimex
                                17, # Rhodnius
                                18, # Triatoma
                                15, # aphid
                                16) # Gryllus
  ) +
  NULL
#ggsave("plots/male_rt/SLAP_PLOTS/SLAPs_TREE.pdf", width = 8, height = 8)


# plot rooted tree without branch lengths
all_slap_circle <- ggtree(treeio::drop.tip(best_tree_bs, "P00727"),
                          branch.length = "none",
                          layout = "circular") %<+% slap_lab +
  #geom_text(aes(label=node), hjust=-.3) +
  #geom_tiplab(aes(label = Protein2)) +
  #theme(legend.position = "") +
  # geom_hilight(node=93, fill = "#D81B60", type="rect") +
  # geom_hilight(node=97, fill = "#1E88E5", type="rect") +
  # geom_hilight(node=158, fill = "#FFC107", type="rect") +

  geom_tippoint(aes(colour = Order, shape = species, fill = Order), size = 5, stroke = 1.5) +

  # geom_cladelabel(node = 122, label = "Drosophila", align = TRUE, offset = 7, barsize = 1.5) +
  # geom_cladelabel(node = 155, label = "Drosophila", align = TRUE, offset = 7, barsize = 1.5) +
  # geom_cladelabel(node = 81, label = "Drosophila", align = TRUE, offset = 7, barsize = 1.5) +

  geom_cladelabel(node = 93, label = "Group 1", align = TRUE, offset = 1, barsize = 1.5, colour = "#D81B60") +
  geom_cladelabel(node = 97, label = "Group 2", align = TRUE, offset = 1, barsize = 1.5, colour = "#1E88E5") +
  geom_cladelabel(node = 158, label = "granny-smith", align = TRUE, offset = 1, barsize = 1.5, colour = "#FFC107") +

  #geom_cladelabel(node = 94, label = "Group 1", align = TRUE, offset = 15, barsize = 1.5) +
  #geom_cladelabel(node = 98, label = "Group 2", align = TRUE, offset = 15, barsize = 1.5) +
  #geom_cladelabel(node = 159, label = "granny-smith", align = TRUE, offset = 15, barsize = 1.5) +
  #scale_colour_manual(values = ) +
  scale_colour_viridis_d(option = "H") +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  scale_shape_manual(guide = "none",
                     values = c(16, # Aedes
                                21, # Drosphila
                                16, # beetle
                                16, # hempiterus
                                21, # Cimex
                                17, # Rhodnius
                                18, # Triatoma
                                15, # aphid
                                16) # Gryllus
                     ) +
  NULL
#all_slap_circle <- ggtree::rotate(all_slap_circle, 101)
all_slap_circle

gheatmap(all_slap_circle, zn_max, offset=0, width=0.1,
         colnames = TRUE, legend_title="binding", color = "black",
         colnames_angle = 45, hjust=1) +
  scale_fill_viridis_d() +
  #scale_fill_manual(values = c("red", "white", "grey")) +
  NULL

p2 <- gheatmap(all_slap_circle, motifmatrix, offset = 0, width = 0.4,
               colnames = TRUE, legend_title="substitution", color = "black",
               colnames_angle = 45, hjust=1) +
  scale_fill_manual(values = c("orange", "red", "white")) +
  NULL

gheatmap(p2, zn_max, offset=5, width=0.1,
         colnames = TRUE, legend_title="binding", color = "black",
         colnames_angle = 45, hjust=1) +
  # scale_fill_manual(values = c("pink", "red",
  #                              "red", "white", "white", "grey")) +
  NULL


# plot only cimex full results
# tree of Cimex orthologs with S-Lap names
cim_slap_tree <- ggtree(slap_trim, branch.length = "none") %<+% slap_lab +
  geom_tiplab(aes(label = Protein2))

cim_slap_motif <- gheatmap(cim_slap_tree, motifmatrix, offset=0.2, width=1,
                           colnames = TRUE, legend_title="substitution", color = "black",
                           colnames_angle = 45, hjust=1) +
  #scale_x_ggtree() +
  #scale_y_continuous(expand=c(0, 0.3)) +
  scale_fill_manual(values = c("lightpink",
                               "red", "white")) +
  NULL
cim_slap_motif

cim_res <- gheatmap(cim_slap_tree, resmatrix, offset=2, width=1,
                     colnames = TRUE, legend_title="substitution", color = "black",
                     colnames_angle = 45, hjust=1) +
  scale_fill_viridis_d(option = "H") +
  NULL
cim_res

cim_res_mat <- gheatmap(cim_res, motifmatrix, offset=9, width=1,
         colnames = TRUE, legend_title="substitution", color = "black",
         colnames_angle = 45, hjust=1) +
  scale_fill_viridis_d(option = "H") +
  NULL
cim_res_mat

gheatmap(cim_res_mat, zn_max, offset=16, width=0.25,
         colnames = TRUE, legend_title="binding", color = "black",
         colnames_angle = 45, hjust=1) +
  scale_fill_viridis_d(option = "H") +
  # scale_fill_manual(values = c("pink", "red",
  #                              "red", "white", "white", "grey")) +
  NULL
#ggsave("plots/male_rt/SLAP_PLOTS/CIM_SLAPs_COMP.pdf", width = 10, height = 4)


# > ranked abundance ####

# data frame of cimex proteins with HOGs and Dmel orthologs
dmel_ranked_abundance <- wide_mns %>% select(1:5) %>% # mean values for each tissue/sperm
  left_join(male_cats_long, by = c("Protein" = "accession")) %>% # add bias category data
  # add HOGs
  left_join(hogs_long %>% filter(species == "Cimex"), by = "Protein") %>%
  left_join(hogs, by = c("HOG", "OG", "Gene Tree Parent Clade")) %>%
  select(-c(`Gene Tree Parent Clade`:species)) %>%
  # add Dmel gene IDs
  left_join(flybase2uniprot %>% drop_na(HOG) %>%
              select(gene_name, HOG) %>%
              distinct(gene_name, .keep_all = TRUE) %>%
              group_by(HOG) %>% summarise(gene_name = paste(gene_name, collapse = ", "))) %>%
  # sort data in descending rank order abundance for each tissue/sperm
  pivot_longer(cols = 2:5) %>%
  group_by(Protein) %>%
  dplyr::slice(which.max(value)) %>%
  ungroup() %>%
  group_by(name) %>%
  arrange(desc(value)) %>%
  mutate(rnk = percent_rank(value)) %>%
  ungroup() %>%

  left_join(flybase2uniprot %>% drop_na(HOG) %>% select(FBgn, HOG), by = "HOG") %>%
  distinct(Protein, .keep_all = TRUE) %>%
  mutate(gene_name = case_when(Protein %in% ncbi_slaps$Protein[ncbi_slaps$SLAP_group == "1"] ~ "S-Lap-g1",
                               Protein %in% ncbi_slaps$Protein[ncbi_slaps$SLAP_group == "2"] ~ "S-Lap-g2",
                               TRUE ~ gene_name),
         dmel_prot = case_when(FBgn %in% wigbySFP$FBgn ~ "sfp",
                               FBgn %in% DmSP3$FBgn ~ "DmSP3"))

# look for Y linked or kl proteins
dmel_ranked_abundance[grepl("Y", x = dmel_ranked_abundance$gene_name), "gene_name"]
dmel_ranked_abundance[grepl("kl", x = dmel_ranked_abundance$gene_name), "gene_name"]
dmel_ranked_abundance[grepl("Arg", x = dmel_ranked_abundance$gene_name), "gene_name"]

# make unique names for all the SLAP orthologs which currently only have which group (1/2) they belong to
cimslap_id <- dmel_ranked_abundance[dmel_ranked_abundance$Protein %in% ncbi_slaps$Protein, ] %>%
  #group_by(gene_name) %>%
  mutate(gn = paste(gene_name, row_number(), sep = "_")) %>% pull(gn)

dmel_ranked_abundance[dmel_ranked_abundance$Protein %in% ncbi_slaps$Protein, "gene_name"] <- cimslap_id

# rename the S-Laps with 1-to-1 orthologs
#dmel_ranked_abundance[dmel_ranked_abundance$Protein %in% RBB$Cimex_ID, ]
dmel_ranked_abundance$gene_name[dmel_ranked_abundance$Protein %in% RBB$Cimex_ID] <- RBB$gene_name[c(4, 2, 3, 1)]

# which category/bias are SLAPs in?
dmel_ranked_abundance[dmel_ranked_abundance$Protein %in% ncbi_slaps$Protein, ] # all in sperm category


### plot ###
dmel_ranked_abundance %>%
  drop_na(category) %>%
  ggplot(aes(x = rnk, y = value)) +
  geom_point() +
  #scale_x_discrete(expand=expansion(add=c(rep(50, 2)))) +
  #scale_colour_manual(values = rev(Hiroshige.4)) +
  labs(x = "", y = "Ranked abundance") +
  theme_bw() +
  theme(legend.position = "",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 15)) +
  facet_wrap(~name, ncol = 4,
             labeller = as_labeller(c(mesa = "Mesadenes",
                                      sfc = "Seminal fluid containers",
                                      sperm = "Sperm",
                                      testis = "Testes"))
             ) +
  #ggsave('plots/male_rt/Ranked_Abundance.pdf', height = 3, width = 12, dpi = 600, useDingbats = FALSE)
  NULL

# top 5 % only
dmel_ranked_abundance %>% #select(Protein, category, name, gene_name, rnk, value)
  filter(rnk >= 0.95) %>%
  drop_na(category) %>%
  ggplot(aes(x = rnk, y = value, colour = category)) +
  geom_point() +
  scale_colour_manual(values = rev(Hiroshige.4)) +
  labs(x = "", y = "Ranked abundance") +
  theme_minimal() +
  theme(legend.position = "",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.text.y = element_blank()) +
  facet_wrap(~ name,
             ncol = 4,
             scales = "free",
             labeller = as_labeller(c(mesa = "Mesadenes",
                                      sfc = "Seminal fluid containers",
                                      sperm = "Sperm",
                                      testis = "Testes"))
             ) +
  ggrepel::geom_text_repel(
    data = dmel_ranked_abundance %>%
      drop_na(category) %>%
      filter(rnk >= 0.95),
    aes(label = gene_name),
    size = 2,
    colour = 'black',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    #max.overlaps = 35,
    seed = 12012026
  ) +
  #coord_cartesian(xlim = c(0.95, 1)) +
  #ggsave('plots/male_rt/Ranked_Abundance_top5.pdf', height = 3, width = 12, dpi = 600, useDingbats = FALSE)
  NULL

### table of Dmel orthologs ###
dmel_ranked_abundance %>% drop_na(gene_name) %>%
  #mutate(gene_name = if_else(is.na(dmel_prot) == TRUE, gene_name, paste(gene_name, dmel_prot, sep = "."))) %>%
  group_by(category) %>%
  arrange(desc(rnk)) %>%
  select(category, gene_name) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = category,
              values_from = gene_name) %>%
  select(-row) %>%
  select(tes, spe, sfp, ags, unbiased = `NA`) %>% #dplyr::slice(1:20) %>% write_csv("output/male_rt/topDmelorths.csv")
  print(n = 50)


# count how many of the top ranked in each group actually have orthologs
dmel_ranked_abundance %>% #drop_na(gene_name) %>%
  #mutate(gene_name = if_else(is.na(dmel_prot) == TRUE, gene_name, paste(gene_name, dmel_prot, sep = "."))) %>%
  group_by(category) %>%
  arrange(desc(rnk)) %>%
  select(category, gene_name) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = category,
              #names_sep = ".",
              values_from = gene_name) %>%
  select(-row) %>%
  select(tes, spe, sfp, ags, unbiased = `NA`) %>%
  dplyr::slice(1:20)



# >> count number of sfps/sperm orthologs ####
# mel counts - sperm proteins
dmel_ranked_abundance %>%
  filter(category == "spe") %>%
  filter(HOG %in% DmSP3$HOG) %>%
  select(Protein, rnk, Drosophila, gene_name) #%>% write_csv("output/male_rt/all_spe_orths.csv")

# mel counts - sfps
dmel_ranked_abundance %>%
  filter(category == "sfp") %>%
  filter(HOG %in% wigbySFP$HOG) %>%
  select(Protein, rnk, Drosophila, gene_name) #%>% write_csv("output/male_rt/all_sfp_orths.csv")

# aedes counts - sperm
dmel_ranked_abundance %>%
  filter(category == "spe") %>%
  filter(HOG %in% aedes_sperm$HOG) %>%
  select(Protein, HOG, Aedes, Dmel_gene = gene_name)

# aedes counts - sfps
dmel_ranked_abundance %>%
  filter(category == "sfp") %>%
  filter(HOG %in% aedes_sfps$HOG) %>%
  select(Protein, HOG, Aedes, Dmel_gene = gene_name)


# evolutionary rates of S-Laps
one_ratio %>% filter(accession %in% sperm_proteome$Protein) %>%
  ggplot(aes(x = omega)) +
  geom_histogram(alpha = .5) +
  geom_point(data = one_ratio %>%
               inner_join(dmel_ranked_abundance %>%
                            filter(HOG %in% cimex_slaps$HOG) %>%
                            select(accession = Protein, gene_name)),
             aes(y = 1),
             colour = "red") +
  ggrepel::geom_text_repel(data = one_ratio %>%
                             inner_join(dmel_ranked_abundance %>%
                                          filter(HOG %in% cimex_slaps$HOG) %>%
                                          select(accession = Protein, gene_name)),
    aes(y = 1, label = gene_name),
    size = 3,
    colour = 'black',
    direction = "y",
    hjust = 0,
    force        = 10,
    nudge_x      = 0.5,
    max.overlaps = Inf,
    seed = 12012026
  ) +
  theme_bw() +
  NULL
#ggsave('plots/male_rt/SLAP_PLOTS/SLAP_RATES.pdf', height = 4, width = 4, dpi = 600, useDingbats = FALSE)



# > Endosymbiont proteins ####
# have a deeper look at all the endosymbiont proteins before filtering them out
endosymbiont_all <- ms_data %>% filter(accession %in% endosymbiont_prots$accession) %>%
  dplyr::select(accession, number_unique_peptides,
                abundance_f1_sample_ags:abundance_f17_sample_testis)
colnames(endosymbiont_all)[-c(1:2)] <- gsub(".*_", "", x = colnames(endosymbiont_all)[-c(1:2)])
colnames(endosymbiont_all)[-c(1:2)] <- paste(colnames(endosymbiont_all)[-c(1:2)],
                                unlist(lapply(rle(colnames(endosymbiont_all))$lengths[-c(1:2)],
                                              function(x) seq(from = 1, to = x, by = 1))), sep = "_")

colnames(endosymbiont_all)[5:7] <- c("ags_3.1", "ags_3.2", "ags_4")
endosymbiont_all %>% filter(number_unique_peptides >= 2) %>% dim

endosymbiont_ms <- endosymbiont_all %>% #filter(number_unique_peptides >= 2) %>%
  rowwise() %>%
  # take the average of the 2 reps for ags_3
  mutate(ags_3 = mean(c(ags_3.1, ags_3.2))) %>%
  #mutate(ags_3 = ags_3.2) %>%
  dplyr::select(-c(ags_3.1, ags_3.2)) %>%
  relocate(ags_3, .before = ags_4) %>%
  dplyr::select(-2) %>% dplyr::rename(ProteinName = accession) %>%
  pivot_longer(cols = 2:17) %>%
  group_by(name) %>% mutate(Run = cur_group_id()) %>% ungroup() %>%
  mutate(PeptideSequence = ProteinName,
         PrecursorCharge = NA,
         FragmentIon = NA,
         ProductCharge = NA,
         IsotopeLabelType = "L",
         Condition = str_sub(name, 1, end = -3),
         BioReplicate = name,
         Intensity = value) %>%
  dplyr::select(-name, -value) %>% relocate(Run, .after = BioReplicate) %>% mutate(across(2:9, as.factor)) %>%
  data.frame()

# >>> MSstats ####
if(!file.exists("output/male_rt_endosymbiont.rds")){

  male_rt_endosymbiont <- dataProcess(raw = endosymbiont_ms,
                                   normalization = 'quantile',
                                   summaryMethod = 'TMP',
                                   censoredInt = "NA",
                                   MBimpute = TRUE,
                                   maxQuantileforCensored = 0.999,
                                   use_log_file = FALSE)

  saveRDS(male_rt_endosymbiont, "output/male_rt_endosymbiont.rds")
} else {
  male_rt_endosymbiont <- readRDS('output/male_rt_endosymbiont.rds')
}

# some proteins are dropped at this processing stage

# >>> Normalised data ####
norm_endosymbiont <- quantification(male_rt_endosymbiont)
colnames(norm_endosymbiont)[-1] <- gsub("^[A-Za-z1-9]+_", "", x = colnames(norm_endosymbiont)[-1])

norm_endosymbiont <- norm_endosymbiont %>%
  mutate(spp = case_when(Protein %in% endosymbiont_prots$accession[grepl("Wolb", endosymbiont_prots$description)] ~ "Wolbachia",
                         Protein %in% endosymbiont_prots$accession[grepl("Rick", endosymbiont_prots$description)] ~ "Rickettsia",
                         TRUE ~ NA))

norm_endosymbiont %>% dplyr::count(spp)

# do tests
# we can reuse the design matrix from above
if(!file.exists("output/male_rt_endo_comparisons.rds")){

  male_rt_endo_comparisons <- groupComparison(contrast.matrix = all_comps, data = male_rt_endosymbiont)

  saveRDS(male_rt_endo_comparisons, "output/male_rt_endo_comparisons.rds")
} else {
  male_rt_endo_comparisons <- readRDS('output/male_rt_endo_comparisons.rds')
}

da_endo <- male_rt_endo_comparisons$ComparisonResult %>%
  mutate(DA = if_else(abs(log2FC) > 1 & adj.pvalue < 0.05, "SD", "NS")) %>%
  drop_na(DA)

da_endo %>%
  ggplot(aes(x = log2FC, y = -log10(pvalue), colour = DA)) +
  geom_point() +
  facet_wrap(~ Label)

# identify proteins most abundant in each tissue.
ags.endo <- da_endo %>% filter(str_detect(Label, "ags"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`ags-mesa` > 1, `ags-sperm` > 1, `ags-testis` > 1)

mes.endo <- da_endo %>% filter(str_detect(Label, "mesa"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`ags-mesa` < -1, `mesa-sperm` > 1, `mesa-testis` > 1)

spe.endo <- da_endo %>% filter(str_detect(Label, "sperm"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`ags-sperm` < -1, `mesa-sperm` < -1, `sperm-testis` > 1)

tes.endo <- da_endo %>% filter(str_detect(Label, "testis"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`ags-testis` < -1, `mesa-testis` < -1, `sperm-testis` < -1)

# proteins more abundant in AGS/MES compared to sperm/testes
sem.endo <- da_endo %>% filter(str_detect(Label, "a"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`ags-sperm` > 1, `ags-testis` > 1, `mesa-sperm` > 1, `mesa-testis` > 1) %>%
  filter(!Protein %in% c(ags.endo$Protein, mes.endo$Protein))

# proteins more abundant in sperm/testis compared to ags/mes, i.e. could be immature sperm etc. in germline
ger.endo <- da_endo %>% filter(str_detect(Label, "sperm|testis"), adj.pvalue < 0.05) %>%
  select(Protein, Label, log2FC) %>%
  pivot_wider(names_from = Label, values_from = log2FC) %>%
  filter(`ags-sperm` < -1, `ags-testis` < -1, `mesa-sperm` < -1, `mesa-testis` < -1) %>%
  filter(!Protein %in% c(spe.endo$Protein, tes.endo$Protein))

tissue_biased_endo <- bind_rows(ags.endo %>% select(Protein) %>% mutate(bias = "ags"),
                                mes.endo %>% select(Protein) %>% mutate(bias = "mes"),
                                spe.endo %>% select(Protein) %>% mutate(bias = "spe"),
                                tes.endo %>% select(Protein) %>% mutate(bias = "tes"),

                                sem.endo %>% select(Protein) %>% mutate(bias = "sem"),
                                ger.endo %>% select(Protein) %>% mutate(bias = "ger"))

# >>> Heatmap ####
hm_endo <- pheatmap:::scale_rows(norm_endosymbiont %>% dplyr::select(2:17) %>%
                                   mutate(across(everything(), ~replace_na(.x, 0)))) %>%
  mutate(Protein = norm_endosymbiont$Protein) %>%
  relocate(Protein, .before = ags_1) %>%
  # this protein has no values
  filter(Protein != "A0A8K1J0Q4") %>%
  mutate(spp = case_when(Protein %in% endosymbiont_prots$accession[grepl("Wolb", endosymbiont_prots$description)] ~ "Wolbachia",
                         Protein %in% endosymbiont_prots$accession[grepl("Rick", endosymbiont_prots$description)] ~ "Rickettsia",
                         TRUE ~ NA)) %>%
  left_join(tissue_biased_endo)


wolbachia_lab <- rowAnnotation(
  spp = hm_endo$spp,
  #bias = hm_endo$bias,
  col = list(bias = setNames(Hiroshige.6,
                             #rainbow(n_distinct(na.omit(bias.frame$bias))),
                             levels(as.factor(hm_endo$bias))),
             spp = c('Rickettsia' = "orange",
                     'Wolbachia' = "turquoise")
  ),
  na_col = NA,
  title = NULL,
  #                          show_legend = FALSE,
  show_annotation_name = FALSE)


# draw a heatmap
#pdf('plots/male_rt/endosymbiont_hm.pdf', height = 6, width = 5)
Heatmap(as.matrix(hm_endo %>%
                    # only Wolbachia proteins
                    filter(Protein %in% endosymbiont_prots$accession[grepl("Wolb", endosymbiont_prots$description)]) %>%
                    select(-Protein, -spp, -bias)),
        #col = viridis::inferno(50),
        na_col = NA,
        heatmap_legend_param = list(title = "log2(intensities)",
                                    direction = "horizontal"),
        #left_annotation = wolbachia_lab,
        #right_annotation = right_lab,
        top_annotation = vmha,
        show_row_names = FALSE,
        show_column_names = FALSE,
        #show_row_dend = FALSE,
        #row_order = order(hm_endo$spp),
        clustering_method_rows = "complete",
        column_split = 4,
        column_gap = unit(0, "mm"),
        row_title = NULL,
        #row_split = 4,
        column_title = NULL)
#dev.off()


# > Alphafold picks ####

# all Wolbachia proteins (after filtering) for AlphaFold
norm_endosymbiont %>% filter(spp == "Wolbachia") %>% select(Protein) #%>% write_csv("output/male_rt/alphafold_stuff/cimex_wolbachia_110.csv")

