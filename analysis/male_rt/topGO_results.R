library(tidyverse)
library(topGO)

# DAVID uniprot ID to entrez_geneID
uniprot2entrez <- read.delim("../Bedbug_RNAseq/DAVID2entrez_geneID.txt") %>%
  mutate(gene_id = paste0("LOC", To)) %>%
  dplyr::rename(accession = From) #%>% left_join(read_csv("output/male_rt/Cimex2Dmel.csv"), by = c("accession" = "Protein_Clec"))

# load GO ids for each gene
go_ids <- read_rds("output/cimex_go_ids.rds") %>%
  left_join(uniprot2entrez, by = c("Accession" = "accession"))

# a small number of genes dont have uniprot IDs associated
go_ids %>% filter(is.na(gene_id))

golist = strsplit(as.character(go_ids$GOID), split = ';')
names(golist) = go_ids$gene_id

go_cats <- c("BP", "CC", "MF")


# > Male proteins ####

# >> male tissue-biased proteins ####
male_cats <- read_rds("output/male_rt/male_cats.rds")

# background proteome - all proteins
go_ids <- read_tsv("../proteomes/uniprot_cimex/uniprotkb_taxonomy_id_79782_2023_08_09.tsv") %>%
  dplyr::select(Accession = Entry, GOID = `Gene Ontology IDs`) %>%
  mutate(GOID = gsub(" ", "", x = GOID))

golist = strsplit(as.character(go_ids$GOID), split = ';')
names(golist) = go_ids$Accession

all_lists <- list(sfp_genes = factor(as.integer(as.logical(go_ids$Accession %in% male_cats$sfp))),
                  ags_genes = factor(as.integer(as.logical(go_ids$Accession %in% male_cats$ags))),
                  spe_genes = factor(as.integer(as.logical(go_ids$Accession %in% male_cats$spe))),
                  tes_genes = factor(as.integer(as.logical(go_ids$Accession %in% male_cats$tes))))

names(all_lists[[1]]) <- go_ids$Accession
names(all_lists[[2]]) <- go_ids$Accession
names(all_lists[[3]]) <- go_ids$Accession
names(all_lists[[4]]) <- go_ids$Accession



GO_results <- list()
k <- 1

for(i in 1:n_distinct(go_cats)) {

  for(j in 1:n_distinct(all_lists)) {

    GOdata = new('topGOdata', ontology = go_cats[i], allGenes = all_lists[[j]], nodeSize = 15,
                 annot = annFUN.gene2GO, gene2GO = golist)

    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

    g2 = GenTable(GOdata, classic = resultFisher, orderBy = "weight", ranksOf = "classic", topNodes = 10)
    # correct p-value for multiple testing
    p.adj = p.adjust(g2$classic, method = "fdr")
    sig_results = data.frame(g2, p.adj = p.adj) %>% mutate(GO.cat = go_cats[i])

    print(sig_results)

    # save results
    #write_csv(sig_results, paste0("output/topGO_results/male_rt/", names(all_lists[j]), "_biased_", go_cats[i], ".csv"))

    # save results to object
    GO_results[[k]] <- data.frame(.id = names(all_lists[j]), sig_results)
    k <- k + 1

  }

}

GO_results_comb <- do.call(rbind, GO_results)

# save results as 1 file
# write_csv(GO_results_comb, "output/topGO_results/male_rt/male_cat_results.csv")

GO_results_comb %>% filter(.id == "sfp_genes", classic < 0.05)
#GO_results_comb %>% filter(.id == "sfp_genes", p.adj < 0.05)

# plot results
facet_names <- c(BP = "Biological Process",
                 CC = "Cellular Component",
                 MF = "Molecular Function")

# look at each tissue separately
GO_results_comb %>%
  #filter(Significant >= 5) %>%
  filter(GO.cat == "BP") %>% #print(n = Inf)
  ggplot(aes(x = reorder(Term, p.adj), y = -log10(p.adj), colour = GO.cat)) +
  geom_segment(aes(x = reorder(Term, p.adj), xend = reorder(Term, p.adj), y = 0, yend = -log10(p.adj)),
               colour = "black") +
  geom_point(size = 7) +
  scale_colour_brewer(palette = 'Dark2') +
  facet_wrap( ~ .id, scales = 'free_x', labeller = as_labeller(facet_names), nrow = 1) +
  labs(y = '-log10(p-value)') +
  theme_bw() +
  theme(legend.position = 'none',# = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  geom_text(aes(label = Significant), colour = 'black', size = 3) +
  #ggsave('plots/topGO_figs/msl_DAP_all_GO.pdf', height = 8, width = 20, dpi = 600, useDingbats = FALSE) +
  NULL




# >> Wolbachia proteins ####
wClec_go_ids <- read_tsv("../uniprotkb_UP000031663_wClec_2025_06_11.tsv") %>%
  dplyr::select(Accession = Entry, GOID = `Gene Ontology IDs`) %>%
  mutate(GOID = gsub(" ", "", x = GOID))

wClec_golist <- strsplit(as.character(wClec_go_ids$GOID), split = ';')
names(wClec_golist) = wClec_go_ids$Accession


wClec_prots <- norm_endosymbiont %>%
  filter(Protein %in% endosymbiont_prots$accession[grepl("Wolb", endosymbiont_prots$description)]) %>%
  pull(Protein)


wClec_prots <- factor(as.integer(as.logical(wClec_go_ids$Accession %in% wClec_prots)))
names(wClec_prots) <- wClec_go_ids$Accession


go_cats <- c("BP", "CC", "MF")

# for(i in 1:n_distinct(go_cats)) {
#
#   GOdata = new('topGOdata', ontology = c("BP", "CC", "MF")[i], allGenes = wClec_prots, nodeSize = 15,
#                annot = annFUN.gene2GO, gene2GO = wClec_golist)
#
#   resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
#
#   g2 = GenTable(GOdata, classic = resultFisher, orderBy = "weight", ranksOf = "classic", topNodes = 10)
#   # correct p-value for multiple testing
#   p.adj = p.adjust(g2$classic, method = "fdr")
#   sig_results = data.frame(g2, p.adj = p.adj) %>% mutate(GO.cat = go_cats[i])
#   print(sig_results)
#   # change this for each gene list
#   #write_csv(sig_results, paste0("output/topGO_results/wClec_110_", go_cats[i], ".csv"))
#
# }

# plot results
facet_names <- c(BP = "Biological Process",
                 CC = "Cellular Component",
                 MF = "Molecular Function")

bind_rows(read_csv("output/topGO_results/wClec_110_BP.csv") %>% mutate(GO = "BP"),
          read_csv("output/topGO_results/wClec_110_CC.csv") %>% mutate(GO = "CC"),
          read_csv("output/topGO_results/wClec_110_MF.csv") %>% mutate(GO = "MF")) %>% #print(n = Inf)
  #filter(Significant >= 5) %>%
  ggplot(aes(x = reorder(Term, p.adj), y = -log10(p.adj), colour = GO)) +
  geom_segment(aes(x = reorder(Term, p.adj), xend = reorder(Term, p.adj), y = 0, yend = -log10(p.adj)),
               colour = "black") +
  geom_point(size = 7) +
  scale_colour_brewer(palette = 'Dark2') +
  facet_wrap(~ GO, scales = 'free_x', ncol = 3, labeller = as_labeller(facet_names)) +
  labs(y = '-log10(p-value)') +
  theme_bw() +
  theme(legend.position = 'none',# = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  geom_text(aes(label = Significant), colour = 'black', size = 3) +
  #ggsave('plots/topGO_figs/msl_DAP_all_GO.pdf', height = 8, width = 20, dpi = 600, useDingbats = FALSE) +
  NULL



# > switch analysis ####
names(golist) = go_ids$Accession

ejac_dap <- read.csv("output/topGO_lists/ejac_DAP_all.csv") %>%
  distinct(effect, Protein) %>% filter(Protein %in% pellet_ejac$Protein)

ejac_dap_split <- split(ejac_dap$Protein, ejac_dap$effect)

# change this for each gene list
ejac_genes <- factor(as.integer(as.logical(go_ids$Accession %in% ejac_dap_split$e)))
names(ejac_genes) <- go_ids$Accession


go_cats <- c("BP", "CC", "MF")

for(i in 1:n_distinct(go_cats)) {

  GOdata = new('topGOdata', ontology = c("BP", "CC", "MF")[i], allGenes = ejac_genes, nodeSize = 15,
               annot = annFUN.gene2GO, gene2GO = golist)

  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

  g2 = GenTable(GOdata, classic = resultFisher, orderBy = "weight", ranksOf = "classic", topNodes = 10)
  # correct p-value for multiple testing
  p.adj = p.adjust(g2$classic, method = "fdr")
  sig_results = data.frame(g2, p.adj = p.adj) %>% mutate(GO.cat = go_cats[i])
  print(sig_results)
  # change this for each gene list
  #write_csv(sig_results, paste0("output/topGO_results/ejac_DAP_all_pellet_e_", go_cats[i], ".csv"))

}


# plot results
facet_names <- c(BP = "Biological Process",
                 CC = "Cellular Component",
                 MF = "Molecular Function")

bind_rows(read_csv("output/topGO_results/pellet_DAP_e_BP.csv") %>% mutate(GO = "BP"),
          read_csv("output/topGO_results/pellet_DAP_e_CC.csv") %>% mutate(GO = "CC"),
          read_csv("output/topGO_results/pellet_DAP_e_MF.csv") %>% mutate(GO = "MF")) %>% #print(n = Inf)
  #filter(Significant >= 5) %>%
  ggplot(aes(x = reorder(Term, p.adj), y = -log10(p.adj), colour = GO)) +
  geom_segment(aes(x = reorder(Term, p.adj), xend = reorder(Term, p.adj), y = 0, yend = -log10(p.adj)),
               colour = "black") +
  geom_point(size = 7) +
  scale_colour_brewer(palette = 'Dark2') +
  facet_wrap(~ GO, scales = 'free_x', ncol = 3, labeller = as_labeller(facet_names)) +
  labs(y = '-log10(p-value)') +
  theme_bw() +
  theme(legend.position = 'none',# = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  geom_text(aes(label = Significant), colour = 'black', size = 3) +
  #ggsave('plots/topGO_figs/msl_DAP_all_GO.pdf', height = 8, width = 20, dpi = 600, useDingbats = FALSE) +
  NULL



all_bp <- bind_rows(read_csv("output/topGO_results/ag_bias_BP.csv") %>% mutate(tissue = "AG"),
                    read_csv("output/topGO_results/testes_bias_BP.csv") %>% mutate(tissue = "testes"),
                    read_csv("output/topGO_results/sperm_bias_BP.csv") %>% mutate(tissue = "sperm"),
                    read_csv("output/topGO_results/mesadenes_bias_BP.csv") %>% mutate(tissue = "mesa"))

level_order <- all_bp$Term[all_bp$Significant >=5]

all_bp %>%
  filter(Significant >= 5) %>%
  ggplot(aes(x = reorder(Term, p.adj), y = -log10(p.adj), colour = tissue)) +
  geom_segment(aes(x = reorder(Term, p.adj), xend = reorder(Term, p.adj), y = 0, yend = -log10(p.adj)),
               colour = "black") +
  geom_point(size = 7) +
  scale_colour_manual(values = Hiroshige.pal, labels = c("SFCs", "Mesadenes", "Sperm", "Testes")) +
  scale_x_discrete(limits = level_order) +
  # facet_wrap(~ tissue, scales = 'free_x', ncol = 4,
  #            labeller = as_labeller(c(testes = "Testes", AG = "Accessory glands",
  #                                     sperm = "Sperm", mesa = "Mesadenes"))) +
  labs(y = '-log10(p-value)') +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  geom_text(aes(label = Significant), colour = 'black', size = 3) +
  #ggsave('plots/topGO_figs/all_msl_GO.pdf', height = 26, width = 8, dpi = 600, useDingbats = FALSE) +
  NULL



library(clusterprofiler)
library(GO.db)

GO_tau_comb <- read_csv("output/rna_seq_results/topGO_results/topGO_combined_tau.csv")
parents <- mget(GO_tau_comb$GO.ID[GO_tau_comb$GO.cat == "BP"], GOBPPARENTS)

library(AnnotationDbi)

# summarise_go_terms <- function(go_ids, ontologies = c("BP","MF","CC"), max_terms = 10) {
#
#   # match ontology argument
#   ontologies <- match.arg(ontologies, several.ok = TRUE)
#
#   # ancestor database
#   ancestor_db <- list(
#     BP = GOBPANCESTOR,
#     MF = GOMFANCESTOR,
#     CC = GOCCANCESTOR
#   )
#
#   root_terms <- c(BP="GO:0008150", MF="GO:0003674", CC="GO:0005575")
#
#   results <- list()
#
#   for(ont in ontologies){
#     # filter GO IDs for this ontology
#     # fallback: we try to keep only IDs in GO.db
#     go_in_ont <- go_ids[sapply(go_ids, function(x) !is.null(ancestor_db[[ont]][[x]]))]
#     if(length(go_in_ont)==0) next
#
#     # get all ancestors (excluding root)
#     all_ancestors <- unlist(lapply(go_in_ont, function(go){
#       anc <- ancestor_db[[ont]][[go]]
#       anc[anc != root_terms[ont]]  # drop root
#     }))
#
#     # count frequency of each ancestor
#     anc_counts <- sort(table(all_ancestors), decreasing = TRUE)
#
#     if(length(anc_counts)==0) next
#
#     # select top N most common ancestors
#     top_n <- head(anc_counts, max_terms)
#
#     go_ids_sel <- names(top_n)
#     counts <- as.integer(top_n)
#
#     # get GO term names
#     term_names <- Term(go_ids_sel)
#
#     # store in data frame
#     df <- data.frame(
#       GO_ID = go_ids_sel,
#       Term = term_names,
#       Count = counts,
#       Ontology = ont,
#       stringsAsFactors = FALSE
#     )
#
#     results[[ont]] <- df
#   }
#
#   # combine all ontologies
#   summary_df <- do.call(rbind, results)
#
#   return(summary_df)
# }
#summarise_go_terms(GO_tau_comb$GO.ID[GO_tau_comb$.id == "m_genes"], max_terms = 15)


summarise_go_terms <- function(go_ids,
                               ontologies = c("BP","MF","CC"),
                               max_terms = 10,
                               min_depth = 3,
                               max_cover = 0.7,
                               drop_generic_terms = NULL) {

  result_list <- list()

  for (ont in ontologies) {

    # Select ancestor database
    parents_db <- switch(ont,
                         BP = GOBPANCESTOR,
                         MF = GOMFANCESTOR,
                         CC = GOCCANCESTOR)

    # Keep only GO IDs present in this ontology
    go_in_ont <- go_ids[sapply(go_ids, function(x) x %in% keys(parents_db))]
    if(length(go_in_ont) == 0) next

    # Get ancestors for each GO ID
    ancestors_list <- lapply(go_in_ont, function(go) {
      anc <- parents_db[[go]]
      if(!is.null(anc)) anc <- anc[anc != "all"]  # remove placeholder
      return(anc)
    })

    all_ancestors <- unlist(ancestors_list)
    if(length(all_ancestors) == 0) next

    # Count frequency
    anc_counts <- sort(table(all_ancestors), decreasing = TRUE)

    # Remove generic terms if supplied
    if(!is.null(drop_generic_terms)) {
      anc_counts <- anc_counts[!names(anc_counts) %in% drop_generic_terms]
    }

    # Remove terms that are shared by too many GO IDs
    anc_counts <- anc_counts[anc_counts <= max_cover * length(go_in_ont)]
    if(length(anc_counts) == 0) next

    # Get GO term names using AnnotationDbi::select
    term_info <- AnnotationDbi::select(GO.db,
                                       keys = names(anc_counts),
                                       columns = c("TERM","ONTOLOGY","GOID"),
                                       keytype = "GOID")

    # Filter by min_depth using term hierarchy levels
    # Depth approximated by the number of ancestors (higher = more specific)
    depth_approx <- sapply(names(anc_counts), function(go) length(parents_db[[go]]))
    keep_terms <- names(anc_counts)[depth_approx >= min_depth]
    if(length(keep_terms) == 0) next

    top_anc <- head(anc_counts[keep_terms], max_terms)

    # Merge with term_info
    df <- term_info[match(names(top_anc), term_info$GOID), c("GOID","TERM")]
    df$Count <- as.integer(top_anc)
    df$Depth <- depth_approx[match(names(top_anc), names(depth_approx))]
    df$Ontology <- ont

    result_list[[ont]] <- df
  }

  # Combine all ontologies
  if(length(result_list) == 0) return(data.frame())
  do.call(rbind, result_list)
}

drop_generic <- c("GO:0008150", "GO:0003674", "GO:0005575",
                  "GO:0009987", "GO:0008152", "GO:0003824")  # super-generic terms

summary_msl <- summarise_go_terms(
  go_ids = GO_tau_comb$GO.ID[GO_tau_comb$.id == "m_genes"],
  max_terms = 10,
  min_depth = 7,
  max_cover = 0.7,
  drop_generic_terms = drop_generic
)

summary_msl





# >> ClusterProfiler Compare Cimex to ... ####

# >>> Sperm proteins
cimexsperm_gos <- go_ids %>% filter(Accession %in% sperm_proteome$Protein)
cimexsperm_golist <- strsplit(as.character(cimexsperm_gos$GOID), split = ';')
names(cimexsperm_golist) <- cimexsfp_gos$Accession

clec_sperm_go <- c(unlist(cimexsperm_golist)) %>% na.omit()
dmel_sperm_go <- c(unlist(DmSP3_golist)) %>% na.omit()

# >>> Sfps
cimexsfp_gos <- go_ids %>% filter(Accession %in% unique(c(sfc.list$Protein,
                                                          mes.list$Protein,
                                                          sfp.list$Protein)),
                                  Accession %in% signal_peps)
cimexsfp_golist <- strsplit(as.character(cimexsfp_gos$GOID), split = ';')
names(cimexsfp_golist) <- cimexsfp_gos$Accession

clec_sfp_go <- c(unlist(cimexsfp_golist)) %>% na.omit()
dmel_sfp_go <- c(unlist(wigbysfp_golist)) %>% na.omit()




# go1 = c("GO:0004022","GO:0004024","GO:0004174")
# go2 = c("GO:0009055","GO:0005515")
# d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
# mgoSim(go1, go2, semData=d, measure="Wang", combine=NULL)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("org.Dm.eg.db")

dmGO.BP <- godata(annoDb = 'org.Dm.eg.db', ont="BP")
dmGO.CC <- godata(annoDb = 'org.Dm.eg.db', ont="CC")
dmGO.MF <- godata(annoDb = 'org.Dm.eg.db', ont="MF")

str(dmGO.BP)

# bp.sfp <- mgoSim(clec_sfp_go, dmel_sfp_go,
#                  semData=dmGO.BP, measure="Wang", combine="BMA")

#write_rds(bp.sfp, "output/male_rt/gosemsim_bp.sfp.rds")
bp.sfp <- read_rds("output/male_rt/gosemsim_bp.sfp.rds")

str(bp.sfp)

# bp.sperm <- mgoSim(clec_sperm_go, dmel_sperm_go,
#                    semData=dmGO.BP, measure="Wang", combine="BMA")
#write_rds(bp.sperm, "output/male_rt/gosemsim_bp.sperm.rds")
bp.sperm <- read_rds("output/male_rt/gosemsim_bp.sperm.rds")


# get a random subset of genes

clec_x = go_ids %>% sample_n(size = n_distinct(cimexsfp_gos$Accession))

clec_y = strsplit(as.character(clec_x$GOID), split = ';')
names(clec_y) = clec_x$Accession
clec_z = c(unlist(clec_y)) %>% na.omit()

dmel_x = flybase2uniprot %>%
  sample_n(size = n_distinct(wigbySFP$FBgn)) %>%
  mutate(GOID = gsub(" ", "", x = `Gene Ontology IDs`))

dmel_y = strsplit(as.character(dmel_x$GOID), split = ';')
names(dmel_y) = dmel_x$FBgn
dmel_z = c(unlist(dmel_y)) %>% na.omit()

bp.rep = mgoSim(clec_z, dmel_z,
                semData=dmGO.BP, measure="Wang", combine="BMA")

bp.rep


# Create a vector to store all 1000 similarity scores
similarity_scores <- numeric(1000)

# # Run the process 1000 times
# for (i in 1:1000) {
#   # Sample clec set
#   clec_x = go_ids %>%
#     sample_n(size = n_distinct(cimexsfp_gos$Accession), replace = FALSE)
#
#   clec_y = strsplit(as.character(clec_x$GOID), split = ';')
#   names(clec_y) = clec_x$Accession
#   clec_z = unlist(clec_y) %>% na.omit()
#
#   # Sample dmel set
#   dmel_x = flybase2uniprot %>%
#     sample_n(size = n_distinct(wigbySFP$FBgn), replace = FALSE) %>%
#     mutate(GOID = gsub(" ", "", x = `Gene Ontology IDs`))
#
#   dmel_y = strsplit(as.character(dmel_x$GOID), split = ';')
#   names(dmel_y) = dmel_x$FBgn
#   dmel_z = unlist(dmel_y) %>% na.omit()
#
#   # Compute semantic similarity (Biological Process)
#   bp.rep = mgoSim(clec_z, dmel_z,
#                    semData = dmGO.BP,
#                    measure = "Wang",
#                    combine = "BMA")
#
#   # Save the result
#   similarity_scores[i] = bp.rep
#   print(i)
# }

#saveRDS(similarity_scores, file = "output/male_rt/similarity_scores_sfps1000.rds")
similarity_scores_sfp <- readRDS("output/male_rt/similarity_scores_sfps1000.rds")

summary(similarity_scores_sfp)

data.frame(draw = similarity_scores_sfp) %>%
  ggplot(aes(x = draw)) +
  #geom_histogram(bins = 100) +
  tidybayes::stat_halfeye(aes(fill = after_stat(x >= bp.sfp)), slab_colour = "black", alpha = .5) +
  geom_vline(xintercept = bp.sfp, lty = 2, colour = "red") +
  scale_fill_manual(values = c("white", "red")) +
  labs(x = "Semantic similarity score", y = "Density") +
  theme_bw() +
  theme(legend.position = "") +
  #annotate(geom = "text", x = 25, y = 0.5, label = paste0(expression(italic('p = ')), p_draws), size = 3) +
  # annotate("text", x = 30, y = 0.85, size = 5,
  #          label = paste0("italic(p)==", round(p_draws, 3)),
  #          parse = TRUE) +
  #ggsave('plots/male_rt/GOSS_sfp.pdf', height = 4, width = 4, dpi = 600, useDingbats = FALSE) +
  NULL


similarity_scores_sperm <- readRDS("output/male_rt/similarity_scores_sperm1000.rds")
summary(similarity_scores_sperm)
hist(similarity_scores_sperm)

data.frame(draw = similarity_scores_sperm) %>%
  ggplot(aes(x = draw)) +
  #geom_histogram(bins = 100) +
  tidybayes::stat_halfeye(aes(fill = after_stat(x >= bp.sperm)), slab_colour = "black", alpha = .5) +
  geom_vline(xintercept = bp.sperm, lty = 2, colour = "red") +
  scale_fill_manual(values = c("white", "red")) +
  labs(x = "Semantic similarity score", y = "Density") +
  theme_bw() +
  theme(legend.position = "") +
  #annotate(geom = "text", x = 25, y = 0.5, label = paste0(expression(italic('p = ')), p_draws), size = 3) +
  # annotate("text", x = 30, y = 0.85, size = 5,
  #          label = paste0("italic(p)==", round(p_draws, 3)),
  #          parse = TRUE) +
  #ggsave('plots/male_rt/GOSS_sfp.pdf', height = 4, width = 4, dpi = 600, useDingbats = FALSE) +
  NULL

similarity_scores_sfp <- readRDS("output/male_rt/similarity_scores_sfp1000_hpc.rds")
summary(similarity_scores_sfp)
hist(similarity_scores_sfp)

data.frame(tiss = rep(c("sfp", "sperm"), each = 1000),
           draw = c(similarity_scores_sfp, similarity_scores_sperm)) %>%
  ggplot(aes(x = tiss, y = draw)) +
  tidybayes::stat_interval(.width = c(.5, .8, .95)) +
  #tidybayes::stat_pointinterval() +
  #tidybayes::stat_pointinterval(.width = c(.66, .95), position = position_nudge(x = -0.2)) +
  scale_color_brewer() +
  labs(x = "", y = "GO semantic similarity score") +
  theme_bw() +
  theme(legend.position = "") +
  annotate("point", x = 1, y = bp.sfp, size = 3, colour = "red") +
  annotate("point", x = 2, y = bp.sperm, size = 3, colour = "red") +
  #ggsave('plots/male_rt/GOSS_sfp-sperm.pdf', height = 4, width = 4.5, dpi = 600, useDingbats = FALSE) +
  NULL


# >>> Dmel SFP/sperm protein GO terms ####

Dmel_go_ids <- flybase2uniprot %>%
  distinct(FBgn, .keep_all = TRUE) %>%
  mutate(GOID = gsub(" ", "", x = `Gene Ontology IDs`))

Dmel_golist <- strsplit(as.character(Dmel_go_ids$GOID), split = ';')
names(Dmel_golist) = Dmel_go_ids$FBgn

wigbysfp_golist


go_cats <- c("BP", "CC", "MF")

for(i in 1:n_distinct(go_cats)) {

  GOdata = new('topGOdata', ontology = c("BP", "CC", "MF")[i], allGenes = Dmel_golist, nodeSize = 15,
               annot = annFUN.gene2GO, gene2GO = wigbysfp_golist)

  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

  g2 = GenTable(GOdata, classic = resultFisher, orderBy = "weight", ranksOf = "classic", topNodes = 10)
  # correct p-value for multiple testing
  p.adj = p.adjust(g2$classic, method = "fdr")
  sig_results = data.frame(g2, p.adj = p.adj) %>% mutate(GO.cat = go_cats[i])
  print(sig_results)
  # change this for each gene list
  #write_csv(sig_results, paste0("output/topGO_results/wClec_110_", go_cats[i], ".csv"))

}

