# Load the data and search for the models
library("inteRmodel")
library("dplyr")
library("tidyr")

transcriptomics <- read.delim("data/host_tx_counts.tsv.gz", check.names = FALSE)
microbiome <- read.delim("data/taxonomic_profiles.tsv.gz", check.names = FALSE,
                         row.names = 1)
meta <- readRDS("data_out/metadata.RDS")

# To find which colum we did this:
# lm <- sapply(meta, function(x){sum(x %in% colnames(microbiome))})
# lm[lm != 0]
microbiome <- microbiome[, colnames(microbiome) %in% meta$`External ID`]

# To find which colum we did this:
# lm <- sapply(meta, function(x){sum(x %in% colnames(transcriptomics))})
# lm[lm != 0]
transcriptomics <- transcriptomics[, colnames(transcriptomics) %in% meta$`External ID`]

# Match samples and order
rnaseq_samples <- colnames(transcriptomics)
names(rnaseq_samples) <- meta$`Participant ID`[match(rnaseq_samples, meta$`External ID`)]
micro_samples <- colnames(microbiome)
names(micro_samples) <- meta$`Participant ID`[match(micro_samples, meta$`External ID`)]

m2 <- meta[meta$`External ID` %in% c(colnames(microbiome), colnames(transcriptomics)), ]

no_match_tx_bp <- m2 %>%
    group_by(site_sub_coll) %>%
    summarize(n = n_distinct(data_type)) %>%
    filter(n != 2) %>%
    pull(site_sub_coll)
m2 <- filter(m2, !site_sub_coll %in% no_match_tx_bp)
microbiome <- microbiome[, colnames(microbiome) %in% m2$`External ID`]
transcriptomics <- transcriptomics[, colnames(transcriptomics) %in% m2$`External ID`]
m3 <- m2 %>%
    select(data_type, site_sub_coll, `External ID`) %>%
    pivot_wider(names_from = data_type, values_from = `External ID`)
m2 <- select(m2, sex, week_num, site_sub_coll,
           visit_num, `Age at diagnosis`, `Age at diagnosis`, diagnosis,
           biopsy_location) %>%
    distinct()
m <- merge(m3, m2, all.x = TRUE, all.y = FALSE)
