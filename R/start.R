library("tidyverse")

hmp2 <- read.csv("data/hmp2_metadata.csv", check.names = FALSE,
                 stringsAsFactors = FALSE)

hmp2b <- hmp2 %>%
  filter(nzchar(biopsy_location),
         data_type %in% c("host_transcriptomics", "biopsy_16S"),
         !is.na(week_num),
         !startsWith(IntervalName, "Stool"),
         biopsy_location != "Non-inflamed")

# See particiipant and visit with data
hmp2b %>%
  ggplot() +
  geom_count(aes(visit_num, `Participant ID`)) +
  theme_minimal()

# Count samples with different time but samepatient
hmp2b %>%
  count(biopsy_location, `Participant ID`, data_type, week_num, sort = TRUE) %>%
  group_by(biopsy_location, `Participant ID`, week_num) %>%
  summarize(both = n_distinct(data_type)) %>%
  arrange(`Participant ID`, week_num, biopsy_location) %>%
  ungroup() %>%
  group_by(`Participant ID`) %>%
  summarize(n  = n_distinct(week_num)) %>%
  count(n)

hmp2b %>%
  ggplot() +
  geom_count(aes(biopsy_location, `Participant ID`, col = data_type)) +
  facet_wrap(~data_type) +
  theme_minimal()

hmp2b %>%
  group_by(`Participant ID`, visit_num) %>%
  mutate(both_types = all(data_type %in% c("biopsy_16S", "host_transcriptomics"))) %>%
  ungroup() %>%
  count(both_types)


# Remove uninformative columns ####
entropy <- function(x) {

  y <- prop.table(table(x, useNA = "ifany"))
  if (length(y) == 1) {
    return(0)
  }
  sum(-y*log(y, base = length(y)))
}
low_entro <- function(x) {
  x < 0.5
}

numeric_filter <- function(x) {
  if (sum(is.na(x)) == length(x)) {
    return(FALSE)
  }
  if (var(x, na.rm = TRUE) == 0) {
    return(FALSE)
  }

  if (sd(x, na.rm = TRUE)/abs(mean(x, na.rm = TRUE)) < 0.25) {
    return(FALSE)
  }

  return(TRUE)
}

mix <- function(x) {
  if (is.numeric(x)) {
    !low_entro(x)
  } else {
    isTRUE(x)
  }
}

meaningful_columns <- hmp2b %>%
  summarise(across(where(is.numeric), .fns = numeric_filter),
            across(where(is.character), .fns = entropy)) %>%
  select(where(mix)) %>%
  colnames()

paired <- select(hmp2b, meaningful_columns)

# Some samples from the same participant and visit
paired %>%
  group_by(`Participant ID`, visit_num) %>%
  summarize(n = n_distinct(biopsy_location)) %>%
  filter(n != 1)

# Some participants from same place different time
paired %>%
  group_by(`Participant ID`, biopsy_location) %>%
  summarize(n = n_distinct(visit_num)) %>%
  filter(n != 1)

p2 <- filter(paired, biopsy_location != "Non-inflamed") %>% count(`Participant ID`)

