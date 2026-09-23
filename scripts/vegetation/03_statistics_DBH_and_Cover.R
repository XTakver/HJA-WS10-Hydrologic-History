# ------------------------------------------------------------ #
# 03_statistics_DBH_and_Cover.R
#
# Purpose
#   Compare tree diameter-at-breast-height (DBH; overall and for
#   the dominant species P. menziesii) and total ground-cover area
#   between hillslopes, by survey year, for the WS10 long-term
#   vegetation plots.
#
# Inputs
#   data/vegetation/processed/WS10_TreeSurvey_AllPlots.csv
#   data/vegetation/processed/WS10_GroundCover_AllPlots.csv
#
# Outputs (written to SD)
#   - dbh_2016_summary.csv          / dbh_2016_ttest.csv
#   - dbh_PSME_2016_summary.csv     / dbh_PSME_2016_ttest.csv
#   - dbh_PSME_byYear_ttests.csv    (1989-present, Holm-adjusted)
#   - dbh_PSME_lm_summary.txt       / dbh_PSME_emmeans_pairwise_by_year.csv
#   - coverage_2016_summary.csv
#   - coverage_byYear_ttests.csv    (Holm-adjusted)
#
# Notes
#   - Welch's t-test (default for t.test()) is used throughout;
#     Holm and Benjamini-Hochberg-adjusted p-values are computed
#     where multiple years are compared.
#
# Created:  2025-08-26
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, glue, here, broom, emmeans)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

# Data Directory
DD <- here('data', 'vegetation', 'processed')

# Stats Directory
SD <- here('output', 'statistics', 'vegetation')

ensure_dir(SD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# All years, all plots. Trees
df_overstory <- read_csv(file.path(DD, 'WS10_TreeSurvey_AllPlots_processed.csv'))
# All years, all plots. Groundcover
df_grndcover.all <- read_csv(file.path(DD, 'WS10_GroundCover_AllPlots_processed.csv'))

## ----------------------------------------- ##
#               DBH Metrics ----
## ----------------------------------------- ##

dbh_data <- df_overstory |>
  filter(!is.na(dbh))

### All Species, 2016 ----

dbh_summary <- dbh_data |>
  filter(year == 2016) |>
  group_by(hillslope) |>
  summarize(n = n(),
            avg_dbh = mean(dbh, na.rm = TRUE),
            sd_dbh = sd(dbh, na.rm = TRUE),
            se_dbh = sd_dbh / sqrt(n),
            ci_95 = qt(0.975, df = n - 1) * se_dbh)          

dbh_2016 <- dbh_data |>
  filter(year == 2016)

t_res <- t.test(dbh ~ hillslope, data = dbh_2016)  # default = Welch t-test
t_res

t_res$estimate   # mean DBH in each hillslope
t_res$p.value    # p-value for difference in means
t_res$conf.int   # 95% CI for the difference in means

diff_means <- diff(t_res$estimate)  # order depends on factor levels
diff_means

write_csv(dbh_summary, file.path(SD, 'dbh_2016_summary.csv'))
write_csv(tidy(t_res), file.path(SD, 'dbh_2016_ttest.csv'))

### PSME Only, 2016 ----

dbh_PSME <- dbh_2016 |>
  filter(tree_spp == 'PSME')

dbh_PSME_summary <- dbh_PSME |>
  group_by(hillslope) |>
  summarize(n = n(),
            avg_dbh = mean(dbh, na.rm = TRUE),
            sd_dbh = sd(dbh, na.rm = TRUE),
            se_dbh = sd_dbh / sqrt(n),
            ci_95 = qt(0.975, df = n - 1) * se_dbh)

t_PSME <- t.test(dbh ~ hillslope, data = dbh_PSME)  # default = Welch t-test
t_PSME

write_csv(dbh_PSME_summary, file.path(SD, 'dbh_PSME_2016_summary.csv'))
write_csv(tidy(t_PSME), file.path(SD, 'dbh_PSME_2016_ttest.csv'))

### PSME, All Years (1989-present) ----

dbh_PSME_all <- dbh_data |>
  filter(tree_spp == 'PSME',
         year >= 1989)

# Check tree # for each year
dbh_PSME_all |>
  group_by(year, hillslope) |>
  summarize(n = n(), .groups = 'drop') |>
  print(n = Inf)

per_year_tests <- dbh_PSME_all |>
  group_by(year) |>
  do(tidy(t.test(dbh ~ hillslope, data = .))) |>
  ungroup() |>
  mutate(p_adj = p.adjust(p.value, method = 'holm'))

per_year_tests

write_csv(per_year_tests, file.path(SD, 'dbh_PSME_byYear_ttests.csv'))

# Linear model (alternative to per-year t-tests)
mod <- lm(dbh ~ hillslope * factor(year), data = dbh_PSME_all)

emm <- emmeans(mod, ~ hillslope | year)
pairs(emm, adjust = 'holm')

writeLines(capture.output(summary(mod)), file.path(SD, 'dbh_PSME_lm_summary.txt'))
write_csv(as.data.frame(pairs(emm, adjust = 'holm')),
          file.path(SD, 'dbh_PSME_emmeans_pairwise_by_year.csv'))

## ----------------------------------------- ##
#         Coverage area estimates ----
## ----------------------------------------- ##

total_CA <- df_grndcover.all |>
  group_by(hillslope, year, plot) |>
  summarize(CA = sum(cover)) |>
  ungroup()

# summary of coverage area
CA_2016 <- total_CA |>
  filter(year == 2016) |>
  group_by(hillslope) |>
  summarize(n = n(),
            avg_CA = mean(CA, na.rm = TRUE),
            sd_CA = sd(CA, na.rm = TRUE),
            se_CA = sd_CA / sqrt(n),
            ci_95 = qt(0.975, df = n - 1) * se_CA)

write_csv(CA_2016, file.path(SD, 'coverage_2016_summary.csv'))

# Check CA for each year
total_CA |>
  group_by(year, hillslope) |>
  summarize(n = n(), .groups = 'drop') |>
  print(n = Inf)

per_year_tests <- total_CA |>
  group_by(year) |>
  do(tidy(t.test(CA ~ hillslope, data = .))) |>
  ungroup() |>
  mutate(p_adj = p.adjust(p.value, method = 'holm'))

per_year_tests

write_csv(per_year_tests, file.path(SD, 'coverage_byYear_ttests.csv'))
