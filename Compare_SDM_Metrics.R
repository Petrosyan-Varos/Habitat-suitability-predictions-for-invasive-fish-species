###############################################################
# Script name: Compare_SDM_Metrics.R
# Purpose: Statistical comparison of SDM model performance
# Author: Varos Petrosyan, Fedor Osipov
# Project: Multi-scale SDM analysis (Armenia)
# R version: 3.6.2+
#
# Description:
# This script performs:
# - Standardization of SDM evaluation metrics (AUC, TSS, CBI)
# - Construction of combined performance indices
# - Model ranking
# - Boxplots and density plots
# - ANOVA, Kruskal–Wallis, pairwise Wilcoxon tests (BH correction)
# - Correlation analysis and heatmap visualization
#
# Input:
# CSV file with columns:
# model, fold, AUC, TSS, CBI
#
# Output:
# - Summary CSV tables
# - PNG plots
# - TXT statistical reports
####  License
# This code is released under the MIT License.
# You are free to use, modify, and redistribute with attribution.
###############################################################

###############################################################

# ============================================================
# 1. Libraries
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(reshape2)
library(multcompView)

# ============================================================
# 2. Paths and input files
# ============================================================

driveName      <- "working_directory" #
MainDirectory  <- paste0(driveName, "Armenia_2025")
PathToSources  <- "Sources_directory" # 
PathToResults  <- "Results_directory" # 

# Choose input file from Compare_SDM.zip
## optional 
Compare_SDM <- paste0(PathToSources, "/", "Parva_Compare_SDM_V2.csv")
## optional 
 Compare_SDM <- paste0(PathToSources, "/", "Caras_Compare_SDM_V2.csv")

# Load data
df <- read_csv(Compare_SDM)

# Basic formatting
df$model <- as.factor(df$model)
df$fold  <- as.factor(df$fold)

# Safety check
stopifnot(all(c("model","fold","AUC","TSS","CBI") %in% names(df)))

# ============================================================
# 3. Metric standardization
# ============================================================

df <- df %>%
  mutate(
    AUC   = pmin(pmax(AUC, 0), 1),
    TSS   = pmax(TSS, -1),
    TSS01 = (TSS + 1) / 2,
    CBI01 = (CBI + 1) / 2
  )

# ============================================================
# 4. Metric grading
# ============================================================

auc_grade <- function(x){
  if(x >= 0.9) return("excellent")
  if(x >= 0.8) return("good")
  if(x >= 0.7) return("satisfactory")
  return("poor")
}

tss_grade <- function(x){
  if(x >= 0.8) return("excellent")
  if(x >= 0.6) return("good")
  if(x >= 0.4) return("satisfactory")
  return("poor")
}

cbi_grade <- function(x){
  v <- (x + 1)/2
  if(v >= 0.8) return("excellent")
  if(v >= 0.6) return("good")
  if(v >= 0.4) return("satisfactory")
  return("poor")
}

grade_map <- c("poor"=2, "satisfactory"=3, "good"=4, "excellent"=5)

df <- df %>%
  mutate(
    AUC_grade = sapply(AUC, auc_grade),
    TSS_grade = sapply(TSS, tss_grade),
    CBI_grade = sapply(CBI, cbi_grade),
    AUC_score = grade_map[AUC_grade],
    TSS_score = grade_map[TSS_grade],
    CBI_score = grade_map[CBI_grade]
  )

# ============================================================
# 5. Combined performance indices
# ============================================================

# (1) Equal weights
df$COMB_equal <- (df$AUC + df$TSS01 + df$CBI01) / 3

# (2) Expert-defined weights
w <- c(AUC = 0.3, TSS01 = 0.3, CBI01 = 0.4)
df$COMB_expert <- df$AUC*w["AUC"] + df$TSS01*w["TSS01"] + df$CBI01*w["CBI01"]

# (3) Inverse-variance weighting
var_tbl <- df %>%
  group_by(model) %>%
  summarise(
    var_auc = var(AUC, na.rm = TRUE),
    var_tss = var(TSS01, na.rm = TRUE),
    var_cbi = var(CBI01, na.rm = TRUE)
  )

eps <- 1e-6

var_tbl <- var_tbl %>%
  mutate(
    w_auc = 1/(var_auc + eps),
    w_tss = 1/(var_tss + eps),
    w_cbi = 1/(var_cbi + eps),
    wsum  = w_auc + w_tss + w_cbi,
    w_auc = w_auc/wsum,
    w_tss = w_tss/wsum,
    w_cbi = w_cbi/wsum
  )

df <- df %>% left_join(var_tbl, by="model")

df$COMB_invvar <- df$AUC*df$w_auc + df$TSS01*df$w_tss + df$CBI01*df$w_cbi

# ============================================================
# 6. Summary table
# ============================================================

summary_tbl <- df %>%
  group_by(model) %>%
  summarise(
    AUC_mean  = mean(AUC),
    TSS_mean  = mean(TSS),
    CBI_mean  = mean(CBI),
    COMB_equal_mean   = mean(COMB_equal),
    COMB_expert_mean  = mean(COMB_expert),
    COMB_invvar_mean  = mean(COMB_invvar)
  ) %>%
  arrange(desc(COMB_equal_mean))

write_csv(summary_tbl,
          paste0(PathToResults, "/", "SDM_Comparison_Summary.csv"))

# ============================================================
# 7. Visualization
# ============================================================

# ---- Boxplots ----
df_melt <- df %>%
  select(model, AUC, TSS, CBI) %>%
  reshape2::melt(id.vars="model")

p1 <- ggplot(df_melt, aes(x=model, y=value, fill=variable)) +
  geom_boxplot() +
  theme_bw() +
  labs(title="Distribution of AUC, TSS and CBI")

ggsave(paste0(PathToResults, "/", "Metrics_Boxplot.png"),
       p1, width=9, height=6, dpi=300)

# ---- Density plots ----
p2 <- ggplot(df_melt, aes(x=value, color=variable)) +
  geom_density(size=1) +
  theme_bw() +
  labs(title="Density distribution of SDM metrics")

ggsave(paste0(PathToResults, "/", "Metrics_Density.png"),
       p2, width=9, height=6, dpi=300)

# ---- Ranking plots ----
p3 <- ggplot(summary_tbl,
             aes(x=reorder(model, COMB_equal_mean),
                 y=COMB_equal_mean)) +
  geom_point(size=4) +
  geom_segment(aes(xend=model, y=0, yend=COMB_equal_mean)) +
  coord_flip() +
  theme_bw() +
  labs(title="Model ranking (COMB_equal)",
       x="Model", y="Score")

ggsave(paste0(PathToResults, "/", "Ranking_COMB_equal.png"),
       p3, width=7, height=5, dpi=300)

p4 <- ggplot(summary_tbl,
             aes(x=reorder(model, COMB_invvar_mean),
                 y=COMB_invvar_mean)) +
  geom_point(size=4) +
  geom_segment(aes(xend=model, y=0, yend=COMB_invvar_mean)) +
  coord_flip() +
  theme_bw() +
  labs(title="Model ranking (COMB_invvar)",
       x="Model", y="Score")

ggsave(paste0(PathToResults, "/", "Ranking_COMB_invvar_mean.png"),
       p4, width=7, height=5, dpi=300)

# ============================================================
# 8. Statistical tests for individual metrics
# ============================================================

# ---- ANOVA ----
anova_auc <- aov(AUC ~ model, data=df)
anova_tss <- aov(TSS ~ model, data=df)
anova_cbi <- aov(CBI ~ model, data=df)

sink(paste0(PathToResults, "/ANOVA_results.txt"))
cat("=== ANOVA AUC ===\n"); print(summary(anova_auc))
cat("\n=== ANOVA TSS ===\n"); print(summary(anova_tss))
cat("\n=== ANOVA CBI ===\n"); print(summary(anova_cbi))
sink()

# ---- Kruskal–Wallis ----
kw_auc <- kruskal.test(AUC ~ model, data=df)
kw_tss <- kruskal.test(TSS ~ model, data=df)
kw_cbi <- kruskal.test(CBI ~ model, data=df)

sink(paste0(PathToResults, "/KW_results.txt"))
print(kw_auc); print(kw_tss); print(kw_cbi)
sink()

# ---- Pairwise Wilcoxon (BH correction) ----
wilc_auc <- pairwise.wilcox.test(df$AUC, df$model, p.adjust.method="BH")
wilc_tss <- pairwise.wilcox.test(df$TSS, df$model, p.adjust.method="BH")
wilc_cbi <- pairwise.wilcox.test(df$CBI, df$model, p.adjust.method="BH")

sink(paste0(PathToResults, "/Wilcoxon_results.txt"))
print(wilc_auc); print(wilc_tss); print(wilc_cbi)
sink()

# ============================================================
# 9. Correlation analysis
# ============================================================

cor_tbl <- df %>%
  select(AUC, TSS, CBI, COMB_equal, COMB_expert, COMB_invvar)

cor_matrix <- cor(cor_tbl, use="pairwise.complete.obs")

write_csv(as.data.frame(cor_matrix),
          paste0(PathToResults, "/", "Correlation_Table.csv"))

# ---- Correlation heatmap ----
cor_melt <- reshape2::melt(cor_matrix)
colnames(cor_melt) <- c("Metric1","Metric2","Correlation")

p <- ggplot(cor_melt, aes(x=Metric1, y=Metric2, fill=Correlation)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0.5, limit=c(0,1), name="r") +
  geom_text(aes(label=round(Correlation, 2)), size=4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Correlation heatmap of SDM metrics")

ggsave(paste0(PathToResults, "/", "Correlation_Heatmap.png"),
       p, width=7, height=5, dpi=300)

# ============================================================
# 10. Statistics for combined indices
# ============================================================

anova_comb_equal  <- aov(COMB_equal ~ model, data=df)
anova_comb_expert <- aov(COMB_expert ~ model, data=df)
anova_comb_invvar <- aov(COMB_invvar ~ model, data=df)

sink(paste0(PathToResults, "/ANOVA_COMB_results.txt"))
cat("=== ANOVA COMB_equal ===\n");  print(summary(anova_comb_equal))
cat("\n=== ANOVA COMB_expert ===\n"); print(summary(anova_comb_expert))
cat("\n=== ANOVA COMB_invvar ===\n"); print(summary(anova_comb_invvar))
sink()

kw_comb_equal  <- kruskal.test(COMB_equal ~ model, data=df)
kw_comb_expert <- kruskal.test(COMB_expert ~ model, data=df)
kw_comb_invvar <- kruskal.test(COMB_invvar ~ model, data=df)

sink(paste0(PathToResults, "/KW_COMB_results.txt"))
print(kw_comb_equal); print(kw_comb_expert); print(kw_comb_invvar)
sink()

wilc_comb_equal  <- pairwise.wilcox.test(df$COMB_equal,  df$model, p.adjust.method="BH")
wilc_comb_expert <- pairwise.wilcox.test(df$COMB_expert, df$model, p.adjust.method="BH")
wilc_comb_invvar <- pairwise.wilcox.test(df$COMB_invvar, df$model, p.adjust.method="BH")

sink(paste0(PathToResults, "/Wilcoxon_COMB_results.txt"))
print(wilc_comb_equal); print(wilc_comb_expert); print(wilc_comb_invvar)
sink()


#######################################################################

get_letters <- function(x, g){
  pw <- pairwise.wilcox.test(x, g, p.adjust.method = "BH")
  letters <- multcompView::multcompLetters(pw$p.value)$Letters
  return(data.frame(model = names(letters), letters = letters))
}


# Convert data to long format (COMB ONLY)

df_comp_long <- df %>%
  select(model, COMB_equal, COMB_expert, COMB_invvar) %>%
  pivot_longer(
    cols = c(COMB_equal, COMB_expert, COMB_invvar),
    names_to = "COMB_type",
    values_to = "COMB_value"
  )

df_comp_long$model <- factor(df_comp_long$model,
                             levels = c("G-SDM", "R-SDM", "ON-SDM", "MN-SDM")
)

df_comp_long$COMB_type <- factor(df_comp_long$COMB_type,
                                 levels = c("COMB_equal", "COMB_expert", "COMB_invvar")
)

## We obtain SIGNIFICANCE LETTERS for each index
letters_tbl <- df_comp_long %>%
  group_by(COMB_type) %>%
  do(get_letters(.$COMB_value, .$model))

## Maximums for placing letters above a boxplot

ymax_tbl <- df_comp_long %>%
  group_by(COMB_type, model) %>%
  summarise(y = max(COMB_value)) %>%
  left_join(letters_tbl, by = c("COMB_type", "model"))

## FIGURE X ? 

model_colors <- c(
  "G-SDM"  = "#F8766D",
  "R-SDM"  = "#7CAE00",
  "ON-SDM" = "#00BFC4",
  "MN-SDM" = "#C77CFF"
)
### ?????? 11
p_Figure_X <- ggplot(df_comp_long, aes(x = model, y = COMB_value)) +
  geom_boxplot(outlier.size = 1.3) +
  facet_wrap(~ COMB_type, scales = "free_y") +
  geom_text(
    data = ymax_tbl,
    aes(x = model, y = y * 1.06, label = letters),
    size = 5,
    fontface = "bold"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Figure X. Composite indices across SDM frameworks",
    x = "SDM framework",
    y = "Composite index value"
  )
########################################################################
p_Figure_X <- ggplot(df_comp_long, aes(x = model, y = COMB_value, fill = model)) +
  geom_boxplot(outlier.size = 1.3) +
  scale_fill_manual(values = model_colors) +
  facet_wrap(~ COMB_type, scales = "free_y") +
  geom_text(
    data = ymax_tbl,
    aes(x = model, y = y * 1.06, label = letters),
    size = 5,
    fontface = "bold"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Figure X. Composite indices across SDM frameworks",
    x = "SDM framework",
    y = "Composite index value"
  )

plot(p_Figure_X)

ggsave(
  paste0(PathToResults, "/Figure_X_Composite_Indices_colored.png"),
  p_Figure_X,
  width = 10, height = 5, dpi = 300
)

# ============================================================
# End of script
# ============================================================
