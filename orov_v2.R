##----- Dependencies -----------------------------------------------------------

# Set system language
Sys.setenv(LANG = "en")

# Packages
library(tidyverse)
library(xlsx)
library(ggpubr)
library(ggh4x)
library(scales)
library(ggbeeswarm)
library(rstatix)
library(CompBioTools)
library(reshape2)
library(patchwork)
library(STRINGdb)
library(RITAN)
library(ggraph)
library(igraph)
library(ggbreak)
library(scales)
library(ppcor)
library(boot)
library(future.apply)

#@ Improving plotting in plot panel
trace(grDevices:::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

trace(grDevices::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)


## Functions

# Shrinking negative values to make HL effect sizes plot better
neg_shrink_trans <- function(factor = 5) {
  trans_new(
    name = paste0("negshrink-", factor),
    transform = function(x) ifelse(x < 0, x / factor, x),
    inverse   = function(x) ifelse(x < 0, x * factor, x)
  )
}

# Boostrapping correlation difference matrix
bootstrap_corr_diff <- function(data_A, data_B, n_boot = 1000) {
  
  # Remove columns with all NAs
  data_A <- data_A[, colSums(!is.na(data_A)) > 0, drop=FALSE]
  data_B <- data_B[, colSums(!is.na(data_B)) > 0, drop=FALSE]
  
  p <- ncol(data_A)
  vars <- colnames(data_A)
  
  diff_boot <- array(0, dim = c(p, p, n_boot))
  
  # Bootstrapping differences
  for (b in 1:n_boot) {
    
    idxA <- sample(1:nrow(data_A), replace = TRUE)
    idxB <- sample(1:nrow(data_B), replace = TRUE)
    
    Ca_b <- suppressWarnings(
      cor(data_A[idxA,], method = "spearman",
          use = "pairwise.complete.obs")
    )
    
    Cb_b <- suppressWarnings(
      cor(data_B[idxB,], method = "spearman",
          use = "pairwise.complete.obs")
    )
    
    diff_boot[,,b] <- Cb_b - Ca_b
  }
  
  # Observed differences
  Ca <- cor(data_A, method = "spearman",
            use = "pairwise.complete.obs")
  Cb <- cor(data_B, method = "spearman",
            use = "pairwise.complete.obs")
  diff_obs <- Cb - Ca
  
  # CI
  CI_low  <- apply(diff_boot, c(1,2), quantile, probs = 0.025)
  CI_high <- apply(diff_boot, c(1,2), quantile, probs = 0.975)
  
  # p-values
  p_mat <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      x <- diff_boot[i,j,]
      obs <- diff_obs[i,j]
      p_mat[i,j] <- mean(sign(x) != sign(obs))
    }
  }
  
  # FDR
  p_fdr <- matrix(
    p.adjust(as.vector(p_mat), method="BH"),
    nrow = p, ncol = p
  )
  
  colnames(diff_obs) <- rownames(diff_obs) <- vars
  colnames(CI_low)   <- rownames(CI_low)   <- vars
  colnames(CI_high)  <- rownames(CI_high)  <- vars
  colnames(p_mat)    <- rownames(p_mat)    <- vars
  colnames(p_fdr)    <- rownames(p_fdr)    <- vars
  
  list(
    diff   = diff_obs,
    CI_low = CI_low,
    CI_high= CI_high,
    p      = p_mat,
    p_fdr  = p_fdr,
    boot   = diff_boot
  )
}

# Partial Spearman function to use for all columns
partial_spearman_fun <- function(data, indices, cols) {
  
  # Resampled data
  d <- data[indices, cols]  
  d <- d[complete.cases(d),]
  #if (nrow(d) < length(cols) + 1) return(NA)
  
  # Run partical cor and get value
  r <- pcor(d, method = "spearman")$estimate[2]
  
  return(r)
}

# Partial Spearman bootstrap function
boot_partial_spearman <- function(data, perm, cols_vars, cols_mult) {
  
  # Start results df
  results <- data.frame()
  
  # For each cytokine column, run partial spearman
  for (col in cols_mult) {
    
    # Setting up columns
    cols <- c(col, cols_vars)
    
    # Bootstrapping
    boot_res <- boot(data = data,
                     statistic = partial_spearman_fun,
                     R = perm, 
                     cols = cols)
    
    # Getting res
    p_val <- mean(boot_res$t * boot_res$t0 <= 0)
    ci_low <- boot.ci(boot_res, type = "bca")$bca[4]
    ci_high <- boot.ci(boot_res, type = "bca")$bca[5]
    results <- rbind(results,
                     data.frame(variable = col,
                                correlation = boot_res$t0,
                                ci_low = ci_low,
                                ci_high = ci_high,
                                p_val = p_val)
    )
    
  }
  
  return(results)
  
}

## -----------------------------------------------------------------------------

## Initial processing

# Reading data
df_orov <- readxl::read_xlsx("OROV_Cytokines.xlsx", sheet = 2)

# Making it long
df_orov_melt <- melt(df_orov,
                     id.vars = c("Sample", "Days", "CT", "Age", "Sex", "Group"),
                     value.name = "Concentration",
                     variable.name = "Cytokine")

## ----- Outliers --------------------------------------------------------------

# Identify the target columns (6 to 18)
target_cols <- names(df_orov)[6:18]

# Create summary table
df_outliers <- df_orov %>%
  mutate(across(all_of(target_cols), ~is.na(.x))) %>%
  pivot_longer(
    cols = all_of(target_cols),
    names_to = "Column",
    values_to = "IsNA"
  ) %>%
  group_by(Sample, Group) %>%
  summarize(
    Outliers = paste(Column[IsNA], collapse = ", "),
    `Outlier Count` = sum(IsNA),
    `Proportion` = sum(IsNA)/length(target_cols),
    .groups = "drop"
  ) %>%
  mutate(Outliers = ifelse(Outliers == "", NA, Outliers))  %>%
  as.data.frame()

df_outliers_anon <- df_outliers
df_outliers_anon$Sample <- paste0("Sample", seq_len(nrow(df_outliers_anon)))
  
#  Marking outliers in original table
df_orov_melt$is_outlier <- NA

for (s in 1:length(df_outliers$Sample)) {
  
  sample <- df_outliers$Sample[s]
  outs <- str_split(df_outliers[s, "Outliers"], ", ")[[1]]

  for (o in outs) {
    
    if (!is.na(o)) {
      df_orov_melt[df_orov_melt$Sample==sample&df_orov_melt$Cytokine==o,"is_outlier"] <- TRUE
    }
    
  }
}

for (i in 1:length(df_orov_melt$is_outlier)) {
  if (is.na(df_orov_melt$is_outlier[i])) {
    df_orov_melt$is_outlier[i] <- FALSE
  }
}

# Checking if counts match
df_orov_melt$is_outlier %>% table()
sum(df_outliers$`Outlier Count`)

## ----- Statistics ------------------------------------------------------------

df_orov_melt$Group <- factor(df_orov_melt$Group, levels = c("HC", "OROV"))

# Wilcoxon with FDR (BH) and Effect Sizes
wilcox_results <- df_orov_melt[df_orov_melt$is_outlier == FALSE, ] %>%
  group_by(Cytokine) %>%
  summarise(
    # Store test object as a list-column
    wt = list(
      tryCatch(
        wilcox.test(Concentration ~ Group, exact = FALSE, conf.int = TRUE),
        error = function(e) NULL
      )
    ),
    
    # p-value
    p_value = if (!is.null(wt[[1]])) wt[[1]]$p.value else NA_real_,
    
    # Hodges–Lehmann median difference and CI
    hl_estimate = if (!is.null(wt[[1]])) wt[[1]]$estimate else NA_real_,
    hl_ci_low  = if (!is.null(wt[[1]])) wt[[1]]$conf.int[1] else NA_real_,
    hl_ci_high = if (!is.null(wt[[1]])) wt[[1]]$conf.int[2] else NA_real_,
    
    # Sample sizes
    n_group1 = sum(Group == unique(Group)[1]),
    n_group2 = sum(Group == unique(Group)[2])) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

# Coding stars
wilcox_results <- pvalueCodes(wilcox_results, "p_adj")

### Saving results
write.xlsx(wilcox_results[,-2], "results/results.xlsx", sheetName = "Wilcoxon", showNA = FALSE)
write.xlsx(df_outliers_anon, "results/results.xlsx", sheetName = "Outliers", showNA = FALSE, append = TRUE)

### Plotting and saving effect sizes
tiff("results/effect_sizes.tiff", res = 300, w = 2500, h = 1500, units = "px")

ggplot(
  wilcox_results,
  aes(x = Cytokine, y = hl_estimate,
      ymin = hl_ci_low, ymax = hl_ci_high)
) +
  geom_pointrange(size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_y_continuous(trans = neg_shrink_trans(factor = 15),
                     breaks = c(-3000,-2000,-1000,0,35,70)) +
  coord_flip() +
  theme_bw() +
  labs(
    y = "Hodges–Lehmann Median Difference",
    x = "Cytokine"
  ) +
  theme_bw() +
  geom_point(aes(x = Inf, y = Inf, colour = "HL Estimate"), size = 2.75, alpha = 0) +
  geom_line(aes(x = Inf, y = Inf, colour = "95% CI"), alpha = 0) +
  scale_colour_manual(values = c("black", "black"),
                      labels = c("95% CI", "HL Estimate"),
                      name = "Legend") +
  guides(colour = guide_legend(reverse = TRUE, override.aes = list(alpha = 1))) +
  theme(text = element_text(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"), 
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold")) +
  annotate("text",
           x = "CCL2 (MCP-1)",
           y = -1913.2599663,
           label = "*",
           size = 6,
           vjust = 0) +
  annotate("text",
           x = "IL-18",
           y = -344.6399680,
           label = "*",
           size = 6,
           vjust = 0) +
  annotate("text",
           x = "TNF-α",
           y = 33.1857605,
           label = "*",
           size = 6,
           vjust = 0)

dev.off()

# Preparing significance for boxplot
stars_df <- df_orov_melt %>%
  filter(is_outlier == FALSE) %>%
  group_by(Cytokine) %>%
  summarise(ymax = max(Concentration, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(y = ymax * 1.05,
         y_bracket = ymax * 1.02,
         label = "*",
         Group_left  = "HC",
         Group_right = "OROV") %>%
  as.data.frame()

# Reordering
wilcox_results$Cytokine <- as.character(wilcox_results$Cytokine)
wilcox_results <- wilcox_results[order(wilcox_results$Cytokine),]
stars_df$Cytokine <- as.character(stars_df$Cytokine)
stars_df <- stars_df[order(stars_df$Cytokine),]
identical(stars_df$Cytokine, wilcox_results$Cytokine)

# Replacing
stars_df$label <- wilcox_results$pcode
stars_df <- stars_df[stars_df$label!="",]

# Setting up order
cytok_order <- unique(df_orov_melt$Cytokine)[order(unique(as.character(df_orov_melt$Cytokine)))]
df_orov_melt$Cytokine <- factor(df_orov_melt$Cytokine, levels = cytok_order)

### Saving boxplot
tiff("results/boxplots.tiff", res = 300, w = 2500, h = 2500, units = "px")

ggplot(df_orov_melt[df_orov_melt$is_outlier==FALSE,], aes(x = Group, y = Concentration)) +
  facet_wrap(~ Cytokine, scales = "free") +
  geom_quasirandom() +
  geom_boxplot(outliers = FALSE, staplewidth = 0.3, width = 0.5, alpha = 0.75) +
  geom_text(
    data = stars_df,
    aes(x = 1.5, y = y, label = label),
    size = 6,
    vjust = 0
    ) +
  geom_segment(
    data = stars_df,
    aes(x = Group_left, xend = Group_right,
        y = y_bracket*1.1, yend = y_bracket*1.1),
    linewidth = 0.4
  ) +
  scale_y_continuous(expand = c(0,0,0.15,0)) +
  labs(x = "", y = "Concentration (pg/mL)") +
  theme_bw() +
  theme(text = element_text(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"), 
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"))

dev.off()

## ----- Fold change analysis --------------------------------------------------

# Medians
median_orov <- df_orov[,5:18] %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), \(x) median(x, na.rm = TRUE))) %>%
  t() %>%
  as.data.frame() 

# Replacing colnames
colnames(median_orov) <- c("HC", "OROV")

# Removing extra line and setting up rownames/datatype
median_orov <- median_orov[-c(1),]
median_orov$Cytokine <- rownames(median_orov)
median_orov[,1:2] <- apply(median_orov[,1:2], 2, as.numeric)

# Calc fold change
median_orov$FoldChange <- median_orov$OROV/median_orov$HC

# Checking order to replace with pcodes
median_orov <- median_orov[wilcox_results$Cytokine,]
identical(as.character(wilcox_results$Cytokine), median_orov$Cytokine)
median_orov$stats <- wilcox_results$pcode

### Ploting and saving
tiff("results/fold_change.tiff", res = 300, w = 2000, h = 1250, units = "px")

ggplot(median_orov, aes(x = Cytokine, y = log2(FoldChange))) +
  geom_bar(stat = "identity", fill = "gray70", color = "black", width = 0.85) +
  geom_text(aes(label = stats), vjust = ifelse(log2(median_orov$FoldChange) < 0, 1.5, 0), size = 4) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0.5)) +
  theme(axis.text.x = element_text(size = 12, angle = -45, hjust = 0), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth = 0.5),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(10, 10, 10, 20)) +
  labs(y = "Fold Change (log2)") +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.25) +
  coord_cartesian(clip = "off")

dev.off()

## ----- Correlation matrices --------------------------------------------------

## Processing data

# Separating patients and controls
df_orov_hc <- df_orov[df_orov$Group=="HC",]
df_orov_pt <- df_orov[df_orov$Group=="OROV",]

# Calculating diff matrix with bootstrap and FDR
bt_corr_diff <- bootstrap_corr_diff(df_orov_hc[,6:18], df_orov_pt[,6:18])

# Preparing for plot
heat_df <- as.data.frame(bt_corr_diff$diff) %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "diff") %>%
  mutate(
    p   = as.vector(bt_corr_diff$p),
    FDR = as.vector(bt_corr_diff$p_fdr),
    CI_high = as.vector(bt_corr_diff$CI_high),
    CI_low = as.vector(bt_corr_diff$CI_low)
  )

# Adjusting levels
levels_order <- unique(heat_df$var1)[order(unique(heat_df$var1))]
heat_df$var1 <- factor(heat_df$var1, levels = levels_order)
heat_df$var2 <- factor(heat_df$var2, levels = levels_order)

# Coding p-values
heat_df <- pvalueCodes(heat_df, "FDR")

### Saving plot
tiff("results/diff_network.tiff", res = 300, w = 2000, h = 1500, units = "px")

# Get significant results onlye
sig <- ((heat_df$CI_low > 0 | heat_df$CI_high < 0) & heat_df$FDR <= 0.05)

ggplot(heat_df, aes(x = var1, y = var2, fill = diff)) +
  geom_tile(color = "white") +
  geom_rect(data = heat_df %>% filter(sig), 
            aes(xmin = as.numeric(var1) - 0.5,
                xmax = as.numeric(var1) + 0.5,
                ymin = as.numeric(var2) - 0.5,
                ymax = as.numeric(var2) + 0.5),
            fill = NA,
            color = "black",
            linewidth = 0.25) +
  geom_text(aes(label = pcode), color = "black", size = 3, lineheight = 0.75) +
  scale_fill_gradient2(low = "#2166AC", high = "#B2182B", mid = "#F7F7F7", 
                       midpoint = 0, limit = c(-2,2), space = "Lab", 
                       name="Δρ") +
  guides(colour = "none") +
  scale_colour_manual(values = c("no" = NA, "yes" = "black")) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.ticks = element_blank(),
        legend.justification = c(0.5, 0), 
        legend.position = "right",
        legend.direction = "vertical", 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  guides(fill = guide_colorbar(barwidth = 1, 
                               barheight = 5,
                               title.position = "top", 
                               title.hjust = 0.5)) +
  coord_fixed() +
  labs(x = "", y = "")

dev.off()


## ----- STRING analysis -------------------------------------------------------

# Cytokine names
cytokine_genes <- read.xlsx("citocinas_genes.xlsx", 1)

# Reading file from string website
net_cytokines <- read_tsv("string_interactions_short.tsv")[,c(1,2,9,10,11,13)]

# Replacing names
net_cytokines$`#node1` <- with(cytokine_genes, name[match(net_cytokines$`#node1`, hgnc_symbol)])
net_cytokines$`node2` <- with(cytokine_genes, name[match(net_cytokines$`node2`, hgnc_symbol)])

# Constructing edges
edges <- data.frame(
  from = net_cytokines$`#node1`,
  to   = net_cytokines$node2,
  Score = net_cytokines$combined_score
)

# Defining breaks in edges
edges$Score_cat <- cut(
  edges$Score,
  breaks = c(-Inf, 0.4, 0.7, 0.9, Inf),
  labels = c("Low", "Medium", "High", "Highest"),
  right = FALSE
)

# Setting up nodes
nodes <- data.frame(
  name = median_orov$Cytokine,
  FoldChange = log2(median_orov$FoldChange)
)

# Creating graph
graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# Create the layout
layout <- create_layout(graph, layout = 'circle')

# Calculate nudge_x and nudge_y to move the labels outside the circle
layout <- layout %>%
  mutate(label_length = nchar(name), 
         nudge_x = x * 0.075,
         nudge_y = y * 0.075)
layout <- layout %>%
  mutate(vjust = ifelse(y > 0, 0, 1),
         hjust = ifelse(x > 0, 0, 1))

# Which cytokines (sig) to colour
which_to_colour <- c("TNF-α", "IL-18", "CCL2 (MCP-1)")

# Plotting
tiff("results/string_network.tiff", res = 300, w = 2000, h = 2000, units = "px")

ggraph(layout) + 
  geom_edge_link(aes(color = Score), linewidth = 1) + 
  geom_node_point(aes(fill = ifelse(name %in% which_to_colour, FoldChange, NA)), 
                  color = "black", size = 10, shape = 21, stroke = 1) + 
  geom_node_text(aes(label = name), 
                 nudge_x = layout$nudge_x,
                 nudge_y = layout$nudge_y,
                 vjust = layout$vjust,
                 hjust = layout$hjust,
                 size = 4, 
                 fontface = "bold") + 
  scale_fill_distiller(palette = "Spectral", 
                       direction = -1, 
                       name = "Fold Change\n(log2)",
                       limits = c(-5, 5)) + 
  theme_graph(background = "white") + 
  scale_edge_color_stepsn(breaks = c(0.4, 0.7, 0.9),
                          colors = c("#ffff45", "#ff2345"),
                          name = "Confidence\nScore (STRING)") + 
  guides(fill = guide_colorbar(barwidth = 6, 
                               barheight = 0.75,
                               title.position = "top", 
                               title.hjust = 0.5),
         edge_colour = guide_coloursteps(barwidth = 6, 
                               barheight = 0.75,
                               title.position = "top", 
                               title.hjust = 0.5)) +
  theme(legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.text = element_text(hjust = 0.5),
        legend.key.width = unit(0.75, "cm"),
        legend.key.height = unit(0.6, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.title.position = "top",
        legend.box.just = "center") +
  coord_fixed(clip = "off")

dev.off()

## ----- CT analysis -----------------------------------------------------------

# Adjusting group 
df_orov_pt$Group2 <- ifelse(df_orov_pt$Days <= 2, "1-2 days", "3-4 days")

# Wilcoxon with FDR (BH) and Effect Sizes
wilcox_results_pt <- df_orov_pt %>%
  summarise(
    
    # Store test object as a list-column
    wt = list(
      tryCatch(wilcox.test(CT ~ Group2, exact = TRUE, conf.int = TRUE),
               error = function(e) NULL)
      ),
    
    # p-value
    p_value = if (!is.null(wt[[1]])) wt[[1]]$p.value else NA_real_,
    
    # Hodges–Lehmann median difference and CI
    hl_estimate = if (!is.null(wt[[1]])) wt[[1]]$estimate else NA_real_,
    hl_ci_low  = if (!is.null(wt[[1]])) wt[[1]]$conf.int[1] else NA_real_,
    hl_ci_high = if (!is.null(wt[[1]])) wt[[1]]$conf.int[2] else NA_real_,
    
    # Sample sizes
    n_group1 = sum(Group2 == unique(Group2)[1]),
    n_group2 = sum(Group2 == unique(Group2)[2])) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Coding p
wilcox_results_pt <- pvalueCodes(wilcox_results_pt, "p_adj")

# Preparing significance for boxplot
stars_df_ct <- df_orov_pt %>%
  summarise(ymax = max(CT, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(y = ymax * 1.05,
    y_bracket = ymax * 1.02,
    label = "*",
    Group_left  = "1-2 days",
    Group_right = "3-4 days") %>%
  as.data.frame()

### Plotting and saving
tiff("results/ct_days.tiff", res = 300, w = 1000, h = 1000, units = "px")

ggplot(df_orov_pt, aes(x = Group2, y = -CT)) +
  geom_quasirandom() +
  geom_boxplot(outliers = FALSE, staplewidth = 0.3, width = 0.5, alpha = 0.75) +
  geom_text(data = stars_df_ct,
            aes(x = 1.5, y = -15, label = label),
            size = 6,
            vjust = 0) +
  geom_segment(data = stars_df_ct,
               aes(x = Group_left, xend = Group_right,
                   y = -15, yend = -15),
               linewidth = 0.4) +
  scale_y_continuous(expand = c(0,0.25,0,0.5)) +
  theme_bw() +
  theme(text = element_text(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold")) +
  labs(x = "Days of symptoms")

dev.off()

## ----- Cytokines by days -----------------------------------------------------

# Coding days
df_orov_melt$Group2 <- ifelse(is.na(df_orov_melt$Days), "HC", ifelse(df_orov_melt$Days <= 2, "1-2 days", "3-4 days"))
df_orov_melt$Group2 <- factor(df_orov_melt$Group2, levels = c("HC", "1-2 days", "3-4 days"))

# Stats
wilcox_results_days <- df_orov_melt[df_orov_melt$is_outlier == FALSE, ] %>%
  group_by(Cytokine) %>%
  nest() %>%
  mutate(comparisons = map(data, ~ {
    
    # Get the 3 groups
    grps <- unique(.x$Group2)
    combs <- combn(grps, 2, simplify = FALSE)
    
    # Perform pairwise tests
      map_df(combs, function(cmb) {
        g1 <- cmb[1]
        g2 <- cmb[2]
        
        df_sub <- .x %>% filter(Group2 %in% c(g1, g2))
        
        wt <- tryCatch(wilcox.test(Concentration ~ Group2,
                                   data = df_sub,
                                   exact = FALSE,
                                   conf.int = TRUE),
                       error = function(e) NULL)
        
        tibble(group1 = g1,
               group2 = g2,
               p_value = if (!is.null(wt)) wt$p.value else NA_real_,
               hl_estimate = if (!is.null(wt)) wt$estimate else NA_real_,
               hl_ci_low  = if (!is.null(wt)) wt$conf.int[1] else NA_real_,
               hl_ci_high = if (!is.null(wt)) wt$conf.int[2] else NA_real_,
               n_group1 = sum(df_sub$Group2 == g1),
               n_group2 = sum(df_sub$Group2 == g2))
      })
    })
  ) %>%
  select(-data) %>%
  unnest(comparisons) %>%
  group_by(Cytokine) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Coding stars
wilcox_results_days <- pvalueCodes(wilcox_results_days, "p_adj")

# Max Y for each cytokine (per panel)
y_lims <- df_orov_melt[df_orov_melt$is_outlier==FALSE,] %>%
  group_by(Cytokine) %>%
  summarise(ymax = max(Concentration, na.rm = TRUE))

# Add bracket positions + significance
stars_df_days <- wilcox_results_days %>%
  left_join(y_lims, by = "Cytokine") %>%
  mutate(
    # Offset the bracket height so multiple comparisons don’t overlap
    bracket_id = match(paste0(group1, group2), unique(paste0(group1, group2))),
    y_bracket = ymax * (1 + 0.12 * bracket_id),
    
    # X positions for bracket
    Group_left  = as.numeric(factor(group1, levels = sort(unique(df_orov_melt$Group2)))),
    Group_right = as.numeric(factor(group2, levels = sort(unique(df_orov_melt$Group2)))),
    
    # Code p-values
    label = case_when(p_adj < 0.001 ~ "***",
                      p_adj < 0.01  ~ "**",
                      p_adj < 0.05  ~ "*",
                      TRUE          ~ "ns"),
    y = y_bracket * 1.05)

# Remove ns
stars_df_days <- stars_df_days[stars_df_days$pcode!="",]

### Plotting and saving
tiff("results/cytokine_days.tiff", res = 300, w = 2600, h = 2500, units = "px")

ggplot(df_orov_melt[df_orov_melt$is_outlier==FALSE,],
       aes(x = Group2, y = Concentration)) +
  facet_wrap(~ Cytokine, scales = "free_y") +
  geom_quasirandom(size = 1.8, alpha = 0.7) +
  geom_boxplot(outliers = FALSE, width = 0.45, alpha = 0.75) +
  # significance stars
  geom_text(
    data = stars_df_days,
    aes(x = (Group_left + Group_right)/2, y = y, label = pcode),
    size = 5, vjust = 0.5
  ) +
  # bracket lines
  geom_segment(
    data = stars_df_days,
    aes(x = Group_left, xend = Group_right,
        y = y_bracket, yend = y_bracket),
    linewidth = 0.4
  ) +
  scale_y_continuous(expand = c(0, 0, 0.15, 0)) +
  labs(x = "", y = "Concentration (pg/mL)") +
  theme_bw() +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14, face = "bold"), 
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

dev.off()

### Saving results table
write.xlsx(wilcox_results_days, "results/results.xlsx", sheetName = "Wilcoxon - Days", showNA = FALSE, append = TRUE)

### ----- CT/days partial correlation ------------------------------------------

# Extracting colnames
cyto_cols <- colnames(df_orov_pt)[6:18]

# Running bootstrap pcor
res_boot <- boot_partial_spearman(df_orov_pt, 1000, c("CT", "Days"), cyto_cols)
res_boot <- adjust_pvalue(res_boot, "p_val", method = "BH")
res_boot <- pvalueCodes(res_boot, "p_val.adj")
res_boot$label <- sprintf("ρ = %.2f, ns\n95%% CI: %.2f – %.2f",
                          res_boot$correlation,
                          res_boot$ci_low,
                          res_boot$ci_high)

# Renaming column to avoid plot repetition of text
res_boot <- res_boot %>%
  rename(Cytokine = variable)

### Saving other table
write.xlsx(res_boot, "results/results.xlsx", sheetName = "CT correlations", showNA = FALSE, append = TRUE)
write.xlsx(heat_df, "results/results.xlsx", sheetName = "Difference Network", showNA = FALSE, append = TRUE)

### Plotting and Saving
tiff("results/cor_cyto_ct.tiff", res = 300, w = 3000, h = 2500, type = "cairo", units = "px")

ggplot(df_orov_melt_pt, aes(x = CT, y = Concentration)) +
  facet_wrap(~ Cytokine, scales = "free") +
  geom_point(size = 1, alpha = 0.5) +
  geom_text(data = res_boot, 
            aes(label = label, x = 16, y = Inf),
            hjust = 0, vjust = 1.1, size = 3) +
  geom_smooth(method = lm, alpha = 0.2, fill = "#FF2345", colour = "#FF2345") +
  theme_bw() +
  theme(text = element_text(colour = "black"),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"), 
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

dev.off()
