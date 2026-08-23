library(tidyverse)

SIMULATION_ID = "lowdim1" # commandArgs(trailingOnly=TRUE)[1]
COMPONENTS = 2 # as.numeric(commandArgs(trailingOnly=TRUE)[2])
RESULT_PATH = paste0("simulation_R1_SUBMISSION_", COMPONENTS, "_components_", SIMULATION_ID)
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

methods    = c("mvHE", "mvREHE", "mvREML", "GEMMA", "HE", "REHE", "REML")
mv_methods = c("mvHE", "mvREHE", "mvREML")

palette = ggsci::pal_aaas("default")(4)
names(palette) = c("mvHE", "mvREHE", "mvREML", "GEMMA")
palette = c(palette, palette[1], palette[2], palette[3])
names(palette) = c("mvHE", "mvREHE", "mvREML", "GEMMA", "HE-diag", "REHE-diag", "REML-diag")

map = c("mvHE", "mvREHE", "mvREML", "GEMMA", "HE-diag", "REHE-diag", "REML-diag")
names(map) = methods

# Encoding (matches the colour version of the main data figure):
#   colour   -> estimator family (mv and diag versions share a colour)
#   linetype -> multivariate/GEMMA (solid) vs diag (dashed)
#   shape    -> regression method (regression figure only)
mvdiag_linetypes = c(
  "mvHE" = "solid", "mvREHE" = "solid", "mvREML" = "solid",
  "GEMMA" = "solid", "HE-diag" = "dashed", "REHE-diag" = "dashed", "REML-diag" = "dashed"
)

reg_shapes = c("tensor" = 16, "lasso" = 17, "ridge" = 15)

files   = list.files(RESULT_PATH, full.names = TRUE)
files   = files[grepl("rds", files)]
results <- lapply(files, function(f) {
  tryCatch(
    readRDS(f),
    error = function(e) NULL  # Returns NULL if reading fails
  )
})

label = c("hat(Sigma)[G]", "hat(Sigma)[C]", "hat(Sigma)[E]")
names(label) = c("Sigma_1", "Sigma_2", "Sigma_0")

q_levels = c("q = 5", "q = 10", "q = 20")

### ── Time ───────────────────────────────────────────────────────────────────

time_df = do.call(rbind, lapply(results, function(x) x$time)) %>%
  group_by(n, q, experiment, method) %>%
  summarize(mean = mean(na.rm = TRUE, time), se = sd(na.rm = TRUE, time) / sqrt(n())) %>%
  filter(gsub("\\.elapsed$", "", method) %in% methods) %>%
  mutate(
    method_label = factor(map[gsub("\\.elapsed$", "", method)], levels = names(palette)),
    facet = factor(paste0("q = ", q), levels = q_levels)
  ) %>%
  filter(!is.na(method_label), experiment == "n")

time_plot = ggplot(time_df %>% filter(n < 10000),
  aes(x = n, y = log10(mean), ymin = log10(mean - 1.96 * se), ymax = log10(mean + 1.96 * se),
      color = method_label, linetype = method_label, group = method)
) +
  facet_wrap(~facet, nrow = 1) +
  geom_line() + geom_errorbar(width = 0.1) +
  theme_bw() +
  labs(color = "Method", linetype = "Method", x = NULL, y = "log10(seconds)") +
  theme(legend.position = "none", strip.background = element_blank()) +
  scale_color_manual(values = palette) +
  scale_linetype_manual(values = mvdiag_linetypes) +
  guides(color = "none", linetype = "none")

### ── Spectral error (Sigma_G) ───────────────────────────────────────────────

spectral_error_df = do.call(rbind, lapply(results, function(x) x$spectral_error)) %>%
  mutate(estimate = label[as.character(estimate)]) %>%
  filter(method %in% methods)

spectral_plot = ggplot(
  spectral_error_df %>%
    filter(experiment == "n",  n < 10000, method %in% mv_methods, estimate == label[1]) %>%
    mutate(
      method_label = factor(map[method], levels = names(palette)),
      facet = factor(paste0("q = ", q), levels = q_levels)
    ) %>%
    group_by(n, facet, method, method_label) %>%
    summarize(mean = mean(na.rm = TRUE, spectral_error), se = sd(na.rm = TRUE, spectral_error) / sqrt(n())),
  aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
      color = method_label, linetype = method_label, group = method_label)
) +
  facet_wrap(~facet, nrow = 1, scales = "free_y") +
  geom_line() + geom_errorbar(width = 0.1) +
  theme_bw() +
  labs(color = "Method", linetype = "Method", x = NULL,
       y = expression("||"*hat(Sigma)[G] - Sigma[G]*"||"[2])) +
  theme(legend.position = "none", strip.background = element_blank()) +
  scale_color_manual(values = palette) +
  scale_linetype_manual(values = mvdiag_linetypes) +
  guides(color = "none", linetype = "none")

### ── h2 error ───────────────────────────────────────────────────────────────

h2_df = do.call(rbind, lapply(results, function(x) x$h2_error)) %>%
  group_by(n, q, experiment, method) %>%
  summarize(mean = mean(na.rm = TRUE, h2_error), se = sd(na.rm = TRUE, h2_error) / sqrt(n())) %>%
  filter(method %in% methods)

h2_plot = ggplot(
  h2_df %>%
    filter(experiment == "n", n < 10000) %>%
    mutate(
      method_label = factor(map[method], levels = names(palette)),
      facet = factor(paste0("q = ", q), levels = q_levels)
    ) %>%
    filter(!is.na(method_label)),
  aes(x = n, y = sqrt(mean), ymin = sqrt(mean - 1.96 * se), ymax = sqrt(mean + 1.96 * se),
      color = method_label, linetype = method_label, group = method)
) +
  facet_wrap(~facet, nrow = 1, scales = "free_y") +
  geom_line() + geom_errorbar(width = 0.1) +
  theme_bw() +
  labs(color = "Method", linetype = "Method", x = "n",
       y = expression("||"*hat(h)^2 - h^2*"||"[2])) +
  theme(legend.position = "bottom", strip.background = element_blank()) +
  scale_color_manual(values = palette) +
  scale_linetype_manual(values = mvdiag_linetypes) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

### ── Supplement figure S7: time | Sigma_G | h2, stacked ───────────────────

patchwork::wrap_plots(list(time_plot, spectral_plot, h2_plot), ncol = 1) +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c"))) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_supplement_lowdim_n.pdf"), height = 7.5, width = 8.5)

### ── Regression figure S12: r2_beta_hat at q = 5 ──────────────────────────

r2_beta_hat = lapply(results, function(x) x$r2_beta_hat)
r2_beta_hat = r2_beta_hat[lengths(sapply(r2_beta_hat, colnames)) == 10]
r2_beta_hat = do.call(rbind, r2_beta_hat)
colnames(r2_beta_hat)[10] = "regression_method"
r2_beta_hat$estimate = factor(label[r2_beta_hat$estimate], levels = label)

r2_beta_hat_per_method = r2_beta_hat %>%
  group_by(n, estimate, method, regression_method) %>%
  summarize(r2 = mean(r2, trim =)) %>%
  mutate(
    method_label = factor(map[method], levels = names(palette)),
    regression_method = factor(regression_method, levels = c("ridge", "lasso"))
  ) %>%
  filter(!is.na(method_label))

ggplot(
  r2_beta_hat_per_method %>% filter(estimate == label[1]),
  aes(x = n, y = r2, color = method_label, linetype = method_label,
      shape = regression_method,
      group = interaction(method, regression_method))
) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = palette) +
  scale_linetype_manual(values = mvdiag_linetypes) +
  scale_shape_manual(values = reg_shapes) +
  labs(color = "Method", linetype = "Method",
       shape = "Regression Method", x = "n",
       y = expression(frac(R[G]^"2*" * "(" * hat(beta)[G] * ")", R[G]^"2*" * "(" * beta[G]^"*" * ")"))) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(-0.1, 1))
ggsave(file.path(FIGURES_PATH, "simulation_figure_supplement_regression_n.pdf"), height = 3, width = 5)

