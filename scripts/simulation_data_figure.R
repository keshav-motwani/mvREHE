library(tidyverse)

SIMULATION_ID = "data"
COMPONENTS = 3 # as.numeric(commandArgs(trailingOnly=TRUE)[1])
RESULT_PATH = paste0("simulation_R1_SUBMISSION_", COMPONENTS, "_components_", SIMULATION_ID)
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

methods = c("mvHE", "mvREHE", "HE", "REHE", "REML")
map = c("mvHE", "mvREHE", "HE-diag", "REHE-diag", "REML-diag")
names(map) = methods

# Colour palette for the estimator family
palette = ggsci::pal_aaas("default")(3)
palette = c(palette[1], palette[2], palette[1], palette[2], palette[3])
names(palette) = map

# Two encodings of the combined figure (each aesthetic encodes exactly one variable):
#   COLOUR version (every other figure follows this one too):
#     colour -> estimator family   linetype -> mv (solid) vs diag (dashed)   shape -> regression method
#   B&W version:
#     colour -> mv (black) vs diag (grey)   shape -> estimator family   linetype -> regression method
method_shapes    = c("mvHE" = 16, "mvREHE" = 17, "HE-diag" = 16, "REHE-diag" = 17, "REML-diag" = 15)
mvdiag_linetypes = c("mvHE" = "solid", "mvREHE" = "solid",
                     "HE-diag" = "dashed", "REHE-diag" = "dashed", "REML-diag" = "dashed")
bw_colors        = c("mvHE" = "black", "mvREHE" = "black",
                     "HE-diag" = "grey70", "REHE-diag" = "grey70", "REML-diag" = "grey70")
reg_shapes       = c("tensor" = 16, "lasso" = 17, "ridge" = 15)
reg_linetypes    = c("tensor" = "solid", "lasso" = "dashed", "ridge" = "dotted")

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]

results = lapply(files, readRDS)

label = c("hat(Sigma)[G]", "hat(Sigma)[C]", "hat(Sigma)[E]")
names(label) = c("Sigma_1", "Sigma_2", "Sigma_0")

### ── Time (combined panels only) ─────────────────────────────────────────────

time_df = do.call(rbind, lapply(results, function(x) x$time)) %>%
  group_by(n, q, Sigma, experiment, method) %>%
  summarize(mean = mean(na.rm = TRUE, time), se = sd(na.rm = TRUE, time) / sqrt(n())) %>%
  filter(method %in% methods) %>%
  filter(experiment == "n") %>%
  mutate(y = log10(mean),
         ymax = log10(mean + 1.96 * se),
         ymin = log10(mean - 1.96 * se),
         method_label = factor(map[method], levels = names(method_shapes)))

# colour: colour = family, linetype = mv/diag
time_col = ggplot(time_df,
       aes(x = n, y = y, ymax = ymax, ymin = ymin,
           color = method_label, linetype = method_label, group = method_label)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "Method", linetype = "Method", y = "log10(seconds)") +
  scale_color_manual(values = palette) +
  scale_linetype_manual(values = mvdiag_linetypes)

# b&w: colour = mv/diag, shape = family
time_bw = ggplot(time_df,
       aes(x = n, y = y, ymax = ymax, ymin = ymin,
           color = method_label, shape = method_label, group = method_label)) +
  geom_line() +
  geom_point() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "Method", shape = "Method", y = "log10(seconds)") +
  scale_color_manual(values = bw_colors) +
  scale_shape_manual(values = method_shapes)

### ── Spectral error → simulation_figure_spectral_error_n.pdf (colour version) ──

spectral_error_df = do.call(rbind, lapply(results, function(x) x$spectral_error)) %>%
  mutate(estimate = label[as.character(estimate)]) %>%
  filter(method %in% methods)  %>%
  filter(experiment == "n" & grepl("mv", method)) %>%
  mutate(facet = factor(estimate, levels = label)) %>%
  group_by(n, facet, method) %>%
  summarize(mean = mean(na.rm = TRUE, spectral_error), se = sd(na.rm = TRUE, spectral_error) / sqrt(n()))

ggplot(spectral_error_df,
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(method, levels = names(method_shapes)),
           linetype = factor(method, levels = names(method_shapes)),
           group = method)) +
  facet_wrap(~facet, scales = "free_y", ncol = COMPONENTS, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  labs(color = "Method", linetype = "Method",
       y = expression("||"*hat(Sigma)[k] - Sigma[k]*"||"[2]), x = "n") +
  theme(legend.position = "bottom", strip.background = element_blank(), strip.placement = "outside") +
  scale_color_manual(values = palette) +
  scale_linetype_manual(values = mvdiag_linetypes)
ggsave(file.path(FIGURES_PATH, "simulation_figure_spectral_error_n.pdf"), height = 3.5, width = 8.5)

# Sigma_G-only panels for the combined figure
spectral_G = spectral_error_df %>% filter(facet == label[1])

spectral_col = ggplot(spectral_G,
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(method, levels = names(method_shapes)),
           linetype = factor(method, levels = names(method_shapes)),
           group = method)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  labs(y = expression("||"*hat(Sigma)[G] - Sigma[G]*"||"[2])) +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside") +
  scale_color_manual(values = palette) +
  scale_linetype_manual(values = mvdiag_linetypes) +
  guides(color = "none", linetype = "none")

spectral_bw = ggplot(spectral_G,
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(method, levels = names(bw_colors)),
           shape = factor(method, levels = names(method_shapes)),
           group = method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  labs(y = expression("||"*hat(Sigma)[G] - Sigma[G]*"||"[2])) +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside") +
  scale_color_manual(values = bw_colors) +
  scale_shape_manual(values = method_shapes) +
  guides(color = "none", shape = "none")

### ── R2 → simulation_figure_r2_n.pdf (colour version) ─────────────────────────

r2_beta_hat = lapply(results, function(x) x$r2_beta_hat)
r2_beta_hat = r2_beta_hat[lengths(sapply(r2_beta_hat, colnames)) == 10]
r2_beta_hat = do.call(rbind, r2_beta_hat)
colnames(r2_beta_hat)[10] = "regression_method"
r2_beta_hat$estimate = factor(label[r2_beta_hat$estimate], levels = label)

r2_beta_hat_per_method = r2_beta_hat %>%
  group_by(n, estimate, method, regression_method) %>%
  summarize(r2 = mean(r2, trim = 0.05)) %>%
  mutate(regression_method = factor(recode(regression_method, matrix = "tensor"),
                                    levels = c("tensor", "lasso", "ridge")))

ggplot(r2_beta_hat_per_method,
       aes(x = n, y = r2,
           color = factor(map[method], levels = names(method_shapes)),
           linetype = factor(map[method], levels = names(method_shapes)),
           shape = regression_method,
           group = interaction(method, regression_method))) +
  facet_wrap(~estimate, scales = "free", labeller = labeller(estimate = label_parsed)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        strip.background = element_blank(), strip.placement = "outside") +
  scale_color_manual(values = palette) +
  scale_linetype_manual(values = mvdiag_linetypes) +
  scale_shape_manual(values = reg_shapes) +
  labs(color = "Method", linetype = "Method", shape = "Regression Method") +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  ylab(expression(
    "Average " *
      frac(R[kj]^"2*" * "(" * hat(beta)[kj] * ")", R[kj]^"2*" * "(" * beta[kj]^"*" * ")") * " across " * j %in% F
  ))
ggsave(file.path(FIGURES_PATH, "simulation_figure_r2_n.pdf"), height = 3.5, width = 8.5)

# Sigma_G-only panels for the combined figure
r2_G = r2_beta_hat_per_method %>% filter(estimate == label[1])

r2_col = ggplot(r2_G,
    aes(x = n, y = r2,
        color = factor(map[method], levels = names(method_shapes)),
        linetype = factor(map[method], levels = names(method_shapes)),
        shape = regression_method,
        group = interaction(method, regression_method))) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = palette) +
  scale_linetype_manual(values = mvdiag_linetypes) +
  scale_shape_manual(values = reg_shapes) +
  ylab(expression(
    "Average " * frac(R[Gj]^"2*" * "(" * hat(beta)[Gj] * ")", R[Gj]^"2*" * "(" * beta[Gj]^"*" * ")") * " across " * j %in% F
  )) +
  guides(color = "none", linetype = "none") +
  labs(shape = "Regression Method")

r2_bw = ggplot(r2_G,
    aes(x = n, y = r2,
        color = factor(map[method], levels = names(bw_colors)),
        shape = factor(map[method], levels = names(method_shapes)),
        linetype = regression_method,
        group = interaction(method, regression_method))) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = bw_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = reg_linetypes) +
  ylab(expression(
    "Average " * frac(R[Gj]^"2*" * "(" * hat(beta)[Gj] * ")", R[Gj]^"2*" * "(" * beta[Gj]^"*" * ")") * " across " * j %in% F
  )) +
  guides(color = "none", shape = "none") +
  labs(linetype = "Regression Method")

### ── Combined figure → colour + B&W versions ─────────────────────────────────

combine = function(a, b, c) {
  patchwork::wrap_plots(a, b, c) +
    patchwork::plot_annotation(tag_levels = list(c("a", "b", "c"))) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# Colour: colour = family, linetype = mv/diag, shape = regression method
combine(time_col, spectral_col, r2_col)
ggsave(file.path(FIGURES_PATH, "simulation_figure_combined_n.pdf"), height = 3.5, width = 10)

# B&W: colour = mv/diag, shape = family, linetype = regression method
combine(time_bw, spectral_bw, r2_bw)
ggsave(file.path(FIGURES_PATH, "simulation_figure_combined_n_bw.pdf"), height = 3.5, width = 10)
