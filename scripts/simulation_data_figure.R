library(tidyverse)

SIMULATION_ID = "data"
COMPONENTS = 3 # as.numeric(commandArgs(trailingOnly=TRUE)[1])
RESULT_PATH = paste0("simulation_FINAL_FINAL_FINAL_ALL_", COMPONENTS, "_components_", SIMULATION_ID)
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

methods = c("mvHE", "mvREHE", "HE", "REHE", "REML")
map = gsub("mv", "", methods)
names(map) = methods

palette = ggsci::pal_aaas("default")(3)
names(palette) = c("HE", "REHE", "REML")

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]

results = lapply(files, readRDS)

time_df = do.call(rbind, lapply(results, function(x) x$time)) %>%
  group_by(n, q, Sigma, experiment, method) %>%
  summarize(mean = mean(na.rm = TRUE, time), se = sd(na.rm = TRUE, time) / sqrt(n())) %>%
  filter(method %in% methods) %>%
  filter(experiment == "n") %>%
  mutate(y = log10(mean),
         ymax = log10(mean + 1.96 * se),
         ymin = log10(mean - 1.96 * se),
         color = factor(map[method], levels = names(palette)),
         linetype = factor(ifelse(grepl("mv", method), "Multivariate", "Univariate")))

time = ggplot(time_df,
       aes(x = n, y = y, ymax = ymax, ymin = ymin, color = color, linetype = linetype)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "Method", linetype = "", y = "log10(seconds)") +
  scale_color_manual(values = palette)

label = c("hat(Sigma)[G]", "hat(Sigma)[C]", "hat(Sigma)[E]")
names(label) = c("Sigma_1", "Sigma_2", "Sigma_0")

Sigma_list_truth = readRDS("data_analysis_3_components_10ROI/fit.rds")$Sigma_hat
indices = which(diag(Sigma_list_truth[[1]]) > 1e-14)
h2_truth = sapply(indices, function(j) Sigma_list_truth[[2]][j, j] / sum(sapply(Sigma_list_truth, function(Sigma) Sigma[j, j])))

h2 = lapply(results, function(result) {
  Sigma_list_estimate = result$output$estimate$Sigma_hat
  h2_estimate = sapply(indices, function(j)
    Sigma_list_estimate[[2]][j, j] / sum(sapply(Sigma_list_estimate, function(Sigma)
      Sigma[j, j])))
  cbind(result$h2_error[, -2], h2_error = sqrt(sum(h2_truth - h2_estimate)^2))
})

h2_df = do.call(rbind, h2) %>%
  group_by(n, q, Sigma, experiment, method) %>%
  summarize(mean = mean(na.rm = TRUE, h2_error), se = sd(na.rm = TRUE, h2_error) / sqrt(n())) %>%
  filter(method %in% methods) %>%
  filter(experiment == "n")

ggplot(h2_df,
       aes(x = n, y = sqrt(mean), ymin = sqrt(mean - 1.96 * se), ymax = sqrt(mean + 1.96 * se),
           color = factor(map[method], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"), group = method)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(color = "Method", linetype = "", y = expression("||"*hat(h)^2 - h^2*"||"[2])) +
  scale_color_manual(values = palette)

ggsave(file.path(FIGURES_PATH, "simulation_figure_h2_n.pdf"), height = 3.5, width = 6)

palette_mv = palette
names(palette_mv) = paste0("mv", names(palette_mv))

spectral_error_df = do.call(rbind, lapply(results, function(x) x$spectral_error)) %>%
  mutate(estimate = label[as.character(estimate)]) %>%
  filter(method %in% methods)  %>%
  filter(experiment == "n" & grepl("mv", method)) %>%
  mutate(facet = estimate) %>%
  mutate(facet = factor(facet, levels = label)) %>%
  group_by(n, facet, method) %>%
  summarize(mean = mean(na.rm = TRUE, spectral_error), se = sd(na.rm = TRUE, spectral_error) / sqrt(n()))

ggplot(spectral_error_df,
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(method, levels = names(palette_mv)))) +
  facet_wrap(~facet, scales = "free_y", ncol = COMPONENTS, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", y = expression("||"*hat(Sigma)[k] - Sigma[k]*"||"[2]), x = "n") +
  theme(legend.position = "bottom", strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(NA, NA)) +
  scale_color_manual(values = palette_mv)
ggsave(file.path(FIGURES_PATH, "simulation_figure_spectral_error_n.pdf"), height = 3.5, width = 8.5)


spectral = ggplot(spectral_error_df %>%
         filter(facet == label[1]),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(method, levels = names(palette_mv)))) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  labs(color = "Method", y = expression("||"*hat(Sigma)[G] - Sigma[G]*"||"[2])) +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(NA, NA)) +
  scale_color_manual(values = palette_mv) +
  guides(linetype="none", color = "none")

pseudoinverse = function(mat) {
  eig = eigen(mat)
  r = sum(eig$val > 1e-12)
  eig$vec[, 1:r] %*% diag(1/eig$val[1:r]) %*% t(eig$vec[, 1:r])
}

cov2cor_NA0 = function(cov) {
  cov_hat = cov2cor(cov)
  cov_hat[diag(cov) < .Machine$double.eps, ] = 0
  cov_hat[, diag(cov) < .Machine$double.eps] = 0
  cov_hat
}

fit = readRDS(file.path("data_analysis_3_components_10ROI/", "fit.rds"))
outcomes = setdiff(1:55, c(1, cumsum(10:1) + 1))
covariates = 56:110

r2_true = sapply(1:length(fit$Sigma_hat), function(k) {
  Sigma = cov2cor_NA0(fit$Sigma_hat[[k]])
  sapply(outcomes, function(outcome) {
    beta = pseudoinverse(Sigma[covariates, covariates]) %*% Sigma[outcome, covariates]
    r2 = 1 - (Sigma[outcome, outcome] -
                2 * Sigma[outcome, covariates] %*% beta +
                t(beta) %*% Sigma[covariates, covariates] %*% beta) /
      Sigma[outcome, outcome]
    r2
  })
})
r2_true = cbind(reshape2::melt(r2_true), outcome = outcomes)
r2_true$estimate = paste0("Sigma_", r2_true$Var2 - 1)
r2_true$r2_true = r2_true$value
r2_true = r2_true[, 4:6]

r2_beta_hat = lapply(results, function(x) x$r2_beta_hat)
r2_beta_hat = r2_beta_hat[lengths(sapply(r2_beta_hat, colnames)) == 11]
r2_beta_hat = do.call(rbind, r2_beta_hat)
r2_beta_hat = left_join(r2_beta_hat, r2_true, by = c("estimate", "outcome")) %>% mutate(r2 = r2/r2_true) 
r2_beta_hat$estimate = factor(label[r2_beta_hat$estimate], levels = label)

r2_beta_hat_per_method = r2_beta_hat %>%
  group_by(n, estimate, method, regression_method, outcome) %>%
  summarize(r2 = mean(r2)) %>% 
  group_by(n, estimate, method, regression_method) %>% 
  summarize(r2 = mean(r2)) %>%
  mutate(regression_method = factor(regression_method, levels = c("matrix", "lasso", "ridge")))

ggplot(r2_beta_hat_per_method,
       aes(
         x = n,
         y = r2,
         color = method,
         shape = regression_method
       )) +
  facet_wrap(~estimate,
             scales = "free",
             labeller = labeller(estimate = label_parsed)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.placement = "outside"
  ) +
  scale_color_manual(values = palette_mv) +
  labs(color = "Method", shape = "Regression Method") +
  ylab(expression(
    "Average " *
      frac(R[kj]^"2*" * "(" * hat(beta)[kj] * ")", R[kj]^"2*" * "(" * beta[kj]^"*" * ")") * " across " * j %in% F
  ))
ggsave(file.path(FIGURES_PATH, "simulation_figure_r2_n.pdf"), height = 3.5, width = 8.5)


r2 = ggplot(
  r2_beta_hat_per_method %>%
    filter(estimate == label[1]),
  aes(
    x = n,
    y = r2,
    color = method,
    shape = regression_method
  )
) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = palette_mv) +
  ylab(expression(
    "Average " * frac(R[Gj]^"2*" * "(" * hat(beta)[Gj] * ")", R[Gj]^"2*" * "(" * beta[Gj]^"*" * ")") * " across " * j %in% F
  )) +
  guides(color = "none") +
  labs(shape = "Regression Method")

patchwork::wrap_plots(time, spectral, r2) +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c"))) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position =  "bottom")

ggsave(file.path(FIGURES_PATH, "simulation_figure_combined_n.pdf"), height = 3.5, width = 10)
