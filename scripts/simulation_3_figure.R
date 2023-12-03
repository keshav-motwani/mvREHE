library(tidyverse)

SIMULATION_ID = 3
RESULT_PATH = paste0("simulation_hcp_results_", SIMULATION_ID)
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

methods = c("mvHE", "mvREHE", "mvHE-smoothed", "mvREHE-smoothed")
palette = ggsci::pal_aaas("default")(length(methods))
names(palette) = methods
options(ggplot2.discrete.colour = palette)
map = methods
names(map) = methods

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]
files = files[grepl("_rInf_", files)]

results = lapply(files, readRDS)

Sigma_labels = c("alpha == 1", "alpha == 2")
names(Sigma_labels) = c("smooth_1", "smooth_2")

for (exp in c("n", "q")) {

  ### Time

  time_df = do.call(rbind, lapply(results, function(x) x$time)) %>%
    group_by(n, q, Sigma, experiment, method) %>%
    summarize(time = mean(time)) %>%
    filter(method %in% methods)

  ggplot(time_df %>%
           filter(experiment == exp) %>%
           mutate(facet = Sigma_labels[Sigma]) %>%
           mutate(facet = factor(facet, levels = Sigma_labels)),
         aes(x = get(exp), y = log10(time), color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"), group = method)) +
    facet_wrap(~facet, scales = "free_y", nrow = 1, labeller = labeller(facet = label_parsed)) +
    geom_line() +
    theme_bw() +
    xlab(exp) +
    labs(color = "Method", linetype = "Method", y = "log10(seconds)") +
    theme(legend.position = "bottom")
  ggsave(file.path(FIGURES_PATH, paste0("simulation_figure_time_", exp, ".pdf")), height = 3, width = 8.5)

  ### Spectral error

  label = c("hat(Sigma)[0]", "hat(Sigma)[1]")
  names(label) = c("Sigma_0", "Sigma_1")

  Sigmas = apply(expand.grid(label, Sigma_labels), 1, function(x) paste0(x[2], "~(", x[1], ")"))

  spectral_error_df = do.call(rbind, lapply(results, function(x) x$spectral_error)) %>%
    mutate(estimate = label[estimate]) %>%
    filter(method %in% methods)

  ggplot(spectral_error_df %>%
           filter(experiment == exp & grepl("mv", method)) %>%
           mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
           mutate(facet = factor(facet, levels = Sigmas)),
         aes(x = as.factor(get(exp)), y = spectral_error, color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)))) +
    facet_wrap(~facet, scales = "free_y", nrow = 2, dir = "v", labeller = labeller(facet = label_parsed)) +
    geom_boxplot(outlier.size = 0) +
    theme_bw() +
    xlab(exp) +
    labs(color = "Method", linetype = "Method", y = expression("||"*hat(Sigma)[k] - Sigma[k]*"||"[2]), x = "n") +
    theme(legend.position = "bottom")
  ggsave(file.path(FIGURES_PATH, paste0("simulation_figure_spectral_error_", exp, ".pdf")), height = 5, width = 8.5)

  ### Squared error

  squared_error_df = do.call(rbind, lapply(results, function(x) x$squared_error)) %>%
    mutate(estimate = label[estimate]) %>%
    filter(method %in% methods)

  ggplot(squared_error_df %>%
           filter(experiment == exp & grepl("mv", method)) %>%
           mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
           mutate(facet = factor(facet, levels = Sigmas)),
         aes(x = as.factor(get(exp)), y = squared_error, color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)))) +
    facet_wrap(~facet, scales = "free_y", nrow = 2, dir = "v", labeller = labeller(facet = label_parsed)) +
    geom_boxplot(outlier.size = 0) +
    theme_bw() +
    xlab(exp) +
    labs(color = "Method", linetype = "Method", y = expression("||"*hat(Sigma)[k] - Sigma[k]*"||"["F"]), x = "n") +
    theme(legend.position = "bottom")
  ggsave(file.path(FIGURES_PATH, paste0("simulation_figure_squared_error_", exp, ".pdf")), height = 5, width = 8.5)

  ### Squared error on diagonal

  diag_squared_error_df = do.call(rbind, lapply(results, function(x) x$diag_squared_error)) %>%
    mutate(estimate = label[estimate]) %>%
    filter(method %in% methods)

  ggplot(diag_squared_error_df %>%
           filter(experiment == exp) %>%
           mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
           mutate(facet = factor(facet, levels = Sigmas)),
         aes(x = as.factor(get(exp)), y = diag_squared_error, color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"))) +
    facet_wrap(~facet, scales = "free_y", nrow = 2, dir = "v", labeller = labeller(facet = label_parsed)) +
    geom_boxplot(outlier.size = 0) +
    theme_bw() +
    xlab(exp) +
    labs(color = "Method", linetype = "Method", y = expression("||diag("*hat(Sigma)[k] - Sigma[k]*")||"[2]), x = "n") +
    theme(legend.position = "bottom")
  ggsave(file.path(FIGURES_PATH, paste0("simulation_figure_diag_squared_error_", exp, ".pdf")), height = 5, width = 8.5)

}
