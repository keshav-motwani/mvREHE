library(tidyverse)

SIMULATION_ID = 2
RESULT_PATH = paste0("final_simulation_hcp_results_", SIMULATION_ID)
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

methods = c("mvREHE", "mvREHE_cvL2", "mvREHE_cvDR", paste0("DR", c(5), "-mvREML"))
palette = ggsci::pal_aaas("default")(length(methods))
names(palette) = methods
options(ggplot2.discrete.colour = palette)
map = methods
names(map) = methods

Sigmas = c("fast", "moderate", "slow")

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]
files = files[grepl("_rInf_", files)]

results = lapply(files, readRDS)

### Time

time_df = do.call(rbind, lapply(results, function(x) x$time)) %>%
  group_by(n, q, Sigma, experiment, method) %>%
  summarize(time = mean(time)) %>%
  filter(method %in% methods)

ggplot(time_df %>%
         filter(experiment == "n") %>%
         mutate(facet = Sigma) %>%
         mutate(facet = factor(facet, levels = Sigmas)),
       aes(x = n, y = log10(time), color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"), group = method)) +
  facet_wrap(~facet, scales = "free_y", nrow = 1) +
  geom_line() +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "Method", y = "log10(seconds)") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_time_n.pdf"), height = 3, width = 8.5)

### Spectral error

label = c("hat(Sigma)[E]", "hat(Sigma)[G]")
names(label) = c("Sigma_0", "Sigma_1")

Sigmas = apply(expand.grid(label, Sigmas), 1, function(x) paste0(x[2], "~(", x[1], ")"))

spectral_error_df = do.call(rbind, lapply(results, function(x) x$spectral_error)) %>%
  mutate(estimate = label[estimate]) %>%
  filter(method %in% methods)

ggplot(spectral_error_df %>%
         filter(experiment == "n" & grepl("mv", method)) %>%
         mutate(facet = paste0(Sigma, "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = Sigmas)),
       aes(x = as.factor(n), y = spectral_error, color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)))) +
  facet_wrap(~facet, scales = "free_y", nrow = 2, dir = "v", labeller = labeller(facet = label_parsed)) +
  geom_boxplot(outlier.size = 0) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "Method", y = expression("||"*hat(Sigma)[k] - Sigma[k]*"||"[2]), x = "n") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_spectral_error_n.pdf"), height = 5, width = 8.5)

### Squared error

squared_error_df = do.call(rbind, lapply(results, function(x) x$squared_error)) %>%
  mutate(estimate = label[estimate]) %>%
  filter(method %in% methods)

ggplot(squared_error_df %>%
         filter(experiment == "n" & grepl("mv", method)) %>%
         mutate(facet = paste0(Sigma, "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = Sigmas)),
       aes(x = as.factor(n), y = squared_error, color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)))) +
  facet_wrap(~facet, scales = "free_y", nrow = 2, dir = "v", labeller = labeller(facet = label_parsed)) +
  geom_boxplot(outlier.size = 0) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "Method", y = expression("||"*hat(Sigma)[k] - Sigma[k]*"||"["F"]), x = "n") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_squared_error_n.pdf"), height = 5, width = 8.5)

### Squared error on diagonal

diag_squared_error_df = do.call(rbind, lapply(results, function(x) x$diag_squared_error)) %>%
  mutate(estimate = label[estimate]) %>%
  filter(method %in% methods)

ggplot(diag_squared_error_df %>%
         filter(experiment == "n") %>%
         mutate(facet = paste0(Sigma, "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = Sigmas)),
       aes(x = as.factor(n), y = diag_squared_error, color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"))) +
  facet_wrap(~facet, scales = "free_y", nrow = 2, dir = "v", labeller = labeller(facet = label_parsed)) +
  geom_boxplot(outlier.size = 0) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "Method", y = expression("||diag("*hat(Sigma)[k] - Sigma[k]*")||"[2]), x = "n") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_diag_squared_error_n.pdf"), height = 5, width = 8.5)

