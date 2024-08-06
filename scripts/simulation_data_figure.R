library(tidyverse)

SIMULATION_ID = "data"
COMPONENTS = as.numeric(commandArgs(trailingOnly=TRUE)[1])
RESULT_PATH = paste0("simulation_", COMPONENTS, "_components_", SIMULATION_ID)
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

methods = c("mvHE", "mvREHE", "HE", "REHE", "REML")
map = gsub("mv", "", methods)
names(map) = methods

palette = ggsci::pal_aaas("default")(3)
names(palette) = c("HE", "REHE", "REML")
options(ggplot2.discrete.colour = palette)

Sigmas = "data_100" # c("data_1.01", "data_10", "data_100", "data_1000", "data_10000")
Sigma_labels = "data" # paste0("cond == ", unlist(lapply(strsplit(Sigmas, "_"), `[[`, 2)))
names(Sigma_labels) = Sigmas

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]

results = lapply(files, readRDS)

### Time

time_df = do.call(rbind, lapply(results, function(x) x$time)) %>%
  group_by(n, q, Sigma, experiment, method) %>%
  summarize(mean = mean(na.rm = TRUE, time), se = sd(na.rm = TRUE, time) / sqrt(n())) %>%
  filter(method %in% methods)

ggplot(time_df %>%
         filter(experiment == "n") %>%
         mutate(facet = Sigma_labels[Sigma]),
       aes(x = n, y = log10(mean), ymax = log10(mean + 1.96 * se), ymin = log10(mean - 1.96 * se),
           color = factor(map[method], levels = names(palette)),
           linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"))) +
  facet_wrap(~facet, scales = "fixed", nrow = 1) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "", y = "log10(seconds)") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside", plot.margin = unit(c(5 / 72.27, 2.5, 5 / 72.27, 2.5), "in"))
ggsave(file.path(FIGURES_PATH, "simulation_figure_time_n.pdf"), height = 3, width = 8.5)

### Squared error on diagonal

label = c("hat(Sigma)[G]", "hat(Sigma)[C]", "hat(Sigma)[E]")
names(label) = c("Sigma_1", "Sigma_2", "Sigma_0")

Sigmas = apply(expand.grid(label, Sigma_labels), 1, function(x) paste0(x[2], "~(", x[1], ")"))

diag_squared_error_df = do.call(rbind, lapply(results, function(x) x$diag_squared_error)) %>%
  mutate(estimate = label[as.character(estimate)]) %>%
  filter(method %in% methods)

ggplot(diag_squared_error_df %>%
         filter(experiment == "n") %>%
         mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = Sigmas)) %>%
         group_by(n, facet, method) %>%
         summarize(mean = mean(na.rm = TRUE, diag_squared_error), se = sd(na.rm = TRUE, diag_squared_error) / sqrt(n())),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(map[method], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"), group = method)) +
  facet_wrap(~facet, scales = "free_y", ncol = COMPONENTS, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", y = expression("E||diag("*hat(Sigma)[k] - Sigma[k]*")||"[2]), x = "n") +   theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(NA, NA))
ggsave(file.path(FIGURES_PATH, "simulation_figure_diag_squared_error_n.pdf"), height = 7.5 * 0.8, width = 8.5)

### h2 error

h2_df = do.call(rbind, lapply(results, function(x) x$h2_error)) %>%
  group_by(n, q, Sigma, experiment, method) %>%
  summarize(mean = mean(na.rm = TRUE, h2_error), se = sd(na.rm = TRUE, h2_error) / sqrt(n())) %>%
  filter(method %in% methods)

ggplot(h2_df %>%
         filter(experiment == "n") %>%
         mutate(facet = Sigma_labels[Sigma]),
       aes(x = n, y = sqrt(mean), ymin = sqrt(mean - 1.96 * se), ymax = sqrt(mean + 1.96 * se),
           color = factor(map[method], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"), group = method)) +
  facet_wrap(~facet, scales = "free_y", nrow = 1) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "", y = expression("E||"*hat(h)^2 - h^2*"||"[2])) +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside", plot.margin = unit(c(5 / 72.27, 2.5, 5 / 72.27, 2.5), "in"))
ggsave(file.path(FIGURES_PATH, "simulation_figure_h2_error_n.pdf"), height = 3, width = 8.5)

### Spectral error

names(palette) = paste0("mv", names(palette))
options(ggplot2.discrete.colour = palette)

spectral_error_df = do.call(rbind, lapply(results, function(x) x$spectral_error)) %>%
  mutate(estimate = label[as.character(estimate)]) %>%
  filter(method %in% methods)

ggplot(spectral_error_df %>%
         filter(experiment == "n" & grepl("mv", method)) %>%
         mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = Sigmas)) %>%
         group_by(n, facet, method) %>%
         summarize(mean = mean(na.rm = TRUE, spectral_error), se = sd(na.rm = TRUE, spectral_error) / sqrt(n())),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(method, levels = names(palette)))) +
  facet_wrap(~facet, scales = "free_y", ncol = COMPONENTS, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", y = expression("E||"*hat(Sigma)[k] - Sigma[k]*"||"[2]), x = "n") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(NA, NA))
ggsave(file.path(FIGURES_PATH, "simulation_figure_spectral_error_n.pdf"), height = 3, width = 8.5)

### Squared error

squared_error_df = do.call(rbind, lapply(results, function(x) x$squared_error)) %>%
  mutate(estimate = label[as.character(estimate)]) %>%
  filter(method %in% methods)

ggplot(squared_error_df %>%
         filter(experiment == "n" & grepl("mv", method)) %>%
         mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = Sigmas)) %>%
         group_by(n, facet, method) %>%
         summarize(mean = mean(na.rm = TRUE, squared_error), se = sd(na.rm = TRUE, squared_error) / sqrt(n())),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(method, levels = names(palette)))) +
  facet_wrap(~facet, scales = "free_y", ncol = COMPONENTS, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", y = expression("E||"*hat(Sigma)[k] - Sigma[k]*"||"[F]), x = "n") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(NA, NA))
ggsave(file.path(FIGURES_PATH, "simulation_figure_squared_error_n.pdf"), height = 7.5 * 0.8, width = 8.5)


### Regression coefficients error

beta_error_df = do.call(rbind, lapply(results, function(x) x$beta_error)) %>%
  mutate(estimate = label[as.character(estimate)]) %>%
  filter(method %in% methods)

ggplot(beta_error_df %>%
         filter(experiment == "n" & grepl("mv", method)) %>%
         mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = Sigmas)) %>%
         group_by(n, facet, method) %>%
         summarize(mean = mean(na.rm = TRUE, beta_error), se = sd(na.rm = TRUE, beta_error) / sqrt(n())),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(map[method], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"), group = method)) +
  facet_wrap(~facet, scales = "free_y", ncol = COMPONENTS, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "", y = expression("E||"*hat(beta)[k] - beta[k]*"||"[2]), x = "n") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(NA, NA))
ggsave(file.path(FIGURES_PATH, "simulation_figure_beta_error_n.pdf"), height = 3, width = 8.5)

### Max principal angle

max_principal_angle_df = do.call(rbind, lapply(results, function(x) x$max_principal_angle)) %>%
  mutate(estimate = label[as.character(estimate)]) %>%
  filter(method %in% methods)

ggplot(max_principal_angle_df %>%
         filter(experiment == "n" & grepl("mv", method) & r == 1) %>%
         mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = Sigmas)) %>%
         group_by(n, facet, method) %>%
         summarize(mean = mean(na.rm = TRUE, value), se = sd(na.rm = TRUE, value) / sqrt(n())),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(map[method], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"), group = method)) +
  facet_wrap(~facet, scales = "free_y", ncol = COMPONENTS, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "", y = "Angle between true and estimated\nloading vector of PC 1", x = "n") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(NA, NA))
ggsave(file.path(FIGURES_PATH, "simulation_figure_max_principal_angle_1_n.pdf"), height = 3, width = 8.5)

ggplot(max_principal_angle_df %>%
         filter(experiment == "n" & grepl("mv", method) & r == 3) %>%
         mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = Sigmas)) %>%
         group_by(n, facet, method) %>%
         summarize(mean = mean(na.rm = TRUE, value), se = sd(na.rm = TRUE, value) / sqrt(n())),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(method, levels = names(palette)))) +
  facet_wrap(~facet, scales = "free_y", ncol = COMPONENTS, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", y = "Principal angle - 3 PCs", x = "n") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(limits = c(NA, NA))
ggsave(file.path(FIGURES_PATH, "simulation_figure_max_principal_angle_3_n.pdf"), height = 7.5 * 0.8, width = 8.5)

