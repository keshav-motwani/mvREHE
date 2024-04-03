library(tidyverse)

SIMULATION_ID = 1
RESULT_PATH = paste0("simulation_hcp_results_", SIMULATION_ID)
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

methods = c("mvHE", "mvREHE", "mvREML", "HE", "REHE", "REML")
palette = ggsci::pal_npg("nrc")(3)
names(palette) = c("HE", "REHE", "REML")
options(ggplot2.discrete.colour = palette)
map = gsub("mv", "", methods)
names(map) = methods

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]

results = lapply(files, readRDS)

### Time

time_df = do.call(rbind, lapply(results, function(x) x$time)) %>%
  group_by(n, q, experiment, method) %>%
  summarize(mean = mean(time), se = sd(time) / sqrt(n())) %>%
  filter(method %in% methods)

ggplot(time_df %>%
         filter(experiment == "n") %>%
         mutate(facet = paste0("q = ", q)) %>%
         mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))),
       aes(x = n, y = log10(mean), ymin = log10(mean - 1.96 * se), ymax = log10(mean + 1.96 * se), color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"), group = method)) +
  facet_wrap(~facet, scales = "free_y", nrow = 1) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "Method", y = "log10(seconds)") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside")
ggsave(file.path(FIGURES_PATH, "simulation_figure_time_n.pdf"), height = 3, width = 8.5)

### Spectral error

label = c("hat(Sigma)[G]", "hat(Sigma)[C]", "hat(Sigma)[E]")
names(label) = c("Sigma_1", "Sigma_2", "Sigma_0")

facets = apply(expand.grid(label, paste0("q == ", c(5, 20, 100))), 1, function(x) paste0(x[2], "~(", x[1], ")"))

spectral_error_df = do.call(rbind, lapply(results, function(x) x$spectral_error)) %>%
  mutate(estimate = label[estimate]) %>%
  filter(method %in% methods)

ggplot(spectral_error_df %>%
         filter(experiment == "n" & grepl("mv", method)) %>%
         mutate(facet = paste0("q == ", q, "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = facets)) %>%
         group_by(n, facet, method) %>%
         summarize(mean = mean(spectral_error), se = sd(spectral_error) / sqrt(n())),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)))) +
  facet_wrap(~facet, scales = "free_y", ncol = 3, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "Method", y = expression("E||"*hat(Sigma)[k] - Sigma[k]*"||"[2]), x = "n") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
ggsave(file.path(FIGURES_PATH, "simulation_figure_spectral_error_n.pdf"), height = 7.5 * 0.8, width = 8.5)

### Squared error

squared_error_df = do.call(rbind, lapply(results, function(x) x$squared_error)) %>%
  mutate(estimate = label[estimate]) %>%
  filter(method %in% methods)

ggplot(squared_error_df %>%
         filter(experiment == "n" & grepl("mv", method)) %>%
         mutate(facet = paste0("q == ", q, "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = facets)) %>%
         group_by(n, facet, method) %>%
         summarize(mean = mean(squared_error), se = sd(squared_error) / sqrt(n())),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)))) +
  facet_wrap(~facet, scales = "free_y", ncol = 3, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "Method", y = expression("E||"*hat(Sigma)[k] - Sigma[k]*"||"[F]), x = "n") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
ggsave(file.path(FIGURES_PATH, "simulation_figure_squared_error_n.pdf"), height = 7.5 * 0.8, width = 8.5)

### Squared error on diagonal

diag_squared_error_df = do.call(rbind, lapply(results, function(x) x$diag_squared_error)) %>%
  mutate(estimate = label[estimate]) %>%
  filter(method %in% methods)

ggplot(diag_squared_error_df %>%
         filter(experiment == "n" & grepl("mv", method)) %>%
         mutate(facet = paste0("q == ", q, "~(", estimate, ")")) %>%
         mutate(facet = factor(facet, levels = facets)) %>%
         group_by(n, facet, method) %>%
         summarize(mean = mean(diag_squared_error), se = sd(diag_squared_error) / sqrt(n())),
       aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
           color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)))) +
  facet_wrap(~facet, scales = "free_y", ncol = 3, dir = "h", labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "Method", y = expression("||diag("*hat(Sigma)[k] - Sigma[k]*")||"[2]), x = "n") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
ggsave(file.path(FIGURES_PATH, "simulation_figure_diag_squared_error_n.pdf"), height = 7.5 * 0.8, width = 8.5)


### h2 error

h2_df = do.call(rbind, lapply(results, function(x) x$h2_error)) %>%
  group_by(n, q, experiment, method) %>%
  summarize(mean = mean(h2_error), se = sd(h2_error) / sqrt(n())) %>%
  filter(method %in% methods)

ggplot(h2_df %>%
         filter(experiment == "n") %>%
         mutate(facet = paste0("q = ", q)) %>%
         mutate(facet = factor(facet, levels = facets)),
       aes(x = n, y = mean, ymin = mean - 1.96 * se, ymax = mean + 1.96 * se,
           color = factor(map[gsub(".elapsed", "", method)], levels = names(palette)), linetype = ifelse(grepl("mv", method), "Multivariate", "Univariate"), group = method)) +
  facet_wrap(~facet, scales = "free_y", nrow = 1) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", linetype = "Method", y = expression(E(hat(h)^2 - h^2)^2)) +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside")
ggsave(file.path(FIGURES_PATH, "simulation_figure_h2_error_n.pdf"), height = 3, width = 8.5)
