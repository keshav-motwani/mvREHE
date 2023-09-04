library(tidyverse)

RESULT_PATH = "simulation_results/"
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]
files = files[grepl("_rInf_", files)]

results = lapply(files, readRDS)

spectral_error_df = do.call(rbind, lapply(results, function(x) x$spectral_error))
time_df = do.call(rbind, lapply(results, function(x) x$time))
truncated = do.call(rbind, lapply(results, function(x) x$truncated))
min_eigenvalue_df = do.call(rbind, lapply(results, function(x) x$min_eigenvalue))

truncated_proportion = truncated %>%
  filter(method == "mvHE") %>%
  group_by(replicate, n, q, experiment, method) %>%
  summarize(truncated = max(truncated)) %>%
  group_by(n, q, experiment, method) %>%
  summarize(proportion = mean(truncated), count = n())
truncated_proportion1 = truncated_proportion
truncated_proportion = do.call(rbind, lapply(unique(truncated$estimate), function(x) truncated_proportion %>% mutate(estimate = x)))

spectral_error_df = left_join(spectral_error_df, truncated %>% mutate(method = NULL) %>% group_by(n, q, replicate, experiment) %>% summarize(truncated = max(truncated)), by = c("n", "q", "replicate", "experiment")) %>%
  mutate(truncated = ifelse(truncated == 0, "not truncated", "truncated"))

time_df = left_join(time_df, truncated %>% mutate(method = NULL) %>% group_by(n, q, replicate, experiment) %>% summarize(truncated = max(truncated)), by = c("n", "q", "replicate", "experiment")) %>%
  mutate(truncated = ifelse(truncated == 0, "not truncated", "truncated"))

methods = c("mvHE", "mvREHE", "mvREHE_L2", "cv_mvREHE_L2", "mvREML", "naive")

spectral_error_df$method = factor(spectral_error_df$method, levels = methods)
time_df$method = factor(time_df$method, levels = methods)
truncated$method = factor(truncated$method, levels = methods)
min_eigenvalue_df$method = factor(min_eigenvalue_df$method, levels = methods)

saveRDS(list(spectral_error = spectral_error_df, time = time_df), file.path(FIGURES_PATH, "data.rds"))

palette = ggsci::pal_npg("nrc")(length(methods))
names(palette) = methods
options(ggplot2.discrete.colour= palette)

### Spectral error

ggplot(spectral_error_df %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = spectral_error, color = method)) +
  facet_wrap(~facet, scales = "free_y", ncol = 2) +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = truncated_proportion %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  ylab("Spectral error") +
  labs(color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_spectral_error_n.pdf"), height = 5.5, width = 7)

ggplot(spectral_error_df %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = spectral_error, color = method)) +
  facet_grid(facet ~ truncated, scales = "free_y") +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = spectral_error_df %>% filter(method == "mvHE", experiment == "n") %>% group_by(n, estimate, q, truncated, method) %>% summarize(count = n()) %>% group_by(n, estimate, q, method) %>% mutate(proportion = count / sum(count)) %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  ylab("Spectral error") +
  labs(color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_spectral_error_truncated_n.pdf"), height = 8, width = 7)

### Min eigenvalue

ggplot(min_eigenvalue_df %>% filter(experiment == "n", method == "mvHE") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = min_eigenvalue, color = method)) +
  facet_wrap(~facet, scales = "free_y", ncol = 2) +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = truncated_proportion %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  ylab("Minimum eigenvalue") +
  labs(color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_min_eigenvalue_n.pdf"), height = 4, width = 7)

### Time

ggplot(time_df %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = log10(time), color = gsub(".elapsed", "", method))) +
  facet_wrap(~facet, scales = "free_y", ncol = 1, dir = "v") +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = truncated_proportion1 %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_time_n.pdf"), height = 5.5, width = 5)

ggplot(time_df %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = log10(time), color = gsub(".elapsed", "", method))) +
  facet_grid(facet ~ truncated, scales = "free_y") +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = time_df  %>% mutate(method = gsub(".elapsed", "", method)) %>% filter(method == "mvHE", experiment == "n") %>% group_by(n, q, truncated, method) %>% summarize(count = n()) %>% group_by(n, q, method) %>% mutate(proportion = count / sum(count)) %>% mutate(facet = paste0("q = ", q)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_time_truncated_n.pdf"), height = 4, width = 7)


# ### Spectral error
#
# ggplot(spectral_error_df %>% filter(experiment == "q"), aes(x = as.factor(q), y = spectral_error, color = method)) +
#   facet_wrap(~paste0("n = ", n, "; ", estimate), scales = "free_y", ncol = 2) +
#   geom_boxplot(outlier.size = 0.25) +
#   geom_text(data = truncated_proportion %>% filter(experiment == "q"), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
#   theme_bw() +
#   xlab("q") +
#   ylab("Spectral error") +
#   labs(color = "Method") +
#   theme(legend.position = "bottom")
# ggsave(file.path(FIGURES_PATH, "simulation_figure_spectral_error_q.pdf"), height = 4, width = 7)
#
# ggplot(spectral_error_df %>% filter(experiment == "q"), aes(x = as.factor(q), y = spectral_error, color = method)) +
#   facet_grid(paste0("n = ", n, "; ", estimate) ~ truncated, scales = "free_y") +
#   geom_boxplot(outlier.size = 0.25) +
#   geom_text(data = spectral_error_df %>% filter(method == "mvHE", experiment == "q") %>% group_by(n, estimate, q, truncated, method) %>% summarize(count = n()) %>% group_by(n, estimate, q, method) %>% mutate(proportion = count / sum(count)), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
#   theme_bw() +
#   xlab("q") +
#   ylab("Spectral error") +
#   labs(color = "Method") +
#   theme(legend.position = "bottom")
# ggsave(file.path(FIGURES_PATH, "simulation_figure_spectral_error_truncated_q.pdf"), height = 8, width = 7)
#
# ### Min eigenvalue
#
# ggplot(min_eigenvalue_df %>% filter(experiment == "q", method == "mvHE"), aes(x = as.factor(q), y = min_eigenvalue, color = method)) +
#   facet_wrap(~paste0("n = ", n, "; ", estimate), scales = "free_y", ncol = 2) +
#   geom_boxplot(outlier.size = 0.25) +
#   geom_text(data = truncated_proportion %>% filter(experiment == "q"), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
#   theme_bw() +
#   xlab("q") +
#   ylab("Minimum eigenvalue") +
#   labs(color = "Method") +
#   theme(legend.position = "bottom")
# ggsave(file.path(FIGURES_PATH, "simulation_figure_min_eigenvalue_q.pdf"), height = 4, width = 7)
#
# ### Time
#
# ggplot(time_df %>% filter(experiment == "q"), aes(x = as.factor(q), y = log10(time), color = gsub(".elapsed", "", method))) +
#   facet_wrap(~paste0("n = ", n), scales = "free_y", ncol = 1, dir = "v") +
#   geom_boxplot(outlier.size = 0.25) +
#   geom_text(data = truncated_proportion1 %>% filter(experiment == "q"), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
#   theme_bw() +
#   xlab("q") +
#   labs(color = "Method") +
#   theme(legend.position = "bottom")
# ggsave(file.path(FIGURES_PATH, "simulation_figure_time_q.pdf"), height = 4, width = 5)
#
# ggplot(time_df %>% filter(experiment == "q"), aes(x = as.factor(q), y = log10(time), color = gsub(".elapsed", "", method))) +
#   facet_grid(paste0("n = ", n) ~ truncated, scales = "free_y") +
#   geom_boxplot(outlier.size = 0.25) +
#   geom_text(data = time_df %>% mutate(method = gsub(".elapsed", "", method)) %>% filter(method == "mvHE", experiment == "q") %>% group_by(n, q, truncated, method) %>% summarize(count = n()) %>% group_by(n, q, method) %>% mutate(proportion = count / sum(count)), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
#   theme_bw() +
#   xlab("q") +
#   labs(color = "Method") +
#   theme(legend.position = "bottom")
# ggsave(file.path(FIGURES_PATH, "simulation_figure_time_truncated_q.pdf"), height = 4, width = 7)
#
#
