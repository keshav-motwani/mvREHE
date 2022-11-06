library(tidyverse)

RESULT_PATH = "results/"
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]

results = lapply(files, readRDS)

error_df = do.call(rbind, lapply(results, function(x) x$spectral_error))
error_df$error = error_df$spectral_error
time_df = do.call(rbind, lapply(results, function(x) x$time))
truncated = do.call(rbind, lapply(results, function(x) x$truncated))
min_eigenvalue_df = do.call(rbind, lapply(results, function(x) x$min_eigenvalue))

truncated_proportion = truncated %>% 
  group_by(replicate, n, q, experiment) %>%
  summarize(truncated = max(truncated)) %>%
  group_by(n, q, experiment) %>%
  summarize(proportion = mean(truncated), count = n()) %>%
  mutate(method = "mvHE")
truncated_proportion1 = truncated_proportion
truncated_proportion = do.call(rbind, lapply(unique(truncated$estimate), function(x) truncated_proportion %>% mutate(estimate = x)))

error_df = left_join(error_df, truncated %>% mutate(method = NULL) %>% group_by(n, q, replicate, experiment) %>% summarize(truncated = max(truncated)), by = c("n", "q", "replicate", "experiment")) %>%
  mutate(truncated = ifelse(truncated == 0, "not truncated", "truncated"))

time_df = left_join(time_df, truncated %>% mutate(method = NULL) %>% group_by(n, q, replicate, experiment) %>% summarize(truncated = max(truncated)), by = c("n", "q", "replicate", "experiment")) %>%
  mutate(truncated = ifelse(truncated == 0, "not truncated", "truncated"))

palette = c("mvHE" = "#006E7F", "mvREHE" = "#EE5007", "mvREHE_init" = "#EE5007", "mvREHE_optim" = "#EE5007", "mvREML" = "#4DBBD5B2", "naive" = "gray")

ggplot(error_df %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = error, color = method)) +
  facet_wrap(~facet, scales = "free_y", ncol = 2) +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = truncated_proportion %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  ylab("Estimation error") +
  labs(color = "Method") +
  scale_color_manual(values = palette) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_error_n.pdf"), height = 4, width = 7)

ggplot(error_df %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = error, color = method)) +
  facet_grid(facet ~ truncated, scales = "free_y") +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = error_df %>% filter(method == "mvHE", experiment == "n") %>% group_by(n, estimate, q, truncated, method) %>% summarize(count = n()) %>% group_by(n, estimate, q, method) %>% mutate(proportion = count / sum(count)) %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  ylab("Estimation error") +
  labs(color = "Method") +
  scale_color_manual(values = palette) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_error_truncated_n.pdf"), height = 8, width = 7)

ggplot(min_eigenvalue_df %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = min_eigenvalue, color = method)) +
  facet_wrap(~facet, scales = "free_y", ncol = 2) +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = truncated_proportion %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q, "; ", estimate)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  ylab("Minimum eigenvalue") +
  labs(color = "Method") +
  scale_color_manual(values = palette["mvHE"]) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_min_eigenvalue_n.pdf"), height = 4, width = 7)

ggplot(time_df %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = log10(time), color = gsub(".elapsed", "", method))) +
  facet_wrap(~facet, scales = "free_y", ncol = 2, dir = "v") +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = truncated_proportion1 %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method") +
  scale_color_manual(values = palette[1:4]) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_time_n.pdf"), height = 2.5, width = 7)

ggplot(time_df %>% filter(experiment == "n") %>% mutate(facet = paste0("q = ", q)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(x = as.factor(n), y = log10(time), color = gsub(".elapsed", "", method))) +
  facet_grid(facet ~ truncated, scales = "free_y") +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = time_df  %>% mutate(method = gsub(".elapsed", "", method)) %>% filter(method == "mvHE", experiment == "n") %>% group_by(n, q, truncated, method) %>% summarize(count = n()) %>% group_by(n, q, method) %>% mutate(proportion = count / sum(count)) %>% mutate(facet = paste0("q = ", q)) %>% mutate(facet = factor(facet, levels = gtools::mixedsort(unique(facet)))), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method") +
  scale_color_manual(values = palette[1:4]) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_time_truncated_n.pdf"), height = 4, width = 7)


ggplot(error_df %>% filter(experiment == "q"), aes(x = as.factor(q), y = error, color = method)) +
  facet_wrap(~paste0("n = ", n, "; ", estimate), scales = "free_y", ncol = 2) +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = truncated_proportion %>% filter(experiment == "q"), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("q") +
  ylab("Estimation error") +
  labs(color = "Method") +
  scale_color_manual(values = palette) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_error_q.pdf"), height = 4, width = 7)

ggplot(error_df %>% filter(experiment == "q"), aes(x = as.factor(q), y = error, color = method)) +
  facet_grid(paste0("n = ", n, "; ", estimate) ~ truncated, scales = "free_y") +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = error_df %>% filter(method == "mvHE", experiment == "q") %>% group_by(n, estimate, q, truncated, method) %>% summarize(count = n()) %>% group_by(n, estimate, q, method) %>% mutate(proportion = count / sum(count)), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("q") +
  ylab("Estimation error") +
  labs(color = "Method") +
  scale_color_manual(values = palette) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_error_truncated_q.pdf"), height = 8, width = 7)

ggplot(min_eigenvalue_df %>% filter(experiment == "q"), aes(x = as.factor(q), y = min_eigenvalue, color = method)) +
  facet_wrap(~paste0("n = ", n, "; ", estimate), scales = "free_y", ncol = 2) +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = truncated_proportion %>% filter(experiment == "q"), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("q") +
  ylab("Minimum eigenvalue") +
  labs(color = "Method") +
  scale_color_manual(values = palette["mvHE"]) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_min_eigenvalue_q.pdf"), height = 4, width = 7)

ggplot(time_df %>% filter(experiment == "q"), aes(x = as.factor(q), y = log10(time), color = gsub(".elapsed", "", method))) +
  facet_wrap(~paste0("n = ", n), scales = "free_y", ncol = 2, dir = "v") +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = truncated_proportion1 %>% filter(experiment == "q"), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("q") +
  labs(color = "Method") +
  scale_color_manual(values = palette[1:4]) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_time_q.pdf"), height = 2.5, width = 7)

ggplot(time_df %>% filter(experiment == "q"), aes(x = as.factor(q), y = log10(time), color = gsub(".elapsed", "", method))) +
  facet_grid(paste0("n = ", n) ~ truncated, scales = "free_y") +
  geom_boxplot(outlier.size = 0.25) +
  geom_text(data = time_df %>% mutate(method = gsub(".elapsed", "", method)) %>% filter(method == "mvHE", experiment == "q") %>% group_by(n, q, truncated, method) %>% summarize(count = n()) %>% group_by(n, q, method) %>% mutate(proportion = count / sum(count)), aes(label = paste0(round(proportion * 100, 2), "%"), y = Inf, vjust = 2), size = 2.5) +
  theme_bw() +
  xlab("q") +
  labs(color = "Method") +
  scale_color_manual(values = palette[1:4]) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGURES_PATH, "simulation_figure_time_truncated_q.pdf"), height = 4, width = 7)

