library(tidyverse)

SIMULATION_ID = 3
RESULT_PATH = paste0("simulation_hcp_results_", SIMULATION_ID)
FIGURES_PATH = file.path(RESULT_PATH, "figures")
dir.create(FIGURES_PATH, recursive = TRUE)

methods = c("mvHE", "mvREHE", "mvHE-smoothed", "mvREHE-smoothed")
palette = ggsci::pal_aaas("default")(6)
palette = palette[c(1, 2, 5, 6)]
names(palette) = methods
options(ggplot2.discrete.colour = palette)

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]

results = lapply(files, readRDS)

Sigma_labels = c("alpha == 1", "alpha == 2")
names(Sigma_labels) = c("smooth_1", "smooth_2")

### Time

time_df = do.call(rbind, lapply(results, function(x) x$time)) %>%
  group_by(n, q, Sigma, experiment, method) %>%
  summarize(mean = mean(time), se = sd(time) / sqrt(n())) %>%
  filter(method %in% methods)

ggplot(time_df %>%
         filter(experiment == "n") %>%
         mutate(facet = Sigma_labels[Sigma]) %>%
         mutate(facet = factor(facet, levels = Sigma_labels)),
       aes(x = n, y = log10(mean), ymin = log10(mean - 1.96 * se), ymax = log10(mean + 1.96 * se), color = factor(gsub(".elapsed", "", method), levels = names(palette)), group = method)) +
  facet_wrap(~facet, scales = "free_y", nrow = 1, labeller = labeller(facet = label_parsed)) +
  geom_line() +
  geom_errorbar(width = 0.1) +
  theme_bw() +
  xlab("n") +
  labs(color = "Method", y = "log10(seconds)") +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.placement = "outside")
ggsave(file.path(FIGURES_PATH, paste0("simulation_figure_time_", "n", ".pdf")), height = 3, width = 6)

label = c("hat(Sigma)[G]", "hat(Sigma)[C]", "hat(Sigma)[E]")
names(label) = c("Sigma_1", "Sigma_2", "Sigma_0")

Sigmas = apply(expand.grid(label, Sigma_labels), 1, function(x) paste0(x[2], "~(", x[1], ")"))

### Squared error

squared_error_df = do.call(rbind, lapply(results, function(x) x$squared_error)) %>%
  mutate(estimate = label[estimate]) %>%
  filter(method %in% methods)

  ggplot(squared_error_df %>%
           filter(experiment == "n" & grepl("mv", method)) %>%
           mutate(facet = paste0(Sigma_labels[Sigma], "~(", estimate, ")")) %>%
           mutate(facet = factor(facet, levels = Sigmas)) %>%
           group_by(n, facet, method) %>%
           summarize(mean = mean(squared_error) / 1000^2, se = sd(squared_error) / sqrt(n()) / 1000^2),
         aes(x = n, y = mean, ymax = mean + 1.96 * se, ymin = mean - 1.96 * se,
             color = factor(method, levels = names(palette)))) +
    facet_wrap(~facet, scales = "free_y", ncol = 3, dir = "h", labeller = labeller(facet = label_parsed)) +
    geom_line() +
    geom_errorbar(width = 0.1) +
    theme_bw() +
    xlab("n") +
    labs(color = "Method", linetype = "Method", y = expression(integral(integral((hat(C)[k](s, t) - C[k](s, t)))^2)*ds*dt), x = "n") +
    theme(legend.position = "bottom") +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    scale_y_continuous(limits = c(0, NA))
  ggsave(file.path(FIGURES_PATH, paste0("simulation_figure_squared_error_", "n", ".pdf")), height = 7.5 * 0.8 * 2.2 / 3, width = 8.5)
