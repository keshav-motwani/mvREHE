library(mvREHE)
library(Matrix)

source("scripts/simulation_setup.R")

COMPONENTS = 3
n = 1000
q = NA
Sigma = "data_100"
SIMULATION_ID = "data"
replicate = 1
DATA_ANALYSIS_RESULT_PATH = paste0("data_analysis_", COMPONENTS, "_components")

output_mvREHE = simulation(COMPONENTS, n, q, Sigma, method = "mvREHE", SIMULATION_ID, replicate, DATA_ANALYSIS_RESULT_PATH)
output_REML = simulation(COMPONENTS, n, q, Sigma, method = "REML", SIMULATION_ID, replicate, DATA_ANALYSIS_RESULT_PATH)

h2 = function(Sigma_list) {
  sapply(1:ncol(Sigma_list[[1]]), function(j) Sigma_list[[2]][j, j] / sum(sapply(Sigma_list, function(Sigma) Sigma[j, j])))
}

result = data.frame(mvREHE = h2(output_mvREHE$estimate$Sigma_hat),
                    REML = h2(output_REML$estimate$Sigma_hat),
                    truth = h2(output_mvREHE$Sigma_list_truth),
                    type = rep(c("Functional", "Structural"), each = 91))

REML = ggplot(result, aes(x = truth, y = REML, color = type)) +
  geom_point() +
  theme_bw() +
  ylab(expression(hat(h)[j]^2~"from REML")) +
  xlab(expression("True"~h[j]^2)) +
  xlim(c(0, 0.4)) +
  ylim(c(0, 0.4)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(color = "Connection Type")

mvREHE = ggplot(result, aes(x = truth, y = mvREHE, color = type)) +
  geom_point() +
  theme_bw() +
  ylab(expression(hat(h)[j]^2~"from mvREHE")) +
  xlab(expression("True"~h[j]^2)) +
  xlim(c(0, 0.4)) +
  ylim(c(0, 0.4)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(color = "Connection Type")

both = ggplot(result, aes(x = REML, y = mvREHE, color = type)) +
  geom_point() +
  theme_bw() +
  ylab(expression(hat(h)[j]^2~"from mvREHE")) +
  xlab(expression(hat(h)[j]^2~"from REML")) +
  xlim(c(0, 0.4)) +
  ylim(c(0, 0.4)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(color = "Connection Type")

patchwork::wrap_plots(list(both, REML, mvREHE)) +
  patchwork::plot_layout(guides = 'collect') &
  theme(legend.position='bottom')

ggsave(file.path(DATA_ANALYSIS_RESULT_PATH, "simulated_h2.pdf"), height = 3.4, width = 8.5)
