require(data.table)
require(ggplot2)

.args <- if (interactive()) c(
  file.path("data", "state_digest.rds"),
  file.path("Figures", "state_digest.rds")
) else commandArgs(trailingOnly = TRUE)

plot.dt <- readRDS(.args[1])

aes(x = X1, ymin=Min, lower=`2.5%`, middle = `50%`, upper = `97.5%`, ymax = Max)) + 
  geom_boxplot(stat="identity")

p <- ggplot(plot.dt[movModel == 0 & moveRate == 1 & ES_Detection == 0]) +
  aes(
    factor(village_size),
    color = label
  ) + facet_grid(vacRate ~ variable) +
  geom_boxplot(
    aes(
      ymin =lo95/village_size,
      lower =lo50/village_size,
      middle = md/village_size,
      upper = hi50/village_size,
      ymax = hi95/village_size,
      group = interaction(cksum, village_size)
    ),
    position = "dodge", stat = "identity"
  ) +
  theme_minimal(base_size = 14) +
  scale_x_discrete(
    "Village Size",
    labels = function(f) gsub("000","K", as.character(f))
  ) +
  coord_cartesian(ylim=c(0.001, 0.6)) +
  scale_y_continuous(
    "% Compartment (logit scale)", trans = "logit",
    breaks = c(0.001, 0.01, 0.1, 0.25, 0.5)
  ) +
  scale_color_discrete(NULL) +
  theme(legend.position = "bottom")

ggsave(tail(.args, 1), p, width = 14, height = 7, dpi = 600, bg = "white")
