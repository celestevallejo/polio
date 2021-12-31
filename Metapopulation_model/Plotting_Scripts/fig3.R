
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.args <- if (interactive()) c(
  file.path("data", "digest.rds"),
  file.path("Figures", "fig3.png")
) else commandArgs(trailingOnly = TRUE)

plot.dt <- readRDS(.args[1])
plot.dt[, icd := ic / sum(ic), by=cksum ]
plot.dt[len != -1, exd := ex / sum(ex), by=cksum ]

keys <- c("movModel", "moveRate", "ES_Detection", "vacRate", "vilModel", "cksum", "label")
p.dt <- melt(plot.dt, id.vars = c(keys, "len"), measure.vars = c("icd", "exd"))
refscen <- expression(movModel == 0 & moveRate == 0 & ES_Detection == 0 & vacRate == 0)

p <- ggplot(p.dt[eval(refscen)][len>0]) +
  facet_grid(
    variable ~ .,
    labeller = labeller(
      variable = c(icd="Intercase Intervals", exd="Extinction Intervals")
    )
#    , scales = "free_y"
  ) +
  aes(
    len/365, value,
    color = label
  ) + # facet_grid(moveRate ~ vacRate) +
  geom_line(data = function(dt) dt[variable == "icd"] ) +
  geom_line(data = function(dt) dt[variable != "icd"], alpha = 0.1) +
  geom_smooth(data = function(dt) dt[variable != "icd"]) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(
    "Time Since Case Observed (years)"
  ) +
  coord_cartesian(ylim=c(0, 0.005)) +
  scale_y_continuous(
    "Density" #trans = "logit",
#    , breaks = seq(-.2, 1, by=0.2)
#    , breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.90, 1)
  ) +
#  scale_linetype_discrete("Patch Arrangement") +
#  scale_color_continuous("Movement Rate", guide = "legend", breaks = plot.dt[, sort(unique(moveRate))]) +
  theme(legend.position = "bottom")

ggsave(tail(.args, 1), p, width = 7, height = 10, dpi = 600, bg = "white")
