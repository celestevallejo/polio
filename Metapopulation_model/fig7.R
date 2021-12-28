
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.args <- if (interactive()) c(
  file.path("data", "digest.rds"),
  file.path("Figures", "fig7.png")
)

plot.dt <- readRDS(.args[1])

p <- ggplot(plot.dt[
  (label %in% c("64Kx1","16Kx4")) & movModel == 0 & ES_Detection == 0 & vacRate == 0
][len >= 0]) +
  aes(
    len/365, 1-pext,
    linetype = label, color = moveRate,
    group = interaction(label, moveRate)
  ) + # facet_grid(moveRate ~ vacRate) +
  geom_line(data = function(dt) dt[label != "64Kx1"]) +
  geom_line(data = function(dt) dt[label == "64Kx1" & moveRate == 0]) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(
    "Time Since Case Observed (years)"
  ) +
  coord_cartesian(ylim=c(0.001, 0.999), xlim = c(0, 3)) +
  scale_y_continuous(
    "P(silent circulation)", #trans = "logit",
    breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.90, 1)
  ) +
  scale_linetype_discrete("Patch Arrangement") +
  scale_color_continuous("Movement Rate", guide = "legend", breaks = plot.dt[, sort(unique(moveRate))]) +
  theme(legend.position = "bottom")

ggsave(tail(.args, 1), p, width = 7, height = 5, dpi = 600, bg = "white")
