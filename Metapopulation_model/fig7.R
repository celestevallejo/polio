
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.args <- if (interactive()) c(
  file.path("data", "digest.rds"),
  file.path("Figures", "fig7.png")
) else commandArgs(trailingOnly = TRUE)

plot.dt <- readRDS(.args[1])

p <- ggplot(plot.dt[
  (label %in% c("16Kx4","4Kx16")) & movModel == 0 & ES_Detection == 0 & vacRate == 0
][len >= 0]) + facet_grid(label ~ .) +
  aes(
    len/365, 1-pext,
    linetype = label == "64Kx1", color = moveRate,
    group = interaction(label, moveRate)
  ) + # facet_grid(moveRate ~ vacRate) +
  geom_line(data = function(dt) dt[label != "64Kx1"]) +
  geom_line(aes(linetype = TRUE, group = NULL), data = plot.dt[
    label == "64Kx1" & moveRate == 0 & movModel == 0 & ES_Detection == 0 & vacRate == 0
  ][len >= 0, .SD, .SDcols = -c("label")]) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(
    "Time Since Case Observed (years)"
  ) +
  coord_cartesian(ylim=c(0.001, 0.999), xlim = c(0, 3)) +
  scale_y_continuous(
    "P(silent circulation)", #trans = "logit",
    breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.90, 1)
  ) +
  scale_linetype_discrete(NULL, breaks = TRUE, labels = "64Kx1") +
  scale_color_continuous(expression(alpha), guide = "legend", breaks = plot.dt[, sort(unique(moveRate))]) +
  theme(
    legend.position = c(0,0), legend.justification = c(0,0),
    legend.spacing = unit(0, "line"), legend.margin = margin(t=-0.5, unit="line")
  )

pinsetylim = plot.dt[
  between(len/365, 3, 4) & (label %in% c("16Kx4", "4Kx16")) & movModel == 0 & ES_Detection == 0 & vacRate == 0
][len >= 0][, ceiling(max(1-pext)*100)/100 ]

pinset16 <- ggplot(plot.dt[
  between(len/365, 3, 4) & (label %in% c("16Kx4")) & movModel == 0 & ES_Detection == 0 & vacRate == 0
][len >= 0]) + #facet_grid(label ~ .) +
  aes(
    len/365, 1-pext,
    linetype = label == "64Kx1", color = moveRate,
    group = interaction(label, moveRate)
  ) + # facet_grid(moveRate ~ vacRate) +
  geom_line(data = function(dt) dt[label != "64Kx1"]) +
  geom_line(aes(linetype = TRUE, group = NULL), data = plot.dt[
    between(len/365, 3, 4) & (label == "64Kx1") & moveRate == 0 & movModel == 0 & ES_Detection == 0 & vacRate == 0
  ][len >= 0, .SD, .SDcols = -c("label")]) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(NULL, breaks = seq(3,4,by=0.5)) +
  coord_cartesian(xlim = c(3, 4), ylim = c(0, pinsetylim), expand = FALSE) +
  scale_y_continuous(NULL) +
  scale_linetype_discrete(NULL, guide = "none") +
  scale_color_continuous(NULL, guide = "none") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.1), plot.margin = margin(r=1, t=1, unit = "line"),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )

pinset4 <- ggplot(plot.dt[
  between(len/365, 3, 4) & (label %in% c("4Kx16")) & movModel == 0 & ES_Detection == 0 & vacRate == 0
][len >= 0]) + #facet_grid(label ~ .) +
  aes(
    len/365, 1-pext,
    linetype = label == "64Kx1", color = moveRate,
    group = interaction(label, moveRate)
  ) + # facet_grid(moveRate ~ vacRate) +
  geom_line(data = function(dt) dt[label != "64Kx1"]) +
  geom_line(aes(linetype = TRUE, group = NULL), data = plot.dt[
    between(len/365, 3, 4) & (label == "64Kx1") & moveRate == 0 & movModel == 0 & ES_Detection == 0 & vacRate == 0
  ][len >= 0, .SD, .SDcols = -c("label")]) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(NULL, breaks = seq(3,4,by=0.5)) +
  coord_cartesian(xlim = c(3, 4), ylim = c(0, pinsetylim), expand = FALSE) +
  scale_y_continuous(NULL) +
  scale_linetype_discrete(NULL, guide = "none") +
  scale_color_continuous(NULL, guide = "none") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.1), plot.margin = margin(r=1, t=1, unit = "line"),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )

resp <- p + inset_element(pinset16, right = 1, top = 1, left = 0.6, bottom = 0.75) + inset_element(pinset4, right = 1, top = 0.5, left = 0.6, bottom = 0.25)

ggsave(tail(.args, 1), resp, width = 7, height = 5, dpi = 600, bg = "white")
