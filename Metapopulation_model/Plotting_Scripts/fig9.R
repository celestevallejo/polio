
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.args <- if (interactive()) c(
  file.path("data", "digest.rds"),
  file.path("Figures", "fig9_alt.png")
) else commandArgs(trailingOnly = TRUE)

plot.dt <- readRDS(.args[1])
keys <- c("movModel", "moveRate", "ES_Detection", "vacRate", "vilModel", "cksum", "label")

ref.dt <- plot.dt[
  label == "1x64K",
  .(len = 0:max(len), pext = approx(len, pext, 0:max(len))$y),
  keyby=keys
]

alt.dt <- plot.dt[
  label != "1x64K" & !(label %like% "&"),
  .(len = 0:max(len), pext = approx(len, pext, 0:max(len))$y),
  keyby=keys
]

del.dt <- alt.dt[
  ref.dt[,.SD,.SDcols = -c("cksum", "label", "vilModel")],
  on=c(setdiff(keys, c("cksum","label","vilModel")), "len"),
  roll = Inf, rollends = TRUE
]

del.dt[, del.p := pext - i.pext ]
del.dt[, OF.p := 1/(i.pext/(1-i.pext) * (1-pext)/pext) ]

refscen <- expression(movModel == 1 & ES_Detection == 0 & vacRate == 0)

p <- ggplot(del.dt[eval(refscen)]) +
  aes(
    len/365, del.p,
    color = label
  ) + # facet_grid(moveRate ~ vacRate) +
  geom_line(data = function(dt) dt[moveRate == 0], alpha = 0.2) +
  geom_line(data = function(dt) dt[moveRate == 0.1]) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(
    "Time Since Case Observed (years)"
  ) +
#  coord_cartesian(ylim=c(-1, 1)) +
  scale_y_continuous(
    expression(Delta*"P(silent circulation)") #trans = "logit",
    , breaks = seq(-.2, 1, by=0.2)
#    , breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.90, 1)
  ) +
  scale_color_discrete(NULL) +
#  scale_linetype_discrete("Patch Arrangement") +
#  scale_color_continuous("Movement Rate", guide = "legend", breaks = plot.dt[, sort(unique(moveRate))]) +
  theme(legend.position = "bottom")

pinset <- ggplot(del.dt[eval(refscen)][between(len/365, 2.25, 3.75)]) +
  aes(
    len/365, del.p,
    color = label
  ) + 
  geom_line(data = function(dt) dt[moveRate == 0], alpha = 0.2) +
  geom_line(data = function(dt) dt[moveRate == 0.1]) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(NULL, breaks = seq(2.5, 3.5, by=.5)) +
  coord_cartesian(xlim=c(2.5, 3.5), ylim = c(-5,5)/100) +
  scale_y_continuous(NULL) +
  scale_color_discrete(guide = "none") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.1), plot.margin = margin(r=1, t=1, unit = "line"),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )

resp <- p + inset_element(pinset, left = 0.6, right = 1, top = 1, bottom = 0.5)

ggsave(tail(.args, 1), resp, width = 7, height = 5, dpi = 600, bg = "white")
