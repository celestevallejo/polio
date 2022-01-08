
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
  require(RColorBrewer)
})

.args <- if (interactive()) c(
  file.path("data", "digest.rds"),
  file.path("Figures", "fig3.png")
) else commandArgs(trailingOnly = TRUE)

plot.dt <- readRDS(.args[1])
plot.dt[, icd := cumsum(ic) / sum(ic), by=cksum ]
plot.dt[len != -1, exd := cumsum(ex) / sum(ex), by=cksum ]

keys <- c("movModel", "moveRate", "ES_Detection", "vacRate", "vilModel", "cksum", "label")
p.dt <- melt(plot.dt, id.vars = c(keys, "len"), measure.vars = c("icd", "exd"))
p.dt$variable = factor(p.dt$variable, levels = unique(p.dt$variable), ordered = TRUE)

refscen <- expression(
  movModel == 0 & moveRate == 0 & ES_Detection == 0 & vacRate == 0 &
  !(label %like% "&")
)

marker <- colorRampPalette(c('orange','red','purple','royalblue'))(5)

p <- ggplot(p.dt[eval(refscen)][len>0][between(value, 0.01, 0.99)]) +
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
  geom_line() +
  geom_text(
    aes(label = label, color = NULL),
    data = function(dt) dt[, .(len = 365/10, value = (max(value)-min(value))*.99+min(value), label = c("A", "B")[as.numeric(variable)]), by=variable],
    size = 8,
    show.legend = FALSE
  ) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(
    "Interval Length (years)"
  ) +
  coord_cartesian(xlim=c(0, 3.5), ylim=c(0.01, 0.99), expand = FALSE) +
  scale_y_continuous(
    "Simulated CDF (logit scale)", trans = "logit"#,
    , breaks = c(0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.90, 0.975, 0.99),
minor_breaks = NULL,
labels = scales::label_percent()
  ) +
  scale_color_manual(
    NULL, breaks = c(rev(unique(del.dt$label)), factor("1x64K")), values = c(rev(marker), "black")
  ) +
  theme(legend.position = c(1, 0.5), legend.justification = c(1, 0), legend.direction = "horizontal")

ggsave(tail(.args, 1), p, width = 6, height = 10, dpi = 600, bg = "white")
