
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(RColorBrewer)
  require(patchwork)
  require(ggrepel)
  require(facetscales)
})

.args <- if (interactive()) c(
  file.path("data", "digest.rds"),
  file.path("Figures", "fig9.png")
) else commandArgs(trailingOnly = TRUE)

refscen <- expression(movModel == 0 & ES_Detection == 0 & vacRate == 0 & (moveRate == 0.1 | moveRate == 0))

plot.dt <- readRDS(.args[1])[eval(refscen)]
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

del.dt <- setkeyv(alt.dt[
  ref.dt[,.SD,.SDcols = -c("cksum", "label", "vilModel")],
  on=c(setdiff(keys, c("cksum","label","vilModel")), "len"),
  roll = Inf, rollends = TRUE
], keys)

del.dt[, p.sc := 1 - pext ]
del.dt[, del.p := pext - i.pext ]
del.dt[, OF.p := log10(((1-i.pext)/i.pext * pext/(1-pext))) ]

plot.src <- rbind(
  melt(del.dt, id.vars = c(key(del.dt), "len"), measure.vars = c("p.sc", "del.p", "OF.p")),
  melt(ref.dt[, p.sc := 1 - pext], id.vars = c(key(del.dt), "len"), measure.vars = c("p.sc"))
)

marker <- colorRampPalette(c('orange','red','purple','royalblue'))(5)

exps <- c("Pr(Silent~Circulation)", "Delta*Pr(SC)", "OR~Pr(SC)~(log~scale)")

plot.src[variable == "p.sc", exp := factor(exps[1], levels = exps, ordered = TRUE) ]
plot.src[variable == "del.p", exp := factor(exps[2], levels = exps, ordered = TRUE) ]
plot.src[variable == "OF.p", exp := factor(exps[3], levels = exps, ordered = TRUE) ]

lbls <- plot.src[moveRate != 0, {
  tarlen <- ((.GRP*0.25)+0.5)*365
  .SD[,{
    tarval <- value[which.max(tarlen<len)]
    .(len=tarlen, value=tarval)
  }, by=exp]
},by=.(vilModel, label, moveRate)]

scales_y <- list(y=list(
  "Pr(Silent~Circulation)" = scale_y_continuous(breaks = seq(0,1,by=0.25), minor_breaks = NULL),
  "Delta*Pr(SC)" = scale_y_continuous(breaks = seq(-.2,1,by=.2), minor_breaks = NULL),
  "OR~Pr(SC)~(log~scale)" = scale_y_continuous(breaks = seq(-2,3,by=1), labels = function(b) 10^b)
), x="fixed")

p <- ggplot(plot.src[!is.na(value) & !is.infinite(value)]) +
  aes(len/365, value, color = label, alpha = factor(moveRate), group = interaction(label, moveRate)) +
  facet_grid_sc(
    exp ~ ., scales = scales_y, switch = "y",
    labeller = label_parsed
  ) +
  geom_line() +
  geom_text(
    aes(label = label, color = NULL, group = NULL, alpha = NULL),
    data = function(dt) dt[, .(len = 365/10, value = (max(value)-min(value))*.90+min(value), label = c("A", "B", "C")[as.numeric(exp)]), by=exp],
    size = 8,
    show.legend = FALSE
  ) +
  geom_blank(aes(color = NULL, alpha = NULL, group = NULL), data = function(dt) dt[,.(len = c(0,3.5)*365, value = 0),by=exp]) +
  geom_label_repel(
    aes(label=label, group = NULL), data=lbls, show.legend = FALSE,
    fill = alpha(c("white"), 0.75),
    #force = 0.5, force_pull = 2,
    box.padding = 0,
    direction = "y"
  ) +
  coord_cartesian(expand = FALSE) +
  scale_color_manual(
    NULL, breaks = c(rev(unique(del.dt$label)), factor("1x64K")), values = c(rev(marker), "black")
  ) +
  scale_x_continuous(
    "Time Since Case Observed (years)"
  ) +
  scale_alpha_manual(NULL, guide = "none", values = c(`0`=0.25, `0.1`=1)) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outer",
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing = unit(1.5, units = "line")
  )

crossbits <- list(
  aes(len/365, value, color = label, alpha = factor(moveRate), group = interaction(label, moveRate)),
  #  geom_line(),
  scale_color_manual(
    NULL, breaks = c(rev(unique(del.dt$label)), factor("1x64K")), values = c(rev(marker), "black")
  ),
  scale_x_continuous(NULL, breaks = seq(2.5, 3.5, by=.5)),
  scale_y_continuous(NULL),
  scale_alpha_manual(NULL, guide = "none", values = c(`0`=0.25, `0.1`=1)),
  coord_cartesian(expand = FALSE),
  theme_minimal(),
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.1),
    # plot.margin = margin(l = 0.1, r = 0.1, unit = "mm"),
    plot.margin = margin(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
)

theme_trim <- list(
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "null"),
    axis.title = element_blank(),
    axis.line = element_blank()
  ),
  geom_blank(aes(color = NULL, alpha = NULL, group = NULL), data = data.table(len = c(2.5, 3.5)*365, value = 0))
)

pinset1 <- ggplot(
  plot.src[!is.na(value) & !is.infinite(value)][variable == "p.sc"][between(len/365, 2.5, 3.5)]
) + crossbits + geom_line() + 
  theme_trim

pinset1axis <- ggplot(
  plot.src[!is.na(value) & !is.infinite(value)][variable == "p.sc"][between(len/365, 2.5, 3.5)]
) + crossbits + geom_blank(data = function(dt) dt[,.(len = rep(len[1],2), value = range(value)),by=.(moveRate, label)]) +theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.length.x = unit(0, "null"),
  axis.title.x = element_blank(),
  axis.line.x = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank()
)

pinset2 <- ggplot(
  plot.src[!is.na(value) & !is.infinite(value)][variable == "del.p"][between(len/365, 2.5, 3.5)]
) + crossbits + geom_line() + theme_trim

pinset2axis <- ggplot(
  plot.src[!is.na(value) & !is.infinite(value)][variable == "del.p"][between(len/365, 2.5, 3.5)]
) + crossbits + geom_blank(data = function(dt) dt[,.(len = rep(len[1],2), value = range(value)),by=.(moveRate, label)]) +theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.length.x = unit(0, "null"),
  axis.title.x = element_blank(),
  axis.line.x = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank()
)

resp <- p +
  inset_element(pinset1, left = 2.5/3.5, right = 1, top = 1, bottom = 0.66+0.33*0.3) +
  inset_element(pinset1axis, left = 2.5/3.5-0.06, right = 2.5/3.5, top = 1, bottom = 0.66+0.33*0.3) +
  inset_element(pinset2, left = 2.5/3.5, right = 1, top = 0.66, bottom = 0.33+0.33*0.3) +
  inset_element(pinset2axis, left = 2.5/3.5-0.05, right = 2.5/3.5, top = 0.66, bottom = 0.33+0.33*0.3)# +

ggsave(tail(.args, 1), resp, width = 7, height = 10, dpi = 600, bg = "white")
