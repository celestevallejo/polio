
require(data.table)
require(ggplot2)

.args <- if (interactive()) c(
	"~/Downloads/qqpolio",
	"data/parameters.rds",
	"Figures/qqpolio.png"
) else commandArgs(trailingOnly = TRUE)

qqfls <- list.files(.args[1], full.names = TRUE)
names(qqfls) <- gsub("^.*_(\\w+)\\.csv$", "\\1", qqfls)

pars <- readRDS(.args[2])

qq.dt <- rbindlist(
	lapply(qqfls, fread, col.names = c("outcome", "I1", "Ir")),
	idcol = "cksum"
)[outcome %in% c("pe","pp")][pars, on=.(cksum), nomatch = 0]

plot.dt <- setkey(rbind(qq.dt[,
	.(N=.N, measure = "I1"), keyby=.(label, outcome = factor(outcome, c("pp","pe"), ordered = TRUE), count = I1)
], qq.dt[,
	.(N=.N, measure = "Ir"), keyby=.(label, outcome = factor(outcome, c("pp","pe"), ordered = TRUE), count = Ir)
]), label, outcome, measure, count)

plot.dt[order(count), q := cumsum(N)/sum(N), by=.(label, outcome, measure)]
plot.dt[, refval := qgeom(q, prob = 1/200) ]

marker <- colorRampPalette(c('orange','red','purple','royalblue'))(5)

pmain <- ggplot(plot.dt[measure == "I1"]) +
	aes(count, refval, color = label) +
	facet_grid(outcome ~ measure, labeller = labeller(
		outcome = c(pe="Extinction", pp="Between Observations")
	)) +
	coord_equal(xlim=c(0,3000), ylim=c(0,3000), expand = FALSE) +
	geom_point(alpha=0.5) + theme_minimal() +
	scale_color_manual(
		NULL, breaks = c(rev(unique(plot.dt$label)), factor("1x64K")), values = c(rev(marker), "black"),
		guide = guide_legend(override.aes = list(alpha = 1))
	) +
	scale_x_continuous("Simulated Quantiles", breaks = seq(0,3000,by=500)) +
	scale_y_continuous("Geometric Quantiles", breaks = seq(0,3000,by=500)) +
	theme(
		legend.position = c(0, 1),
		legend.justification = c(0, 1),
		panel.spacing = unit(1, units = "line")
	)

ggsave(tail(.args, 1), pmain, width = 5, height = 10, dpi = 600, bg = "white")

psi <- ggplot(plot.dt[measure == "Ir"]) +
	aes(count, refval, color = label) +
	facet_grid(outcome ~ measure, labeller = labeller(
		outcome = c(pe="Extinction", pp="Between Observations")
	)) +
	coord_equal(xlim=c(0,3000), ylim=c(0,3000), expand = FALSE) +
	geom_point(alpha=0.5) + theme_minimal() +
	scale_color_manual(
		NULL, breaks = c(rev(unique(plot.dt$label)), factor("1x64K")), values = c(rev(marker), "black"),
		guide = guide_legend(override.aes = list(alpha = 1))
	) +
	scale_x_continuous("Simulated Quantiles", breaks = seq(0,3000,by=500)) +
	scale_y_continuous("Geometric Quantiles", breaks = seq(0,3000,by=500)) +
	theme(
		legend.position = c(0, 1),
		legend.justification = c(0, 1),
		panel.spacing = unit(1, units = "line")
	)

ggsave(gsub("\\.","_alt.",tail(.args, 1)), psi, width = 5, height = 10, dpi = 600, bg = "white")
