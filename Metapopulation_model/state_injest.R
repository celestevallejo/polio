
require(data.table)

.args <- if (interactive()) c(
  file.path("data", "parameters.txt"),
  file.path("data", "output"),
  file.path("data", "state_digest.rds")
) else commandArgs(trailingOnly = TRUE)

dt <- fread(.args[1])[,.(
  cksum, movModel, moveRate, ES_Detection, vacRate, vilPop
)]

dt[, label :=  {
  pops <- rle(strsplit(gsub("[{}]","",gsub("000([},])","K\\1",vilPop)),",")[[1]])
  paste(sprintf("%sx%i", pops$values, pops$lengths), collapse=" & ")
}, by=vilPop]

dt[, vilModel := .GRP, by=label]
dt$vilPop <- NULL

outputs <- grep("^.+_\\d+_.+\\.out$", list.files(.args[2], "states.out$", full.names = TRUE), invert = TRUE, value = TRUE)

exthash <- function(fl) { gsub("^([^_]+)_.+$","\\1",basename(fl)) }
injest <- function(fl) {
  res <- melt(
    fread(fl)[, serial := as.integer(as.hexmode(serial)) ],
    id.vars = c("serial", "village_id", "village_size")
  )
  qtile <- res[, {
    qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(qs) <- c("lo95", "lo50", "md", "hi50", "hi95")
    as.list(qs)
  }, by=.(village_size, variable)]
  qtile[, cksum := exthash(fl) ]
  qtile
}

#' ~20Gb with full data set
plot.dt <- rbindlist(lapply(outputs, injest))[dt, on=.(cksum) ]

saveRDS(plot.dt, tail(.args, 1))