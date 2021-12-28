
require(data.table)
require(ggplot2)

.args <- if (interactive()) c(
  file.path("data", "parameters.txt"),
  file.path("data", "output"),
  file.path("data", "digest.rds")
)

dt <- fread(.args[1])[,.(
  cksum, movModel, moveRate, ES_Detection, vacRate, vilPop
)]

dt[, label :=  {
  pops <- rle(strsplit(gsub("[{}]","",gsub("000([},])","K\\1",vilPop)),",")[[1]])
  paste(sprintf("%sx%i", pops$values, pops$lengths), collapse=" & ")
}, by=vilPop]

dt[, vilModel := .GRP, by=label]
dt$vilPop <- NULL

outputs <- grep("^.+_\\d+_.+\\.out$", list.files(.args[2], "cts.out$", full.names = TRUE), invert = TRUE, value = TRUE)

exthash <- function(fl) { gsub("^([^_]+)_.+$","\\1",basename(fl)) }
#' TODO show cross sample diversity?
injest <- function(fl) {
  res <- dcast(
    fread(fl), serial + len ~ type, value.var = "ct", fun.aggregate = sum
  )[, .(ex = sum(ex), ic = sum(ic)), keyby = len]
  res[, cksum := exthash(fl) ]
  # res[len > 0, c("cex", "cic") := .(cumsum(ex), cumsum(ic)) ]
  # res[, pext := cex/(cex[.N]+cic[.N]-c(0, cex[-.N])-c(0, cic[-.N])) ]
  res
}

#' ~20Gb with full data set
plot.dt <- rbindlist(lapply(outputs, injest))[dt, on=.(cksum) ]
plot.dt[len >= 0, c("cex", "cic") := .(cumsum(ex), sum(ic) - cumsum(c(0, ic[-.N]))), by = cksum ]
plot.dt[, pext := cex/(cex + cic), by=cksum ]

saveRDS(plot.dt, tail(.args, 1))
