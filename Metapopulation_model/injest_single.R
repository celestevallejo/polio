
require(data.table)

.args <- if (interactive()) c(
  "~/Downloads/all_parameters.txt",
  "~/Downloads/morepolio/tmp",
  file.path("data", "singledigest.rds")
)

dt <- fread(.args[1])[,.(
  cksum, movModel, moveRate, ES_Detection, vacRate, vilPop
)][, c("label", "patches", "vilsize") :=  {
  pops <- rle(strsplit(gsub("[{}]","",gsub("000([},])","K\\1",vilPop)),",")[[1]])
  .(
    paste(sprintf("%ix%s", pops$lengths, pops$values), collapse=" & "),
    sum(pops$lengths),
    as.integer(gsub("K","",pops$values[1]))
  )
}, by=vilPop][patches == 1]

unilabels <- dt[
  order(vilsize),
  unique(label)
]

dt[
  order(patches),
  flabel := factor(label, levels = unilabels, ordered = TRUE)
]

dt[, vilModel := .GRP, keyby = flabel]
dt$label <- dt$flabel
dt$vilPop <- dt$patches <- dt$flabel <- NULL

merge.fls <- list.files(.args[2], ".+_\\d+_.+\\.out", full.names = TRUE)

uniqs <- unique(gsub("(.+)_\\d+_.+\\.out", "\\1", merge.fls))

big <- list()

for (uid in uniqs) {
  tars <- grep(uid, merge.fls, value = TRUE)[1:20]
  big[[basename(uid)]] <- setkey(rbindlist(
      lapply(tars, function(fl) {
        ret <- fread(fl)
        rng <- gsub(".+_(\\d+)_.+\\.out", "\\1", fl)
        ret[, seedid := as.numeric(rng) ]
        ret[, serial := as.integer(as.hexmode(serial)) ]
        ret
      })
    ), seedid, serial, type, len)
}

tot.dt <- rbindlist(lapply(names(big), function(hsh) {
  intct <- big[[hsh]]
  intct[, newserial := .GRP - 1, keyby=.(seedid, serial)]
  intct$serial <- intct$newserial
  intct$seedid <- intct$newserial <- NULL
  dcast(
    intct[serial < 10000], serial + len ~ type, value.var = "ct", fun.aggregate = sum
  )[, cksum := hsh ][len >= 0]
}))[, .(ex = sum(ex), ic = sum(ic)), keyby = .(len, cksum)]

#' ~20Gb with full data set
plot.dt <- tot.dt[dt, on=.(cksum), nomatch = 0]
plot.dt[len >= 0, c("cex", "cic") := .(cumsum(ex), sum(ic) - cumsum(c(0, ic[-.N]))), by = cksum ]
plot.dt[, pext := cex/(cex + cic), by=cksum ]

saveRDS(plot.dt, tail(.args, 1))
