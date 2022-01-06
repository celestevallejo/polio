
require(data.table)

.args <- if (interactive()) c(
  file.path("data", "parameters.txt"),
  file.path("data", "parameters.rds")
)

dt <- fread(.args[1])[,.(
  cksum, movModel, moveRate, ES_Detection, vacRate, vilPop
)]

dt[, c("label", "patches") :=  {
  pops <- rle(strsplit(gsub("[{}]","",gsub("000([},])","K\\1",vilPop)),",")[[1]])
  .(
    paste(sprintf("%ix%s", pops$lengths, pops$values), collapse=" & "),
    sum(pops$lengths)
  )
}, by=vilPop]

unilabels <- dt[
  order(patches),
  unique(label)
]

dt[
  order(patches),
  flabel := factor(label, levels = unilabels, ordered = TRUE)
]

dt[, vilModel := .GRP, keyby = flabel]
dt$label <- dt$flabel
dt$vilPop <- dt$patches <- dt$flabel <- NULL

saveRDS(dt, tail(.args, 1))