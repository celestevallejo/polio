
-include local.makefile

RSCRIPT = $(shell which Rscript)

R = $(strip ${RSCRIPT} $^ $(1) $@)

SRCDIR ?= ~/Downloads/tmp
DATADIR ?= data
FIGDIR ?= Figures

data/output: | ${SRCDIR}
	ln -s $| $@

${DATADIR}/digest.rds: injest.R data/parameters.txt data/output
	$(call R)

${DATADIR}/state_digest.rds: state_injest.R data/parameters.txt data/output
	$(call R)

${FIGDIR}/fig7.png: fig7.R ${DATADIR}/digest.rds:
	$(call R)

${FIGDIR}/fig11.png: fig11.R ${DATADIR}/state_digest.rds:
	$(call R)

# figure 11: initial distribution of states
# uses one scenario, for multiple village sizes
figure11.png: starting_conditions.R
