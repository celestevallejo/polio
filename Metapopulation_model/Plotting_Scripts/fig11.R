require(data.table)
require(ggplot2)

.args <- if (interactive()) c(
  file.path("data", "state_digest.rds"),
  file.path("Figures", "fig11.png")
) else commandArgs(trailingOnly = TRUE)

plot.dt <- readRDS(.args[1])

#run model to find equilibrium values
#differential algebraic equation (ODE + constraint)
library(deSolve)

dmodel <- function(y, pars) with(as.list(c(y, pars)), {
  trans = beta*(I1+kappa*IR)
  
  I1_inf = S*trans
  IR_inf = kappa*P*trans
  
  death <- -mu*y
  birth <- c((1-vacc_rate)*mu, 0, vacc_rate*mu, 0, 0, 0)
  inf1 <- c(-I1_inf, I1_inf, 0, 0, 0, 0)
  rinf <- c(0, 0, 0, 0, -IR_inf, IR_inf)
  rec <- c(0, -gamma*I1, 0, gamma*I1+gamma/kappa*IR, 0, -gamma/kappa*IR)
  wane <- c(0, 0, 0, -omega*R, omega*R, 0)
  return(death + birth + inf1 + rinf + rec + wane)
})

model <- function(t,y,dy,pars) return(list(dy - dmodel(y, pars)))

find_equilibrium <- function(vacc_rate){
  #input: total population size, vaccination rate
  #output: proportion of individuals in each compartment in the order
  #S,I1,V,R,P,IR
  
  #model parameters
  pars <- c(mu = 0.02, beta = 135, gamma = 13, omega = 0.2, kappa = 0.4179,
            vacc_rate=vacc_rate)
  
  #initial starting conditions for compartments
  yini <-c(S=0.90-vacc_rate, I1=0.1, V=vacc_rate, R=0, P=0, IR=0)
  
  # initial starting conditions for differential equation
  dyini <- dmodel(yini, pars)
  
  #length of time to run the model
  times = seq(0,10000,0.1)
  
  #run the model
  out<-daspk(y=yini,dy=dyini,res=model,times=times,parms=pars)
  return(as.list(out[dim(out)[1],-1]))
}

eq.dt <- plot.dt[
  movModel == 0 & moveRate == 0 & ES_Detection == 0 & village_size == 32000
][,
  find_equilibrium(vacc_rate = vacRate),
  by=.(vacRate)
]

marker <- colorRampPalette(c('orange','red','purple','royalblue'))(5)

p <- ggplot(plot.dt[movModel == 0 & moveRate == 1 & ES_Detection == 0]) +
  aes(
    factor(village_size),
    color = label
  ) + facet_grid(vacRate ~ variable) +
  geom_boxplot(
    aes(
      ymin =lo95/village_size,
      lower =lo50/village_size,
      middle = md/village_size,
      upper = hi50/village_size,
      ymax = hi95/village_size,
      group = interaction(cksum, village_size)
    ),
    position = "dodge", stat = "identity"
  ) +
  geom_hline(
    aes(x=NULL, yintercept = value, color = NULL),
    data = melt(eq.dt, id.vars = "vacRate"),
    color = "red", show.legend = FALSE
  ) +
  theme_minimal(base_size = 14) +
  scale_x_discrete(
    "Village Size",
    labels = function(f) gsub("000","K", as.character(f))
  ) +
  coord_cartesian(ylim=c(0.0001, 0.6)) +
  scale_y_continuous(
    "% Compartment (logit scale)", trans = "logit",
    breaks = c(0.0001, 0.001, 0.01, 0.1, 0.25, 0.5)
  ) +
  scale_color_discrete(NULL) +
  theme(legend.position = "bottom")

ggsave(tail(.args, 1), p, width = 14, height = 7, dpi = 600, bg = "white")
