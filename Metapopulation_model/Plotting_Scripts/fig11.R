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

model <- function(t,y,dy,pars){
  with(as.list(c(y,dy,pars)),{
    
    trans = beta*(y[2]+kappa*y[6])/pop_size
    
    I1_inf = y[1]*trans
    Ir_inf = kappa*y[5]*trans
    
    s_eq = dy[1] - mu*pop_size + vacc_rate*y[1] + I1_inf + mu*y[1]
    i1_eq = dy[2] - I1_inf + (gamma + mu)*y[2]
    v_eq = dy[3] - vacc_rate*y[1] + mu*y[3]
    r_eq = dy[4] - gamma*y[2] - (gamma/kappa)*y[6] + (omega + mu)*y[4]
    p_eq = dy[5] - omega*y[4] + Ir_inf + mu*y[5]
    ir_eq = dy[6] - Ir_inf + ((gamma/kappa) + mu)*y[6]
    constraint = pop_size - (y[1]+y[2]+y[3]+y[4]+y[5]+y[6])
    
    list(c(s_eq,i1_eq,v_eq,r_eq,p_eq,ir_eq))
    
  })
}

find_equilibrium <- function(pop_size, vacc_rate){
  #input: total population size, vaccination rate
  #output: proportion of individuals in each compartment in the order
  #S,I1,V,R,P,Ir
  
  #model parameters
  pars <- c(mu = 0.02, beta = 135, gamma = 13, omega = 0.2, kappa = 0.4179,
            pop_size = pop_size,vacc_rate=vacc_rate)
  
  #initial starting conditions for compartments
  yini <-c(S=(pop_size - 1),I1=1,V=0,R=0,P=0,Ir=0)
  
  #initial starting conditions for differential equation
  #set to zero b/c want equilibrium
  dyini <- c(0,0,0,0,0,0)
  
  #length of time to run the model
  times = seq(0,10000,0.1)
  
  #run the model
  out<-daspk(y=yini,dy=dyini,res=model,times=times,parms=pars)
  return(as.list(pmax(out[dim(out)[1],-1]/pop_size, 0)))
}

eq.dt <- melt(plot.dt[
  movModel == 0 & moveRate == 0 & ES_Detection == 0 & village_size == 32000
][,
  find_equilibrium(pop_size = village_size, vacc_rate = vacRate),
  by=.(village_size, vacRate)
], id.vars = c("village_size", "vacRate"))

marker <- colorRampPalette(c('orange','red','purple','royalblue'))(5)

p <- ggplot(plot.dt[movModel == 0 & moveRate == 0 & ES_Detection == 0]) +
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
    data = eq.dt,
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
