#run model to find equilibrium values

#differential algebraic equation (ODE + constraint)


library(deSolve)

model <- function(t,y,dy,pars){
  with(as.list(c(y,dy,pars)),{
    
    trans = beta*(y[2]+kappa*y[6])
    
    I1_inf = y[1]*trans
    Ir_inf = kappa*y[5]*trans
    
    s_eq = dy[1] - mu + vacc_rate*y[1] + I1_inf + mu*y[1]
    i1_eq = dy[2] - I1_inf + (gamma + mu)*y[2]
    v_eq = dy[3] - vacc_rate*y[1] + mu*y[3]
    r_eq = dy[4] - gamma*y[2] - (gamma/kappa)*y[6] + (omega + mu)*y[4]
    p_eq = dy[5] - omega*y[4] + Ir_inf + mu*y[5]
    ir_eq = dy[6] - Ir_inf + ((gamma/kappa) + mu)*y[6]
    constraint = 1 - (y[1]+y[2]+y[3]+y[4]+y[5]+y[6])
    
    list(c(s_eq,i1_eq,v_eq,r_eq,p_eq,ir_eq))
    
  })
}

find_equilibrium <-function(vacc_rate){
  #input: vaccination rate
  #output: proportion of individuals in each compartment in the order
  #S,I1,V,R,P,Ir
  
  #model parameters
  pars <- c(mu = 0.02, beta = 135, gamma = 13, omega = 0.2, kappa = 0.4179,
            vacc_rate=vacc_rate)
  
  #initial starting conditions for compartments
  yini <-c(0.99,0.1,0,0,0,0)
  #yini <-c((pop_size - 500),10,0,0,0,0)
  
  #initial starting conditions for differential equation
  #set to zero b/c want equilibrium
  dyini <- c(0,0,0,0,0,0)
  
  #length of time to run the model
  times = seq(0,10000,0.1)
  
  #run the model
  out<-daspk(y=yini,dy=dyini,res=model,times=times,parms=pars)

  return(data.frame(S = out[length(times),2],
                    I1 = out[length(times),3],
                    V = out[length(times),4],
                    R = out[length(times),5],
                    P = out[length(times),6],
                    Ir = out[length(times),7]))
}






