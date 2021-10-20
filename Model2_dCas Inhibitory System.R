library(tidyverse)
library(ggplot2)
library(deSolve)


times <- seq(from=0,to=150,by=0.01)


dcas.system.model <- function (t, x, params) {
  require(deSolve)
  SgRNA <- x[1]
  dcpf1 <- x[2]
  SgRNA.dcpf1 <- x[3]
  DNA <- x[4]
  DNA.SgRNA.dcpf1 <-x[5]
  GFP <-x[6]
  
  #SgRNA + dcpf1 <-> SgRNA.dcpf1
  #SgRNA.dcpf1 + DNA <-> DNA.SgRNA.dcpf1
  #DNA-> DNA+ GFP 
  #GFP -> N/A
  
  ## now extract the parameters

  k9 <- params["k9"]
  k10 <- params["k10"]
  k11 <- params["k11"]
  k12 <- params["k12"]
  k13 <- params["k13"]
  k14 <- params["k14"]
  
  ## now code the model equations
      dSgRNAdt <- k10*x[3]-k9*x[1]*x[2]
      ddcpf1dt <- k10*x[3]-k9*x[1]*x[2]
      dSgRNA.dcpf1dt <- k9*x[1]*x[2]-k10*x[3]+k12*x[5]-k11*x[3]*x[4]
      dDNAdt <- -k11*x[3]*x[4] + k12*x[5]
      dDNA.SgRNA.dcpf1dt <- k11*x[3]*x[4]-k12*x[5]
      dGFPdt <-k13*x[4]-k14*x[6]
    
      dxdt <- c(dSgRNAdt,ddcpf1dt,dSgRNA.dcpf1dt,dDNAdt,dDNA.SgRNA.dcpf1dt,dGFPdt)
      ## return result as a list
      ode <- list(dxdt)

}
xstart <- c(SgRNA=50,dcpf1=50,SgRNA.dcpf1=0,DNA=100,DNA.SgRNA.dcpf1=0,GFP=0)
parms <- c(k9=0.001,k10=200,k11=8000,k12=0,k13=4,k14=1)

ode <- ode(
  dcas.system.model,
  y=xstart,
  times=times,
  parms=parms
)
  
ode <- as.data.frame(ode)
SgRNA_data <- as.data.frame(ode$SgRNA)
GFP_data <- as.data.frame(ode$GFP)
GFP_SgRNA <- as.data.frame(c(GFP_data, SgRNA_data))

ode_data <- gather(ode, substance,value,-time) 

ggplot(data=ode_data, aes(x=time,y=value,color=substance))+
  geom_line(size=1)+
  labs(title = "dCas Inhibitory System", x='time (s)',y='concentration (¦ÌM)') +
  theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(0.5, "inches")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 2))
  


