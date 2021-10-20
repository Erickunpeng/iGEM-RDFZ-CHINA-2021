library(tidyverse)
library(ggplot2)
library(deSolve)

#Set up
maxtime = 700
interval  = 1
v <-(0:100)
GFP_change <- vector()
SgRNA_dcpf1_change <-vector()
times <- seq(from=0,to=maxtime,by=interval)

# loop: plot the GFP when equilibrium with the change in SgDNA/dCpf1 input
for (i in v) {

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
xstart <- c(SgRNA=i,dcpf1=i,SgRNA.dcpf1=0,DNA=100,DNA.SgRNA.dcpf1=0,GFP=0)
parms <- c(k9=0.001,k10=200,k11=8000,k12=0,k13=4,k14=1)

ode(
  dcas.system.model,
  y=xstart,
  times=times,
  parms=parms
) %>%
  as.data.frame() -> out


GFP_change<-as.data.frame(append(GFP_change,out[maxtime/interval+1,7]))
SgRNA_dcpf1_change <- as.data.frame(append (SgRNA_dcpf1_change,i))
data <- as.data.frame(c(GFP_change, SgRNA_dcpf1_change))
names(data)[1:2] <- c("GFP_change", "SgRNA_dcpf1_change")



ggplot(data = data, aes(x = SgRNA_dcpf1_change, y = GFP_change)) +
  geom_line(color="#6495ED", size=1, linetype=1) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(title = "GFP concentration when equilibrium with the change of SgRNR/dCpf1 input",
       x='SgRNA/dCpf1 concentration(¦ÌM)',
       y='GFP concentration(¦ÌM)') +
  theme(plot.title = element_text(hjust = 0.5))



