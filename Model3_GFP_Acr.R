library(tidyverse)
library(ggplot2)
library(deSolve)

#Set up
maxtime = 200
interval  = 1
v <-(0:500)
GFP_change <- vector()
Acr_change <-vector()
times <- seq(from=0,to=maxtime,by=interval)

# loop: plot the GFP when equilibrium with the change in SgDNA/dCpf1 input
for (i in v){

dcas.system.model <- function (t, x, params) {
  require(deSolve)
  SgRNA <- x[1]
  dcpf1 <- x[2]
  SgRNA.dcpf1 <- x[3]
  DNA <- x[4]
  DNA.SgRNA.dcpf1 <-x[5]
  GFP <-x[6]
  Acr <- x[7]
  
  #SgRNA + dcpf1 <-> SgRNA.dcpf1
  #SgRNA.dcpf1 + DNA <-> DNA.SgRNA.dcpf1(k11/k12)
  #DNA-> DNA+ GFP 
  #GFP -> N/A
  #vc2->vc2 + Acr (k15)
  #Acr + SgRNA.dcpf1 -> SgRNA.dcpf1* (k16)
  #Acr -> N/A (k17)
  
  
  ## now extract the parameters
  
  k9 <- params["k9"]
  k10 <- params["k10"]
  k11 <- params["k11"]
  k12 <- params["k12"]
  k13 <- params["k13"]
  k14 <- params["k14"]
  k15 <- params["k15"]
  k16 <- params["k16"]
  k17 <- params["k17"]
  
  ## now code the model equations
  dSgRNAdt <- k10*x[3]-k9*x[1]*x[2]
  ddcpf1dt <- k10*x[3]-k9*x[1]*x[2]
  dSgRNA.dcpf1dt <- k9*x[1]*x[2]-k10*x[3]+k12*x[5]-k11*x[3]*x[4]- k16*x[3]*x[7] 
  dDNAdt <- -k11*x[3]*x[4] + k12*x[5]
  dDNA.SgRNA.dcpf1dt <- k11*x[3]*x[4]-k12*x[5]
  dGFPdt <-k13*x[4]-k14*x[6]
  dAcr <- k15*i - k16*x[3]*x[7]-k17*x[7]
  
  
  dxdt <- c(dSgRNAdt,ddcpf1dt,dSgRNA.dcpf1dt,dDNAdt,dDNA.SgRNA.dcpf1dt,dGFPdt,dAcr)
  ## return result as a list
  ode <- list(dxdt)
  
}
xstart <- c(SgRNA=100,dcpf1=100,SgRNA.dcpf1=0,DNA=100,DNA.SgRNA.dcpf1=0,GFP=0,Acr=i)
parms <- c(k9=0.001,k10=200,k11=8000,k12=0,k13=4,k14=1,k15=4,k16=350,k17= 0.0001)

ode(
  dcas.system.model,
  y=xstart,
  times=times,
  parms=parms
) %>%
  as.data.frame() -> out
out %>% gather(substance,value,-time) %>% ggplot(aes(x=time,y=value,color=substance))+geom_line(size=1)+theme_classic()+
  labs(x='time (s)',y='concentration')

GFP_change<-append(GFP_change,out[maxtime/interval+1,7])
Acr_change <-append (Acr_change,i)

}

GFP_change <- as.data.frame(append(GFP_change,out[maxtime/interval+1,7]))
Acr_change <- as.data.frame(append (Acr_change,i))
GFP_Acr <- as.data.frame(c(GFP_change, Acr_change))
names(GFP_Acr)[1:2] <- c("GFP_change", "Acr_change")

ggplot(data = GFP_Acr, aes(x = Acr_change, y = GFP_change)) + 
  geom_line(color="#6495ED", size=1, linetype=1) +
  labs(title = "Effect of Acr input on downstream protein expression (GFP)", 
       x='Acr concentration(¦ÌM)', y='GFP concentration(¦ÌM)') +
  theme(plot.title = element_text(hjust = 0.5))
