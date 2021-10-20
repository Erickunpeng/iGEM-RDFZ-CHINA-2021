
%% Main funtions 

%   #DSF (x5) + RpfC (x3) <-> DSF~RpfC (x6) (k1,k2)
%   #DSF~RpfC (x6) <-> DSF~RpfC(P) (x7) (k3,k4)
%   #DSF~RpfC(P) (x7) + RpfG (x1) <-> DSF~RpfC~RpfG(P) (x10) (k5,k6)
%   #DSF~RpfC~RpfG(P) (x10) <-> DSF~RpfC (x6) + RpfG(P) (x2) (k7,k8)
%   #N/A -> cGMP (x4) (k9)
%   #cGMP (x4) -> N/A (k10)
%   #cGMP (x4) + vc2 (x8) <-> cGMP~vc2 (x9) (k11,k12)
%   #RpfG(P) (x2) + cGMP~vc2 (x9) <-> RpfG(x1) + vc2 + cGMP (k13,k14)

%% Total value 


DSF_total =100;
RpfG_total =100;
RpfC_total=100;
cdigmp_vc2=100;

%% Constant assigning 

k1=42;
k2=k1*20;   
k3=9*10^3;
k4=k3*10^-3;
k5=2;
k6=k5*200; 
k7=1.7*10^2;  
k8=5*10^-4;
k9=9.2*10^1; 
k10=8.5*10^-5;
k11=2*10^2; 
k12=3*10^4;
k13=2000; 
k14=1*10^-4;

%% Initial value

RpfG=RpfG_total;
RpfGP=0;
RpfC=RpfC_total;
cdigmp=0;
DSF=DSF_total;
DSF_RpfC=0;
DSF_RpfCP=0;
vc2 = 0;
cdigmp_vc2 = cdigmp_vc2;
DSF_RpfC_RpfGP = 0;

%% ODE


tspan =(0:2);
y0=[RpfG;RpfGP;RpfC;cdigmp;DSF;DSF_RpfC;DSF_RpfCP;vc2;cdigmp_vc2;DSF_RpfC_RpfGP];
[t,y]=ode45(@(t,y)odeRDFZ(t,y,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14),tspan,y0);

%% Plot 

plot(t,y);
plot(t,y(:,1),'b-', t,y(:,2),'r-', t,y(:,3),'k-',t,y(:,4),'-',t,y(:,5),'-',t,y(:,6),'-',t,y(:,7),'-',t,y(:,8),'-',t,y(:,9),'-',t,y(:,10),'-');
legend('RpfG','RpfGP','RpfC','cdigmp','DSF','DSF·RpfC','DSF·RpfCP','vc2','cdigmp·vc2','DSF·RpfC·RpfGP');
xlabel('time(s)')
ylabel('concentration(μM)');



