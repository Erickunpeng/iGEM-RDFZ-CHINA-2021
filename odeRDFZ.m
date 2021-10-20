function dydt = odeRDFZ(~,y,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14)
dydt=[
%RpfG
- k5*y(7)*y(1) + k6*y(10) + k13*y(2)*y(9) - k14*y(4)*y(1)*y(8);
%RpfGP
k7*y(10) - k8*y(6)*y(2) - k13*y(2)*y(9) + k14*y(4)*y(1)*y(8);
%RpfC
- k1*y(3)*y(5) + k2*y(6);
%cdigmp
k9 - k10*y(4) - k11*y(4)*y(8) + k12*y(9) - k13*y(2)*y(9) + k14*y(4)*y(1)*y(8);
%DSF
- k1*y(3)*y(5) + k2*y(6);
%DSF.RpfC 
k1*y(3)*y(5) - k2*y(6) - k3*y(6)+ k4*y(7) + k7*y(10) - k8*y(6)*y(2);
%DSF.RpfCP
k3*y(6)- k4*y(7) - k5*y(7)*y(1) + k6*y(10); 
%vc2
- k11*y(4)*y(8) + k12*y(9)+k13*y(2)*y(9)-k14*y(4)*y(1)*y(8);
%cdigmp.vc2
k11*y(4)*y(8) - k12*y(9)+k14*y(4)*y(1)*y(8)-k13*y(2)*y(9);
%DSF.RpfC.RpfGP
k5*y(7)*y(1) - k6*y(10) - k7*y(10) + k8*y(6)*y(2);

    ];
end



