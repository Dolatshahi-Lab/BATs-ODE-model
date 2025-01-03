function dy = ddefun_MA_nogen_ligandvariable_Jul_nodel(t,y,p,MN2,ML2,MP2,MD2,CN2,CL2,CP2,CD2,TN2,TL2,TP2,TD2) % equation being solved
    
    dy=zeros(27,1);
   

%%MCF7
%dTumor/dt y(1)
dy(1) = p.Matp*y(1)-(p.Matp*y(1)^2)/p.Mmax...
       - p.aKill_N*y(2)*y(1)...
       - p.aKill_L*y(3)*y(1)...
       - p.aKill_P*y(4)*y(1)...
       - p.aKill_D*y(5)*y(1);

%dBAT_N/dt y(5)
dy(2) = ppval(MN2,t);

%dBAT_L/dt y(8)
dy(3) = ppval(ML2,t);

%dBAT_P/dt y(11)
dy(4) = ppval(MP2,t);

%dBAT_D/dt y(11)
dy(5) = ppval(MD2,t);

dy(6) = p.aKill_N*y(2)*y(1);
dy(7) = p.aKill_L*y(3)*y(1);
dy(8) = p.aKill_P*y(4)*y(1);
dy(9) = p.aKill_D*y(5)*y(1);

%%CAMA1
%dTumor/dt y(1)
dy(10) = p.Catp*y(10)-(p.Catp*y(10)^2)/p.Cmax...
       - p.aKill_N*y(11)*y(10)...
       - p.aKill_L*y(12)*y(10)...
       - p.aKill_P*y(13)*y(10)...
       - p.aKill_D*y(14)*y(10);

%dBAT_N/dt y(5)
dy(11) = ppval(CN2,t);

%dBAT_L/dt y(8)
dy(12) = ppval(CL2,t);

%dBAT_P/dt y(11)
dy(13) = ppval(CP2,t);

%dBAT_D/dt y(11)
dy(14) = ppval(CD2,t);

dy(15) = p.aKill_N*y(11)*y(10);
dy(16) = p.aKill_L*y(12)*y(10);
dy(17) = p.aKill_P*y(13)*y(10);
dy(18) = p.aKill_D*y(14)*y(10);


%%T47D
%dTumor/dt y(1)
dy(19) = p.Tatp*y(19)-(p.Tatp*y(19)^2)/p.Tmax...
       - p.aKill_N*y(20)*y(19)...
       - p.aKill_L*y(21)*y(19)...
       - p.aKill_P*y(22)*y(19)...
       - p.aKill_D*y(23)*y(19);

%dBAT_N/dt y(5)
dy(20) = ppval(TN2,t);

%dBAT_L/dt y(8)
dy(21) = ppval(TL2,t);

%dBAT_P/dt y(11)
dy(22) = ppval(TP2,t);

%dBAT_D/dt y(11)
dy(23) = ppval(TD2,t);

dy(24) = p.aKill_N*y(20)*y(19);
dy(25) = p.aKill_L*y(21)*y(19);
dy(26) = p.aKill_P*y(22)*y(19);
dy(27) = p.aKill_D*y(23)*y(19);


end 