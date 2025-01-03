function dy = ddefun_MA_nogen_ligandvariable_Jul_nodel(t,y,p,NC_2,LC_2,PC_2,DC_2,ND_2,LD_2,PD_2,DD_2) % equation being solved
    
    dy=zeros(18,1);
   
% a = 14.9021;
% b = 9.92e-5;
% c = 5.24e4;

%dTumor/dt y(1)
dy(1) = p.atp_C*y(1)-(p.atp_C*y(1)^2)/p.Cmax...
       - p.aKill_N*y(2)*y(1)...
       - p.aKill_L*y(3)*y(1)...
       - p.aKill_P*y(4)*y(1)...
       - p.aKill_D*y(5)*y(1);

%dBAT_N/dt y(5)
dy(2) = ppval(NC_2,t);

%dBAT_L/dt y(8)
dy(3) = ppval(LC_2,t);

%dBAT_P/dt y(11)
dy(4) = ppval(PC_2,t);

%dBAT_D/dt y(11)
dy(5) = ppval(DC_2,t);

dy(6) = p.aKill_N*y(2)*y(1);
dy(7) = p.aKill_L*y(3)*y(1);
dy(8) = p.aKill_P*y(4)*y(1);
dy(9) = p.aKill_D*y(5)*y(1);

%dTumor/dt y(1)
dy(10) = p.atp_D*y(10)-(p.atp_D*y(10)^2)/p.Dmax...
       - p.aKill_N*y(11)*y(10)...
       - p.aKill_L*y(12)*y(10)...
       - p.aKill_P*y(13)*y(10)...
       - p.aKill_D*y(14)*y(10);

%dBAT_N/dt y(5)
dy(11) = ppval(ND_2,t);

%dBAT_L/dt y(8)
dy(12) = ppval(LD_2,t);

%dBAT_P/dt y(11)
dy(13) = ppval(PD_2,t);

%dBAT_D/dt y(11)
dy(14) = ppval(DD_2,t);

dy(15) = p.aKill_N*y(11)*y(10);
dy(16) = p.aKill_L*y(12)*y(10);
dy(17) = p.aKill_P*y(13)*y(10);
dy(18) = p.aKill_D*y(14)*y(10);


end 