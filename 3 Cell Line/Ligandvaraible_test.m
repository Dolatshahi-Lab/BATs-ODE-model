% clear
expdeath=table2array(readtable('Mar24_3Cancer_tumor.xlsx','Sheet','Total_counts','VariableNamingRule','preserve'));
exp_death=expdeath(1:386,:);
exp_death_fit=exp_death([1,17,49,97,145,193,241,289,337,386],:);
%Phenotype data
cell_prog_time = [0 4 24 48 72 96];

%MCF7 splines
BATs_numbers=table2array(readtable('BATs_numbers_C.xlsx','Sheet','MCF7','VariableNamingRule','preserve'));

            MNBAT_fit = BATs_numbers(:,2);
            MLBAT_fit = BATs_numbers(:,4);
            MPBAT_fit = BATs_numbers(:,3);
            MDBAT_fit = BATs_numbers(:,5);
            N1 = csaps(cell_prog_time,MNBAT_fit,.05);
            L1 = csaps(cell_prog_time,MLBAT_fit,0.01,[],[10 10 1 10 10 0]);
            P1 = csaps(cell_prog_time,MPBAT_fit,.05);
            D1 = csaps(cell_prog_time,MDBAT_fit,.05);
            MN2 = fnder(N1);
            ML2 = fnder(L1);
            MP2 = fnder(P1);
            MD2 = fnder(D1);

%CAMA1 splines
BATs_numbers=table2array(readtable('BATs_numbers_C.xlsx','Sheet','CAMA1','VariableNamingRule','preserve'));

            CNBAT_fit = BATs_numbers(:,2);
            CLBAT_fit = BATs_numbers(:,4);
            CPBAT_fit = BATs_numbers(:,3);
            CDBAT_fit = BATs_numbers(:,5);
            N1 = csaps(cell_prog_time,CNBAT_fit,.05);
            L1 = csaps(cell_prog_time,CLBAT_fit,.05);
            P1 = csaps(cell_prog_time,CPBAT_fit,.05);
            D1 = csaps(cell_prog_time,CDBAT_fit,.05);
            CN2 = fnder(N1);
            CL2 = fnder(L1);
            CP2 = fnder(P1);
            CD2 = fnder(D1);

%T47D splines
BATs_numbers=table2array(readtable('BATs_numbers_C.xlsx','Sheet','T47D','VariableNamingRule','preserve'));

            TNBAT_fit = BATs_numbers(:,2);
            TLBAT_fit = BATs_numbers(:,4);
            TPBAT_fit = BATs_numbers(:,3);
            TDBAT_fit = BATs_numbers(:,5);
            N1 = csaps(cell_prog_time,TNBAT_fit,.05);
            L1 = csaps(cell_prog_time,TLBAT_fit,.05);
            P1 = csaps(cell_prog_time,TPBAT_fit,.05);
            D1 = csaps(cell_prog_time,TDBAT_fit,.05);
            TN2 = fnder(N1);
            TL2 = fnder(L1);
            TP2 = fnder(P1);
            TD2 = fnder(D1);

modelX=0:1:96;
MCF7=passedtotalModelRuns(:,:,1);
CAMA1=passedtotalModelRuns(:,:,2);
T47D=passedtotalModelRuns(:,:,3);

MCF7err=MCF7([1,5,13,25,37,49,61,73,85,97],:);
MCF7err=(MCF7err-exp_death_fit(:,2)).^2;
MCF7err=transpose(sum(MCF7err,1));
MCF7prop=MCF7err/max(MCF7err);
[~,Merr]=min(MCF7err(:,1));

CAMA1err=CAMA1([1,5,13,25,37,49,61,73,85,97],:);
CAMA1err=(CAMA1err-exp_death_fit(:,3)).^2;
CAMA1err=transpose(sum(CAMA1err,1));
CAMA1prop=CAMA1err/max(CAMA1err);
[~,Cerr]=min(CAMA1err(:,1));

T47Derr=T47D([1,5,13,25,37,49,61,73,85,97],:);
T47Derr=(T47Derr-exp_death_fit(:,4)).^2;
T47Derr=transpose(sum(T47Derr,1));
T47Dprop=T47Derr/max(T47Derr);
[~,Terr]=min(T47Derr(:,1));

error=MCF7err+CAMA1err+T47Derr;
[~,minerror]=min(error(:,1));

error=MCF7prop+CAMA1prop+T47Dprop;
[~,minprop]=min(error(:,1));

%Tumor growth parameters
p.Matp = 0.1287;%doubling time
p.Mmax = 1.4;%max tumor cells
p.Catp = 0.06;%doubling time
p.Cmax = 3.2;%max tumor cells
p.Tatp = 0.16;%doubling time
p.Tmax = 1.2;%max tumor cells

n=964;
% n=6370; Donor 3 optimal
n=minerror;
% n=Merr;
% n=minprop;
% BAT killing rates
p.aKill_N = passedtotalParams(n,1);%0.9;
p.aKill_P = passedtotalParams(n,2);%0.9;
p.aKill_L = passedtotalParams(n,3);%0.9;
p.aKill_D = passedtotalParams(n,4);%0.9;
% p.aKill_N = 2.3280e-07;%passedtotalParams_1(n,1);%0.9;
% p.aKill_P = 3.4691e-07;%passedtotalParams_1(n,2);%0.9;
% p.aKill_L = 1.1684e-06;%passedtotalParams_1(n,3);%0.9;
% p.aKill_D = 6.8321e-07;%passedtotalParams_1(n,4);%0.9;

%Setup
MnT=exp_death(1,2);%#Tumor cells
CnT=exp_death(1,3);%#Tumor cells
TnT=exp_death(1,4);%#Tumor cells
tspan = [0 96];

% %2:1 E:T
% ETR=0.45;%E:T ratio
% dose=89804;%nT*ETR;

%Natural
history_2 = [MnT;MNBAT_fit(1);MLBAT_fit(1);MPBAT_fit(1);MDBAT_fit(1);0;0;0;0;CnT;CNBAT_fit(1);CLBAT_fit(1);CPBAT_fit(1);CDBAT_fit(1);0;0;0;0;TnT;TNBAT_fit(1);TLBAT_fit(1);TPBAT_fit(1);TDBAT_fit(1);0;0;0;0];
sol_2 = ode15s(@(t,y) ddefun_MA_nogen_ligandvariable_Jul_nodel(t,y,p,MN2,ML2,MP2,MD2,CN2,CL2,CP2,CD2,TN2,TL2,TP2,TD2), tspan, history_2);


%Counts - MCF7
figure(1)
plot(sol_2.x,sol_2.y(1,:),'linewidth',2)
hold on
plot(exp_death(:,1),exp_death(:,2),'--','linewidth',2,Color=[0.4660 0.6740 0.1880])
ylabel('MCF7 Cell count [Conductivity]')
legend('Tumor','Tumor - Exp')
title("Tumor Cell Counts")
%ylim([0 5e4])
xlim([0 72])
hold off
xlabel('Time [h]')
hold off

figure(2)
plot(sol_2.x,sol_2.y(2,:),'linewidth',2,'color',[0 0.4470 0.7410])
hold on
plot(cell_prog_time,MNBAT_fit,'--o','linewidth',2,'color',[0 0.4470 0.7410])
ylabel('MCF7 Cell count [Total Cells]')
legend('NBATs Total','N - EXP')
title("NBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')

figure(3)
plot(sol_2.x,sol_2.y(3,:),'linewidth',2,'color',[0.85 0.3250 0.0980])
hold on
plot(cell_prog_time,MLBAT_fit,'--o','linewidth',2,'color',[0.85 0.3250 0.0980])
ylabel('MCF7 Cell count [Total Cells]')
legend('LBAT total','L - EXP')
title("LBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')


figure(4)
plot(sol_2.x,sol_2.y(4,:),'linewidth',2,'color',[0.9290 0.6940 0.1250])
hold on
plot(cell_prog_time,MPBAT_fit,'--o','linewidth',2,'color',[0.9290 0.6940 0.1250])
ylabel('MCF7 Cell count [Total Cells]')
legend('TBAT total','T - EXP')
title("TBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')
% 
figure(5)
plot(sol_2.x,sol_2.y(5,:),'linewidth',2,'color',[0.4940 0.1840 0.5560])
hold on
plot(cell_prog_time,MDBAT_fit,'--o','linewidth',2,'color',[0.4940 0.1840 0.5560])
ylabel('MCF7 Cell count [Total Cells]')
legend('DBAT total','D - EXP')
title("DBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')

figure(6)
plot(sol_2.x,sol_2.y([6,7,8,9],:),'linewidth',2)
hold on
ylabel('MCF7 Cell count [Total Cells]')
legend('NBATs','LBATs','TBATs','DBATs')
title("Total cell deaths by population")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')

%Counts - CAMA1
figure(7)
plot(sol_2.x,sol_2.y(10,:),'linewidth',2)
hold on
plot(exp_death(:,1),exp_death(:,3),'--','linewidth',2,Color=[0.4660 0.6740 0.1880])
ylabel('CAMA1 Cell count [Conductivity]')
legend('Tumor','Tumor - Exp')
title("Tumor Cell Counts")
%ylim([0 5e4])
xlim([0 72])
hold off
xlabel('Time [h]')
hold off

figure(8)
plot(sol_2.x,sol_2.y(11,:),'linewidth',2,'color',[0 0.4470 0.7410])
hold on
plot(cell_prog_time,CNBAT_fit,'--o','linewidth',2,'color',[0 0.4470 0.7410])
ylabel('CAMA1 Cell count [Total Cells]')
legend('NBATs Total','N - EXP')
title("NBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')

figure(9)
plot(sol_2.x,sol_2.y(12,:),'linewidth',2,'color',[0.85 0.3250 0.0980])
hold on
plot(cell_prog_time,CLBAT_fit,'--o','linewidth',2,'color',[0.85 0.3250 0.0980])
ylabel('CAMA1 Cell count [Total Cells]')
legend('LBAT total','L - EXP')
title("LBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')


figure(10)
plot(sol_2.x,sol_2.y(13,:),'linewidth',2,'color',[0.9290 0.6940 0.1250])
hold on
plot(cell_prog_time,CPBAT_fit,'--o','linewidth',2,'color',[0.9290 0.6940 0.1250])
ylabel('CAMA1 Cell count [Total Cells]')
legend('TBAT total','T - EXP')
title("TBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')
% 
figure(11)
plot(sol_2.x,sol_2.y(14,:),'linewidth',2,'color',[0.4940 0.1840 0.5560])
hold on
plot(cell_prog_time,CDBAT_fit,'--o','linewidth',2,'color',[0.4940 0.1840 0.5560])
ylabel('CAMA1 Cell count [Total Cells]')
legend('DBAT total','D - EXP')
title("DBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')

figure(12)
plot(sol_2.x,sol_2.y([15,16,17,18],:),'linewidth',2)
hold on
ylabel('CAMA1 Cell count [Total Cells]')
legend('NBATs','LBATs','TBATs','DBATs')
title("Total cell deaths by population")
%ylim([0 5e4])
xlim([0 72])
hold off
xlabel('Time [h]')


%Counts - T47D
figure(13)
plot(sol_2.x,sol_2.y(19,:),'linewidth',2)
hold on
plot(exp_death(:,1),exp_death(:,4),'--','linewidth',2,Color=[0.4660 0.6740 0.1880])
ylabel('T47D Cell count [Conductivity]')
legend('Tumor','Tumor - Exp')
title("Tumor Cell Counts")
%ylim([0 5e4])
xlim([0 72])
hold off
xlabel('Time [h]')
hold off

figure(14)
plot(sol_2.x,sol_2.y(20,:),'linewidth',2,'color',[0 0.4470 0.7410])
hold on
plot(cell_prog_time,TNBAT_fit,'--o','linewidth',2,'color',[0 0.4470 0.7410])
ylabel('T47D Cell count [Total Cells]')
legend('NBATs Total','N - EXP')
title("NBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')

figure(15)
plot(sol_2.x,sol_2.y(21,:),'linewidth',2,'color',[0.85 0.3250 0.0980])
hold on
plot(cell_prog_time,TLBAT_fit,'--o','linewidth',2,'color',[0.85 0.3250 0.0980])
ylabel('T47D Cell count [Total Cells]')
legend('LBAT total','L - EXP')
title("LBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')


figure(16)
plot(sol_2.x,sol_2.y(22,:),'linewidth',2,'color',[0.9290 0.6940 0.1250])
hold on
plot(cell_prog_time,TPBAT_fit,'--o','linewidth',2,'color',[0.9290 0.6940 0.1250])
ylabel('T47D Cell count [Total Cells]')
legend('TBAT total','T - EXP')
title("TBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')
% 
figure(17)
plot(sol_2.x,sol_2.y(23,:),'linewidth',2,'color',[0.4940 0.1840 0.5560])
hold on
plot(cell_prog_time,TDBAT_fit,'--o','linewidth',2,'color',[0.4940 0.1840 0.5560])
ylabel('T47D Cell count [Total Cells]')
legend('DBAT total','D - EXP')
title("DBAT Cell Counts")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')

figure(18)
plot(sol_2.x,sol_2.y([24,25,26,27],:),'linewidth',2)
hold on
ylabel('T47D Cell count [Total Cells]')
legend('NBATs','LBATs','TBATs','DBATs')
title("Total cell deaths by population")
xlim([0 72])
%ylim([0 5e4])
hold off
xlabel('Time [h]')


% 
% 
Max = max(passedtotalParams);
Min = min(passedtotalParams);
% 
% 
Max_T = max(transpose(MCF7));
Min_T = min(transpose(MCF7));

figure(22)
plot(exp_death(:,1),exp_death(:,2),'-.o','linewidth',2,Color=[0 0 0])
hold on
plot(sol_2.x,sol_2.y(1,:),'linewidth',2,Color=[0 0.4470 0.7410])
plot(modelX,Min_T,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Max_T,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
ylabel('Cell count [Conductivity]')
legend('Tumor - Exp','Tumor - best','Tumor - range')
title("MCF7 Tumor Cell Counts")
%ylim([0 5e4])
xlim([0 72])
hold off
xlabel('Time [h]')

Max_T = max(transpose(CAMA1));
Min_T = min(transpose(CAMA1));

figure(23)
plot(exp_death(:,1),exp_death(:,3),'-.o','linewidth',2,Color=[0 0 0])
hold on
plot(sol_2.x,sol_2.y(10,:),'linewidth',2,Color=[0 0.4470 0.7410])
plot(modelX,Min_T,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Max_T,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
ylabel('Cell count [Conductivity]')
legend('Tumor - Exp','Tumor - best','Tumor - range')
title("CAMA1 Tumor Cell Counts")
%ylim([0 5e4])
xlim([0 72])
hold off
xlabel('Time [h]')

Max_T = max(transpose(T47D));
Min_T = min(transpose(T47D));

figure(24)
plot(exp_death(:,1),exp_death(:,4),'-.o','linewidth',2,Color=[0 0 0])
hold on
plot(sol_2.x,sol_2.y(19,:),'linewidth',2,Color=[0 0.4470 0.7410])
plot(modelX,Min_T,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Max_T,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
ylabel('Cell count [Conductivity]')
legend('Tumor - Exp','Tumor - best','Tumor - range')
title("T47D Tumor Cell Counts")
%ylim([0 5e4])
xlim([0 72])
hold off
xlabel('Time [h]')
% 
% 
% 
% 
% 
% MCF7
% for i = 1:length(passedtotalParams)
% for i = 1:3320
% n=i;
% % 
% p.aKill_N = passedtotalParams(n,1);%0.9;
% p.aKill_P = passedtotalParams(n,2);%0.9;
% p.aKill_L = passedtotalParams(n,3);%0.9;
% p.aKill_D = passedtotalParams(n,4);%0.9;
% 
% %Setup
% 
% 
% %Natural
% history_2 = [MnT;MNBAT_fit(1);MLBAT_fit(1);MPBAT_fit(1);MDBAT_fit(1);0;0;0;0;CnT;CNBAT_fit(1);CLBAT_fit(1);CPBAT_fit(1);CDBAT_fit(1);0;0;0;0;TnT;TNBAT_fit(1);TLBAT_fit(1);TPBAT_fit(1);TDBAT_fit(1);0;0;0;0];
% sol_range = ode(ODEFcn=@(t,y) ddefun_MA_nogen_ligandvariable_Jul_nodel(t,y,p,MN2,ML2,MP2,MD2,CN2,CL2,CP2,CD2,TN2,TL2,TP2,TD2), InitialValue=history_2);
% tspan = linspace(0,96,97);
% S = solve(sol_range,tspan);
% MN_elim(n,:) = S.Solution(6,:);
% ML_elim(n,:) = S.Solution(7,:);
% MP_elim(n,:) = S.Solution(8,:);
% MD_elim(n,:) = S.Solution(9,:);
% CN_elim(n,:) = S.Solution(15,:);
% CL_elim(n,:) = S.Solution(16,:);
% CP_elim(n,:) = S.Solution(17,:);
% CD_elim(n,:) = S.Solution(18,:);
% TN_elim(n,:) = S.Solution(24,:);
% TL_elim(n,:) = S.Solution(25,:);
% TP_elim(n,:) = S.Solution(26,:);
% TD_elim(n,:) = S.Solution(27,:);
% end
% 
%MCF7 range 
Max_Nelim = max((MN_elim));
Min_Nelim = min((MN_elim));

Max_Lelim = max((ML_elim));
Min_Lelim = min((ML_elim));

Max_Pelim = max((MP_elim));
Min_Pelim = min((MP_elim));

Max_Delim = max((MD_elim));
Min_Delim = min((MD_elim));

figure(25)
plot(sol_2.x,sol_2.y([6,7,8,9],:),'linewidth',2)
hold on
plot(modelX,Min_Nelim,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Min_Lelim,'linewidth',2,'LineStyle',':',Color=[0.85 0.3250 0.0980])
plot(modelX,Min_Pelim,'linewidth',2,'LineStyle',':',Color=[0.9290 0.6940 0.1250])
plot(modelX,Min_Delim,'linewidth',2,'LineStyle',':',Color=[0.4940 0.1840 0.5560])
plot(modelX,Max_Nelim,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Max_Lelim,'linewidth',2,'LineStyle',':',Color=[0.85 0.3250 0.0980])
plot(modelX,Max_Pelim,'linewidth',2,'LineStyle',':',Color=[0.9290 0.6940 0.1250])
plot(modelX,Max_Delim,'linewidth',2,'LineStyle',':',Color=[0.4940 0.1840 0.5560])
ylabel('Cell count [Conductivity]')
legend('NBATs','LBATs','TBATs','DBATs')
title("MCF7 Total cell deaths by population")
xlim([0 72])
hold off
xlabel('Time [h]')


%CAMA1 range 
Max_Nelim = max((CN_elim));
Min_Nelim = min((CN_elim));

Max_Lelim = max((CL_elim));
Min_Lelim = min((CL_elim));

Max_Pelim = max((CP_elim));
Min_Pelim = min((CP_elim));

Max_Delim = max((CD_elim));
Min_Delim = min((CD_elim));

figure(26)
plot(sol_2.x,sol_2.y([15,16,17,18],:),'linewidth',2)
hold on
plot(modelX,Min_Nelim,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Min_Lelim,'linewidth',2,'LineStyle',':',Color=[0.85 0.3250 0.0980])
plot(modelX,Min_Pelim,'linewidth',2,'LineStyle',':',Color=[0.9290 0.6940 0.1250])
plot(modelX,Min_Delim,'linewidth',2,'LineStyle',':',Color=[0.4940 0.1840 0.5560])
plot(modelX,Max_Nelim,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Max_Lelim,'linewidth',2,'LineStyle',':',Color=[0.85 0.3250 0.0980])
plot(modelX,Max_Pelim,'linewidth',2,'LineStyle',':',Color=[0.9290 0.6940 0.1250])
plot(modelX,Max_Delim,'linewidth',2,'LineStyle',':',Color=[0.4940 0.1840 0.5560])
ylabel('Cell count [Conductivity]')
legend('NBATs','LBATs','TBATs','DBATs')
title("CAMA1 Total cell deaths by population")
xlim([0 72])
hold off
xlabel('Time [h]')


%T47D range 
Max_Nelim = max((TN_elim));
Min_Nelim = min((TN_elim));

Max_Lelim = max((TL_elim));
Min_Lelim = min((TL_elim));

Max_Pelim = max((TP_elim));
Min_Pelim = min((TP_elim));

Max_Delim = max((TD_elim));
Min_Delim = min((TD_elim));

figure(21)
plot(sol_2.x,sol_2.y([24,25,26,27],:),'linewidth',2)
hold on
plot(modelX,Min_Nelim,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Min_Lelim,'linewidth',2,'LineStyle',':',Color=[0.85 0.3250 0.0980])
plot(modelX,Min_Pelim,'linewidth',2,'LineStyle',':',Color=[0.9290 0.6940 0.1250])
plot(modelX,Min_Delim,'linewidth',2,'LineStyle',':',Color=[0.4940 0.1840 0.5560])
plot(modelX,Max_Nelim,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Max_Lelim,'linewidth',2,'LineStyle',':',Color=[0.85 0.3250 0.0980])
plot(modelX,Max_Pelim,'linewidth',2,'LineStyle',':',Color=[0.9290 0.6940 0.1250])
plot(modelX,Max_Delim,'linewidth',2,'LineStyle',':',Color=[0.4940 0.1840 0.5560])
ylabel('Cell count [Conductivity]')
legend('NBATs','LBATs','TBATs','DBATs')
title("T47D Total cell deaths by population")
xlim([0 72])
hold off
xlabel('Time [h]')
% 
% % 
passCriteria=2;

LogicalArrayCondition1 = modelOutput(:,:) >= expData(:,1) / passCriteria;
LogicalArrayCondition2 = modelOutput(:,:) <= expData(:,2) * passCriteria;
% LogicalArrayCondition3 = modelOutput(2,:) <= expData(2,2) * 1.5;
% LogicalArrayCondition4 = modelOutput(2,:) >= expData(2,1) / 1.5;
% LogicalArrayCondition5 = modelOutput(9,:) <= expData(9,2) * 1.5;
% LogicalArrayCondition6 = modelOutput(9,:) >= expData(9,1) / 1.25;
% LogicalArraytot = LogicalArrayCondition1 & LogicalArrayCondition2 & LogicalArrayCondition3 & LogicalArrayCondition4 & LogicalArrayCondition5 & LogicalArrayCondition6;
LogicalArraytot = LogicalArrayCondition1 & LogicalArrayCondition2;% & LogicalArrayCondition5;

critsum=sum(LogicalArraytot,2);


min(passedtotalParams(:,1))
min(passedtotalParams(:,2))
min(passedtotalParams(:,3))
min(passedtotalParams(:,4))

max(passedtotalParams(:,1))
max(passedtotalParams(:,2))
max(passedtotalParams(:,3))
max(passedtotalParams(:,4))
