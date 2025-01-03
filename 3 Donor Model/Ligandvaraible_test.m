% clear
expdeath_D=table2array(readtable('Mar24_3Donor_tumor.xlsx','Sheet','Donor 1','VariableNamingRule','preserve'));
exp_death_D=expdeath_D(1:583,:);
exp_death_fit_D=exp_death_D([1,17,49,97,145,193,241,289,337,386],:);
expdeath_C=expdeath_D;
exp_death_C=exp_death_D;
exp_death_fit_C=exp_death_fit_D;
cell_prog_time = [0 4 24 48 72 96];

%Phenotype data
cell_prog_time = [0 4 24 48 72 96];

BATs_numbers=table2array(readtable('BATs_numbers_D.xlsx','Sheet','Donor_1','VariableNamingRule','preserve'));

CNBAT_fit = BATs_numbers(:,2);
CLBAT_fit = BATs_numbers(:,4);
CPBAT_fit = BATs_numbers(:,3);
CDBAT_fit = BATs_numbers(:,5);

N1 = csaps(cell_prog_time,CNBAT_fit,0.2,[],[1 0.04 .1 10 10 10]);
L1 = csaps(cell_prog_time,CLBAT_fit,0.05);
P1 = csaps(cell_prog_time,CPBAT_fit,0.2,[],[1 0.04 .1 10 10 10]);
D1 = csaps(cell_prog_time,CDBAT_fit,0.05);

NC_2 = fnder(N1);
LC_2 = fnder(L1);
PC_2 = fnder(P1);
DC_2 = fnder(D1);

BATs_numbers=table2array(readtable('BATs_numbers_D.xlsx','Sheet','Donor_1','VariableNamingRule','preserve'));

DNBAT_fit = BATs_numbers(:,2);
DLBAT_fit = BATs_numbers(:,4);
DPBAT_fit = BATs_numbers(:,3);
DDBAT_fit = BATs_numbers(:,5);

N1 = csaps(cell_prog_time,CNBAT_fit,0.2,[],[1 0.04 .1 10 10 10]);
L1 = csaps(cell_prog_time,DLBAT_fit,0.05);
P1 = csaps(cell_prog_time,DPBAT_fit,0.05);
D1 = csaps(cell_prog_time,DDBAT_fit,0.05);

ND_2 = fnder(N1);
LD_2 = fnder(L1);
PD_2 = fnder(P1);
DD_2 = fnder(D1);


modelX=0:1:96;%transpose(CalibrationOutputArray{1, end}.odeSettings.tspan);

TumorC=passedtotalModelRuns(:,:,1);
TumorD=passedtotalModelRuns(:,:,2);

tumorerr=TumorC([1,5,13,25,37,49,61,73,85,97],:);
tumorerr=(tumorerr-exp_death_fit_C(:,3)).^2;
% tumorerr=(tumorerr-cyto_fit_C(:,3)).^2;
tumorerrC=transpose(sum(tumorerr,1));
tumorprop=tumorerrC/max(tumorerrC);
[~,minerrorC]=min(tumorerrC(:,1));

tumorerr=TumorD([1,5,13,25,37,49,61,73,85,97],:);
tumorerr=(tumorerr-exp_death_fit_D(:,3)).^2;
% tumorerr=(tumorerr-cyto_fit_D(:,3)).^2;
tumorerrD=transpose(sum(tumorerr,1));
tumorprop=tumorerrD/max(tumorerrD);
[~,minerrorD]=min(tumorerrD(:,1));

totalerror=tumorerrC+tumorerrD;
[~,minerror]=min(totalerror(:,1));
[~,maxerror]=max(totalerror(:,1));

%Tumor growth parameters
p.atp_C = 0.0748;%doubling time
p.atp_D = 0.0748;%doubling time
% p.atp_C = 0.1694;%doubling time
% p.atp_D = 0.0748;%doubling time
p.Cmax = 7.6;%max tumor cells
p.Dmax = 7.6;%max tumor cells

n=1;
n=minerrorD;
% n=minerror;
% n=maxerror;
% n=minlagerror;
% n=minprop;
% BAT killing rates
p.aKill_N =passedtotalParams(n,1);%0.9;
p.aKill_P =passedtotalParams(n,2);%0.9;
p.aKill_L =passedtotalParams(n,3);%0.9;
p.aKill_D =passedtotalParams(n,4);%0.9;
% p.aKill_N = 1.23e-8;%passedtotalParams_1(n,1);%0.9;
% p.aKill_P = 1.85e-8;%passedtotalParams_1(n,2);%0.9;
% p.aKill_L = 7.26e-7*1.1;%passedtotalParams_1(n,3);%0.9;
% p.aKill_D = 6.54e-7;%passedtotalParams_1(n,4);%0.9;
% 
% %Setup
CnT=exp_death_C(1,3);%#Tumor cells
DnT=exp_death_D(1,3);%#Tumor cells
tspan = [0 96];

%2:1 E:T
ETR=0.45;%E:T ratio
dose=89804;%nT*ETR;

%Natural
history_2 = [CnT;CNBAT_fit(1);CLBAT_fit(1);CPBAT_fit(1);CDBAT_fit(1);0;0;0;0;DnT;DNBAT_fit(1);DLBAT_fit(1);DPBAT_fit(1);DDBAT_fit(1);0;0;0;0];
sol_2 = ode15s(@(t,y) ddefun_MA_nogen_ligandvariable_Jul_nodel(t,y,p,NC_2,LC_2,PC_2,DC_2,ND_2,LD_2,PD_2,DD_2), tspan, history_2);

% cytotox_C = (1-(sol_2.y(1,:)./sol_2.y(19,:)))*100;
% cytotox_D = (1-(sol_2.y(10,:)./sol_2.y(20,:)))*100;

%Counts
% figure(1)
% plot(sol_2.x,sol_2.y(1,:),'linewidth',2)
% hold on
% plot(exp_death_C(:,1),exp_death_C(:,3),'--','linewidth',2,Color=[0.4660 0.6740 0.1880])
% ylabel('Cell count [Total Cells]')
% legend('Tumor','Tumor CTRL','Tumor - Exp')
% title("Tumor Cell Counts")
% %ylim([0 5e4])
% xlim([0 100])
% hold off
% xlabel('Time [h]')
% hold off
% 
figure(2)
plot(sol_2.x,sol_2.y(2,:),'linewidth',2,'color',[0 0.4470 0.7410])
hold on
plot(cell_prog_time,CNBAT_fit,'--o','linewidth',2,'color',[0 0.4470 0.7410])
ylabel('Cell count [Total Cells]')
legend('NBATs Total','N - EXP')
title("NBAT Cell Counts")
%ylim([0 5e4])
hold off
xlabel('Time [h]')
% 
% figure(3)
% plot(sol_2.x,sol_2.y(3,:),'linewidth',2,'color',[0.85 0.3250 0.0980])
% hold on
% plot(cell_prog_time,CLBAT_fit,'--o','linewidth',2,'color',[0.85 0.3250 0.0980])
% ylabel('Cell count [Total Cells]')
% legend('LBAT total','L - EXP')
% title("LBAT Cell Counts")
% %ylim([0 5e4])
% hold off
% xlabel('Time [h]')
% 
% 
figure(4)
plot(sol_2.x,sol_2.y(4,:),'linewidth',2,'color',[0.9290 0.6940 0.1250])
hold on
plot(cell_prog_time,CPBAT_fit,'--o','linewidth',2,'color',[0.9290 0.6940 0.1250])
ylabel('Cell count [Total Cells]')
legend('PBAT total','P - EXP')
title("PBAT Cell Counts")
%ylim([0 5e4])
hold off
xlabel('Time [h]')
% 
% figure(5)
% plot(sol_2.x,sol_2.y(5,:),'linewidth',2,'color',[0.4940 0.1840 0.5560])
% hold on
% plot(cell_prog_time,CDBAT_fit,'--o','linewidth',2,'color',[0.4940 0.1840 0.5560])
% ylabel('Cell count [Total Cells]')
% legend('DBAT total','D - EXP')
% title("DBAT Cell Counts")
% %ylim([0 5e4])
% hold off
% xlabel('Time [h]')
% 
% 
figure(6)
plot(sol_2.x,sol_2.y([6,7,8,9],:),'linewidth',2)
hold on
ylabel('Cell count [Total Cells]')
legend('NBATs','LBATs','TBATs','DBATs')
title("Donor 1: Total cell deaths by population")
%ylim([0 5e4])
hold off
xlabel('Time [h]')


figure(7)
plot(sol_2.x,sol_2.y(10,:),'linewidth',2)
hold on
plot(exp_death_D(:,1),exp_death_D(:,3),'--','linewidth',2,Color=[0.4660 0.6740 0.1880])
ylabel('Cell count [Total Cells]')
legend('Tumor','Tumor CTRL','Tumor - Exp')
title("Tumor Cell Counts")
%ylim([0 5e4])
xlim([0 100])
hold off
xlabel('Time [h]')
hold off

figure(8)
plot(sol_2.x,sol_2.y(11,:),'linewidth',2,'color',[0 0.4470 0.7410])
hold on
plot(cell_prog_time,DNBAT_fit,'--o','linewidth',2,'color',[0 0.4470 0.7410])
ylabel('Cell count [Total Cells]')
legend('NBATs Total','N - EXP')
title("NBAT Cell Counts")
%ylim([0 5e4])
hold off
xlabel('Time [h]')

% figure(9)
% plot(sol_2.x,sol_2.y(12,:),'linewidth',2,'color',[0.85 0.3250 0.0980])
% hold on
% plot(cell_prog_time,DLBAT_fit,'--o','linewidth',2,'color',[0.85 0.3250 0.0980])
% ylabel('Cell count [Total Cells]')
% legend('LBAT total','L - EXP')
% title("LBAT Cell Counts")
% %ylim([0 5e4])
% hold off
% xlabel('Time [h]')
% 
% 
figure(10)
plot(sol_2.x,sol_2.y(13,:),'linewidth',2,'color',[0.9290 0.6940 0.1250])
hold on
plot(cell_prog_time,DPBAT_fit,'--o','linewidth',2,'color',[0.9290 0.6940 0.1250])
ylabel('Cell count [Total Cells]')
legend('PBAT total','P - EXP')
title("PBAT Cell Counts")
%ylim([0 5e4])
hold off
xlabel('Time [h]')
% % 
% figure(11)
% plot(sol_2.x,sol_2.y(14,:),'linewidth',2,'color',[0.4940 0.1840 0.5560])
% hold on
% plot(cell_prog_time,DDBAT_fit,'--o','linewidth',2,'color',[0.4940 0.1840 0.5560])
% ylabel('Cell count [Total Cells]')
% legend('DBAT total','D - EXP')
% title("DBAT Cell Counts")
% %ylim([0 5e4])
% hold off
% xlabel('Time [h]')


figure(12)
plot(sol_2.x,sol_2.y([15,16,17,18],:),'linewidth',2)
hold on
ylabel('Cell count [Total Cells]')
legend('NBATs','LBATs','TBATs','DBATs')
title("Donor 1: Total cell deaths by population")
xlim([0 72])
hold off
xlabel('Time [h]')

% 
% figure(15)
% plot(sol_2.x,cytotox_C,'--o','linewidth',2,Color=[0.4660 0.6740 0.1880])
% hold on
% plot(cyto_C(:,1),cyto_C(:,3),'--','linewidth',2,Color=[0.4660 0.6740 0.1880])
% plot(sol_2.x,cytotox_D,'--o','linewidth',2,Color=[0.4940 0.1840 0.5560])
% plot(cyto_D(:,1),cyto_D(:,3),'--','linewidth',2,Color=[0.4940 0.1840 0.5560])
% ylabel('Cell count [Total Cells]')
% legend('Cyto model','Cyto - Exp')
% title("Tumor Cell Counts")
% %ylim([0 5e4])
% xlim([0 100])
% hold off
% xlabel('Time [h]')
% hold off



Max = max(passedtotalParams);
Min = min(passedtotalParams);


Max_TC = max(transpose(TumorC));
Min_TC = min(transpose(TumorC));



figure(13)
plot(exp_death_C(:,1),exp_death_C(:,3),'-.o','linewidth',2,Color=[0 0 0])
hold on
plot(sol_2.x,sol_2.y(1,:),'linewidth',2,Color=[0 0.4470 0.7410])
plot(modelX,Min_TC,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Max_TC,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
ylabel('Cell count [Total Cells]')
legend('Tumor - Exp','Tumor - best','Tumor - range')
title("Donor 1 Tumor Cell Counts")
%ylim([0 5e4])
xlim([0 72])
hold off
xlabel('Time [h]')

Max_TD = max(transpose(TumorD));
Min_TD = min(transpose(TumorD));

figure(14)
plot(exp_death_D(:,1),exp_death_D(:,3),'-.o','linewidth',2,Color=[0 0 0])
hold on
plot(sol_2.x,sol_2.y(10,:),'linewidth',2,Color=[0 0.4470 0.7410])
plot(modelX,Min_TD,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
plot(modelX,Max_TD,'linewidth',2,'LineStyle',':',Color=[0 0.4470 0.7410])
ylabel('Cell count [Total Cells]')
legend('Tumor - Exp','Tumor - best','Tumor - range')
title("Donor 1 Tumor Cell Counts")
%ylim([0 5e4])
xlim([0 100])
hold off
xlabel('Time [h]')


