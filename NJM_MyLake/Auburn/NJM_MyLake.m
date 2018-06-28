%Nick Messina MyLake script 
clear all %clears memory space before running
close all
path(path,'C:\Users\Nick\Desktop\NJM_MyLake\air_sea') %path for air-sea toolbox
path(path,'C:\Users\Nick\Desktop\NJM_MyLake\v12') %path for MyLake model code
m_start=[2005,1,1]; %year,month,day
m_stop=[2014,12,31];

initfile='C:\Users\Nick\Desktop\NJM_MyLake\Auburn\NJM_init_v12.xls';%initial conditions
parafile='C:\Users\Nick\Desktop\NJM_MyLake\Auburn\NJM_para_v12.xls';%parameters
inputfile='C:\Users\Nick\Desktop\NJM_MyLake\Auburn\NJM_input_PGS84_00_v12.xls';%input

[zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,Qzt_sed,lambdazt,...
        P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt]...
           = solvemodel_v12(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake'); 
       
figure(1)
plot(Tzt(:,2429),zz);
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('Temperature (C)');
ylabel('Depth (m)');
title('August 26,2011');
xlim([0 25]);

figure(2)
plot(Tzt(:,2432),zz);
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('Temperature (C)');
ylabel('Depth (m)');
title('August 29,2011');
xlim([0 25]);

figure(3)
plot(Pzt(:,2429),zz);
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('Dissolved Inorganic Phosphorus (ppb)');
ylabel('Depth (m)');
title('August 26,2011');
xlim([0 10]);

figure(4)
plot(Pzt(:,2432),zz);
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('Dissolved Inorganic Phosphorus (ppb)');
ylabel('Depth (m)');
title('August 29,2011');
xlim([0 10]);

[obs_temp]=xlsread('C:\Users\Nick\Desktop\NJM_MyLake\Auburn\Observations\temp.xlsx');
figure(5)
plot(Tzt(1:19,3034),zz(1:19),'DisplayName','Modelled');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('Temperature (C)');
ylabel('Depth (m)');
title('April 22,2013')
xlim([0 25]);
hold on
plot(obs_temp(:,1),zz(1:19),'DisplayName','Observed');
legend('show');

figure(6)
plot(Tzt(1:19,3113),zz(1:19),'DisplayName','Modelled');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('Temperature (C)');
ylabel('Depth (m)');
title('July 10,2013')
xlim([0 30]);
hold on
plot(obs_temp(:,2),zz(1:19),'DisplayName','Observed');
legend('show');