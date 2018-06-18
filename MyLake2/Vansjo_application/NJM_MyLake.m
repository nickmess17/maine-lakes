%MyLake script 
clear all
close all
m_start=[2007,5,1]; 
m_stop=[2013,12,31];

%initfile='C:\Users\Nick\Desktop\MYLake\LA_init_v12.xls';
%parafile='C:\Users\Nick\Desktop\MYLake\LA_para_v12.xls';
%inputfile='C:\Users\Nick\Desktop\MYLake\LA_input_2005_2014_v12.xls';
initfile='C:\Users\Nick\Desktop\MyLake2\Vansjo_application\NJM_init_v12.xls';
parafile='C:\Users\Nick\Desktop\MyLake2\Vansjo_application\NJM_para_v12.xls';
inputfile='C:\Users\Nick\Desktop\MyLake2\Vansjo_application\NJM_input_PGS84_00_v12.xls';

[zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,Qzt_sed,lambdazt,...
        P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt]...
           = solvemodel_v12(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake'); 
figure(2)
clf
plot([0 25],[0 25],':b', TempObs(:,3), TempMod, '.r');
set(gca,'fontsize',9);
ylabel('Modelled temperature (^oC)');
xlabel('Measured temperature (^oC)');
axis([0 25 0 25]);
axis square;
title([lake ' ' datestr(datenum(m_start),28) '--' datestr(datenum(m_stop),28)]); 
grid on;