m_start=[1994,5,1]; 
m_stop=[2000,12,31];

initfile='C:\Users\Nick\Desktop\MyLake2\Vansjo_application\NJM_init_v12.xls';
parafile='C:\Users\Nick\Desktop\MyLake2\Vansjo_application\NJM_para_v12.xls';
inputfile='C:\Users\Nick\Desktop\MyLake2\Vansjo_application\NJM_input_PGS84_00_v12.xls';

[zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,Qzt_sed,lambdazt,...
        P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt]...
           = solvemodel_v12(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake'); 