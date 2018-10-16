%Nick Messina MyLake script 
function [SS_temp]=NJM_MyLake(par)%clc, clear, close all%clears figures and data before running
path(path,'C:\Users\Nick\Desktop\NJM_MyLake\air_sea') %path for air-sea toolbox
path(path,'C:\Users\Nick\Desktop\NJM_MyLake\v12') %path for MyLake model code

lake='Lake Auburn';
year=2006;
m_start=[2006,3,31]; %year,month,day
m_stop=[2014,12,31];

% for i=1:72
% filename=sprintf('data%03d.xlsx',i);
% sheet='lake';
% xlswrite(fullfile('C:\Users\Nick\Desktop\NJM_MyLake\Auburn',filename),sheet);
% xlswrite(filename,par(1),sheet,'B4');%open water diffusion scaling coefficient
% xlswrite(filename,par(2),sheet,'B7');%wind sheltering coefficient
% xlswrite(filename,par(3),sheet,'B18');%scaling of inflow volume
% xlswrite(filename,par(4),sheet,'B19');%scaling of inflow temperature
% xlswrite(filename,par(5),sheet,'B26');%non PAR light attenuation coefficient
% xlswrite(filename,par(6),sheet,'B27');%PAR light attenuation coefficient
% end

filename='C:\Users\Nick\Desktop\NJM_MyLake\Auburn\NJM_para_v12.xls';
sheet='lake';
xlswrite(filename,par(1),sheet,'B4');%open water diffusion scaling coefficient
xlswrite(filename,par(2),sheet,'B7');%wind sheltering coefficient
xlswrite(filename,par(3),sheet,'B18');%scaling of inflow volume
xlswrite(filename,par(4),sheet,'B19');%scaling of inflow temperature
xlswrite(filename,par(5),sheet,'B26');%non PAR light attenuation coefficient
xlswrite(filename,par(6),sheet,'B27');%PAR light attenuation coefficient

initfile='C:\Users\Nick\Desktop\NJM_MyLake\Auburn\NJM_init_v12.xls';%initial profile
parafile=filename;%parameters
% parafile='C:\Users\Nick\Desktop\NJM_MyLake\Auburn\NJM_para_v12.xls';
inputfile='C:\Users\Nick\Desktop\NJM_MyLake\Auburn\NJM_input_PGS84_00_v12.xls';%daily inputs

[zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,Qzt_sed,lambdazt,...
        P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt]...
           = solvemodel_v12(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake'); 

tt_mod=tt-datenum(year,1,1);%time now scaled so that it begins from the 1 january of the "year" (=0)

%Temperature profile observations(tt_mod, z, T)
[obs_temp]=xlsread('C:\Users\Nick\Desktop\NJM_MyLake\Auburn\Observations\MyLake_temp_profiles.xlsx');%year,month,day,depth,temperature
date_of_observation=[obs_temp(:,1),obs_temp(:,2),obs_temp(:,3)];
TempObs=[datenum(date_of_observation) - datenum(year,1,1), obs_temp(:,4), obs_temp(:,5)];

%=align temperature observations with model results
initial=[1;find(diff(TempObs(:,1))~=0)+1];
finish=[find(diff(TempObs(:,1))~=0); length(TempObs)];
for i=1:length(initial)
operation=find(tt_mod==TempObs(initial(i),1));
    if (isempty(operation)==0)
    TempMod(initial(i):finish(i))=interp1(zz,Tzt(:,operation),TempObs(initial(i):finish(i),2));
    else
    TempMod(initial(i):finish(i))=NaN;    
    end    
end

%calculate sum squares for temperature
ss=zeros(1,numel(TempMod));
for i=1:numel(TempMod)
    ss(1,i)=((TempMod(1,i)-TempObs(i,3)).^2);
end
SS_temp=nansum(ss);
disp(SS_temp)

%calculate TP
TPzt=Chlzt+Pzt+PPzt;

%TP profile observations
[obs_TP]=xlsread('C:\Users\Nick\Desktop\NJM_MyLake\Auburn\Observations\MyLake_TP_profiles.xlsx');%year,month,day,depth,temperature
date_of_observation_TP=[obs_TP(:,1),obs_TP(:,2),obs_TP(:,3)];
TPobs=[datenum(date_of_observation_TP) - datenum(year,1,1), obs_TP(:,4), obs_TP(:,5)];

%=align TP observations with model results
initialP=[1;find(diff(TPobs(:,1))~=0)+1];
finishP=[find(diff(TPobs(:,1))~=0); length(TPobs)];
for i=1:length(initialP)
operationP=find(tt_mod==TPobs(initialP(i),1));
    if (isempty(operationP)==0)
    TPmod(initialP(i):finishP(i))=interp1(zz,TPzt(:,operationP),TPobs(initialP(i):finishP(i),2));
    else
    TPmod(initialP(i):finishP(i))=NaN;    
    end    
end

%prepare for plotting
zlim = [0 max(zz)];
tlim = [min(tt_mod) max(tt_mod)];
end

% figure(1)
% clf
% fign=9;
% post_operation=find(isnan(TempMod(initial))==0);
% post_initial=initial(post_operation);
% post_finish=finish(post_operation);
% for i = 1:min(fign,length(post_initial))
%    subplot(3,3,i);
%    operation=find(tt_mod==TempObs(post_initial(i),1));
%    plot(Tzt(:,operation),zz,'-b',TempObs(post_initial(i):post_finish(i),3),TempObs(post_initial(i):post_finish(i),2),'.r');
%    axis([0 30 zlim])
%    axis ij
%    title(['Lake Auburn ' datestr(tt_mod(operation)+datenum(year,1,1)),],'fontsize',8); 
%    set(gca,'FontSize',8) 
% end
% subplot(337)
% ylabel('Depth (m)','fontsize',8)
% xlabel('Temperature (^{o}C)','fontsize',8)
% subplot(331)
% ylabel('Depth (m)','fontsize',8)
% legend('Modelled','Observed');
% subplot(338)
% xlabel('Temperature (^{o}C)','fontsize',8)
% subplot(339)
% xlabel('Temperature (^{o}C)','fontsize',8)
% 
% figure(2)
% clf
% for i = fign+1:min(18,length(post_initial))
%    subplot(3,3,i-fign);
%    set(gca,'FontSize',8) 
%    operation=find(tt_mod==TempObs(post_initial(i),1));
%    plot(Tzt(:,operation),zz,'-b',TempObs(post_initial(i):post_finish(i),3),TempObs(post_initial(i):post_finish(i),2),'.r');
%    axis([0 30 zlim])
%    axis ij
%    title(['Lake Auburn ' datestr(tt_mod(operation)+datenum(year,1,1)),],'fontsize',8); 
% end
%    subplot(337)
%    ylabel('Depth (m)','fontsize',8)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%     subplot(331)
%    ylabel('Depth (m)','fontsize',8)
%    legend('Modelled','Observed');
%     subplot(338)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%      subplot(339)
%    xlabel('Temperature (^{o}C)','fontsize',8)
% 
% figure(3)
% clf
% plot([0 25],[0 25],':b', TempObs(:,3), TempMod, '.r');
% set(gca,'fontsize',9);
% ylabel('Modelled temperature (^oC)');
% xlabel('Measured temperature (^oC)');
% axis([0 25 0 25]);
% axis square;
% title([lake ' ' datestr(datenum(m_start),28) '--' datestr(datenum(m_stop),28)]); 
% grid on;
% mdl=fitlm(TempObs(:,3),TempMod);
% 
%  tlims=[datenum([2013,1,1]) datenum(m_stop)];
% 
% figure(4)
% clf
% subplot(511)
% inx=find(round(TempObs(:,2))<2);
% H=plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),Tzt(1,:),'-');
%  set(gca,'ylim',[0 30]);
%  ylabel('^oC','fontsize',9)
%  title('Temperature  0-1m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
% set(gca,'xlim',tlims)
%  
% subplot(512)
% inx=find((round(TempObs(:,2))==3));
% plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),Tzt(3,:),'-');
%  set(gca,'ylim',[0 30]);
%  ylabel('^oC','fontsize',9)
%  title('Temperature  2-3m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
%  set(gca,'xlim',tlims)
%  
% subplot(513)
% inx=find((round(TempObs(:,2))==5));
% plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),Tzt(5,:),'-');
%  set(gca,'ylim',[0 30]);
%  ylabel('^oC','fontsize',9)
%  title('Temperature  4-5m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
%  set(gca,'xlim',tlims)
%  
% subplot(514)
% inx=find((round(TempObs(:,2))==11));
% plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),Tzt(11,:),'-');
%  set(gca,'ylim',[0 25]);
%  ylabel('^oC','fontsize',9)
%  title('Temperature  10-11m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
%  set(gca,'xlim',tlims)
%  
%  subplot(515)
% inx=find((round(TempObs(:,2))==30)|(round(TempObs(:,2))==31));
% plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),Tzt(31,:),'-');
%  set(gca,'ylim',[0 25]);
%  ylabel('^oC','fontsize',9)
%  xlabel('year','fontsize',9)
%  title('Temperature 30-31m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xlim',tlims)
 
%  figure(5)
%  clf
% subplot(311)
% inxP=find(round(TPobs(:,2))<2);
% H=plot(TPobs(inxP,1)+datenum(year,1,1),TPobs(inxP,3),'r+','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),TPzt(1,:),'-');
%  set(gca,'ylim',[0 55]);
%  ylabel('ppb','fontsize',9)
%  title('TP  0-1m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
% %set(gca,'xlim',tlims)
% 
% subplot(312)
% inxP=find((round(TPobs(:,2))==15));
% plot(TPobs(inxP,1)+datenum(year,1,1),TPobs(inxP,3),'r+','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),TPzt(15,:),'-');
%  set(gca,'ylim',[0 55]);
%  ylabel('ppb','fontsize',9)
%  title('TP 14-15m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
%  %set(gca,'xlim',tlims)
%  
%   subplot(313)
% inxP=find((round(TPobs(:,2))==24)|(round(TPobs(:,2))==25));
% plot(TPobs(inxP,1)+datenum(year,1,1),TPobs(inxP,3),'r+','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),TPzt(25,:),'-');
%  set(gca,'ylim',[0 55]);
%  ylabel('ppb','fontsize',9)
%  xlabel('year','fontsize',9)
%  title('TP 30-31m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  %set(gca,'xlim',tlims)