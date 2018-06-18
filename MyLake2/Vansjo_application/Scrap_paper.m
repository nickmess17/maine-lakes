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

figure(3)
clf
fign=6;
inxpos=find(isnan(TempMod(alku))==0);
posalku=alku(inxpos);
posloppu=loppu(inxpos);

for i = 1:min(fign,length(posalku))
   subplot(3,3,i);
   inxt=find(tt_mod==TempObs(posalku(i),1));
   plot(Tzt(:,inxt),zz,'-b',TempObs(posalku(i):posloppu(i),3),TempObs(posalku(i):posloppu(i),2),'.r');
   axis([0 25 zlim])
   axis ij
   title(['Vansjø-Stfj. ' datestr(tt_mod(inxt)+datenum(year,1,1)),],'fontsize',8); 
   set(gca,'FontSize',8) 
end;

   subplot(334)
   ylabel('Depth (m)','fontsize',8)
   xlabel('Temperature (^{o}C)','fontsize',8)
   
   subplot(331)
   ylabel('Depth (m)','fontsize',8)
    
    subplot(335)
   xlabel('Temperature (^{o}C)','fontsize',8)
   
     subplot(336)
   xlabel('Temperature (^{o}C)','fontsize',8)

figure(4)
clf
for i = fign+1:min(12,length(posalku))
   subplot(3,3,i-fign);
   set(gca,'FontSize',8) 
   inxt=find(tt_mod==TempObs(posalku(i),1));
   plot(Tzt(:,inxt),zz,'-b',TempObs(posalku(i):posloppu(i),3),TempObs(posalku(i):posloppu(i),2),'.r');
   axis([0 25 zlim])
   axis ij
   title(['Vansjø-Stfj. ' datestr(tt_mod(inxt)+datenum(year,1,1)),],'fontsize',8); 
end;

   subplot(334)
   ylabel('Depth (m)','fontsize',8)
   xlabel('Temperature (^{o}C)','fontsize',8)
   
    subplot(331)
   ylabel('Depth (m)','fontsize',8)
   
    subplot(335)
   xlabel('Temperature (^{o}C)','fontsize',8)
   
     subplot(336)
   xlabel('Temperature (^{o}C)','fontsize',8)

figure(8)
clf
zinx=find(zz<4); %depth layer considered
kk=0;
dum=datevec(tt_mod+datenum(year,1,1));
yrs=dum(:,1);
months=dum(:,2);

monthly_qual=NaN*ones(12*(max(yrs)-min(yrs)+1),12);
F_OM=1e+6*0.012;    %mass fraction [mg kg-1] of P of dry organic matter (assuming 50% of C, and Redfield ratio)
for i=min(yrs):max(yrs)
   for j=1:12
    kk=kk+1;   
    tinx=find((yrs==i)&(months==j)); 

    monthly_qual(kk,1)=i; %year
    monthly_qual(kk,2)=j; %month
    monthly_qual(kk,3)=15; %day
    monthly_qual(kk,4)=datenum(i,j,15); %dtenumber
    monthly_qual(kk,5)=mean(mean(Pzt(zinx,tinx)+PPzt(zinx,tinx)+DOPzt(zinx,tinx)+Chlzt(zinx,tinx))); %TotP
    monthly_qual(kk,6)=mean(mean(Chlzt(zinx,tinx))); %Chla
    monthly_qual(kk,7)=mean(mean(Szt(zinx,tinx)+Chlzt(zinx,tinx)./F_OM)); %SS
    monthly_qual(kk,8)=mean(mean(Chlzt(zinx,tinx))); %Chla
    monthly_qual(kk,9)=mean(mean(Pzt(zinx,tinx))); %P
    monthly_qual(kk,10)=mean(mean(PPzt(zinx,tinx))); %PP
    monthly_qual(kk,11)=mean(mean(DOPzt(zinx,tinx))); %DOP
       zdum=-9.6 + 3.6*log10(Pzt(zinx,tinx)+PPzt(zinx,tinx)+DOPzt(zinx,tinx)+Chlzt(zinx,tinx)) + 0.23*Tzt(zinx,tinx);
       if(isempty(zdum)==0)
        monthly_qual(kk,12)=max(mean(1./(1+exp(-zdum)))); %max(P[>10% cyanobacteria])
       end
   end
 end
 
 %synchronize observations to model simulations
         TP_Chla_Mod=NaN*ones(length(Obs_TP_Chla),2); 
        for i=1:length(Obs_TP_Chla)
            inx=find(datenum(Obs_TP_Chla(i,1:3))==monthly_qual(:,4));
            if(isempty(inx)==0)
                TP_Chla_Mod(i,1)=monthly_qual(inx,5); %1) TP
                TP_Chla_Mod(i,2)=monthly_qual(inx,6); %2) Chla
            end
        end
 
 yearly_qual=NaN*ones((max(yrs)-min(yrs)+1),6);
 kk=0;
for i=min(yrs):max(yrs)
    kk=kk+1;   
    tinx=find((monthly_qual(:,1)==i)&(monthly_qual(:,2)>5)&(monthly_qual(:,2)<10)); 
    yearly_qual(kk,1)=i; %year
    yearly_qual(kk,2)=mean(monthly_qual(tinx,5)); %TotP
    yearly_qual(kk,3)=mean(monthly_qual(tinx,6)); %Chla
    yearly_qual(kk,4)=mean(monthly_qual(tinx,7)); %SS
    yearly_qual(kk,5)=mean(monthly_qual(tinx,10)); %PIP
    yearly_qual(kk,6)=max(monthly_qual(tinx,12)); %P[cyano>10%]
 end
 
subplot(311)
 plot(monthly_qual(:,4),monthly_qual(:,5),'.-')
 hold on
 plot(datenum(Obs_TP_Chla(:,1:3)),Obs_TP_Chla(:,4),'r.')
 set(gca,'ylim',[0 70]);
 ylabel('TotP  (\mug/l)')
 title('Monthly mean (0-4m)')
 datetick('x','mmmyy')
 grid on
 
 subplot(312)
 plot(monthly_qual(:,4),monthly_qual(:,6),'.-')
 hold on
 plot(datenum(Obs_TP_Chla(:,1:3)),Obs_TP_Chla(:,5),'r.')
 set(gca,'ylim',[0 20]);
 ylabel('Chl a  (\mug/l)')
 title('Monthly mean (0-4m)')
 datetick('x','mmmyy')
 grid on
 
 subplot(313)
 plot(monthly_qual(:,4),monthly_qual(:,9),'.-')
 hold on
 plot(datenum(Obs_PO4P(:,1:3)),Obs_PO4P(:,4),'r.') %mg/m3
 set(gca,'ylim',[0 40]);
 ylabel('PO4 (\mug/l)')
 title('Monthly mean (0-4m)')
 datetick('x','mmmyy')
 grid on
 

 figure(15)
 clf
 subplot(211)
 plot(tt_mod, cumsum(MixStat(13,:)),'r')
 hold on
 plot(tt_mod, cumsum(MixStat(14,:)),'b')
 plot(tt_mod, cumsum(MixStat(15,:)),'c')
 plot(tt_mod, (MixStat(18,:)),'m')
 plot(tt_mod, cumsum(MixStat(13,:)-MixStat(14,:)-MixStat(15,:)+MixStat(16,:)+MixStat(17,:))-MixStat(18,:),'k')
 
 datetick('x','mmm');
 H=legend('Inflow of P', 'Outflow of P', 'Sedimentation of P', 'Change in lake P', 'P-Balance');
 set(H,'fontsize',8);
 grid on;
  ylabel('kg')
  
 dum=datevec(tt_mod+datenum(year,1,1));
 yrs=dum(:,1);
 kk=0;
 Intern=NaN*ones(length(max(yrs)-min(yrs))+1,5);
 for i=min(yrs):max(yrs)
 inx=find(yrs==i);
 kk=kk+1;
 Intern(kk,1)=sum(MixStat(16,inx));
 Intern(kk,2)=sum(MixStat(17,inx));
 Intern(kk,3)=sum(MixStat(19,inx)); %Net flux out of sediment
 Intern(kk,4)=sum(MixStat(13,inx)); %Inflow
 Intern(kk,5)=sum(MixStat(14,inx)); %Outflow
 Intern(kk,6)=sum(MixStat(20,inx)); %Algae available inflow
 end    
 
 subplot(212)
 plot([min(yrs)+1:max(yrs)], -1e-3*Intern(2:end,3),'k','linewidth',2)
 hold on
 plot([min(yrs)+1:max(yrs)], 1e-3*Intern(2:end,4),'r','linewidth',2)
 plot([min(yrs)+1:max(yrs)], 1e-3*Intern(2:end,6),'r--','linewidth',2)
 plot([min(yrs)+1:max(yrs)], 1e-3*Intern(2:end,5),'b','linewidth',2)

 H=legend('Net TP flux to sed.','TP inflow','P04 inflow','TP outflow');
 title('P budget','fontweight','bold')
 set(H,'fontsize',8);
 ylabel('tons/year')
 grid on;
 
 
 %=Figures for modelpaper
 tlims=[datenum([1984,12,15]) datenum([2001,1,15])];

figure(22)
clf
subplot(311)
inx=find(round(TempObs(:,2))<2);
H=plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),Tzt(1,:),'-');
 set(gca,'ylim',[0 25]);
 ylabel('^oC','fontsize',9)
 title('Temperature  0-1m','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xticklabel',[]);
set(gca,'xlim',tlims)
 
subplot(312)
inx=find((round(TempObs(:,2))==10)|(round(TempObs(:,2))==11));
plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),Tzt(11,:),'-');
 set(gca,'ylim',[0 25]);
 ylabel('^oC','fontsize',9)
 title('Temperature  10-11m','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xticklabel',[]);
 set(gca,'xlim',tlims)
 
 subplot(313)
inx=find((round(TempObs(:,2))==30)|(round(TempObs(:,2))==31));
plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),Tzt(31,:),'-');
 set(gca,'ylim',[0 25]);
 ylabel('^oC','fontsize',9)
 xlabel('year','fontsize',9)
 title('Temperature 30-31m','fontweight','bold')
 datetick('x','yy')
 grid on
 set(gca,'fontsize',9)
 set(gca,'xlim',tlims)
 %===============
 
 disp(['Sum of P sinks: '  num2str(round(sum(MixStat(14,:)+MixStat(15,:)))) ' kg']); 
 disp(['Sum of P sources: '  num2str(round( sum(sum(Intern)) + sum(MixStat(13,:)) )) ' kg']);
 
 disp(['Average summer season TotP, Chla, and SS: ' num2str(mean(yearly_qual(:,2:4)))])