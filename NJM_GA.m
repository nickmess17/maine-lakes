function [SS_temp] = NJM_GA(par)
% par(1)=0.0241;%open water diffusion scaling coefficient
% par(2)=0.93;%wind sheltering coefficient
% par(3)=1.66;%scaling of inflow volume
% par(4)=2.5;%non PAR light attenuation coefficient
% par(5)=1.05;%PAR light attenuation coefficient
f=@NJM_MyLake;
SS_temp=NJM_MyLake(par);
end

