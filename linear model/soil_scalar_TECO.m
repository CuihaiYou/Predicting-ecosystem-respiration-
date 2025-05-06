% this function is for Tem and Mois response function
% using the TECO response function from ALei
function y=tau_tem_moist_TECO(Tem_soil,RH_soil)
mscut=0.2;
temp=Tem_soil;
moist_teco=RH_soil./100;
tmp=2.^((temp-10)./10);
if moist_teco>mscut
    moisture1=1;
else
moisture1=1.0-5*(mscut-moist_teco); % 5*moist when mscut=0.2
end
y=tmp.*moisture1;
  
  