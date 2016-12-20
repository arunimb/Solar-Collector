clear all
clc

Tw=100+273;
T1=273;
Ta=30+273;

l=0.001;
m_dot=0.288;
%------Parabola Characteristic-------%
Lt=3;
Girr=684;
Gf=Girr*Lt*l;
%-----end----%
%---Tube characteristics-----%
id=0.01;
od=0.015;
K=385;

UAcond=2*pi*385*l/(log(od/id));
ab=0.7;
emi=0.7;
Gabs=ab*Gf;
%---End-------%
%----Thermal Fluid Characteristic------%
rho=1060;
v=m_dot*4/(rho*pi*(id)^2);

%-----end-----%
%------Ambient Conditions-----%
ha=60;
UAconv_o=ha*pi*od*l;
d=2;
counter=0;
qtotal=0;
while Tw<(200+273)
    T=Tw-273;
    
    cp=0.002414*T+5.9591e-6*T^2-2.9879e-8*T^3+4.4172e-11*T^4+1.498;
   cp=cp*1000;
%    k=-8.19477*10e-5*T-1.92257e-7*T^2+2.5034e-11*T^3-7.2974e-15*T^4+0.137;
    k=0.1;
   rho=-0.90797*T+0.00078116*T^2-2.367e-6*T^3+1083.25;
   gamma=exp(544.149/(T+114.3)-2.59578);
   vis=gamma*rho*1e-6;
   
   Pr=cp*vis/k;
   Re=rho*v*id/vis;
   
   Nu=0.023*Re^0.8*Pr^0.33;
   
   Htf=Nu*k/id;
   UAconv_i=pi*id*l*Htf;
   
while abs(d)>=0.001
    
    cond_ratio=UAcond/UAconv_i;
    T2=(cond_ratio*T1+Tw)/(1+cond_ratio);
    d1=UAcond*(T1-T2);
    d2=-1*Gabs+emi*(5.67e-8)*(T1^4-Ta^4)*pi*od*l+UAconv_o*(T1-Ta);
    d=d1+d2;
    T1=T1+0.001;
end
   q=UAcond*(T1-T2);
   T=Tw-273;
   cp=0.002414*T+5.9591e-6*T^2-2.9879e-8*T^3+4.4172e-11*T^4+1.498;
   cp=cp*1000;
   dT=q/(m_dot*cp);
   Tw=Tw+dT;
   counter=counter+1;
   
  qtotal=qtotal+q;
end
Length=counter*l
 qtotal
 Length*Girr
eta=qtotal/(Length*Girr*Lt)