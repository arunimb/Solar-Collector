function [qtotal,etat]= general(ff)

T1=30+273;
T6=100+273;
l=0.001;
Lt=3;
G=ff;
Gin=G*l*Lt;
%-----glass parameter------%
t=0.9;
abg=0.04;
eg=abg;
idg=0.02;
odg=0.021;
Kg=4;
Hg=100;
%-----copper tube parameter----%
abc=0.70;
ec=abc;
idc=0.01;
odc=0.015;
Kc=385;


%------View Factors-----%
Aoc=pi*odc*l;
Aic=pi*idc*l;
Aog=pi*odg*l;
Aig=pi*idg*l;

F34=Aoc/Aig;
F33=1-F34;
%----Thermal Fluid Characteristic------%
m_dot=0.288;
rho=1060;
v=m_dot*4/(rho*pi*(idc)^2);
%----constants----%
si=5.67e-8;
Cg=2*pi*Kg*l/(log(odg/idg));
Cc=2*pi*Kc*l/(log(odc/idc));
e1=2;
T2=T1;
pr2=0.01;
pr1=0.1;
b=0;
T4=T6;
counter=0;
Length=0;
qt=0;
while Length<(76.6500)
    
   T=T6-273;
   cp=0.002414*T+5.9591e-6*T^2-2.9879e-8*T^3+4.4172e-11*T^4+1.498;
   cp=cp*1000;
   %k=-8.19477*10e-5*T-1.92257e-7*T^2+2.5034e-11*T^3-7.2974e-15*T^4+0.137;
   k=0.1;
   rho=-0.90797*T+0.00078116*T^2-2.367e-6*T^3+1083.25;
   gamma=exp(544.149/(T+114.3)-2.59578);
   vis=gamma*rho*1e-6;   
   Pr=cp*vis/k;
   Re=rho*v*idc/vis;  
   Nu=0.023*Re^0.8*Pr^0.33; 
   Hc=Nu*k/idc;
 
while abs(e1)>0.01
    
    T3=(-1*abg*Gin/(Cg))+T2+(eg*si*(T2^4-T1^4)*Aog/Cg)+(Hg*Aog*(T2-T1)/Cg);
    e2=2;
    c=0;
    T4=T6+pr2;
    while abs(e2)>0.001
        x=(Hc*Aic)/Cc;       
        T5=(x*T6+T4)/(1+x);
        e2=(F34*abc*si*T3^4*Aig)+abc*t*Gin-Cc*(T4-T5)-ec*si*T4^4*Aoc;%+F34*(1-t-abg)*abc*si*(T4^4)*Aoc;
        T4=T4+pr2;
   end   
    T4=T4-pr2;
    e1=+Cg*(T2-T3)+abg*ec*si*T4^4*Aoc-F34*eg*si*T3^4*Aig+F33*eg*si*T3^4*abg*Aig+abg*(1-abc)*t*Gin+abg*(1-abc)*F34*si*(T3^4)*Aig;
    T2=T2+pr1;
end
   q=Cc*(T4-T5);
   dT=q/(m_dot*cp);
   T6=T6+dT;
   counter=counter+1;
   Length=l*counter;
   qt=qt+q;
end
qtotal=qt;

etat=qtotal/(Lt*Length*G);
end
