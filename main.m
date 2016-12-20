
clear all
clc

qreq=8.5601e+04;

n=120; % day of a year
phi=0.731293;
Is=1367;
I=Is*(1+0.034*cos(2*pi*n/265.25));
I=I/2;
delta=23.45*pi/180*sin(2*pi*(284+n)/36.25);
sunset_angle=acos(-1*tan(phi)*tan(delta));
N=2*sunset_angle*180/(15*pi);
alpha=N/(180);
x=linspace(0,180,10);
q=ones(1,10);
h=alpha.*x;
eta=zeros(1,10);
for i=1:10
    
    Ir=I*sind(x(i));
    if Ir>=12
    [qtol,efficiency]= general(Ir);
    q(i)=qtol;
    eta(i)=efficiency;
    end
    if I<=12
        q(i)=0;
        eta(i)=0;
    end
end
figure(1);
plot(h,q)
title('Power Collected');
xlabel('Time of day');
ylabel('Power Collected');


figure(2)
qneed=qreq*ones(1,10);
qneed=qneed-q;
plot(h,qneed)
title('Auxiliary Power Required');
xlabel('Time of day');
ylabel('Auxiliary Power Required');


figure(3)
plot(h,eta)
title('Efficiency of Solar Collector');
xlabel('Time of day');
ylabel('Efficiency of Solar Collector');
axis([0 h(10) 0.55 0.65])