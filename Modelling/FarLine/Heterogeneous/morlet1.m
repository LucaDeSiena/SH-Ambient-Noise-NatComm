function [f,nw]                 =   morlet1(fb,dt,E,rhoSource)
t0                              =   6/5/fb;
nw                              =   2.5*t0/dt;
lb                              =   -2;
ub                              =   2;
% lb                              =   -4;
% ub                              =   4;
[f,~]                           =   morlet(lb,ub,nw);
% f                               =...
%     iomega(sqrt(E*dt/rhoSource)*f/sum(abs(cumtrapz(f))),1,1,2);
