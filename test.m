
close all; clc; clear all
%%

param=[1708/2557,-1,451/593,-1,335/849,769/789,1,1,385/804,-1,1,543/994];

syms x_1 x_2

Xj=sigmoidal(0,[x_1,x_2],param);
J=jacobian([Xj(1);Xj(2)],[x_1,x_2]);
J_e=subs(J,{x_1,x_2},{0.5,2.5});
J_e=double(J_e);
ev=eig(J_e);


function dx = sigmoidal(t,x,param) % 12 parameters 

    x1=x(1);
    x2=x(2);
    
    % Sigmoid x1: a1=param(1), b1=param(2), c1=param(3).
    s1=param(1)+(1-param(1))./(1+exp((-4*param(2)/param(3))*(x1-param(3))/(1-param(1))));
    % Sigmoid x2: a2=param(4), b2=param(5), c2=param(6).
    s2=param(4)+(1-param(4))./(1+exp((-4*param(5)/param(6))*(x2-param(6))/(1-param(4))));
    % Input
    xI=0;
    
    % W11=param(9), W12=param(10)
    dx1=-x1./param(7)+param(9).*s1+param(10).*s2+xI;
    
    % W21=param(11), W22=param(12)    
    dx2=-x2./param(8)+param(11).*s1+param(12).*s2+xI;
    
    dx=[dx1 dx2]';
    
end
    