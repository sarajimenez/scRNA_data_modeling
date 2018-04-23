function dx = sigmoidal_s(t,x,param) % 12 parameters 

    x1=x(1);
    x2=x(2);
    
    % Sigmoid x1: a1=param(1), b1=param(2), c1=param(3).
    s1=param(1)+(1-param(1))./(1+exp((-4*param(2))*(x1-param(3))/(1-param(1))));
    % Sigmoid x2: a2=param(4), b2=param(5), c2=param(6).
    s2=param(4)+(1-param(4))./(1+exp((-4*param(5))*(x2-param(6))/(1-param(4))));
    
%     % Sigmoid x1: a1=param(1), b1=param(2), c1=param(3).
%     s1=param(1)+1./(1+exp(-param(2)*(x1-param(3))));
%     % Sigmoid x2: a2=param(4), b2=param(5), c2=param(6).
%     s2=param(4)+1./(1+exp(-param(5)*(x2-param(6))));
    
    % Input
    xI=0;
    
    % W11=param(9), W12=param(10)
    dx1=-x1./param(7)+param(9).*s1+param(10).*s2+xI;
    
    % W21=param(11), W22=param(12)    
    dx2=-x2./param(8)+param(11).*s1+param(12).*s2+xI;
    
    dx=[dx1 dx2]';
    
end