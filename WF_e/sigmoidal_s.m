function dx = sigmoidal_s(t,x,param) % 12 parameters 

    x1 = x(1);
    x2 = x(2);
    
    % Sigmoid x1: a1=param(1), b1=param(2), c1=param(3).
    s1 = param(1) + 1/(1+exp(-param(3)*(x1-param(2)/param(3))));
    % Sigmoid x2: a2=param(4), b2=param(5), c2=param(6).
    s2 = param(4) + 1/(1+exp(-param(6)*(x2-param(5)/param(6))));
    
    % Input
    xI = 0;
    
    % W11=param(7), W12=param(8)
    dx1 = -x1./1.0+param(7).*s1+param(8).*s2+xI;
    
    % W21=param(9), W22=param(10)    
    dx2 = -x2./1.0+param(9).*s1+param(10).*s2+xI;
    
    dx = [dx1 dx2]';
    
end