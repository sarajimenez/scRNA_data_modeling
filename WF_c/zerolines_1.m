function x1_nl = zerolines_1(x2,param) % 12 parameters 

    x1(1) = param(14);

    % Sigmoid x1: a1=param(1), b1=param(2), c1=param(3).
    s1 = param(1)+(1-param(1))./(1+exp((-4*param(2)/param(3))*(x1-param(3))/(1-param(1))));
    % Sigmoid x2: a2=param(4), b2=param(5), c2=param(6).
    s2 = param(4)+(1-param(4))./(1+exp((-4*param(5)/param(6))*(x2-param(6))/(1-param(4))));
    
    % W11=param(9), W12=param(10)
    x1 = param(7).*(param(9).*s1+param(10).*s2);
    
    x1_nl = x1;
    
end