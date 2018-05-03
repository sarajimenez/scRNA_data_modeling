function x2_nl = zerolines_2(x1,param) % 12 parameters 

    x2(1) = param(13);

    % Sigmoid x1: a1=param(1), b1=param(2), c1=param(3).
    s1 = param(1)+(1-param(1))./(1+exp((-4*param(2)/param(3))*(x1-param(3))/(1-param(1))));
    % Sigmoid x2: a2=param(4), b2=param(5), c2=param(6).
    s2 = param(4)+(1-param(4))./(1+exp((-4*param(5)/param(6))*(x2-param(6))/(1-param(4))));
    
    % W21=param(11), W22=param(12)    
    x2 = param(8).*(param(11).*s1+param(12).*s2);
    
    x2_nl = x2;
    
end
