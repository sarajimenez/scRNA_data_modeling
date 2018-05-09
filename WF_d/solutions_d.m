%{ 
Model for decoding data 
Model of solutions: Sigmoidal 
Following: Daniels, B. C., & Nemenman, I. (2015). Nature Communications, 6, 1?8. https://doi.org/10.1038/ncomms9133 
%}

function cf = solutions_d(param)

% Non-linear constraint related with the stability of the fixed points 
[F1, F2, F3, G] = nlc_d(param);

% Constraints related with the relation of the parameters "Shape of the sigmoid"

H1 = exp(param(9)*param(10));
H2 = exp(param(11)*param(12));
H3 = exp(-param(10)*param(11));
H4 = exp(-param(9)*param(12));
H5 = param(2)/param(3)-1;
H6 = param(5)/param(6)-1;
    
% Cost function
cf = 1*F1(1,1) + 1*F1(2,1) + 1*F2(1,1) + 1*F2(2,1) + 1*F3(1,1) + 1*F3(2,1) + 10*G(1,1) + 10*G(2,1) + 10*H1 + 10*H2 + 10*H3 + 10*H4 + 10*H5 + 10*H6;

end
