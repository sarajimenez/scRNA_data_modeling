%{ 
Model for decoding data 
Model of solutions: Sigmoidal 
Following: Wang et al. 2010. Biophysical Journal Volume 99 July 2010 29?39. https://doi.org/10.1016/j.bpj.2010.03.058
%}

function [cf]  = solutions_SH(param)

% Non-linear constraint related with the stability of the fixed points 
[F1, F2, F3, G] = nlc_SH(param);
    
% Cost function
cf =  1*F1(1,1) + 1*F1(2,1) + 1*F2(1,1) + 1*F2(2,1) + 1*F3(1,1) + 1*F3(2,1) + 10*G(1,1) + 10*G(2,1) ...
    + 1*F1(3,1) + 1*F1(4,1) + 1*F2(3,1) + 1*F2(4,1) + 1*F3(3,1) + 1*F3(4,1) + 10*G(3,1) + 10*G(4,1);


end
