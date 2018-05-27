%{ 
Model for decoding data 
Model of solutions: Sigmoidal 
Following: Daniels, B. C., & Nemenman, I. (2015). Nature Communications, 6, 1?8. https://doi.org/10.1038/ncomms9133 
%}

function cf = solutions_e(param)

% Non-linear constraint related with the stability of the fixed points 
[F1, F2, F3, G] = nlc_e(param);

assignin('base','F1',F1) % Save the jacobian in the workspace
assignin('base','F2',F2) % Save the jacobian in the workspace
assignin('base','F3',F3) % Save the jacobian in the workspace
assignin('base','G',G) % Save the jacobian in the workspace


% Cost function
cf = 1*F1(1,1) + 1*F1(2,1) + 1*F2(1,1) + 1*F2(2,1) + 1*F3(1,1) + 1*F3(2,1) + 10*G(1,1) + 10*G(1,2) + 10*G(2,1) + 10*G(2,2);

end
