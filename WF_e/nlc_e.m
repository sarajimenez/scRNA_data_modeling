%% Non-linear constraint related with the jacobian 

function [F1, F2, F3] = nlc_e(param)

global s1 
global s2

x_1 = sym('x_1');
x_2 = sym('x_2');

% Functions
Xj=sigmoidal_s(0,[x_1,x_2],param);
% Jacobian of the function wrt x_1 and x_2
J=jacobian([Xj(1);Xj(2)],[x_1,x_2]);
% Evaluate the jacobian in the "desired states"
J_e=subs(J,{x_1,x_2},{s1,s2}); % Here we modify the states  
J_e=double(J_e);
% Eigen values 
ev=eig(J_e);

% The determinant of the jacobian evaluated at the desired state must be
% positive 

D=ev(1)*ev(2);

if D>0
    F1=0;
else
    F1=1;
end

% The trace of the jacobian evaluated at the desired state must be negative

T=ev(1)+ev(2);

if T<0
    F2=0;
else
    F2=1;
end

% The eigenvalues must be real then the term inside the square root must be
% positive

SRT=T^2-4*D;

if SRT>0
    F3=0;
else
    F3=1;
end

end