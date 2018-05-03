%% Non-linear constraint related with the jacobian 

function [F1, F2, F3] = nlc_b(param)

x_1 = sym('x_1');
x_2 = sym('x_2');

% Functions
Xj = sigmoidal_s(0,[x_1,x_2],param);

% Jacobian of the function wrt x_1 and x_2
J = jacobian([Xj(1);Xj(2)],[x_1,x_2]);

% Pre-location
ev = zeros(2,3);
D = zeros(3,1);
T = zeros(3,1);
SRT = zeros(3,1);
F1 = zeros(3,1);
F2 = zeros(3,1);
F3 = zeros(3,1);

s = [0 1; 0.5 0.5; 1 0]; % Desired states (each row is a state)

for j = 1:3
       
    % Evaluate the jacobian in the "desired states"
    J_e = subs(J,{x_1,x_2},{s(j,1),s(j,2)}); % Here we modify the states  
    J_e = double(J_e);
    
    % Eigen values 
    ev(:,j) = eig(J_e);
    
    % Determinant 
    D(j,1) = ev(1,j).*ev(2,j);
    
    F1(j,1) = exp(-D(j,1));
    
    % Trace
    T(j,1) = ev(1,j)+ev(2,j);
    
    F2(j,1) = exp(T(j,1));
    
    % Term inside the square root
    SRT(j,1) = T(j,1)^2-4.*D(j,1);
    
    F3(j,1) = exp(-SRT(j,1));
        
end

assignin('base','ev',ev) % Save the eigenvalues in the workspace

% The determinant of the jacobian evaluated at the desired state must be
% positive 

% if D > 0
%     F1 = 0;
% else
%     F1 = 1;
% end

%F1 = sum(F1);

% The trace of the jacobian evaluated at the desired state must be negative

% if T < 0
%     F2 = 0;
% else
%     F2 = 1;
% end

%F2 = sum(F2);

% The eigenvalues must be real then the term inside the square root must be
% positive

% if SRT > 0
%     F3 = 0;
% else
%     F3 = 1;
% end

%F3 = sum(F3);

end