%% Non-linear constraint related with the jacobian 

function [F1, F2, F3] = nlc_a(param)

x_1 = sym('x_1');
x_2 = sym('x_2');

% Functions
Xj = sigmoidal_s(0,[x_1,x_2],param);

% Jacobian of the function wrt x_1 and x_2
J = jacobian([Xj(1);Xj(2)],[x_1,x_2]);

% Pre-location
ev = zeros(2,2);
D = zeros(2,1);
T = zeros(2,1);
SRT = zeros(2,1);
F1 = zeros(2,1);
F2 = zeros(2,1);
F3 = zeros(2,1);

s = [0 1.5; 1.5 0]; % Desired states (each row is a state)

for j = 1:2
       
    % Evaluate the jacobian in the "desired states"
    J_e = subs(J,{x_1,x_2},{s(j,1),s(j,2)}); % Here we modify the states  
    J_e = double(J_e);
    
    % Eigen values 
    ev(:,j) = eig(J_e);
    
    % Determinant 
    D(j) = ev(1,j).*ev(2,j);
    
    F1(j) = exp(-D(j));
    
    % Trace
    T(j) = ev(1,j)+ev(2,j);
    
    F2(j) = exp(T(j));
    
    % Term inside the square root
    SRT(j) = T(j)^2-4.*D(j);
    
    F3(j) = exp(-SRT(j));
        
end

assignin('base','ev',ev) % Save the eigenvalues in the workspace

end