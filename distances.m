clear all; close all; clc

%Euclidean distances between the fixed points 

S1 = [0.4722, 1.0205];
S2 = [0.4830, 1.0126];
S3 = [0.4733, 1.0292];
S4 = [0.4820, 1.0130];

S = [S1; S2; S3; S4];

D = zeros(4,4);

for i=1:4
    
    for j=1:4
        
        D(i,j) = sqrt((S(i,1)-S(j,1))^2+(S(i,2)-S(j,2))^2);
        
    end
    
end

