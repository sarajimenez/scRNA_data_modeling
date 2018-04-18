clc; clear; close all; 

%% Computing the vector field for the "decoded system"

% Set of parameters 
param1 = [0.2450,-0.8541,-0.8033,-0.6813,0.6444,0.4517,0.9396,0.7381,0.7446,-0.3172,0.6456,0.9565]; % S1
param2 = [0.1460,0.4952,0.1653,-0.8734,0.7707,0.2729,0.9545,0.9975,0.9203,-0.4175,0.0632,0.9760]; % S2
param3 = [-0.3250,-0.9606,-0.0596,-0.3903,-0.5258,0.4042,1.0000,0.9826,0.5589,0.3473,0.8453,-0.8078]; % S3
param4 = [0.6777,0.4688,0.8196,-0.6336,0.8570,0.3318,0.9468,0.9882,0.8312,-0.0795,0.1405,0.9476]; % S4

solutions = [param1; param2; param3; param4];

% Fixed points

xe = [0.4609, 1.0408; ... 
      0.7564, 0.3608; ...
      0.4629, 1.0500; ...
      0.6338, -0.4808];

% Initial conditions equally distributed 
x1=linspace(0,3,20);
x2=linspace(0,3,20);

[x1_s,x2_s]=meshgrid(x1,x2);

% Pre-location
u_s=zeros(size(x1_s));
v_s=zeros(size(x2_s));

t=0;

%for j = 1:4
    for i = 1:numel(x1_s)
    
        Xprime_s = sigmoidal_s(t,[x1_s(i);x2_s(i)],solutions(1,:));
        u_s(i) = Xprime_s(1);
        v_s(i) = Xprime_s(2);    

    end
    
    figure(3)
    quiver(x2_s,x1_s,v_s,u_s),xlabel('x2'),ylabel('x1'),title('Vector field of the decoded system'),axis tight equal;
                
    hold on
    plot(xe(1,2),xe(1,1),'.k','MarkerSize',20) % Fixed point 
    text(xe(1,2)+0.1,xe(1,1)+(0.05)*1,num2str(1))     
    
%end




%% Plotting solutions on the vector field of the "decoded system"

hold on
    
for x10 = [0 3.0]
    for x20 = [0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
        [t, S] = ode45(@sigmoidal_s,[0,5],[x10,x20],[],param1); 
        plot(S(:,2),S(:,1),'r')
    end
end

for x20 = [0 3.0]
    for x10 = [0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
        [t, S] = ode45(@sigmoidal_s,[0,5],[x10,x20],[],param1); 
        plot(S(:,2),S(:,1),'r')
    end
end 
% 
% hold on
%     
% for x10 = [-1 3.0]
%     for x20 = [-1 -0.5 0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
%         [t, S] = ode45(@sigmoidal_s,[0,10],[x10,x20],[],param2); 
%         plot(S(:,2),S(:,1),'g')
%     end
% end
% 
% hold on
%     
% for x20 = [-1 3.0]
%     for x10 = [-1 -0.5 0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
%         [t, S] = ode45(@sigmoidal_s,[0,10],[x10,x20],[],param2); 
%         plot(S(:,2),S(:,1),'g')
%     end
% end
%
% hold on
% 
% for x10 = [0 3.0]
%     for x20 = [0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
%         [t, S] = ode45(@sigmoidal_s,[0,5],[x10,x20],[],param3); 
%         plot(S(:,2),S(:,1),'b')
%     end
% end
% 
% hold on
% 
% for x20 = [0 3.0]
%     for x10 = [0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
%         [t, S] = ode45(@sigmoidal_s,[0,5],[x10,x20],[],param3); 
%         plot(S(:,2),S(:,1),'b')
%     end
% end
% 
% hold on
% 
% for x10 = [-1 3.0]
%     for x20 = [-1 -0.5 0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
%         [t, S] = ode45(@sigmoidal_s,[0,15],[x10,x20],[],param1); 
%         plot(S(:,2),S(:,1),'m')
%     end
% end
% 
% for x20 = [-1 3.0]
%     for x10 = [-1 -0.5 0 0.3 0.5 1.0 1.5 2.0 2.5 3.0]
%         [t, S] = ode45(@sigmoidal_s,[0,15],[x10,x20],[],param1); 
%         plot(S(:,2),S(:,1),'m')
%     end
% end
