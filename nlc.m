%% Non-linear constraint related with the jacobian 

function ev = nlc(param)

x_1 = sym('x_1');
x_2 = sym('x_2');

Xj=sigmoidal_s(0,[x_1,x_2],param);
J=jacobian([Xj(1);Xj(2)],[x_1,x_2]);
J_e=subs(J,{x_1,x_2},{0.5,2.5});
J_e=double(J_e);
ev=eig(J_e);

end