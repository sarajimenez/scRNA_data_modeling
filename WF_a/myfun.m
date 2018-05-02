%{
Model for generating data
Following: Wang et al. 2010. Biophysical Journal Volume 99 July 2010 29?39. https://doi.org/10.1016/j.bpj.2010.03.058
%}

function dx = myfun(t,x)

% Parameters of the model for generating data

qa = 1;
qb = 1;
ka = 1;
kb = 1;
S = 0.5;
n = 4;

% pa = 1.2; % High autoactivation of xa
% pb = 1.2; % High expression of xb

pa = 0.6; % Low autoactivation of xa
pb = 0.6; % Low autoactivation of xb

xa = x(1);
xb = x(2);

dxa = (pa.*xa.^n)./(S^n+xa.^n)+(qa*S^n)./(S^n+xb.^n)-ka.*xa;
dxb = (pb.*xb.^n)./(S^n+xb.^n)+(qb*S^n)./(S^n+xa.^n)-kb.*xb;

dx = [dxa dxb]';

end