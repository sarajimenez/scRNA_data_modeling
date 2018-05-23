function dx = sigmoidal_s(t,x,param) % 12 parameters 

qa = param(3);
qb = param(4);
ka = 1;
kb = 1;
S = param(5);
n = param(6);
pa = param(1);
pb = param(2);

xa = x(1);
xb = x(2);

dxa = (pa.*xa.^n)./(S^n+xa.^n)+(qa*S^n)./(S^n+xb.^n)-ka.*xa;
dxb = (pb.*xb.^n)./(S^n+xb.^n)+(qb*S^n)./(S^n+xa.^n)-kb.*xb;

dx = [dxa dxb]';
    
end