function dx = sigmoidal_s(t,x,param) % 12 parameters 

qa = param(1);
qb = param(2);
ka = param(3);
kb = param(4);
S = param(5);
n = param(6);
pa = param(7);
pb = param(8);

xa = x(1);
xb = x(2);

dxa = (pa.*xa.^n)./(S^n+xa.^n)+(qa*S^n)./(S^n+xb.^n)-ka.*xa;
dxb = (pb.*xb.^n)./(S^n+xb.^n)+(qb*S^n)./(S^n+xa.^n)-kb.*xb;

dx = [dxa dxb]';
    
end