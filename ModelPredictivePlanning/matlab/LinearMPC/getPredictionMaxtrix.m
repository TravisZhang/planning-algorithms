function [Tp,Tv,Ta,Tj,Bp,Bv,Ba,Bj] = getPredictionMaxtrix(K,dt,p0,v0,a0,j0)
Ta = zeros(K);
Tv = zeros(K);
Tp = zeros(K);
Tj = eye(K);
Ba = zeros(K,1);
Bv = zeros(K,1);
Bp = zeros(K,1);
Bj = zeros(K,1);

for i = 1:K
    for j = 1:i
        Ta(i,j) = dt;
        Tv(i,j) = (1/2+(i-j))*dt^2;
        Tp(i,j) = (1/6+1/2*(i-j+1)*(i-j))*dt^3;
    end
    Ba(i) = a0;
    Bv(i) = v0+i*a0*dt;
    Bp(i) = p0+i*v0*dt+1/2*i^2*a0*dt^2;
    Bj(i) = j0;
end

end

