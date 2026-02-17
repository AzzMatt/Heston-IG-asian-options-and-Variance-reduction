%Computing the path of the variance as in Glasserman book
function [V,Ic]=Glasserman_pathV(d,dt,V_0,alpha,beta,sigma)
%from pag 124 Glasserman
dparam=4*alpha*beta/sigma^2;
k=20;%number used for the integral
dt_small=dt/k;
V=zeros(k*d+1,1);
Ic=zeros(d,1);
V(1)=V_0;
% computation of the path of V(t)
if dparam>1
    for ii=1:k*d
        c=sigma^2*(1-exp(-alpha*dt_small))/(4*alpha);
        lambda=V(ii)*exp(-alpha*dt_small)/c;
        Z=randn(1);
        X=chi2rnd(dparam-1);
        V(ii+1)=c*((Z+sqrt(lambda))^2+X);
    end
end
if dparam<=1
    for ii=1:k*d
        c=sigma^2*(1-exp(-alpha*dt_small))/(4*alpha);
        lambda=V(ii)*exp(-alpha*dt_small)/c;
        N=poissrnd(lambda/2);
        X=chi2rnd(dparam+2*N);
        V(ii+1)=c*X;

    end
end
% computing the integral with routine intration technique
for j=1:d
    idx_start=(j-1)*k+1;
    idx_end=j*k+1;
    V_seg=V(idx_start:idx_end);
    for m = 1:k
        Ic(j) = Ic(j) + (V_seg(m) + V_seg(m+1))/2 * dt_small;
    end
end
% taking back V(t) to its d+1 values
V=V(1:k:end);


end