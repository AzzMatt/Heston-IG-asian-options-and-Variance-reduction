%Computing the price using the variance reduction technique
function [price,conf_int]=Varreduction(dt,T,d,n,S0,K,V_0,alpha,beta,sigma,rho,r)
g_0=(S0/d)*sum(exp((1:d).*(r*dt)))-K;
L=log(K/S0);
% Grid for numerical integration with trapezoidal rule
w_max = 100;  % upper limit (empirically set)
dw = 0.1;     % integration step
w = -w_max:dw:w_max;
w(w==0) = [];  % remove w=0 to avoid division by zero

% Vectorized computation of the integrand
integrand_vals = zeros(size(w));
for i = 1:length(w)
    integrand_vals(i) = imag(exp(-1i*w(i)*L)*g(w(i),d,S0,K,dt,V_0,alpha,beta,sigma,rho,r))/w(i);
end

% trapezoidal rule
mu_W = g_0/2 + trapz(w, integrand_vals)/(2*pi);

% Using the built in function for the integral, slower with higher values
% of d
%integrand=@(w) arrayfun(@(ww) imag(exp(-1i*ww*L)*g(ww,d,S0,K,dt,V_0,alpha,beta,sigma,rho,r))/ww,w);
%mu_W=g_0/2+integral(integrand,-Inf,Inf)/(2*pi);
Y=zeros(n,1);

for ii=1:n

        [Vpath,Ic]=Glasserman_pathV(d,dt,V_0,alpha,beta,sigma);
        Z=randn(d,1); 
        [condexp,~,~,~,~,~]=conditionalexp(Z,Vpath,Ic,d,S0,K,alpha,beta,sigma,rho,r,dt);
        
        
        Y(ii)=exp(-r*T)*(condexp+mu_W);

end

price=mean(Y);
conf_int=norminv(1-0.05/2)*std(Y)/sqrt(n);

end