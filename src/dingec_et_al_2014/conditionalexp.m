%Computing conditional expectation and quantities useful later
function [ris,a,kv,s,sigmaquad,b]=conditionalexp(Z,Vpath,Ic,d,S0,K,alpha,beta,sigma,rho,r,dt)
% creating the vector for efficiency
Gamma=zeros(d,1);
Lambda=zeros(d,1);
Theta=zeros(d,1);
coeff=0;
coeff2=0;
coeff3=0;

for j=1:d
    Gamma(j)=(r-alpha*beta*rho/sigma)*dt+(rho/sigma)*(Vpath(j+1)-Vpath(j))+(alpha*rho/sigma-0.5)*Ic(j);
    Lambda(j)=sqrt((1-rho^2)*Ic(j));
    coeff=coeff+(d-j+1)^2*Lambda(j)^2;
    coeff2=coeff2+Gamma(j)*(d-j+1);
    coeff3=coeff3+Lambda(j)^2*(d-j+1)^2;
end
for jj=1:d
    Theta(jj)=(d-jj+1)*Lambda(jj)/sqrt(coeff);
end
mu=log(S0)+coeff2/d;
sigmaquad=coeff3/d^2;
kv=(log(K)-mu)/sqrt(sigmaquad);
xi=Z-Theta*(Theta'*Z);
s=zeros(d+1,1);
a=zeros(d+1,1);
s(1)=S0;
for ii=2:d+1
    a(ii)=a(ii-1)+Lambda(ii-1)*Theta(ii-1);
    s(ii)=s(ii-1)*exp(Gamma(ii-1)+Lambda(ii-1)*xi(ii-1)); 
end
%setting Netwon's method
pfun=@(x) sum(exp(a(2:d+1).*x).*s(2:d+1))/d-K;
Dpfun=@(x) sum(s(2:d+1).*a(2:d+1).*exp(a(2:d+1).*x))/d;
x0=kv-pfun(kv)/Dpfun(kv); %point set as in the paper
b=Newton(pfun,Dpfun,x0,1e-8,100);
vec=zeros(d,1);
for jj=1:d
    vec(jj)=s(jj+1)*exp(a(jj+1)^2/2)*(normcdf(kv-a(jj+1))-normcdf(b-a(jj+1)));
end
ris=sum(vec)/d-K*(normcdf(kv)-normcdf(b));
end