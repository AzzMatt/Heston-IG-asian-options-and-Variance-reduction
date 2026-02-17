%Computing the g function
function g=g(w,d,S0,K,dt,V_0,alpha,beta,sigma,rho,r)
% creating vectors for more efficiency
a=zeros(d,1);
b=zeros(d,1);
x=zeros(d,1);
for j=1:d
    a(j)=w*(d-j+1)/d-1i;
    b(j)=w*(d-j+1)/d;
end
for ii=1:d-1
    v=b;
    v(1:ii)=a(1:ii);
    x(ii)=charfun(v,dt,d,V_0,alpha,beta,sigma,rho,r);
end
x(d)=charfun(a,dt,d,V_0,alpha,beta,sigma,rho,r);
y=charfun(b,dt,d,V_0,alpha,beta,sigma,rho,r);
% final value from the formula
g=S0/d*sum(x)-K*y;
end