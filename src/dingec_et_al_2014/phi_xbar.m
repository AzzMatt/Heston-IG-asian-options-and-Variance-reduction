%Computing the formula for the characteristic function that we need
function phi_xbar=phi_xbar(w1,dt,d,V_0,alpha,beta,sigma,rho,r)
u=zeros(d,1);
for j=1:d
    u(j)=w1*(d-j+1)/d;
end
phi_xbar=charfun(u,dt,d,V_0,alpha,beta,sigma,rho,r);
end