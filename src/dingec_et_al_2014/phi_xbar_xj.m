%Computing the formula for the characteristic function that we need
function phi_xbar_xj=phi_xbar_xj(w1,w2,nj,dt,d,V_0,alpha,beta,sigma,rho,r)
u=zeros(d,1);
for ii=1:nj
    u(ii)=w2+w1*(d-ii+1)/d;
end
for k=nj+1:d
    u(k)=w1*(d-k+1)/d;
end
phi_xbar_xj=charfun(u,dt,d,V_0,alpha,beta,sigma,rho,r);

end