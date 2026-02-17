function ris=charf(u1,dt,d,V_0,alpha,beta,sigma,rho,r)

if ~isvector(u1) || length(u1) ~= d
    error('charf:dimensionError','u1 must be a vector of length d.');
end

u2=zeros(d,1);
for k=d-1:-1:1
    u2(k)=-1i*D(dt,u1(k+1),u2(k+1),alpha,sigma,rho);
end
sommare=zeros(d,1);
for ii=1:d
sommare(ii)=C(dt,r,u1(ii),u2(ii),alpha,beta,sigma,rho);
end

ris=exp(sum(sommare)+D(dt,u1(1),u2(1),alpha,sigma,rho)*V_0);% questo sarebbe il risultato della funzione caratteristica


function  b=b(w1,alpha,sigma,rho)
b=1i*w1*sigma*rho-alpha;
end

function h=h(w1,sigma,alpha,rho)
h=sqrt(sigma^2*1i*w1*(1i*w1-1)-b(w1,alpha,sigma,rho)^2);
end

function gamma=gamma(w1,w2,sigma,alpha,rho,dt)
gamma=tan(h(w1,sigma,alpha,rho)*dt/2+atan((b(w1,alpha,sigma,rho)+1i*w2*sigma^2)/2));
end


function C=C(dt,r,w1,w2,alpha,beta,sigma,rho)
C=1i*w1*r*dt+alpha*beta/sigma^2*(log(h(w1,sigma,alpha,rho)^2*(1+gamma(w1,w2,sigma,alpha,rho,dt)^2)/(sigma^2*(2*b(w1,alpha,sigma,rho)*1i*w2-sigma^2*w2^2+1i*w1*(1i*w1-1))))-b(w1,alpha,sigma,rho)*dt);
end

function D=D(dt,w1,w2,alpha,sigma,rho)
D=(gamma(w1,w2,sigma,alpha,rho,dt)*h(w1,sigma,alpha,rho)-b(w1,alpha,sigma,rho))/sigma^2;
end
end