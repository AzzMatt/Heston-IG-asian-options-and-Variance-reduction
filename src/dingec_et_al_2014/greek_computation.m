%Computing the greeks
function [Delta_greek,Gamma_greek]=greek_computation(dt,T,d,n,S0,K,V_0,alpha,beta,sigma,rho,r)

L=log(K/S0);
gammaval=sum(exp((1:d)*r*dt))/d;
C1=@(w) sum(arrayfun(@(qq) phi_xbar_xj(w,-1i,qq,dt,d,V_0,alpha,beta,sigma,rho,r),1:d));
C2=@(w) K*phi_xbar(w,dt,d,V_0,alpha,beta,sigma,rho,r);

%Integrand functions
integrand1 = @(w) arrayfun(@(ww) imag(exp(-1i*ww*L)./ww .* (((1i*ww+1)/d).*C1(ww) - (1i*ww)/S0 .* C2(ww))), w);
integrand2 = @(w) arrayfun(@(ww) imag(exp(-1i*ww*L) .* (((1i-ww)/(S0*d)).*C1(ww) + ((1i+ww)/S0^2).*C2(ww))), w);

%Integral using built in function
mu_W_prime=gammaval/2+(integral(integrand1,-Inf,Inf,'AbsTol', 1e-6, 'RelTol', 1e-6))/(2*pi);
mu_W_second=(integral(integrand2,-Inf,Inf,'AbsTol', 1e-6, 'RelTol', 1e-6))/(2*pi);


f_prime=zeros(n,1);
f_second=zeros(n,1);
for ii=1:n

        [Vpath,Ic]=Glasserman_pathV(d,dt,V_0,alpha,beta,sigma);
        Z=randn(d,1);
        [~,a,kv,s,sigmaquad,b]=conditionalexp(Z,Vpath,Ic,d,S0,K,alpha,beta,sigma,rho,r,dt);
        ai=a(2:end);
        si=s(2:end);
        first_term_f_prime=sum(si.*exp(ai.^2./2).*(normcdf(kv-ai)-normcdf(b-ai)))/d;
        second_term_f_prime=normpdf(kv)/sqrt(sigmaquad)*(sum(si.*exp(ai.*kv))/d-K);
        f_prime(ii)=(first_term_f_prime-second_term_f_prime)/S0;
        
        first_term_f_second=K^2*normpdf(b)/(sum(ai.*si.*exp(ai.*b))/d);
        second_term_f_second=normpdf(kv)/sigmaquad*((1/d)*sum(ai.*si.*exp(ai.*kv))-2*sqrt(sigmaquad)*(1/d)*sum(si.*exp(ai.*kv))+((1/d)*sum(si.*exp(ai.*kv))-K)*(sqrt(sigmaquad)-kv));
        f_second(ii)=(first_term_f_second+second_term_f_second)/S0^2;
        
      
end

Delta_greek=exp(-r*T)*(mean(f_prime)+mu_W_prime);
Gamma_greek=exp(-r*T)*(mean(f_second)+mu_W_second);
end