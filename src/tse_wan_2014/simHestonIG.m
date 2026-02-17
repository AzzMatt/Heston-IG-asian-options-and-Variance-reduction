%Computing the price of the option
function estimation=simHestonIG(N,M,sigma,k,theta,V0,rho,r,T,X0,K)
dt=T/N;
% given by proposition 2.1
delta=4*k*theta/(sigma^2);
vu=delta/2-1;
% formulae given by propoposition 3.1 and equation 3
n_dt=4*k*exp(-k*dt)/(sigma^2*(1-exp(-k*dt)));
C_z=2*k*(sigma^2*sinh(k*dt/2))^(-1);
C_1=coth(k*dt/2);
C_2=csch(k*dt/2)^2;
E_X2=delta*sigma^2*(-2+k*dt*C_1)/(4*k^2);
sgX2=delta*sigma^4*(-8+2*k*dt*C_1+(k^2)*(dt^2)*C_2)/(8*(k^4));
E_Z=4*E_X2/delta;
sgZ=4*sgX2/delta;
V=zeros(N+1,M);
X=zeros(N+1,M);
V(1,:)=V0;
X(1,:)=X0;
% Preallocation of the gaussian values
Gaussian=randn(N,M);




for i=1:N
    coeff=2*exp(-k*dt)/n_dt;
    m_p=poissrnd(V(i,:).*(n_dt/2));
    V(i+1,:)=coeff.*gamrnd(m_p+delta/2,1); %computing the next value of the variance
    z=sqrt(V(i,:).*V(i+1,:)).*C_z;
    % computation for the moment matching Ic, formulae of proposition 3.1
    E_nu=z.*besseli(vu+1,z)./(2.*besseli(vu,z));
    E_nu2=z.^2.*besseli(vu+2,z)./(4.*besseli(vu,z))+E_nu;
    E_X1=(V(i,:)+V(i+1,:)).*(C_1/k-dt*C_2/2);
    sgX1=(V(i,:)+V(i+1,:)).*(sigma^2*C_1/k^3+sigma^2*dt*C_2/(2*k^2)-(sigma^2)*(dt^2)*C_1*C_2/(2*k));
    E_Ic=E_X1+E_X2+E_nu.*E_Z;
    Var_Ic=sgX1+sgX2+E_nu.*sgZ+(E_nu2-E_nu.^2).*(E_Z^2);
    % computation of Ic
    Ic=IG(E_Ic,(E_Ic.^3)./Var_Ic);
    % computing the log ratio
    log_ratio=r*dt-0.5.*Ic+(rho/sigma).*(V(i+1,:)-V(i,:)-k*theta*dt+k.*Ic)+sqrt(Ic.*(1-rho^2)).*Gaussian(i,:);
    X(i+1,:)=exp(log_ratio).*X(i,:);
    
end

% price for the European Call option
S_T=X(N+1,:);
payoffs=max(S_T-K,0);
estimation=exp(-r*T)*mean(payoffs);

%Price for Asian Call option

% Na_values=[2, 4, 8, 16];
% estimated_price=zeros(size(Na_values));
% for j=1:length(Na_values)
%     Na=Na_values(j);
% 
%     sums=X(N/Na+1,:);
%     for ii=2:Na
%         sums = sums + X(ii * (N/Na) + 1, :);
%     end
%     sums=sums./Na;
%     payoffs=max(sums-K,0);
%     estimated_price(j)=exp(-r*T).*mean(payoffs);

end