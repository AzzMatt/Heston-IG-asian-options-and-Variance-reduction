%Simultating the variance path using IG, from the paper of Tse and Wan
function [V,Ic]=pathsim(d,dt,alpha,beta,sigma,V_0)

delta=4*alpha*beta/(sigma^2);
vu=delta/2-1;

n_dt=4*alpha*exp(-alpha*dt)/(sigma^2*(1-exp(-alpha*dt)));
C_z=2*alpha*(sigma^2*sinh(alpha*dt/2))^(-1);
C_1=coth(alpha*dt/2);
C_2=csch(alpha*dt/2)^2;
E_X2=delta*sigma^2*(-2+alpha*dt*C_1)/(4*alpha^2);
sgX2=delta*sigma^4*(-8+2*alpha*dt*C_1+(alpha^2)*(dt^2)*C_2)/(8*(alpha^4));
E_Z=4*E_X2/delta;
sgZ=4*sgX2/delta;
V=zeros(d+1,1);

V(1)=V_0;
Ic=zeros(d,1);




for i=1:d
    coeff=2*exp(-alpha*dt)/n_dt;
    m_p=poissrnd(V(i)*(n_dt/2));
    V(i+1)=coeff*gamrnd(m_p+delta/2,1);
    z=sqrt(V(i,:)*V(i+1,:))*C_z;
    E_nu=z.*besseli(vu+1,z)/(2*besseli(vu,z));
    E_nu2=z.^2.*besseli(vu+2,z)/(4*besseli(vu,z))+E_nu;
    E_X1=(V(i)+V(i+1))*(C_1/alpha-dt*C_2/2);
    sgX1=(V(i)+V(i+1))*(sigma^2*C_1/alpha^3+sigma^2*dt*C_2/(2*alpha^2)-(sigma^2)*(dt^2)*C_1*C_2/(2*alpha));
    E_Ic=E_X1+E_X2+E_nu*E_Z;
    Var_Ic=sgX1+sgX2+E_nu*sgZ+(E_nu2-E_nu^2)*(E_Z^2);
    Ic(i)=IG(E_Ic,(E_Ic^3)/Var_Ic);
end