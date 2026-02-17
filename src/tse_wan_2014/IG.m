%generation of IG, Algorithm 3 of the paper
function ris=IG(m,s)  
% parameter: m is the mean, s the shape 
ris=zeros(size(m)); 
% generating a vectors N,U from multivariate normal and uniform distribution  
N=randn(size(m));  
U=rand(size(m));
x=1+N.^2./(2.*s./m)-sqrt((4.*s./m).*N.^2+(N.^4))./(2.*s./m);
mask=((U.*(1+x))>1);
ris(mask)=m(mask)./x(mask);
ris(~mask)=m(~mask).*x(~mask);
end