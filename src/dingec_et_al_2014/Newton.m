%Newton's method taken from Matematica Numerica (book)
function zero=Newton(fun,Dfun,x0,tol,nmax)
%pag228 of the book Matematica Numerica
err=tol+1;
iter=0;
fx0=fun(x0);
while iter<nmax && err>tol
    dfx0=Dfun(x0);
    if dfx0==0
        fprintf('Stop because the derivative is zero\n');
        zero=[];
        return
    end
    x=x0-fx0/dfx0;
    err=abs(x-x0);x0=x;fx0=fun(x0);iter=iter+1;
end
zero=x0;
end