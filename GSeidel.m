function [T,C,x0,count,E] = GSeidel(A,b,k,tol)
    D=diag(diag(A));
    L=-tril(A,-1);
    U=-(A-D+L);
    T=inv(D-L)*(U);
    C=inv(D-L)*b;
    x0=zeros(size(A,1),1);
    count=0;
    E=1000;
    while E>tol && count<k   
        xn= T*x0+C;
        E=norm(x0-xn);
        x0=xn;
        count=count+1;
    end
end