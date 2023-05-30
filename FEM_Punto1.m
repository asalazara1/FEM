function [uh , Err, A,b,xr, V,D] = FEM_Punto1(u,n,x)
    syms w z positive integer
    f = -diff(diff(u));
    dx = 1/(2*n);
    N = 2*n-1;
    H = 0:dx:1;
    y = (x-w.*dx)/(2*dx);
    v = piecewise(mod(z,2)==0,piecewise((-1<=subs(y,w,z)) & (subs(y,w,z)<0),...
        (1+subs(y,w,z))*(1+2*subs(y,w,z)),(0<=subs(y,w,z))&(subs(y,w,z)<=1),...
        (1-subs(y,w,z))*(1-2*subs(y,w,z)),0), (mod(z,2)==1),...
        piecewise((-1/2 <= subs(y,w,z))&(subs(y,w,z)<=1/2),...
        1-4*subs(y,w,z).^2,0));
    I = 1:N;

    V = subs(v,z,I);
    D = diff(V,x);
    
      
    A = zeros(N,N);
    b = zeros(N,1);
    for k = 1:N
        b(k) = eval(int(f*V(k),0,1));
        for j = k:min(k+2,N)
            A(k,j) = eval(int(D(k)*D(j),0,1));
            A(j,k) = A(k,j);
        end
    end
    
    try
        xr = A\b;
    catch
        [T,C,xr,count,E] = GSeidel(A,b,100000,1e-6); 
    end
   
    uh=V*xr;
    
    Err=sqrt(eval(int((u-uh).^2,0,1)));
end