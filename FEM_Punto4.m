function [uh , Err, A,b,xr, V,D,DD] = FEM_Punto4(u,n,x)
    syms w z r positive integer
    f = diff(diff(diff(u)))-diff(diff(u));
    dx = 1/(2*n);
    N = 2*n-1;
    H = 0:dx:1;
    y = (x-w.*dx)/(2*dx);
    phi = piecewise((-1<=subs(y,w,z)) & (subs(y,w,z)<0),...
        (1+subs(y,w,z)).^2*(1-2*subs(y,w,z)),(0<=subs(y,w,z))&(subs(y,w,z)<=1),...
        (1-subs(y,w,z)).^2*(1+2*subs(y,w,z)),0);
    psi = piecewise(z==0,0,z==N+1,0,piecewise((-1 <= subs(y,w,z))&(subs(y,w,z)<=0),...
        subs(y,w,z).*(1+subs(y,w,z)).^2,(0 < subs(y,w,z))&...
        (subs(y,w,z)<=1),subs(y,w,z).*(1-subs(y,w,z)).^2,0));

    v = piecewise(mod(r,2)==0,subs(phi,z,r),mod(r,2)==1,...
        -1/(dx) * (subs(psi,z,r+1)- subs(psi,z,r-1)));
    I = 1:N;

    V = subs(v,r,I);
    D = diff(V,x);
    DD = diff(D,x);
      
    A = zeros(N,N);
    b = zeros(N,1);
    for k = 1:N
        b(k) = eval(int(f*V(k),0,1));
        for j = max(1,k-4):min(N,k+4)
            A(k,j) = eval(int((DD(k)+D(k))*D(j),0,1));
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