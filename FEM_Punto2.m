function [U, Time, A ,xr, V,D] = FEM_Punto2(u,n,x,t,T,dt)
    syms w z positive integer
    f = -diff(diff(u,x),x) + u + diff(u,t);
    u0 = subs(u,t,0);
    dx = 1/(2*n);
    N = 2*n-1;
    H = dx:dx:1-dx;
    y = (x-w.*dx)/(2*dx);
    v = piecewise(mod(z,2)==0,piecewise((-1<=subs(y,w,z)) & (subs(y,w,z)<0),...
        (1+subs(y,w,z))*(1+2*subs(y,w,z)),(0<=subs(y,w,z))&(subs(y,w,z)<=1),...
        (1-subs(y,w,z))*(1-2*subs(y,w,z)),0), (mod(z,2)==1),...
        piecewise((-1/2 <= subs(y,w,z))&(subs(y,w,z)<=1/2),...
        1-4*subs(y,w,z).^2,0));
    I = 1:N;
    
    V = subs(v,z,I);
    D = diff(V,x);
    
    R = zeros(N,N);
    M = zeros(N,N);
    F = zeros(N,1);
    F = sym(F);
    for k = 1:N
        F(k) = int(f*V(k),x,0,1);
        for j = k:min(k+2,N)
            R(k,j) = eval(int(D(k)*D(j),x,0,1));
            R(j,k) = R(k,j);
            M(k,j) = eval(int(V(k)*V(j),x,0,1));
            M(j,k) = M(k,j);
        end
    end

    U = {u0};
    a0 = subs(u0,x,H)';
    A(1,:) = a0';
    Time = [0];
    c = 2;
    tp=dt;
    
    K = M/dt + R + M;

    while tp<=T
        b = eval(subs(F,t,tp))+1/dt*M*a0;
        try
            xr = K\b;
        catch
            [TG,CG,xr,count,E] = GSeidel(K,b,100000,1e-6); 
        end
        A(c,:) = xr';
        Time = [Time tp]; 
        tp = tp+dt;   
        uh = V*xr;
        U{c} = uh;
    
        c=c+1;
        a0 = xr;
    end
    %Err=sqrt(eval(int((u-Uh).^2,0,1)));
end