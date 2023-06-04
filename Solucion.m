%% Punto 1
syms x t

u1 = (x.^2-x.^3);
N = [2 3 4 5 10 20 35 55 75 100];
E12 = zeros(10,1);
T1 = zeros(10,1);
U1 = zeros(10,1);
U1 = sym(U1);
tic
parfor i = 1:length(N)
    [uh , Err, A,b,xr, V,D] = FEM_Punto1(u1,N(i),x);
    E12(i) = Err;
    U1(i) = uh;
end
toc

%% Evaluacion 1
tic
E1inf = zeros(10,1);
E12 = zeros(10,1);

parfor i = 1:10
    E12(i) = L2(u1-U1(i),x,0,1);
    E1inf(i) = Linfty(u1-U1(i),1000,x);
end
toc

%% Punto 2

u2 = (x.^2-x.^3).*exp(-t);
n = 25;
T = 1;
dt = 0.05;
tic
[U2, Time, A ,xr, V,D] = FEM_Punto2(u2,n,x,t,T,dt);
T2  = toc;
T2
%% Evaluación 2
tic
E22 = [];
E2inf = [];
parfor i = 1:length(Time)
    r = subs(u2,t,Time(i))-U2{i};
    E22(i) = L2(r,x,0,1);
    E2inf(i) = Linfty(r,1000,x);
end
toc

%% Punto 3
u3 = sin(2*pi*x-1/2*pi)+1;
E3 = zeros(10,1);
T3 = zeros(10,1);
U3 = zeros(10,1);
U3 = sym(U3);

tic
parfor i = 1:length(N)
    
    [uh , Err, A,b,xr, V,D] = FEM_Punto3(u3,N(i),x);
    E3(i) = Err;
    U3(i) = uh;
end
toc
%% Evaluacion 3
tic
E3inf = zeros(10,1);
E32 = zeros(10,1);

parfor i = 1:10
    E32(i) = L2(u3-U3(i),x,0,1);
    E3inf(i) = Linfty(u3-U3(i),1000,x);
end
toc
%% Punto 4

E4 = zeros(10,1);
U4 = zeros(10,1);
U4 = sym(U4);

tic
parfor i = 1:length(N)
    
    [uh , Err, A,b,xr, V,D,DD] = FEM_Punto4(u3,N(i),x);
    E4(i) = Err;
    U4(i) = uh;
end
toc

%% Evaluacion 4

tic
E4inf = zeros(10,1);
E42 = zeros(10,1);

parfor i = 1:10
    E42(i) = L2(u3-U4(i),x,0,1);
    E4inf(i) = Linfty(u3-U4(i),1000,x);
end
toc

%% Funciones
function norm_infty = Linfty(u,n,x)
   X_prueba = 0:1/n:1;
   r = abs(subs(u,x,X_prueba));
   norm_infty = double(max(r));
end

function norm_2 = L2(u,x,linf,lsup)
    norm_2 = sqrt(double(int(u.^2,x,linf,lsup)));
    if isnan(norm_2)
        r = @(y) double(subs(u.^2,x,y));
        norm_2 = sqrt(integral(r,0,1));
        norm_2 = sqrt(real(norm_2).^2+imag(norm_2).^2);
    end
    
end