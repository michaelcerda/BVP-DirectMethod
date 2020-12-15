
% Solving BVP (Boundary Value Problem) using Direct method

clear all;  clc; close all;
x0=0; xf=2; n=50; 

x = linspace(x0, xf, n);
dx = x(2)-x(1);

%% Central differencing
%   (2+a(x)*dx)*T(j+1) + (2-a(x)*dx)*T(j-1)+ (b(x)*2(dt^2)-4)*T(j) = 2*f(x)*dx^2  <---- Central difference Eq.
%   4*T(N-1) + (2*b(x)*(dt^2)-4)*T(N) = 2*f(x)*dx^2  <---- Boundary Condition
% set up the AX=B equations
A3=zeros(n);
B3=zeros(n,1);
a = @(x) -((x+3)./(x+1));
b = @(x) (x+3)./((x+1).^2);
f = @(x) 2.*(x+1)+ 3.*b(x);

A3(1,1)=1;  B3(1)=5;   % <--- x1 = 0
for j=2:n-1
    A3(j,j-1)= 2-a(x(j))*dx;
    A3(j,j) = b(x(j))*2*(dx^2)-4;
    A3(j,j+1) = 2+a(x(j))*dx;
    B3(j)= 2*f(x(j))*dx^2;
end
A3(n,n-1) = 4;                                       %<------Boundary Condition
A3(n,n)= 2*b(x(n))*dx^2 -4; B3(n)= 2*f(x(n))*dx^2;   %<------Boundary Condition
 
X3 = A3\B3;


%% Plot Results of Finite Difference

% plot(x,X3,'-*')
% legend('Central difference')
X3 = A3\B3;

a = @(x) -((x+3)/(x+1));
b = @(x) (x+3)/((x+1)^2);
f = @(x) 2*(x+1)+ 3*b(x);

f1 = @(x,z) [z(2); f(x)-a(x)*z(2)- b(x)*z(1)];
x0=0; xf=2; n=50; 

% solve the ODEs with a guess for T'(0)

xp01=-0.1;  % guess a value for x'(0)
xg0=[5;xp01];
[x1,T1] = SysODE_RK4th_subinterval(f1,n,x0,xf,xg0); 
% T1 is solution for the case with x'(0)=xp01 

% plot(x1,T1(1,:), '+')
% hold on;


% solve the ODEs with another guess for T'(0)

xp02=0.1;  % guess another value for x'(0)
xg1=[5;xp02];
[x2,T2] = SysODE_RK4th_subinterval(f1,n,x0,xf,xg1);
% x2 is solution for the case with x'(0)=xp02 

% plot(x2,T2(1,:),'x')
% hold on

% calculate c1 and c2

T1L=T1(2,end); T2L=T2(2,end);
TL = 0;  %the boundary condition at x=L

c1=(TL-T2L)/(T1L-T2L);
c2=(T1L-TL)/(T1L-T2L);

%% Plot Comparison
T3 = c1*T1(1,:)+c2*T2(1,:);
plot(x,X3,'*',x1,T3,'-')
legend('Direct Method-Insulated','Shooting Method-Insulated')
xlabel('x');ylabel('T');
% % the actual x'(0) 
% xp01*c1+xp02*c2

