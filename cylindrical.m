function [t,h,ddh] = cylindrical(mu,L,H0,H1,b,d,t0)
%
% Function 'cylindrical' returns the motion of COM on a cylindrical ramp
%
% Input arguments:
% mu: friction coefficient
% L:  radius of pendulum
% H0: min of COM 
% H1: max of COM
% b: pumping momentum
% d: timing parameter
%
% Output arguments:
% t: angle 
% h: height of COM
% ddh: acceleration of COM
% 
% Quick usage:
% >> [t,h,ddh] = cylindrical(0.02,3,0.8,0.9,14.74,-0.025,-pi/3);
% 

p.mu = mu;
p.L  = L;
p.H0 = H0;
p.H1 = H1;
p.b  = b;
p.d  = d;

dt = 1e-3;
t = t0:dt:pi;
K = kinetic(t,0,p);
i = find(K<0,1);
if isempty(i)
    t1 = pi;
else
    t1 = interp1([K(i-1),K(i)],[t(i-1),t(i)],0);
end

dt = (t1-t0)*1e-3;
t = t0:dt:t1;
h = h_sgmd(t,p);
h1 = dh_sgmd(t,p);
h2 = ddh_sgmd(t,p);
K = kinetic(t,0,p);
ddh = 2*h2.*K + h1.*(f_frc(t,p).*K+k_frc(t,p));

subplot(2,1,1)
plot(t,h,'LineWidth',1);
ylabel('height');
subplot(2,1,2)
plot(t,ddh,'LineWidth',1);
ylabel('acceleration');


function K = kinetic(x,K0,p)
n = length(x);
K = zeros(1,n-1);
K(1) = expf(x(1),x(2),p)*K0+A(x(1),x(2),p);
for i=1:n-2
    K(i+1) = K(i)*expf(x(i+1),x(i+2),p) + A(x(i+1),x(i+2),p);
end
K = [K0,K];

function y = A(t1,t2,p)
y = (itgrd(t1,t2,p)+4*itgrd((t1+t2)/2,t2,p)+itgrd(t2,t2,p))*(t2-t1)/6;

function y = itgrd(x,b,p)
y = k_frc(x,p).*expf(x,b,p);

function y = expf(a,b,p)
z = integral(@(x)f_frc(x,p),a,b);
y = exp(z);

function y = f_frc(x,p)
r0 = p.L-h_sgmd(x,p);
r1 = -dh_sgmd(x,p);
r2 = -ddh_sgmd(x,p);
y = -(4*r1 + 2*p.mu*r0 -2*p.mu*r2)./(r0-p.mu*r1);

function y = k_frc(x,p)
r0 = p.L-h_sgmd(x,p);
r1 = -dh_sgmd(x,p);
y = -9.8*(sin(x) + p.mu*cos(x))./(r0-p.mu*r1);

function y = h_sgmd(x,p)
z = exp(-p.b*(x-p.d));
y = (p.H1-p.H0)./(1+z) + p.H0;

function y = dh_sgmd(x,p)
z = exp(-p.b*(x-p.d));
y = (p.H1-p.H0)*p.b*z./(1+z).^2;

function y = ddh_sgmd(x,p)
z = exp(-p.b*(x-p.d));
y = -(p.H1-p.H0)*p.b^2*z.*(1-z)./(1+z).^3;