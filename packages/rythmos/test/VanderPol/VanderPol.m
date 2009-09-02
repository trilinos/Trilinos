%
% Exact discrete solution to Van der Pol equation with Backward Euler
%

function [x0n1,x1n1] = VanderPol(N)
format long e;
% IC:
x0n = 2.0;
x1n = 0.0;
x0n1 = [x0n]; % IC
x1n1 = [x1n]; % IC
epsilon = 0.5;
h = 0.1;
for i=[1:N]
  t = 0.0 + i*h;
  fn = forcing_term(t,epsilon);
  [x0n1_temp,x1n1_temp] = solveCubic(x0n,x1n,epsilon,h,fn);
  x0n1 = [x0n1,x0n1_temp];
  x1n1 = [x1n1,x1n1_temp];
  x0n = x0n1_temp;
  x1n = x1n1_temp;
end;



function [f] = forcing_term(t,epsilon)
x0 = 2*cos(t)+epsilon*(0.75*sin(t)-0.25*sin(3*t));
x1 = -2*sin(t)+epsilon*(0.75*cos(t)-0.75*cos(3*t));
f = -2*cos(t)+epsilon*(-0.75*sin(t)+(9/4)*sin(3*t))-epsilon*(1-x0*x0)*x1+x0;

function [x0n1,x1n1] = solveCubic(x0n,x1n,epsilon,h,fn)
a3 = epsilon/h;
a2 = -epsilon/h*x0n;
a1 = 1/(h*h)+1-epsilon/h;
a0 = epsilon/h*x0n-(1/(h*h))*x0n-(1/h)*x1n-fn;
poly=[a3,a2,a1,a0];
R = roots(poly);
x0n1 = R(3);
x1n1 = (x0n1-x0n)/h;

