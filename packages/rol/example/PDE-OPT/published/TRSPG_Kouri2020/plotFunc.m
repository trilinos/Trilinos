
function plotFunc

n = 100;

x0 = round(rand(n,1));% max(0,min(1,randn(n,1) + 1));
while (sum(x0)>n*0.5) 
  x0 = round(rand(n,1));% max(0,min(1,randn(n,1) + 1));
end
%x0 = 0.5*x0/sum(x0);
x  = x0 + 500*randn(n,1);
l  = zeros(n,1);
u  = ones(n,1);
A  = ones(1,n)/n;
b  = 0.5;
H  = speye(n);
g  = x;
del = 6;

m  = 1000;
T  = linspace(0,1,m);
F  = zeros(m,1);
for i=1:m
  F(i) = func(T(i),x,x0,l,u,A,b,H)-del;
end

[T0,T1] = regula_falsi(@(t)(func(t,x,x0,l,u,A,b,H)-del),0,1);
[S0,S1] = ridders(@(t)(func(t,x,x0,l,u,A,b,H)-del),0,1);
[U0,U1] = brents(@(t)(func(t,x,x0,l,u,A,b,H)-del),0,1);

figure,
plot(T,F,'b','linewidth',3), hold on
plot(T,0*T,'r--','linewidth',3)
plot(T1,0*T1,'kx','MarkerSize',10,'LineWidth',3), hold off
set(gca,'FontSize',15)
title('Regula Falsi','FontSize',20)
print('-depsc2','regula-falsi.eps');

pause

figure,
plot(T,F,'b','linewidth',3), hold on
plot(T,0*T,'r--','linewidth',3)
plot(S1,0*S1,'kx','MarkerSize',10,'LineWidth',3), hold off
set(gca,'FontSize',15)
title('Ridders','FontSize',20)
print('-depsc2','ridders.eps');

pause

figure,
plot(T,F,'b','linewidth',3), hold on
plot(T,0*T,'r--','linewidth',3)
plot(U1,0*U1,'kx','MarkerSize',10,'LineWidth',3), hold off
set(gca,'FontSize',15)
title('Brents','FontSize',20)
print('-depsc2','brents.eps');

fprintf('Regula Falsi Iters: %d\n',length(T0))
fprintf('Ridders Iters:      %d\n',length(S0))
fprintf('Brents Iters:       %d\n',length(U0))

end

function f = func(t,x,x0,l,u,A,b,H)

zt = t*x + (1-t)*x0;
z  = quadprog(H,-zt,A,b,[],[],l,u,zt);
f  = norm(z - x0);

end

function [T0,T1] = regula_falsi(f,t0,t1)

f0 = f(t0);
f1 = f(t1);
T0 = [t0];
T1 = [t1];

for i=1:40
  tnew = (t0*f1-t1*f0)/(f1-f0);
  fnew = f(tnew);
  if (fnew < 0)
    t0 = tnew;
    f0 = fnew;
  else
    t1 = tnew;
    f1 = fnew;
  end
  T0 = [T0;t0];
  T1 = [T1;t1];
end

end

function [T0,T1] = ridders(f,t0,t1)

f0 = f(t0);
f1 = f(t1);
T0 = [t0];
T1 = [t1];

for i=1:40
  tmid = 0.5*(t0+t1);
  fmid = f(tmid);
  if (fmid == 0), break; end

  tnew = tmid-(tmid-t0)*fmid/sqrt(fmid*fmid-f0*f1);
  fnew = f(tnew);
  if (fnew == 0) break;
  elseif (fnew < 0)
    if (fmid < 0)
      if (tnew < tmid)
        t0 = tmid;
        f0 = fmid;
      else
        t0 = tnew;
        f0 = fnew;
      end
    else
      t0 = tnew;
      f0 = fnew;
      t1 = tmid;
      f1 = fmid;
    end
  else
    if (fmid > 0)
      if (tnew > tmid)
        t1 = tmid;
        f1 = fmid;
      else
        t1 = tnew;
        f1 = fnew;
      end
    else
     t0 = tmid;
     f0 = fmid
     t1 = tnew;
     f1 = fnew;
    end
  end
  T0 = [T0;t0];
  T1 = [T1;t1];
end

end


function [T0,T1] = brents(f,t0,t1)
tol0 = sqrt(eps)*1e-2;

f0 = f(t0);
f1 = f(t1);
T0 = [t0];
T1 = [t1];

tc = t0;
fc = f0;
e  = t1 - t0;
d  = e;

while(true)
  if (abs(fc) < abs(f1))
    t0 = t1; t1 = tc; tc = t0;
    f0 = f1; f1 = fc; fc = f0;
  end

  tol = 2*eps*abs(t1) + tol0;
  m = 0.5*(tc - t1);

  if (abs(m) <= tol || f1 == 0) break; end

  if (abs(e) < tol || abs(f0) <= abs(f1))
    e = m; d = e;
  else
    s = f1 / f0;
    if (t0 == tc)
      p = 2 * m * s; q = 1 - s;
    else
      q = f0 / fc; r = f1 / fc;
      p = s * (2 * m * q * (q - r) - (t1 - t0) * (r - 1));
      q = (q -  1) * (r - 1) * (s - 1);
    end

    if (p > 0)
      q = -q;
    else
      p = -p;
    end

    s = e; e = d;

    if (2 * p < 3 * m * q - abs(tol *  q) && p < abs(0.5 * s * q))
      d = p / q;
    else
      e = m; d = e;
    end
  end

  t0 = t1;
  f0 = f1;

  if (abs(d) > tol)
    t1 = t1 + d;
  elseif (m > 0)
    t1 = t1 + tol;
  else
    t1 = t1 - tol;
  end

  f1 = f(t1);

  if ((f1 > 0 && fc > 0) || (f1 <= 0 && fc <=0))
    tc = t0;
    fc = f0;
    e  = t1 - t0;
    d  = e;
  end
  T0 = [T0;t0];
  T1 = [T1;t1];
end
end

