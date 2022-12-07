
function sample(M, full, useQMC)
if nargin < 3
  useQMC = true;
end
if nargin < 2
  full = true;
end

M2 = floor(M/2);

r  = 1.5;
R  = 0.5;

% Vertices of the split octagon
a = 2*r/(1+sqrt(2));
vx1 = [-0.5*(r+0.5*a); -0.5*a; 0.5*a; r; r; 0.5*(r+0.5*a)];
vy1 = [-0.5*(r+0.5*a); -r; -r; -0.5*a; 0.5*a; 0.5*(r+0.5*a)];
vx2 = [0.5*(r+0.5*a); 0.5*a; -0.5*a; -r; -r; -0.5*(r+0.5*a)];
vy2 = [0.5*(r+0.5*a); r; r; 0.5*a; -0.5*a; -0.5*(r+0.5*a)];

% Use Quasi Monte Carlo to sample
S0 = 2*r*rand(20*M,2)-r;
if (useQMC)
  rng default
  skip = 1e3;
  leap = 1e2;
  if (full)
    skip = 5e3;
    leap = 1e3;
  end
  p  = haltonset(2,'Skip',skip,'Leap',leap);
  p  = scramble(p,'RR2');
  S0 = 2*r*net(p,20*M)-r;
end

in = inpolygon(S0(:,1),S0(:,2),vx1,vy1);
S  = S0(in,:);
if (full) 
  S = S(1:M2,:);
else
  [m,n] = size(S);
  S0    = [];
  for (i = 1:m)
    xnorm = norm(S(i,:));
    if (xnorm > 0.5)
      S0 = [S0;S(i,:)];
    end
    if (length(S0(:,1)) == M2)
      break;
    end
  end
  S = S0;
end

S = [S;[S(:,2),S(:,1)]];

%% Vertices of octagon
%a  = 2*r/(1+sqrt(2));
%vx = [-0.5*a; 0.5*a; r; r; 0.5*a; -0.5*a; -r; -r];
%vy = [-r; -r; -0.5*a; 0.5*a; r; r; 0.5*a; -0.5*a];
%
%% Use Quasi Monte Carlo to sample
%rng default
%skip = 1e3;
%leap = 1e2;
%if (full)
%  skip = 5e3;
%  leap = 1e3;
%end
%p  = haltonset(2,'Skip',skip,'Leap',leap);
%p  = scramble(p,'RR2');
%S0 = 2*r*net(p,20*M)-r;
%
%in = inpolygon(S0(:,1),S0(:,2),vx,vy);
%S  = S0(in,:);
%if (full) 
%  S = S(1:M,:);
%else
%  [m,n] = size(S);
%  S0    = [];
%  for (i = 1:m)
%    xnorm = norm(S(i,:));
%    if (xnorm > 0.5)
%      S0 = [S0;S(i,:)];
%    end
%    if (length(S0(:,1)) == M)
%      break;
%    end
%  end
%  S = S0;
%end

[m,n] = size(S);
if (m ~= M)
  fprintf('There are not enough samples in S!\n')
end

figure, plot(S(:,1),S(:,2),'.')
axis('equal','square')

pname = ['hole_points_',int2str(M),'.txt'];
wname = ['hole_weights_',int2str(M),'.txt'];
if (full)
  pname = ['points_',int2str(M),'.txt'];
  wname = ['weights_',int2str(M),'.txt'];
end

fileP = fopen(pname,'w');
fileW = fopen(wname,'w');
for (i = 1:M)
  fprintf(fileP,'%1.16f     %1.16f\n',S(i,1),S(i,2));
  fprintf(fileW,'%1.16f\n',1/M);
end
fclose(fileP);
fclose(fileW);

return
