
function flowMovie(name, type, T, xlimit, ylimit)
%% Parse data from input files
tic
fprintf('Parsing input data.\n')
[Ux, Uy, P, adj, nodes] = parseData(name);
toc
%% Compute quantity to plot
tic
fprintf('Computing quantity for plotting.\n')
if (strcmp(type,'mag'))
  data = sqrt(Ux.^2 + Uy.^2);
elseif (strcmp(type,'vort'))
  data = vorticity(nodes, Ux, Uy);
elseif (strcmp(type,'velx'))
  data = Ux;
elseif (strcmp(type,'vely'))
  data = Uy;
elseif (strcmp(type,'pres'))
  data = P;
else
  fprintf('Unknown type %s -- Terminating.\n', type);
  return;
end
toc
%% Capture movie
tic
fprintf('Capturing movie.\n')
trisurfMovie(adj, nodes, data, [name,'_',type], T, xlimit, ylimit);
toc
