%Prototype function to test evolution strength-of-connection
%
% reference: Algorithm 1, "Evolution Measure", p. 724 of
%            "A new perspective on strength measures in algebraic multigrid"
%            Luke Olson, Jacob Schroder, and Ray Tuminaro
%            Numerical Linear Algebra with Applications, Vol. 17, p 713-733, 2010
%            doi 10.1002/nla.669
%
% Notes:  Dofs per node is hardwired to 1
%
function [Graph, DofsPerNode, Filtering, S] = evolutionSoC(A, B)

  %TODO Should not be hardwired to 1.
  DofsPerNode = 1;
  % TODO get k from Parameter List
  % TODO get theta from Parameter List
  % TODO use lambda max estimate if it's available

  Graph.edges = logical(A);
  Graph.boundaryNodes = int32([]);

  % temporary matrix recording SoC measures
  S = sparse(A);

  D = diag(A);
  DinvA = diag(D)\A;
  lambdaMax = eigs(DinvA,1,'LM');
  if exist('k') ~= 1
    k = max(floor(lambdaMax),1);
  end
  if exist('theta') ~= 1
    theta = 4.0;
  end

  Filtering = false; % If any connections are dropped, then Filtering must be true
  [m,n] = size(A);
  I = speye(m,n);
  deltaFn = zeros(m,1);
  e = zeros(m,1);

  % iterate over rows of A
  numDropped = 0;
  for ii=1:m
    e(ii) = 1;
    deltaFn(ii) = 1;
    z = deltaFn;  % FIXME deep copy
    for jj=1:k
      %z = (I - 1/(k*lambdaMax)*DinvA)*z;
      z = z - 1/(k*lambdaMax)*(DinvA*z);
    end

    % Omega is a mock aggregate that contains all possible connections (due to A's sparsity pattern) to root ii
    Omega = find(A(ii,:));
    % record index of root node
    rootIndex = find(Omega==ii);
    %solve constrained minimization problem
    Q = diag(D(Omega,:));
    Bbar = B(Omega,:);
    ebar = e(Omega,:);
    zbar = z(Omega,:);
    q = Q(rootIndex,rootIndex);
    dd=size(ebar,2);
    ee=size(Bbar,2);
    blockA = [2*Bbar'*Q*Bbar   q*Bbar'*ebar     ;...  %FIXME : missing Q(ii,ii) from 2,2 term
                ebar'* Bbar    zeros(dd,ee)     ];
    blockRhs = [2*Bbar'*Q*zbar ; ...
              zbar(rootIndex)];
    
    x = blockA \ blockRhs;
    ztilde = Bbar * x(1:end-1);

    %calculate SoC measure for each nonroot point in mock aggregate (each off-diagonal in row ii's stencil)
    WEAK=intmax;
    soc = zeros(size(Omega));
    for kk=1:length(Omega)
      if ztilde(kk)/zbar(kk) < 0
        soc(kk) = WEAK;
      else
        soc(kk) = abs( (zbar(kk) - ztilde(kk)) / zbar(kk) );
      end
    end

    %find minimum over all nonroot points
    soc(rootIndex) = WEAK;
    minOverOmega = min(soc);
    %check SoC for each nonroot point
    nnz=length(Omega);
    for kk=1:length(Omega)
      S(ii,Omega(kk)) = soc(kk);
      if (soc(kk) > theta * minOverOmega) && (Omega(kk) ~= ii)
        %Omega(kk)'s connection to ii is weak, but never drop diagonal
        Graph.edges(ii,Omega(kk)) = false;
        numDropped = numDropped + 1;
        Filtering = true;
        nnz = nnz-1;
      else
        Graph.edges(ii,Omega(kk)) = true;
      end
    end
    S(ii,ii) = 1;
    if nnz == 1
      Graph.boundaryNodes(end+1) = ii;
    end
    e(ii) = 0;
    deltaFn(ii) = 0;
  end %for ii=1:m

  fprintf('Evolution Soc: dropped %d connections\n', numDropped);
  fprintf('Evolution Soc: detected %d boundary nodes\n', length(Graph.boundaryNodes));
  
end 
