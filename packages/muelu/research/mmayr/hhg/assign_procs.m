function [ nodesToProc, nodes ] = assign_procs(nodesToRegion, nodes)

  nNodes = size(nodesToRegion,1);
  
  nodesToProc = cell(nNodes,1);
  
  parfor i = 1 : nNodes
    regions = nodesToRegion{i};
    region = regions(1);
    nodesToProc{i} = region;
  end
  
  procs = zeros(size(nodes,1),1);
  for i = 1 : nNodes
    procs(find(nodes(:,1) == i)) = nodesToProc{i};
  end
  nodes = [nodes procs];
  
end
