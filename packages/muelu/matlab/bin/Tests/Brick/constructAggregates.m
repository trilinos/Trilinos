function agg = constructAggregates(nVertices, nAggregates, vertexToAggID, rootNodes, aggSizes)
  agg = struct('nVertices', nVertices, 'nAggregates', nAggregates, 'vertexToAggID', vertexToAggID, 'rootNodes', rootNodes, 'aggSizes', aggSizes);
end
