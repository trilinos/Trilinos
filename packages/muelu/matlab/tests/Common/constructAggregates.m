function agg = constructAggregates(nVertices, nAggregates, vertexToAggID, rootNodes, aggSizes)
  agg = struct('nVertices', int32(nVertices), 'nAggregates', int32(nAggregates), 'vertexToAggID', int32(vertexToAggID), 'rootNodes', int32(rootNodes), 'aggSizes', int32(aggSizes));
end
