%A function implementing a very simple aggregation scheme.
%Triplets of nodes with consecutive IDs are grouped together.
%Should simulate brick with some set of parameters?
function agg = simpleAggregation(A)
        [m, n] = size(A);
	nVerts = m; %number of rows -> number of nodes
	nAggs = nVerts / 3;
	vertToAgg = int32(zeros([nVerts, 1]));
	for i = 1:nVerts
		vertToAgg(i) = int32((i - 2) / 3);
	end
	rootNodes = int32(zeros([nAggs, 1]));
	for i = 1:nAggs
		rootNodes(i) = int32(i * 3 - 2);
	end
	aggSizes = int32(zeros([nAggs, 1]));
	for i = 1:nAggs
		aggSizes(i) = int32(3);
	end
	agg = constructAggregates(nVerts, nAggs, vertToAgg, rootNodes, aggSizes);
end
