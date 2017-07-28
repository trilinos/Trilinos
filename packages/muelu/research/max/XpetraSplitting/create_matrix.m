function [nodes, A] = create_matrix(nregion_nodes_x, nregion_nodes_y, nintervals_x, nintervals_y)

    [nodes, nodesToRegion, nodesOwnership, nodes_neighbourhood] = create_regions( nregion_nodes_x, nregion_nodes_y, nintervals_x, nintervals_y );
    
    nTotal_nodes = size(nodesToRegion,1);
    
    A = sparse(nTotal_nodes, nTotal_nodes);
    
    for i = 1: nTotal_nodes
        
        cols = nodes_neighbourhood{i};
        
        for j = 1:length(cols)
            
            A(i, cols(j)) = nodesOwnership(cols(j));
            
        end
        
    end
    
end
