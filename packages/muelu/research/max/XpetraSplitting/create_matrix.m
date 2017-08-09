function [nodes, A] = create_matrix(nregion_nodes_x, nregion_nodes_y, nintervals_x, nintervals_y)

%	INPUT: 
%		nregions_nodes_x : number of nodes along the x direction for each region	
%		nregions_nodes_y : number of nodes along the y direction for each region	
%		nintervals_x : number of intervals splitting the domain along the x direction
%		nintervals_y : number of intervals splitting the domain along the y direction
%
%	OUTPUT:
%		nodes: global node labels with region index associated with it
%		A: composite matrix
% 		
%		A(i,j) = k if node j belongs to regions to region k

    [nodes, nodesToRegion, ~, nodes_neighbourhood] = create_regions( nregion_nodes_x, nregion_nodes_y, nintervals_x, nintervals_y );
    
    nTotal_nodes = size(nodesToRegion,1);
    
    A = sparse(nTotal_nodes, nTotal_nodes);
    
    for i = 1: nTotal_nodes
        
        cols = nodes_neighbourhood{i};
        
        for j = 1:length(cols)
            
            %A(i, cols(j)) = nodesOwnership(cols(j));
            if( i==cols(j) )
                A(i, cols(j)) = 8;
            else
                A(i, cols(j)) = -1;
            end
            
        end
        
    end
    
end
