function [global_indices, nodesToRegion, nodesOwnership, nodes_neighbourhood] = create_regions( nregion_nodes_x, nregion_nodes_y, nintervals_x, nintervals_y )

    initial_nodes_info = zeros( nregion_nodes_x*nregion_nodes_y*nintervals_x*nintervals_y );
    ntotal_nodes_x = nregion_nodes_x * nintervals_x - (nintervals_x-1);
    ntotal_nodes_y = nregion_nodes_y * nintervals_y - (nintervals_y-1);
    ntotal_nodes = ntotal_nodes_x * ntotal_nodes_y;
       
    % detects the global indices of region 1
    global_indices_region1 = zeros(nregion_nodes_x*nregion_nodes_y,1);
    
    node_idex = 1;
    
    for i = 1 : nregion_nodes_x
        
        for j = 1 : nregion_nodes_y
           
            global_indices_region1( (i-1)*nregion_nodes_x + j ) = j + (i-1) * ntotal_nodes_x;
            
        end
        
    end
    
    global_indices = [];
    
    for i = 1 : nintervals_x
        
        for j = 1 : nintervals_y
           
            region_index = (i-1)*nintervals_y + j;
            
            global_indices = [global_indices; [global_indices_region1 + (i-1) * (nregion_nodes_x-1) * ntotal_nodes_y + (j-1) * (nregion_nodes_y-1), region_index * ones(nregion_nodes_x*nregion_nodes_y,1) ] ];
            
        end
        
    end
    
    nodesToRegion = cell(ntotal_nodes, 1);
    
    for i = 1 : ntotal_nodes
       
        find_result = global_indices(:,1)==i;
        regions = global_indices(find_result,2);
        nodesToRegion{i} = regions;
        
    end
    
    nodesOwnership = zeros(ntotal_nodes, 1);
    
    for i = 1 : ntotal_nodes
       
        nodesOwnership(i) = min(nodesToRegion{i});
        
    end
    
    nodes_neighbourhood = cell(ntotal_nodes, 1);
    
    for i = 1 : ntotal_nodes
        
        neighbourhood = [];
        regions =  nodesToRegion{i};
        
        for region = 1 : length(regions)
                   
            find_result = global_indices(:,2)==regions(region);
            regional_nodes = global_indices(find_result,1);

            neighbours = [ i-ntotal_nodes_y, i-ntotal_nodes_y-1, i-ntotal_nodes_y+1, i, i-1, i+1, i+ntotal_nodes_y, i+ntotal_nodes_y-1, i+ntotal_nodes_y+1 ];
            neighbours = intersect( neighbours, regional_nodes );

            if(~isempty(neighbours))
                    neighbourhood = [ neighbourhood  neighbours];
            end

        end
        
        nodes_neighbourhood{i} = unique(neighbourhood);
        
    end

end







