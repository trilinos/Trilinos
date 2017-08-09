function [global_indices, nodesToRegion, nodesOwnership, nodes_neighbourhood] = create_regions( nregion_nodes_x, nregion_nodes_y, nintervals_x, nintervals_y )

%	INPUT: 
%		nregions_nodes_x : number of nodes along the x direction for each region	
%		nregions_nodes_y : number of nodes along the y direction for each region	
%		nintervals_x : number of intervals splitting the domain along the x direction
%		nintervals_y : number of intervals splitting the domain along the y direction
%
%	OUTPUT:
%		global_indices: global node labels with region index associated with it
%		nodesToRegion: cell structure that contains the list of regions a node belongs to
% 		nodesOwnership: structure that associates a globally idnexed node with the region that owns it
% 				(even if a node sits in more than one region, only a single region can own it)
%		nodes_neighbourhood: for each node, this structure provides the list of neighbors (itself included) with distance equal to 1

%	number of nodes along the x direction of the main
    ntotal_nodes_x = nregion_nodes_x * nintervals_x - (nintervals_x-1);

%	number of nodes along the y direction of the domain
    ntotal_nodes_y = nregion_nodes_y * nintervals_y - (nintervals_y-1);

%	total number of nodes that discretize the domain
    ntotal_nodes = ntotal_nodes_x * ntotal_nodes_y;
       
    % detects the global indices ONLY for nodes in region 1
    global_indices_region1 = zeros(nregion_nodes_x*nregion_nodes_y,1);
    
    for i = 1 : nregion_nodes_x
        
        for j = 1 : nregion_nodes_y
           
            global_indices_region1( (i-1)*nregion_nodes_y + j ) = j + (i-1) * ntotal_nodes_y;
            
        end
        
    end
  

%   use the global indices of nodes in region 1 to reconstruct the global indices of the nodes
%   in every other region  
    global_indices = [];
    
    for i = 1 : nintervals_x
        
        for j = 1 : nintervals_y
           
            region_index = (i-1)*nintervals_y + j;
            
            global_indices = [global_indices; [global_indices_region1 + (i-1) * (nregion_nodes_x-1) * ntotal_nodes_y + (j-1) * (nregion_nodes_y-1), region_index * ones(nregion_nodes_x*nregion_nodes_y,1) ] ];
            
        end
        
    end
    
    nodesToRegion = cell(ntotal_nodes, 1);
    
    parfor i = 1 : ntotal_nodes
       
        find_result = global_indices(:,1)==i;
        regions = global_indices(find_result,2);
        nodesToRegion{i} = regions;
        
    end
    
    nodesOwnership = zeros(ntotal_nodes, 1);
    
%     for i = 1 : ntotal_nodes
%        
%         nodesOwnership(i) = min(nodesToRegion{i});
%         
%     end
    
    nodes_neighbourhood = cell(ntotal_nodes, 1);
    
    parfor i = 1 : ntotal_nodes
        
        neighbourhood = [];
        regions =  nodesToRegion{i};
        
        for region = 1 : length(regions)
                   
            find_result = global_indices(:,2)==regions(region);
            regional_nodes = global_indices(find_result,1);

            neighbours = [ i-ntotal_nodes_y, i-ntotal_nodes_y-1, i-ntotal_nodes_y+1, i, i-1, i+1, i+ntotal_nodes_y, i+ntotal_nodes_y-1, i+ntotal_nodes_y+1 ];
            neighbours = intersect( neighbours, regional_nodes );

            if(~isempty(neighbours))
                    neighbourhood = [neighbourhood; neighbours];
            end

        end
        
        nodes_neighbourhood{i} = unique(neighbourhood);
        
    end

end







