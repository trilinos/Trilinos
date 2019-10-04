#if HAVE_CUBIT

#include <percept/structured/PGeomAssocStructured.hpp>
#include <PGeom.hpp>
#include <PGeomRangeSearch.hpp>
#include "CubitVector.hpp"
#include <iostream>

const std::string PGeomAssocStructured::direction_string[3] = {{"I"}, {"J"}, {"K"}};

void PGeomAssocStructured::find_mesh_association(std::shared_ptr<PGeom> pgeom,
                                                 std::shared_ptr<percept::BlockStructuredGrid> grid,
                                                 const bool local_tol,
                                                 const double associate_tolerance,
                                                 const int print_level)
{
    mBlockAssociation.resize(grid->m_sblocks.size());

    find_node_association(pgeom, grid, local_tol, associate_tolerance);
    find_edge_association(pgeom, grid);
    find_face_association(pgeom, grid);

    if (print_level) print_block_geometry_association();
}

void PGeomAssocStructured::find_node_association(std::shared_ptr<PGeom> pgeom,
                                                 std::shared_ptr<percept::BlockStructuredGrid> grid,
                                                 bool local_tol,
                                                 double associate_tolerance)
{
  std::shared_ptr<PGeomRangeSearch> range_search(new PGeomRangeSearch());

  bool surf_tree_successful = range_search->create_surface_tree();
  bool curve_tree_successful = range_search->create_curve_tree();
  bool vertex_tree_successful = range_search->create_vertex_tree();

  if (!surf_tree_successful || !curve_tree_successful || !vertex_tree_successful)
    return;


  for (unsigned iblock=0; iblock < grid->m_sblocks.size(); ++iblock)
  {
    std::shared_ptr<percept::StructuredBlock> sgi = grid->m_sblocks[iblock];

    unsigned i_min = sgi->m_sizes.node_min[0];
    unsigned j_min = sgi->m_sizes.node_min[1];
    unsigned k_min = sgi->m_sizes.node_min[2];
    unsigned i_max = sgi->m_sizes.node_max[0];
    unsigned j_max = sgi->m_sizes.node_max[1];
    unsigned k_max = sgi->m_sizes.node_max[2];


    // find association for nodes on i_min block face
    std::array<unsigned, 3> range_min = {{i_min, j_min, k_min}};
    std::array<unsigned, 3> range_max = {{i_min, j_max, k_max}};
    create_node_association_data(pgeom, range_search, sgi, range_min, range_max, local_tol,
                                 associate_tolerance, iblock);

    // find association for nodes on i_max block face
    range_min[0] = i_max;
    range_max[0] = i_max;
    create_node_association_data(pgeom, range_search, sgi, range_min, range_max, local_tol,
                                 associate_tolerance, iblock);

    // find association for nodes on j_min block face - skip i_min and i_max
    range_min[0] = i_min + 1;
    range_min[1] = j_min;
    range_max[0] = i_max - 1;
    range_max[1] = j_min;
    create_node_association_data(pgeom, range_search, sgi, range_min, range_max, local_tol,
                                 associate_tolerance, iblock);

    // find association for nodes on j_max block face - skip i_min and i_max
    range_min[1] = j_max;
    range_max[1] = j_max;
    create_node_association_data(pgeom, range_search, sgi, range_min, range_max, local_tol,
                                 associate_tolerance, iblock);

    // find association for nodes on k_min block face - skip i_min, i_max, j_min, j_max
    range_min[1] = j_min + 1;
    range_min[2] = k_min;
    range_max[1] = j_max - 1;
    range_max[2] = k_min;
    create_node_association_data(pgeom, range_search, sgi, range_min, range_max, local_tol,
                                 associate_tolerance, iblock);

    // find association for nodes on k_min block face - skip i_min, i_max, j_min, j_max
    range_min[2] = k_max;
    range_max[2] = k_max;
    create_node_association_data(pgeom, range_search, sgi, range_min, range_max, local_tol,
                                 associate_tolerance, iblock);
  }
}


void PGeomAssocStructured::create_node_association_data(std::shared_ptr<PGeom> pgeom,
                                                        std::shared_ptr<PGeomRangeSearch> range_search,
                                                        std::shared_ptr<percept::StructuredBlock> sgi,
                                                        std::array<unsigned, 3> range_min,
                                                        std::array<unsigned, 3> range_max,
                                                        bool local_tol,
                                                        double associate_tolerance,
                                                        unsigned block_num)
{
  std::multimap< int, int > surface_to_node_index;
  std::multimap< int, int > curve_to_node_index;
  std::multimap< int, int > vertex_to_node_index;

  std::map<int, NodeAssociationData> &node_map = mBlockAssociation[block_num].node_to_geom_map;
  std::multimap<int, int> &surf_map = mBlockAssociation[block_num].surface_to_node_map;
  std::multimap<int, int> &curve_map = mBlockAssociation[block_num].curve_to_node_map;
  std::multimap<int, int> &vertex_map = mBlockAssociation[block_num].vertex_to_node_map;

  std::array<unsigned,3> indx{{0,0,0}};

  //TODO - do we need to get this data if we have local_offset use the data internally
  //     - can we get rid of the access ordering and always access in the same way
  const unsigned A0 = sgi->m_access_ordering[0], A1 = sgi->m_access_ordering[1], A2 = sgi->m_access_ordering[2];
  unsigned Asizes[3] = {sgi->m_sizes.node_size[A0], sgi->m_sizes.node_size[A1], sgi->m_sizes.node_size[A2]};


  for (indx[2] = range_min[2]; indx[2] <= range_max[2]; ++indx[2]) // K
    {
      for (indx[1] = range_min[1]; indx[1] <= range_max[1]; ++indx[1]) // J
        {
          for (indx[0] = range_min[0]; indx[0] <= range_max[0]; ++indx[0]) // I
            {
              int local_node_index = sgi->local_offset(indx[0], indx[1], indx[2], Asizes);

              NodeAssociationData nd;
              nd.edge_len_tol = associate_tolerance;
              if (local_tol)
                nd.edge_len_tol = structured_node_local_tolerance(sgi, indx, associate_tolerance);

              node_map[local_node_index] = nd;

              CubitVector node_pt(&sgi->m_sgrid_coords(indx[0], indx[1], indx[2], 0));
              CubitBox node_box(node_pt);

              std::vector<int> close_surfs, close_curves, close_vertices;

              range_search->find_close_surfs(node_box, nd.edge_len_tol, close_surfs);
              range_search->find_close_curves(node_box, nd.edge_len_tol, close_curves);
              range_search->find_close_vertices(node_box, nd.edge_len_tol, close_vertices);

              for(int cur_surf : close_surfs)
              {
                surface_to_node_index.insert( std::multimap<int, int>::value_type( cur_surf, local_node_index));
              }

              for(int cur_curve : close_curves )
              {
                curve_to_node_index.insert( std::multimap<int, int>::value_type( cur_curve, local_node_index));
              }

              for(int cur_vert : close_vertices )
              {
                vertex_to_node_index.insert( std::multimap<int, int>::value_type( cur_vert, local_node_index));
              }

            }
        }
    }

  get_vertices_close_to_nodes(sgi, pgeom, vertex_to_node_index, node_map, vertex_map);
  get_curves_close_to_nodes(sgi, pgeom, curve_to_node_index, node_map, curve_map);
  get_surfaces_close_to_nodes(sgi, pgeom, surface_to_node_index, node_map, surf_map);
}

// TODO - shouldn't this function move to the StructuredBlock class?
CubitVector PGeomAssocStructured::node_coordinates(std::shared_ptr<percept::StructuredBlock> sgi,
                                                   int node_index)
{
    std::array<unsigned,3> indx;
    sgi->multi_dim_indices_from_local_offset(node_index, indx);
    CubitVector node_pt(&sgi->m_sgrid_coords(indx[0], indx[1], indx[2], 0));
    return node_pt;
}

void PGeomAssocStructured::gather_node_coordinates(std::shared_ptr<percept::StructuredBlock> sgi,
                                                   const std::multimap<int,int>::const_iterator lower,
                                                   const std::multimap<int,int>::const_iterator upper,
                                                   std::vector<CubitVector> &node_positions)
{
    std::multimap< int, int >::const_iterator iter;

    for(iter = lower; iter!=upper; ++iter)
    {
      int node_index = iter->second;
      CubitVector node_point = node_coordinates(sgi, node_index);
      node_positions.push_back( node_point );
    }
}


void PGeomAssocStructured::get_surfaces_close_to_nodes(std::shared_ptr<percept::StructuredBlock> sgi,
                                                       std::shared_ptr<PGeom> pgeom,
                                                       const std::multimap< int, int > &candidate_surf_map,
                                                       std::map<int, NodeAssociationData> &node_map,
                                                       std::multimap< int, int > &surface_map)
{
  std::multimap< int, int >::const_iterator iter, lower_iter, upper_iter, key_iter, tmp_iter;
  iter = candidate_surf_map.begin();

  while( iter != candidate_surf_map.end() )
  {
      lower_iter = iter;

      int surface_id = lower_iter->first;
      upper_iter = candidate_surf_map.upper_bound( surface_id );
      iter = upper_iter;

      std::vector<CubitVector> node_positions;
      gather_node_coordinates(sgi, lower_iter, upper_iter, node_positions);

      std::vector<CubitVector> closest_positions;
      pgeom->get_surf_closest_points_trimmed( surface_id, node_positions, closest_positions );

    //put the surface and distance into the node map if within tolerance
    key_iter = lower_iter;
    for(size_t i=0; i<closest_positions.size(); ++i )
    {
        tmp_iter = key_iter;
        key_iter++;

        int node_index = tmp_iter->second;
        NodeAssociationData &nd = node_map[ node_index ];
        double dist = closest_positions[i].distance_between( node_positions[i] );
        if (dist < nd.edge_len_tol)
        {
            nd.add_surface(surface_id, dist);
            surface_map.insert(std::multimap<int,int>::value_type(surface_id, node_index));
        }
    }
  }
}


// TODO - this is too much copy/paste fromt the function for surfaces
// refactor!
void PGeomAssocStructured::get_curves_close_to_nodes(std::shared_ptr<percept::StructuredBlock> sgi,
                                                         std::shared_ptr<PGeom> pgeom,
                                                         const std::multimap< int, int > &candidate_curve_map,
                                                         std::map<int, NodeAssociationData> &node_map,
                                                         std::multimap< int, int > &curve_map)
{
  std::multimap< int, int >::const_iterator iter, lower_iter, upper_iter, key_iter, tmp_iter;
  iter = candidate_curve_map.begin();

  while( iter != candidate_curve_map.end() )
  {
      lower_iter = iter;

      int curve_id = iter->first;
      upper_iter = candidate_curve_map.upper_bound( curve_id );
      iter = upper_iter;

      std::vector<CubitVector> node_positions;
      gather_node_coordinates(sgi, lower_iter, upper_iter, node_positions);

      std::vector<CubitVector> closest_positions;
      bool points_projected = pgeom->get_curv_closest_points_trimmed( curve_id, node_positions, closest_positions );

    if (points_projected)
    {
      //put the curve and distance into the node map
      key_iter = lower_iter;
      for(size_t i=0; i<closest_positions.size(); i++)
      {
          tmp_iter = key_iter;
          key_iter++;

          int node_index = tmp_iter->second;
          //PRINT_INFO("node index %d\n", node_index );
          NodeAssociationData &nd = node_map[ node_index ];
          double dist = closest_positions[i].distance_between( node_positions[i] );
          if (dist < nd.edge_len_tol)
          {
              nd.add_curve(curve_id, dist);
              curve_map.insert(std::multimap<int, int>::value_type(curve_id, node_index));
          }
      }
    }
  }
}

void PGeomAssocStructured::get_vertices_close_to_nodes(std::shared_ptr<percept::StructuredBlock> sgi,
                                                           std::shared_ptr<PGeom> pgeom,
                                                           const std::multimap<int, int> &candidate_vertex_map,
                                                           std::map<int, NodeAssociationData> &node_map,
                                                           std::multimap< int, int > &vertex_map)
{
  //get the closest points for each curve
  std::multimap< int, int >::const_iterator iter, lower_iter, upper_iter, key_iter, tmp_iter;
  iter = candidate_vertex_map.begin();

  while( iter != candidate_vertex_map.end() )
  {
      lower_iter = iter;

      int vertex_id = lower_iter->first;
      upper_iter = candidate_vertex_map.upper_bound( vertex_id );
      iter = upper_iter;

      CubitVector vertex_pt = pgeom->get_vert_coord(vertex_id);

      for(key_iter = lower_iter; key_iter!=upper_iter; )
      {
          tmp_iter = key_iter;
          key_iter++;

          int node_index = tmp_iter->second;
          CubitVector node_pt = node_coordinates(sgi, node_index);

          double dist = vertex_pt.distance_between(node_pt);

          NodeAssociationData &nd = node_map[ node_index ];
          if (dist < nd.edge_len_tol)
          {
              nd.add_vertex(vertex_id, dist);
              vertex_map.insert(std::multimap<int, int>::value_type(vertex_id, node_index));
          }
      }
  }
}


void PGeomAssocStructured::structured_node_edge_lengths(std::shared_ptr<percept::StructuredBlock> sgi,
                                               std::array<unsigned,3> indx,
                                               std::vector<double> &edge_lengths)
{
    //TODO - this handles a single block, but we also need to handle edges that go into an adjacent block

  CubitVector node_pt(&sgi->m_sgrid_coords(indx[0], indx[1], indx[2], 0));

  // minus direction lengths
  for (int i=0; i<3; i++)
  {
      if (indx[i] > sgi->m_sizes.node_min[i])
      {
          std::array<unsigned,3> minus_indx = indx;
          minus_indx[i]--;
          CubitVector minus_pt(&sgi->m_sgrid_coords(minus_indx[0], minus_indx[1], minus_indx[2], 0));
          edge_lengths.push_back(node_pt.distance_between(minus_pt));
      }
  }

  // plus direction lengths
  for (int i=0; i<3; i++)
  {
      if (indx[i] < sgi->m_sizes.node_max[i])
      {
          std::array<unsigned,3> plus_indx = indx;
          plus_indx[i]++;
          CubitVector minus_pt(&sgi->m_sgrid_coords(plus_indx[0], plus_indx[1], plus_indx[2], 0));
          edge_lengths.push_back(node_pt.distance_between(minus_pt));
      }
  }
}


double PGeomAssocStructured::structured_node_local_tolerance(std::shared_ptr<percept::StructuredBlock> sgi,
                                                  std::array<unsigned,3> indx,
                                                  double default_tol)
{
    std::vector<double> edge_lengths;
    structured_node_edge_lengths(sgi, indx, edge_lengths);

    std::sort( edge_lengths.begin(), edge_lengths.end() );

    double shortest_edge_length = edge_lengths[0];

    double local_tol = 0;
    double fraction_of_edge_length = 0.25;

    if( shortest_edge_length > 0.0)
      local_tol = shortest_edge_length*fraction_of_edge_length;

    if( default_tol > local_tol )
      local_tol = default_tol;

    return local_tol;
}


void PGeomAssocStructured::find_edge_association(std::shared_ptr<PGeom> pgeom,
                                                 std::shared_ptr<percept::BlockStructuredGrid> grid)
{
  for (unsigned iblock=0; iblock < grid->m_sblocks.size(); ++iblock)
  {
    std::shared_ptr<percept::StructuredBlock> sgi = grid->m_sblocks[iblock];

    unsigned i_min = sgi->m_sizes.cell_min[0];
    unsigned j_min = sgi->m_sizes.cell_min[1];
    unsigned k_min = sgi->m_sizes.cell_min[2];
    unsigned i_max = sgi->m_sizes.cell_max[0];
    unsigned j_max = sgi->m_sizes.cell_max[1];
    unsigned k_max = sgi->m_sizes.cell_max[2];

    // I direction edges
    // J_MIN face
    std::array<unsigned, 3> range_min = {{i_min, j_min, k_min}};
    std::array<unsigned, 3> range_max = {{i_max, j_min, k_max+1}};
    create_edge_association_data(sgi, range_min, range_max, I_DIRECTION, iblock);
    // J_MAX face
    range_min[1] = j_max+1;
    range_max[1] = j_max+1;
    create_edge_association_data(sgi, range_min, range_max, I_DIRECTION, iblock);
    // K_MIN face
    range_min = {{i_min, j_min+1, k_min}};
    range_max = {{i_max, j_max, k_min}};
    create_edge_association_data(sgi, range_min, range_max, I_DIRECTION, iblock);
    // K_MAX face
    range_min[2] = k_max+1;
    range_max[2] = k_max+1;
    create_edge_association_data(sgi, range_min, range_max, I_DIRECTION, iblock);

    // J direction edges
    // I_MIN face
    range_min = {{i_min, j_min, k_min}};
    range_max = {{i_min, j_max, k_max+1}};
    create_edge_association_data(sgi, range_min, range_max, J_DIRECTION, iblock);
    // I_MAX face
    range_min[0] = i_max+1;
    range_max[0] = i_max+1;
    create_edge_association_data(sgi, range_min, range_max, J_DIRECTION, iblock);
    // K_MIN face
    range_min = {{i_min+1, j_min, k_min}};
    range_max = {{i_max, j_max, k_min}};
    create_edge_association_data(sgi, range_min, range_max, J_DIRECTION, iblock);
    // K_MAX face
    range_min[2] = k_max+1;
    range_max[2] = k_max+1;
    create_edge_association_data(sgi, range_min, range_max, J_DIRECTION, iblock);

    // K direction edges
    // I_MIN face
    range_min = {{i_min, j_min, k_min}};
    range_max = {{i_min, j_max+1, k_max}};
    create_edge_association_data(sgi, range_min, range_max, K_DIRECTION, iblock);
    // I_MAX face
    range_min[0] = i_max+1;
    range_max[0] = i_max+1;
    create_edge_association_data(sgi, range_min, range_max, K_DIRECTION, iblock);
    // J_MIN face
    range_min = {{i_min+1, j_min, k_min}};
    range_max = {{i_max, j_min, k_max}};
    create_edge_association_data(sgi, range_min, range_max, K_DIRECTION, iblock);
    // J_MAX face
    range_min[1] = j_max+1;
    range_max[1] = j_max+1;
    create_edge_association_data(sgi, range_min, range_max, K_DIRECTION, iblock);

  }
}

void PGeomAssocStructured::create_edge_association_data(std::shared_ptr<percept::StructuredBlock> sgi,
                                                        std::array<unsigned, 3> range_min,
                                                        std::array<unsigned, 3> range_max,
                                                        enum IJK_DIRECTION edge_direction,
                                                        unsigned block_num)
{
    std::map<int, NodeAssociationData> &node_map = mBlockAssociation[block_num].node_to_geom_map;

    std::map<int, EdgeAssociationData> &edge_map = mBlockAssociation[block_num].edge_to_geom_map[edge_direction];
    std::multimap<int, int> &curve_map = mBlockAssociation[block_num].curve_to_edge_map[edge_direction];
    std::multimap<int, int> &surf_map = mBlockAssociation[block_num].surface_to_edge_map[edge_direction];

    // create an array with the correct i,j,k increment to get the second node on the edge
    std::array<unsigned, 3> edge_delta{{0, 0, 0}};
    edge_delta[edge_direction] = 1;

    //TODO - do we need to get this data if we make have local_offset use the data internally
    //     - can we get rid of the access ordering and always access in the same way
    const unsigned A0 = sgi->m_access_ordering[0], A1 = sgi->m_access_ordering[1], A2 = sgi->m_access_ordering[2];
    unsigned Asizes[3] = {sgi->m_sizes.node_size[A0], sgi->m_sizes.node_size[A1], sgi->m_sizes.node_size[A2]};

    std::array<unsigned,3> indx{{0,0,0}};
    for (indx[2] = range_min[2]; indx[2] <= range_max[2]; ++indx[2]) // K
    {
        for (indx[1] = range_min[1]; indx[1] <= range_max[1]; ++indx[1]) // J
        {
            for (indx[0] = range_min[0]; indx[0] <= range_max[0]; ++indx[0]) // I
            {
                std::vector<int> edge_nodes(2);  // assumes linear elements

                //TODO - local_offset should be a member function and use the block's internal data rather than passing it in
                edge_nodes[0] = sgi->local_offset(indx[0], indx[1], indx[2], Asizes);
                edge_nodes[1] = sgi->local_offset(indx[0]+edge_delta[0],
                                                  indx[1]+edge_delta[1],
                                                  indx[2]+edge_delta[2],
                                                  Asizes);

                std::vector<int> common_curves;
                find_common_curves(edge_nodes, node_map, common_curves);

                std::vector<int> common_surfs;
                find_common_surfaces(edge_nodes, node_map, common_surfs);

                if (common_curves.size() || common_surfs.size())
                {
                    int edge_id = edge_nodes[0];

                    EdgeAssociationData ed;

                    if (common_curves.size())
                    {
                        ed.close_curves = common_curves;

                        for (int i_curve : common_curves)
                            curve_map.insert(std::multimap<int, int>::value_type(i_curve, edge_id));
                    }

                    if (common_surfs.size())
                    {
                        ed.close_surfs = common_surfs;

                        for (int i_surf : common_surfs)
                            surf_map.insert(std::multimap<int, int>::value_type(i_surf, edge_id));
                    }

                    edge_map[edge_id] = ed;
                }
            }
        }
    }
}

// TODO - combine functions find_common_curves and find_common_surfaces
void PGeomAssocStructured::find_common_curves(std::vector<int> const &elem_nodes,
                                              std::map<int, NodeAssociationData> &node_map,
                                              std::vector<int> &common_curves)
{
    int master_node = elem_nodes[0];
    NodeAssociationData master_nd = node_map[master_node];

    for(size_t j=0; j<master_nd.close_curves.size(); j++)
    {
        // For every curve that the master node is on within
        // tolerance see if the other nodes are on the same
        // curve within tolerance.
//TODO - should we sort the node data before this so we can always take the curves in order ?????
        int cur_master_curve = master_nd.close_curves[j].geometry_id;

        // now see if the other nodes have this curve.
        unsigned num_matches = 1;
        for(size_t k=1; k<elem_nodes.size(); k++)
        {
            int cur_node = elem_nodes[k];
            if (find_matching_geometry(cur_master_curve, node_map[cur_node].close_curves))
            {
                num_matches++;
            }
        }
        if(num_matches == elem_nodes.size())
            common_curves.push_back(cur_master_curve);
    }
}


void PGeomAssocStructured::find_common_surfaces(std::vector<int> const &elem_nodes,
                                                std::map<int, NodeAssociationData> &node_map,
                                                std::vector<int> &common_surfs)
{
    int master_node = elem_nodes[0];
    NodeAssociationData master_nd = node_map[master_node];

    for(size_t j=0; j<master_nd.close_surfs.size(); j++)
    {
        // For every surface that the master node is on within
        // tolerance see if the other nodes are on the same
        // surface within tolerance.
        int cur_master_surf = master_nd.close_surfs[j].geometry_id;

        // now see if the other nodes have this surface
        unsigned num_matches = 1;
        for(size_t k=1; k<elem_nodes.size(); k++)
        {
            int cur_node = elem_nodes[k];
            if (find_matching_geometry(cur_master_surf, node_map[cur_node].close_surfs))
            {
                num_matches++;
            }
        }

        if(num_matches == elem_nodes.size())
            common_surfs.push_back(cur_master_surf);
    }
}

bool PGeomAssocStructured::find_matching_geometry(int geom_id_to_match, std::vector<CloseGeometry> &close_geometry)
{
    for(size_t m=0; m<close_geometry.size(); m++ )
    {
        if(close_geometry[m].geometry_id == geom_id_to_match)
        {
            return true;
        }
    }
    return false;
}

void PGeomAssocStructured::find_face_association(std::shared_ptr<PGeom> pgeom,
                                                 std::shared_ptr<percept::BlockStructuredGrid> grid)
{
  for (unsigned iblock=0; iblock < grid->m_sblocks.size(); ++iblock)
  {
    std::shared_ptr<percept::StructuredBlock> sgi = grid->m_sblocks[iblock];

    unsigned i_min = sgi->m_sizes.cell_min[0];
    unsigned j_min = sgi->m_sizes.cell_min[1];
    unsigned k_min = sgi->m_sizes.cell_min[2];
    unsigned i_max = sgi->m_sizes.cell_max[0];
    unsigned j_max = sgi->m_sizes.cell_max[1];
    unsigned k_max = sgi->m_sizes.cell_max[2];

    // I direction faces
    // I_MIN face
    std::array<unsigned, 3> range_min = {{i_min, j_min, k_min}};
    std::array<unsigned, 3> range_max = {{i_min, j_max, k_max}};
    create_face_association_data(sgi, range_min, range_max, I_DIRECTION, iblock);

    // I_MAX face
    range_min[0] = i_max+1;
    range_max[0] = i_max+1;
    create_face_association_data(sgi, range_min, range_max, I_DIRECTION, iblock);

    // J direction faces
    // J_MIN face
    range_min = {{i_min, j_min, k_min}};
    range_max = {{i_max, j_min, k_max}};
    create_face_association_data(sgi, range_min, range_max, J_DIRECTION, iblock);

    // J_MAX face
    range_min[1] = j_max+1;
    range_max[1] = j_max+1;
    create_face_association_data(sgi, range_min, range_max, J_DIRECTION, iblock);

    // K_MIN face
    range_min = {{i_min, j_min, k_min}};
    range_max = {{i_max, j_max, k_min}};
    create_face_association_data(sgi, range_min, range_max, K_DIRECTION, iblock);

    // K_MAX face
    range_min[2] = k_max+1;
    range_max[2] = k_max+1;
    create_face_association_data(sgi, range_min, range_max, K_DIRECTION, iblock);

  }
}

void PGeomAssocStructured::create_face_association_data(std::shared_ptr<percept::StructuredBlock> sgi,
                                                        std::array<unsigned, 3> range_min,
                                                        std::array<unsigned, 3> range_max,
                                                        enum IJK_DIRECTION face_direction,
                                                        unsigned block_num)
{
    std::map<int, NodeAssociationData> &node_map = mBlockAssociation[block_num].node_to_geom_map;
    std::map<int, FaceAssociationData> &face_map = mBlockAssociation[block_num].face_to_geom_map[face_direction];
    std::multimap<int, int> &surf_map = mBlockAssociation[block_num].surface_to_face_map[face_direction];

    std::array<uint64_t,3> indx{{0,0,0}};

    for (indx[2] = range_min[2]; indx[2] <= range_max[2]; ++indx[2]) // K
    {
        for (indx[1] = range_min[1]; indx[1] <= range_max[1]; ++indx[1]) // J
        {
            for (indx[0] = range_min[0]; indx[0] <= range_max[0]; ++indx[0]) // I
            {
                std::vector<int> face_nodes(4);  // assumes linear elements

                get_face_nodes(sgi, indx[0], indx[1], indx[2], face_direction, face_nodes);

                std::vector<int> common_surfs;
                find_common_surfaces(face_nodes, node_map, common_surfs);

                if (common_surfs.size())
                {
                    int face_id = face_nodes[0];

                    FaceAssociationData fd;
                    fd.close_surfs = common_surfs;

                    face_map[face_id] = fd;

                    for (int i_surf : common_surfs)
                        surf_map.insert(std::multimap<int, int>::value_type(i_surf, face_id));
                }

            }
        }
    }
}

void PGeomAssocStructured::get_face_nodes(std::shared_ptr<percept::StructuredBlock> sgi,
                                          uint64_t i,
                                          uint64_t j,
                                          uint64_t k,
                                          enum IJK_DIRECTION face_direction,
                                          std::vector<int> &face_nodes)
{
    //TODO - do we need to get this data if we make have local_offset use the data internally
    //     - can we get rid of the access ordering and always access in the same way
    const unsigned A0 = sgi->m_access_ordering[0], A1 = sgi->m_access_ordering[1], A2 = sgi->m_access_ordering[2];
    unsigned Asizes[3] = {sgi->m_sizes.node_size[A0], sgi->m_sizes.node_size[A1], sgi->m_sizes.node_size[A2]};


//    - gather the four face nodes
//    - set a delta (i,j,k triplet) for each of the two adjacent nodes of the quad face
//      and a delta for the opposite node of the face
//    - the i, j, k direction is used to get the face with constant i, j, or k index respectively.
//    - use the i,j,k computed from the deltas to get the local offset from the structured block

//
//    i,j,k
//    + adj_delta2
//     _______________i,j,k + opp_delta
//    |               |
//    |               |
//    |               |
//    |               |
//    |_______________|
//    i,j,k      i,j,k + adj_delta1
    std::array<unsigned, 3> adj_delta1{{0, 0, 0}};
    std::array<unsigned, 3> adj_delta2{{0, 0, 0}};
    std::array<unsigned, 3> opp_delta {{0, 0, 0}};

    switch(face_direction)
    {
        case I_DIRECTION:
            adj_delta1 = {{0, 1, 0}};
            opp_delta = {{0, 1, 1}};
            adj_delta2 = {{0, 0, 1}};
            break;
        case J_DIRECTION:
            adj_delta1 = {{1, 0, 0}};
            opp_delta = {{1, 0, 1}};
            adj_delta2 = {{0, 0, 1}};
            break;
        case K_DIRECTION:
            adj_delta1 = {{1, 0, 0}};
            opp_delta = {{1, 1, 0}};
            adj_delta2 = {{0, 1, 0}};
            break;
        default:
            return;
    }

    //TODO - local_offset should use the block's internal (Asizes) data rather than passing it in

    face_nodes[0] = sgi->local_offset(i, j, k, Asizes);
    face_nodes[1] = sgi->local_offset(i + adj_delta1[0],
                                      j + adj_delta1[1],
                                      k + adj_delta1[2],
                                      Asizes);
    face_nodes[2] = sgi->local_offset(i + opp_delta[0],
                                      j + opp_delta[1],
                                      k + opp_delta[2],
                                      Asizes);
    face_nodes[3] = sgi->local_offset(i + adj_delta2[0],
                                      j + adj_delta2[1],
                                      k + adj_delta2[2],
                                      Asizes);
    return;
}

void PGeomAssocStructured::make_refined_mesh_association(std::shared_ptr<percept::BlockStructuredGrid> input_grid,
                                                         std::shared_ptr<percept::BlockStructuredGrid> output_grid)
{
    mRefinedNodeToGeomMaps.resize(input_grid->m_sblocks.size());

    make_refined_edge_mesh_association(input_grid, output_grid, I_DIRECTION);
    make_refined_edge_mesh_association(input_grid, output_grid, J_DIRECTION);
    make_refined_edge_mesh_association(input_grid, output_grid, K_DIRECTION);

    make_refined_face_mesh_association(input_grid, output_grid, I_DIRECTION);
    make_refined_face_mesh_association(input_grid, output_grid, J_DIRECTION);
    make_refined_face_mesh_association(input_grid, output_grid, K_DIRECTION);
}

void PGeomAssocStructured::make_refined_edge_mesh_association(std::shared_ptr<percept::BlockStructuredGrid> input_grid,
                                                              std::shared_ptr<percept::BlockStructuredGrid> output_grid,
                                                              enum IJK_DIRECTION edge_direction)
{
//
// for now assumes mesh is being doubled
//

    for (unsigned iblock=0; iblock < input_grid->m_sblocks.size(); ++iblock)
    {
        std::shared_ptr<percept::StructuredBlock> sg_in = input_grid->m_sblocks[iblock];
        std::shared_ptr<percept::StructuredBlock> sg_out = output_grid->m_sblocks[iblock];

        std::map<int, EdgeAssociationData> &block_edge_map = mBlockAssociation[iblock].edge_to_geom_map[edge_direction];
        std::map<uint64_t, NodeAssociatedGeom> &refine_node_map = mRefinedNodeToGeomMaps[iblock];


        std::map< int, EdgeAssociationData >::const_iterator edge_data_iter;

        for(edge_data_iter = block_edge_map.begin();
            edge_data_iter != block_edge_map.end();
            edge_data_iter++)
        {
          int edge_id = edge_data_iter->first;

          uint64_t refined_edge_node = refined_edge_node_id(edge_id, edge_direction, sg_in, sg_out);

          NodeAssociatedGeom node_geom;

          EdgeAssociationData const &ed = edge_data_iter->second;
          if (ed.close_curves.size())
          {
            node_geom.dimension = 1;
            node_geom.geom_ids = ed.close_curves;
            refine_node_map[refined_edge_node] = node_geom;
          }
          else if (ed.close_surfs.size())
          {
            node_geom.dimension = 2;
            node_geom.geom_ids = ed.close_surfs;
            refine_node_map[refined_edge_node] = node_geom;
          }
          else
          {
              std::ostringstream oss;
              oss << "ERROR: edge " << edge_id << " " << direction_string[edge_direction] << " is not associated with geometry";
              throw std::runtime_error(oss.str());
          }
        }
    }
}

void PGeomAssocStructured::make_refined_face_mesh_association(std::shared_ptr<percept::BlockStructuredGrid> input_grid,
                                                              std::shared_ptr<percept::BlockStructuredGrid> output_grid,
                                                              enum IJK_DIRECTION face_direction)
{
//
// for now assumes mesh is being doubled
//

    for (unsigned iblock=0; iblock < input_grid->m_sblocks.size(); ++iblock)
    {
        std::shared_ptr<percept::StructuredBlock> sg_in = input_grid->m_sblocks[iblock];
        std::shared_ptr<percept::StructuredBlock> sg_out = output_grid->m_sblocks[iblock];

        std::map<int, FaceAssociationData> &block_face_map = mBlockAssociation[iblock].face_to_geom_map[face_direction];

        std::map<uint64_t, NodeAssociatedGeom> &refine_node_map = mRefinedNodeToGeomMaps[iblock];

        std::map< int, FaceAssociationData >::const_iterator face_data_iter;

        for(face_data_iter = block_face_map.begin();
            face_data_iter != block_face_map.end();
            face_data_iter++)
        {
          int face_id = face_data_iter->first;

          uint64_t refined_face_node = refined_face_node_id(face_id, face_direction, sg_in, sg_out);

          NodeAssociatedGeom node_geom;

          FaceAssociationData const &fd = face_data_iter->second;
          if (fd.close_surfs.size())
          {
            node_geom.dimension = 2;
            node_geom.geom_ids = fd.close_surfs;
            refine_node_map[refined_face_node] = node_geom;
          }
          else
          {
              std::ostringstream oss;
              oss << "ERROR: face " << face_id << " " << direction_string[face_direction] << " is not associated with geometry";
              throw std::runtime_error(oss.str());
          }
        }

    }

}


uint64_t PGeomAssocStructured::refined_edge_node_id(int edge_id,
                                                    enum IJK_DIRECTION direction,
                                                    std::shared_ptr<percept::StructuredBlock> sg_in,
                                                    std::shared_ptr<percept::StructuredBlock> sg_out)
{
    //TODO - maybe this should be in the StructuredBlock class???

    std::array<unsigned,3> indx;
    sg_in->multi_dim_indices_from_local_offset(edge_id, indx);

    // mesh is doubled in size
    for (unsigned &ijk : indx)
    {
        ijk *= 2;
    }

    indx[direction] += 1;

    unsigned sizes[3] = {sg_out->m_sizes.node_size[0],
                         sg_out->m_sizes.node_size[1],
                         sg_out->m_sizes.node_size[2]};
    uint64_t refined_id = sg_out->local_offset(indx[0], indx[1], indx[2], sizes);

    return refined_id;
}

uint64_t PGeomAssocStructured::refined_face_node_id(int start_node,
                                                    enum IJK_DIRECTION direction,
                                                    std::shared_ptr<percept::StructuredBlock> sg_in,
                                                    std::shared_ptr<percept::StructuredBlock> sg_out)
{
    //TODO - maybe this should be in the StructuredBlock class???

    std::array<unsigned,3> indx;
    sg_in->multi_dim_indices_from_local_offset(start_node, indx);

    // mesh is doubled in size
    for (unsigned &ijk : indx)
    {
        ijk *= 2;
    }

    switch (direction)
    {
        case I_DIRECTION: // constant i face
            indx[1] += 1;
            indx[2] += 1;
            break;

        case J_DIRECTION: // constant j face
            indx[0] += 1;
            indx[2] += 1;
            break;

        case K_DIRECTION: // constant k face
            indx[0] += 1;
            indx[1] += 1;
            break;
    }

    unsigned sizes[3] = {sg_out->m_sizes.node_size[0],
                          sg_out->m_sizes.node_size[1],
                          sg_out->m_sizes.node_size[2]};
    uint64_t refined_id = sg_out->local_offset(indx[0], indx[1], indx[2], sizes);

    return refined_id;
}

std::map<uint64_t, NodeAssociatedGeom> & PGeomAssocStructured::get_refined_node_map(unsigned iblock)
{
    return mRefinedNodeToGeomMaps[iblock];
}



void PGeomAssocStructured::print_block_geometry_association() const
{
    for (unsigned iblock=0; iblock < mBlockAssociation.size(); ++iblock)
    {
        std::cout << "Block " << iblock << " Association Data" << std::endl;
        const BlockAssociationData &block_assoc = mBlockAssociation[iblock];

        print_node_association(block_assoc.node_to_geom_map);

        for (int i=0; i<3; i++)
        {
            print_edge_association_map(block_assoc.edge_to_geom_map[i], direction_string[i]);
        }

        for (int i=0; i<3; i++)
        {
            print_face_association_map(block_assoc.face_to_geom_map[i], direction_string[i]);
        }

        print_geometry_association("Vertex", " node", "", block_assoc.vertex_to_node_map);

        print_geometry_association("Curve", " node", "", block_assoc.curve_to_node_map);
        for (int i=0; i<3; i++)
        {
            print_geometry_association("Curve", " edge", direction_string[i], block_assoc.curve_to_edge_map[i]);
        }

        print_geometry_association("Surface", " node", "", block_assoc.surface_to_node_map);
        for (int i=0; i<3; i++)
        {
            print_geometry_association("Surface", " edge", direction_string[i], block_assoc.surface_to_edge_map[i]);
        }
        for (int i=0; i<3; i++)
        {
            print_geometry_association("Surface", " face", direction_string[i], block_assoc.surface_to_face_map[i]);
        }
    }
}

void PGeomAssocStructured::print_node_association(const std::map<int, NodeAssociationData> &node_to_geom_map) const
{
    std::map< int, NodeAssociationData >::const_iterator node_data_iter;

    for(node_data_iter = node_to_geom_map.begin();
        node_data_iter != node_to_geom_map.end();
        node_data_iter++)
    {
      int node_id = node_data_iter->first;
      NodeAssociationData const &nd = node_data_iter->second;

      std::cout << "Node " << node_id << std::endl;

      print_close_geometry_ids(" vertices", nd.close_verts);
      print_close_geometry_ids(" curves", nd.close_curves);
      print_close_geometry_ids(" surfaces", nd.close_surfs);
    }
}

void PGeomAssocStructured::print_edge_association_map(const std::map<int, EdgeAssociationData> &edge_map,
                                                      const std::string &direction_string) const
{
    std::map< int, EdgeAssociationData >::const_iterator edge_data_iter;

    for(edge_data_iter = edge_map.begin();
        edge_data_iter != edge_map.end();
        edge_data_iter++)
    {
      int edge_id = edge_data_iter->first;
      EdgeAssociationData const &ed = edge_data_iter->second;

      std::cout << "Edge: " << edge_id << " " << direction_string << std::endl;

      print_ids(" curves", ed.close_curves);
      print_ids(" surfaces", ed.close_surfs);
    }
}

void PGeomAssocStructured::print_face_association_map(const std::map<int, FaceAssociationData> &face_map,
                                                      const std::string &direction_string) const
{
  std::map< int, FaceAssociationData >::const_iterator face_data_iter;

    for(face_data_iter = face_map.begin();
        face_data_iter != face_map.end();
        face_data_iter++)
    {
      int face_id = face_data_iter->first;
      FaceAssociationData const &fd = face_data_iter->second;

      std::cout << "Face: " << face_id << " " << direction_string << std::endl;
      print_ids(" surfaces", fd.close_surfs);

    }
}

void PGeomAssocStructured::print_geometry_association(const std::string &key_name,
                                                      const std::string &value_name,
                                                      const std::string &direction_string,
                                                      const std::multimap< int, int > &geometry_map)
{
    std::multimap< int, int >::const_iterator iter, upper_iter;

    iter = geometry_map.begin();
    while (iter != geometry_map.end())
    {
        int geom_id = iter->first;
        std::cout << key_name << " " << geom_id << " " << direction_string << std::endl;

        upper_iter = geometry_map.upper_bound(geom_id);
        print_multimap_ids(value_name, iter, upper_iter);
        iter = upper_iter;
    }
}

void PGeomAssocStructured::print_ids(const std::string &label,
                                     const std::vector<int> &id_list)
{
    std::cout << label;

    for (int i : id_list)
    {
        std::cout << " " << i;
    }
    std::cout << std::endl;
}

void PGeomAssocStructured::print_multimap_ids(const std::string &label,
                                              const std::multimap<int,int>::const_iterator lower,
                                              const std::multimap<int,int>::const_iterator upper)
{
    std::multimap< int, int >::const_iterator iter;

    std::cout << label;
    for (iter=lower ; iter!=upper; iter++)
    {
        std::cout << " " << iter->second;
    }
    std::cout << std::endl;
}

void PGeomAssocStructured::print_close_geometry_ids(const std::string &label,
                                                    const std::vector<struct CloseGeometry> &close_geom)
{
    std::cout << label;
    for (struct CloseGeometry close : close_geom)
    {
        std::cout << " " << close.geometry_id;// << " " << close.distance << " ";
    }
    std::cout << std::endl;
}

void PGeomAssocStructured::check_cylinder_surface_projection(std::shared_ptr<PGeom> pgeom,
                                                             std::shared_ptr<percept::BlockStructuredGrid> grid,
                                                             int surface_id,
                                                             double radius,
                                                             double tolerance) const
{
    std::cout << "Testing node points on cylinder (surface " << surface_id << ") for correct xy plane radius " << radius << std::endl;
    std::map<uint64_t, NodeAssociatedGeom>::const_iterator iter;

    for (unsigned iblock=0; iblock < mBlockAssociation.size(); ++iblock)
    {
        std::shared_ptr<percept::StructuredBlock> sgi = grid->m_sblocks[iblock];
        const std::map<uint64_t, NodeAssociatedGeom> &node_to_geom_map = mRefinedNodeToGeomMaps[iblock];

        for (iter=node_to_geom_map.begin(); iter!=node_to_geom_map.end(); iter++)
        {
            const NodeAssociatedGeom &geom = iter->second;

            if (2 == geom.dimension)
            {
                if (std::find(geom.geom_ids.begin(), geom.geom_ids.end(), surface_id) != geom.geom_ids.end())
                {
                    int node_id = iter->first;
                    CubitVector node_point = node_coordinates(sgi, node_id);
                    double xy_radius = sqrt(node_point.x()*node_point.x() +
                                            node_point.y()*node_point.y());
                    if (fabs(radius - xy_radius) > tolerance)
                    {
                        std::ostringstream oss;
                        oss << "ERROR: projected point does not lie on cylinder.  Desired radius=" << radius << " Actual radius=" << xy_radius;
                        throw std::runtime_error(oss.str());
                    }
                }
            }
        }
    }
}

#endif
