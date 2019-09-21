// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "TopologyVerifier.hpp"
#include <iostream>
#include "PerceptMesh.hpp"

  namespace percept
  {
    using namespace interface_table;

    TopologyVerifier::TopologyVerifier()
    {
      build_invalid_edge_sets();
    }

    void TopologyVerifier::build_invalid_edge_sets()
    {
      // build tables for edge checking
      const CellTopologyData *topos[NELT];

      m_invalid_edge_set = std::vector<invalid_edge_set_type>(NELT);

      for (unsigned i = 0; i < NELT; i++)
        {
          topos[i] = ShardsInterfaceTable::s_elemInfo[i].cellTopoData;
        }
      //topos[interface_table::shards_Wedge_6] = shards::getCellTopologyData<shards::Wedge<6> >();
      // etc. FIXME

      for (unsigned ielt = 0; ielt < NELT; ielt++)
        {
          if (topos[ielt])
            {
              // generate all possible edges in the element, including diagonals, etc.
              unsigned vertex_count = topos[ielt]->vertex_count;
              if (vertex_count > 1)
                {
                  for (unsigned iv = 0; iv < vertex_count-1; iv++)
                    {
                      for (unsigned jv = iv+1; jv < vertex_count; jv++)
                        {
                          MyEdge<unsigned> edge(iv, jv);
                          m_invalid_edge_set[ielt].insert(edge);
                        }
                    }

                  // remove the actual valid edges so only invalid remain
                  unsigned edge_count = topos[ielt]->edge_count;
                  for (unsigned ie = 0; ie < edge_count; ie++)
                    {
                      MyEdge<unsigned> edge(topos[ielt]->edge[ie].node[0], topos[ielt]->edge[ie].node[1]);
                      m_invalid_edge_set[ielt].erase(edge);
                    }
                }
            }
        }

    }

    /// return true if topology is bad
    bool TopologyVerifier::isTopologyBad(stk::mesh::BulkData& bulkData, stk::mesh::Entity elem)
    {
      const CellTopologyData * const top = stk::mesh::get_cell_topology(bulkData.bucket(elem).topology()).getCellTopologyData();

      const MyPairIterRelation elem_nodes(bulkData, elem, stk::topology::NODE_RANK );

#if 0
      std::cout << "top->node_count = " << top->node_count << "\n";
      std::cout << "elem_nodes.size() = " << elem_nodes.size() << "\n";
      std::cout.flush();
#endif

      for (unsigned j = 0; j < elem_nodes.size()-1; j++)
        {
          stk::mesh::Entity node = elem_nodes[ j ].entity();

          for ( unsigned i = j+1 ; i < top->node_count ; ++i ) {
            {
              stk::mesh::Entity nodei = elem_nodes[ i ].entity();
              if (bulkData.identifier(node) == bulkData.identifier(nodei))
                {
                  return true;
                }
            }
          }
        }
      return false;
    }



    /**
     * Algorithm:
     *   1. for each element, loop over its edges, for each edge's node, loop over elements attached to node
     *   1a.   for each neighboring element, check if current edge is invalid
     */
    bool TopologyVerifier::isTopologyBad(stk::mesh::BulkData& bulk) //, stk::mesh::Part& mesh_part )
    {
      const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

      stk::mesh::Field<double, stk::mesh::Cartesian> *coord_field =
        meta.get_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");

      //mesh::Selector select_owned( meta_data.locally_owned_part() );

      const stk::mesh::BucketVector & buckets = bulk.buckets(stk::topology::ELEMENT_RANK );

      for ( stk::mesh::BucketVector::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
        {
          // if ( select_owned( **ik ) ) {...

          const stk::mesh::Bucket & bucket = **ik ;

          // Number of elems in this bucket of elems and elem field data
          const unsigned number_elems = bucket.size();

          double * elem_node_data = stk::mesh::field_data( *coord_field , bucket, 0 );
          //double * elem_centroid_data = field_data( elem_centroid_field , bucket.begin() );

          // FIXME
          if (0) { elem_node_data[0]++;}

#if 1
          //const CellTopologyData * const bucket_cell_topo = m_eMesh.get_cell_topology(bucket);
          const CellTopologyData * const bucket_cell_topo = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();
          int bucket_shardsId = ShardsInterfaceTable::s_singleton.lookupShardsId(bucket_cell_topo->name);
#endif

          //if (0) { std::cout << bucket_cell_topo->name; }
          if (0) { std::cout << "bucket_shardsId= " << bucket_shardsId << " name= " << bucket_cell_topo->name <<  std::endl; }

          if (0) { std::cout << "number_elems= " << number_elems << std::endl;}

          for ( unsigned i = 0 ; i < number_elems ; ++i)
            {
              stk::mesh::Entity elem = bucket[i] ;
              bool isDuplicateNode = isTopologyBad(bulk, elem);
              if (isDuplicateNode)
                {
                  std::cout << "duplicate node found: elem = " << elem << std::endl;
                  return true;
                }
              if (0) std::cout << "elemOfBucket= " << elem << std::endl;
              const MyPairIterRelation elem_nodes(bulk, elem, stk::topology::NODE_RANK );

              const CellTopologyData * const cell_topo = stk::mesh::get_cell_topology(bulk.bucket(elem).topology()).getCellTopologyData();

              int shardsId = ShardsInterfaceTable::s_singleton.lookupShardsId(cell_topo->name);
              if (0) { std::cout << "shardsId= " << shardsId << " name= " << cell_topo->name <<  std::endl; }

              for (unsigned iedgeOrd = 0; iedgeOrd < cell_topo->edge_count; iedgeOrd++)
                {
                  //const CellTopologyData_Subcell& edge =

                  unsigned in0 = cell_topo->edge[iedgeOrd].node[0];
                  unsigned in1 = cell_topo->edge[iedgeOrd].node[1];

                  MyEdge<unsigned> potential_bad_edge(bulk.identifier(elem_nodes[in0].entity()), bulk.identifier(elem_nodes[in1].entity()));
                  //if (potential_bad_edge.getId0() == 3 && potential_bad_edge.getId1() == 6)
                  if (0) std::cout << "potential_bad_edge: " << potential_bad_edge.getId0() << " " << potential_bad_edge.getId1() << std::endl;
                }
              for (unsigned iedgeOrd = 0; iedgeOrd < cell_topo->edge_count; iedgeOrd++)
                {
                  //const CellTopologyData_Subcell& edge =

                  unsigned in0 = cell_topo->edge[iedgeOrd].node[0];
                  unsigned in1 = cell_topo->edge[iedgeOrd].node[1];

                  MyEdge<unsigned> potential_bad_edge(bulk.identifier(elem_nodes[in0].entity()), bulk.identifier(elem_nodes[in1].entity() ));
                  //if (potential_bad_edge.getId0() == 3 && potential_bad_edge.getId1() == 6)
                  if (0) std::cout << "potential_bad_edge: " << potential_bad_edge.getId0() << " " << potential_bad_edge.getId1() << std::endl;

                  if (0 && MyEdge<unsigned>(3,6) == potential_bad_edge)
                    {
                      std::cout << "bad edge" << std::endl;
                    }

                  for (unsigned inodeOnPotBadEdge = 0; inodeOnPotBadEdge < 2; inodeOnPotBadEdge++)
                    {
                      unsigned inodeOnPotBadEdgeInElem = cell_topo->edge[iedgeOrd].node[inodeOnPotBadEdge];

                      const MyPairIterRelation node_elems(bulk, elem_nodes[inodeOnPotBadEdgeInElem].entity(), stk::topology::ELEMENT_RANK );
                      unsigned num_elems_on_node = node_elems.size();

                      for (unsigned iele = 0; iele < num_elems_on_node; iele++)
                        {
                          stk::mesh::Entity elemOnNode = node_elems[iele].entity();
                          const MyPairIterRelation elemOnNode_nodes(bulk, elemOnNode, stk::topology::NODE_RANK );

                          //const CellTopologyData * const local_cell_topo = m_eMesh.get_cell_topology(elemOnNode);
                          const CellTopologyData * const local_cell_topo = stk::mesh::get_cell_topology(bulk.bucket(elemOnNode).topology()).getCellTopologyData();
                          int local_shardsId = ShardsInterfaceTable::s_singleton.lookupShardsId(local_cell_topo->name);
                          //if (1) { std::cout << "shardsId= " << shardsId << " name= " << cell_topo->name <<  std::endl; }

                          if (0) std::cout << "elemOnNode= " << elemOnNode << std::endl;
                          unsigned num_invalid_edges = m_invalid_edge_set[local_shardsId].size();
                          if (0) { std::cout << num_invalid_edges; }
                          for (invalid_edge_set_type::iterator inv_edge = m_invalid_edge_set[local_shardsId].begin();
                               inv_edge != m_invalid_edge_set[local_shardsId].end(); inv_edge++)
                            {
                              MyEdge<unsigned> globalIdInvEdge( bulk.identifier(elemOnNode_nodes[ (*inv_edge).getId0()].entity()),
                                                                bulk.identifier(elemOnNode_nodes[ (*inv_edge).getId1()].entity()) );

                              if (0) std::cout << "globalIdInvEdge: " << globalIdInvEdge.getId0() << " " << globalIdInvEdge.getId1() << std::endl;
                              if(potential_bad_edge == globalIdInvEdge)
                                {
                                  return true;
                                }
                            }
                        }
                    }
                }


            }
        }
      return false;
    }


  }//namespace percept
