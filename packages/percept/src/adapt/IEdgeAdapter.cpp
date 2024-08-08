// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <adapt/IEdgeAdapter.hpp>

  namespace percept {

    void IEdgeAdapter::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
                                            vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      int spatialDimension = m_eMesh.get_spatial_dim();

      CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

      stk::mesh::FieldBase* coordField = m_eMesh.get_coordinates_field();

      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;
          stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              numSubDimNeededEntities = cell_topo_data->edge_count;
            }
          else if (needed_entity_rank == m_eMesh.face_rank())
            {
              numSubDimNeededEntities = cell_topo_data->side_count;
              throw std::runtime_error("IEdgeAdapter::apply can't use IEdgeAdapter for RefinerPatterns that require face nodes");
            }
          else if (needed_entity_rank == stk::topology::ELEMENT_RANK)
            {
              numSubDimNeededEntities = 1;
              throw std::runtime_error("IEdgeAdapter::apply can't use IEdgeAdapter for RefinerPatterns that require volume nodes");
            }

          // see how many edges are already marked
          std::vector<int> edge_marks(numSubDimNeededEntities,0);
          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
                  if (!is_empty)
                    {
                      edge_marks[iSubDimOrd] = 1;
                    }
                }
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              if (needed_entity_rank == m_eMesh.edge_rank())
                {
                  stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 = static_cast<double*>(stk::mesh::field_data( *coordField , node0 ));
                  double * const coord1 = static_cast<double*>(stk::mesh::field_data( *coordField , node1 ));
                  double  dcoord0[3] = {coord0[0],coord0[1], (spatialDimension==2?0:coord0[2])};
                  double  dcoord1[3] = {coord1[0],coord1[1], (spatialDimension==2?0:coord1[2])};

                  int markInfo = markEdge(element, iSubDimOrd, node0, node1, dcoord0, dcoord1, &edge_marks);

                  bool needNodes = (DO_REFINE & markInfo);
                    {
                      (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, needNodes,bucket_topo_data);
                    }
                }

            } // iSubDimOrd
        } // ineed_ent
    }

    int IEdgeAdapter::markCountRefinedEdges(const stk::mesh::Entity element)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
      int spatialDimension = m_eMesh.get_spatial_dim();
      CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

      stk::mesh::FieldBase* coordField = m_eMesh.get_coordinates_field();

      unsigned numSubDimNeededEntities = 0;
      numSubDimNeededEntities = cell_topo_data->edge_count;

      int ref_count=0;
      for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
        {
          stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
          stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
          double * const coord0 = static_cast<double*>(stk::mesh::field_data( *coordField , node0 ));
          double * const coord1 = static_cast<double*>(stk::mesh::field_data( *coordField , node1 ));
          double  dcoord0[3] = {coord0[0],coord0[1], (spatialDimension==2?0:coord0[2])};
          double  dcoord1[3] = {coord1[0],coord1[1], (spatialDimension==2?0:coord1[2])};

          int markInfo = markEdge(element, iSubDimOrd, node0, node1, dcoord0, dcoord1, 0);
          bool do_ref = markInfo & DO_REFINE;
          if (do_ref)
            {
              ++ref_count;
            }
        }
      return ref_count;
    }

    int IEdgeAdapter::markUnrefine(const stk::mesh::Entity element)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

      int spatialDimension = m_eMesh.get_spatial_dim();
      stk::mesh::FieldBase* coordField = m_eMesh.get_coordinates_field();

      unsigned numSubDimNeededEntities = 0;
      numSubDimNeededEntities = cell_topo_data->edge_count;

      bool unrefAllEdges = true;
      for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
        {
          stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
          stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
          double * const coord0 = static_cast<double*>(stk::mesh::field_data( *coordField , node0 ));
          double * const coord1 = static_cast<double*>(stk::mesh::field_data( *coordField , node1 ));
          double  dcoord0[3] = {coord0[0],coord0[1], (spatialDimension==2?0:coord0[2])};
          double  dcoord1[3] = {coord1[0],coord1[1], (spatialDimension==2?0:coord1[2])};

          int markInfo = markEdge(element, iSubDimOrd, node0, node1, dcoord0, dcoord1, 0);
          bool do_unref = markInfo & DO_UNREFINE;
          if (!do_unref)
            {
              unrefAllEdges = false;
              break;
            }
        }
      if (unrefAllEdges)
        return DO_UNREFINE;
      else
        return DO_NOTHING;
    }

    void IEdgeAdapter::buildUnrefineList(ElementUnrefineCollection& elements_to_unref)
    {
      //ElementUnrefineCollection elements_to_unref(*m_eMesh.get_bulk_data());
      elements_to_unref.clear();

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];

                // FIXME
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error is isParentElement)
                //const bool check_for_family_tree = false;
                //bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
                //bool hasFamilyTree = m_eMesh.hasFamilyTree(element);

                if (0)
                  {
                    const bool check_for_family_tree = false;
                    bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);

                    if (isParent)
                      continue;
                  }

                //if (hasFamilyTree && !isParent)
                  {
                    const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

                    if (elem_nodes.size())
                      {
                        int markInfo = markUnrefine(element);
                        if (markInfo & DO_UNREFINE)
                          {
                            elements_to_unref.insert(element);
                            //std::cout << "tmp unref element id= " << m_eMesh.identifier(element) << std::endl;
                            //m_eMesh.print_entity(std::cout, element);
                          }
                      }
                  }
              }
          }
        }

      //return elements_to_unref;
    }

  }

