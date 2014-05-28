#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <unit_tests/TestLocalRefinerTri2.hpp>

//#define STK_PERCEPT_HAS_GEOMETRY
#undef STK_PERCEPT_HAS_GEOMETRY
#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <stk_adapt/geometry/GeometryKernelOpenNURBS.hpp>
#include <stk_adapt/geometry/MeshGeometry.hpp>
#include <stk_adapt/geometry/GeometryFactory.hpp>
#endif

// FIXME
// #include <stk_mesh/baseImpl/EntityImpl.hpp>
// #include <stk_mesh/base/Entity.hpp>
// FIXME

namespace stk {
  namespace adapt {


    TestLocalRefinerTri2::TestLocalRefinerTri2(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field, bool diagonals) : 
      Refiner(eMesh, bp, proc_rank_field), m_diagonals(diagonals)
    {
    }

    void TestLocalRefinerTri2::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks)
    {
      const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                
      CellTopology cell_topo(cell_topo_data);
      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

      VectorFieldType* coordField = m_eMesh.get_coordinates_field();

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
            }
          else if (needed_entity_rank == m_eMesh.element_rank())
            {
              numSubDimNeededEntities = 1;
            }

          // see how many edges are already marked
          int num_marked=0;
          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
                  if (!is_empty) ++num_marked;
                }
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
              //SubDimCell_SDSEntityType subDimEntity;
              //getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
              //bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
              //if(1||!is_empty)

              if (needed_entity_rank == m_eMesh.edge_rank())
                {
                  stk::mesh::Entity & node0 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk::mesh::Entity & node1 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 = stk::mesh::field_data( *coordField , node0 );
                  double * const coord1 = stk::mesh::field_data( *coordField , node1 );

                  // only refine diagonals of the split quads
                  if (m_diagonals)
                    {
                  if ( std::abs(coord0[0]-coord1[0]) > 1.e-3 && std::abs(coord0[1]-coord1[1]) > 1.e-3 )
                    {
                      (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true);
                    }
                    }
                  else
                    {
                      if ( std::abs(coord0[0]-coord1[0]) < 1.e-3 && std::abs(coord0[1]-coord1[1]) > 1.e-3 )
                        {
                          (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true);
                        }
                    }
                }

            } // iSubDimOrd
        } // ineed_ent
    }


  } // namespace adapt
} // namespace stk
