#include <stk_adapt/NodeRegistry.hpp>

#if 0 && NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE
namespace Teuchos
{
  //template <class T> int hashCode(const T& x);

  template <> int hashCode(const SDCell_HashTable_Key& x)
  {
    return (int)x.getHash();
  }

}
#endif

namespace stk {
  namespace adapt {

#if !NODE_REGISTRY_MAP_ACCESSORS_INLINED
    SubDimCellData& NodeRegistry::getFromMap(SubDimCell_EntityId& subDimEntity)
    {
      //ftest(subDimEntity);
      return m_cell_2_data_map[subDimEntity];
    }
    void NodeRegistry::putInMap(SubDimCell_EntityId& subDimEntity, SubDimCellData& data)
    {
      //m_cell_2_data_map.insert(std::subDimEntity, data);
      //SubDimCellData& dataInMap = m_cell_2_data_map[subDimEntity];
      //dataInMap = data;
      //ftest(subDimEntity);
      m_cell_2_data_map[subDimEntity] = data;
    }
#endif

      /// fill
      ///    @param subDimEntity with the stk::mesh::EntityId's of
      ///    the ordinal @param iSubDimOrd sub-dimensional entity of
      ///    @param element of rank
      ///    @param needed_entity_rank
      ///
    void NodeRegistry::
      noInline_getSubDimEntity(SubDimCell_SDSEntityType& subDimEntity, const stk::mesh::Entity& element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        subDimEntity.clear();
        // in the case of elements, we don't share any nodes so we just make a map of element id to node
        if (needed_entity_rank == m_eMesh.element_rank())
          {
            subDimEntity.insert( const_cast<stk::mesh::Entity*>(&element) );
            //!!subDimEntity.insert(element.identifier());
            return;
          }

        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);

        //CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

        const unsigned *  inodes = 0;
        unsigned nSubDimNodes = 0;
        static const unsigned edge_nodes_2[2] = {0,1};
        static const unsigned face_nodes_3[3] = {0,1,2};
        static const unsigned face_nodes_4[4] = {0,1,2,3};

        // special case for faces in 3D
        if (needed_entity_rank == m_eMesh.face_rank() && needed_entity_rank == element.entity_rank())
          {
            nSubDimNodes = cell_topo_data->vertex_count;

            // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
            if (nSubDimNodes ==3 )
              inodes = face_nodes_3;
            else
              inodes = face_nodes_4;

          }
        // special case for edges in 2D
        else if (needed_entity_rank == m_eMesh.edge_rank() && needed_entity_rank == element.entity_rank())
          {
            nSubDimNodes = cell_topo_data->vertex_count;

            if (nSubDimNodes == 2 )
              {
                inodes = edge_nodes_2;
              }
            else
              {
                throw std::runtime_error("NodeRegistry bad for edges");
              }
          }
        else if (needed_entity_rank == m_eMesh.edge_rank())
          {
            inodes = cell_topo_data->edge[iSubDimOrd].node;
            nSubDimNodes = 2;
          }
        else if (needed_entity_rank == m_eMesh.face_rank())
          {
            nSubDimNodes = cell_topo_data->side[iSubDimOrd].topology->vertex_count;
            // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
            inodes = cell_topo_data->side[iSubDimOrd].node;
          }

        //subDimEntity.reserve(nSubDimNodes);
        for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
          {
            subDimEntity.insert( elem_nodes[inodes[jnode]].entity() );
            //!!subDimEntity.insert(elem_nodes[inodes[jnode]].entity()->identifier());
          }

      }


  }
}
