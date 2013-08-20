#ifndef stk_adapt_NodeRegistryDef_hpp
#define stk_adapt_NodeRegistryDef_hpp

#include <stk_adapt/NodeRegistry.hpp>

namespace stk {
  namespace adapt {

    NodeIdsOnSubDimEntityType* NodeRegistry::getNewNodesOnSubDimEntity(const stk::mesh::Entity element,  stk::mesh::EntityRank& needed_entity_rank, 
                                                                       unsigned iSubDimOrd)
    {
      EXCEPTWATCH;
      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
      static  SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;

      if (is_empty)
        {
          if (0)
            {
              const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
              CellTopology cell_topo(cell_topo_data);

              std::cout << "NodeRegistry::getNewNodesOnSubDimEntity: no node found, cell_topo = " << cell_topo.getName()
                        << "\n subDimEntity= " << subDimEntity 
                        << "\n element= " << element 
                        << "\n eMesh.entity_rank(element) = " << m_eMesh.entity_rank(element)
                        << "\n needed_entity_rank= " << needed_entity_rank
                        << "\n iSubDimOrd= " << iSubDimOrd << std::endl;
              throw std::runtime_error("NodeRegistry::getNewNodesOnSubDimEntity: no node found");

            }
          return 0;
        }
      NodeIdsOnSubDimEntityType& nodeId = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
      return &nodeId;
    }

    bool NodeRegistry::is_empty( const stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
    {
      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
      static SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      //SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;
      return is_empty;
    }


  }
}
#endif
