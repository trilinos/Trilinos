// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_NodeRegistryDef_hpp
#define adapt_NodeRegistryDef_hpp

#include <adapt/NodeRegistry.hpp>

  namespace percept {

    NodeIdsOnSubDimEntityType* NodeRegistry::getNewNodesOnSubDimEntity(const stk::mesh::Entity element,  stk::mesh::EntityRank& needed_entity_rank,
                                                                       unsigned iSubDimOrd)
    {
      EXCEPTWATCH;
      SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
      getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
      static  SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnerId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnerId = (nodeId_elementOwnerId_ptr ? *nodeId_elementOwnerId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnerId_ptr == 0;

      bool debug = false;
      if (debug)
        {
          const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
          CellTopology cell_topo(cell_topo_data);

          std::cout << "NodeRegistry::getNewNodesOnSubDimEntity: cell_topo = " << cell_topo.getName()
                    << "\n subDimEntity= " << subDimEntity
                    << "\n element= " << m_eMesh.id(element)
                    << "\n eMesh.entity_rank(element) = " << m_eMesh.entity_rank(element)
                    << "\n needed_entity_rank= " << needed_entity_rank
                    << "\n iSubDimOrd= " << iSubDimOrd
                    << "\n nodeId_elementOwnerId.size= " << std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId).size()
                    << "\n node_valid= " << m_eMesh.is_valid(std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId)[0])
                    << "\n node_id= " << std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId).m_entity_id_vector[0]
                    << "\n is_empty= " << is_empty
                    << "\n owningElementKey= " << (is_empty ? stk::mesh::EntityKey(m_eMesh.node_rank(),0) : std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnerId))
                    << "\n owningElementSubDimRank= " << (is_empty ? 0 : std::get<SDC_DATA_OWNING_SUBDIM_RANK>(nodeId_elementOwnerId))
                    << "\n owningElementSubDimOrd= " << (is_empty ? 0 : std::get<SDC_DATA_OWNING_SUBDIM_ORDINAL>(nodeId_elementOwnerId))
                    << std::endl;
        }

      if (is_empty)
        {
          return 0;
        }
      NodeIdsOnSubDimEntityType& nodeId = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
      if (debug)
        {
          std::cout << "NodeRegistry::getNewNodesOnSubDimEntity: nodeId.size= " << nodeId.size() << " node= " << nodeId[0] << " " << nodeId.m_entity_id_vector[0] << std::endl;
        }
      return &nodeId;
    }

    bool NodeRegistry::is_empty( const stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
    {
      SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
      getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
      static SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnerId_ptr = getFromMapPtr(subDimEntity);
      bool is_empty = nodeId_elementOwnerId_ptr == 0;
      return is_empty;
    }


  }

#endif
