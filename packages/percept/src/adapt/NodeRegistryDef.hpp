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

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;

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
                    << "\n nodeId_elementOwnderId.size= " << nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().size()
                    << "\n node_valid= " << m_eMesh.is_valid(nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>()[0])
                    << "\n node_id= " << nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().m_entity_id_vector[0]
                    << "\n is_empty= " << is_empty
                    << "\n owningElementKey= " << (is_empty ? stk::mesh::EntityKey(m_eMesh.node_rank(),0) : nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>())
                    << "\n owningElementSubDimRank= " << (is_empty ? 0 : nodeId_elementOwnderId.get<SDC_DATA_OWNING_SUBDIM_RANK>())
                    << "\n owningElementSubDimOrd= " << (is_empty ? 0 : nodeId_elementOwnderId.get<SDC_DATA_OWNING_SUBDIM_ORDINAL>())
                    << std::endl;
        }

      if (is_empty)
        {
          return 0;
        }
      NodeIdsOnSubDimEntityType& nodeId = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
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

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      //SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;
      return is_empty;
    }


  }

#endif
