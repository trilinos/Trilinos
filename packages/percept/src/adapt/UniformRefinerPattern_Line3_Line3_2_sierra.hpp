// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Line3_Line3_2_sierra_hpp
#define adapt_UniformRefinerPattern_Line3_Line3_2_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::Line<3>, shards::Line<3>, 2, SierraPort > : public URP<shards::Line<3>, shards::Line<3>  >
    {
    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Line<3>, shards::Line<3>  >(eMesh)
      {
        m_primaryEntityRank = m_eMesh.edge_rank();
        if (m_eMesh.get_spatial_dim() == 1)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[0].second = 3u;
      }

      virtual unsigned getNumNewElemPerElem() { return 2; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool, 
                                        proc_rank_field);
      }

    };

  }

#endif
