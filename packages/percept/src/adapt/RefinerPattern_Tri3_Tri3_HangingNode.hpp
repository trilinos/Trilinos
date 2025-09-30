// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Tri3_Tri3_HangingNode_sierra_hpp
#define adapt_RefinerPattern_Tri3_Tri3_HangingNode_sierra_hpp


#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Line2_Line2_N.hpp"
#include "RefinerPattern_Tri3_Tri3_N.hpp"

  namespace percept {

    extern bool s_do_transition_break;

    struct TriHangingNode {};

    template <>
    class RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1, TriHangingNode > : public URP<shards::Triangle<3>, shards::Triangle<3>  >
    {

      RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;
      RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 >  * m_transition_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Triangle<3>, shards::Triangle<3>  >(eMesh),
                                                                                                    m_edge_breaker (0), m_transition_breaker(0)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        if (m_eMesh.get_spatial_dim() == 2)
          {
            m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;
          }
        m_transition_breaker = new RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 >  (eMesh, block_names);

      }

      ~RefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
        if (m_transition_breaker) delete m_transition_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
      {
        EXCEPTWATCH;
        bp.resize(2);

        if (eMesh.get_spatial_dim() == 2)
          {
            bp[0] = this;
            bp[1] = m_edge_breaker;
          }
        else if (eMesh.get_spatial_dim() == 3)
          {
            bp.resize(0);
          }

      }

      void setSubPatternsForSetNeededParts( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
      {
        EXCEPTWATCH;
        bp.resize(2);

        if (eMesh.get_spatial_dim() == 2)
          {
            bp[0] = this;
            bp[1] = m_edge_breaker;
          }
        else if (eMesh.get_spatial_dim() == 3)
          {
            bp.resize(0);
          }
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 4; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 3; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }
        if (num_edges_marked != 3)
          {
            if (s_do_transition_break)
              m_transition_breaker->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
            return;
          }
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool,
                                        proc_rank_field);
      }

    };

  }

#endif
