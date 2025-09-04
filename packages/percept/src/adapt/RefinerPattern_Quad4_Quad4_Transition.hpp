// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Quad4_Quad4_Transition_hpp
#define adapt_RefinerPattern_Quad4_Quad4_Transition_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Line2_Line2_N.hpp"
#include "RefinerPattern_Quad4_Het_N.hpp"

  namespace percept {
#define LOCAL_DEBUG_PP 0

    struct QuadTransition {};

    template <>
    class RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition > :
      public URP<shards::Quadrilateral<4>,shards::Quadrilateral<4>  >
    {
    private:
      RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;
    public:
      typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 8, QuadHet > QuadTransitionRefinePatternType;
      QuadTransitionRefinePatternType *m_transition_breaker;

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType(), bool use_only_tris=true, bool avoid_centroid_node=false) :
        URP<shards::Quadrilateral<4>, shards::Quadrilateral<4>  >(eMesh)
      {
        m_primaryEntityRank = eMesh.face_rank();
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        m_transition_breaker = new QuadTransitionRefinePatternType(eMesh, block_names, use_only_tris, avoid_centroid_node) ;

        bool merge=false;
        mergeOrAddParts(m_transition_breaker, this, merge);

        if (m_eMesh.get_spatial_dim() == 2)
          {
            m_mark_centroid_always = true;
            m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;
          }
        else
          m_edge_breaker = 0;

        if (LOCAL_DEBUG_PP)
          {
            std::cout << "tmp Quad_Tri2 printParts trans breaker= \n" ;
            printParts(m_transition_breaker);

            std::cout << "tmp Quad_Tri2 printParts this= \n" ;
            printParts(this);
          }
      }
      ~RefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
        if (m_transition_breaker) delete m_transition_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        bp.resize(0);
        bp.push_back(this);
        for (unsigned ii=0; ii < m_transition_breaker->m_bp_exported.size(); ++ii)
          {
            bp.push_back(m_transition_breaker->m_bp_exported[ii]);
          }
        if (m_edge_breaker)
          bp.push_back(m_edge_breaker);
      }

      void setSubPatternsForSetNeededParts( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        bp.resize(0);
        bp.push_back(this);
        bp.push_back(m_transition_breaker);
        if (m_edge_breaker)
          bp.push_back(m_edge_breaker);
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(2);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank());
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 8; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 4; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }
        unsigned num_faces_marked = 0;
        stk::mesh::EntityRank rank = (eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : eMesh.face_rank());
        if ( new_sub_entity_nodes[rank].size() )
          {
            int iface = 0;
            num_faces_marked = new_sub_entity_nodes[rank][iface].size();
          }
        if (num_edges_marked == 4 && num_faces_marked == 1)
          {

            genericRefine_createNewElements(eMesh, nodeRegistry,
                                            element, new_sub_entity_nodes, element_pool, ft_element_pool,
                                            proc_rank_field);
          }
        else
          {
            if (s_do_transition_break && m_transition_breaker)
              m_transition_breaker->createNewElements(eMesh, nodeRegistry, element,
                                                      new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
          }
      }

    };

  }
#undef LOCAL_DEBUG_PP
#endif
