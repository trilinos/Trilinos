// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Hex8_Hex8_Transition_hpp
#define adapt_RefinerPattern_Hex8_Hex8_Transition_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Quad4_Quad4_Transition.hpp"
#include "RefinerPattern_Hex8_Het_N.hpp"

#define LOCAL_DEBUG 0

  namespace percept {

    struct HexTransition {};

    template <>
    class RefinerPattern<shards::Hexahedron<8>, shards::Hexahedron<8>, -1, HexTransition > :
      public URP<shards::Hexahedron<8>,shards::Hexahedron<8>  >
    {
    private:
      //RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;
      typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition >  FaceBreakerType;
      FaceBreakerType * m_face_breaker;
      typedef RefinerPattern<shards::Hexahedron<8>, shards::Hexahedron<8>, 24, HexHet > HexTransitionRefinePatternType;
      HexTransitionRefinePatternType *m_transition_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :
        URP<shards::Hexahedron<8>, shards::Hexahedron<8>  >(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        m_transition_breaker = new HexTransitionRefinePatternType(eMesh, block_names) ;

        bool merge=false;
        mergeOrAddParts(m_transition_breaker, this, merge);

        m_mark_centroid_always = true;

        m_face_breaker =  new FaceBreakerType(eMesh, block_names) ;

        if (LOCAL_DEBUG)
          {
            std::cout << "tmp Hex_ printParts this= \n" ;
            UniformRefinerPatternBase::printParts(this);
          }
      }
      ~RefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
        if (m_transition_breaker) delete m_transition_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        //m_transition_breaker->setSubPatterns(bp, eMesh);
        bp.resize(0);
        bp.push_back(this);
        for (unsigned ii=0; ii < m_transition_breaker->m_bp_exported.size(); ++ii)
          {
            bp.push_back(m_transition_breaker->m_bp_exported[ii]);
          }
        bp.push_back(m_face_breaker);
      }


      void setSubPatternsForSetNeededParts( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        //m_transition_breaker->setSubPatterns(bp, eMesh);
        bp.resize(0);
        bp.push_back(this);
        bp.push_back(m_transition_breaker);
        bp.push_back(m_face_breaker);
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(3);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = m_eMesh.face_rank();
        needed_entities[2].first = stk::topology::ELEMENT_RANK;
        setToOne(needed_entities);
      }

      // FIXME - too many here...
      virtual unsigned getNumNewElemPerElem() override { return 48; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 12; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }
        unsigned num_faces_marked = 0;
        stk::mesh::EntityRank rank = eMesh.face_rank();
        for (int iface = 0; iface < 6; ++iface)
          {
            unsigned num_nodes_on_face = new_sub_entity_nodes[rank][iface].size();
            if (num_nodes_on_face)
              ++num_faces_marked;
          }
        int num_centroid_marked = 0;
        if (new_sub_entity_nodes[eMesh.element_rank()][0].size())
          ++num_centroid_marked;

        if (0)
          {
            std::cout << "tmp srk b1 elem = " << eMesh.identifier(element) << " topo= " << eMesh.bucket(element).topology()
                      << " num_edges_marked= " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                      << " num_centroid_marked= " << num_centroid_marked << std::endl;
          }
        if (num_edges_marked == 0)
          return;
        //if (num_edges_marked == 12 && num_faces_marked == 6 && num_centroid_marked == 1)
        if (num_edges_marked == 12)
          {
            genericRefine_createNewElements(eMesh, nodeRegistry,
                                            element, new_sub_entity_nodes, element_pool, ft_element_pool,
                                            proc_rank_field);
          }
        else
          {
            if (s_do_transition_break && m_transition_breaker)
              {
                //m_transition_breaker->m_ep_begin = m_ep_begin;
                //m_transition_breaker->m_ep_end = m_ep_end;

                m_transition_breaker->createNewElements(eMesh, nodeRegistry, element,
                                                        new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
              }
          }
      }

    };

  }
#undef LOCAL_DEBUG

#endif
