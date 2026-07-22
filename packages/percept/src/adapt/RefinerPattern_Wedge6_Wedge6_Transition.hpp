// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Wedge6_Wedge6_Transition_hpp
#define adapt_RefinerPattern_Wedge6_Wedge6_Transition_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Quad4_Quad4_Transition.hpp"
#include "RefinerPattern_Tri3_Tri3_HangingNode.hpp"
#include "RefinerPattern_Wedge6_Het_N.hpp"

  namespace percept {

    struct WedgeTransition {};

    template <>
    class RefinerPattern<shards::Wedge<6>, shards::Wedge<6>, -1, WedgeTransition > :
      public URP<shards::Wedge<6>,shards::Wedge<6>  >
    {
    private:
      typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition >  FaceBreakerType;
      FaceBreakerType * m_face_breaker;
      typedef RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1, TriHangingNode >  TriFaceBreakerType;
      TriFaceBreakerType * m_face_breaker_tri;
      typedef RefinerPattern<shards::Wedge<6>, shards::Wedge<6>, 20, WedgeHet > WedgeTransitionRefinePatternType;
      WedgeTransitionRefinePatternType *m_transition_breaker;
      typedef UniformRefinerPattern<shards::Wedge<6>, shards::Wedge<6>, 8, SierraPort > UniformWedgeRefinePatternType;
      UniformWedgeRefinePatternType *m_uniform_refiner;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :
        URP<shards::Wedge<6>, shards::Wedge<6>  >(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        m_transition_breaker = new WedgeTransitionRefinePatternType(eMesh, block_names) ;
        m_uniform_refiner = new UniformWedgeRefinePatternType(eMesh, block_names) ;

        bool merge=false;
        mergeOrAddParts(m_transition_breaker, this, merge);

        m_mark_centroid_always = true;

        m_face_breaker =  new FaceBreakerType(eMesh, block_names) ;
        m_face_breaker_tri = new TriFaceBreakerType(eMesh, block_names) ;
        if (0)
          {
            std::cout << "tmp Wedge_ printParts this= \n" ;
            printParts(this);
          }
      }
      ~RefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
        if (m_face_breaker_tri) delete m_face_breaker_tri;
        if (m_transition_breaker) delete m_transition_breaker;
        if (m_uniform_refiner) delete m_uniform_refiner;
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
        bp.push_back(m_face_breaker_tri);
      }

      void setSubPatternsForSetNeededParts( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        //m_transition_breaker->setSubPatterns(bp, eMesh);
        bp.resize(0);
        bp.push_back(this);
        bp.push_back(m_transition_breaker);
        bp.push_back(m_face_breaker);
        bp.push_back(m_face_breaker_tri);
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(3);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = m_eMesh.face_rank();
        needed_entities[2].first = stk::topology::ELEMENT_RANK;
        setToOne(needed_entities);

        int faces[5] = {1,1,1,0,0};
        needed_entities[1].third.resize(5);
        needed_entities[1].third.assign(faces,faces+5);
      }

      // FIXME - too many here...
      virtual unsigned getNumNewElemPerElem() override { return 32; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 9; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }
        unsigned num_faces_marked = 0;
        stk::mesh::EntityRank rank = eMesh.face_rank();
        for (int iface = 0; iface < 5; ++iface)
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

        //if (num_edges_marked == 9 && num_faces_marked == 5)
        if (num_edges_marked == 9)
          {
            m_uniform_refiner->createNewElements(eMesh, nodeRegistry, element,
                                                 new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
          }
        else
          {
            if (s_do_transition_break && m_transition_breaker)
              {
                m_transition_breaker->createNewElements(eMesh, nodeRegistry, element,
                                                        new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
              }
          }
      }

    };

  }

#endif
