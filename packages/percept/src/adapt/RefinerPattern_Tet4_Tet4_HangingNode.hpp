// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Tet4_Tet4_HangingNode_sierra_hpp
#define adapt_RefinerPattern_Tet4_Tet4_HangingNode_sierra_hpp


#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Line2_Line2_N.hpp"
#include "RefinerPattern_Tri3_Tri3_HangingNode.hpp"
#include "RefinerPattern_Tet4_Tet4_N.hpp"

  namespace percept {

    extern bool s_do_transition_break;

    struct TetHangingNode {};

    template <>
    class RefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, -1, TetHangingNode > : public URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >
    {

      RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;
      //RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > * m_face_breaker;
      RefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,-1, TriHangingNode > *m_face_breaker;

    public:
      RefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, -1 >  * m_transition_breaker;

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >(eMesh),
                                                                                                    m_edge_breaker (0), m_face_breaker(0), m_transition_breaker(0)
      {
        m_primaryEntityRank = m_eMesh.element_rank();
        setNeededParts(eMesh, block_names, true);

        if (0)
          {
            std::cout << "tmp TetHangingNode=\n" ;
            printParts(this, true);
          }

        Elem::StdMeshObjTopologies::bootstrap();

          {
            m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;
            //m_face_breaker =  new RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > (eMesh, block_names) ;
            m_face_breaker = new RefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,-1, TriHangingNode > (eMesh, block_names);
            m_transition_breaker = new RefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, -1 >  (eMesh, block_names);

            // repeat calls to setNeededParts to catch newly created parts
            bool skipConvertedParts = false;
            m_transition_breaker->setNeededParts(eMesh, block_names, true, skipConvertedParts);

          }

      }

      ~RefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
        if (m_face_breaker) delete m_face_breaker;
        if (m_transition_breaker) delete m_transition_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp.resize(3); //= std::vector<UniformRefinerPatternBase *>(3u, 0);

        bp[0] = this;
        bp[1] = m_face_breaker;
        bp[2] = m_edge_breaker;
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        if (!m_mark_centroid_always)
          {
            needed_entities.resize(1);
            needed_entities[0].first = m_eMesh.edge_rank();
            //needed_entities[1].first = m_eMesh.face_rank();
            setToOne(needed_entities);
          }
        else
          {
            needed_entities.resize(2);
            needed_entities[0].first = m_eMesh.edge_rank();
            needed_entities[1].first = m_eMesh.face_rank();
            setToOne(needed_entities);
          }
      }

      virtual unsigned getNumNewElemPerElem() { return 8; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        m_transition_breaker->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
      }

    };

  }

#endif
