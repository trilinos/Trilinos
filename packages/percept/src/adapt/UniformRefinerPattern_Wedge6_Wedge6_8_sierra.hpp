// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Wedge6_Wedge6_8_sierra_hpp
#define adapt_UniformRefinerPattern_Wedge6_Wedge6_8_sierra_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Quad4_Quad4_4_sierra.hpp"
#include "UniformRefinerPattern_Tri3_Tri3_4_sierra.hpp"

#include <percept/PerceptBoostArray.hpp>

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::Wedge<6>, shards::Wedge<6>, 8, SierraPort > : public URP<shards::Wedge<6>,shards::Wedge<6>  >
    {

#define FACE_BREAKER_W6_W6_8 1
#if FACE_BREAKER_W6_W6_8
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > * m_face_breaker;
      UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > * m_face_breaker_tri;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Wedge<6>, shards::Wedge<6>  >(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if FACE_BREAKER_W6_W6_8

        m_face_breaker =  new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > (eMesh, block_names) ;
        m_face_breaker_tri = new UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > (eMesh, block_names);
#endif

      }

      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
        if (m_face_breaker_tri) delete m_face_breaker_tri;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;
        bp.resize(3);

        bp[0] = this;
#if FACE_BREAKER_W6_W6_8
        bp[1] = m_face_breaker;
        bp[2] = m_face_breaker_tri;
#endif
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(3);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = m_eMesh.face_rank();
        setToOne(needed_entities);

        int faces[5] = {1,1,1,0,0};
        needed_entities[1].third.resize(5);
        needed_entities[1].third.assign(faces,faces+5);
      }

      virtual unsigned getNumNewElemPerElem() override { return 8; }


      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
#if 0
        static bool s_not_printed = true;

        if (s_not_printed)
          {
            s_not_printed = false;
            printRefinementTopoX_Table();
          }
#endif
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool,
                                        proc_rank_field);
      }

    };

  }

#endif
