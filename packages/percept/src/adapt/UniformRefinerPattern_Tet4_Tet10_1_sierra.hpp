// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Tet4_Tet10_1_sierra_hpp
#define adapt_UniformRefinerPattern_Tet4_Tet10_1_sierra_hpp


#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <percept/PerceptBoostArray.hpp>

#define FACE_BREAKER_TET4_TET10_1 1
#if FACE_BREAKER_TET4_TET10_1
#include "UniformRefinerPattern_Tri3_Tri6_1_sierra.hpp"
#endif

  namespace percept {

    template <>
    class UniformRefinerPattern< shards::Tetrahedron<4>, shards::Tetrahedron<10>, 1, SierraPort > : public URP<shards::Tetrahedron<4> , shards::Tetrahedron<10> >
    {

#if FACE_BREAKER_TET4_TET10_1
      UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<6>, 1, SierraPort > * m_subDim_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Tetrahedron<4> , shards::Tetrahedron<10> >(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, false);
        Elem::StdMeshObjTopologies::bootstrap();
#if FACE_BREAKER_TET4_TET10_1

        m_subDim_breaker =  new UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<6>, 1, SierraPort > (eMesh, block_names) ;
#endif

      }

      ~UniformRefinerPattern()
      {
        if (m_subDim_breaker) delete m_subDim_breaker;
      }


      virtual void doBreak() override {}

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;
        bp.resize(2);

        bp[0] = this;
#if FACE_BREAKER_TET4_TET10_1
        bp[1] = m_subDim_breaker;
#endif

      }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 1; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        genericEnrich_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool,
                                        proc_rank_field);
      }

    };

  } // namespace percept

#endif
