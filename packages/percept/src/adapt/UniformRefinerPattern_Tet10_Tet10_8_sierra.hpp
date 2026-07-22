// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Tet10_Tet10_8_sierra_hpp
#define adapt_UniformRefinerPattern_Tet10_Tet10_8_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#define FACE_BREAKER_T10_T10 1
#if FACE_BREAKER_T10_T10
#include "UniformRefinerPattern_Tri6_Tri6_4_sierra.hpp"
#endif

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::Tetrahedron<10>, shards::Tetrahedron<10>, 8, SierraPort > : public URP<shards::Tetrahedron<10>, shards::Tetrahedron<10>  >
    {

      //Teuchos::RCP< UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort >  > m_face_breaker;
      //UniformRefinerPatternBase * m_face_breaker;
#if FACE_BREAKER_T10_T10
      UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Tetrahedron<10>, shards::Tetrahedron<10>  >(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if FACE_BREAKER_T10_T10

        m_face_breaker =  new UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > (eMesh, block_names);

#endif

      }
      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;
#if FACE_BREAKER_T10_T10
        bp.resize(2);
#else
        bp.resize(1);
#endif

        bp[0] = this;
#if FACE_BREAKER_T10_T10
        bp[1] = m_face_breaker;
#endif
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(3);
        // 4 vertices
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u); // 18
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank(), 3u); // 12
        needed_entities[2] = NeededEntityType(stk::topology::ELEMENT_RANK, 1u); // 1
        //setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 8; }

      //typedef std::array<unsigned, ToTopology::node_count > refined_element_type;


      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool, 
                                        proc_rank_field);
      }
    };

  }

#endif
