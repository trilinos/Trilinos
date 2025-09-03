// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Quad9_Quad9_4_sierra_hpp
#define adapt_UniformRefinerPattern_Quad9_Quad9_4_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#define EDGE_BREAKER_Q9_Q9 1
#if EDGE_BREAKER_Q9_Q9
#include "UniformRefinerPattern_Line3_Line3_2_sierra.hpp"
#endif

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::Quadrilateral<9>, shards::Quadrilateral<9>, 4, SierraPort > : public URP<shards::Quadrilateral<9>,shards::Quadrilateral<9>  >
    {

      //Teuchos::RCP< UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort >  > m_edge_breaker;
      //UniformRefinerPatternBase * m_edge_breaker;
#if EDGE_BREAKER_Q9_Q9
      UniformRefinerPattern<shards::Line<3>, shards::Line<3>, 2, SierraPort > * m_edge_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Quadrilateral<9>, shards::Quadrilateral<9>  >(eMesh)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if EDGE_BREAKER_Q9_Q9
        if (m_eMesh.get_spatial_dim() == 2)
          m_edge_breaker =  new UniformRefinerPattern<shards::Line<3>, shards::Line<3>, 2, SierraPort > (eMesh, block_names) ;
        else
          m_edge_breaker = 0;
#endif

      }
      ~UniformRefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
      {
        EXCEPTWATCH;
#if EDGE_BREAKER_Q9_Q9
        bp.resize(2);
#else
        bp.resize(1);
#endif

        if (eMesh.get_spatial_dim() == 2)
          {
            bp[0] = this;
#if EDGE_BREAKER_Q9_Q9
            bp[1] = m_edge_breaker;
#endif
          }
        else if (eMesh.get_spatial_dim() == 3)
          {
            bp.resize(0);
            // FIXME
            //             std::cout << "ERROR" ;
            //             exit(1);
          }

      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(2);
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u);
        needed_entities[1] = NeededEntityType( (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank()), 9u);
        //setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 4; }

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
