// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_ShellTri6_ShellTri6_4_sierra_hpp
#define adapt_UniformRefinerPattern_ShellTri6_ShellTri6_4_sierra_hpp


#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#define EDGE_TR6_TR6_4_BREAKER 1
#if EDGE_TR6_TR6_4_BREAKER

#include "UniformRefinerPattern_ShellLine3_ShellLine3_2_sierra.hpp"
#include "UniformRefinerPattern_Tri6_Tri6_4_sierra.hpp"
#endif

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::ShellTriangle<6>, shards::ShellTriangle<6>, 4, SierraPort > :
      public URP<shards::ShellTriangle<6>,shards::ShellTriangle<6>  >
    {


#if EDGE_TR6_TR6_4_BREAKER
      UniformRefinerPattern<shards::ShellLine<3>, shards::ShellLine<3>, 2, SierraPort > * m_edge_breaker;
      UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :
        URP<shards::ShellTriangle<6>, shards::ShellTriangle<6>  >(eMesh)
      {

        if (m_eMesh.get_spatial_dim() != 3)
          {
            throw std::runtime_error("can't refine shell elements in 2D");
          }
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if EDGE_TR6_TR6_4_BREAKER

        //m_edge_breaker = Teuchos::rcp( new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) );
        if (m_eMesh.get_spatial_dim() == 3)
          {
            m_edge_breaker = new UniformRefinerPattern<shards::ShellLine<3>, shards::ShellLine<3>, 2, SierraPort > (eMesh, block_names) ;
            m_face_breaker = new UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > (eMesh, block_names) ;
          }
#endif

      }

      ~UniformRefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
        if (m_face_breaker) delete m_face_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp.resize(0);

        if (eMesh.get_spatial_dim() == 3)
          {
            bp.push_back(this);
#if EDGE_TR6_TR6_4_BREAKER
            bp.push_back( m_face_breaker);
            bp.push_back( m_edge_breaker);
#endif
          }
        else if (eMesh.get_spatial_dim() != 3)
          {
            // FIXME
             std::cout << "ERROR" ;
             throw std::runtime_error("ERROR in shell quad class");
          }

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        if (m_eMesh.get_spatial_dim() == 2)
          {
            throw std::runtime_error("ERROR in shell quad class fillNeededEntities");
          }

        needed_entities.resize(2);
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u);
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank() , 3u);
      }

      virtual unsigned getNumNewElemPerElem() { return 4; }

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
