// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Quad4_Quad9_1_sierra_hpp
#define adapt_UniformRefinerPattern_Quad4_Quad9_1_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <percept/PerceptBoostArray.hpp>

#define EDGE_BREAKER_Q4_Q9_1 1
#if EDGE_BREAKER_Q4_Q9_1
#include "UniformRefinerPattern_Line2_Line3_1_sierra.hpp"
#endif

  namespace percept {

    template <>
    class UniformRefinerPattern< shards::Quadrilateral<4>, shards::Quadrilateral<9>, 1, SierraPort > : public URP<shards::Quadrilateral<4> , shards::Quadrilateral<9> >
    {

#if EDGE_BREAKER_Q4_Q9_1
      UniformRefinerPattern<shards::Line<2>, shards::Line<3>, 1, SierraPort > * m_edge_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Quadrilateral<4> , shards::Quadrilateral<9> >(eMesh)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, false);
        Elem::StdMeshObjTopologies::bootstrap();
#if EDGE_BREAKER_Q4_Q9_1

        m_edge_breaker =  new UniformRefinerPattern<shards::Line<2>, shards::Line<3>, 1, SierraPort > (eMesh, block_names) ;
#endif

      }

      ~UniformRefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
      }

      virtual void doBreak() {}

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp.resize(2);

        if (eMesh.get_spatial_dim() == 2)
          {
            bp[0] = this;
#if EDGE_BREAKER_Q4_Q9_1
            bp[1] = m_edge_breaker;
#endif
          }
        else if (eMesh.get_spatial_dim() == 3)
          {
            bp.resize(0);
          }

      }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(2);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank());
        setToOne(needed_entities);

      }

      virtual unsigned getNumNewElemPerElem() { return 1; }

      virtual StringStringMap fixSurfaceAndEdgeSetNamesMap()
      {
        StringStringMap str_map;
        //str_map["hex8"] = "tet4";
        str_map["quad4"] = "quad9";
        return str_map;
      }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        genericEnrich_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool, 
                                        proc_rank_field);
        // write a diagram of the refinement pattern as a vtk file, or a latex/tikz/pgf file
#define WRITE_DIAGRAM 0
#if WRITE_DIAGRAM

#endif
      }
    };

  }

#endif
