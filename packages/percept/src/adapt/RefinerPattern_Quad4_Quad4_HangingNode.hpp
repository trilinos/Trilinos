// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Quad4_Quad4_HangingNode_hpp
#define adapt_RefinerPattern_Quad4_Quad4_HangingNode_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Line2_Line2_N.hpp"
#include "RefinerPattern_Quad4_Het_N.hpp"

  namespace percept {

    struct QuadHangingNode {};

    template <>
    class RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadHangingNode > :
      public URP<shards::Quadrilateral<4>,shards::Quadrilateral<4>  >
    {
    protected:

      RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :
        URP<shards::Quadrilateral<4>, shards::Quadrilateral<4>  >(eMesh)
      {
        m_mark_centroid_always = true;

        m_primaryEntityRank = eMesh.face_rank();
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        if (m_eMesh.get_spatial_dim() == 2)
          m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;
        else
          m_edge_breaker = 0;
      }

      ~RefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
      {
        EXCEPTWATCH;
        bp.resize(2); // = std::vector<UniformRefinerPatternBase *>(2u, 0);

        if (eMesh.get_spatial_dim() == 2)
          {
            bp[0] = this;
            bp[1] = m_edge_breaker;
          }
        else if (eMesh.get_spatial_dim() == 3)
          {
            bp.resize(0);
          }
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
            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }
        unsigned num_faces_marked = 0;
        stk::mesh::EntityRank rank = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank());
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
      }

    };

  }

#endif
