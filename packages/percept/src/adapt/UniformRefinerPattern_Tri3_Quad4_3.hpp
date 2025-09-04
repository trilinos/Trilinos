// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Tri3_Quad4_3_hpp
#define adapt_UniformRefinerPattern_Tri3_Quad4_3_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"

  namespace percept {

    /**
     *      Face #0:
     *
     *           2
     *           o
     *          / \
     *         /   \
     *        /  6  \
     *     5 o-._ _.-o 4
     *      /    *    \
     *     /     |     \
     *    /      |      \
     *   o-------o-------o
     *  0        3        1
     *
     *
     * Centroid node: 14
     *
     *   CHILD 4-Node Quadrilateral Object Node Maps:
     * |
     * | static const UInt child[3][4] = { {0, 3, 6, 5}, {1, 4, 6, 3}, {2, 5, 6, 4} };
     */

    template <>
    class UniformRefinerPattern< shards::Triangle<3>, shards::Quadrilateral<4>, 3, Specialization > : public URP< shards::Triangle<3>, shards::Quadrilateral<4> >
    {
      UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > * m_edge_breaker;
    public:

      virtual bool edgeMarkIsEnough() override { return false; }

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP< shards::Triangle<3>, shards::Quadrilateral<4>  >(eMesh)
       {
         EXCEPTWATCH;
         m_primaryEntityRank = eMesh.face_rank();
         if (m_eMesh.get_spatial_dim() == 2)
           m_primaryEntityRank = stk::topology::ELEMENT_RANK;

         setNeededParts(eMesh, block_names, false);

        Elem::StdMeshObjTopologies::bootstrap();

        if (m_eMesh.get_spatial_dim() == 2)
          m_edge_breaker =  new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) ;
        else
          m_edge_breaker = 0;

       }
      ~UniformRefinerPattern()
      {
#if EDGE_BREAKER_Q4_T3_6_S
        if (m_edge_breaker) delete m_edge_breaker;
#endif
      }

      virtual void doBreak() override {}

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
      {
        EXCEPTWATCH;
        bp.resize(2);

        if (eMesh.get_spatial_dim() == 2)
          {
            bp[0] = this;
            bp[1] = m_edge_breaker;
          }
        else if (eMesh.get_spatial_dim() == 3)
          {
            bp.resize(0);
            // FIXME
            //             std::cout << "ERROR" ;
            //             exit(1);
          }

      }
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(2);
        needed_entities[0].first = m_eMesh.edge_rank(); // edges have 2 nodes
        needed_entities[1].first = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK :  m_eMesh.face_rank());
        setToOne(needed_entities);
      }


      virtual unsigned getNumNewElemPerElem() override { return 3; }

      virtual StringStringMap fixSurfaceAndEdgeSetNamesMap() override
      {
        StringStringMap str_map;
        str_map["tet4"] = "hex8";
        str_map["tri3"] = "quad4";
        return str_map;
      }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        static stk::mesh::EntityId elems[3][4];

        CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

        add_parts = m_toParts;

        stk::mesh::EntityRank my_rank = m_primaryEntityRank;

        nodeRegistry.prolongateCoords(element, my_rank, 0u);
        nodeRegistry.addToExistingParts(element, my_rank, 0u);
        nodeRegistry.prolongateFields(element, my_rank, 0u);

#define CENTROID_N NN(m_primaryEntityRank, 0)

        static const unsigned child[3][4] = { {0, 3, 6, 5}, {1, 4, 6, 3}, {2, 5, 6, 4} };

        for (unsigned i_child = 0; i_child < 3; i_child++)
          {
            for (unsigned jnode = 0; jnode < 4; jnode++)
              {
                unsigned kc = 0;
                unsigned jc = child[i_child][jnode];
                if (jc <= 2)
                  {
                    kc = VERT_N(jc);
                  }
                else if (jc <= 5)
                  {
                    kc = EDGE_N(jc - 3);
                  }
                else
                  {
                    kc = CENTROID_N;
                  }
                elems[i_child][jnode] = kc;
              }
          }

#undef CENTROID_N

        bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

        for (unsigned ielem=0; ielem < 3; ielem++)
          {
            stk::mesh::Entity newElement = stk::mesh::Entity();
            if (!use_declare_element_side)
              newElement = *element_pool;

            stk::mesh::Entity nodes[4];
            for (unsigned jnode=0; jnode < 4; ++jnode)
              {
                nodes[jnode] = createOrGetNode(nodeRegistry, eMesh, elems[ielem][jnode]);
              }

            create_side_element(eMesh, use_declare_element_side, nodes, 4, newElement);

            if (proc_rank_field && m_eMesh.entity_rank(element) == stk::topology::ELEMENT_RANK)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            change_entity_parts(eMesh, element, newElement);

            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            ft_element_pool++;
            if (!use_declare_element_side)
              element_pool++;

          }

      }

    };

  }

#endif
