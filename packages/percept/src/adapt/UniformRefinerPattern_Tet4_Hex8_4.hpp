// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Tet4_Hex8_4_hpp
#define adapt_UniformRefinerPattern_Tet4_Hex8_4_hpp

#include "UniformRefinerPattern.hpp"
#include "UniformRefinerPattern_Tri3_Quad4_3.hpp"

  namespace percept {

    /**
     *
     *                         PARENT 4-Node Tetrahedron Object Nodes
     *              3
     *              o
     *             /|\
     *            / | \       (PARENT) 4-Node Tetrahedron Object
     *           /  |  \               Edge Node Map:
     *          /   |   \
     *         /    |    \    { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} };
     *      0 o-----|-----o 2
     *         \    |    /
     *          \   |   /
     *           \  |  /
     *            \ | /
     *             \|/
     *              o
     *              1
     *
     *     After refinement (new nodes = *):
     *
     *              3
     *              o
     *             /|\  .12 (0,2,3)
     *            / | \
     *         7 *  |  * 9
     *          /10 |   \
     *         / . 6| .11\
     *      0 o----*|-----o 2
     *         \    *8   /        Face nodes:
     *          \   |.13/         
     *         4 *  |  * 5        10 {0,1,3}
     *            \ | /           11 {1,2,3}
     *             \|/            12 {0,2,3}
     *              o             13 {0,1,2} 
     *              1
     *
     *      Face #3:               Face #0:              Face #1:                Face #2:
     *
     *           1                      3                     3                       3
     *           o                      o                     o                       o
     *          / \                    / \                   / \                     / \
     *         / B \                  /   \                 /   \                   / B \      B = Base of Hex
     *        /  13 \                /  10 \               /  11 \                 /  12 \
     *     4 o-._ _.-o 5          7 o-._ _.-o 8         8 o-._ _.-o 9           9 o-._ _.-o 7
     *      /    *    \            /    *    \           /    *    \             /    *    \
     *     /     |     \          /  B  |     \         /     | B   \           /     |     \
     *    /      |      \        /      |      \       /      |      \         /      |      \
     *   o-------o-------o      o-------o-------o     o-------o-------o       o-------o-------o
     *  0        6        2    0        4        1   1        5        2     2        6        0
     *
     *
     * Centroid node: 14
     *
     *   CHILD 8-Node Hexahedron Object Node Maps:
     * |
     * | static const UInt child_0[] = { 0, 7, 10, 4,  6, 12, 14, 13 };
     * | static const UInt child_1[] = { 1, 5, 13, 4,  8, 11, 14, 10 };
     * | static const UInt child_2[] = { 2, 5, 11, 9,  6, 13, 14, 12 };
     * | static const UInt child_3[] = { 3, 7, 12, 9,  8, 10, 14, 11 };
     * |
     *
     */


    template <>
    class UniformRefinerPattern< shards::Tetrahedron<4>, shards::Hexahedron<8>, 4 > : public URP<shards::Tetrahedron<4>, shards::Hexahedron<8>  >
    {
    private:
      UniformRefinerPattern< shards::Triangle<3>, shards::Quadrilateral<4>, 3, Specialization > *m_face_breaker;

    public:
      virtual bool edgeMarkIsEnough() override { return false; }

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP< shards::Tetrahedron<4>, shards::Hexahedron<8> >(eMesh)
       {
         EXCEPTWATCH;

         m_primaryEntityRank = stk::topology::ELEMENT_RANK;

         setNeededParts(eMesh, block_names, false);

         m_face_breaker = new UniformRefinerPattern< shards::Triangle<3>, shards::Quadrilateral<4>, 3, Specialization > (eMesh, block_names);
       }
      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }


      virtual void doBreak() override {}

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(3);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = m_eMesh.face_rank();
        needed_entities[2].first = m_eMesh.element_rank();
        setToOne(needed_entities);
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;
        bp.resize(2);
        bp[0] = this;
        bp[1] = m_face_breaker;
      }

      virtual unsigned getNumNewElemPerElem() override { return 4u; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        EXCEPTWATCH;
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        static stk::mesh::EntityId elems[4][8];

        shards::CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);


        // FIXME - maybe the computation of node coorinates should go in the calling code?
        double tmp_x[3];
        for (unsigned i_face = 0; i_face < 4; i_face++)
          {
            double * pts[3] = {FACE_COORD(i_face, 0), FACE_COORD(i_face, 1), FACE_COORD(i_face, 2)};

            double * mp = getCentroid(pts, 3, eMesh.get_spatial_dim(), tmp_x);
            (void)mp;
#if 0
            std::cout << "pts = \n"

                      << pts[0][0] << " " << pts[0][1] << " " << pts[0][2] << " \n"
                      << pts[1][0] << " " << pts[1][1] << " " << pts[1][2] << " \n"
                      << pts[2][0] << " " << pts[2][1] << " " << pts[2][2] << " \n"
                      << pts[3][0] << " " << pts[3][1] << " " << pts[3][2] << " \n"
                      << " centroid = "
                      << mp[0] << " " << mp[1] << " " << mp[2] <<  std::endl;
#endif

            //stk::mesh::Entity new_node =eMesh.createOrGetNode(FACE_N(i_face), mp);

            //std::cout << "P[" << eMesh.get_rank() << "] tmp srk createOrGetNode 1 face node = " << FACE_N(i_face) << std::endl;
            //stk::mesh::Entity node = createOrGetNode(nodeRegistry, eMesh, FACE_N(i_face));
            eMesh.createOrGetNode(FACE_N(i_face), mp);
            nodeRegistry.addToExistingParts(element, m_eMesh.face_rank(), i_face);
            nodeRegistry.prolongateFields(element, m_eMesh.face_rank(), i_face);

          }

        nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);
        nodeRegistry.addToExistingParts(element, stk::topology::ELEMENT_RANK, 0u);
        nodeRegistry.prolongateFields(element, stk::topology::ELEMENT_RANK, 0u);

        //#define C 14

// new_sub_entity_nodes[i][j]
#define CENTROID_N NN(stk::topology::ELEMENT_RANK, 0)

        static const unsigned child_0[] = { 0, 7, 10, 4,  6, 12, 14, 13 };
        static const unsigned child_1[] = { 1, 5, 13, 4,  8, 11, 14, 10 };
        static const unsigned child_2[] = { 2, 5, 11, 9,  6, 13, 14, 12 };
        static const unsigned child_3[] = { 3, 7, 12, 9,  8, 10, 14, 11 };

        static const unsigned *child[4] = {child_0, child_1, child_2, child_3 };

        for (unsigned i_child = 0; i_child < 4; i_child++)
          {
            for (unsigned jnode = 0; jnode < 8; jnode++)
              {
                unsigned kc = 0;
                unsigned jc = child[i_child][jnode];
                if (jc <= 3) 
                  {
                    kc = VERT_N(jc);
                  }
                else if (jc <= 9)
                  {
                    kc = EDGE_N(jc - 4);
                  }
                else if (jc <= 13)
                  {
                    kc = FACE_N(jc - 10);
                  }
                else
                  {
                    kc = CENTROID_N;
                  }
                elems[i_child][jnode] = kc;

                //std::cout << "P[" << eMesh.get_rank() << "] tmp srk createOrGetNode id= " << kc << " jc= " << jc << std::endl;
                //createOrGetNode(nodeRegistry, eMesh, kc);

              }
          }

#undef CENTROID_N


        for (unsigned ielem=0; ielem < 4; ielem++)
          {
            stk::mesh::Entity newElement = *element_pool;

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            {
              EXCEPTWATCH;
              change_entity_parts(eMesh, element, newElement);
            }

            {
              if (!elems[ielem][0])
                {
                  std::cout << "P[" << eMesh.get_rank() << " nid = 0 << " << std::endl;
                  exit(1);
                }
            }
            for (unsigned jnode=0; jnode < 8; jnode++)
              {
                //std::cout << "P[" << eMesh.get_rank() << "] tmp srk createOrGetNode 2" << std::endl;
                //eMesh.get_bulk_data()->declare_relation(newElement, createOrGetNode(nodeRegistry, eMesh, elems[ielem][jnode]), jnode);
                eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem][jnode]), jnode);
              }

            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            ft_element_pool++;
            element_pool++;

          }

      }

    };

  }

#endif
