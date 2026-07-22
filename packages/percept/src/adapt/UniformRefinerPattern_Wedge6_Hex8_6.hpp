// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Wedge6_Hex8_6_hpp
#define adapt_UniformRefinerPattern_Wedge6_Hex8_6_hpp

#include "UniformRefinerPattern.hpp"
#include "UniformRefinerPattern_Tri3_Quad4_3.hpp"
#include "UniformRefinerPattern_Quad4_Quad4_4_sierra.hpp"

  namespace percept {

    /**
     *  Wedge Element with 6 Nodes
     *
     *                                   PARENT Linear 6-Node Wedge Nodes
     *                       5           
     *                     . o
     *                    . / \
     *                   . /   \         Face_Quad_4_3D()  0-1-4-3
     *                  . /     \        Face_Quad_4_3D()  1-4-5-2
     *                 . /       \       Face_Quad_4_3D()  0-2-5-3
     *              2 . o---------o 4    Face_Tri_3_3D()   0-1-2
     *               o . 3      .        Face_Tri_3_3D()   3-4-5
     *              /.\        .
     *             /.  \      .
     *            /.    \    .
     *           /.      \  .
     *          o---------o
     *         0           1
     *
     *---
     *                                   after refinement
     *
     *
     *  Centroid: 20
     *
     *                        5
     *                      . o
     *                     . / \
     *                    . /   \
     *                   . /  19 \
     *                  . o-._ _.-o 13
     *   face: 17      . /14  *    \              {0, 6, 18, 8,  9, 15, 20, 17 }
     *   0253  \   11 * /     |     \             {1, 7, 18, 6,  10, 16, 20, 15 }
     *          \    . /      |12    \            {2, 8, 18, 7,  11, 17, 20, 16 }
     *           \  . o-------o-------o           {3, 14, 19, 12, 9, 17, 20, 15 }
     *             . . 3            .  4          {4, 12, 19, 13, 10, 15, 20, 16 }
     *          2 . * 9        16  .              {5, 13, 19, 14, 11, 16, 20, 17 }
     *           o .     X 20  /  .
     *          /.\             * 10
     *         /.  \    +15    .
     *        /. 18 \         .
     *     8 o-._ _.-o 7     .
     *      /.   *    \     .
     *     /.    |     \   .
     *    /.     |      \ .
     *   o-------o-------o
     *  0        6        1
     *
     *
     *
     *  Face #0 0-1-4-3       Face #1 1-4-5-2       Face #2 0-2-5-3
     *
     *  3        12       4   2        11       5   3        14       5
     *   o-------*-------o     o-------*-------o     o-------*-------o
     *   |       |       |     |       |       |     |       |       |
     *   |       |       |     |       |       |     |       |       |
     *   |       |       |     |       |       |     |       |       |
     *   |     15|       |     |     16|       |     |     17|       |
     *  9*-------*-------*10  7*-------*-------*13  9*-------*-------*11
     *   |       |       |     |       |       |     |       |       |
     *   |       |       |     |       |       |     |       |       |
     *   |       |       |     |       |       |     |       |       |
     *   |       |       |     |       |       |     |       |       |
     *   o-------*-------o     o-------*-------o     o-------*-------o
     *  0        6        1   1       10        4   0        8        2
     *
     *
     *      Face #3:               Face #4:
     *
     *           2                      5
     *           o                      o
     *          / \                    / \
     *         /   \                  /   \
     *        /  18 \                /  19 \
     *     8 o-._ _.-o 7         14 o-._ _.-o 13
     *      /    *    \            /    *    \
     *     /     |     \          /     |     \
     *    /      |      \        /      |      \
     *   o-------o-------o      o-------o-------o
     *  0        6        2    3        12       4
     *
     *
     *
     * Centroid node: 20
     *
     *
     * Hex table:
     * unsigned child[6][8] =  {
     *   {0, 6, 18, 8,  9, 15, 20, 17 },
     *   {1, 7, 18, 6,  10, 16, 20, 15 },
     *   {2, 8, 18, 7,  11, 17, 20, 16 },
     *   {3, 14, 19, 12, 9, 17, 20, 15 },
     *   {4, 12, 19, 13, 10, 15, 20, 16 },
     *   {5, 13, 19, 14, 11, 16, 20, 17 }};
     *
     *   Refined Wedge6 Edge node tables:  NOTE: nodes are not in same order as edges (12-14 reversed with 9-11)
     *
     * unsigned edge_fixup_table[14] = {0,0,0,0, 0,0,0,0,0, 12,13,14, 9, 10, 11};
     *
     * |
     * | static const UInt edge_0[] = { 0, 1,  6 };
     * | static const UInt edge_1[] = { 1, 2,  7 };
     * | static const UInt edge_2[] = { 2, 0,  8 };
     * | static const UInt edge_3[] = { 3, 4, 12 };
     * | static const UInt edge_4[] = { 4, 5, 13 };
     * | static const UInt edge_5[] = { 5, 3, 14 };
     * | static const UInt edge_6[] = { 0, 3,  9 };
     * | static const UInt edge_7[] = { 1, 4, 10 };
     * | static const UInt edge_8[] = { 2, 5, 11 };
     *
     */



    template <>
    class UniformRefinerPattern< shards::Wedge<6>, shards::Hexahedron<8>, 6 > : public URP<shards::Wedge<6>, shards::Hexahedron<8>  >
    {
    private:
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > * m_face_breaker;
      UniformRefinerPattern< shards::Triangle<3>, shards::Quadrilateral<4>, 3, Specialization > * m_face_breaker_tri;

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP< shards::Wedge<6>, shards::Hexahedron<8> >(eMesh)
       {
         EXCEPTWATCH;

         m_primaryEntityRank = stk::topology::ELEMENT_RANK;

         setNeededParts(eMesh, block_names, false);

         m_face_breaker_tri = new UniformRefinerPattern< shards::Triangle<3>, shards::Quadrilateral<4>, 3, Specialization > (eMesh, block_names);
         m_face_breaker = new UniformRefinerPattern< shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > (eMesh, block_names);
       }

      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
        if (m_face_breaker_tri) delete m_face_breaker_tri;
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
        bp.resize(3);
        bp[0] = this;
        bp[1] = m_face_breaker;
        bp[2] = m_face_breaker_tri;
      }

      virtual unsigned getNumNewElemPerElem() override { return 6u; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        EXCEPTWATCH;
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        static stk::mesh::EntityId elems[6][8];

        shards::CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);


// new_sub_entity_nodes[i][j]
#define CENTROID_N NN(stk::topology::ELEMENT_RANK, 0)

#if 0
        // FIXME - maybe the computation of node coorinates should go in the calling code?
        double tmp_x[3];

        for (unsigned i_edge = 0; i_edge < 9; i_edge++)
          {
            double * pts[2] = {EDGE_COORD(i_edge, 0), EDGE_COORD(i_edge, 1)};

            double * mp = getCentroid(pts, 2, eMesh.get_spatial_dim(), tmp_x);
            (void)mp;
#if 1
            std::cout << "i_edge= " << i_edge << " pts = \n"

                      << pts[0][0] << " " << pts[0][1] << " " << pts[0][2] << " \n"
                      << pts[1][0] << " " << pts[1][1] << " " << pts[1][2] << " \n"
                      << " edge = "
                      << mp[0] << " " << mp[1] << " " << mp[2] <<  std::endl;
#endif

            //stk::mesh::Entity new_node =eMesh.createOrGetNode(EDGE_N(i_edge), mp);

            //std::cout << "P[" << eMesh.get_rank() << "] tmp srk createOrGetNode 1 edge node = " << EDGE_N(i_edge) << std::endl;
            //stk::mesh::Entity node = createOrGetNode(nodeRegistry, eMesh, EDGE_N(i_edge));
            eMesh.createOrGetNode(EDGE_N(i_edge), mp);
            nodeRegistry.addToExistingParts(element, m_eMesh.edge_rank(), i_edge);
            nodeRegistry.prolongateFields(element, m_eMesh.edge_rank(), i_edge);
          }

        for (unsigned i_face = 0; i_face < 3; i_face++)
          {
            double * pts[4] = {FACE_COORD(i_face, 0), FACE_COORD(i_face, 1), FACE_COORD(i_face, 2), FACE_COORD(i_face, 3)};

            double * mp = getCentroid(pts, 4, eMesh.get_spatial_dim(), tmp_x);
            (void)mp;
#if 1
            std::cout << "i_face= " << i_face << " pts = \n"

                      << pts[0][0] << " " << pts[0][1] << " " << pts[0][2] << " \n"
                      << pts[1][0] << " " << pts[1][1] << " " << pts[1][2] << " \n"
                      << pts[2][0] << " " << pts[2][1] << " " << pts[2][2] << " \n"
                      << pts[3][0] << " " << pts[3][1] << " " << pts[3][2] << " \n"
                      << " quad face = "
                      << mp[0] << " " << mp[1] << " " << mp[2] <<  std::endl;
#endif

            //stk::mesh::Entity new_node =eMesh.createOrGetNode(FACE_N(i_face), mp);

            //std::cout << "P[" << eMesh.get_rank() << "] tmp srk createOrGetNode 1 face node = " << FACE_N(i_face) << std::endl;
            //stk::mesh::Entity node = createOrGetNode(nodeRegistry, eMesh, FACE_N(i_face));
            eMesh.createOrGetNode(FACE_N(i_face), mp);
            nodeRegistry.addToExistingParts(element, m_eMesh.face_rank(), i_face);
            nodeRegistry.prolongateFields(element, m_eMesh.face_rank(), i_face);
          }

        for (unsigned j_face = 0; j_face < 2; j_face++)
          {
            unsigned i_face = j_face + 3u;

            double * pts[3] = {FACE_COORD(i_face, 0), FACE_COORD(i_face, 1), FACE_COORD(i_face, 2)};

            double * mp = getCentroid(pts, 3, eMesh.get_spatial_dim(), tmp_x);
            (void)mp;
#if 1
            std::cout << "i_face= " << i_face << " pts = \n"

                      << pts[0][0] << " " << pts[0][1] << " " << pts[0][2] << " \n"
                      << pts[1][0] << " " << pts[1][1] << " " << pts[1][2] << " \n"
                      << pts[2][0] << " " << pts[2][1] << " " << pts[2][2] << " \n"
                      << " tri face = "
                      << mp[0] << " " << mp[1] << " " << mp[2] <<  std::endl;
#endif

            //stk::mesh::Entity new_node =eMesh.createOrGetNode(FACE_N(i_face), mp);

            //std::cout << "P[" << eMesh.get_rank() << "] tmp srk createOrGetNode 1 face node = " << FACE_N(i_face) << std::endl;
            //stk::mesh::Entity node = createOrGetNode(nodeRegistry, eMesh, FACE_N(i_face));
            eMesh.createOrGetNode(FACE_N(i_face), mp);
            nodeRegistry.addToExistingParts(element, m_eMesh.face_rank(), i_face);
            nodeRegistry.prolongateFields(element, m_eMesh.face_rank(), i_face);
          }

        // centroid
        double * pts[6] = {VERT_COORD(0), VERT_COORD( 1), VERT_COORD( 2), VERT_COORD(3), VERT_COORD(4), VERT_COORD(5) };
        double * mp = getCentroid(pts, 6, eMesh.get_spatial_dim(), tmp_x);
#if 1
            std::cout << "CENTROID pts = \n"

                      << pts[0][0] << " " << pts[0][1] << " " << pts[0][2] << " \n"
                      << pts[1][0] << " " << pts[1][1] << " " << pts[1][2] << " \n"
                      << pts[2][0] << " " << pts[2][1] << " " << pts[2][2] << " \n"
                      << pts[3][0] << " " << pts[3][1] << " " << pts[3][2] << " \n"
                      << pts[4][0] << " " << pts[4][1] << " " << pts[4][2] << " \n"
                      << pts[5][0] << " " << pts[5][1] << " " << pts[5][2] << " \n"
                      << " centroid = "
                      << mp[0] << " " << mp[1] << " " << mp[2] <<  std::endl;
#endif

        eMesh.createOrGetNode(CENTROID_N, mp);

        nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);
        nodeRegistry.addToExistingParts(element, stk::topology::ELEMENT_RANK, 0u);
        nodeRegistry.prolongateFields(element, stk::topology::ELEMENT_RANK, 0u);

#endif



        static const unsigned edge_fixup_table[15] = {0,0,0,0, 0,0, 6,7,8, 12,13,14, 9, 10, 11};

        static const unsigned child[6][8] =  {
          {0, 6, 18, 8,  9, 15, 20, 17 },
          {1, 7, 18, 6,  10, 16, 20, 15 },
          {2, 8, 18, 7,  11, 17, 20, 16 },
          {3, 14, 19, 12, 9, 17, 20, 15 },
          {4, 12, 19, 13, 10, 15, 20, 16 },
          {5, 13, 19, 14, 11, 16, 20, 17 }};
        
        for (unsigned i_child = 0; i_child < 6; i_child++)
          {
            //std::cout << "\nP[" << eMesh.get_rank() << "] tmp srk createOrGetNode i_child= " << i_child << std::endl;
            for (unsigned jnode = 0; jnode < 8; jnode++)
              {
                unsigned kc = 0;
                unsigned jc = child[i_child][jnode];
                if (jc <= 5) 
                  {
                    kc = VERT_N(jc);
                  }
                else if (6 <= jc && jc <= 14)
                  {
                    unsigned jce = edge_fixup_table[jc];
                    kc = EDGE_N(jce - 6);
                  }
                else if (15 <= jc && jc <= 19)
                  {
                    kc = FACE_N(jc - 15);
                  }
                else
                  {
                    kc = CENTROID_N;
                  }
                elems[i_child][jnode] = kc;

                if (0)
                  {
                    stk::mesh::Entity node = createOrGetNode(nodeRegistry, eMesh, kc);
                    VERIFY_OP_ON(m_eMesh.is_valid(node), ==, true, " hmmm");
                    double *coord = static_cast<double*>(stk::mesh::field_data( *eMesh.get_coordinates_field() , node ));
                    std::cout << "P[" << eMesh.get_rank() << "] tmp srk createOrGetNode id= " << kc 
                              << " coord = " << coord[0] << " " << coord[1] << " " << coord[2]
                              << " jc= " << jc << std::endl;
                  }

              }
          }

#undef CENTROID_N


        for (unsigned ielem=0; ielem < 6; ielem++)
          {
            stk::mesh::Entity newElement = *element_pool;

            if (proc_rank_field)
              {
                double *fdata = static_cast<double*>(stk::mesh::field_data( *proc_rank_field , newElement ));
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            change_entity_parts(eMesh, element, newElement);

            {
              if (!elems[ielem][0])
                {
                  std::cout << "P[" << eMesh.get_rank() << " nid = 0 << " << std::endl;
                  throw std::logic_error("nid = 0 " );
                }
            }
            for (unsigned jnode=0; jnode < 8; jnode++)
              {
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
