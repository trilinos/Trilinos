// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_TriangulateTri_hpp
#define percept_TriangulateTri_hpp

#include <percept/PerceptMesh.hpp>
#include <percept/PerceptBoostArray.hpp>
#include <adapt/CompareCoordinates.hpp>

namespace percept {


  class TriangulateTri {

  public:
    typedef std::array<int, 3> tri_tuple_type_local;
    typedef std::array<stk::mesh::EntityId, 3> tri_tuple_type;

      /**
       *
       *   Convention: input is the element's nodes and the marks on the edges.  Output is an array
       *     of "elements" defined as local id's of nodes forming those elements, where {0,1,2} represent
       *     the original vertices and {3,4,5} are the edges:
       *
       *                2
       *                o
       *               / \
       *              /   \
       *             /     \
       *          5 *       * 4
       *           /         \
       *          /           \
       *         /             \
       *        o-------*-------o
       *       0        3        1
       */

      // Note: this will form the basis of triangulating faces in 3D, so it is generalized to a
      //   generic method.
      // Note: code below doesn't orient the face except for a rotation - we need a polarity flip check as
      //   well for the general, 3D face case
      //

#define T_VERT_N(i) (i)
#define T_EDGE_N(i) ((i)+3)

      static void triangulate_face(PerceptMesh& eMesh, stk::mesh::Entity elem_nodes[3], unsigned edge_marks[3],
                                   vector<tri_tuple_type_local>& elems)
      {
        elems.resize(0);

        const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< shards::Triangle<3> >();

        shards::CellTopology cell_topo(cell_topo_data);
        //const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK); /NLM
        stk::mesh::FieldBase* coordField = eMesh.get_coordinates_field();

        std::array<double *,3> node_coords;
        for (int i=0; i<3; ++i)
          {
             node_coords[i] = static_cast<double*>(stk::mesh::field_data( *coordField , elem_nodes[i] ));
          }
        const std::array<int, 3> node_rank = get_rank_of_nodes_based_on_coordinates( node_coords );

        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 3; iedge++)
          {
            unsigned num_nodes_on_edge = edge_marks[iedge];
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }

        //std::cout << "tmp RefinerPattern_Tri3_Tri3_N::num_edges_marked= " << num_edges_marked << std::endl;
        if (num_edges_marked == 3)
          {
            elems.resize(4);

            elems[0] = { T_VERT_N(0),    T_EDGE_N(0), T_EDGE_N(2) };
            elems[1] = { T_VERT_N(1),    T_EDGE_N(1), T_EDGE_N(0) };
            elems[2] = { T_VERT_N(2),    T_EDGE_N(2), T_EDGE_N(1) };
            elems[3] = { T_EDGE_N(0),    T_EDGE_N(1), T_EDGE_N(2) };
          }
        else if (num_edges_marked == 2)
          {
            /**
             *
             *   case 1: jedge == max length edge
             *
             *                i2
             *                o
             *               /|\
             *              / | \
             *             /  |  \
             *            /   |   * jedgep
             *           /    |  / \
             *          /     | /   \
             *         /      |/     \
             *        o-------*-------o
             *       i0      jedge     i1
             *
             *
             *   case 2: jedge+1 == max length edge
             *
             *                i2
             *                o
             *               / \
             *              /   \
             *             /     \
             *            /     _.* jedgep
             *           /   _.* / \
             *          / _.*   /   \
             *         /.*     /     \
             *        o-------*-------o
             *       i0      jedge     i1
             *
             */

            elems.resize(3);

            // find first of two marked edges in sequence (get in "standard" orientation), and longest marked edge
            int jedge = -1;
            int jedge_max_edge = -1;
            double max_edge_length = -1.0;
            unsigned id_diff_0 = 0u;
            unsigned id_diff_1 = 0u;
            for (int iedge = 0; iedge < 3; iedge++)
              {

                unsigned num_nodes_on_edge = edge_marks[iedge];
                unsigned num_nodes_on_edge_p = edge_marks[(iedge+1)%3];
                if (num_nodes_on_edge && num_nodes_on_edge_p)
                  {
                    jedge = iedge;
                  }

                if (num_nodes_on_edge)
                  {
                    stk::mesh::Entity node_0 = elem_nodes[cell_topo_data->edge[iedge].node[0]];
                    stk::mesh::Entity node_1 = elem_nodes[cell_topo_data->edge[iedge].node[1]];
                    int rank_node_0 = node_rank[cell_topo_data->edge[iedge].node[0]];
                    int rank_node_1 = node_rank[cell_topo_data->edge[iedge].node[1]];

                    //bool reverse = false;
                    // ensure edge_len is computed identically, independent of edge orientation
                    if (node_rank[cell_topo_data->edge[iedge].node[0]] > node_rank[cell_topo_data->edge[iedge].node[1]])
                      {
                        //reverse = true;
                        stk::mesh::Entity node_temp = node_0;
                        node_0 = node_1;
                        node_1 = node_temp;
                        const int rank_temp = rank_node_0;
                        rank_node_0 = rank_node_1;
                        rank_node_1 = rank_temp;
                      }

                    double * const coord_0 = static_cast<double*>(stk::mesh::field_data( *coordField , node_0 ));
                    double * const coord_1 = static_cast<double*>(stk::mesh::field_data( *coordField , node_1 ));
                    double edge_len_squared = 0.0;

                    edge_len_squared =
                      (coord_0[0] - coord_1[0])*(coord_0[0] - coord_1[0])+
                      (coord_0[1] - coord_1[1])*(coord_0[1] - coord_1[1])+
                      (eMesh.get_spatial_dim() == 2 ? 0 :
                       (coord_0[2] - coord_1[2])*(coord_0[2] - coord_1[2]) );

                    if (edge_len_squared > max_edge_length)
                      {
                        id_diff_0 = rank_node_0;
                        id_diff_1 = rank_node_1;
                        max_edge_length = edge_len_squared;
                        jedge_max_edge = iedge;
                      }
                    // intentional floating-point comparison (tie-break)
                    else if (edge_len_squared == max_edge_length)
                      {
                        unsigned loc_id_diff_0 = rank_node_0;
                        unsigned loc_id_diff_1 = rank_node_1;
                        bool lexical_less = false;
                        if (loc_id_diff_0 < id_diff_0)
                          {
                            lexical_less = true;
                          }
                        else if (loc_id_diff_0 == id_diff_0 && loc_id_diff_1 < id_diff_1)
                          {
                            lexical_less = true;
                          }
                        if (!lexical_less)
                          {
                            max_edge_length = edge_len_squared;
                            jedge_max_edge = iedge;
                          }
                      }
                  }
              }

            if (jedge < 0 || jedge_max_edge < 0)
              {
                std::cout << "jedge = " << jedge << " jedge_max_edge = " << jedge_max_edge << std::endl;
                throw std::runtime_error("RefinerPattern_Tri3_Tri3_N jedge < 0");
              }

            //stk::mesh::Entity node0 = *elem_nodes[iii].entity();
            int i0 = cell_topo_data->edge[jedge].node[0];
            if (i0 != jedge)
              {
                std::cout << "i0 = " << i0 << " jedge= " << jedge << std::endl;
                throw std::runtime_error("RefinerPattern_Tri3_Tri3_N i0 != jedge");
              }

            int i1 = (i0+1)%3;
            int i2 = (i0+2)%3;
            int jedgep = (jedge+1)%3;
            if (jedge_max_edge == jedge)
              {
                elems[0] = { T_VERT_N(i0),    T_EDGE_N(jedge),         T_VERT_N(i2)       };
                elems[1] = { T_EDGE_N(jedge), T_EDGE_N( jedgep ),      T_VERT_N(i2)       };
                elems[2] = { T_EDGE_N(jedge), T_VERT_N(i1),            T_EDGE_N( jedgep ) };
              }
            else
              {
                elems[0] = { T_VERT_N(i0),    T_EDGE_N(jedge),         T_EDGE_N( jedgep ) };
                elems[1] = { T_VERT_N(i0),    T_EDGE_N( jedgep ),      T_VERT_N(i2)       };
                elems[2] = { T_EDGE_N(jedge), T_VERT_N(i1),            T_EDGE_N( jedgep ) };
              }
          }
        else if (num_edges_marked == 1)
          {
            elems.resize(2);
            for (int iedge = 0; iedge < 3; iedge++)
              {
                unsigned num_nodes_on_edge = edge_marks[iedge];
                if (num_nodes_on_edge)
                  {
                    elems[0] = {T_VERT_N(iedge), T_EDGE_N(iedge), T_VERT_N((iedge+2)%3) };
                    elems[1] = {T_EDGE_N(iedge), T_VERT_N((iedge+1)%3), T_VERT_N((iedge+2)%3) };
                    break;
                  }
              }
          }
        else if (num_edges_marked == 0)
          {
#if 1
            if (allow_single_refine)
              {
                // this allows each level to be at the same hierarchical level by having a single parent to single child
                elems.resize(1);
                elems[0] = {T_VERT_N(0), T_VERT_N(1), T_VERT_N(2) };
              }
#else
            if (elems.size() != 0)
              {
                std::cout << "tmp num_edges_marked= 0 " << elems.size() << std::endl;
                throw std::logic_error("hmmmmmmmmmmm");
              }

            return;
#endif
          }

      }
  };

#undef T_VERT_N
#undef T_EDGE_N
}

#endif
