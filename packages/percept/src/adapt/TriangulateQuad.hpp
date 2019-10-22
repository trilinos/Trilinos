// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_TriangulateQuad_hpp
#define percept_TriangulateQuad_hpp

#include <percept/PerceptMesh.hpp>
#include <percept/PerceptBoostArray.hpp>

namespace percept {


  typedef std::array<unsigned, 3> quad_to_tri_tuple_type_local;
  typedef std::array<stk::mesh::EntityId, 3 > quad_to_tri_tuple_type;

  typedef std::array<unsigned, 4> quad_to_quad_tuple_type_local;
  typedef std::array<stk::mesh::EntityId, 4> quad_to_quad_tuple_type;

  /**
   *
   *   Convention: input is the element's nodes and the marks on the edges.  Output is an array
   *     of "elements" defined as local id's of nodes forming those elements, where {0,1,2,3} represent
   *     the original vertices and {4,5,6,7} are the edges, and {8} is the centroid:
   *
   *        3       6
   *        o-------*-------o 2
   *        |               |
   *        |               |
   *        |        8      |
   *     7  *       *       * 5
   *        |               |
   *        |               |
   *        |               |
   *        o-------*-------o
   *       0        4        1
   */


#define Q_VERT_N(i) (i)
#define Q_EDGE_N(i) ((i)+4)
#define Q_CENTROID_N (8)

  class TriangulateQuad {
    bool m_use_only_tris;
    bool m_avoid_centroid_node;
  public:

    TriangulateQuad(bool use_only_tris=true, bool avoid_centroid_node=false) : m_use_only_tris(use_only_tris), m_avoid_centroid_node(avoid_centroid_node) {}

    void triangulate_quad_face(unsigned edge_marks[4],
                               vector<quad_to_tri_tuple_type_local>& elems,
                               vector<quad_to_quad_tuple_type_local>& elems_quad
                               )
    {
      elems.resize(0);
      elems_quad.resize(0);

      const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< shards::Quadrilateral<4> >();

      shards::CellTopology cell_topo(cell_topo_data);

      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < 4; iedge++)
        {
          unsigned num_nodes_on_edge = edge_marks[iedge];
          if (num_nodes_on_edge)
            {
              ++num_edges_marked;
            }
        }

      //std::cout << "tmp RefinerPattern_Tri3_Tri3_N::num_edges_marked= " << num_edges_marked << std::endl;
      if (num_edges_marked == 0)
        {
          return;
        }

      if (num_edges_marked == 4)
        {
          /**
           *
           *
           *                centroid
           *   je+3 o-------o
           *        |       |
           *        |       |
           *        |       |
           *        |       |
           *        o-------*-------o
           *       i0   je = i0+4    i1
           */

          for (unsigned iedge = 0; iedge < 4; ++iedge)
            {
              unsigned i0 = iedge;
              //unsigned i1 = (iedge+1) % 4;
              unsigned jedge = iedge + 4;
              unsigned kedge = (iedge + 3) % 4 + 4;
              elems_quad.push_back({i0, jedge, Q_CENTROID_N, kedge});
            }
          return;
        }

      if (m_use_only_tris)
        {
          /**
           *
           *   case 1: edge is marked
           *
           *                centroid
           *                o
           *               /|\
           *              / | \
           *             /  |  \
           *            /   |   \
           *           /    |    \
           *          /     |     \
           *         /      |      \
           *        o-------*-------o
           *       i0      i0+4     i1
           *
           *
           *   case 2: edge is not marked
           *
           *                centroid
           *                o
           *               / \
           *              /   \
           *             /     \
           *            /       \
           *           /         \
           *          /           \
           *         /             \
           *        o---------------o
           *       i0      i0+4     i1
           *
           */

          elems.resize(0);
          for (unsigned iedge = 0; iedge < 4; ++iedge)
            {
              unsigned i0 = iedge;
              unsigned i1 = (iedge+1) % 4;
              unsigned jedge = iedge + 4;
              if (edge_marks[iedge])
                {
                  elems.push_back({i0, jedge, Q_CENTROID_N});
                  elems.push_back({jedge, i1, Q_CENTROID_N});
                }
              else
                {
                  elems.push_back({i0, i1, Q_CENTROID_N});
                }
            }
        }
      else
        {
          // else if (num_edges_marked == 4)
          //   {
          //     // error?
          //   }
          if (num_edges_marked == 1)
            {
              /**
               *       i3
               *        o---------------o i2
               *        |\             /|
               *        | \           / |
               *        |  \         /  |
               *        |   \       /   |
               *        |    \     /    |
               *        |     \   /     |
               *        |      \ /      |
               *        o-------*-------o
               *       i0      i0+4     i1
               */
              unsigned iedge = 0;
              for (unsigned jedge = 0; jedge < 4; jedge++)
                {
                  if (edge_marks[jedge])
                    {
                      iedge = jedge;
                      break;
                    }
                }
              unsigned i0 = iedge;
              unsigned i1 = (iedge+1) % 4;
              unsigned i2 = (iedge+2) % 4;
              unsigned i3 = (iedge+3) % 4;
              unsigned jedge = iedge + 4;
              elems.push_back({i0, jedge, i3});
              elems.push_back({jedge, i2, i3});
              elems.push_back({jedge, i1, i2});
            }
          else if (num_edges_marked == 2)
            {

              unsigned iedge0 = 0, iedge1 = 0;
              bool adj = false;
              for (unsigned jedge = 0; jedge < 4; jedge++)
                {
                  if (edge_marks[jedge] && edge_marks[(jedge+1) % 4])
                    {
                      adj = true;
                      iedge0 = jedge;
                      iedge1 = (jedge+1) % 4;
                      break;
                    }
                  if (edge_marks[jedge] && edge_marks[(jedge+2) % 4])
                    {
                      adj = false;
                      iedge0 = jedge;
                      iedge1 = (jedge+2) % 4;
                      break;
                    }
                }

              unsigned i0 = iedge0;
              unsigned i1 = (iedge0+1) % 4;
              unsigned i2 = (iedge0+2) % 4;
              unsigned i3 = (iedge0+3) % 4;
              unsigned jedge0 = iedge0 + 4;
              unsigned jedge1 = iedge1 + 4;

              if (adj)
                {
                  /**
                   *       i3
                   *        o---------------o i2
                   *        |\ \            |
                   *        | \   \         |
                   *        |  \     \      |
                   *        |   \        \  |
                   *        |    \       / -* i1+4  jedge1
                   *        |     \    /    |
                   *        |      \ /      |
                   *        o-------*-------o
                   *       i0      i0+4     i1
                   *             jedge0
                   */

                  if (m_avoid_centroid_node)
                    {
                      elems.push_back({i0, jedge0, i3});
                      elems.push_back({jedge0, i1, jedge1});
                      elems.push_back({jedge1, i2, i3});
                      elems.push_back({jedge0, jedge1, i3});
                    }
                  else
                    {
                      /**
                       *       i3
                       *        o---------------o i2
                       *        |\              |
                       *        |  \            |
                       *        |    \          |
                       *        |      \        |
                       *        |       *-------* i1+4
                       *        |       |       |
                       *        |       |       |
                       *        o-------*-------o
                       *       i0      i0+4     i1
                       */
                      elems_quad.push_back({i0, jedge0, Q_CENTROID_N, i3});
                      elems_quad.push_back({jedge0, i1, jedge1, Q_CENTROID_N});
                      elems_quad.push_back({jedge1, i2, i3, Q_CENTROID_N});
                    }
                }
              else
                {
                  /**
                   *
                   *       i3       i2+4
                   *        o-------*-------o i2
                   *        |       |       |
                   *        |       |       |
                   *        |       |       |
                   *        |       |       |
                   *        |       |       |
                   *        |       |       |
                   *        |       |       |
                   *        o-------*-------o
                   *       i0      i0+4     i1
                   */

                  elems_quad.push_back({i0, jedge0, jedge1, i3});
                  elems_quad.push_back({jedge0, i1, i2, jedge1});
                }
            }
          else if (num_edges_marked == 3)
            {
              /**
               *
               *               jedge2
               *       i3       i2+4
               *        o-------*-------o i2
               *        |       | \     |
               *        |       |   \   |
               *        |       |     \ |
               *        |       |     / * i1+4, jedge1
               *        |       |    /  |
               *        |       |  /    |
               *        |       |/      |
               *        o-------*-------o
               *       i0      i0+4     i1
               *               jedge0
               */
              unsigned iedge0 = 0, iedge1 = 0, iedge2 = 0;
              for (unsigned jedge = 0; jedge < 4; jedge++)
                {
                  if (edge_marks[jedge] && edge_marks[(jedge+1) % 4] && edge_marks[(jedge+2) % 4])
                    {
                      iedge0 = jedge;
                      iedge1 = (jedge+1) % 4;
                      iedge2 = (jedge+2) % 4;
                      break;
                    }
                }

              unsigned i0 = iedge0;
              unsigned i1 = (iedge0+1) % 4;
              unsigned i2 = (iedge0+2) % 4;
              unsigned i3 = (iedge0+3) % 4;
              unsigned jedge0 = iedge0 + 4;
              unsigned jedge1 = iedge1 + 4;
              unsigned jedge2 = iedge2 + 4;

              elems_quad.push_back({i0, jedge0, jedge2, i3});
              elems.push_back({jedge0, i1, jedge1});
              elems.push_back({jedge0, jedge1, jedge2});
              elems.push_back({jedge1, i2, jedge2});
            }
          else if (num_edges_marked == 4)
            {
              // this shouldn't occur here, see above, just keeping it here for history and possible 
              //   future use
              bool allow_4_edge_marks = true;
              if (!allow_4_edge_marks) return;
              /**
               *
               *               jedge2
               *       i3       i2+4
               *        o-------*-------o i2
               *        |     / | \     |
               *        |   /   |   \   |
               *        | /     |     \ |
               *  i3+4 *|       |     / * i1+4, jedge1
               * jedge3 | \     |    /  |
               *        |  \    |  /    |
               *        |    \  |/      |
               *        o-------*-------o
               *       i0      i0+4     i1
               *               jedge0
               */
              unsigned iedge0 = 0, iedge1 = 1, iedge2 = 2, iedge3 = 3;

              unsigned i0 = 0;
              unsigned i1 = 1;
              unsigned i2 = 2;
              unsigned i3 = 3;
              unsigned jedge0 = iedge0 + 4;
              unsigned jedge1 = iedge1 + 4;
              unsigned jedge2 = iedge2 + 4;
              unsigned jedge3 = iedge3 + 4;

              elems.push_back({jedge0, i1, jedge1});
              elems.push_back({jedge0, jedge1, jedge2});
              elems.push_back({jedge1, i2, jedge2});

              elems.push_back({jedge2, i3, jedge3});
              elems.push_back({jedge0, jedge2, jedge3});
              elems.push_back({i0, jedge0, jedge3});
            }

        }

    }
  };

#undef Q_VERT_N
#undef Q_EDGE_N
#undef Q_CENTROID_N

}

#endif
