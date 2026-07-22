// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_DiscretizeWedge_hpp
#define percept_DiscretizeWedge_hpp

#include <percept/PerceptMesh.hpp>
#include <percept/PerceptBoostArray.hpp>
#include <adapt/TriangulateQuad.hpp>
#include <adapt/TriangulateTri.hpp>

namespace percept {

  typedef std::array<unsigned, 4> wedge_to_tet_tuple_type_local;
  typedef std::array<unsigned, 5> wedge_to_pyr_tuple_type_local;
  typedef std::array<unsigned, 6> wedge_to_wedge_tuple_type_local;
  typedef std::array<stk::mesh::EntityId, 4> wedge_to_tet_tuple_type;
  typedef std::array<stk::mesh::EntityId, 5> wedge_to_pyr_tuple_type;
  typedef std::array<stk::mesh::EntityId, 6> wedge_to_wedge_tuple_type;


  /**
   *
   *
   *   Convention: input is the element's nodes and the marks on the edges.  Output is an array
   *     of "elements" defined as local id's of nodes forming those elements, where {0,..,5} represent
   *     the original vertices and {6..14} are the edges, {15..17} are quadratic face centroids
   *     and {18} is the centroid used to discretize the wedge
   *     (these are numberings used by Shards for the extra/quadratic nodes, and so we map these
   *      numbers back to the local edge/face/centroid using tables; for example, quadratic node 12
   *      corresponds to edge {3,4} which is not the 7th edge (12-6+1  = 7), but is the 4th in
   *      Shards' convention - see table "edges" and edge_map_rev below)
   *
   *
   *   typedef
   *     MakeTypeList< IndexList< 0 , 1 ,   6 > ,
   *                   IndexList< 1 , 2 ,   7 > ,
   *                   IndexList< 2 , 0 ,   8 > ,
   *                   IndexList< 3 , 4 ,  12 > ,
   *                   IndexList< 4 , 5 ,  13 > ,
   *                   IndexList< 5 , 3 ,  14 > ,
   *                   IndexList< 0 , 3 ,   9 > ,
   *                   IndexList< 1 , 4 ,  10 > ,
   *                   IndexList< 2 , 5 ,  11 >
   *     >::type WedgeEdgeNodeMap ;
   *
   *   typedef
   *     MakeTypeList< IndexList< 0 , 1 , 4 , 3 ,   6 , 10 , 12 ,  9 ,  15 > ,
   *                   IndexList< 1 , 2 , 5 , 4 ,   7 , 11 , 13 , 10 ,  16 > ,
   *                   IndexList< 0 , 3 , 5 , 2 ,   9 , 14 , 11 ,  8 ,  17 > ,
   *                   IndexList< 0 , 2 , 1 ,       8 ,  7 ,  6 > ,
   *                   IndexList< 3 , 4 , 5 ,      12 , 13 , 14 >
   *     >::type WedgeFaceNodeMap ;
   *
   *
   *  Wedge Element with 6 Nodes
   *  Pentahedral Element with 6 Nodes
   *  3-Sided Prism Element with 6 Nodes
   *
   *                                   PARENT Linear 6-Node Wedge Nodes
   *                       2           (SPACE_DIM = 3!)
   *                     . o
   *                    . / \
   *            /      . /   \         Face_Quad_4_3D()  0-1-4-3
   *           /      . /     \        Face_Quad_4_3D()  1-2-5-4 (error in orig code - 1-4-5-2)
   *          /      . /       \       Face_Quad_4_3D()  0-3-5-2 (error in orig code - 0-2-5-3)
   *         /    5 . o---------o 1    Face_Tri_3_3D()   0-2-1 (error was 0-1-2)
   *        /      o . 0      .        Face_Tri_3_3D()   3-4-5
   *       /      /.\        .
   *      /      /.  \      .
   *    |/_     /.    \    .
   *           /.      \  .
   *    Z     o---------o
   *         3           4
   *
   * -------------------------------------------------------------------------------
   *
   *                                   PARENT Linear 6-Node Wedge Nodes
   *                       2           (SPACE_DIM = 3!)
   *                     . o
   *                    . / \
   *            /      . 2   1         Face_Quad_4_3D()  0-1-4-3 edges = {0,7,3,6}
   *           /      8 /     \        Face_Quad_4_3D()  1-2-5-4 edges = {1,8,4,7}
   *          /      . /       \       Face_Quad_4_3D()  0-3-5-2 edges = {6,5,8,2}
   *         /    5 . o----0----o 1    Face_Tri_3_3D()   0-2-1   edges = {0,2,1}
   *        /      o . 0      .        Face_Tri_3_3D()   3-4-5   edges = {3,4,5}
   *       /      /6\        .
   *      /      5.  4      7
   *    |/_     /.    \    .
   *           /.      \  .
   *    Z     o----3----o
   *         3           4
   *
   * -------------------------------------------------------------------------------
   *                          5
   *                          o
   *                         / \
   *                     14 /   \
   *                       o     o 13
   *                      /       \
   *            Z      3 /         \
   *          _         o-----o-----o
   *           /|            12      4
   *          /       11
   *         /         o
   *        /         . .
   *       /         .   .
   *      /         .     .
   *     /         .       .
   *              .         .
   *           9 o...........o 10
   *
   *             2
   *             o
   *            / \
   *           /   \
   *        8 o     o 7
   *         /       \
   *        /         \
   *       o-----o-----o
   *      0      6      1
   *
   *
   *   4-wedge refinement: (only 1st 6 edges marked)
   *     {{0, 6, 8, 3, 12, 14},
   *      {6, 1, 7, 12, 4, 13},
   *      {6, 7, 8, 12, 13, 14},
   *      {8, 7, 2, 14, 13, 5}
   *     }
   *   2-wedge refinement: (only 2 edges marked (6,12), (7,13) or (8,14)
   *     {{{0, 6, 2, 3, 12, 5},
   *      {6, 1, 2, 12, 4, 5}},
   *     {{1, 7, 0, 4, 13, 3},
   *      {7, 2, 0, 13, 5, 3}},
   *     {{2, 8, 1, 5, 14, 4},
   *      {8, 0, 1, 14, 3, 4}
   *     }
   *     }
   *
   *   CHILD Wedge6 node tables;
   * |
   * | static const UInt child_0[] = {  0, 6, 8,   9, 15, 17 };
   * | static const UInt child_1[] = {  6, 1, 7,  15, 10, 16 };
   * | static const UInt child_2[] = {  8, 7, 2,  17, 16, 11 };
   * | static const UInt child_3[] = {  7, 8, 6,  16, 17, 15 };
   * | static const UInt child_4[] = {  9, 15, 17,   3, 12, 14 };
   * | static const UInt child_5[] = { 15, 10, 16,  12, 4, 13 };
   * | static const UInt child_6[] = { 17, 16, 11,  14, 13, 5 };
   * | static const UInt child_7[] = { 16, 17, 15,  13, 14, 12 };
   *
   *   Refined Wedge6 Edge node tables:
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
   *   Refined Wedge6 Face node tables:
   * |
   * | static const UInt face_0[] =  { 0, 1, 4, 3,  6, 10, 12, 9, 15    };
   * | static const UInt face_1[] =  { 1, 2, 5, 4,  7, 11, 13, 10, 16    };
   * | static const UInt face_2[] =  { 0, 3, 5, 2,  9, 14, 11, 8, 17    };
   * | static const UInt face_3[] =  { 0, 2, 1,    8, 7, 6, -1, -1, -1};
   * | static const UInt face_4[] =  { 3, 4, 5,   12, 13, 14, -1, -1, -1};
   *
   *
   *
   *                          5
   *                          o
   *                         / \
   *                     14 /   \
   *                       o     o 13
   *          Z           /       \
   *        _          3 /         \
   *         /|         o-----o-----o
   *        /                12      4
   *       /          11                        errors here... FIXME
   *      /            o              Face_Quad_8_3D() 0-6-1-10-4-12-3-9
   *     /            . .             Face_Quad_8_3D() 1-10-4-13-5-11-2-7
   *    /         17 .   . 16         Face_Quad_8_3D() 0-8-2-11-5-14-3-9
   *                o     o           Face_Tri_6_3D()  0-6-1-7-2-8
   *               .       .          Face_Tri_6_3D()  3-12-4-13-5-14
   *              .         .
   *           9 o.....o.....o 10
   *                  15
   *             2
   *             o
   *            / \
   *           /   \
   *        8 o     o 7
   *         /       \
   *        /         \
   *       o-----o-----o
   *      0      6      1
   *

   *
   *
   */

  class DiscretizeWedge {
  public:
    static const bool s_debug = false;

    static const int nedges = 9;
    static const int nfaces = 5;
    static const int nquad_faces = 3;  // these come first
    static const int ntri_faces = 2;

    static const int edge_offset = 6;
    static const int face_offset = 15;

    static const int centroid = 18;
    static const int mapped_centroid = 18;

    static bool is_v(int i) { return i < edge_offset; }
    static bool is_e(int i) { return (edge_offset <= i && i < edge_offset + nedges); }
    static bool is_f(int i) { return (face_offset <= i && i < face_offset + nquad_faces); }
    static bool is_c(int i) { return i == centroid; }

    static bool is_mapped_v(int i) { return is_v(i); }
    static bool is_mapped_e(int i) { return is_e(i); }
    static bool is_mapped_f(int i) { return is_f(i); }
    static bool is_mapped_c(int i) { return i == mapped_centroid; }

    static int edges_of_face(int face, int edge)
    {
      static const int t_edges_of_face[5][4] = {
        {0,7,3,6},
        {1,8,4,7},
        {6,5,8,2},
        {0,2,1,-1},
        {3,4,5,-1}
      };
      return t_edges_of_face[face][edge];
    }

    static int map_to_wedge(int hp, const int *edge_map_rev, const int *face_map_rev)
    {
      int rv = hp;
      if (is_v(hp))
        return hp;
      if (is_e(hp))
        rv = edge_offset + edge_map_rev[ hp - edge_offset ];
      else if (is_f(hp))
        rv = face_offset + face_map_rev[ hp - face_offset ];
      else if (is_c(hp))
        rv = mapped_centroid;
      else
        throw std::runtime_error("bad map");
      return rv;
    }

    static void discretize(PerceptMesh& eMesh, stk::mesh::Entity wedge,
                           unsigned edge_marks[nedges], unsigned face_marks[nfaces],
                           vector<wedge_to_tet_tuple_type_local>& elems_tet,
                           vector<wedge_to_pyr_tuple_type_local>& elems_pyr,
                           vector<wedge_to_wedge_tuple_type_local>& elems_wedge,
                           bool allow_special_refine = false,
                           bool m_debug = false
                           )
    {
      const percept::MyPairIterRelation elem_nodes (eMesh, wedge, stk::topology::NODE_RANK);

      // Shards edges {node0, node1, quadratic_node}
      static const int edges[][3] = {
        { 0, 1,  6 },
        { 1, 2,  7 },
        { 2, 0,  8 },
        { 3, 4, 12 },
        { 4, 5, 13 },
        { 5, 3, 14 },
        { 0, 3,  9 },
        { 1, 4, 10 },
        { 2, 5, 11 }
      };

      static const int edge_map[] =     {6, 7, 8,    12, 13, 14,   9, 10, 11};
      // map of quadratic node to index in edges table
      static const int edge_map_rev[] = {0, 1,  2,  6, 7, 8,  3, 4, 5}; // edge_map_rev[edge_map[i] - 6] = i;

      // similar to edges: {node0...node3, edge0(quadraticNode),...edge4, faceQuadraticNode}
      static const int faces[][9] =  {
        { 0, 1, 4, 3,  6, 10, 12, 9, 15    },
        { 1, 2, 5, 4,  7, 11, 13, 10, 16    },
        { 0, 3, 5, 2,  9, 14, 11, 8, 17    },
        { 0, 2, 1,    8, 7, 6, -1, -1, -1},
        { 3, 4, 5,   12, 13, 14, -1, -1, -1}
      };
      static const int tri_faces_not_outward_normal[][6] =  {
        { 0, 1, 2, 6, 7, 8},
        { 3, 4, 5, 12, 13, 14}
      };
      static const int face_map[] =     { 15, 16, 17, -1, -1 };
      static const int face_map_rev[] = {  0,  1,  2,  -1, -1 };   // face_map_rev[face_map[i] - 15] = i;

      (void)edge_map;
      (void)edges;
      (void)face_map;

      elems_tet.resize(0);
      elems_pyr.resize(0);
      elems_wedge.resize(0);

      const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< shards::Wedge<6> >();
      shards::CellTopology cell_topo(cell_topo_data);

      bool special_iso = true;
      int special_te = -1;
      int special_te_case_2 = -1;
      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < nedges; iedge++)
        {
          unsigned num_nodes_on_edge = edge_marks[iedge];
          if (iedge < 6 && !num_nodes_on_edge) special_iso=false;
          if (iedge >= 6 && num_nodes_on_edge) special_iso=false;
          if (num_nodes_on_edge)
            {
              ++num_edges_marked;
            }
        }
      unsigned num_faces_marked=0;
      for (int iface = 0; iface < 3; iface++)
        {
          unsigned num_nodes_on_face = face_marks[iface];
          if (num_nodes_on_face)
            {
              if (iface < 3) special_iso = false;
              ++num_faces_marked;
            }
        }
      if (num_edges_marked == 2) // && num_faces_marked == 0)
        {
          special_te = -1;
          for (unsigned ii=0; ii < 3; ii++)
            {
              if (edge_marks[ii] && edge_marks[ii+3])
                {
                  special_te = ii;
                  break;
                }
            }
        }
      if (num_edges_marked == 4) // && num_faces_marked == 0)
        {
          special_te_case_2 = -1;
          int nfound = 0;
          for (unsigned ii=0; ii < 3; ii++)
            {
              if (edge_marks[ii] && edge_marks[ii+3])
                {
                  ++nfound;
                  if (special_te_case_2 < 0)
                    special_te_case_2 = 0;
                  special_te_case_2 += 1 << (ii);
                }
            }
          if (nfound != 2)
            special_te_case_2 = -1;
        }

      //VERIFY_OP_ON(special_te_case_2, == , -1, "bad");
      if (m_debug) //  || special_te_case_2 >= 0)
        {
          std::vector<unsigned> em(&edge_marks[0], &edge_marks[0] + nedges);
          std::vector<unsigned> fm(&face_marks[0], &face_marks[0] + nfaces);
          std::cout << "P[" << eMesh.get_rank() << "] tmp DiscretizeWedge::num_edges_marked for wedge[ " << eMesh.identifier(wedge)
                    << "] = " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                    << " allow_special_refine= " << allow_special_refine
                    << " special_iso= " << special_iso
                    << " special_te= " << special_te
                    << " special_te_case_2= " << special_te_case_2
                    << " em= " << em
                    << "\n fm= " << fm
                    << std::endl;
        }

      if (num_edges_marked == 0 && num_faces_marked == 0)
        {
          return;
        }

      if ((m_debug) && allow_special_refine && (special_te_case_2 || special_te || special_iso))
        {
          std::vector<unsigned> em(&edge_marks[0], &edge_marks[0] + nedges);
          std::cout << "P[" << eMesh.get_rank() << "] tmp DiscretizeWedge::num_edges_marked for wedge[ " << eMesh.identifier(wedge)
                    << "] = " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                    << " allow_special_refine= " << allow_special_refine
                    << " special_iso= " << special_iso
                    << " special_te= " << special_te
                    << " special_te_case_2= " << special_te_case_2
                    << " em= " << em
                    << std::endl;
        }

      if (allow_special_refine && special_iso)
        {
          //4-wedge refinement: (only 1st 6 edges marked)
          static const unsigned T[][6] =
            {{0, 6, 8, 3, 12, 14},
             {6, 1, 7, 12, 4, 13},
             {6, 7, 8, 12, 13, 14},
             {8, 7, 2, 14, 13, 5}
            };
          wedge_to_wedge_tuple_type_local hw;
          for (unsigned ii=0; ii < 4; ++ii)
            {
              for (unsigned j=0; j < 6; ++j)
                {
                  hw[j] = map_to_wedge(T[ii][j], edge_map_rev, face_map_rev);
                }
              elems_wedge.push_back(hw);
            }
          return;
        }

      if (allow_special_refine && special_te >= 0)
        {
          //2-wedge refinement: (only 2 edges marked (6,12), (7,13) or (8,14)
          static const unsigned T[][2][6] = {
            {{0, 6, 2, 3, 12, 5},
             {6, 1, 2, 12, 4, 5}},
            {{1, 7, 0, 4, 13, 3},
             {7, 2, 0, 13, 5, 3}},
            {{2, 8, 1, 5, 14, 4},
             {8, 0, 1, 14, 3, 4}
            }
          };
          int J = special_te;
          wedge_to_wedge_tuple_type_local hw;
          for (unsigned ii=0; ii < 2; ++ii)
            {
              for (unsigned j=0; j < 6; ++j)
                {
                  hw[j] = map_to_wedge(T[J][ii][j], edge_map_rev, face_map_rev);
                }
              elems_wedge.push_back(hw);
            }
          return;
        }

      if (allow_special_refine && special_te_case_2 >= 0)
        {
          // 4-wedge refinement where two edges are marked (an their vertical (offset) edge is also marked)
          std::vector<tri_tuple_type_local> elems_tri[2];
          unsigned tri_local_nodes[2][6];
          for (unsigned iface=3; iface < 5; ++iface)
            {
              unsigned tri_edge_marks[3] = {0};
              //unsigned num_tri_edge_marks = 0;
              //stk::mesh::Entity tri_local_entities[3] = {stk::mesh::Entity()};
              stk::mesh::Entity tri_local_entities_non_outward_normal[3] = {stk::mesh::Entity()};
              for (unsigned j = 0; j < 6; j++)
                {
                  tri_local_nodes[iface-3][j] = tri_faces_not_outward_normal[iface-3][j];
                }
              for (unsigned j = 0; j < 3; j++)
                {
                  //tri_local_entities[j] = elem_nodes[faces[iface][j]].entity();
                  tri_local_entities_non_outward_normal[j] = elem_nodes[tri_faces_not_outward_normal[iface-3][j]].entity();
                }
              for (unsigned j = 0; j < 3; j++)
                {
                  int q_edge = tri_faces_not_outward_normal[iface - 3][j+3];  // edge of the face (quadratic node numbering)
                  int l_edge = edge_map_rev[q_edge - DiscretizeWedge::edge_offset];  // extract Shards index of edge
                  VERIFY_OP_ON(((0 <= l_edge) && (l_edge < nedges)), == , true, "l_edge");
                  tri_edge_marks[j] = edge_marks[l_edge];
                  //if (edge_marks[l_edge])
                  //  ++num_tri_edge_marks;
                }

              TriangulateTri tt;
              tt.triangulate_face(eMesh, tri_local_entities_non_outward_normal, tri_edge_marks, elems_tri[iface-3]);
              // since TriangulateTri uses node coordinate positions to choose a diagonal,
              //   we can have the situation where the two opposing faces are not triangulated
              //   the same way, in which case we default to the method of introducing a centroid
              //   node and joining to the faces
              if (0 && iface == 3)
                {
                  for (unsigned ii=0; ii < elems_tri[0].size(); ++ii)
                    {
                      tri_tuple_type_local tl;
                      tl[0] = elems_tri[0][ii][0];
                      tl[1] = elems_tri[0][ii][2];
                      tl[2] = elems_tri[0][ii][1];
                      elems_tri[0][ii] = tl;
                    }
                }
              if (iface == 4)
                {
                  VERIFY_OP_ON(elems_tri[0].size(), ==, elems_tri[1].size(), "hmm");
                  for (unsigned ii=0; ii < elems_tri[0].size(); ++ii)
                    {
                      if (0)
                        std::cout << "P[" << eMesh.get_rank() << "] tmp DiscretizeWedge::checking special_te_case_2 "
                                  << " elems_tri[0][ii]= " << elems_tri[0][ii]
                                  << " elems_tri[1][ii]= " << elems_tri[1][ii]  << std::endl;

                      if (elems_tri[0][ii] != elems_tri[1][ii])
                        {
                          if (m_debug)
                            {
                              std::ostringstream ostr;
                              for (unsigned jj=0; jj < elems_tri[0].size(); ++jj)
                                {
                                  ostr
                                    << " elems_tri[0][jj]= " << elems_tri[0][jj]
                                    << " elems_tri[1][jj]= " << elems_tri[1][jj]  << std::endl;
                                }

                              std::cout << "P[" << eMesh.get_rank() << "] tmp DiscretizeWedge::turn off special_te_case_2 "
                                        << " elems_tri[0][ii]= " << elems_tri[0][ii]
                                        << " elems_tri[1][ii]= " << elems_tri[1][ii]  << "\n" << ostr.str() << std::endl;
                            }
                          special_te_case_2 = -1;
                          break;
                        }
                    }
                }
            }
          if (special_te_case_2 >= 0)
            {
              wedge_to_wedge_tuple_type_local hw;
              for (unsigned ii=0; ii < elems_tri[0].size(); ++ii)
                {
                  hw[0]   = map_to_wedge(tri_local_nodes[0][elems_tri[0][ii][0]], edge_map_rev, face_map_rev);
                  hw[1]   = map_to_wedge(tri_local_nodes[0][elems_tri[0][ii][1]], edge_map_rev, face_map_rev);
                  hw[2]   = map_to_wedge(tri_local_nodes[0][elems_tri[0][ii][2]], edge_map_rev, face_map_rev);

                  hw[3]   = map_to_wedge(tri_local_nodes[1][elems_tri[1][ii][0]], edge_map_rev, face_map_rev);
                  hw[4]   = map_to_wedge(tri_local_nodes[1][elems_tri[1][ii][1]], edge_map_rev, face_map_rev);
                  hw[5]   = map_to_wedge(tri_local_nodes[1][elems_tri[1][ii][2]], edge_map_rev, face_map_rev);

                  elems_wedge.push_back(hw);
                }
              return;
            }
        }

      //if ( allow_special_refine && (special_te_case_2 || special_te || special_iso))
      if (0 && allow_special_refine )
        {
          std::vector<unsigned> fm(&face_marks[0], &face_marks[0] + nfaces);
          std::vector<unsigned> em(&edge_marks[0], &edge_marks[0] + nedges);
          std::cout << "P[" << eMesh.get_rank() << "] tmp DiscretizeWedge::num_edges_marked for wedge[ " << eMesh.identifier(wedge)
                    << "] = " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                    << " allow_special_refine= " << allow_special_refine
                    << " special_iso= " << special_iso
                    << " special_te= " << special_te
                    << " special_te_case_2= " << special_te_case_2
                    << " em= " << em
                    << " \nfm= " << fm
                    << std::endl;

          std::set<stk::mesh::Entity> neighbors;
          eMesh.get_node_neighbors(wedge, neighbors);

          {
            std::set<stk::mesh::Entity>::iterator it;
            for (it = neighbors.begin(); it != neighbors.end(); ++it)
              {
                stk::mesh::Entity pyr = *it;
                int face_0=0, face_1=0;
                bool ifn = eMesh.is_face_neighbor(wedge, pyr, &face_0, &face_1);
                if (ifn && face_0 < 3)
                  {
                    int *refine_field_pyr = stk::mesh::field_data( *eMesh.m_refine_field, pyr );
                    std::cout << " found pyramid: " << " across face: " << face_0 << " topo: " << eMesh.bucket(pyr).topology() << " refine_field_pyr= " << refine_field_pyr[0] << std::endl;
                  }
              }
          }

          eMesh.dump_vtk("err."+toString(eMesh.get_rank())+".vtk", false, &neighbors);
          neighbors.clear();
          neighbors.insert(wedge);
          eMesh.dump_vtk("werr."+toString(eMesh.get_rank())+".vtk", false, &neighbors);


        }

      //VERIFY_OP_ON(allow_special_refine, == , false, "bad pattern - please contact developer");

      bool simple_discretization = true; // for future expansion using other patterns
      if (simple_discretization)
        {
          // for each face, triangulate face, then join to centroid

          // quad faces
          for (unsigned iface=0; iface < 3; ++iface)
            {
              unsigned quad_edge_marks[4] = {0};
              unsigned num_quad_edge_marks = 0;
              unsigned quad_local_nodes[9] = {0};
              for (unsigned j = 0; j < 9; j++)
                {
                  quad_local_nodes[j] = faces[iface][j];
                }
              for (unsigned j = 0; j < 4; j++)
                {
                  int q_edge = faces[iface][j+4];  // edge of the face (quadratic node numbering)
                  int l_edge = edge_map_rev[q_edge - DiscretizeWedge::edge_offset];  // extract Shards index of edge
                  VERIFY_OP_ON(((0 <= l_edge) && (l_edge < nedges)), == , true, "l_edge");
                  quad_edge_marks[j] = edge_marks[l_edge];
                  if (edge_marks[l_edge])
                    ++num_quad_edge_marks;
                }

              std::vector<quad_to_tri_tuple_type_local> elems_tri;
              std::vector<quad_to_quad_tuple_type_local> elems_quad;
              bool use_only_tris = false;
              bool avoid_centroid_node = true;
              TriangulateQuad tq(use_only_tris, avoid_centroid_node);
              tq.triangulate_quad_face(quad_edge_marks,
                                       elems_tri, elems_quad);

              // if no edges are marked, we still need to create a prism joining the centroid to the quad face
              if (elems_tri.size() == 0 && elems_quad.size() == 0)
                {
                  VERIFY_OP_ON(num_quad_edge_marks, ==, 0, "hmm");
                  elems_quad.push_back({0,1,2,3});
                }
              if (m_debug)
                {
                  std::vector<unsigned> qem(&quad_edge_marks[0], &quad_edge_marks[0]+4);
                  //std::vector<unsigned> em(&edge_marks[0], &edge_marks[0]+12);
                  //std::vector<unsigned> fm(&face_marks[0], &face_marks[0]+6);
                  std::cout << "tmp srk TriangulateQuad[iface= " << iface << "]: num_quad_edge_marks= " << num_quad_edge_marks
                            << " elems_tri.size()= " << elems_tri.size() << " elems_quad.size= " << elems_quad.size()
                            << " Wedge num_edges_marked= " << num_edges_marked << " Wedge num_faces_marked= " << num_faces_marked
                            << "\n quad face edge marks= " << qem
                            << std::endl;
                }
              for (unsigned iquad=0; iquad < elems_quad.size(); ++iquad)
                {
                  if (m_debug)
                    {
                      std::cout << "tmp srk TriangulateQuad[iface= " << iface << "] elems_quad[" << iquad << "] = " << elems_quad[iquad];
                    }

                  wedge_to_pyr_tuple_type_local hp;
                  // reverse order to get proper volume - not here since faces have proper outward-pointing normals
                  hp[3] = elems_quad[iquad][0];
                  hp[2] = elems_quad[iquad][1];
                  hp[1] = elems_quad[iquad][2];
                  hp[0] = elems_quad[iquad][3];
                  for (unsigned ii=0; ii < 4; ++ii)
                    {
                      int hpii = quad_local_nodes[hp[ii]];
                      if (m_debug) std::cout << " , " << hpii;
                      hp[ii] = map_to_wedge(hpii, edge_map_rev, face_map_rev);
                    }
                  hp[4] = mapped_centroid;  //centroid;
                  if (m_debug) std::cout << "\n hp= " << hp << std::endl;
                  elems_pyr.push_back(hp);
                }
              for (unsigned itri=0; itri < elems_tri.size(); ++itri)
                {
                  if (m_debug)
                    {
                      std::cout << "tmp srk TriangulateQuad[iface= " << iface << "] elems_tri[" << itri << "] = " << elems_tri[itri];
                    }
                  wedge_to_tet_tuple_type_local ht;
                  // reverse order to get proper volume
                  ht[2] = elems_tri[itri][0];
                  ht[1] = elems_tri[itri][1];
                  ht[0] = elems_tri[itri][2];
                  for (unsigned ii=0; ii < 3; ++ii)
                    {
                      int htii = quad_local_nodes[ht[ii]];
                      if (m_debug) std::cout << " , " << htii;
                      ht[ii] = map_to_wedge(htii, edge_map_rev, face_map_rev);
                    }
                  ht[3] = mapped_centroid;  //centroid;
                  if (m_debug) std::cout << "\n ht= " << ht << std::endl;
                  elems_tet.push_back(ht);
                }
            }

          // triangle faces
          for (unsigned iface=3; iface < 5; ++iface)
            {
              unsigned tri_edge_marks[3] = {0};
              unsigned num_tri_edge_marks = 0;
              unsigned tri_local_nodes[6] = {0};
              stk::mesh::Entity tri_local_entities[3] = {stk::mesh::Entity()};
              for (unsigned j = 0; j < 6; j++)
                {
                  tri_local_nodes[j] = faces[iface][j];
                }
              for (unsigned j = 0; j < 3; j++)
                {
                  tri_local_entities[j] = elem_nodes[faces[iface][j]].entity();
                }
              for (unsigned j = 0; j < 3; j++)
                {
                  int q_edge = faces[iface][j+3];  // edge of the face (quadratic node numbering)
                  int l_edge = edge_map_rev[q_edge - DiscretizeWedge::edge_offset];  // extract Shards index of edge
                  VERIFY_OP_ON(((0 <= l_edge) && (l_edge < nedges)), == , true, "l_edge");
                  tri_edge_marks[j] = edge_marks[l_edge];
                  if (edge_marks[l_edge])
                    ++num_tri_edge_marks;
                }

              std::vector<tri_tuple_type_local> elems_tri;
              TriangulateTri tt;
              tt.triangulate_face(eMesh, tri_local_entities, tri_edge_marks, elems_tri);

              // if no edges are marked, we still need to create a tet joining the centroid to the tri face
              if (elems_tri.size() == 0)
                {
                  VERIFY_OP_ON(num_tri_edge_marks, ==, 0, "hmm");
                  elems_tri.push_back({0,1,2});
                }

              if (m_debug)
                {
                  std::vector<unsigned> qem(&tri_edge_marks[0], &tri_edge_marks[0]+3);
                  std::cout << "tmp srk TriangulateTri[iface= " << iface << "] num_tri_edge_marks= " << num_tri_edge_marks
                            << " elems_tri.size()= " << elems_tri.size() << " elems_tri.size= " << elems_tri.size()
                            << " num_edges_marked= " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                            << "\n tri face edge marks= " << qem
                            << std::endl;
                }
              for (unsigned itri=0; itri < elems_tri.size(); ++itri)
                {
                  if (m_debug)
                    {
                      std::cout << "tmp srk TriangulateTri[ " << iface << "] elems_tri[" << itri << "] = " << elems_tri[itri];
                    }
                  wedge_to_tet_tuple_type_local ht;
                  // reverse order to get proper volume
                  ht[2] = elems_tri[itri][0];
                  ht[1] = elems_tri[itri][1];
                  ht[0] = elems_tri[itri][2];
                  for (unsigned ii=0; ii < 3; ++ii)
                    {
                      int htii = tri_local_nodes[ht[ii]];
                      if (m_debug) std::cout << " , " << htii;
                      ht[ii] = map_to_wedge(htii, edge_map_rev, face_map_rev);
                    }
                  ht[3] = mapped_centroid;  //centroid;
                  if (m_debug) std::cout << "\n ht= " << ht << std::endl;
                  elems_tet.push_back(ht);
                }
            }

        }
    }
  };

}

#endif
