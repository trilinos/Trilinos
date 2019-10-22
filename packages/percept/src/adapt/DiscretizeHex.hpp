// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_DiscretizeHex_hpp
#define percept_DiscretizeHex_hpp

#include <percept/PerceptMesh.hpp>
#include <percept/PerceptBoostArray.hpp>
#include <adapt/TriangulateQuad.hpp>

namespace percept {

  typedef std::array<unsigned, 4> hex_to_tet_tuple_type_local;
  typedef std::array<unsigned, 5> hex_to_pyr_tuple_type_local;
  typedef std::array<stk::mesh::EntityId, 4> hex_to_tet_tuple_type;
  typedef std::array<stk::mesh::EntityId, 5> hex_to_pyr_tuple_type;

  typedef std::array<unsigned, 8> hex_to_hex_tuple_type_local;
  typedef std::array<stk::mesh::EntityId, 8> hex_to_hex_tuple_type;

  /**
   *
   *   Convention: input is the element's nodes and the marks on the edges.  Output is an array
   *     of "elements" defined as local id's of nodes forming those elements, where {0,..,7} represent
   *     the original vertices and {8..19} are the edges, {21..26} are faces and {20} is the centroid
   *     (these are numberings used by Shards for the extra/quadratic nodes, and so we map these
   *      numbers back to the local edge/face/centroid using tables; for example, quadratic node 12
   *      corresponds to edge {0,4} which is not the 5th edge (12-8 + 1 = 5), but is the 7th in
   *      Shards' convention - see table "edges" and edge_map_rev below)
   *
   *   Linear 8-Node Hexahedron node locations.
   *
   *          7                    6
   *           o------------------o
   *          /|                 /|
   *         / |                / |
   *        /  |               /  |
   *       /   |              /   |
   *      /    |             /    |
   *     /     |            /     |
   *  4 /      |         5 /      |
   *   o------------------o       |
   *   |       |          |       |
   *   |     3 o----------|-------o 2
   *   |      /           |      /
   *   |     /            |     /
   *   |    /             |    /
   *   |   /              |   /
   *   |  /               |  /
   *   | /                | /
   *   |/                 |/
   *   o------------------o
   *  0                    1
   *
   *
   *
   *           7         18         6
   *            o--------o---------o
   *           /|       /|        /|
   *          / |      / |       / |
   *         /  |   22   |      /  |
   *      19o ------ T ----- 17o   |
   *       /| 15o --/-- B26 - /|---o14
   *      / |  /|  /     |   / |  /|
   *   4 /  | / | 16     |  /  | / |             20 centroid
   *    o---------o--------o 5 |/  |             21 bottom
   *    | 23L   | |     10 |   R24 |             22 top
   *    |  /| 3 o-|-------o|--/|---o 2           23 left
   *    | / |  /  |      / | / |  /              24 right
   *    |/  | /  25     /  |/  | /               25 front
   *  12o-------- F ---/---o13 |/                26 back
   *    |   o 11 -| B21--- |-- o 9
   *    |  /      |  /     |  /
   *    | /       | /      | /
   *    |/        |/       |/
   *    o---------o--------o
   *   0          8         1
   *
   *
   *           H         S          G
   *            o--------o---------o
   *           /|       /|        /|
   *          / |      / |       / |
   *         /  |     /  |      /  |
   *      T o ------ V ----- R o   |
   *       /| P o --/--  X -- /|---o O
   *      / |  /|  /     |   / |  /|
   *   E /  | / | Q      | F/  | / |
   *    o---------o--------o   |/  |
   *    |   Y   | |      K |   Z   |
   *    |  /| D o-|-------o|--/|---o C
   *    | / |  /  |      / | / |  /
   *    |/  | /   |     /  |/  | /
   *  M o-------- W ---/---o N |/
   *    |   o L  -|   U -- |-- o J
   *    |  /      |  /     |  /
   *    | /       | /      | /
   *    |/        |/       |/
   *    o---------o--------o
   *   A          I         B
   *
   *
   *
   */

  class DiscretizeHex {
    static const bool m_debug = false;
  public:
    static const int nedges = 12;
    static const int nfaces = 6;
    static const int edge_offset = 8;
    static const int face_offset = 21;

    static const int centroid = 20;
    static const int mapped_centroid = 20;

    static bool is_v(int i) { return i < 8; }
    static bool is_e(int i) { return (edge_offset <= i && i < edge_offset + nedges); }
    static bool is_f(int i) { return (face_offset <= i && i < face_offset + nfaces); }
    static bool is_c(int i) { return i == centroid; }

    static bool is_mapped_v(int i) { return is_v(i); }
    static bool is_mapped_e(int i) { return is_e(i); }
    static bool is_mapped_f(int i) { return is_f(i); }
    static bool is_mapped_c(int i) { return i == mapped_centroid; }

    static int map_to_hex(int hp, const int *edge_map_rev, const int *face_map_rev)
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

    static void discretize(
                    unsigned edge_marks[nedges], unsigned face_marks[nfaces],
                    vector<hex_to_tet_tuple_type_local>& elems_tet,
                    vector<hex_to_pyr_tuple_type_local>& elems_pyr,
                    vector<hex_to_hex_tuple_type_local>& elems_hex )
    {
      // Shards edges {node0, node1, quadratic_node}
      static const int edges[][3] = {
        { 0 , 1 ,   8 } ,
        { 1 , 2 ,   9 } ,
        { 2 , 3 ,  10 } ,
        { 3 , 0 ,  11 } ,
        { 4 , 5 ,  16 } ,
        { 5 , 6 ,  17 } ,
        { 6 , 7 ,  18 } ,
        { 7 , 4 ,  19 } ,
        { 0 , 4 ,  12 } ,
        { 1 , 5 ,  13 } ,
        { 2 , 6 ,  14 } ,
        { 3 , 7 ,  15 } };
      static const int edge_map[] =     {8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
      // map of quadratic node to index in edges table
      static const int edge_map_rev[] = {0, 1,  2,  3,  8,  9, 10, 11,  4,  5,  6,  7}; // edge_map_rev[edge_map[i] - 8] = i;


      // similar to edges: {node0...node3, edge0(quadraticNode),...edge4, faceQuadraticNode}
      static const int faces[][9] =  {
        { 0, 1, 5, 4,   8, 13, 16, 12,   25 } ,
        { 1, 2, 6, 5,   9, 14, 17, 13,   24 } ,
        { 2, 3, 7, 6,  10, 15, 18, 14,   26 } ,
        { 0, 4, 7, 3,  12, 19, 15, 11,   23 } ,
        { 0, 3, 2, 1,  11, 10,  9,  8,   21 } ,
        { 4, 5, 6, 7,  16, 17, 18, 19,   22 }
      };
      static const int face_map[] =     { 25, 24, 26, 23, 21, 22 };
      static const int face_map_rev[] = {  4,  5,  3,  1,  0,  2 };   // face_map_rev[face_map[i] - 21] = i;

      (void)edge_map;
      (void)edges;
      (void)face_map;

      elems_tet.resize(0);
      elems_pyr.resize(0);
      elems_hex.resize(0);

      const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< shards::Hexahedron<8> >();
      shards::CellTopology cell_topo(cell_topo_data);

      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < nedges; iedge++)
        {
          unsigned num_nodes_on_edge = edge_marks[iedge];
          if (num_nodes_on_edge)
            {
              ++num_edges_marked;
            }
        }
      unsigned num_faces_marked=0;
      for (int iface = 0; iface < nfaces; iface++)
        {
          unsigned num_nodes_on_face = face_marks[iface];
          if (num_nodes_on_face)
            {
              ++num_faces_marked;
            }
        }

      //std::cout << "tmp RefinerPattern_Tri3_Tri3_N::num_edges_marked= " << num_edges_marked << std::endl;
      if (num_edges_marked == 0 && num_faces_marked == 0)
        {
          return;
        }

      bool simple_discretization = true; // for future expansion using other patterns
      if (simple_discretization)
        {
          // for each face, triangulate face, then join to centroid

          for (int iface=0; iface < nfaces; ++iface)
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
                  int l_edge = edge_map_rev[q_edge - DiscretizeHex::edge_offset];  // extract Shards index of edge
                  VERIFY_OP_ON(((0 <= l_edge) && (l_edge < nedges)), == , true, "l_edge");
                  quad_edge_marks[j] = edge_marks[l_edge];
                  if (edge_marks[l_edge])
                    ++num_quad_edge_marks;
                }

              std::vector<quad_to_tri_tuple_type_local> elems_tri;
              std::vector<quad_to_quad_tuple_type_local> elems_quad;
              TriangulateQuad tq(false,true);
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
                  //std::vector<unsigned> em(&edge_marks[0], &edge_marks[0]+nedges);
                  //std::vector<unsigned> fm(&face_marks[0], &face_marks[0]+nfaces);
                  std::cout << "tmp srk num_quad_edge_marks= " << num_quad_edge_marks
                            << " elems_tri.size()= " << elems_tri.size() << " elems_quad.size= " << elems_quad.size()
                            << " num_edges_marked= " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                            << "\n qem= " << qem
                            << std::endl;
                }
              for (unsigned iquad=0; iquad < elems_quad.size(); ++iquad)
                {
                  if (m_debug)
                    {
                      std::cout << "tmp srk elems_quad[" << iquad << "] = " << elems_quad[iquad];
                    }

                  hex_to_pyr_tuple_type_local hp;
                  // reverse order to get proper volume
                  hp[3] = elems_quad[iquad][0];
                  hp[2] = elems_quad[iquad][1];
                  hp[1] = elems_quad[iquad][2];
                  hp[0] = elems_quad[iquad][3];
                  for (unsigned ii=0; ii < 4; ++ii)
                    {
                      int hpii = quad_local_nodes[hp[ii]];
                      if (m_debug) std::cout << " , " << hpii;
                      hp[ii] = map_to_hex(hpii, edge_map_rev, face_map_rev);
                    }
                  hp[4] = mapped_centroid;  //centroid;
                  if (m_debug) std::cout << "\n hp= " << hp << std::endl;
                  elems_pyr.push_back(hp);
                }
              for (unsigned itri=0; itri < elems_tri.size(); ++itri)
                {
                  if (m_debug)
                    {
                      std::cout << "tmp srk elems_tri[" << itri << "] = " << elems_tri[itri];
                    }
                  hex_to_tet_tuple_type_local ht;
                  // reverse order to get proper volume
                  ht[2] = elems_tri[itri][0];
                  ht[1] = elems_tri[itri][1];
                  ht[0] = elems_tri[itri][2];
                  for (unsigned ii=0; ii < 3; ++ii)
                    {
                      int htii = quad_local_nodes[ht[ii]];
                      if (m_debug) std::cout << " , " << htii;
                      ht[ii] = map_to_hex(htii, edge_map_rev, face_map_rev);
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
