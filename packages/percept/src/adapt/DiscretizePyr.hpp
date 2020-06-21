// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_DiscretizePyr_hpp
#define percept_DiscretizePyr_hpp

#include <percept/PerceptMesh.hpp>
#include <percept/PerceptBoostArray.hpp>
#include <adapt/TriangulateQuad.hpp>
#include <adapt/TriangulateTri.hpp>

namespace percept {

  typedef std::array<unsigned, 4> pyr_to_tet_tuple_type_local;
  typedef std::array<unsigned, 5> pyr_to_pyr_tuple_type_local;
  typedef std::array<stk::mesh::EntityId, 4> pyr_to_tet_tuple_type;
  typedef std::array<stk::mesh::EntityId, 5> pyr_to_pyr_tuple_type;


  /**
   *
   *   Convention: input is the element's nodes and the marks on the edges.  Output is an array
   *     of "elements" defined as local id's of nodes forming those elements, where {0,..,4} represent
   *     the original vertices and {5..12} are the edges, {13} are faces and {14} is the centroid used
   *     to join to faces to discretize the pyramid
   *     (these are numberings used by Shards for the extra/quadratic nodes, and so we map these
   *      numbers back to the local edge/face/centroid using tables; for example, quadratic node 9
   *      corresponds to edge {0,4} see table "edges" and edge_map_rev below)
   *
   *
   *   typedef
   *     MakeTypeList< IndexList< 0 , 1 ,   5 > ,
   *                   IndexList< 1 , 2 ,   6 > ,
   *                   IndexList< 2 , 3 ,   7 > ,
   *                   IndexList< 3 , 0 ,   8 > ,
   *                   IndexList< 0 , 4 ,   9 > ,
   *                   IndexList< 1 , 4 ,  10 > ,
   *                   IndexList< 2 , 4 ,  11 > ,
   *                   IndexList< 3 , 4 ,  12 > >::type
   *     PyramidEdgeNodeMap ;
   *
   *   typedef
   *     MakeTypeList< IndexList< 0, 1, 4,     5, 10,  9 > ,
   *                   IndexList< 1, 2, 4,     6, 11, 10 > ,
   *                   IndexList< 2, 3, 4,     7, 12, 11 > ,
   *                   IndexList< 3, 0, 4,     8,  9, 12 > ,
   *                   IndexList< 0, 3, 2, 1,  8,  7,  6,  5,  13 > >::type
   *     PyramidFaceNodeMap ;
   *
   *
   *   Linear 5-Node Pyramid node locations.
   *
   *  4  - top vertex of pyramid
   *  13 - centroid of quad face
   *
   *
   *                       7
   *            3 o--------o---------o 2
   *             / \ 12          /  /
   *            /   o       o 11   /
   *           /     \    /       /
   *          /      o 4         /
   *         /     /  \         /
   *      8 o  9 /  o \       o 6
   *       /   o   13  o 10   /
   *      /  /          \    /
   *     / /            \   /
   *    //               \ /
   *   o--------o---------o
   *  0         5          1
   *
   *
   *
   */

  class DiscretizePyr {
    static const bool m_debug = false;
  public:
    static const int nnodes = 5;
    static const int nedges = 8;
    static const int nfaces = 5;
    static const int nquad_faces = 1; // last
    static const int ntri_faces = 4;

    static const int edge_offset = 5;
    static const int face_offset = 13;  // where the quad face lives
    /// following is the starting point for tri faces in the pyramid which are the first 4
    //    so to get the 5'th, quad face, subtract face_offset_tri... 13 - 9 = 4 == 5'th face
    static const int face_offset_tri = face_offset - ntri_faces;  // == 9

    static const int centroid = 14;
    static const int mapped_centroid = 14;

    static bool is_v(int i) { return i < edge_offset; }
    static bool is_e(int i) { return (edge_offset <= i && i < edge_offset + nedges); }
    static bool is_f(int i) { return (face_offset <= i && i < face_offset + nquad_faces); }
    static bool is_c(int i) { return i == centroid; }

    static bool is_mapped_v(int i) { return is_v(i); }
    static bool is_mapped_e(int i) { return is_e(i); }
    static bool is_mapped_f(int i) { return is_f(i); }
    static bool is_mapped_c(int i) { return i == mapped_centroid; }

    static int map_to_pyr(int hp, const int *edge_map_rev, const int *face_map_rev)
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

    static void discretize(PerceptMesh& eMesh, stk::mesh::Entity pyramid,
                    unsigned edge_marks[nedges], unsigned face_marks[nfaces],
                    vector<pyr_to_tet_tuple_type_local>& elems_tet,
                    vector<pyr_to_pyr_tuple_type_local>& elems_pyr)
    {
      const percept::MyPairIterRelation elem_nodes (eMesh, pyramid, stk::topology::NODE_RANK);

      // Shards edges {node0, node1, quadratic_node}
      static const int edges[][3] = {
        { 0 , 1 ,   5 } ,
        { 1 , 2 ,   6 } ,
        { 2 , 3 ,   7 } ,
        { 3 , 0 ,   8 } ,
        { 0 , 4 ,   9 } ,
        { 1 , 4 ,  10 } ,
        { 2 , 4 ,  11 } ,
        { 3 , 4 ,  12 } };

      static const int edge_map[] =     {5, 6, 7, 8, 9, 10, 11, 12};
      // map of quadratic node to index in edges table
      static const int edge_map_rev[] = {0, 1,  2,  3,  4,  5,  6,  7}; // edge_map_rev[edge_map[i] - 5] = i;

      // similar to edges: {node0...node3, edge0(quadraticNode),...edge4, faceQuadraticNode}
      static const int faces[][9] =  {
        { 0, 1, 4,        5, 10,  9, -1, -1, -1 } ,
        { 1, 2, 4,        6, 11, 10, -1, -1, -1 } ,
        { 2, 3, 4,        7, 12, 11, -1, -1, -1 } ,
        { 3, 0, 4,        8,  9, 12, -1, -1, -1 } ,
        { 0, 3, 2,  1,    8,  7,  6,  5,  13 }
      };
      static const int face_map[] =     { -1, -1, -1, -1, 13 };
      static const int face_map_rev[] = {  0,  1,  2,  3,  4 };   // face_map_rev[face_map[i] - 13] = i;

      (void)edge_map;
      (void)edges;
      (void)face_map;

      elems_tet.resize(0);
      elems_pyr.resize(0);

      const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< shards::Pyramid<5> >();
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
      {
        int iface = 4;
        unsigned num_nodes_on_face = face_marks[iface];
        if (num_nodes_on_face)
          {
            ++num_faces_marked;
          }
      }

      if (m_debug) std::cout << "tmp DiscretizePyr::num_edges_marked for pyramid[ " << eMesh.identifier(pyramid)
                             << "] = " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                             << std::endl;
      if (num_edges_marked == 0 && num_faces_marked == 0)
        {
          return;
        }

      bool simple_discretization = true; // for future expansion using other patterns
      if (simple_discretization)
        {
          // for each face, triangulate face, then join to centroid

          for (unsigned iface=4; iface <= 4; ++iface)
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
                  int l_edge = edge_map_rev[q_edge - DiscretizePyr::edge_offset];  // extract Shards index of edge
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
                  //std::vector<unsigned> em(&edge_marks[0], &edge_marks[0]+12);
                  //std::vector<unsigned> fm(&face_marks[0], &face_marks[0]+6);
                  std::cout << "tmp srk TriangulateQuad[iface= " << iface << "]: num_quad_edge_marks= " << num_quad_edge_marks
                            << " elems_tri.size()= " << elems_tri.size() << " elems_quad.size= " << elems_quad.size()
                            << " Pyr num_edges_marked= " << num_edges_marked << " Pyr num_faces_marked= " << num_faces_marked
                            << "\n quad face edge marks= " << qem
                            << std::endl;
                }
              for (unsigned iquad=0; iquad < elems_quad.size(); ++iquad)
                {
                  bool ldb=false;
                  if (ldb || m_debug)
                    {
                      std::cout << "tmp srk TriangulateQuad[iface= " << iface << "] elems_quad[" << iquad << "] = " << elems_quad[iquad];
                    }

                  pyr_to_pyr_tuple_type_local hp;
                  // reverse order to get proper volume
                  hp[3] = elems_quad[iquad][0];
                  hp[2] = elems_quad[iquad][1];
                  hp[1] = elems_quad[iquad][2];
                  hp[0] = elems_quad[iquad][3];
                  for (unsigned ii=0; ii < 4; ++ii)
                    {
                      int hpii = quad_local_nodes[hp[ii]];
                      if (ldb || m_debug) std::cout << " , " << hpii;
                      hp[ii] = map_to_pyr(hpii, edge_map_rev, face_map_rev);
                    }
                  hp[4] = mapped_centroid;  //centroid;
                  if (ldb || m_debug) std::cout << "\n hp= " << hp << std::endl;
                  elems_pyr.push_back(hp);
                }
              for (unsigned itri=0; itri < elems_tri.size(); ++itri)
                {
                  if (m_debug)
                    {
                      std::cout << "tmp srk TriangulateQuad[iface= " << iface << "] elems_tri[" << itri << "] = " << elems_tri[itri];
                    }
                  pyr_to_tet_tuple_type_local ht;
                  // reverse order to get proper volume
                  ht[2] = elems_tri[itri][0];
                  ht[1] = elems_tri[itri][1];
                  ht[0] = elems_tri[itri][2];
                  for (unsigned ii=0; ii < 3; ++ii)
                    {
                      int htii = quad_local_nodes[ht[ii]];
                      if (m_debug) std::cout << " , " << htii;
                      ht[ii] = map_to_pyr(htii, edge_map_rev, face_map_rev);
                    }
                  ht[3] = mapped_centroid;  //centroid;
                  if (m_debug) std::cout << "\n ht= " << ht << std::endl;
                  elems_tet.push_back(ht);
                }
            }

          // triangle faces
          for (unsigned iface=0; iface < 4; ++iface)
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
                  int l_edge = edge_map_rev[q_edge - DiscretizePyr::edge_offset];  // extract Shards index of edge
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
                  pyr_to_tet_tuple_type_local ht;
                  // reverse order to get proper volume
                  ht[2] = elems_tri[itri][0];
                  ht[1] = elems_tri[itri][1];
                  ht[0] = elems_tri[itri][2];
                  for (unsigned ii=0; ii < 3; ++ii)
                    {
                      int htii = tri_local_nodes[ht[ii]];
                      if (m_debug) std::cout << " , " << htii;
                      ht[ii] = map_to_pyr(htii, edge_map_rev, face_map_rev);
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
