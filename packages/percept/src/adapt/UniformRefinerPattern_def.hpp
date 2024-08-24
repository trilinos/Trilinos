// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_def_hpp
#define adapt_UniformRefinerPattern_def_hpp

#include <adapt/UniformRefinerPattern.hpp>

  namespace percept {

    /// utility methods for converting Sierra tables to new format (which groups DOF's on sub-entities)
    template<typename FromTopology,  typename ToTopology >
    bool
    URP<FromTopology, ToTopology>::
    on_parent_vertex(unsigned childNodeIdx)
    {
      return (childNodeIdx < FromTopology::vertex_count);
    }

    template<typename FromTopology,  typename ToTopology >
    bool
    URP<FromTopology, ToTopology>::
    on_parent_edge(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo)
    {
      unsigned num_child_nodes = ref_topo.num_child_nodes();

      shards::CellTopology cell_topo ( shards::getCellTopologyData< FromTopology >() );
      unsigned n_edges = cell_topo.getEdgeCount();

      for (unsigned i_edge = 0; i_edge < n_edges; i_edge++)
        {
          const UInt *edge_nodes = ref_topo.edge_node(i_edge);

          unsigned i_ord = 0;
          for (unsigned i_edge_n = 0; edge_nodes[i_edge_n] != END_UINT_ARRAY; i_edge_n++)
            {
              unsigned j_e_node = edge_nodes[i_edge_n];
              if (childNodeIdx == j_e_node)
                {
                  return true;
                }
              if (j_e_node < num_child_nodes)
                ++i_ord;
            }
        }
      return false;
    }

    template<typename FromTopology,  typename ToTopology >
    bool
    URP<FromTopology, ToTopology>::
    on_parent_face(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo)
    {
      unsigned num_child_nodes = ref_topo.num_child_nodes();

      shards::CellTopology cell_topo ( shards::getCellTopologyData< FromTopology >() );
      if (cell_topo.getDimension() == 2)
        {
          return true; // by definition
        }

      unsigned n_faces = cell_topo.getFaceCount();
      if (n_faces == 0) n_faces = 1; // 2D face has one "face"

      for (unsigned i_face = 0; i_face < n_faces; i_face++)
        {
          const UInt *face_nodes = ref_topo.face_node(i_face);

          unsigned i_ord = 0;
          for (unsigned i_face_n = 0; face_nodes[i_face_n] != END_UINT_ARRAY; i_face_n++)
            {
              unsigned j_e_node = face_nodes[i_face_n];
              if (childNodeIdx == j_e_node)
                {
                  return true;
                }
              if (j_e_node < num_child_nodes)
                ++i_ord;
            }
        }
      return false;
    }

    template<typename FromTopology,  typename ToTopology >
    bool
    URP<FromTopology, ToTopology>::
    on_parent_edge_interior(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo, unsigned& i_edge, unsigned& i_ord, unsigned& n_ord)
    {
      if (on_parent_vertex(childNodeIdx))
        return false;

      unsigned num_child_nodes = ref_topo.num_child_nodes();
      shards::CellTopology cell_topo ( shards::getCellTopologyData< FromTopology >() );
      unsigned n_edges = cell_topo.getEdgeCount();
      if (n_edges == 0) n_edges = 1;

      for ( i_edge = 0; i_edge < n_edges; i_edge++)
        {
          const UInt *edge_nodes = ref_topo.edge_node(i_edge);

          n_ord = 0;
          int n_edge_n = 0;
          for (unsigned i_edge_n = 0; edge_nodes[i_edge_n] != END_UINT_ARRAY; i_edge_n++)
            {
              if (on_parent_vertex(edge_nodes[i_edge_n]))
                continue;
              if (edge_nodes[i_edge_n] < num_child_nodes)
                ++n_ord;
              ++n_edge_n;
            }

          i_ord = 0;
          for (unsigned i_edge_n = 0; edge_nodes[i_edge_n] != END_UINT_ARRAY; i_edge_n++)
            // go in reverse to put mid node at the end
            //for (int i_edge_n = n_edge_n-1; i_edge_n >= 0; i_edge_n--)
            {
              if (on_parent_vertex(edge_nodes[i_edge_n]))
                continue;
              unsigned j_e_node = edge_nodes[i_edge_n];
              if (childNodeIdx == j_e_node)
                {
                  if (i_ord == 0)
                    i_ord = n_ord-1;
                  else
                    --i_ord;

                  return true;
                }
              if (j_e_node < num_child_nodes)
                ++i_ord;
            }
        }
      return false;
    }

    /** SPECIAL CASE ALERT (see below for triangle faces)
     *
     *  To prepare for future truly hierarchical elements, we want to number the nodes on the interior of each face
     *  in a manner that mimics the parent element.  For example, nodes {8, 17-20, and 21-24} below are not numbered
     *  consistent with the parent quadratic element.  The consistent numbering would be {21-24, 17-20, 8} to correspond
     *  with parent's {0-8} numbering:
     *
     * After refinement:
     *
     *  3    14    6   13     2   CHILD 9-Node Quadrilateral Object Nodes
     *   o----*----o----*----o    (new nodes = *)
     *   |         |         |
     *   |   24    |    23   |
     * 15*    *    *19  *    *12
     *   |         |         |
     *   |        8|    18   |
     * 7 o----*----o----*----o 5
     *   |   20    |         |
     *   |         |         |
     * 16*    *  17*    *    *11
     *   |   21    |   22    |
     *   |         |         |
     *   o----*----o----*----o
     *  0     9    4   10     1
     *
     *
     * After refinement:
     *
     *  3    14    6   13     2   CHILD 8-Node Quadrilateral Object Nodes
     *   o----*----o----*----o    (new nodes = *)
     *   |         |         |
     *   |         |         |
     * 15*         *19       *12    This case is so similar to the full quadratic case we just re-use the node numbering and
     *   |         |         |         go ahead and generate 9 nodes per face, which is 4 more than we need, but we drop them later.
     *   |        8|    18   |
     * 7 o----*----o----*----o 5
     *   |   20    |         |
     *   |         |         |
     * 16*       17*         *11
     *   |         |         |
     *   |         |         |
     *   o----*----o----*----o
     *  0     9    4   10     1
     *
     * The way this meshes with hierarchical elements is to imagine the face being
     *  part of a p-element with local polynomial degree 4x3 (p=4 in local-x of the face, p=3 in local-y).  In that
     *  case, we can think of the nodes laid out in a rectangular pattern as:
     *
     *      o---o---o
     *      |   |   |
     *      o---o---o
     *
     *  or in the parent's hierarchical element (e.g. node 5 is the hiearchical edge-bubble DOF), and numbered
     *  in a lexicographical (x-first, then y) ordering:
     *
     *  3        6        2
     *   o-------o-------o
     *   |               |
     *   |   11  12  13  |
     *   |    o--o--o    |
     * 7 o    |  |  |    o 5
     *   |    o--o--o    |
     *   |   8   9   10  |
     *   |               |
     *   o-------o-------o
     *  0        4        1
     *
     *
     *  Of course, this is only symbolic of the hiearchical DOF's at the center of the face, shown as Lagrange-type DOF's for
     *  exposition only.  In reality, we denote just with a + at the center:
     *
     *            6
     *   3        _       2
     *    o---------------o
     *    |               |
     *    |               |
     *    |         8     |
     * 7 ||       +       || 5
     *    |               |
     *    |               |
     *    |               |
     *    o---------------o
     *   0        -        1
     *            4
     *
     *  So, we renumber the nodes on the face interior as: {21-24, 17-20, 8} - this is used below to choose face DOF ordinals
     *  in building the tables in printRefinementTopoX_Table
     *
     */

    template<typename FromTopology,  typename ToTopology >
    unsigned
    URP<FromTopology, ToTopology>::
    renumber_quad_face_interior_nodes(unsigned original_node)
    {
      static int face_interior_inverse_map[] = { -1, /* 0 */
                                                 -1, -2, -3, -4, -5, -6, -7,
                                                 8,  /* 8 */
                                                 -9, -10,
                                                 -11, -12, -13, -14, -15, -16,
                                                 4, 5, 6, 7, // -17, -18, -19, -20
                                                 0, 1, 2, 3 }; //-21, -22, -23, -24};

      /*
        static int face_interior_map[] = {21, 22, 23, 24,
        17, 18, 19, 20,
        8 };
      */
      if (original_node >= 25) throw std::logic_error("renumber_quad_face_interior_nodes 1");
      int val = face_interior_inverse_map[original_node];
      if (val < 0) throw std::logic_error("renumber_quad_face_interior_nodes 2");
      return (unsigned)val;
    }

    // not used (yet)
    template<typename FromTopology,  typename ToTopology >
    unsigned
    URP<FromTopology, ToTopology>::
    renumber_quad_face_interior_nodes_quad8(unsigned original_node)
    {
      static int face_interior_inverse_map[] = { -1, /* 0 */
                                                 -1, -2, -3, -4, -5, -6, -7,
                                                 4,  /* 8 */
                                                 -9, -10,
                                                 -11, -12, -13, -14, -15, -16,
                                                 0, 1, 2, 3 // -17, -18, -19, -20
      };

      /*
        static int face_interior_map[] = {21, 22, 23, 24,
        17, 18, 19, 20,
        8 };
      */
      if (original_node >= 21) throw std::logic_error("renumber_quad_face_interior_nodes_quad8 1");
      int val = face_interior_inverse_map[original_node];
      if (val < 0) throw std::logic_error("renumber_quad_face_interior_nodes_quad8 2");
      return (unsigned)val;
    }

    /*--------------------------------------------------------------------*/
    /**
     *           2            PARENT 6-Node Triangle Object Nodes
     *           o
     *          / \
     *         /   \          (PARENT) 6-Node Triangle Object Edge Node Map:
     *        /     \
     *     5 o       o 4      { {0, 1, 3}, {1, 2, 4}, {2, 0, 5} };
     *      /         \
     *     /           \
     *    /             \
     *   o-------o-------o
     *  0        3        1
     *
     *   After refinement:
     *
     *           2            CHILD 6-Node Triangle Object Nodes
     *           o                  (new nodes = *)
     *          / \
     *      10 *   * 9
     *        / 14  \
     *     5 o---*---o 4
     *      / \     / \
     *  11 * 12*   *13 * 8
     *    /     \ /     \
     *   o---*---o---*---o
     *  0    6   3   7    1
     *
     * | CHILD 6-Node Triangle Object Node Maps:
     * |
     * | static const UInt child_0[] = { 0, 3, 5, 6, 12, 11  };
     * | static const UInt child_1[] = { 3, 1, 4, 7, 8, 13   };
     * | static const UInt child_2[] = { 5, 4, 2, 14, 9, 10  };
     * | static const UInt child_3[] = { 4, 5, 3, 14, 12, 13 };
     * |
     *
     *  Refined 6-Node Triangle Object PERMUTATION Node Maps:
     *
     *  Rotation  Polarity
     *     0          1       { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };
     *     0          0       { 0, 2, 1, 5, 4, 3, 11, 10, 9, 8, 7, 6, 12, 14, 13 };
     *     1          1       { 2, 0, 1, 5, 3, 4, 10, 11, 6, 7, 8, 9, 14, 12, 13 };
     *     1          0       { 2, 1, 0, 4, 3, 5, 9, 8, 7, 6, 11, 10, 14, 13, 12 };
     *     2          1       { 1, 2, 0, 4, 5, 3, 8, 9, 10, 11, 6, 7, 13, 14, 12 };
     *     2          0       { 1, 0, 2, 3, 5, 4, 7, 6, 11, 10, 9, 8  13, 12, 14 };
     *
     **/

    template<typename FromTopology,  typename ToTopology >
    bool
    URP<FromTopology, ToTopology>::
    on_parent_face_interior(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo, unsigned& i_face, unsigned& i_ord, unsigned& n_ord)
    {
      if (on_parent_edge(childNodeIdx, ref_topo))
        return false;

      static bool doRenumber = true;

      unsigned num_child_nodes = ref_topo.num_child_nodes();
      shards::CellTopology cell_topo ( shards::getCellTopologyData< FromTopology >() );

      int topoDim = getTopoDim(cell_topo);
      //unsigned cell_topo_key = fromTopoKey;


      if (topoDim == 2)
        {
          i_face = 0;
          n_ord = 0;
          int n_face_n = 0;
          for (unsigned i_face_n = 0; i_face_n < num_child_nodes; i_face_n++)
            {
              if (on_parent_edge(i_face_n, ref_topo))
                continue;
              ++n_ord;
              ++n_face_n;
            }

          if (fromTopoKey == topo_key_quad8 || fromTopoKey == topo_key_shellquad8)
            {
              n_ord = 9;
              std::cout << "n_ord = " << n_ord << " for cell_topo= " << cell_topo.getName() << std::endl;
            }

          i_ord = 0;
          for (unsigned i_face_n = 0; i_face_n < num_child_nodes; i_face_n++)
            {
              if (on_parent_edge(i_face_n, ref_topo))
                continue;
              if (i_face_n == childNodeIdx)
                {
                  if (fromTopoKey == topo_key_quad9 || fromTopoKey == topo_key_quad8 || fromTopoKey == topo_key_shellquad8)
                    {
                      if (doRenumber)
                        {
                          i_ord = renumber_quad_face_interior_nodes(i_face_n);
                        }
                    }
                  return true;
                }
              ++i_ord;
            }
          return false;
        }  // cell dim == 2

      unsigned n_faces = cell_topo.getFaceCount();

      for ( i_face = 0; i_face < n_faces; i_face++)
        {
          const UInt *face_nodes = ref_topo.face_node(i_face);

          n_ord = 0;
          int n_face_n = 0;
          for (unsigned i_face_n = 0; face_nodes[i_face_n] != END_UINT_ARRAY; i_face_n++)
            {
              if (on_parent_edge(face_nodes[i_face_n], ref_topo))
                continue;
              if (face_nodes[i_face_n] < num_child_nodes)
                ++n_ord;
              ++n_face_n;
            }

          if (fromTopoKey == topo_key_hex20)
            {
              n_ord = 9;
            }
          if (fromTopoKey == topo_key_wedge15)
            {
              if (i_face <= 2)  // quad faces
                n_ord = 9;
            }
          if (fromTopoKey == topo_key_pyramid13)
            {
              if (i_face == 4)  // quad face
                n_ord = 9;
            }

          i_ord = 0;
          unsigned fnl=0;
          for (unsigned i_face_n = 0; face_nodes[i_face_n] != END_UINT_ARRAY; i_face_n++)
            {
              ++fnl;
            }
          for (unsigned i_face_n = 0; face_nodes[i_face_n] != END_UINT_ARRAY; i_face_n++)
            {
              if (on_parent_edge(face_nodes[i_face_n], ref_topo))
                continue;
              unsigned j_e_node = face_nodes[i_face_n];
              if (childNodeIdx == j_e_node)
                {
                  if (doRenumber &&

                      ( (fromTopoKey == topo_key_hex27) ||
                        (fromTopoKey == topo_key_hex20) ||
                        (fromTopoKey == topo_key_wedge15 && i_face <= 2) ||
                        (fromTopoKey == topo_key_pyramid13 && i_face == 4) )
                      )
                    {
                      i_ord = renumber_quad_face_interior_nodes(i_face_n);
                      //std::cout << "tmp childNodeIdx= " << childNodeIdx << " i_ord= " << i_ord << " i_face_n= " << i_face_n << " fnl= " << fnl <<  std::endl;
                    }

                  return true;
                }
              if (j_e_node < num_child_nodes)
                ++i_ord;
            }
        }

      return false;
    }

    template<typename FromTopology,  typename ToTopology >
    bool
    URP<FromTopology, ToTopology>::
    on_parent_volume_interior(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo, unsigned& i_volume, unsigned& i_ord, unsigned& n_ord)
    {
      if (on_parent_face(childNodeIdx, ref_topo))
        return false;

      unsigned num_child_nodes = ref_topo.num_child_nodes();
      shards::CellTopology cell_topo ( shards::getCellTopologyData< FromTopology >() );

      i_volume = 0;

      n_ord = 0;
      for (unsigned i_volume_n = 0; i_volume_n < num_child_nodes; i_volume_n++)
        {
          if (on_parent_face(i_volume_n, ref_topo))
            continue;

          ++n_ord;
        }

      i_ord = 0;
      for (unsigned i_volume_n = 0; i_volume_n < num_child_nodes; i_volume_n++)
        {
          if (on_parent_face(i_volume_n, ref_topo))
            continue;

          if (childNodeIdx == i_volume_n)
            {
              return true;
            }
          ++i_ord;
        }

      return false;
    }


    /// utility to help convert Sierra tables - this method takes the index of the child node and finds it and adds
    ///   the associated info to the ref_topo_x tables containing the new/additional refinement table information
    template<typename FromTopology,  typename ToTopology >
    void
    URP<FromTopology, ToTopology>::
    findRefinedCellTopoInfo(unsigned childNodeIdx,
                            const Elem::RefinementTopology& ref_topo,
                            RefTopoX_arr ref_topo_x,  // assumed good for the vertices
                            unsigned& rank_of_subcell,
                            unsigned& ordinal_of_subcell,
                            unsigned& ordinal_of_node_on_subcell,
                            unsigned& num_node_on_subcell)
    {
      shards::CellTopology cell_topo ( shards::getCellTopologyData< FromTopology >() );

      //bool found = false;
      bool on_parent_edge = false;
      bool on_parent_face = false;
      bool on_parent_volume = false;
      if (on_parent_vertex(childNodeIdx))
        {
          rank_of_subcell = 0;
          ordinal_of_subcell = childNodeIdx;
          ordinal_of_node_on_subcell = 0;
          num_node_on_subcell = 1;
          return;
        }



      if ( (on_parent_edge = on_parent_edge_interior(childNodeIdx, ref_topo, ordinal_of_subcell, ordinal_of_node_on_subcell, num_node_on_subcell)))
        {
          rank_of_subcell = 1;
          // SPECIAL CASE
          if (cell_topo.getKey() == base_s_beam_2_key ||
              cell_topo.getKey() == base_s_beam_3_key)
            {
              //rank_of_subcell = 3;   // wrong - see UniformRefinerPattern_Beam2_Beam2_2_sierra.hpp and ...Beam3...
              rank_of_subcell = 1;   // wrong - see UniformRefinerPattern_Beam2_Beam2_2_sierra.hpp and ...Beam3...
            }
          return;
        }

      if ( (on_parent_face = on_parent_face_interior(childNodeIdx, ref_topo, ordinal_of_subcell, ordinal_of_node_on_subcell, num_node_on_subcell)))
        {
          rank_of_subcell = 2;
          return;
        }

      // FIXME
      if (cell_topo.getDimension() == 2)
        {
        }

      if ( (on_parent_volume = on_parent_volume_interior(childNodeIdx, ref_topo, ordinal_of_subcell, ordinal_of_node_on_subcell, num_node_on_subcell)))
        {
          rank_of_subcell = 3;
          return;
        }
      throw std::logic_error("findRefinedCellTopoInfo:: hmmmm");
    }

    template<typename FromTopology,  typename ToTopology >
    void
    URP<FromTopology, ToTopology>::
    findRefinedCellParamCoords(const Elem::RefinementTopology& ref_topo,
                               RefTopoX_arr ref_topo_x)
    {
      shards::CellTopology parent_cell_topo ( shards::getCellTopologyData< FromTopology >() );

      unsigned num_child = ref_topo.num_child();

      for (unsigned i_child = 0; i_child < num_child; i_child++)
        {
          shards::CellTopology cell_topo = ref_topo.child_cell_topology(i_child);

          const unsigned *child_nodes = ref_topo.child_node(i_child);

          unsigned n_edges = cell_topo.getEdgeCount();
          if (n_edges == 0) n_edges = 1; // 1D edge has one "edge"
          unsigned n_faces = cell_topo.getFaceCount();
          if (parent_cell_topo.getDimension() > 1 && n_faces == 0) n_faces = 1; // 2D face has one "face"

          for (unsigned i_edge = 0; i_edge < n_edges; i_edge++)
            {
              // FIXME for 2d
              shards::CellTopology edge_topo = parent_cell_topo.getDimension()==1? parent_cell_topo : shards::CellTopology(cell_topo.getCellTopologyData( 1, i_edge));

              if (edge_topo.getNodeCount() == 3)
                {
                  unsigned i0 = 0;
                  unsigned i1 = 0;
                  unsigned i2 = 0;

                  if (parent_cell_topo.getDimension() == 1)
                    {
                      i0 = child_nodes[0];
                      i1 = child_nodes[1];
                      i2 = child_nodes[2];
                    }
                  else
                    {
                      i0 = child_nodes[cell_topo.getCellTopologyData()->edge[i_edge].node[0]];
                      i1 = child_nodes[cell_topo.getCellTopologyData()->edge[i_edge].node[1]];
                      i2 = child_nodes[cell_topo.getCellTopologyData()->edge[i_edge].node[2]];
                    }

                  double *param_coord = ref_topo_x[i2].parametric_coordinates;
                  param_coord[0] = (ref_topo_x[i0].parametric_coordinates[0]+ref_topo_x[i1].parametric_coordinates[0])/2.;
                  param_coord[1] = (ref_topo_x[i0].parametric_coordinates[1]+ref_topo_x[i1].parametric_coordinates[1])/2.;
                  param_coord[2] = (ref_topo_x[i0].parametric_coordinates[2]+ref_topo_x[i1].parametric_coordinates[2])/2.;

                  if (0)
                    std::cout<<"param_coord in findRefinedCellParamCoords edge= " << i2 << " "
                             << param_coord[0] << " "
                             << param_coord[1] << " "
                             << param_coord[2] << std::endl;
                }
            }

          for (unsigned i_face = 0; i_face < n_faces; i_face++)
            {
              // FIXME for 2d
              shards::CellTopology face_topo = cell_topo.getDimension()==2 ? cell_topo : shards::CellTopology(cell_topo.getCellTopologyData( 2, i_face));

              // skip triangle faces
              if (face_topo.getVertexCount() == 3)
                continue;

              // NOTE: if this is a serendipity 8-node face, it has no interior node - only 9-noded quad faces have an interior node
              if (face_topo.getNodeCount() == 9)
                {
                  unsigned i0 = cell_topo.getDimension()==2 ? 8 : cell_topo.getCellTopologyData()->side[i_face].node[8];
                  i0 = child_nodes[i0];

                  double *param_coord = ref_topo_x[i0].parametric_coordinates;
                  param_coord[0] = 0.0;
                  param_coord[1] = 0.0;
                  param_coord[2] = 0.0;
                  for (unsigned i_face_n=0; i_face_n < 4; i_face_n++)
                    {
                      unsigned i1 = cell_topo.getDimension()==2 ? i_face_n : cell_topo.getCellTopologyData()->side[i_face].node[i_face_n];
                      i1 = child_nodes[i1];
                      param_coord[0] += ref_topo_x[i1].parametric_coordinates[0]/4.;
                      param_coord[1] += ref_topo_x[i1].parametric_coordinates[1]/4.;
                      param_coord[2] += ref_topo_x[i1].parametric_coordinates[2]/4.;
                    }
                  if (0)
                    std::cout<<"param_coord in findRefinedCellParamCoords face= " << i0 << " "
                             << param_coord[0] << " "
                             << param_coord[1] << " "
                             << param_coord[2] << std::endl;

                }
            }

          if (cell_topo.getDimension() == 3 && toTopoKey == topo_key_hex27)
            {
              unsigned i0 = child_nodes[centroid_node]; // Has to be a Hex27 to have an interior node
              double *param_coord = ref_topo_x[i0].parametric_coordinates;

              for (unsigned ix=0; ix < 3; ix++)
                {
                  param_coord[ix] = 0.0;
                }
              for (unsigned k_node = 0; k_node < 8; k_node++)
                {
                  unsigned l_node = child_nodes[k_node];
                  for (unsigned ix=0; ix < 3; ix++)
                    {
                      param_coord[ix] += ref_topo_x[l_node].parametric_coordinates[ix]/8.0;
                    }
                }
              if (1)
                std::cout<<"param_coord in findRefinedCellParamCoords vol= "
                         << param_coord[0] << " "
                         << param_coord[1] << " "
                         << param_coord[2] << std::endl;

            }
        }
    }


    /// continuing in the convert tables theme, this helps to find the new nodes' parametric coordinates
    template<typename FromTopology,  typename ToTopology >
    void
    URP<FromTopology, ToTopology>::
    findRefinedCellParamCoordsLinear(const Elem::RefinementTopology& ref_topo,
                                     RefTopoX_arr ref_topo_x  // assumed good for the vertices
                                     )
    {
      double param_coord[3];
      shards::CellTopology cell_topo ( shards::getCellTopologyData< FromTopology >() );

      unsigned num_child_nodes = ref_topo.num_child_nodes();
      for (unsigned childNodeIdx = 0; childNodeIdx < num_child_nodes; childNodeIdx++)
        {
          bool found = false;
          //bool on_edge = false;
          bool on_vertex = false;
          if (childNodeIdx < FromTopology::vertex_count)
            {
              param_coord[0] = ref_topo_x[childNodeIdx].parametric_coordinates[0];
              param_coord[1] = ref_topo_x[childNodeIdx].parametric_coordinates[1];
              param_coord[2] = ref_topo_x[childNodeIdx].parametric_coordinates[2];
              found = true;
              on_vertex = true;
            }

          if (!on_vertex)
            {
              unsigned n_edges = cell_topo.getEdgeCount();
              if (n_edges == 0) n_edges = 1; // 1D face has one "edge"
              unsigned n_faces = cell_topo.getFaceCount();
              if (cell_topo.getDimension() > 1 && n_faces == 0) n_faces = 1; // 2D face has one "face"
              //unsigned n_sides = cell_topo.getSideCount();

              // check for shell line elements
              int topoDim = getTopoDim(cell_topo);
              if (topoDim == 1) n_faces = 0;

              for (unsigned i_edge = 0; i_edge < n_edges; i_edge++)
                {
                  const UInt *edge_nodes = ref_topo.edge_node(i_edge);

                  if (childNodeIdx == edge_nodes[2])  // FIXME
                    {
                      //on_edge = true;
                      found = true;
                      param_coord[0] = (ref_topo_x[edge_nodes[0]].parametric_coordinates[0]+ref_topo_x[edge_nodes[1]].parametric_coordinates[0])/2.;
                      param_coord[1] = (ref_topo_x[edge_nodes[0]].parametric_coordinates[1]+ref_topo_x[edge_nodes[1]].parametric_coordinates[1])/2.;
                      param_coord[2] = (ref_topo_x[edge_nodes[0]].parametric_coordinates[2]+ref_topo_x[edge_nodes[1]].parametric_coordinates[2])/2.;
                      break;
                    }
                }

              if (!found)
                {
                  for (unsigned i_face = 0; i_face < n_faces; i_face++)
                    {

                      // FIXME for 2d
                      shards::CellTopology face_topo = cell_topo.getDimension()==2 ? cell_topo : shards::CellTopology(cell_topo.getCellTopologyData( 2, i_face));

                      // skip triangle faces
                      if (face_topo.getVertexCount() == 3)
                        continue;

                      const UInt *face_nodes = ref_topo.face_node(i_face);

                      unsigned j_node_end = 0;
#if 0
                      bool lfnd = false;
                      for (unsigned j_node = 0; j_node < 9; j_node++)
                        {
                          if (face_nodes[j_node] == END_UINT_ARRAY)
                            {
                              lfnd = true;
                              j_node_end = j_node;
                              break;
                            }
                        }
                      if (!lfnd)
                        {
                          throw std::logic_error("findRefinedCellParamCoordsLinear logic err # 0");
                        }
#endif
                      j_node_end = 9;  //!#

                      for (unsigned j_node = 0; j_node < j_node_end; j_node++) // FIXME
                        {
                          if (cell_topo.getDimension() != 2 && face_nodes[j_node] == END_UINT_ARRAY)
                            {
                              throw std::logic_error("findRefinedCellParamCoordsLinear logic err # 1");
                            }
                          unsigned fn =  cell_topo.getDimension()==2 ? j_node : face_nodes[j_node];

                          if (childNodeIdx == fn)
                            {
                              found = true;

                              for (unsigned ix=0; ix < 3; ix++)
                                {
                                  param_coord[ix] = 0.0;
                                }

                              for (unsigned k_node = 0; k_node < 4; k_node++)
                                {
                                  if (cell_topo.getDimension() != 2 && face_nodes[k_node] == END_UINT_ARRAY)
                                    {
                                      throw std::logic_error("findRefinedCellParamCoordsLinear logic err # 2");
                                    }
                                  unsigned fnk = cell_topo.getDimension()==2 ? k_node : face_nodes[k_node];
                                  for (unsigned ix=0; ix < 3; ix++)
                                    {
                                      param_coord[ix] += ref_topo_x[fnk].parametric_coordinates[ix]/4.0;
                                    }
                                }

                              break;
                            }
                        }
                    }
                }
            }

          if (!found)
            {
              found = true;

              for (unsigned ix=0; ix < 3; ix++)
                {
                  param_coord[ix] = 0.0;
                }
              unsigned nvert = FromTopology::vertex_count;
              double dnvert = (double)nvert;
              for (unsigned k_node = 0; k_node < nvert; k_node++)
                {
                  for (unsigned ix=0; ix < 3; ix++)
                    {
                      param_coord[ix] += ref_topo_x[k_node].parametric_coordinates[ix]/dnvert;
                    }
                }
            }

          for (unsigned ix=0; ix < 3; ix++)
            {
              ref_topo_x[childNodeIdx].parametric_coordinates[ix] = param_coord[ix];
            }
          if (0)
            std::cout<<"tmp param_coord in findRefinedCellParamCoordsLinear= " << childNodeIdx << " "
                     << param_coord[0] << " "
                     << param_coord[1] << " "
                     << param_coord[2] << std::endl;
        }

    }

    /// this is called one time (during code development) to generate and print a table of the extra refinement info

#define DEBUG_PRINT_REF_TOPO_X 1
    template<typename FromTopology,  typename ToTopology >
    void
    URP<FromTopology, ToTopology>::
    printRefinementTopoX_Table(std::ostream& out )
    {
      const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< ToTopology >();

      shards::CellTopology cell_topo(cell_topo_data);

      if (DEBUG_PRINT_REF_TOPO_X)
        {
          std::cout << "toTopoKey: " << toTopoKey << " topo_key_quad8      = " << topo_key_quad8 << " cell_topo= " << cell_topo.getName() << std::endl;
          std::cout << "toTopoKey: " << toTopoKey << " topo_key_shellquad8 = " << topo_key_shellquad8 << " cell_topo= " << cell_topo.getName() << std::endl;
        }

      unsigned n_edges = cell_topo_data->edge_count;
      unsigned n_faces = cell_topo.getFaceCount();
      if (n_faces == 0) n_faces = 1; // 2D face has one "face"
      unsigned n_sides = cell_topo.getSideCount();
      if (DEBUG_PRINT_REF_TOPO_X)  std::cout << "tmp  n_edges= " << n_edges << " n_faces= " << n_faces << " n_sides= " << n_sides << std::endl;

      Elem::CellTopology elem_celltopo = Elem::getCellTopology< FromTopology >();
      const Elem::RefinementTopology* ref_topo_p = Elem::getRefinementTopology(elem_celltopo);
      if (!ref_topo_p)
        throw std::runtime_error("printRefinementTopoX_Table:: error, no refinement topology found");
      const Elem::RefinementTopology& ref_topo = *ref_topo_p;

      unsigned num_child = ref_topo.num_child();
      unsigned num_child_nodes = ref_topo.num_child_nodes();

      if (DEBUG_PRINT_REF_TOPO_X) std::cout << "tmp num_child_nodes= " << num_child_nodes << " num_child= " << num_child << std::endl;

      typedef Elem::StdMeshObjTopologies::RefTopoX RefTopoX;

      RefTopoX_arr ref_topo_x = new Elem::StdMeshObjTopologies::RefinementTopologyExtraEntry[num_child_nodes];

      std::string ct_name = cell_topo.getName();
      Util::replace(ct_name, "_", "<");

      bool useIntrepid = true;
      if (useIntrepid)
      { 
        Kokkos::DynRankView<double,Kokkos::HostSpace> param_coord("param_coord",cell_topo.getDimension());
        for (unsigned iNode = 0; iNode < FromTopology::node_count; iNode++)
          {              
            Intrepid2::CellTools<Kokkos::HostSpace>::getReferenceNode(param_coord, cell_topo, iNode);
            ref_topo_x[iNode].parametric_coordinates[0] = param_coord(0);
            ref_topo_x[iNode].parametric_coordinates[1] = param_coord(1);
            ref_topo_x[iNode].parametric_coordinates[2] = param_coord(2);
            if (0) std::cout<<"tmp param_coord= "
                            << param_coord(0) << " "
                            << param_coord(1) << " "
                            << param_coord(2) << std::endl;
          }
      }

      findRefinedCellParamCoordsLinear(ref_topo, ref_topo_x);
      findRefinedCellParamCoords(ref_topo, ref_topo_x);

      out << "\n template<> RefTopoX RefinementTopologyExtra< shards:: " << ct_name << ">  > :: refinement_topology = {" << std::endl;
      for (unsigned childNodeIdx = 0; childNodeIdx < num_child_nodes; childNodeIdx++)
        {
          //*  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}

          unsigned rank_of_subcell            = 0;
          unsigned ordinal_of_subcell         = 0;
          unsigned ordinal_of_node_on_subcell = 0;
          unsigned num_node_on_subcell        = 0;

          findRefinedCellTopoInfo(childNodeIdx, ref_topo, &ref_topo_x[0], rank_of_subcell,
                                  ordinal_of_subcell, ordinal_of_node_on_subcell, num_node_on_subcell);

          double *param_coord = ref_topo_x[childNodeIdx].parametric_coordinates;
          out << "{\t" << childNodeIdx << ",\t" << rank_of_subcell << ",\t" << ordinal_of_subcell << ",\t"
              << ordinal_of_node_on_subcell << ",\t" << num_node_on_subcell
              << ",\t{" << param_coord[0] << ",\t" << param_coord[1] << ",\t" << param_coord[2] << "} }"
              << (childNodeIdx==(num_child_nodes-1)? " " : ",") << std::endl;
        }
      out << "\n};" << std::endl;
      delete[] ref_topo_x;
    }

  }

#endif
