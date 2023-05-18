// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// FIXME - cleanup all the comments with New ref info
/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

/**
 * NOTE: the text of the RefTopoX tables herein are generated from adapt/UniformRefinerPattern::printRefinementTopoX_Table
 *   This need only be done once at code development time, or if Intrepid2 changes their definitions of parametric coordinates.
 *   Note: these tables could be generated each time the code is run by using a 'bootstrap' method and by changing the types
 *     of RefTopoX to be a pointer to an array of RefinementTopologyExtra entries that can be allocated and filled on the fly
 *     in the bootstrap process.  However, we liked the idea of being able to view the tables, that's why they are generated
 *     and pasted in below.
 */

#include <stdexcept>

#include <Shards_BasicTopologies.hpp>

#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <adapt/sierra_element/MeshObjTopology.hpp>
#include <adapt/sierra_element/RefinementTopology.hpp>

#include <adapt/sierra_element/percept_code_types.hpp>

#include <percept/Util.hpp>


//using sierra::Diag::Tracespec;
//#define RuntimeError() std::runtime_exception

  using namespace percept;

  namespace percept {
    namespace Elem {

  namespace StdMeshObjTopologies {

    enum { NOT_ELEMENT = 0 , PARTICLE = 1, ROD = 2, SHELL = 3 };

    enum { EUA = END_UINT_ARRAY };
    typedef const Elem::MeshObjTopology * const_top_ptr ;

    static const_top_ptr node();
    static const_top_ptr point();
    static const_top_ptr line( UInt eclass, UInt nnode);
    static const_top_ptr tri(  UInt eclass, UInt nnode);
    static const_top_ptr tri4( UInt eclass);
    static const_top_ptr quad( UInt eclass, UInt nnode);
    static const_top_ptr tet(               UInt nnode);
    static const_top_ptr hex(               UInt nnode);
    static const_top_ptr wedge(             UInt nnode);
    static const_top_ptr pyramid(           UInt nnode);

    const MeshObjTopology *
    node()
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology node3D(shards::getCellTopologyData<shards::Node>());

      {
        static const_top_ptr node3D_child[] = { & node3D };

        static UInt   child_node[] = { 0 , EUA};
        static UInt * child_node_table[] = { child_node };

        static RefinementTopology node3D_refinement(&node3D, 1, node3D_child, 1, child_node_table, 0, NULL, 0, NULL, 0, NULL, NULL, false);
      }

      return &node3D;
    }


    const MeshObjTopology *
    point()
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology part3D(shards::getCellTopologyData<shards::Particle>());

      {
        static const_top_ptr part3D_child[] = { & part3D };

        static UInt   child_node[] = { 0 , EUA};
        static UInt * child_node_table[] = { child_node };

        static RefinementTopology part3D_refinement(&part3D, 1, part3D_child, 1, child_node_table, 0, NULL, 0, NULL, 0, NULL, NULL, false);
        //    static RefinementTopology part3D_refinement(shards::getCellTopologyData<shards::Particle>(), 1, part3D_child, 1, child_node_table, 0, NULL, 0, NULL, 0, NULL, NULL, false);
      }

      return &part3D;
    }


    /*--------------------------------------------------------------------*/
    /**
     *   0               1    PARENT Linear Edge Element Nodes (SPACE_DIM = 1!)
     *   o---------------o
     *
     *
     *   After refinement:
     *
     *
     *   0       2       1    CHILD Linear Edge Element Nodes (new nodes = *)
     *   o-------*-------o
     *
     *                      | CHILD Linear Edge Node Maps (global node numbers!)
     *   0       1          |
     *   o-------o          |
     *      E#1             | Element (or edge) 0: childNodeMap[0] = { 0, 2 };
     *                      |
     *           0       1  |
     *           o-------o  |
     *              E#2     | Element (or edge) 1: childNodeMap[1] = { 2, 1 };
     *
     *
     *  Refined Linear Edge (or Linear Bar element) PERMUTATION Node Maps:
     *
     *  Polarity = 1  { 0, 1; 2 }
     *  Polarity = 0  { 1, 0; 2 }
     *
     **/

    /*--------------------------------------------------------------------*/
    /**
     *   0       2       1    PARENT 3-Node Line Object Nodes
     *   o-------o-------o
     *
     *
     *   After refinement:
     *
     *
     *   0   3   2   4   1    CHILD Objects (new nodes = *)
     *   o---*---o---*---o
     *
     *
     *                      | CHILD Line Node Maps (global node numbers!)
     *   0   2   1          |
     *   o---o---o          |
     *      E#1             | Object  (or line) 0: childNodeMap[0] = { 0, 2, 3 };
     *                      |
     *           0   2   1  |
     *           o---o---o  |
     *              E#2     | Object  (or line) 1: childNodeMap[1] = { 2, 1, 4 };
     *
     *
     *  Refined 3-Node Line Object PERMUTATION Node Maps:
     *
     *  Polarity = 1  { 0, 1, 2; 3, 4 }
     *  Polarity = 0  { 1, 0, 2; 4, 3 }
     *
     **/

    /*  New ref topo info Line2
     *  ------------------
     *
     *  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
     */

#include "GeneratedRefinementTable.hpp"

    /*--------------------------------------------------------------------*/
    // Line topologies with 2 or 3 nodes
    //   eclass ==  NOT_ELEMENT   = >  edge
    //   eclass ==  ROD
    //   eclass ==  SHELL
    const MeshObjTopology *
    line(
         UInt                  eclass,
         UInt                  nnode)
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology edge_2(shards::getCellTopologyData<shards::Line<2> >());
      static MeshObjTopology edge_3(shards::getCellTopologyData<shards::Line<3> >());

      static MeshObjTopology rod_2(shards::getCellTopologyData<shards::Beam<2> >());
      static MeshObjTopology rod_3(shards::getCellTopologyData<shards::Beam<3> >());

      static MeshObjTopology shell_2(shards::getCellTopologyData<shards::ShellLine<2> >());
      static MeshObjTopology shell_3(shards::getCellTopologyData<shards::ShellLine<3> >());

      static bool first = true ;

      if ( first ) {

        first = false ;

        { // Initialize num-child and child-topology
          static const_top_ptr edge_2_child[] = { &edge_2 , &edge_2 };
          static const_top_ptr edge_3_child[] = { &edge_3 , &edge_3 };

          static const_top_ptr rod_2_child[] = { &rod_2 , &rod_2 };
          static const_top_ptr rod_3_child[] = { &rod_3 , &rod_3 };

          static const_top_ptr shell_2_child[] = { &shell_2 , &shell_2 };
          static const_top_ptr shell_3_child[] = { &shell_3 , &shell_3 };

          static const UInt   child_0[] = { 0 , 2 , 3 , EUA};
          static const UInt   child_1[] = { 2 , 1 , 4 , EUA};
          static const UInt * child_node_table[] = { child_0 , child_1 };

          // Initialize edge permutations
          static const UInt   perm_P1[] = { 0, 1, 2,  3, 4 , EUA}; // Vertices + children
          static const UInt   perm_P0[] = { 1, 0, 2,  4, 3 , EUA}; // Vertices + children
          static const UInt * perm_table[2] = { NULL , NULL }; // Polarity only

          if ( perm_table[0] ==  NULL ) {
            perm_table[ 0 ] = perm_P1 ;
            perm_table[ 1 ] = perm_P0 ;
          }

          // Initialize rod-edge and shell-edge

          static const UInt edge_0[] = { 0 , 1 , 2 ,  3 , 4 , EUA};
          static const UInt edge_1[] = { 1 , 0 , 2 ,  4 , 3 , EUA};
          static const UInt * edge_table[] = { edge_0 , edge_1 };

//           static RefinementTopology edge_2_refinement(&edge_2, 2, edge_2_child, 3, child_node_table, 0, NULL, 0, NULL, 2, perm_table, NULL, true);
//           static RefinementTopology edge_3_refinement(&edge_3, 2, edge_3_child, 5, child_node_table, 0, NULL, 0, NULL, 2, perm_table, NULL, true);

          static RefinementTopology edge_2_refinement(&edge_2, 2, edge_2_child, 3, child_node_table, 2, edge_table, 0, NULL, 2, perm_table, NULL, true);
          static RefinementTopology edge_3_refinement(&edge_3, 2, edge_3_child, 5, child_node_table, 2, edge_table, 0, NULL, 2, perm_table, NULL, true);

          static RefinementTopology  rod_2_refinement(&rod_2,  2, rod_2_child,  3, child_node_table, 2, edge_table, 0, NULL, 2, perm_table, NULL, true);
          static RefinementTopology  rod_3_refinement(&rod_3,  2, rod_3_child,  5, child_node_table, 2, edge_table, 0, NULL, 2, perm_table, NULL, true);

          static RefinementTopology shell_2_refinement(&shell_2, 2, shell_2_child, 3, child_node_table, 2, edge_table, 0, NULL, 2, perm_table, NULL, true);
          static RefinementTopology shell_3_refinement(&shell_3, 2, shell_3_child, 5, child_node_table, 2, edge_table, 0, NULL, 2, perm_table, NULL, true);
        }
      }

      MeshObjTopology * top = NULL ;

      switch( ( eclass << 8 ) | ( nnode << 4 ) ) {
      case 0x0020 : top = & edge_2 ; break ;
      case 0x0030 : top = & edge_3 ; break ;
      case 0x0220 : top = & rod_2 ; break ;
      case 0x0230 : top = & rod_3 ; break ;
      case 0x0320 : top = & shell_2 ; break ;
      case 0x0330 : top = & shell_3 ; break ;
      default :
        //throw RuntimeError() << "Invalid eclass and nnode specified" << std::endl;// << StackTrace;
        throw std::runtime_error( "Invalid eclass and nnode specified") ; // << std::endl;// << StackTrace;
      }

      return top ;
    }


    /*--------------------------------------------------------------------*/
    /**
     *  3                 2   PARENT Linear 4-Node Quadrilateral Element Nodes
     *   o---------------o    (SPACE_DIM = 2!)
     *   |               |
     *   |               |
     *   |               |
     *   |               |
     *   |               |    (PARENT) Linear 4-Node Quadrilateral
     *   |               |             Element Edge Node Map:
     *   |               |
     *   |               |    { {0, 1}, {1, 2}, {2, 3} {3, 0} };
     *   |               |
     *   o---------------o
     *  0                 1
     *
     *   After refinement:
     *
     *  3        6        2   CHILD Linear 4-Node Quadrilateral Element Nodes
     *   o-------*-------o    (SPACE_DIM = 2!) (new nodes = *)
     *   |       |       |
     *   |       |       |
     *   |       |       |
     *   |      8|       |
     *  7*-------*-------*5
     *   |       |       |
     *   |       |       |
     *   |       |       |
     *   |       |       |
     *   o-------*-------o
     *  0        4        1 | CHILD Linear 4-Node Quadrilateral Element Node Maps:
     *                      |
     *                      |
     *                      | Element 0: childNodeMap[0] = { 0, 4, 8, 7 }
     *                      | Element 1: childNodeMap[1] = { 4, 1, 5, 8 }
     *                      | Element 2: childNodeMap[2] = { 8, 5, 2, 6 }
     *                      | Element 3: childNodeMap[3] = { 7, 8, 6, 3 }
     *                      |
     *
     *  New ref topo info Quad4
     *  ------------------
     *
     *  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
     * {
     *  {0,   0,         0,    0,         1,             {0.0, 0.0, 0.0} },
     *  {1,   0,         1,    0,         1,             {1.0, 0.0, 0.0} },
     *  {2,   0,         2,    0,         1,             {1.0, 1.0, 0.0} },
     *  {3,   0,         3,    0,         1,             {0.0, 1.0, 0.0} },
     *  {4,   1,         0,    0,         1,             {0.5, 0.0, 0.0} },
     *  {5,   1,         1,    0,         1,             {1.0, 0.5, 0.0} },
     *  {6,   1,         2,    0,         1,             {0.5, 1.0, 0.0} },
     *  {7,   1,         3,    0,         1,             {0.0, 0.5, 0.0} },
     *  {8,   2,         0,    0,         1,             {0.5, 0.5, 0.0} }
     * }
     *
     *  Refined Linear 4-Node Quadrilateral Element PERMUTATION Node Maps:
     *
     *  Rotation  Polarity
     *     0          1      { 0, 1, 2, 3; 4, 5, 6, 7, 8 }
     *     0          0      { 0, 3, 2, 1; 7, 6, 5, 4, 8 }
     *     1          1      { 3, 0, 1, 2; 7, 4, 5, 6, 8 }
     *     1          0      { 3, 2, 1, 0; 6, 5, 4, 7, 8 }
     *     2          1      { 2, 3, 0, 1; 6, 7, 4, 5, 8 }
     *     2          0      { 2, 1, 0, 3; 5, 4, 7, 6, 8 }
     *     3          1      { 1, 2, 3, 0; 5, 6, 7, 4, 8 }
     *     3          0      { 1, 0, 3, 2; 4, 7, 6, 5, 8 }
     *
     */
    /*--------------------------------------------------------------------*/
    /**
     *  3        6        2   PARENT 9-Node Quadrilateral Object Nodes
     *   o-------o-------o
     *   |               |
     *   |               |
     *   |       8       |
     * 7 o       o       o 5  (PARENT) 9-Node Quadrilateral Object's
     *   |               |             Edge Node Map:
     *   |               |
     *   |               |    { {0, 1, 4}, {1, 2, 5}, {2, 3, 6} {3, 0, 7} };
     *   o-------o-------o
     *  0        4        1
     *
     *
     *   After refinement:
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
     *   CHILD 9-Node Quadrilateral Object Node Maps:
     *                       |
     *                       |
     *                       | Object  0: childNodeMap[0] = { 0, 4, 8, 7, 9, 17, 20, 16;  21 }
     *                       | Object  1: childNodeMap[1] = { 4, 1, 5, 8, 10, 11, 18, 17; 22 }
     *                       | Object  2: childNodeMap[2] = { 8, 5, 2, 6, 18, 12, 13, 19; 23 }
     *                       | Object  3: childNodeMap[3] = { 7, 8, 6, 3, 20, 19, 14, 15; 24 }
     *                       |
     *  New ref topo info Quad9
     *  ------------------
     *
     *  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
     * {
     *  {0,   0,         0,      0,       1,             {0.0, 0.0, 0.0} },
     *  {1,   0,         1,      0,       1,             {1.0, 0.0, 0.0} },
     *  {2,   0,         2,      0,       1,             {1.0, 1.0, 0.0} },
     *  {3,   0,         3,      0,       1,             {0.0, 1.0, 0.0} },
     *  {4,   1,         0,      0,       1,             {0.5, 0.0, 0.0} },
     *  {5,   1,         1,      0,       1,             {1.0, 0.5, 0.0} },
     *  {6,   1,         2,      0,       1,             {0.5, 1.0, 0.0} },
     *  {7,   1,         3,      0,       1,             {0.0, 0.5, 0.0} },
     *  {8,   2,         0,      8,       9,             {0.5, 0.5, 0.0} },
     *
     *  {9,   1,         0,      1,       3,             {0.25, 0.00, 0.00} },
     *  {10,  1,         0,      2,       3,             {0.75, 0.00, 0.00} },
     *  {11,  1,         1,      1,       3,             {1.00, 0.25, 0.00} },
     *  {12,  1,         1,      2,       3,             {1.00, 0.75, 0.00} },
     *  {13,  1,         2,      1,       3,             {0.75, 1.00, 0.00} },
     *  {14,  1,         2,      2,       3,             {0.25, 1.00, 0.00} },
     *  {15,  1,         3,      1,       3,             {0.00, 0.75, 0.00} },
     *  {16,  1,         3,      2,       3,             {0.00, 0.25, 0.00} }

     *  {17,  2,         0,      4,       9,             {0.50, 0.25, 0.00} },
     *  {18,  2,         0,      5,       9,             {0.75, 0.50, 0.00} },
     *  {19,  2,         0,      6,       9,             {0.50, 0.75, 0.00} },
     *  {20,  2,         0,      7,       9,             {0.25, 0.50, 0.00} },
     *
     *  {21,  2,         0,      0,       9,             {0.25, 0.25, 0.00} },
     *  {22,  2,         0,      1,       9,             {0.75, 0.25, 0.00} },
     *  {23,  2,         0,      2,       9,             {0.75, 0.75, 0.00} },
     *  {24,  2,         0,      3,       9,             {0.25, 0.75, 0.00} },
     *
     *
     *
     * }
     *
     *  Refined 9-Node Quadrilateral Object PERMUTATION Node Maps:
     *
     *  Rotation  Polarity
     *     0          1      { 0, 1, 2, 3, 4, 5, 6, 7; 8, 9, 10, 11, 12, 13, 14, 15, 16,
     *                         17, 18, 19, 20, 21, 22, 23, 24 }
     *     0          0      { 0, 3, 2, 1, 7, 6, 5, 4; 8, 16, 15, 14, 13, 12, 11, 10, 9,
     *                         20, 19, 18, 17, 21, 24, 23, 22 }
     *     1          1      { 3, 0, 1, 2, 7, 4, 5, 6; 8, 15, 16, 9, 10, 11, 12, 13, 14,
     *                         20, 17, 18, 19, 24, 21, 22, 23 }
     *     1          0      { 3, 2, 1, 0, 6, 5, 4, 7; 8, 14, 13, 12, 11, 10, 9, 16, 15,
     *                         19, 18, 17, 20, 24, 23, 22, 21 }
     *     1          1      { 1, 2, 3, 0, 5, 6, 7, 4; 8, 11, 12, 13, 14, 15, 16, 9, 10,
     *                         18, 19, 20, 17, 22, 23, 24, 21 }
     *     1          0      { 1, 0, 3, 2, 4, 7, 6, 5; 8, 10, 9, 16, 15, 14, 13, 12, 11,
     *                         17, 20, 19, 18, 22, 21, 24, 23 }
     *     3          1      { 2, 3, 0, 1, 6, 7, 4, 5; 8, 13, 14, 15, 16, 9, 10, 11, 12,
     *                         19, 20, 17, 18, 23, 24, 21, 22 }
     *     3          0      { 2, 1, 0, 3, 5, 4, 7, 6; 8, 12, 11, 10, 9, 16, 15, 14, 13,
     *                         18, 17, 20, 19, 23, 22, 21, 24 }
     *
     **/
    /*--------------------------------------------------------------------*/
    // Quadrilateral with 4, 8, or 9 nodes.  Face, shell, or 2D solid.
    //   eclass ==  NOT_ELEMENT  = > face in 3D
    //   eclass ==  SHELL        = > element in 3D
    //   eclass ==  SOLID        = > element in 2D

    const MeshObjTopology *
    quad(UInt eclass, UInt nnode)
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology face_4(shards::getCellTopologyData<shards::Quadrilateral<4> >());
      static MeshObjTopology face_8(shards::getCellTopologyData<shards::Quadrilateral<8> >());
      static MeshObjTopology face_9(shards::getCellTopologyData<shards::Quadrilateral<9> >());

      static MeshObjTopology shell_4(shards::getCellTopologyData<shards::ShellQuadrilateral<4> >());
      static MeshObjTopology shell_8(shards::getCellTopologyData<shards::ShellQuadrilateral<8> >());
      static MeshObjTopology shell_9(shards::getCellTopologyData<shards::ShellQuadrilateral<9> >());

      static bool first = true ;

      if ( first ) {

        first = false ;

        { // Child topologies and child nodes

          static const_top_ptr
            face_4_child[] = { &face_4, &face_4, &face_4, &face_4 };
          static const_top_ptr
            face_8_child[] = { &face_8, &face_8, &face_8, &face_8 };
          static const_top_ptr
            face_9_child[] = { &face_9, &face_9, &face_9, &face_9 };

          static const_top_ptr
            shell_4_child[] = {&shell_4, &shell_4, &shell_4, &shell_4};
          static const_top_ptr
            shell_8_child[] = {&shell_8, &shell_8, &shell_8, &shell_8};
          static const_top_ptr
            shell_9_child[] = {&shell_9, &shell_9, &shell_9, &shell_9};

          static const UInt child_0[] = { 0, 4, 8, 7,  9, 17, 20, 16, 21 , EUA};
          static const UInt child_1[] = { 4, 1, 5, 8, 10, 11, 18, 17, 22 , EUA};
          static const UInt child_2[] = { 8, 5, 2, 6, 18, 12, 13, 19, 23 , EUA};
          static const UInt child_3[] = { 7, 8, 6, 3, 20, 19, 14, 15, 24 , EUA};

          static const UInt * child_node_table[] = {child_0, child_1, child_2, child_3};

          // Face permutations
          // Permutation tables including children [ 2 * number_of_vertices ]
          // Five groups of nodes:
          // a) vertices
          // b) outer edge mid-points + centroid
          // c) outer edge quarter points
          // d) inner edge mid-points
          // e) child centroid

          static const UInt perm_P1_R0[] = {  0,    1,    2,    3,
                                              4,    5,    6,    7,    8,
                                              9, 10, 11, 12, 13, 14, 15, 16,
                                              17,   18,   19,   20,
                                              21,   22,   23,   24 , EUA};

          static const UInt perm_P1_R1[] = { 3,    0,    1,    2,
                                             7,    4,    5,    6,    8,
                                             15, 16, 9, 10, 11, 12, 13, 14,
                                             20,   17,   18,   19,
                                             24,   21,   22,   23 , EUA};

          static const UInt perm_P1_R2[] = { 2,    3,    0,    1,
                                             6,    7,    4,    5,    8,
                                             13, 14, 15, 16, 9, 10, 11, 12,
                                             19,   20,   17,   18,
                                             23,   24,   21,   22 , EUA};

          static const UInt perm_P1_R3[] = { 1,    2,    3,    0,
                                             5,    6,    7,    4,    8,
                                             11, 12, 13, 14, 15, 16, 9, 10,
                                             18,   19,   20,   17,
                                             22,   23,   24,   21 , EUA};

          static const UInt perm_P0_R0[] = {  0,    3,    2,    1,
                                              7,    6,    5,    4,    8,
                                              16, 15, 14, 13, 12, 11, 10, 9,
                                              20,   19,   18,   17,
                                              21,   24,   23,   22 , EUA};

          static const UInt perm_P0_R1[] = { 3,    2,    1,    0,
                                             6,    5,    4,    7,    8,
                                             14, 13, 12, 11, 10, 9, 16, 15,
                                             19,   18,   17,   20,
                                             24,   23,   22,   21 , EUA};

          static const UInt perm_P0_R2[] = { 2,    1,    0,    3,
                                             5,    4,    7,    6,    8,
                                             12, 11, 10, 9, 16, 15, 14, 13,
                                             18,   17,   20,   19,
                                             23,   22,   21,   24 , EUA};

          static const UInt perm_P0_R3[] = { 1,    0,    3,    2,
                                             4,    7,    6,    5,    8,
                                             10, 9, 16, 15, 14, 13, 12, 11,
                                             17,   20,   19,   18,
                                             22,   21,   24,   23 , EUA};


          static const UInt * perm_table[] = { NULL , NULL , NULL , NULL , NULL , NULL , NULL , NULL };

          perm_table[ 0 ] = perm_P1_R0 ;
          perm_table[ 1 ] = perm_P1_R1 ;
          perm_table[ 2 ] = perm_P1_R2 ;
          perm_table[ 3 ] = perm_P1_R3 ;
          perm_table[ 4 ] = perm_P0_R0 ;
          perm_table[ 5 ] = perm_P0_R1 ;
          perm_table[ 6 ] = perm_P0_R2 ;
          perm_table[ 7 ] = perm_P0_R3 ;

          // Face edge permutations
          // Permutation tables [ 2 * number_of_vertices ]

          static const UInt edge_perm_P1_R0[] = {  0, 1, 2, 3 , EUA};
          static const UInt edge_perm_P1_R1[] = {  3, 0, 1, 2 , EUA};
          static const UInt edge_perm_P1_R2[] = {  2, 3, 0, 1 , EUA};
          static const UInt edge_perm_P1_R3[] = {  1, 2, 3, 0 , EUA};
          static const UInt edge_perm_P0_R0[] = {  3, 2, 1, 0 , EUA};
          static const UInt edge_perm_P0_R1[] = {  2, 1, 0, 3 , EUA};
          static const UInt edge_perm_P0_R2[] = {  1, 0, 3, 2 , EUA};
          static const UInt edge_perm_P0_R3[] = {  0, 3, 2, 1 , EUA};

          static const UInt * edge_perm_table[] = { NULL , NULL , NULL , NULL ,
                                                    NULL , NULL , NULL , NULL };

          edge_perm_table[ 0 ] = edge_perm_P1_R0 ;
          edge_perm_table[ 1 ] = edge_perm_P1_R1 ;
          edge_perm_table[ 2 ] = edge_perm_P1_R2 ;
          edge_perm_table[ 3 ] = edge_perm_P1_R3 ;
          edge_perm_table[ 4 ] = edge_perm_P0_R0 ;
          edge_perm_table[ 5 ] = edge_perm_P0_R1 ;
          edge_perm_table[ 6 ] = edge_perm_P0_R2 ;
          edge_perm_table[ 7 ] = edge_perm_P0_R3 ;

          // Edge topology and node tables including edges' child-nodes
          static const UInt   edge_0[] = { 0, 1,  4,   9, 10 , EUA};
          static const UInt   edge_1[] = { 1, 2,  5,  11, 12 , EUA};
          static const UInt   edge_2[] = { 2, 3,  6,  13, 14 , EUA};
          static const UInt   edge_3[] = { 3, 0,  7,  15, 16 , EUA};
          static const UInt * edge_table[] = { edge_0 , edge_1 , edge_2 , edge_3 };

          // 3D Shell has two faces

          static const UInt face_0[] = {  0,    1,    2,    3,
                                          4,    5,    6,    7,    8,
                                          9, 10, 11, 12, 13, 14, 15, 16,
                                          17,   18,   19,   20,
                                          21,   22,   23,   24 , EUA};

          static const UInt face_1[] = {  0,    3,    2,    1,
                                          7,    6,    5,    4,    8,
                                          16, 15, 14, 13, 12, 11, 10, 9,
                                          20,   19,   18,   17,
                                          21,   24,   23,   22 , EUA};

          static const UInt * face_table[] = { face_0 , face_1 };

          static RefinementTopology face_4_refinement(&face_4, 4, face_4_child, 9, child_node_table, 4, edge_table, 0, NULL, 8, perm_table, edge_perm_table, true);
          static RefinementTopology face_8_refinement(&face_8, 4, face_8_child, 21, child_node_table, 4, edge_table, 0, NULL, 8, perm_table, edge_perm_table, true);
          static RefinementTopology face_9_refinement(&face_9, 4, face_9_child, 25, child_node_table, 4, edge_table, 0, NULL, 8, perm_table, edge_perm_table, true);

          static RefinementTopology shell_4_refinement(&shell_4, 4, shell_4_child, 9, child_node_table, 4, edge_table, 2, face_table, 0, NULL, NULL, true);
          static RefinementTopology shell_8_refinement(&shell_8, 4, shell_8_child, 21, child_node_table, 4, edge_table, 2, face_table, 0, NULL, NULL, true);
          static RefinementTopology shell_9_refinement(&shell_9, 4, shell_9_child, 25, child_node_table, 4, edge_table, 2, face_table, 0, NULL, NULL, true);
        }
      }

      //--------------------------------------------------------------------

      MeshObjTopology * top = NULL ;

      switch( ( eclass << 8 ) | ( nnode << 4 ) ) {
      case 0x0040 : top = &face_4 ; break ;
      case 0x0080 : top = &face_8 ; break ;
      case 0x0090 : top = &face_9 ; break ;
      case 0x0340 : top = &shell_4 ; break ;
      case 0x0380 : top = &shell_8 ; break ;
      case 0x0390 : top = &shell_9 ; break ;
      default:
        //throw RuntimeError() << "Invalid eclass and nnode specified" << std::endl;// << StackTrace;
        throw std::runtime_error( "Invalid eclass and nnode specified") ; // << std::endl;// << StackTrace;
      }

      return top ;
    }

    /*  New ref topo info Quad4
     *  ------------------
     *
     *  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
     */
#if 0
    template<>
    RefTopoX
    RefinementTopologyExtra< shards::Quadrilateral<4> >::refinement_topology =
      {
        {0,   0,         0,    0,         1,             {0.0, 0.0, 0.0} },
        {1,   0,         1,    0,         1,             {1.0, 0.0, 0.0} },
        {2,   0,         2,    0,         1,             {1.0, 1.0, 0.0} },
        {3,   0,         3,    0,         1,             {0.0, 1.0, 0.0} },
        {4,   1,         0,    0,         1,             {0.5, 0.0, 0.0} },
        {5,   1,         1,    0,         1,             {1.0, 0.5, 0.0} },
        {6,   1,         2,    0,         1,             {0.5, 1.0, 0.0} },
        {7,   1,         3,    0,         1,             {0.0, 0.5, 0.0} },
        {8,   2,         0,    0,         1,             {0.5, 0.5, 0.0} }
      };
#endif
    /* New ref topo info Quad9
     * ------------------
     *
     * {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
     */

#if 0
    template<> RefTopoX RefinementTopologyExtra< shards:: Quadrilateral<9>  > :: refinement_topology = {
      {	0,	0,	0,	0,	1,	{-1,	-1,	0} },
      {	1,	0,	1,	0,	1,	{1,	-1,	0} },
      {	2,	0,	2,	0,	1,	{1,	1,	0} },
      {	3,	0,	3,	0,	1,	{-1,	1,	0} },
      {	4,	1,	0,	0,	3,	{0,	-1,	0} },
      {	5,	1,	1,	0,	3,	{1,	0,	0} },
      {	6,	1,	2,	0,	3,	{0,	1,	0} },
      {	7,	1,	3,	0,	3,	{-1,	0,	0} },
      {	8,	2,	0,	0,	9,	{0,	0,	0} },
      {	9,	1,	0,	1,	3,	{-0.5,	-1,	0} },
      {	10,	1,	0,	2,	3,	{0.5,	-1,	0} },
      {	11,	1,	1,	1,	3,	{1,	-0.5,	0} },
      {	12,	1,	1,	2,	3,	{1,	0.5,	0} },
      {	13,	1,	2,	1,	3,	{0.5,	1,	0} },
      {	14,	1,	2,	2,	3,	{-0.5,	1,	0} },
      {	15,	1,	3,	1,	3,	{-1,	0.5,	0} },
      {	16,	1,	3,	2,	3,	{-1,	-0.5,	0} },
      {	17,	2,	0,	1,	9,	{0,	-0.5,	0} },
      {	18,	2,	0,	2,	9,	{0.5,	0,	0} },
      {	19,	2,	0,	3,	9,	{0,	0.5,	0} },
      {	20,	2,	0,	4,	9,	{-0.5,	0,	0} },
      {	21,	2,	0,	5,	9,	{-0.5,	-0.5,	0} },
      {	22,	2,	0,	6,	9,	{0.5,	-0.5,	0} },
      {	23,	2,	0,	7,	9,	{0.5,	0.5,	0} },
      {	24,	2,	0,	8,	9,	{-0.5,	0.5,	0} }

    };
#endif

    /*--------------------------------------------------------------------*/
    /**
     *           2            PARENT Linear 3-Node Triangle Element Nodes
     *           o            (SPACE_DIM = 2!)
     *          / \
     *         /   \          (PARENT) Linear 3-Node Triangle Edge Node Map:
     *        /     \
     *       /       \        { {0, 1}, {1, 2}, {2, 0} };
     *      /         \
     *     /           \
     *    /             \
     *   o---------------o
     *  0                 1
     *
     *   After refinement:
     *
     *           2            CHILD Linear 3-Node Triangle Element Nodes
     *           o                  (new nodes = *)
     *          / \
     *         /   \
     *        /     \
     *     5 *-------* 4
     *      / \     / \
     *     /   \   /   \
     *    /     \ /     \
     *   o-------*-------o
     *  0        3        1
     *
     * | CHILD Linear 3-Node Triangle Element Node Maps:
     * |
     * |
     * | static const UInt child_0[] = { 0, 3, 5 };
     * | static const UInt child_1[] = { 3, 1, 4 };
     * | static const UInt child_2[] = { 5, 4, 2 };
     * | static const UInt child_3[] = { 4, 5, 3 };
     * |
     */

    /*  New ref topo info Tri3
     *  ------------------
     *
     *  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
     */
#if 0
    template<>
    RefTopoX
    RefinementTopologyExtra< shards::Triangle<3> >::refinement_topology =
      {
        {0,   0,         0,    0,         1,             {0.0, 0.0, 0.0} },
        {1,   0,         1,    0,         1,             {1.0, 0.0, 0.0} },
        {2,   0,         2,    0,         1,             {0.0, 1.0, 0.0} },

        {3,   1,         0,    0,         1,             {0.5, 0.0, 0.0} },
        {4,   1,         1,    0,         1,             {0.5, 0.5, 0.0} },
        {5,   1,         2,    0,         1,             {0.0, 0.5, 0.0} },

      };
#endif

    /*  Refined Linear 3-Node Triangle Element PERMUTATION Node Maps:
     *
     *  Rotation  Polarity
     *     0          1       { 0, 1, 2; 3, 4, 5 }
     *     0          0       { 0, 2, 1; 5, 4, 3 }
     *     1          1       { 2, 0, 1; 5, 3, 4 }
     *     1          0       { 2, 1, 0; 4, 3, 5 }
     *     2          1       { 1, 2, 0; 4, 5, 3 }
     *     2          0       { 1, 0, 2; 3, 5, 4 }
     *
     */

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

    /*--------------------------------------------------------------------*/
    // Triangle with 3 or 6 nodes.  Face, shell, or 2D solid.
    //   eclass ==  NOT_ELEMENT  = > sdim ==  3 && mdim ==  2
    //   eclass ==  SOLID        = > sdim ==  2 && mdim ==  2
    //   eclass ==  SHELL        = > sdim ==  3 && mdim ==  3
    const MeshObjTopology *
    tri(UInt eclass, UInt nnode)
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology face_3(shards::getCellTopologyData<shards::Triangle<3> >());
      static MeshObjTopology face_6(shards::getCellTopologyData<shards::Triangle<6> >());
      static MeshObjTopology shell_3(shards::getCellTopologyData<shards::ShellTriangle<3> >());
      static MeshObjTopology shell_6(shards::getCellTopologyData<shards::ShellTriangle<6> >());

      static bool first = true ;

      if ( first ) {

        first = false ;

        { // Child topologies and child nodes

          static const_top_ptr
            face_3_child[] = { &face_3, &face_3, &face_3, &face_3 };
          static const_top_ptr
            face_6_child[] = { &face_6, &face_6, &face_6, &face_6 };

          static const_top_ptr
            shell_3_child[] = {&shell_3, &shell_3, &shell_3, &shell_3};
          static const_top_ptr
            shell_6_child[] = {&shell_6, &shell_6, &shell_6, &shell_6};

          static const UInt child_0[] = { 0, 3, 5,  6, 12, 11 , EUA};
          static const UInt child_1[] = { 3, 1, 4,  7,  8, 13 , EUA};
          static const UInt child_2[] = { 5, 4, 2, 14,  9, 10 , EUA};
          static const UInt child_3[] = { 4, 5, 3, 14, 12, 13 , EUA};

          static const UInt * child_node_table[4]  =
            { child_0, child_1, child_2, child_3 };


          // Permutation tables including children [ 2 * number_of_vertices ]

          static const UInt perm_P1_R0[] = {  0, 1, 2,
                                              3, 4, 5,
                                              6, 7, 8, 9, 10, 11,
                                              12, 13, 14 , EUA};

          static const UInt perm_P1_R1[] = {  2, 0, 1,
                                              5, 3, 4,
                                              10, 11, 6, 7, 8, 9,
                                              14, 12, 13 , EUA};

          static const UInt perm_P1_R2[] = {  1, 2, 0,
                                              4, 5, 3,
                                              8, 9, 10, 11, 6, 7,
                                              13, 14, 12 , EUA};
          static const UInt perm_P0_R0[] = {  0, 2, 1,
                                              5, 4, 3,
                                              11, 10, 9, 8, 7, 6,
                                              12, 14, 13 , EUA};

          static const UInt perm_P0_R1[] = {  2, 1, 0,
                                              4, 3, 5,
                                              9, 8, 7, 6, 11, 10,
                                              14, 13, 12 , EUA};

          static const UInt perm_P0_R2[] = {  1, 0, 2,
                                              3, 5, 4,
                                              7, 6, 11, 10, 9, 8,
                                              13, 12, 14 , EUA};


          static const UInt * perm_table[] = { NULL, NULL, NULL, NULL, NULL, NULL };

          perm_table[ 0 ] = perm_P1_R0 ;
          perm_table[ 1 ] = perm_P1_R1 ;
          perm_table[ 2 ] = perm_P1_R2 ;
          perm_table[ 3 ] = perm_P0_R0 ;
          perm_table[ 4 ] = perm_P0_R1 ;
          perm_table[ 5 ] = perm_P0_R2 ;

          // Edge permutation tables [ 2 * number_of_vertices ]

          static const UInt edge_perm_P1_R0[] = { 0, 1, 2 , EUA};
          static const UInt edge_perm_P1_R1[] = { 2, 0, 1 , EUA};
          static const UInt edge_perm_P1_R2[] = { 1, 2, 0 , EUA};
          static const UInt edge_perm_P0_R0[] = { 2, 1, 0 , EUA};
          static const UInt edge_perm_P0_R1[] = { 1, 0, 2 , EUA};
          static const UInt edge_perm_P0_R2[] = { 0, 2, 1 , EUA};

          static const UInt * edge_perm_table[]  = { NULL, NULL, NULL, NULL, NULL, NULL };

          edge_perm_table[ 0 ] = edge_perm_P1_R0 ;
          edge_perm_table[ 1 ] = edge_perm_P1_R1 ;
          edge_perm_table[ 2 ] = edge_perm_P1_R2 ;
          edge_perm_table[ 3 ] = edge_perm_P0_R0 ;
          edge_perm_table[ 4 ] = edge_perm_P0_R1 ;
          edge_perm_table[ 5 ] = edge_perm_P0_R2 ;

          // Edge topology and node tables including edges' child-nodes
          static const UInt   edge_0[] = { 0, 1,  3,   6,  7 , EUA};
          static const UInt   edge_1[] = { 1, 2,  4,   8,  9 , EUA};
          static const UInt   edge_2[] = { 2, 0,  5,  10, 11 , EUA};
          static const UInt * edge_table[3] = { edge_0 , edge_1 , edge_2 };

          // 3D shells have two faces
          static const UInt face_0[] = {  0,  1,  2,
                                          3,  4,  5,
                                          6, 7, 8, 9, 10, 11,
                                          12, 13, 14 , EUA};

          static const UInt face_1[] = {  0,  2,  1,
                                          5,  4,  3,
                                          11, 10, 9, 8, 7, 6,
                                          12, 14, 13, EUA};

          static const UInt * face_table[] = { face_0 , face_1 };

          static RefinementTopology face_3_refinement(&face_3, 4, face_3_child, 6, child_node_table, 3, edge_table, 0, NULL, 6, perm_table, edge_perm_table, true);
          static RefinementTopology face_6_refinement(&face_6, 4, face_6_child, 15, child_node_table, 3, edge_table, 0, NULL, 6, perm_table, edge_perm_table, true);
          static RefinementTopology shell_3_refinement(&shell_3, 4, shell_3_child, 6, child_node_table, 3, edge_table, 2, face_table, 0, NULL, NULL, true);
          static RefinementTopology shell_6_refinement(&shell_6, 4, shell_6_child, 15, child_node_table, 3, edge_table, 2, face_table, 0, NULL, NULL, true);
        }
      }

      MeshObjTopology * top = NULL ;

      switch( ( eclass << 8 ) | ( nnode << 4 ) ) {
      case 0x0030 : top = &face_3 ; break ;
      case 0x0060 : top = &face_6 ; break ;
      case 0x0330 : top = &shell_3 ; break ;
      case 0x0360 : top = &shell_6 ; break ;
      case 0x0430 : top = &face_3 ; break ;
      case 0x0460 : top = &face_6 ; break ;
      default:
        //throw RuntimeError() << "Invalid eclass and nnode specified" << std::endl ;// << StackTrace;
        throw std::runtime_error( "Invalid eclass and nnode specified") ; // << std::endl ;// << StackTrace;
      }

      return top ;
    }


    /*-------------------------------------------------------------------*/
    /**
     *           2            PARENT 4-Node Triangle Object Nodes
     *           o
     *          / \
     *         /   \          (PARENT) 4-Node Triangle Object Edge Node Map:
     *        /     \
     *       /       \        { {0, 1}, {1, 2}, {2, 0} };
     *      /    o    \
     *     /     3     \
     *    /             \
     *   o---------------o
     *  0                 1
     *
     *
     *   After refinement:
     *
     *           2            CHILD 4-Node Triangle Object Nodes
     *           o                  (new nodes = *)
     *          / \
     *         / 9 \
     *        /  *  \
     *     6 *-------* 5
     *      / \  o  / \
     *     / 7 \ 3 / 8 \
     *    /  *  \ /  *  \
     *   o-------*-------o
     *  0        4        1
     *
     * | CHILD 4-Node Triangle Object Node Maps:
     * |
     * | static const UInt child_0[] = { 0, 4, 6, 7 };
     * | static const UInt child_1[] = { 4, 1, 5, 8 };
     * | static const UInt child_2[] = { 6, 5, 2, 9 };
     * | static const UInt child_3[] = { 5, 6, 4, 3 };
     * |
     *
     *  4-Node Triangle Object PERMUTATION Node Maps:
     *
     *  Original Parent 4-Node Triangle Object:
     *  Rotation  Polarity
     *     0          1       { 0, 1, 2, 3 }
     *     0          0       { 0, 2, 1, 3 }
     *     1          1       { 2, 0, 1, 3 }
     *     1          0       { 2, 1, 0, 3 }
     *     2          1       { 1, 2, 0, 3 }
     *     2          0       { 1, 0, 2, 3 }
     *
     *  After Refinement and using child node numbering:
     *  Rotation  Polarity
     *
     *     0          1       { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
     *     0          0       { 0, 2, 1, 3, 6, 5, 4, 7, 9, 8 };
     *     1          1       { 2, 0, 1, 3, 6, 4, 5, 9, 7, 8 };
     *     1          0       { 2, 1, 0, 3, 6, 5, 4, 9, 8, 7 };
     *     2          1       { 1, 2, 0, 3, 5, 6, 4, 8, 9, 7 };
     *     2          0       { 1, 0, 2, 3, 4, 6, 5, 8, 7, 9 };
     *
     */

    /*-------------------------------------------------------------------*/
    // Triangle with 4 nodes.  Face, or 2D solid.
    //   eclass ==  NOT_ELEMENT
    //   eclass ==  SOLID
    const MeshObjTopology *
    tri4( UInt /* eclass */)
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology face_4(shards::getCellTopologyData<shards::Triangle<4> >());

      static bool first = true ;

      if ( first ) {

        first = false ;

        { // Child topologies and child nodes

          static const_top_ptr face_4_child[]  =
            { &face_4 , &face_4 , &face_4 , &face_4 };

          static const UInt child_0[] = { 0, 4, 6, 7 , EUA};
          static const UInt child_1[] = { 4, 1, 5, 8 , EUA};
          static const UInt child_2[] = { 6, 5, 2, 9 , EUA};
          static const UInt child_3[] = { 5, 6, 4, 3 , EUA};

          static const UInt * child_node_table[] = {child_0, child_1, child_2, child_3};

          // Permutation tables including children [ 2 * number_of_vertices ]

          static const UInt perm_P1_R0[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 , EUA};
          static const UInt perm_P1_R1[] = { 2, 0, 1, 3, 6, 4, 5, 9, 7, 8 , EUA};
          static const UInt perm_P1_R2[] = { 1, 2, 0, 3, 5, 6, 4, 8, 9, 7 , EUA};
          static const UInt perm_P0_R0[] = { 0, 2, 1, 3, 6, 5, 4, 7, 9, 8 , EUA};
          static const UInt perm_P0_R1[] = { 2, 1, 0, 3, 5, 4, 6, 9, 8, 7 , EUA};
          static const UInt perm_P0_R2[] = { 1, 0, 2, 3, 4, 6, 5, 8, 7, 9 , EUA};

          static const UInt * perm_table[] = { NULL, NULL, NULL, NULL, NULL, NULL };

          perm_table[ 0 ] = perm_P1_R0 ;
          perm_table[ 1 ] = perm_P1_R1 ;
          perm_table[ 2 ] = perm_P1_R2 ;
          perm_table[ 3 ] = perm_P0_R0 ;
          perm_table[ 4 ] = perm_P0_R1 ;
          perm_table[ 5 ] = perm_P0_R2 ;

          // Edge permutation tables [ 2 * number_of_vertices ]

          static const UInt edge_perm_P1_R0[] = { 0, 1, 2 , EUA};
          static const UInt edge_perm_P1_R1[] = { 2, 0, 1 , EUA};
          static const UInt edge_perm_P1_R2[] = { 1, 2, 0 , EUA};
          static const UInt edge_perm_P0_R0[] = { 2, 1, 0 , EUA};
          static const UInt edge_perm_P0_R1[] = { 1, 0, 2 , EUA};
          static const UInt edge_perm_P0_R2[] = { 0, 2, 1 , EUA};

          static const UInt * edge_perm_table[]  = { NULL, NULL, NULL, NULL, NULL, NULL };

          edge_perm_table[ 0 ] = edge_perm_P1_R0 ;
          edge_perm_table[ 1 ] = edge_perm_P1_R1 ;
          edge_perm_table[ 2 ] = edge_perm_P1_R2 ;
          edge_perm_table[ 3 ] = edge_perm_P0_R0 ;
          edge_perm_table[ 4 ] = edge_perm_P0_R1 ;
          edge_perm_table[ 5 ] = edge_perm_P0_R2 ;

          // Edge topology and node tables including edges' child-nodes

          static const UInt   edge_0[] = { 0, 1,  4 , EUA};
          static const UInt   edge_1[] = { 1, 2,  5 , EUA};
          static const UInt   edge_2[] = { 2, 0,  6 , EUA};
          static const UInt * edge_table[] = { edge_0 , edge_1 , edge_2 };

          static RefinementTopology face_4_refinement(&face_4, 4, face_4_child, 10, child_node_table, 3, edge_table, 0, NULL, 6, perm_table, edge_perm_table, false);
        }
      }

      MeshObjTopology * const top = &face_4 ;

      return top ;
    }

    /*--------------------------------------------------------------------*/
    /**
     *                                  PARENT Linear 8-Node Hexahedron Nodes
     *         7                    6   (SPACE_DIM = 3!)
     *          o------------------o
     *         /|                 /|
     *        / |                / |
     *       /  |               /  |
     *      /   |              /   |
     *     /    |             /    |
     *    /     |            /     |
     * 4 /      |         5 /      |
     *  o------------------o       |
     *  |       |          |       |
     *  |     3 o----------|-------o 2
     *  |      /           |      /
     *  |     /            |     /
     *  |    /             |    /
     *  |   /              |   /
     *  |  /               |  /
     *  | /                | /
     *  |/                 |/
     *  o------------------o
     * 0                    1
     *
     *                        (PARENT) Linear 8-Node Hexahedron
     *                                 3D Element Edge Node Map:
     *
     *                              { {0, 1}, {1, 2}, {2, 3}, {3, 0},
     *                                {4, 5}, {5, 6}, {6, 7}, {7, 4},
     *                                {0, 4}, {1, 5}, {2, 6}, {3, 7}  };
     *
     *                                 3D Element Face Node Map:
     *
     *                              { {0, 1, 5, 4}, {1, 2, 6, 5}, { 2, 3, 7, 6}, { 0, 4, 7, 3},  { 0, 3, 2, 1}, { 4, 5, 6, 7} };
     * Shards face list info:
     *
     *  typedef
     *   MakeTypeList< IndexList< 0, 1, 5, 4,   8, 13, 16, 12,   25 > ,
     *                 IndexList< 1, 2, 6, 5,   9, 14, 17, 13,   24 > ,
     *                 IndexList< 2, 3, 7, 6,  10, 15, 18, 14,   26 > ,
     *                 IndexList< 0, 4, 7, 3,  12, 19, 15, 11,   23 > ,
     *                 IndexList< 0, 3, 2, 1,  11, 10,  9,  8,   21 > ,
     *                 IndexList< 4, 5, 6, 7,  16, 17, 18, 19,   22 > >::type
     *   HexahedronFaceNodeMap ;
     *
     *
     *   After refinement (new nodes = *):
     *
     *
     *          7         18         6
     *           o---------*--------o
     *          /|        /|       /|
     *         / |       / |      / |
     *     19 /  |    22/  |     /  |
     *       *---------*--------*17 |
     *      /| 15*----/|---*---/|---*14
     *     / |  /|   / |  /|26/ |  /|                                 |   (PARENT) Linear 8-Node Hexahedron
     *  4 /  | / |16/  | / | /  | / |                                 |            3D Element Edge Node to mid-edge quadratic node map
     *   o---------*--20----o5  |/  |                                 |
     *   | 23*---|-|---*- 10|---*24 |                                 |         { 8, 9, 10, 11,
     *   |  /|  3o-|--/|---*|--/|---o 2                               |           16, 17, 18, 19,
     *   | / |  /  | / |  / | / |  /                                  |           12, 13, 14, 15 }
     *   |/  | / 25|/  | /  |/  | /                                   |
     * 12*---------*--------*13 |/                                    |            Face to mid-face quadratic node map
     *   | 11*-----|---*----|---* 9                                   |         { 25, 24, 26, 23, 21, 22 }
     *   |  /      |  /21   |  /                                      |           0,   1,  2,  3,  4,  5
     *   | /       | /      | /
     *   |/        |/       |/
     *   o---------*--------o
     *  0          8         1
     *
     *
     *   CHILD Linear 8-Node Hexahedron 3D Element Node Maps:
     * |
     * |  static const UInt child_0[] = {  0, 8, 21, 11, 12, 25, 20, 23 };
     * |  static const UInt child_1[] = {  8, 1, 9, 21, 25, 13, 24, 20 };
     * |  static const UInt child_2[] = { 21, 9, 2, 10, 20, 24, 14, 26 };
     * |  static const UInt child_3[] = { 11, 21, 10, 3, 23, 20, 26, 15 };
     * |  static const UInt child_4[] = { 12, 25, 20, 23, 4, 16, 22, 19 };
     * |  static const UInt child_5[] = { 25, 13, 24, 20, 16, 5, 17, 22 };
     * |  static const UInt child_6[] = { 20, 24, 14, 26, 22, 17, 6, 18 };
     * |  static const UInt child_7[] = { 23, 20, 26, 15, 19, 22, 18, 7 };
     * |
     *
     */

    /*  New ref topo info Hex8
     *  ----------------------
     *
     *  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
     */
#if 0
    template<>
    RefTopoX
    RefinementTopologyExtra< shards::Hexahedron<8> >::refinement_topology =
      {
        {0,   0,         0,      0,       1,             {0.0, 0.0, 0.0} },
        {1,   0,         1,      0,       1,             {1.0, 0.0, 0.0} },
        {2,   0,         2,      0,       1,             {1.0, 1.0, 0.0} },
        {3,   0,         3,      0,       1,             {0.0, 1.0, 0.0} },
        {4,   0,         4,      0,       1,             {0.0, 0.0, 1.0} },
        {5,   0,         5,      0,       1,             {1.0, 0.0, 1.0} },
        {6,   0,         6,      0,       1,             {1.0, 1.0, 1.0} },
        {7,   0,         7,      0,       1,             {0.0, 1.0, 1.0} },

        {8,   1,         0,      0,       1,             {0.5, 0.0, 0.0} },
        {9,   1,         1,      0,       1,             {1.0, 0.5, 0.0} },
        {10,  1,         2,      0,       1,             {0.5, 1.0, 0.0} },
        {11,  1,         3,      0,       1,             {0.0, 0.5, 0.0} },

        {12,  1,         8,      0,       1,             {0.0, 0.0, 0.5} },
        {13,  1,         9,      0,       1,             {1.0, 0.0, 0.5} },
        {14,  1,         10,     0,       1,             {1.0, 1.0, 0.5} },
        {15,  1,         11,     0,       1,             {0.0, 1.0, 0.5} },

        {16,  1,         4,      0,       1,             {0.5, 0.0, 1.0} },
        {17,  1,         5,      0,       1,             {1.0, 0.5, 1.0} },
        {18,  1,         6,      0,       1,             {0.5, 1.0, 1.0} },
        {19,  1,         7,      0,       1,             {0.0, 0.5, 1.0} },

        {20,  3,         0,      0,       1,             {0.5, 0.5, 0.5} },

        {21,  2,         4,      0,       1,             {0.5, 0.5, 0.0} },
        {22,  2,         5,      0,       1,             {0.5, 0.5, 1.0} },
        {23,  2,         3,      0,       1,             {0.0, 0.5, 0.5} },
        {24,  2,         1,      0,       1,             {1.0, 0.5, 0.5} },
        {25,  2,         0,      0,       1,             {0.5, 0.0, 0.5} },
        {26,  2,         2,      0,       1,             {0.5, 1.0, 0.5} }

      };
#endif

#if 0
    template<>
    RefTopoX
    RefinementTopologyExtra< shards::Hexahedron<27> >::refinement_topology =
      {
        {0,   0,         0,      0,       1,             {0.0, 0.0, 0.0} }
      };
#endif

    /*--------------------------------------------------------------------------*/
    /**
     *                                   PARENT Quadratic 20-Node Hexahedron Nodes
     *          7         18         6   (SPACE_DIM = 3!)
     *           o--------o---------o
     *          /|                 /|
     *         / |                / |
     *        /  |               /  |
     *     19o   |            17o   |
     *      /  15o             /    o14
     *     /     |            /     |
     *  4 /      | 16        /      |
     *   o---------o--------o 5     |
     *   |       |       10 |       |
     *   |     3 o-------o--|-------o 2
     *   |      /           |      /
     *   |     /            |     /
     * 12o    /             o13  /
     *   |   o11            |   o9
     *   |  /               |  /
     *   | /                | /
     *   |/                 |/
     *   o---------o--------o
     *  0          8         1
     *
     *   PARENT Quadratic 20-Node Hexahedron 3D Element Edge Node Map:
     *
     * |
     * | static const UInt edge_0[]  = { 0, 1,  8 };
     * | static const UInt edge_1[]  = { 1, 2,  9 };
     * | static const UInt edge_2[]  = { 2, 3, 10 };
     * | static const UInt edge_3[]  = { 3, 0, 11 };
     * | static const UInt edge_4[]  = { 4, 5, 16 };
     * | static const UInt edge_5[]  = { 5, 6, 17 };
     * | static const UInt edge_6[]  = { 6, 7, 18 };
     * | static const UInt edge_7[]  = { 7, 4, 19 };
     * | static const UInt edge_8[]  = { 0, 4, 12 };
     * | static const UInt edge_9[]  = { 1, 5, 13 };
     * | static const UInt edge_10[] = { 2, 6, 14 };
     * | static const UInt edge_11[] = { 3, 7, 15 };
     * |
     *
     *   CHILD Quadratic 20-Node Hexahedron 3D Element Node Maps:
     * |
     * |  // Child node tables use Two groups of nodes:
     * |  // a) vertices
     * |  // b) outer edge mid-points
     * |
     * |  static const UInt child_0[] = { 0, 8, 21, 11, 12, 25, 20, 23,
     * |                                  27, 60, 67, 34, 35, 59, 79, 74, 51, 75, 77, 58 };
     * |
     * |  static const UInt child_1[] = { 8, 1, 9, 21, 25, 13, 24, 20,
     * |                                  28, 29, 68, 60, 59, 36, 69, 79, 52, 53, 78, 75 };
     * |
     * |  static const UInt child_2[] = { 21, 9, 2, 10, 20, 24, 14, 26,
     * |                                  68, 30, 31, 61, 79, 69, 37, 62, 78, 54, 55, 76 };
     * |
     * |  static const UInt child_3[] = { 11, 21, 10, 3, 23, 20, 26, 15,
     * |                                  67, 61, 32, 33, 74, 79, 62, 38, 77, 76, 56, 57 };
     * |
     * |  static const UInt child_4[] = { 12, 25, 20, 23, 4, 16, 22, 19,
     * |                                  51, 75, 77, 58, 39, 66, 80, 73, 43, 65, 72, 50 };
     * |
     * |  static const UInt child_5[] = { 25, 13, 24, 20, 16, 5, 17, 22,
     * |                                  52, 53, 78, 75, 66, 40, 70, 80, 44, 45, 71, 65 };
     * |
     * |  static const UInt child_6[] = { 20, 24, 14, 26, 22, 17, 6, 18,
     * |                                  78, 54, 55, 76, 80, 70, 41, 63, 71, 46, 47, 64 };
     * |
     * |  static const UInt child_7[] = { 23, 20, 26, 15, 19, 22, 18, 7,
     * |                                  77, 76, 56, 57, 73, 80, 63, 42, 72, 64, 48, 49 };
     * |
     *
     */

    /*--------------------------------------------------------------------------*/
    /**
     *                                   PARENT Quadratic 27-Node Hexahedron Nodes
     *          7         18         6   (SPACE_DIM = 3!)
     *           o--------o---------o
     *          /|                 /|
     *         / |                / |
     *        /  |               /  |
     *     19o   |            17o   |
     *      /  15o             /    o14
     *     /     |            /     |
     *  4 /      | 16        /      |
     *   o---------o--------o 5     |
     *   |       |       10 |       |
     *   |     3 o-------o--|-------o 2
     *   |      /           |      /
     *   |     /            |     /
     * 12o    /             o13  /
     *   |   o11            |   o9
     *   |  /               |  /
     *   | /                | /
     *   |/                 |/
     *   o---------o--------o
     *  0          8         1
     *
     *
     *           x--------x---------x
     *          /|                 /|
     *         / |                / |
     *        /  |   22          /  |
     *       x   |    o         x   |
     *      /    x       o26   /    x     (Node #20 is at centroid of element)
     *     /     |            /     |
     *    /      |           /      |     "2D surface" containing nodes 0, 8, 1, 13, 5, 16, 4, 12 has
     *   x---------x--------x       |      node 25 at center....
     *   | 23o   |          |   o24 |
     *   |       x-------x--|-------x
     *   |      /           |      /
     *   |     /  25        |     /
     *   x    /    o        x    /
     *   |   x        o21   |   x
     *   |  /               |  /
     *   | /                | /
     *   |/                 |/
     *   x---------x--------x
     *
     *
     *   PARENT Quadratic 27-Node Hexahedron 3D Element Edge Node Map:
     *
     * |
     * | static const UInt edge_0[]  = { 0, 1,  8 };
     * | static const UInt edge_1[]  = { 1, 2,  9 };
     * | static const UInt edge_2[]  = { 2, 3, 10 };
     * | static const UInt edge_3[]  = { 3, 0, 11 };
     * | static const UInt edge_4[]  = { 4, 5, 16 };
     * | static const UInt edge_5[]  = { 5, 6, 17 };
     * | static const UInt edge_6[]  = { 6, 7, 18 };
     * | static const UInt edge_7[]  = { 7, 4, 19 };
     * | static const UInt edge_8[]  = { 0, 4, 12 };
     * | static const UInt edge_9[]  = { 1, 5, 13 };
     * | static const UInt edge_10[] = { 2, 6, 14 };
     * | static const UInt edge_11[] = { 3, 7, 15 };
     * |
     *
     *   Refined 27-Node Hexahedron Edge node tables:
     * |
     * | static const UInt edge_0[] = { 0, 1,  8, 27, 28 };
     * | static const UInt edge_1[] = { 1, 2,  9, 29, 30 };
     * | static const UInt edge_2[] = { 2, 3, 10, 31, 32 };
     * | static const UInt edge_3[] = { 3, 0, 11, 33, 34 };
     * | static const UInt edge_4[] = { 4, 5, 16, 43, 44 };
     * | static const UInt edge_5[] = { 5, 6, 17, 45, 46 };
     * | static const UInt edge_6[] = { 6, 7, 18, 47, 48 };
     * | static const UInt edge_7[] = { 7, 4, 19, 49, 50 };
     * | static const UInt edge_8[] = { 0, 4, 12, 35, 39 };
     * | static const UInt edge_9[] = { 1, 5, 13, 36, 40 };
     * | static const UInt edge_10[] =  { 2, 6, 14, 37, 41 };
     * | static const UInt edge_11[] =  { 3, 7, 15, 38, 42 };
     * |
     *
     *   CHILD 27-Node Hexahedron 3D Element Node Maps:
     * |
     * | // Child node tables use Four groups of nodes:
     * | // a) vertices
     * | // b) outer edge mid-points
     * | // c) centroid
     * | // d) mid-face points
     * |
     * |  static const UInt child_0[] = { 0, 8, 21, 11, 12, 25, 20, 23,
     * |                                  27, 60, 67, 34, 35, 59, 79, 74, 51, 75, 77, 58,
     * |                                  81,   89, 117, 97, 113, 105, 121  };
     * |
     * |  static const UInt child_1[] = { 8, 1, 9, 21, 25, 13, 24, 20,
     * |                                  28, 29, 68, 60, 59, 36, 69, 79, 52, 53, 78, 75,
     * |                                  82,   92, 118, 113, 101, 106, 122 };
     * |
     * |  static const UInt child_2[] = { 21, 9, 2, 10, 20, 24, 14, 26,
     * |                                  68, 30, 31, 61, 79, 69, 37, 62, 78, 54, 55, 76,
     * |                                  83,   91, 119, 114, 102, 122, 109 };
     * |
     * |  static const UInt child_3[] = { 11, 21, 10, 3, 23, 20, 26, 15,
     * |                                  67, 61, 32, 33, 74, 79, 62, 38, 77, 76, 56, 57,
     * |                                  84,   90, 120, 100, 114, 121, 110 };
     * |
     * |  static const UInt child_4[] = { 12, 25, 20, 23, 4, 16, 22, 19,
     * |                                  51, 75, 77, 58, 39, 66, 80, 73, 43, 65, 72, 50,
     * |                                  85,   117, 93, 98, 116, 108, 124  };
     * |
     * |  static const UInt child_5[] = { 25, 13, 24, 20, 16, 5, 17, 22,
     * |                                  52, 53, 78, 75, 66, 40, 70, 80, 44, 45, 71, 65,
     * |                                  86,   118, 94, 116, 104, 107, 123 };
     * |
     * |  static const UInt child_6[] = { 20, 24, 14, 26, 22, 17, 6, 18,
     * |                                  78, 54, 55, 76, 80, 70, 41, 63, 71, 46, 47, 64,
     * |                                  87,   119, 95, 115, 103, 123, 112 };
     * |
     * |  static const UInt child_7[] = { 23, 20, 26, 15, 19, 22, 18, 7,
     * |                                  77, 76, 56, 57, 73, 80, 63, 42, 72, 64, 48, 49,
     * |                                  88,   120, 96, 99, 115, 124, 111  };
     * |
     * |
     * |  Refined Hexagonal Element "Exterior" Faces
     * |  (Local Face node numbering for 'Hierarchical/Consistent' Hex objects)
     * |
     * |   3    14    6   13     2
     * |    o----*----o----*----o
     * |    |         |         |
     * |    |   24    |    23   |
     * |  15*    *    *19  *    *12
     * |    |         |         |
     * |    |        8|    18   |
     * |  7 o----*----o----*----o 5
     * |    |   20    |         |
     * |    |         |         |
     * |  16*    *  17*    *    *11
     * |    |   21    |   22    |
     * |    |         |         |
     * |    o----*----o----*----o
     * |   0     9    4   10     1
     * |
     * |
     * | Hexagonal object face topology child-nodes:
     * |
     * |  // Face node tables use Six groups of nodes:
     * |  // a) vertices                 (Local nodes: 0-1-2-3               )
     * |  // b) edge mid-points          (Local nodes: 4-5-6-7               )
     * |  // c) centroid                 (Local node : 8                     )
     * |  // d) edge quater points       (Local nodes: 9-10-11-12-13-14-15-16)
     * |  // e) interior edge mid-points (Local nodes: 17-18-19-20           )
     * |  // f) mid-quadrant points      (Local nodes: 21-22-23-24           )
     * |
     * |  static const UInt face_0[] = { 0, 1, 5, 4,   8, 13, 16, 12,  25,
     * |                                 27, 28, 36, 40, 44, 43, 39, 35,
     * |                                 59, 52, 66, 51,  105, 106, 107, 108 };
     * |
     * |  static const UInt face_1[] = { 1, 2, 6, 5,   9, 14, 17, 13,  24,
     * |                                 29, 30, 37, 41, 46, 45, 40, 36,
     * |                                 69, 54, 70, 53,  101, 102, 103, 104 };
     * |
     * |  static const UInt face_2[] = { 2, 3, 7, 6,  10, 15, 18, 14,  26,
     * |                                 31, 32, 38, 42, 48, 47, 41, 37,
     * |                                 62, 56, 63, 55,  109, 110, 111, 112 };
     * |
     * |  static const UInt face_3[] = { 0, 4, 7, 3, 12, 19, 15, 11,  23,
     * |                                 35, 39, 50, 49, 42, 38, 33, 34,
     * |                                 58, 73, 57, 74  97, 98, 99, 100   };
     * |
     * |  static const UInt face_4[] = { 0, 3, 2, 1,   11, 10, 9, 8   21,
     * |                                 34, 33, 32, 31, 30, 29, 28, 27,
     * |                                 67, 61, 68, 60,  89, 90, 91, 92   };
     * |
     * |  static const UInt face_5[] = { 4, 5, 6, 7,  16, 17, 18, 19,  22,
     * |                                 43, 44, 45, 46, 47, 48, 49, 50,
     * |                                 65, 71, 64, 72,  93, 94, 95, 96   };
     *
     */

    /*--------------------------------------------------------------------*/
    // Hexahedrons 8, 20, or 27 nodes (3D solid).
    const MeshObjTopology * hex( UInt nnode)
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology hex8(shards::getCellTopologyData<shards::Hexahedron<8> >());
      static MeshObjTopology hex20(shards::getCellTopologyData<shards::Hexahedron<20> >());
      static MeshObjTopology hex27(shards::getCellTopologyData<shards::Hexahedron<27> >());

      static bool first = true ;

      if ( first ) {

        first = false ;

        { // Child topology and nodes

          static const_top_ptr hex8_child[]  =
            { &hex8 , &hex8 , &hex8 , &hex8 , &hex8 , &hex8 , &hex8 , &hex8 };

          static const_top_ptr hex20_child[]  =
            { &hex20, &hex20, &hex20, &hex20, &hex20, &hex20, &hex20, &hex20 };

          static const_top_ptr hex27_child[]  =
            { &hex27, &hex27, &hex27, &hex27, &hex27, &hex27, &hex27, &hex27 };

          // Child node tables use Four groups of nodes:
          // a) vertices
          // b) outer edge mid-points
          // c) centroid
          // d) mid-face points

          static const UInt child_0[] = { 0, 8, 21, 11, 12, 25, 20, 23,
                                          27, 60, 67, 34, 35, 59, 79, 74, 51, 75, 77, 58,
                                          81,   89, 117, 97, 113, 105, 121  , EUA};

          static const UInt child_1[] = { 8, 1, 9, 21, 25, 13, 24, 20,
                                          28, 29, 68, 60, 59, 36, 69, 79, 52, 53, 78, 75,
                                          82,   92, 118, 113, 101, 106, 122 , EUA};

          static const UInt child_2[] = { 21, 9, 2, 10, 20, 24, 14, 26,
                                          68, 30, 31, 61, 79, 69, 37, 62, 78, 54, 55, 76,
                                          83,   91, 119, 114, 102, 122, 109 , EUA};

          static const UInt child_3[] = { 11, 21, 10, 3, 23, 20, 26, 15,
                                          67, 61, 32, 33, 74, 79, 62, 38, 77, 76, 56, 57,
                                          84,   90, 120, 100, 114, 121, 110 , EUA};

          static const UInt child_4[] = { 12, 25, 20, 23, 4, 16, 22, 19,
                                          51, 75, 77, 58, 39, 66, 80, 73, 43, 65, 72, 50,
                                          85,   117, 93, 98, 116, 108, 124  , EUA};

          static const UInt child_5[] = { 25, 13, 24, 20, 16, 5, 17, 22,
                                          52, 53, 78, 75, 66, 40, 70, 80, 44, 45, 71, 65,
                                          86,   118, 94, 116, 104, 107, 123 , EUA};

          static const UInt child_6[] = { 20, 24, 14, 26, 22, 17, 6, 18,
                                          78, 54, 55, 76, 80, 70, 41, 63, 71, 46, 47, 64,
                                          87,   119, 95, 115, 103, 123, 112 , EUA};

          static const UInt child_7[] = { 23, 20, 26, 15, 19, 22, 18, 7,
                                          77, 76, 56, 57, 73, 80, 63, 42, 72, 64, 48, 49,
                                          88,   120, 96, 99, 115, 124, 111  , EUA};

          static const UInt * child_node_table[]  =
            { child_0 , child_1 , child_2 , child_3 ,
              child_4 , child_5 , child_6 , child_7  };

          // for hex20 elements:
          static const UInt hex20_child_0[] = { 0, 8, 21, 11, 12, 25, 20, 23,
                                                27, 60, 67, 34, 35, 59, 79, 74, 51, 75, 77, 58 , EUA };

          static const UInt hex20_child_1[] = { 8, 1, 9, 21, 25, 13, 24, 20,
                                                28, 29, 68, 60, 59, 36, 69, 79, 52, 53, 78, 75 , EUA };

          static const UInt hex20_child_2[] = { 21, 9, 2, 10, 20, 24, 14, 26,
                                                68, 30, 31, 61, 79, 69, 37, 62, 78, 54, 55, 76 , EUA };

          static const UInt hex20_child_3[] = { 11, 21, 10, 3, 23, 20, 26, 15,
                                                67, 61, 32, 33, 74, 79, 62, 38, 77, 76, 56, 57 , EUA };

          static const UInt hex20_child_4[] = { 12, 25, 20, 23, 4, 16, 22, 19,
                                                51, 75, 77, 58, 39, 66, 80, 73, 43, 65, 72, 50 , EUA };

          static const UInt hex20_child_5[] = { 25, 13, 24, 20, 16, 5, 17, 22,
                                                52, 53, 78, 75, 66, 40, 70, 80, 44, 45, 71, 65 , EUA };

          static const UInt hex20_child_6[] = { 20, 24, 14, 26, 22, 17, 6, 18,
                                                78, 54, 55, 76, 80, 70, 41, 63, 71, 46, 47, 64 , EUA };

          static const UInt hex20_child_7[] = { 23, 20, 26, 15, 19, 22, 18, 7,
                                                77, 76, 56, 57, 73, 80, 63, 42, 72, 64, 48, 49 , EUA };


          static const UInt * hex20_child_node_table[]  =
            { hex20_child_0 , hex20_child_1 , hex20_child_2 , hex20_child_3 ,
              hex20_child_4 , hex20_child_5 , hex20_child_6 , hex20_child_7  };



          // Edge topology and node tables including edges' child-nodes

          static const UInt edge_0[] = { 0, 1,  8, 27, 28 , EUA};
          static const UInt edge_1[] = { 1, 2,  9, 29, 30 , EUA};
          static const UInt edge_2[] = { 2, 3, 10, 31, 32 , EUA};
          static const UInt edge_3[] = { 3, 0, 11, 33, 34 , EUA};
          static const UInt edge_4[] = { 4, 5, 16, 43, 44 , EUA};
          static const UInt edge_5[] = { 5, 6, 17, 45, 46 , EUA};
          static const UInt edge_6[] = { 6, 7, 18, 47, 48 , EUA};
          static const UInt edge_7[] = { 7, 4, 19, 49, 50 , EUA};
          static const UInt edge_8[] = { 0, 4, 12, 35, 39 , EUA};
          static const UInt edge_9[] = { 1, 5, 13, 36, 40 , EUA};
          static const UInt edge_10[] =  { 2, 6, 14, 37, 41 , EUA};
          static const UInt edge_11[] =  { 3, 7, 15, 38, 42 , EUA};

          static const UInt * edge_table[]  =
            { edge_0 , edge_1 , edge_2 , edge_3 , edge_4  , edge_5 ,
              edge_6 , edge_7 , edge_8 , edge_9 , edge_10 , edge_11  };

          // Face topology and edge tables including faces' child-nodes

          // Each face will have at 25 nodes to cover all cases of refined
          // Hexahedron faces: Quad_4_3D, Quad_8_3D, or Quad_9_3D.
          // Face node tables use Six groups of nodes:
          // a) vertices                 (Local nodes: 0-1-2-3               )
          // b) edge mid-points          (Local nodes: 4-5-6-7               )
          // c) centroid                 (Local node : 8                     )
          // d) edge quater points       (Local nodes: 9-10-11-12-13-14-15-16)
          // e) interior edge mid-points (Local nodes: 17-18-19-20           )
          // f) mid-quadrant points      (Local nodes: 21-22-23-24           )

          static const UInt face_0[] = { 0, 1, 5, 4,   8, 13, 16, 12,  25,
                                         27, 28, 36, 40, 44, 43, 39, 35,
                                         59, 52, 66, 51,  105, 106, 107, 108 , EUA};

          static const UInt face_1[] = { 1, 2, 6, 5,   9, 14, 17, 13,  24,
                                         29, 30, 37, 41, 46, 45, 40, 36,
                                         69, 54, 70, 53,  101, 102, 103, 104 , EUA};

          static const UInt face_2[] = { 2, 3, 7, 6,  10, 15, 18, 14,  26,
                                         31, 32, 38, 42, 48, 47, 41, 37,
                                         62, 56, 63, 55,  109, 110, 111, 112 , EUA};

          static const UInt face_3[] = { 0, 4, 7, 3, 12, 19, 15, 11,  23,
                                         35, 39, 50, 49, 42, 38, 33, 34,
                                         58, 73, 57, 74, 97, 98, 99, 100   , EUA};

          static const UInt face_4[] = { 0, 3, 2, 1,   11, 10, 9, 8,   21,
                                         34, 33, 32, 31, 30, 29, 28, 27,
                                         67, 61, 68, 60,  89, 90, 91, 92    , EUA};

          static const UInt face_5[] = { 4, 5, 6, 7,  16, 17, 18, 19,  22,
                                         43, 44, 45, 46, 47, 48, 49, 50,
                                         65, 71, 64, 72,   93, 94, 95, 96    , EUA};

          static const UInt * face_table[]  =
            { face_0 , face_1 , face_2 , face_3 , face_4 , face_5 };


          // for hex20 elements
          static const UInt hex20_face_0[] = { 0, 1, 5, 4,   8, 13, 16, 12,  25,
                                         27, 28, 36, 40, 44, 43, 39, 35,
                                         59, 52, 66, 51, EUA};

          static const UInt hex20_face_1[] = { 1, 2, 6, 5,   9, 14, 17, 13,  24,
                                         29, 30, 37, 41, 46, 45, 40, 36,
                                         69, 54, 70, 53,  EUA};

          static const UInt hex20_face_2[] = { 2, 3, 7, 6,  10, 15, 18, 14,  26,
                                         31, 32, 38, 42, 48, 47, 41, 37,
                                         62, 56, 63, 55,  EUA};

          static const UInt hex20_face_3[] = { 0, 4, 7, 3, 12, 19, 15, 11,  23,
                                         35, 39, 50, 49, 42, 38, 33, 34,
                                         58, 73, 57, 74, EUA};

          static const UInt hex20_face_4[] = { 0, 3, 2, 1,   11, 10, 9, 8,   21,
                                         34, 33, 32, 31, 30, 29, 28, 27,
                                         67, 61, 68, 60,  EUA};

          static const UInt hex20_face_5[] = { 4, 5, 6, 7,  16, 17, 18, 19,  22,
                                         43, 44, 45, 46, 47, 48, 49, 50,
                                         65, 71, 64, 72,  EUA};

          static const UInt * hex20_face_table[]  =
            { hex20_face_0 , hex20_face_1 , hex20_face_2 , hex20_face_3 , hex20_face_4 , hex20_face_5 };


          static RefinementTopology hex8_refinement(&hex8, 8, hex8_child, 27, child_node_table, 12, edge_table, 6, face_table, 0, NULL, NULL, true);
          static RefinementTopology hex20_refinement(&hex20, 8, hex20_child, 81, hex20_child_node_table, 12, edge_table, 6, hex20_face_table, 0, NULL, NULL, true);
          static RefinementTopology hex27_refinement(&hex27, 8, hex27_child, 125, child_node_table, 12, edge_table, 6, face_table, 0, NULL, NULL, true);
        }
      }

      MeshObjTopology * top = NULL ;

      switch( nnode ) {
      case  8 : top = & hex8  ; break ;
      case 20 : top = & hex20 ; break ;
      case 27 : top = & hex27 ; break ;
      default :
        //throw RuntimeError() << "Invalid nnode specified" << std::endl ;// << StackTrace;
        throw std::runtime_error( "Invalid nnode specified") ; // << std::endl ;// << StackTrace;
      }

      return top ;
    }


    /*---------------------------------------------------------------------*/
    /**
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
     *             /|\
     *            / | \
     *         7 *  |  * 9
     *          /   |   \
     *         /   6|    \
     *      0 o----*|-----o 2
     *         \    *8   /
     *          \   |   /
     *         4 *  |  * 5
     *            \ | /
     *             \|/
     *              o
     *              1
     *
     *   CHILD 4-Node Tetrahedron 3D Object Node Maps:
     * |
     * | static const UInt child_0[] = { 0, 4, 6, 7 };  // srkenno 091410 fixed (used to be {0, 4, 8, 7} )
     * | static const UInt child_1[] = { 4, 1, 5, 8 };
     * | static const UInt child_2[] = { 6, 5, 2, 9 };
     * | static const UInt child_3[] = { 7, 8, 9, 3 };
     * | static const UInt child_4[] = { 8, 7, 6, 4 };
     * | static const UInt child_5[] = { 6, 9, 8, 5 };
     * | static const UInt child_6[] = { 9, 8, 7, 6 };
     * | static const UInt child_7[] = { 5, 6, 4, 8 };
     * |
     *
     **/

    /*  New ref topo info Tet4
     *  ----------------------
     *
     *  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
     */
#if 0
    template<>
    RefTopoX
    RefinementTopologyExtra< shards::Tetrahedron<4> >::refinement_topology =
      {
        {0,   0,         0,      0,       1,             {0.0, 0.0, 0.0} },
        {1,   0,         1,      0,       1,             {1.0, 0.0, 0.0} },
        {2,   0,         2,      0,       1,             {0.0, 1.0, 0.0} },
        {3,   0,         3,      0,       1,             {0.0, 0.0, 1.0} },
        {4,   1,         0,      0,       1,             {0.5, 0.0, 0.0} },
        {5,   1,         1,      0,       1,             {0.5, 0.5, 0.0} },
        {6,   1,         2,      0,       1,             {0.0, 0.5, 0.0} },
        {7,   1,         3,      0,       1,             {0.0, 0.0, 0.5} },
        {8,   1,         4,      0,       1,             {0.5, 0.0, 0.5} },
        {9,   1,         5,      0,       1,             {0.0, 0.5, 0.5} }
      };

#endif

    /*--------------------------------------------------------------------------*/
    /**
     *                         PARENT 10-Node Tetrahedron Object Nodes
     *              3
     *              o
     *             /|\
     *            / | \
     *         7 o  |  o 9    (PARENT) 10-Node Tetrahedron Object
     *          /   |   \              Edge Node Map:
     *         /   6|    \
     *      0 o----o|-----o 2        { {0, 1, 4}, {1, 2, 5}, {2, 0, 6},
     *         \    o8   /             {0, 3, 7}, {1, 3, 8}, {2, 3, 9} };
     *          \   |   /
     *         4 o  |  o 5
     *            \ | /
     *             \|/
     *              o
     *              1
     *
     *  After refinement (new nodes = *):
     *
     *                3
     *                o
     *               /|\
     *              * | *
     *             /  |  \
     *          7 o   |   o 9
     *           /    *    \
     *          *     |     *
     *         /     6|      \
     *      0 o---*--o|--*----o 2
     *         \      o8     /
     *          *     |     *
     *           \    |    /
     *          4 o   |   o 5
     *             \  *  /
     *              * | *
     *               \|/
     *                o
     *                1
     *
     * |  // Child edge node tables
     * |
     * |  static const UInt   edge_0[]  = { 0, 1,  4,  10, 11 };
     * |  static const UInt   edge_1[]  = { 1, 2,  5,  12, 13 };
     * |  static const UInt   edge_2[]  = { 2, 0,  6,  14, 15 };
     * |  static const UInt   edge_3[]  = { 0, 3,  7,  16, 19 };
     * |  static const UInt   edge_4[]  = { 1, 3,  8,  17, 20 };
     * |  static const UInt   edge_5[]  = { 2, 3,  9,  18, 21 };
     * |
     * |  // Child face node (cfn) tables:
     * |  // Local Face (LF) LF0 uses [0:5], LF1 uses [6:11], LF2 uses [12:17], LF3 uses [18:23]
     * |
     * |  static const UInt cfn_0[] = {0, 4, 7, 10, 27, 16, 4, 6, 7, 31, 25, 27, 0, 7, 6, 16, 25, 15, 0, 6, 4, 15, 31, 10 };
     * |  static const UInt cfn_1[] = {4, 1, 8, 11, 17, 34, 1, 5, 8, 12, 26, 17, 4, 8, 5, 34, 26, 20, 4, 5, 1, 30, 12, 11 };
     * |  static const UInt cfn_2[] = {6, 5, 9, 24, 28, 33, 5, 2, 9, 13, 18, 28, 6, 9, 2, 33, 18, 14, 6, 2, 5, 14, 13, 24 };
     * |  static const UInt cfn_3[] = {7, 8, 3, 23, 20, 19, 8, 9, 3, 32, 21, 20, 7, 3, 9, 19, 21, 29, 7, 9, 8, 29, 32, 33 };
     * |  static const UInt cfn_4[] = {8, 7, 4, 23, 27, 34, 7, 6, 4, 25, 31, 27, 8, 4, 6, 34, 31, 22, 8, 6, 7, 22, 25, 23 };
     * |  static const UInt cfn_5[] = {6, 9, 5, 33, 28, 24, 9, 8, 5, 32, 26, 28, 6, 5, 8, 24, 26, 22, 6, 8, 9, 22, 32, 33 };
     * |  static const UInt cfn_6[] = {9, 8, 6, 32, 22, 33, 8, 7, 6, 23, 25, 22, 9, 6, 7, 33, 25, 29, 9, 7, 8, 29, 23, 32 };
     * |  static const UInt cfn_7[] = {5, 6, 8, 24, 22, 26, 6, 4, 8, 31, 34, 22, 5, 8, 4, 26, 34, 30, 5, 4, 6, 30, 31, 24 };
     * |
     *
     *
     *
     *      Face #0:               Face #1:               Face #2:
     *
     *           3                      3                      3
     *           o                      o                      o
     *          / \                    / \                    / \
     *       19*   *20              20*   *21              21*   *19
     *        / 23  \                / 32  \                / 29  \    <
     *     7 o---*---o 8          8 o---*---o 9          9 o---*---o 7  \
     *      / \     / \            / \     / \            / \     / \    \
     *   16* 27*   *34 *17      17* 26*   *28 *18      18* 33*   *25 *16  |
     *    /     \ /     \        /     \ /     \        /     \ /     \   |
     *   o---*---o---*---o      o---*---o---*---o      o---*---o---*---o  #
     *  0   10   4   11   1    1   12   5   13   2    2   14   6   15   0
     *   #                      #
     *    \                      \
     *     \                      \
     *      -->                    -->
     *
     *      Face #3:
     *
     *           1
     *           o
     *          / \
     *       11*   *12
     *        / 30  \
     *     4 o---*---o 5
     *      / \     / \
     *   10* 31*   *24 *13
     *    /     \ /     \
     *   o---*---o---*---o
     *  0   15   6   14   2
     *   #
     *    \
     *     \
     *      -->
     *
     *   Various Interior Faces of Children:
     *
     *        7              9        7   25    6     8   26    5
     *        o              o         o---*---o       o---*---o
     *       / \            / \         \     /         \     /
     *    23*   *25      32*   *33     27*   *31       34*   *30
     *     /  22 \        / 22  \         \ /             \ /
     *   8o---*---o10   8o---*---o6        o               o
     *     \     /        \     /          4               4
     *    26*   *24      34*   *31
     *       \ /            \ /
     *        o              o
     *        5              4
     *
     *        9              9
     *        o              o
     *       / \            / \
     *    33*   *28      29*   *32
     *     /     \        /     \
     *    o---*---o      o---*---o
     *   6   24    5    7    23   8
     *
     *
     *   Edge node tables for refined 10-Node Tetrahedrons:
     * |
     * | static const UInt   edge_0[]  = { 0, 1,  4,  10, 11 };
     * | static const UInt   edge_1[]  = { 1, 2,  5,  12, 13 };
     * | static const UInt   edge_2[]  = { 2, 0,  6,  14, 15 };
     * | static const UInt   edge_3[]  = { 0, 3,  7,  16, 19 };
     * | static const UInt   edge_4[]  = { 1, 3,  8,  17, 20 };
     * | static const UInt   edge_5[]  = { 2, 3,  9,  18, 21 };
     * |
     *
     *   Face node tables for refined 10-Node Tetrahedrons:
     * |
     * | static const UInt face_0[] = { 0, 1, 3, 4, 8, 7, 10, 11, 17, 20, 19, 16, 27, 34, 23 };
     * | static const UInt face_1[] = { 1, 2, 3, 5, 9, 8, 12, 13, 18, 21, 20, 17, 26, 28, 32 };
     * | static const UInt face_2[] = { 0, 3, 2, 7, 9, 6, 16, 19, 21, 18, 14, 15, 25, 29, 33 };
     * | static const UInt face_3[] = { 0, 2, 1, 6, 5, 4, 15, 14, 13, 12, 11, 10, 31, 24, 30 };
     * |
     *
     *   CHILD 10-Node Tetrahedron Object Node Maps:
     *
     * |
     * | static const UInt child_0[] = { 0, 4, 8, 7, 10, 31, 15, 16, 27, 25 };
     * | static const UInt child_1[] = { 4, 1, 5, 8, 11, 12, 30, 34, 17, 26 };
     * | static const UInt child_2[] = { 6, 5, 2, 9, 24, 13, 14, 33, 28, 18 };
     * | static const UInt child_3[] = { 7, 8, 9, 3, 23, 32, 29, 19, 20, 21 };
     * | static const UInt child_4[] = { 8, 7, 6, 4, 23, 25, 22, 34, 27, 31 };
     * | static const UInt child_5[] = { 6, 9, 8, 5, 33, 32, 22, 24, 28, 26 };
     * | static const UInt child_6[] = { 9, 8, 7, 6, 32, 23, 29, 33, 22, 25 };
     * | static const UInt child_7[] = { 5, 6, 4, 8, 24, 31, 30, 26, 22, 34 };
     * |
     *
     *
     **/

    /*--------------------------------------------------------------------------*/
    /**
     *                          PARENT Semi-Linear 8-Node Tetrahedron Nodes
     *              3           (SPACE_DIM = 3!)
     *              o
     *             /|\
     *            / | \
     *           /  |  \        "Front faces"
     *          /   |   \        Node 4 on 2D surface containing nodes 0, 1, 3
     *         /    |    \       Node 5 on 2D surface containing nodes 1, 2, 3
     *      0 o-----|-----o 2
     *         \ 4  |  + /      (PARENT) Semi-Linear 8-Node Tetrahedron
     *          \ + | 5 /                3D Element Edge Node Map:
     *           \  |  /
     *            \ | /        { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} }
     *             \|/
     *              o
     *              1
     *
     *
     *              3
     *              o
     *             /|\         "Back faces (note mid-face-node does not follow face ordering!)"
     *            / | \         Node 7 on 2D surface containing nodes 0, 3, 2
     *           /  |  \        Node 6 on 2D surface containing nodes 0, 2, 1
     *          /   |+  \
     *         /    |7   \
     *      0 o-----|-----o 2
     *         \   6|    /
     *          \  +|   /
     *           \  |  /
     *            \ | /
     *             \|/
     *              o
     *              1
     *
     *  After refinement (new nodes = *):
     *
     *              3
     *              o
     *             /|\
     *            / | \
     *         11*  |  *13
     *          /   |   \
     *         /  10|    \
     *      0 o----*|-----o 2
     *         \    *12  /
     *          \   |   /
     *         8 *  |  * 9
     *            \ | /
     *             \|/
     *              o
     *              1
     *
     * |  // Child edge node tables
     * |
     * |  static const UInt   edge_0[]  = { 0, 1,   8 };
     * |  static const UInt   edge_1[]  = { 1, 2,   9 };
     * |  static const UInt   edge_2[]  = { 2, 0,  10 };
     * |  static const UInt   edge_3[]  = { 0, 3,  11 };
     * |  static const UInt   edge_4[]  = { 1, 3,  12 };
     * |  static const UInt   edge_5[]  = { 2, 3,  13 };
     * |
     *
     * |
     * |  // Child face node (cfn) tables:
     * |  // Local Face (LF)
     * |    LF0 uses [0:3], LF1 uses [4:7], LF2 uses [8:11], LF3 uses [12:15]
     * |
     * |  static const UInt cfn_0[]  =
     * |    { 0, 8, 11, 14,   8, 10, 11, 30,  0, 11, 10, 17,  0, 10, 8, 16   };
     * |  static const UInt cfn_1[]  =
     * |    { 8, 1, 12, 18,   1, 9, 12, 15,   8, 12, 9, 31,   8, 9, 1, 24    };
     * |  static const UInt cfn_2[]  =
     * |    { 10, 9, 13, 32,  9, 2, 13, 19,   10, 13, 2, 25,  10, 2, 9, 20   };
     * |  static const UInt cfn_3[]  =
     * |    { 11, 12, 3, 22,  12, 13, 3, 23,  11, 3, 13, 21,  11, 13, 12, 33 };
     * |  static const UInt cfn_4[]  =
     * |    { 12, 11, 8, 4,   11, 10, 8, 30,  12, 8, 10, 29,  12, 10, 11, 26 };
     * |  static const UInt cfn_5[]  =
     * |    { 10, 13, 9, 32,  13, 12, 9, 5,   10, 9, 12, 27,  10, 12, 13, 28 };
     * |  static const UInt cfn_6[]  =
     * |    { 13, 12, 10, 28, 12, 11, 10, 26, 13, 10, 11, 7,  13, 11, 12, 33 };
     * |  static const UInt cfn_7[]  =
     * |    { 9, 10, 12, 27,  10, 8, 12, 29,  9, 12, 8, 31,   9, 8, 10, 6    };
     * |
     *
     *
     *      Face #0:               Face #1:               Face #2:
     *
     *           3                      3                      3
     *           o                      o                      o
     *          / \                    / \                    / \
     *         / 16\                  / 19\                  / 24\
     *        /  *  \                /  *  \                /  *  \    <
     *     11*-------*12          12*-------*13          13*-------*11  \
     *      / \  o  / \            / \  o  / \            / \  o  / \    \
     *     /14 \ 4 /15 \          /17 \ 5 /18 \          /25 \ 7 / 23\    |
     *    /  *  \ /  *  \        /  *  \ /  *  \        /  *  \ /  *  \   |
     *   o-------*-------o      o-------*-------o      o-------*-------o  #
     *  0        8        1    1        9        2    2        10       0
     *   #                      #
     *    \                      \
     *     \                      \
     *      -->                    -->
     *
     *      Face #3:
     *
     *           1
     *           o
     *          / \
     *         /22 \
     *        /  *  \
     *      8*-------*9
     *      / \  o  / \
     *     /20 \ 6 / 21\
     *    /  *  \ /  *  \
     *   o-------*-------o
     *  0       10        2
     *   #
     *    \
     *     \
     *      -->
     *
     *   Various Interior Faces of Children:
     *
     *       11              13       11       10     12        9
     *        *              *         *-------*       *-------*
     *       / \            / \         \  *  /         \  *  /
     *      /26 \          /28 \         \30 /           \31 /
     *     /  *  \        /  *  \         \ /             \ /
     *  12*-------*10  10*-------*12       *               *
     *     \  *  /        \  *  /          8               8
     *      \27 /          \29 /
     *       \ /            \ /
     *        *              *
     *        9              8
     *
     *       13              13
     *        *              *
     *       / \            / \
     *      /32 \          /33 \
     *     /  *  \        /  *  \
     *    *-------*      *-------*
     *  10         9    11       12
     *
     * | // tet8 face node tables
     * |
     * | static const UInt    t8_face_0[] = { 0, 1, 3,  4,   8, 12, 11,  14, 15, 16 };
     * | static const UInt    t8_face_1[] = { 1, 2, 3,  5,   9, 13, 12,  17, 18, 19 };
     * | static const UInt    t8_face_2[] = { 0, 3, 2,  7,  11, 13, 10,  23, 24, 25 };
     * | static const UInt    t8_face_3[] = { 0, 2, 1,  6,  10, 9, 8,  20, 21, 22 };
     * |
     *
     *   CHILD 8-Node Tetrahedron Object Node Maps:
     *
     * |
     * | static const UInt child_0[] = {  0, 8, 10, 11,  14, 30, 23, 20 };
     * | static const UInt child_1[] = {  8, 1, 9, 12,  15, 17, 31, 22 };
     * | static const UInt child_2[] = { 10, 9, 2, 13,  32, 18, 25, 21 };
     * | static const UInt child_3[] = { 11, 12, 13, 3,  16, 19, 24, 33 };
     * | static const UInt child_4[] = { 12, 11, 10, 8,   4, 30, 29, 26 };
     * | static const UInt child_5[] = { 10, 13, 12, 9,  32, 5, 27, 28 };
     * | static const UInt child_6[] = { 13, 12, 11, 10,  28, 26, 7, 33 };
     * | static const UInt child_7[] = {  9, 10, 8, 12,  27, 29, 31, 6 };
     * |
     *
     **/

    /*--------------------------------------------------------------------*/
    // Tetrahedron with 4, 8, or 10 Nodes. Object is a 3D solid.
    const MeshObjTopology * tet( UInt nnode )
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology tet4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
      static MeshObjTopology tet8(shards::getCellTopologyData<shards::Tetrahedron<8> >());
      static MeshObjTopology tet10(shards::getCellTopologyData<shards::Tetrahedron<10> >());

      static bool first = true ;

      if ( first ) {

        first = false ;

        { // Tet4 and Tet10 Child and child node specifications

          static const_top_ptr tet4_child[]  =
            { &tet4 , &tet4 , &tet4 , &tet4 , &tet4 , &tet4 , &tet4 , &tet4 };

          static const_top_ptr tet10_child[]  =
            { &tet10, &tet10, &tet10, &tet10, &tet10, &tet10, &tet10, &tet10 };

          static const UInt   tet4_child_0[] = { 0, 4, 6, 7,   10, 31, 15, 16, 27, 25 , EUA};
          static const UInt   tet4_child_1[] = { 4, 1, 5, 8,   11, 12, 30, 34, 17, 26 , EUA};
          static const UInt   tet4_child_2[] = { 6, 5, 2, 9,   24, 13, 14, 33, 28, 18 , EUA};
          static const UInt   tet4_child_3[] = { 7, 8, 9, 3,   23, 32, 29, 19, 20, 21 , EUA};
          static const UInt   tet4_child_4[] = { 8, 7, 6, 4,   23, 25, 22, 34, 27, 31 , EUA};
          static const UInt   tet4_child_5[] = { 6, 9, 8, 5,   33, 32, 22, 24, 28, 26 , EUA};
          static const UInt   tet4_child_6[] = { 9, 8, 7, 6,   32, 23, 29, 33, 22, 25 , EUA};
          static const UInt   tet4_child_7[] = { 5, 6, 4, 8,   24, 31, 30, 26, 22, 34 , EUA};
          static const UInt * tet4_child_node_table[]  =
            { tet4_child_0 , tet4_child_1 , tet4_child_2 , tet4_child_3 ,
              tet4_child_4 , tet4_child_5 , tet4_child_6 , tet4_child_7 };



          // Tet8 Child and child node specifications

          static const_top_ptr tet8_child[]  =
            { &tet8 , &tet8 , &tet8 , &tet8 , &tet8 , &tet8 , &tet8 , &tet8 };

          static const UInt   tet8_child_0[] = {  0, 8, 10, 11,  14, 30, 20, 23 , EUA};
          static const UInt   tet8_child_1[] = {  8, 1, 9, 12,  15, 17, 22, 31 , EUA};
          static const UInt   tet8_child_2[] = { 10, 9, 2, 13,  32, 18, 21, 25 , EUA};
          static const UInt   tet8_child_3[] = { 11, 12, 13, 3,  16, 19, 33, 24 , EUA};
          static const UInt   tet8_child_4[] = { 12, 11, 10, 8,   4, 30, 26, 29 , EUA};
          static const UInt   tet8_child_5[] = { 10, 13, 12, 9,  32, 5, 28, 27 , EUA};
          static const UInt   tet8_child_6[] = { 13, 12, 11, 10,  28, 26, 33, 7 , EUA};
          static const UInt   tet8_child_7[] = {  9, 10, 8, 12,  27, 29, 6, 31 , EUA};
          static const UInt * tet8_child_node_table[]  =
            { tet8_child_0 , tet8_child_1 , tet8_child_2 , tet8_child_3 ,
              tet8_child_4 , tet8_child_5 , tet8_child_6 , tet8_child_7 };

          // Edge topology and node tables including edges' child-nodes

          // Edge child nodes for both Tet4 and Tet10 edges

          static const UInt   edge_0[]  = { 0, 1,  4,  10, 11 , EUA};
          static const UInt   edge_1[]  = { 1, 2,  5,  12, 13 , EUA};
          static const UInt   edge_2[]  = { 2, 0,  6,  14, 15 , EUA};
          static const UInt   edge_3[]  = { 0, 3,  7,  16, 19 , EUA};
          static const UInt   edge_4[]  = { 1, 3,  8,  17, 20 , EUA};
          static const UInt   edge_5[]  = { 2, 3,  9,  18, 21 , EUA};
          static const UInt * tet4_edge_table[]  =
            { edge_0 , edge_1 , edge_2 , edge_3 , edge_4 , edge_5 };

          // Edge child nodes for Tet8 edges

          static const UInt   tet8_edge_0[]  = { 0, 1,   8 , EUA};
          static const UInt   tet8_edge_1[]  = { 1, 2,   9 , EUA};
          static const UInt   tet8_edge_2[]  = { 2, 0,  10 , EUA};
          static const UInt   tet8_edge_3[]  = { 0, 3,  11 , EUA};
          static const UInt   tet8_edge_4[]  = { 1, 3,  12 , EUA};
          static const UInt   tet8_edge_5[]  = { 2, 3,  13 , EUA};

          static const UInt * tet8_edge_table[]  =
            { tet8_edge_0, tet8_edge_1, tet8_edge_2, tet8_edge_3, tet8_edge_4, tet8_edge_5 };

          // Face topology, nodes, and edge tables

          // Face node tables use Four groups of nodes:
          // a) vertices
          // b) outer edge mid-points
          // c) outer edge quarter points
          // d) inner edge mid-points

          static const UInt tet4_face_0[] = {  0,  1,  3,
                                               4,  8,  7,
                                               10, 11, 17, 20, 19, 16,
                                               27, 34, 23 , EUA};

          static const UInt tet4_face_1[] = {  1,  2,  3,
                                               5,  9,  8,
                                               12, 13, 18, 21, 20, 17,
                                               26, 28, 32 , EUA};

          static const UInt tet4_face_2[] = {  0,  3,  2,
                                               7,  9,  6,
                                               16, 19, 21, 18, 14, 15,
                                               25, 29, 33 , EUA};

          static const UInt tet4_face_3[] = {  0,  2,  1,
                                               6,  5,  4,
                                               15, 14, 13, 12, 11, 10,
                                               31, 24, 30 , EUA};

          static const UInt * tet4_face_table[4] = { tet4_face_0, tet4_face_1, tet4_face_2, tet4_face_3 };

          // tet8 face node tables

          static const UInt   tet8_face_0[] = { 0, 1, 3,  4,   8, 12, 11,  14, 15, 16 , EUA};
          static const UInt   tet8_face_1[] = { 1, 2, 3,  5,   9, 13, 12,  17, 18, 19 , EUA};
          static const UInt   tet8_face_2[] = { 0, 3, 2,  7,  11, 13, 10,  23, 24, 25 , EUA};
          static const UInt   tet8_face_3[] = { 0, 2, 1,  6,  10, 9, 8,  20, 21, 22 , EUA};
          static const UInt * tet8_face_table[4]  =
            { tet8_face_0, tet8_face_1, tet8_face_2, tet8_face_3 };

          static RefinementTopology tet4_refinement(&tet4, 8, tet4_child, 10, tet4_child_node_table, 6, tet4_edge_table, 4, tet4_face_table, 0, NULL, NULL, true);
          static RefinementTopology tet10_refinement(&tet10, 8, tet10_child, 35, tet4_child_node_table, 6, tet4_edge_table, 4, tet4_face_table, 0, NULL, NULL, true);
          static RefinementTopology tet8_refinement(&tet8, 8, tet8_child, 34, tet8_child_node_table, 6, tet8_edge_table, 4, tet8_face_table, 0, NULL, NULL, true);
        }
      }

      MeshObjTopology * top = NULL ;

      switch( nnode ) {
      case  4 : top = &tet4 ; break ;
      case  8 : top = &tet8 ; break ;
      case 10 : top = &tet10 ; break ;
      default:
        //throw RuntimeError() << "Invalid nnode specified" << std::endl ;// << StackTrace;
        throw std::runtime_error( "Invalid nnode specified") ; // << std::endl ;// << StackTrace;
      }

      return top ;
    }



    /*--------------------------------------------------------------------*/
    /**
     *  Wedge Element with 6 Nodes
     *  Pentahedral Element with 6 Nodes
     *  3-Sided Prism Element with 6 Nodes
     *
     *                                   PARENT Linear 6-Node Wedge Nodes
     *                       5           (SPACE_DIM = 3!)
     *                     . o
     *                    . / \
     *                   . /   \         Face_Quad_4_3D()  0-1-4-3
     *                  . /     \        Face_Quad_4_3D()  1-2-5-4 (error old was 1-4-5-2)
     *                 . /       \       Face_Quad_4_3D()  0-2-5-3
     *              2 . o---------o 4    Face_Tri_3_3D()   0-2-1 (error old was 0-1-2 )
     *               o . 3      .        Face_Tri_3_3D()   3-4-5
     *              /.\        .
     *             /.  \      .
     *            /.    \    .
     *           /.      \  .
     *          o---------o
     *         0           1
     *
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
     */

    /*  New ref topo info Wedge6
     *  ----------------------
     *
     *  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
     */
#if 0
    template<> RefTopoX RefinementTopologyExtra< shards:: Wedge<6>  > :: refinement_topology = {
      {	0,	0,	0,	0,	1,	{0,	0,	0} },
      {	1,	0,	1,	0,	1,	{1,	0,	0} },
      {	2,	0,	2,	0,	1,	{0,	1,	0} },
      {	3,	0,	3,	0,	1,	{0,	0,	1} },
      {	4,	0,	4,	0,	1,	{1,	0,	1} },
      {	5,	0,	5,	0,	1,	{0,	1,	1} },
      {	6,	1,	0,	0,	1,	{0.5,	0,	0} },
      {	7,	1,	1,	0,	1,	{0.5,	0.5,	0} },
      {	8,	1,	2,	0,	1,	{0,	0.5,	0} },
      {	9,	1,	6,	0,	1,	{0,	0,	0.5} },
      {	10,	1,	7,	0,	1,	{1,	0,	0.5} },
      {	11,	1,	8,	0,	1,	{0,	1,	0.5} },
      {	12,	1,	3,	0,	1,	{0.5,	0,	1} },
      {	13,	1,	4,	0,	1,	{0.5,	0.5,	1} },
      {	14,	1,	5,	0,	1,	{0,	0.5,	1} },
      {	15,	2,	0,	0,	1,	{0.5,	0,	0.5} },
      {	16,	2,	1,	0,	1,	{0.5,	0.5,	0.5} },
      {	17,	2,	2,	0,	1,	{0,	0.5,	0.5} }

    };
#endif

    /*--------------------------------------------------------------------*/
    /**
     *  Wedge Element with 15 Nodes
     *  Pentahedral Element with 15 Nodes
     *  3-Sided Prism Element with 15 Nodes
     *
     *                          PARENT Quadratic 15-Node Wedge Nodes
     *                          (SPACE_DIM = 3!)
     *
     *                          5
     *                          o
     *                         / \
     *                     14 /   \
     *                       o     o 13
     *                      /       \
     *                   3 /         \
     *                    o-----o-----o
     *                         12      4
     *                  11
     *                   o              Face_Quad_8_3D() 0-6-1-10-4-12-3-9
     *                  . .             Face_Quad_8_3D() 1-10-4-13-5-11-2-7
     *                 .   .            Face_Quad_8_3D() 0-8-2-11-5-14-3-9
     *             17 .     .16         Face_Tri_6_3D()  0-6-1-7-2-8
     *               .       .          Face_Tri_6_3D()  3-12-4-13-5-14
     *              .         .
     *           9 o...........o 10
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
     *   |
     *   | After Refinement:
     *   |
     *
     *          Face #0                     Face #1                      Face #2
     *
     *  3    30   12   31     4      4    32   13    33    5      5    34   14   35     3
     *   o----*----o----*----o        o----*----o----*----o        o----*----o----*----o
     *   |         |         |        |         |         |        |         |         |
     *   |         |         |        |         |         |        |         |         |
     * 51*         *52       *53    53*         *54       *55    55*         *56       *51
     *   |         |         |        |         |         |        |         |         |
     *   |       15|    25   |        |       16|    27   |        |       17|   29    |
     * 9 o----*----o----*----o 10  10 o----*----o----*----o 11  11 o----*----o----*----o 9
     *   |   24    |         |        |   26    |         |        |   28    |         |
     *   |         |         |        |         |         |        |         |         |
     * 45*       46*         *47    47*       48*       49*      49*       50*         *45  ^
     *   |         |         |        |         |         |        |         |         |    |
     *   |         |         |        |         |         |        |         |         |    |
     *   o----*----o----*----o        o----*----o----*----o        o----*----o----*----o    #
     *  0    18    6   19     1      1    20    7   21     2      2    22    8   23     0
     *   #                            #
     *    \                            \
     *     \                            \
     *      -->                          -->
     *
     *
     *             Face #4                          Face #5
     *
     *   2     21     7    20      1      5     33    13     32     4
     *    o-----*-----o-----*-----o        o-----*-----o-----*-----o
     *     \         / \         /          \         / \         /
     *      \       /   \       /            \       /   \       /
     *     22*   36*     *38   *19          34*   42*     *44   *31
     *        \   /       \   /                \   /       \   /
     *         \ /   37    \ /                  \ /   43    \ /
     *        8 o-----*-----o 6               14 o-----*-----o 12
     *           \         /                      \         /
     *            \       /                        \       /
     *       ^   23*     *18                      35*     *30    ^
     *        \     \   /                            \   /      /
     *         \     \ /                              \ /      /
     *          #     o                                o      #
     *                0                                3
     *
     *
     *
     *   CHILD Wedge15 node tables:
     * |
     * | static const UInt child_0[] = {  0, 6, 8, 9, 15, 17, 18, 37, 28, 45, 46, 50, 24, 40, 29 };
     * | static const UInt child_1[] = {  6, 1, 7, 15, 10, 16, 19, 20, 38, 46, 47, 48, 25, 26, 41 };
     * | static const UInt child_2[] = {  8, 7, 2, 17, 16, 11, 36, 21, 22, 50, 48, 49, 39, 27, 28 };
     * | static const UInt child_3[] = {  9, 15, 17, 3, 12, 14, 24, 40, 29, 51, 52, 56, 30, 43, 35 };
     * | static const UInt child_4[] = { 15, 10, 16, 12, 4, 13, 25, 26, 41, 52, 53, 54, 31, 32, 44 };
     * | static const UInt child_5[] = { 17, 16, 11, 14, 13, 5, 39, 27, 28, 56, 54, 55, 42, 33, 34 };
     * | static const UInt child_6[] = {  7, 8, 6, 16, 17, 15, 36, 37, 38, 48, 50, 46, 39, 40, 41 };
     * | static const UInt child_7[] = { 16, 17, 15, 13, 14, 12, 39, 40, 41, 54, 56, 52, 42, 43, 44 };
     * |
     *
     *   Refined Wedge15 Edge node tables:
     * |
     * | static const UInt edge_0[] = { 0, 1,  6, 18, 19 };
     * | static const UInt edge_1[] = { 1, 2,  7, 20, 21 };
     * | static const UInt edge_2[] = { 2, 0,  8, 22, 23 };
     * | static const UInt edge_3[] = { 3, 4, 12, 30, 31 };
     * | static const UInt edge_4[] = { 4, 5, 13, 32, 33 };
     * | static const UInt edge_5[] = { 5, 3, 14, 34, 35 };
     * | static const UInt edge_6[] = { 0, 3,  9, 45, 51 };
     * | static const UInt edge_7[] = { 1, 4, 10, 47, 53 };
     * | static const UInt edge_8[] = { 2, 5, 11, 49, 55 };
     *
     *   Refined Wedge15 Face node tables:
     * |
     * | static const UInt face_0[] = {0, 1, 4, 3, 6, 10, 12, 9,  15, 18, 19, 47, 53, 31, 30, 51, 45, 46, 25, 52, 24};
     * | static const UInt face_1[] = {1, 2, 5, 4, 7, 11, 13, 10, 16, 20, 21, 49, 55, 33, 32, 53, 47, 48, 27, 54, 26};
     * | static const UInt face_2[] = {0, 3, 5, 2, 9, 14, 11, 8,  17, 45, 51, 35, 34, 55, 49, 22, 23, 29, 56, 28, 50};
     * | static const UInt face_3[] = {0, 2, 1, 8, 7, 6,    23, 22, 21, 20, 19, 18, 37, 36, 38 };
     * | static const UInt face_4[] = {3, 4, 5, 12, 13, 14, 30, 31, 32, 33, 34, 35, 43, 44, 42 };
     * |
     *
     **/

    /*--------------------------------------------------------------------*/
    // Pentahedrons with 6 or 15 nodes (3D solid).
    // Wedges with 6 or 15 nodes (3D solid).
    //   eclass ==  SOLID        = > sdim ==  3 && mdim ==  3
    const MeshObjTopology * wedge( UInt nnode )
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology wedge6(shards::getCellTopologyData<shards::Wedge<6> >());
      static MeshObjTopology wedge15(shards::getCellTopologyData<shards::Wedge<15> >());
      static MeshObjTopology wedge18(shards::getCellTopologyData<shards::Wedge<18> >());

      static bool first = true ;

      if ( first ) {

        first = false ;

        { // Children and child nodes
          static const_top_ptr wedge6_child[]  =
            {&wedge6, &wedge6, &wedge6, &wedge6,
             &wedge6, &wedge6, &wedge6, &wedge6 };

          static const_top_ptr wedge15_child[]  =
            {&wedge15, &wedge15, &wedge15, &wedge15,
             &wedge15, &wedge15, &wedge15, &wedge15};

          static const UInt   child_0[]  =
            {  0, 6, 8, 9, 15, 17, 18, 37, 23, 45, 46, 50, 24, 40, 29 , EUA};
          static const UInt   child_1[]  =
            {  6, 1, 7, 15, 10, 16, 19, 20, 38, 46, 47, 48, 25, 26, 41 , EUA};
          static const UInt   child_2[]  =
            {  8, 7, 2, 17, 16, 11, 36, 21, 22, 50, 48, 49, 39, 27, 28 , EUA};
          static const UInt   child_3[]  =
            {  9, 15, 17, 3, 12, 14, 24, 40, 29, 51, 52, 56, 30, 43, 35 , EUA};
          static const UInt   child_4[]  =
            { 15, 10, 16, 12, 4, 13, 25, 26, 41, 52, 53, 54, 31, 32, 44 , EUA};
          static const UInt   child_5[]  =
            { 17, 16, 11, 14, 13, 5, 39, 27, 28, 56, 54, 55, 42, 33, 34 , EUA};
          static const UInt   child_6[]  =
            {  7, 8, 6, 16, 17, 15, 36, 37, 38, 48, 50, 46, 39, 40, 41 , EUA};
          static const UInt   child_7[]  =
            { 16, 17, 15, 13, 14, 12, 39, 40, 41, 54, 56, 52, 42, 43, 44 , EUA};
          static const UInt * child_node_table[]  =
            { child_0 , child_1 , child_2 , child_3 ,
              child_4 , child_5 , child_6 , child_7  };



          // Edge topology and node tables including edges' child-nodes

          static const UInt   edge_0[] = { 0, 1,  6, 18, 19 , EUA};
          static const UInt   edge_1[] = { 1, 2,  7, 20, 21 , EUA};
          static const UInt   edge_2[] = { 2, 0,  8, 22, 23 , EUA};
          static const UInt   edge_3[] = { 3, 4, 12, 30, 31 , EUA};
          static const UInt   edge_4[] = { 4, 5, 13, 32, 33 , EUA};
          static const UInt   edge_5[] = { 5, 3, 14, 34, 35 , EUA};
          static const UInt   edge_6[] = { 0, 3,  9, 45, 51 , EUA};
          static const UInt   edge_7[] = { 1, 4, 10, 47, 53 , EUA};
          static const UInt   edge_8[] = { 2, 5, 11, 49, 55 , EUA};
          static const UInt * edge_table[]  =
            { edge_0 , edge_1 , edge_2 ,
              edge_3 , edge_4 , edge_5 ,
              edge_6 , edge_7 , edge_8 };

          // Face topology, nodes tables, and edge tables

          static const UInt face_0[] = {0, 1, 4, 3, 6, 10, 12, 9,  15, 18, 19, 47, 53, 31, 30, 51, 45, 46, 25, 52, 24, EUA};
          static const UInt face_1[] = {1, 2, 5, 4, 7, 11, 13, 10, 16, 20, 21, 49, 55, 33, 32, 53, 47, 48, 27, 54, 26, EUA};
          static const UInt face_2[] = {0, 3, 5, 2, 9, 14, 11, 8,  17, 45, 51, 35, 34, 55, 49, 22, 23, 29, 56, 28, 50, EUA};
          static const UInt face_3[] = {0, 2, 1, 8, 7, 6,    23, 22, 21, 20, 19, 18, 37, 36, 38 , EUA};
          static const UInt face_4[] = {3, 4, 5, 12, 13, 14, 30, 31, 32, 33, 34, 35, 43, 44, 42 , EUA};

          static const UInt * face_table[]  =
            { face_0 , face_1 , face_2 , face_3 , face_4 };

          static RefinementTopology wedge6_refinement(&wedge6, 8, wedge6_child, 18, child_node_table, 9, edge_table, 5, face_table, 0, NULL, NULL, false);
          static RefinementTopology wedge15_refinement(&wedge15, 8, wedge15_child, 57, child_node_table, 9, edge_table, 5, face_table, 0, NULL, NULL, false);
        }
      }

      VERIFY_TRUE( nnode ==  6 || nnode ==  15 );

      MeshObjTopology * const top = nnode ==  6 ? &wedge6 : &wedge15 ;

      return top ;
    }


    /**
     *                               PARENT Linear: 5-Node Pyramid Nodes
     *
     *     3                   2
     *      o-----------------o
     *      |\               /|
     *      | \             / |
     *      |  \           /  |
     *      |   \         /   |
     *      |    \       /    |
     *      |     \     /     |
     *      |      \   /      |
     *      |       \ /       |
     *      |        o 4      |
     *      |       / \       |
     *      |      /   \      |
     *      |     /     \     |
     *      |    /       \    |
     *      |   /         \   |
     *      |  /           \  |
     *      | /             \ |
     *      |/               \|
     *      o-----------------o
     *     0                   1
     *
     *                               CHILDREN 5-Node Pyramid Nodes
     *                               6 new Pyramids, 4 new Tets
     *
     *      3--------7--------2
     *      |\      / \      /|
     *      | \    /   \    / |
     *      |  \  /     \  /  |          static const UInt   child_0[] = {  0, 5, 13, 8, 9, EUA };
     *      |   \/       \/   |          static const UInt   child_1[] = {  5, 1, 6, 13, 10, EUA };
     *      |   12-------11   |          static const UInt   child_2[] = {  13, 6, 2, 7, 11, EUA };
     *      |  /| \     /| \  |          static const UInt   child_3[] = {  8, 13, 7, 3, 12, EUA };
     *      | / |  \   / |  \ |          static const UInt   child_4[] = {  9, 12, 11, 10, 13, EUA };
     *      |/  |   \ /  |   \|          static const UInt   child_5[] = {  9, 10, 11, 12, 4, EUA };
     *      8   |   4,13 |    6          static const UInt   child_6[] = {  13, 10, 9, 5, EUA };
     *      |\  |   / \  |   /|          static const UInt   child_7[] = {  13, 11, 10, 6, EUA };
     *      | \ |  /   \ |  / |          static const UInt   child_8[] = {  13, 12, 11, 7, EUA };
     *      |  \| /     \| /  |          static const UInt   child_9[] = {  13, 9, 12, 8, EUA };
     *      |   9/______10/   |
     *      |   /\       /\   |
     *      |  /  \     /  \  |
     *      | /    \   /    \ |
     *      |/      \ /      \|
     *      o--------5--------o
     *     0                   1
     *
     *                               PARENT Quadratic 13- & 14-Node Pyramid Nodes
     *
     *     3                   2
     *      o--------7--------o
     *      |\               /|
     *      | \             / |
     *      |  \           /  |
     *      |   \         /   |
     *      |   12       11   |
     *      |     \     /     |
     *      |      \   /      |
     *      |       \ /       |
     *      8        o 4,13   6  Node 13 is centroid node of the Quad face
     *      |       / \       |
     *      |      /   \      |
     *      |     /     \     |
     *      |    9      10    |
     *      |   /         \   |
     *      |  /           \  |
     *      | /             \ |
     *      |/               \|
     *      o--------5--------o
     *     0                   1
     *
     *
     *
     *                               CHILDREN 13, and 14-Node Pyramid Nodes
     *                               6 new Pyramids, 4 new Tets
     *
     *      3---19---7---18---2
     *      |\      / \      /|          Note: the 99's are filled in programatically by the RefinementTopoolgy code.
     *      | 42  43   44  45 |
     *      |  \  /     \  /  |          static const UInt   children[][14] = {{  0, 5, 13, 8, 9,    99,99,99,99, 99,99,99,99, 99, EUA},
     *     20   \/       \/  17           {  5, 1, 6, 13, 10,   99,99,99,99, 99,99,99,99, 99, EUA},
     *      |   12--41---11   |           {  13, 6, 2, 7, 11,   99,99,99,99, 99,99,99,99, 99, EUA},
     *      |  /| \     /| \  |           {  8, 13, 7, 3, 12,   99,99,99,99, 99,99,99,99, 99, EUA},
     *      |37 | 38  39 | 40 |           {  9, 12, 11, 10, 13, 99,99,99,99, 99,99,99,99, 99, EUA},
     *      |/  |   \ /  |   \|           {  9, 10, 11, 12, 4,  99,99,99,99, 99,99,99,99, 99, EUA},
     *      8  35   4,13 36   6           {  13, 10, 9, 5,      99,99,99,99, 99,99,99,99, 99, EUA},
     *      |\  |   / \  |   /|           {  13, 11, 10, 6,     99,99,99,99, 99,99,99,99, 99, EUA},
     *      |31 | 32  33 |  34|           {  13, 12, 11, 7,     99,99,99,99, 99,99,99,99, 99, EUA},
     *      |  \| /     \| /  |           {  13, 9, 12, 8,      99,99,99,99, 99,99,99,99, 99, EUA} };
     *      |   9/__30__10/   |
     *      21  /\       /\  16
     *      |  /  \     /  \  |
     *      | 26   27  28  29 |          Edges (note: these may not be in edge-orientation order, but not necessary, just used to build child tables)
     *      |/      \ /      \|
     *      o---14---5---15---o          static const UInt edges[][3] = {
     *     0                   1            {0,5,14},{5,1,15},{1,6,16},{6,2,17},{2,7,18},{7,3,19},{3,8,20},{8,0,21},
     *                                      {5,13,22},{8,13,23},{13,6,24},{13,7,25},
     *                                      {0,9,26},{5,9,27},{5,10,28},{1,10,29},{9,10,30},{8,9,31},{9,4,32},{10,4,33},{6,10,34}
     *                                      {9,12,35},{10,11,36},{8,12,37},{4,12,38},{4,11,39},{6,11,40},{12,11,41},{3,12,42},{7,12,43},{7,11,44},{2,11,45} };
     *
     *
     *    >>OLD, incorrect<<         Additional nodes on edges from the 4 Tets
     *                               Left here for reference - there are holes in the numbering, but that's ok, these are really only
     *      3--------7--------2      labels, anyways, so we can live with the holes.
     *      |\      /|\      /|
     *      | \    / | \    / |
     *      |  \  /  |  \  /  |         add_edges[] = {
     *      |   \/   53  \/   |            {13,5,22},{13,9,47},{13,10,48},
     *      |   12   |   11   |            {13,8,23},{13,6,24},{13,12,51},{13,11,52},{13,7 ,53} };
     *      |  /  \  |  /  \  |
     *      | /   51 | 52   \ |
     *      |/      \|/      \|
     *      8---49--4,13--50--6
     *      |\      /|\      /|
     *      | \   47 | 48   / |
     *      |  \  /  |  \  /  |
     *      |   9/  22  10/   |
     *      |   /\   |   /\   |
     *      |  /  \  |  /  \  |
     *      | /    \ | /    \ |
     *      |/      \|/      \|
     *      o--------5--------o
     *     0                   1
     *
     *    renum 1: 53->25, 50->24, 46->22, 49->23
     *    Holes:  45,*,47,48,*,*,51,52,*,54,55,56,57
     *    New:    45, ,46,47,    48,49,  50,51,52,53
     *    Renumber: 47->46,48->47,49->
     *
     *    >> NEW, corrected <<
     *
     *                               Additional nodes on edges from the 4 Tets
     *
     *      3--------7--------2
     *      |\      /|\      /|
     *      | \    / | \    / |
     *      |  \  /  |  \  /  |         add_edges[] = {
     *      |   \/   25  \/   |            {13,5,22},{13,9,46},{13,10,47},
     *      |   12   |   11   |            {13,8,23},{13,6,24},{13,12,48},{13,11,49},{13,7 ,25} };
     *      |  /  \  |  /  \  |
     *      | /   48 | 49   \ |
     *      |/      \|/      \|
     *      8---23--4,13--24--6
     *      |\      /|\      /|
     *      | \   46 | 47   / |
     *      |  \  /  |  \  /  |
     *      |   9/  22  10/   |
     *      |   /\   |   /\   |
     *      |  /  \  |  /  \  |
     *      | /    \ | /    \ |
     *      |/      \|/      \|
     *      o--------5--------o
     *     0                   1
     *
     *
     * |    Bottom Quadrilateral face
     * |
     * |    3---19----7---18----2        Faces (quad faces on bottom):
     * |    |         |         |
     * |    |         |         |        static const UInt face_centroids[][5] = {
     * |    20  52   25   53   17          {0,5,13,8,22},{5,1,6,13,46},{8,13,7,3,47},{13,6,2,7,23} };
     * |    |         |         |
     * |    |         |         |
     * |    8----23--13----24---6        Original Edges:
     * |    |         |         |
     * |    |         |         |        static const UInt original_edges[][5] = {
     * |    21  50   22   51   16          {0,1,5,14,15},{1,2,6,16,17},{2,3,7,18,19},{3,0,8,20,21},
     * |    |         |         |          {0,4,9,26,32},{1,4,10,29,33},{2,4,11,45,39},{3,4,12,42,38} };
     * |    |         |
     * |    0---14----5----15-- 1        Original Faces:
     * |
     *                                   static const UInt original_quad_faces[][] = {};
     *
     *  quad_face_nodes = {0,3,2,1, 8,7,6,5, 13, 21,20,19,18,17,16,15,14, 23,25,24,22, 50, 52, 53, 51}
     *
     */
    template<unsigned M>
    static void print_table(unsigned N, UInt table[][M], std::string name)
    {
      std::cout << "static UInt " << name << "[" << N << "][" << M << "] = {\n";
      for (unsigned ic=0; ic < N; ic++)
        {
          UInt * child = table[ic];
          for (unsigned jc=0; jc < M; jc++)
            {
              bool end = (child[jc] == EUA || jc== M-1);
              std::cout << (jc==0?"{":"") << (child[jc] == EUA?"EUA":std::to_string(child[jc])) << (end?(ic==(N-1)?"}\n":"},\n"):", ");
              if (end) break;
            }
        }
      std::cout << "};" << std::endl;
    }

    // Pyramids with 5 or 13 nodes (3D solid).
    //   eclass ==  SOLID        = > sdim ==  3 && mdim ==  3
    const MeshObjTopology *
    pyramid(
            UInt                  nnode )
    {
      /* %TRACE[NONE]% */  /* %TRACE% */
      static MeshObjTopology pyramid5(shards::getCellTopologyData<shards::Pyramid<5> >());
      static MeshObjTopology pyramid13(shards::getCellTopologyData<shards::Pyramid<13> >());
      static MeshObjTopology tet4(shards::getCellTopologyData<shards::Tetrahedron<4> >());
      static MeshObjTopology tet10(shards::getCellTopologyData<shards::Tetrahedron<10> >());

      static bool first = true ;

      if ( first ) {

        first = false ;

        { // Children and child nodes

          static const_top_ptr pyramid5_child[]  =
            {&pyramid5, &pyramid5, &pyramid5, &pyramid5, &pyramid5, &pyramid5,
             &tet4, &tet4, &tet4, &tet4};

          static const_top_ptr pyramid13_child[]  =
            {&pyramid13, &pyramid13, &pyramid13, &pyramid13, &pyramid13, &pyramid13,
             &tet10, &tet10, &tet10, &tet10
            };

          //           static const UInt   child_0[] = {  0, 5, 13, 8, 9 , EUA};
          //           static const UInt   child_1[] = {  5, 1, 6, 13, 10 , EUA};
          //           static const UInt   child_2[] = {  13, 6, 2, 7, 11 , EUA};
          //           static const UInt   child_3[] = {  8, 13, 7, 3, 12 , EUA};
          //           static const UInt   child_4[] = {  9, 12, 11, 10, 13 , EUA};
          //           static const UInt   child_5[] = {  9, 10, 11, 12, 4 , EUA};
          //           static const UInt   child_6[] = {  13, 10, 9, 5 , EUA};
          //           static const UInt   child_7[] = {  13, 11, 10, 6 , EUA};
          //           static const UInt   child_8[] = {  13, 12, 11, 7 , EUA};
          //           static const UInt   child_9[] = {  13, 9, 12, 8 , EUA};

          //Note: the 99's are filled in programatically

          typedef const UInt * const *ChildNodeTableConst;

          static UInt child_node_table_0[][15] = {
            {  0, 5, 13, 8, 9,    99,99,99,99, 99,99,99,99, 99, EUA},
            {  5, 1, 6, 13, 10,   99,99,99,99, 99,99,99,99, 99, EUA},
            {  13, 6, 2, 7, 11,   99,99,99,99, 99,99,99,99, 99, EUA},
            {  8, 13, 7, 3, 12,   99,99,99,99, 99,99,99,99, 99, EUA},
            {  9, 12, 11, 10, 13, 99,99,99,99, 99,99,99,99, 99, EUA},
            {  9, 10, 11, 12, 4,  99,99,99,99, 99,99,99,99, 99, EUA},
            {  13, 10, 9, 5,      99,99,99,99, 99,99,99,99, 99, EUA},
            {  13, 11, 10, 6,     99,99,99,99, 99,99,99,99, 99, EUA},
            {  13, 12, 11, 7,     99,99,99,99, 99,99,99,99, 99, EUA},
            {  13, 9, 12, 8,      99,99,99,99, 99,99,99,99, 99, EUA}
          };
          (void)child_node_table_0;

          // Generated from above table - see below
          static UInt cnode[10][15] = {
            {0, 5, 13, 8, 9, 14, 22, 23, 21, 26, 27, 46, 31, EUA},
            {5, 1, 6, 13, 10, 15, 16, 24, 22, 28, 29, 34, 47, EUA},
            {13, 6, 2, 7, 11, 24, 17, 18, 25, 49, 40, 45, 44, EUA},
            {8, 13, 7, 3, 12, 23, 25, 19, 20, 37, 48, 43, 42, EUA},
            {9, 12, 11, 10, 13, 35, 41, 36, 30, 46, 48, 49, 47, EUA},
            {9, 10, 11, 12, 4, 30, 36, 41, 35, 32, 33, 39, 38, EUA},
            {13, 10, 9, 5, 47, 30, 46, 22, 28, 27, EUA},
            {13, 11, 10, 6, 49, 36, 47, 24, 40, 34, EUA},
            {13, 12, 11, 7, 48, 41, 49, 25, 43, 44, EUA},
            {13, 9, 12, 8, 46, 35, 48, 23, 31, 37, EUA}
          };

          static const UInt * child_node_table[]  =
            {cnode[0], cnode[1], cnode[2], cnode[3], cnode[4], cnode[5], cnode[6], cnode[7], cnode[8], cnode[9] };

          // Edges (note: these may not be in edge-orientation order, but not necessary, just used to build child tables)

          static const UInt edges[][3] = { // total 40 (32 + 8 additional tet edges)
            {0,5,14},{5,1,15},{1,6,16},{6,2,17},{2,7,18},{7,3,19},{3,8,20},{8,0,21},   // outer edges = 8
            {5,13,22},{8,13,23},{13,6,24},{13,7,25},  // quad face edges = 4
            {0,9,26},{5,9,27},{5,10,28},{1,10,29},{9,10,30},{8,9,31},{9,4,32},{10,4,33},{6,10,34},  // tri-face edges = 20
            {9,12,35},{10,11,36},{8,12,37},{4,12,38},{4,11,39},{6,11,40},{12,11,41},{3,12,42},{7,12,43},{7,11,44},{2,11,45},

            // additional edges from tets
            {13,5,22},{13,9,46},{13,10,47},
            {13,8,23},{13,6,24},{13,12,48},{13,11,49},{13,7 ,25}

          };
          //static int n_edge_table = 40;
          static UInt n_edge_table = sizeof(edges)/(3*sizeof(UInt));
          if (n_edge_table != 40) throw std::logic_error("bad table");

          static bool did_compute_child_nodes = true; // toggle this to regenerate the child_node_table
          if (!did_compute_child_nodes)
            {
              did_compute_child_nodes = true;
              for (unsigned ic=0; ic < 10; ic++)
                {
                  UInt * child = child_node_table_0[ic];
                  const CellTopologyData * topo_data = shards::getCellTopologyData<shards::Pyramid<5> >();
                  if (ic > 5) {
                    topo_data = shards::getCellTopologyData<shards::Tetrahedron<4> >();
                  }
                  shards::CellTopology topo(topo_data);
                  int nv = topo_data->vertex_count;
                  unsigned n_edge = topo.getEdgeCount();
                  for (unsigned i_edge = 0; i_edge < n_edge; i_edge++)
                    {
                      int i0 = child[topo_data->edge[i_edge].node[0]];
                      int i1 = child[topo_data->edge[i_edge].node[1]];

                      bool found = false;
                      for (unsigned j_edge = 0; j_edge < n_edge_table; j_edge++)
                        {
                          int j0 = edges[j_edge][0];
                          int j1 = edges[j_edge][1];
                          int quadratic_edge_node = edges[j_edge][2];
                          //std::cout << " i0= " << i0 << " i1= " << i1 << " j0= " << j0 << " j1= " << j1 << std::endl;
                          if ( (i0 == j0 && i1 == j1) || (i0 == j1 && i1 == j0))
                            {
                              found = true;
                              child[nv  + i_edge] = quadratic_edge_node;
                            }
                        }
                      if (!found) throw std::logic_error("bad pyramid edge table");
                    }
                  child[nv+n_edge] = EUA;
                }
              print_table(10, child_node_table_0, "cnode");

            }

          // Edge topology and node tables including edges' child-nodes

          static const UInt edge_t[][6] = {
            {0, 1, 5, 14, 15, EUA},
            {1, 2, 6, 16, 17, EUA},
            {2, 3, 7, 18, 19, EUA},
            {3, 0, 8, 20, 21, EUA},
            {0, 4, 9, 26, 32, EUA},
            {1, 4, 10, 29, 33, EUA},
            {2, 4, 11, 45, 39, EUA},
            {3, 4, 12, 42, 38, EUA}
          };

          static const UInt * edge_table[]  = { edge_t[0], edge_t[1], edge_t[2], edge_t[3], edge_t[4], edge_t[5], edge_t[6], edge_t[7] };

          // Face topology, nodes tables, and edge tables

          static const UInt face_0[] = { 0, 1, 4,     5,  10,  9, 14,15,29,33,32,26,27,28,30,EUA};
          static const UInt face_1[] = { 1, 2, 4,     6,  11, 10, 16,17,45,39,33,29,34,40,36,EUA};
          static const UInt face_2[] = { 2, 3, 4,     7,  12, 11, 18,19,42,38,39,45,44,43,41,EUA};
          static const UInt face_3[] = { 3, 0, 4,     8,   9,  12, 20,21,26,32,38,42,31,35,37,EUA};
          //static const UInt face_4[] = { 0, 3, 2, 1,  8,   7,  6,  5, 13 ,14,15,16,17,18,19,20,21, EUA};
          static const UInt face_4[] = {0,3,2,1, 8,7,6,5, 13, 21,20,19,18,17,16,15,14, 23,25,24,22, 50, 52, 53, 51, EUA};
          static const UInt * face_table[]  =
            { face_0 , face_1 , face_2 , face_3 , face_4 };

          static RefinementTopology pyramid5_refinement(&pyramid5, 10, pyramid5_child, 14, child_node_table, 8, edge_table, 5, face_table, 0, NULL, NULL, false);
          static RefinementTopology pyramid13_refinement(&pyramid13, 10, pyramid13_child, 50,child_node_table, 8, edge_table, 5, face_table, 0, NULL, NULL, false);
        }
      }

      VERIFY_TRUE( nnode ==  5 || nnode ==  13 );

      MeshObjTopology * const top = nnode ==  5 ? &pyramid5 : &pyramid13 ;

      return top ;
    }


    void
    bootstrap()
    {
      /* %TRACE[ON]% */  /* %TRACE% */
      static bool first_call = true ;

      if ( first_call ) { // Construct topologies
        first_call = false ;

        // Nodes
        node();

        // Edges
        line(NOT_ELEMENT, 3);

        // Faces
        tri( NOT_ELEMENT, 3);
        tri4(NOT_ELEMENT);
        tri( NOT_ELEMENT, 6);
        quad(NOT_ELEMENT, 4);
        quad(NOT_ELEMENT, 8);
        quad(NOT_ELEMENT, 9);

        // Elements
        point();
        line(ROD, 2);
        line(ROD, 3);
        line(SHELL, 2);
        line(SHELL, 3);
        tri( SHELL, 3);
        tri( SHELL, 6);
        tri4(SHELL);
        quad(SHELL, 4);
        quad(SHELL, 9);

        hex( 8);
        hex(20);
        hex(27);
        tet( 4);
        tet( 8);
        tet(10);
        wedge( 6);
        wedge(15);
        pyramid( 5);
        pyramid(13);
      }
    }

  } // namespace StdMeshObjTopologies
} // namespace Elem
} // namespace percept

