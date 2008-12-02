/*------------------------------------------------------------------------*/
/*                  shards : Shared Discretization Tools                  */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/

#ifndef Shards_BasicTopologies_hpp
#define Shards_BasicTopologies_hpp

#include <Shards_CellTopologyTraits.hpp>

namespace shards {

/** \addtogroup shards_package_cell_topology
 *  \{
 */
//----------------------------------------------------------------------
 
/** \brief Topological traits: Dimension = 0, Vertices = 0, Nodes = 0. */
struct Node : public CellTopologyTraits<0,0,0>
{
#ifndef DOXYGEN_COMPILE
  typedef Node base ;
#endif /* DOXYGEN_COMPILE */
};
 
/** \brief  Singleton for Node topology */
template<> const CellTopologyData * getCellTopologyData< Node >();

//----------------------------------------------------------------------

/** \brief Topological traits: Dimension = 1, Vertices = 1, Nodes = 1. */
struct Particle : public CellTopologyTraits<1,1,1>
{
#ifndef DOXYGEN_COMPILE
  typedef Particle base ;
#endif /* DOXYGEN_COMPILE */
};

/** \brief  Singleton for Particle topology */
template<> const CellTopologyData * getCellTopologyData< Particle >();

//----------------------------------------------------------------------

/** \brief Topological traits: Dimension = 1, Vertices = 2, Nodes = 2 or 3.
 *
 *  A line has local node ordinals as follows: <br>
 *  [0]--------[2]--------[1]  ---> Positive direction
 */
template< unsigned NodeCount = 2 > struct Line {};

/** \brief  Singleton for line topology with two nodes.  */
template<> const CellTopologyData * getCellTopologyData< Line<2> >();

/** \brief  Singleton for line topology with three nodes.  */
template<> const CellTopologyData * getCellTopologyData< Line<3> >();

//----------------------------------------------------------------------

/** \brief  Topological traits: Dimension = 2, Edges = 1, Vertices = 2,
 *          and Nodes = 2 or 3.
 *
 *  \see shards::Line
 */
template< unsigned NodeCount = 2 > struct Beam {};

/** \brief  Singleton for beam topology with two nodes.  */
template<> const CellTopologyData * getCellTopologyData< Beam<2> >();

/** \brief  Singleton for beam topology with three nodes.  */
template<> const CellTopologyData * getCellTopologyData< Beam<3> >();

//----------------------------------------------------------------------

/** \brief  Topological traits: Dimension = 2, Edges = 2, Vertices = 2,
 *          and Nodes = 2 or 3.
 *
 *  \see shards::Line
 */
template< unsigned NodeCount = 2 > struct ShellLine {};

/** \brief  Singleton for shell-line topology with two nodes.  */
template<> const CellTopologyData * getCellTopologyData< ShellLine<2> >();

/** \brief  Singleton for shell-line topology with three nodes.  */
template<> const CellTopologyData * getCellTopologyData< ShellLine<3> >();

//----------------------------------------------------------------------

/** \brief  Topological traits: Dimension = 2, Edges = 3, Vertices = 3,
 *          and Nodes = 3 or 6.
 *
 * <PRE>
 *  Triangle node numbering, up to 6 nodes:
 *
 *                  2
 *                  o
 *                 / \
 *                /   \
 *               /     \
 * Edge #2    5 o       o 4   Edge #1
 *             /         \
 *            /           \
 *           /             \
 *          o-------o-------o
 *         0        3        1
 *
 *                Edge #0
 *
 * </PRE>
 */
template< unsigned NodeCount = 3 > struct Triangle {};

/** \brief  Return CellTopologyData singleton for the Triangle<3> */
template<> const CellTopologyData * getCellTopologyData< Triangle<3> >();

/**  \brief  Return CellTopologyData singleton for the Triangle<6> */
template<> const CellTopologyData * getCellTopologyData< Triangle<6> >();

/**  \brief  Return CellTopologyData singleton for the Triangle<4> (Face of Tet8) */
template<> const CellTopologyData * getCellTopologyData< Triangle<4> >();

//----------------------------------------------------------------------

/** \brief  Topological traits: Dimension = 3, Sides = 2, Edges = 3,
 *          Vertices = 3, and Nodes = 3 or 6.
 *
 *  \see shards::Triangle
 */
template< unsigned NodeCount = 3 > struct ShellTriangle {};

/** \brief  Return CellTopologyData singleton for the ShellTriangle<3> */
template<> const CellTopologyData * getCellTopologyData< ShellTriangle<3> >();

/** \brief  Return CellTopologyData singleton for the ShellTriangle<6> */
template<> const CellTopologyData * getCellTopologyData< ShellTriangle<6> >();

//----------------------------------------------------------------------
/** \brief Topological traits: Dimension = 2, Edges = 4, Vertices = 4,
 *         and Nodes = 4, 8, or 9.
 *
 * Conventional numbering quadrilateral with up to nine-nodes
 *
 *                 Edge #2
 *
 *            3        6        2
 *             o-------o-------o
 *             |               |
 *             |               |
 *             |       8       |
 *  Edge #3  7 o       o       o 5  Edge #1
 *             |               |
 *             |               |
 *             |               |
 *             o-------o-------o
 *            0        4        1
 *
 *                  Edge #0
 *
 */
template< unsigned NodeCount = 4 > struct Quadrilateral {};

/** \brief  Return CellTopologyData singleton for the Quadrilateral<4> */
template<> const CellTopologyData * getCellTopologyData< Quadrilateral<4> >();

/** \brief  Return CellTopologyData singleton for the Quadrilateral<8> */
template<> const CellTopologyData * getCellTopologyData< Quadrilateral<8> >();

/** \brief  Return CellTopologyData singleton for the Quadrilateral<9> */
template<> const CellTopologyData * getCellTopologyData< Quadrilateral<9> >();

//----------------------------------------------------------------------

/** \brief Topological traits: Dimension = 2, Sides = 2, Edges = 4,
 *         Vertices = 4, and Nodes = 4, 8, or 9.
 *
 *  \see shards::Quadrilateral
 */
template< unsigned NodeCount = 4 > struct ShellQuadrilateral {};

/** \brief  Return CellTopologyData singleton for the ShellQuadrilateral<4> */
template<> const CellTopologyData * getCellTopologyData< ShellQuadrilateral<4> >();

/** \brief  Return CellTopologyData singleton for the ShellQuadrilateral<8> */
template<> const CellTopologyData * getCellTopologyData< ShellQuadrilateral<8> >();

/** \brief  Return CellTopologyData singleton for the ShellQuadrilateral<9> */
template<> const CellTopologyData * getCellTopologyData< ShellQuadrilateral<9> >();

//----------------------------------------------------------------------

/** \brief Topological traits: Dimension = 3, Sides = 4, Edges = 6,
 *         Vertices = 4, and Nodes = 4 or 10.
 */
template< unsigned NodeCount = 4 > struct Tetrahedron ;

/** \brief  Return CellTopologyData singleton for the Tetrahedron<4> */
template<> const CellTopologyData * getCellTopologyData< Tetrahedron<4> >();

/** \brief  Return CellTopologyData singleton for the Tetrahedron<10> */
template<> const CellTopologyData * getCellTopologyData< Tetrahedron<10> >();

/** \brief  Return CellTopologyData singleton for the Tetrahedron<8> */
template<> const CellTopologyData * getCellTopologyData< Tetrahedron<8> >();

//----------------------------------------------------------------------

/** \brief Topological traits: Dimension = 3, Sides = 5, Edges = 8,
 *         Vertices = 5, and Nodes = 5, 13, or 14.
 */
template< unsigned NodeCount = 5 > struct Pyramid {};

/**  \brief  Return CellTopologyData singleton for the Pyramid<5> */
template<> const CellTopologyData * getCellTopologyData< Pyramid<5> >();

/** \brief  Return CellTopologyData singleton for the Pyramid<13> */
template<> const CellTopologyData * getCellTopologyData< Pyramid<13> >();

/** \brief  Return CellTopologyData singleton for the Pyramid<14> */
template<> const CellTopologyData * getCellTopologyData< Pyramid<14> >();

//----------------------------------------------------------------------

/** \brief Topological traits: Dimension = 3, Sides = 5, Edges = 9,
 *         Vertices = 6, and Nodes = 6, 15, or 18.
 */
template< unsigned NodeCount = 6 > struct Wedge {};

/** \brief  Return CellTopologyData singleton for the Wedge<6> */
template<> const CellTopologyData * getCellTopologyData< Wedge<6> >();

/** \brief  Return CellTopologyData singleton for the Wedge<15> */
template<> const CellTopologyData * getCellTopologyData< Wedge<15> >();

/** \brief  Return CellTopologyData singleton for the Wedge<18> */
template<> const CellTopologyData * getCellTopologyData< Wedge<18> >();

//----------------------------------------------------------------------

/** \brief Topological traits: Dimension = 3, Sides = 6, Edges = 12,
 *         Vertices = 8, and Nodes = 8, 20, or 27.
 *
 *  <PRE>
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
 *   Quadratic 20-Node Hexahedron node locations:
 *
 *           7         18         6
 *            o--------o---------o
 *           /|                 /|
 *          / |                / |
 *         /  |               /  |
 *      19o   |            17o   |
 *       /  15o             /    o14
 *      /     |            /     |
 *   4 /      | 16        /      |
 *    o---------o--------o 5     |
 *    |       |       10 |       |
 *    |     3 o-------o--|-------o 2
 *    |      /           |      /
 *    |     /            |     /
 *  12o    /             o13  /
 *    |   o11            |   o9
 *    |  /               |  /
 *    | /                | /
 *    |/                 |/
 *    o---------o--------o
 *   0          8         1
 *
 *
 *   Quadratic 27-Node Hexahedron additional node locations:
 *
 *
 *            x--------x---------x
 *           /|                 /|
 *          / |                / |
 *         /  |   22          /  |
 *        x   |    o         x   |
 *       /    x       o26   /    x     Node #20 is at centroid of element
 *      /     |            /     |
 *     /      |           /      |     "2D surface" containing nodes
 *    x---------x--------x       |      0,1,5,4 has node 25 at center....
 *    | 23o   |          |   o24 |
 *    |       x-------x--|-------x
 *    |      /           |      /
 *    |     /  25        |     /
 *    x    /    o        x    /
 *    |   x        o21   |   x
 *    |  /               |  /
 *    | /                | /
 *    |/                 |/
 *    x---------x--------x
 * </PRE>
 */
template< unsigned NodeCount = 8 > struct Hexahedron {};

/** \brief  Return CellTopologyData singleton for the Hexahedron<8> */
template<> const CellTopologyData * getCellTopologyData< Hexahedron<8> >();

/** \brief  Return CellTopologyData singleton for the Hexahedron<20> */
template<> const CellTopologyData * getCellTopologyData< Hexahedron<20> >();

/** \brief  Return CellTopologyData singleton for the Hexahedron<27> */
template<> const CellTopologyData * getCellTopologyData< Hexahedron<27> >();

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

template<> struct Line<2> : public CellTopologyTraits<1,2,2>
{ typedef Line<2> base ; };

template<> struct Line<3> : public CellTopologyTraits<1,2,3>
{ typedef Line<2> base ; };

// Beam is a line with one edge:

typedef
  MakeTypeList< IndexList< 0 , 1 , 2 > >::type BeamEdgeNodeMap ;

template<> struct Beam<2> : public
  CellTopologyTraits< 2 , 2 , 2 ,
                      MakeTypeList< Line<2> >::type ,
                      BeamEdgeNodeMap >
{ typedef Beam<2> base ; };

template<> struct Beam<3> : public
  CellTopologyTraits< 2 , 2 , 3 ,
                      MakeTypeList< Line<3> >::type ,
                      BeamEdgeNodeMap >
{ typedef Beam<2> base ; };

// Shell-line has two edges:

typedef
  MakeTypeList< IndexList< 0 , 1 , 2 > ,
                IndexList< 1 , 0 , 2 > >::type
    ShellLineEdgeNodeMap ;

template<> struct ShellLine<2> : public
  CellTopologyTraits< 2 , 2 , 2 ,
                      MakeTypeList< Line<2> , Line<2> >::type ,
                      ShellLineEdgeNodeMap >
{ typedef ShellLine<2> base ; };

template<> struct ShellLine<3> : public
  CellTopologyTraits< 2 , 2 , 3 ,
                      MakeTypeList< Line<3> , Line<3> >::type ,
                      ShellLineEdgeNodeMap >
{ typedef ShellLine<2> base ; };

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 , 3 > ,
                IndexList< 1 , 2 , 4 > ,
                IndexList< 2 , 0 , 5 > >::type
    TriangleEdgeNodeMap ;

template<> struct Triangle<3> : public
  CellTopologyTraits< 2 , 3 , 3 ,
                            MakeTypeList< Line<2>  ,
                                          Line<2>  ,
                                          Line<2>  >::type ,
                            TriangleEdgeNodeMap >
{ typedef Triangle<3> base ; };

template<> struct Triangle<6> : public
  CellTopologyTraits< 2 , 3 , 6 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            TriangleEdgeNodeMap >
{ typedef Triangle<3> base ; };

template<> struct Triangle<4> : public
  CellTopologyTraits< 2 , 3 , 4 ,
                            MakeTypeList< Line<2>  ,
                                          Line<2>  ,
                                          Line<2>  >::type ,
                            TriangleEdgeNodeMap >
{ typedef Triangle<3> base ; };

//------------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0,1,2,  3,4,5 > ,
                IndexList< 0,2,1,  5,4,3 >
    >::type ShellTriangleFaceNodeMap ;

template<> struct ShellTriangle<3> : public
  CellTopologyTraits< 3 , 3 , 3 ,
                            MakeTypeList< Line<2>  ,
                                          Line<2>  ,
                                          Line<2>  >::type ,
                            TriangleEdgeNodeMap ,
                            MakeTypeList< Triangle<3>  ,
                                          Triangle<3>  >::type ,
                            ShellTriangleFaceNodeMap >
{ typedef ShellTriangle<3> base ; };

template<> struct ShellTriangle<6> : public
  CellTopologyTraits< 3 , 3 , 6 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            TriangleEdgeNodeMap ,
                            MakeTypeList< Triangle<6>  ,
                                          Triangle<6>  >::type ,
                            ShellTriangleFaceNodeMap >
{ typedef ShellTriangle<3> base ; };

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 ,  4 > ,
                IndexList< 1 , 2 ,  5 > ,
                IndexList< 2 , 3 ,  6 > ,
                IndexList< 3 , 0 ,  7 > >::type
  QuadrilateralEdgeNodeMap ;

template<> struct Quadrilateral<4> : public
  CellTopologyTraits< 2 , 4 , 4 ,
                            MakeTypeList< Line<2>  ,
                                          Line<2>  ,
                                          Line<2>  ,
                                          Line<2>  >::type ,
                            QuadrilateralEdgeNodeMap >
{ typedef Quadrilateral<4> base ; };

template<> struct Quadrilateral<8> : public
  CellTopologyTraits< 2 , 4 , 8 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            QuadrilateralEdgeNodeMap >
{ typedef Quadrilateral<4> base ; };

template<> struct Quadrilateral<9> : public
  CellTopologyTraits< 2 , 4 , 9 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            QuadrilateralEdgeNodeMap >
{ typedef Quadrilateral<4> base ; };

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0,1,2,3,  4,5,6,7,  8 > ,
                IndexList< 0,3,2,1,  7,6,5,4,  8 > >::type
    ShellQuadrilateralFaceNodeMap ;

template<> struct ShellQuadrilateral<4> : public
  CellTopologyTraits< 3 , 4 , 4 ,
                            MakeTypeList< Line<2>  ,
                                          Line<2>  ,
                                          Line<2>  ,
                                          Line<2>  >::type ,
                            QuadrilateralEdgeNodeMap ,
                            MakeTypeList< Quadrilateral<4>  ,
                                          Quadrilateral<4>  >::type ,
                            ShellQuadrilateralFaceNodeMap >
{ typedef ShellQuadrilateral<4> base ; };

template<> struct ShellQuadrilateral<8> : public
  CellTopologyTraits< 3 , 4 , 8 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            QuadrilateralEdgeNodeMap ,
                            MakeTypeList< Quadrilateral<8>  ,
                                          Quadrilateral<8>  >::type ,
                            ShellQuadrilateralFaceNodeMap >
{ typedef ShellQuadrilateral<4> base ; };

template<> struct ShellQuadrilateral<9> : public
  CellTopologyTraits< 3 , 4 , 9 ,
                            MakeTypeList< Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  ,
                                          Line<3>  >::type ,
                            QuadrilateralEdgeNodeMap ,
                            MakeTypeList< Quadrilateral<9>  ,
                                          Quadrilateral<9>  >::type ,
                            ShellQuadrilateralFaceNodeMap >
{ typedef ShellQuadrilateral<4> base ; };

//------------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 , 4 > ,
                IndexList< 1 , 2 , 5 > ,
                IndexList< 2 , 0 , 6 > ,
                IndexList< 0 , 3 , 7 > ,
                IndexList< 1 , 3 , 8 > ,
                IndexList< 2 , 3 , 9 > >::type TetrahedronEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0 , 1 , 3 ,   4 , 8 , 7 > ,
                IndexList< 1 , 2 , 3 ,   5 , 9 , 8 > ,
                IndexList< 0 , 3 , 2 ,   7 , 9 , 6 > ,
                IndexList< 0 , 2 , 1 ,   6 , 5 , 4 >
  >::type TetrahedronSideNodeMap ;

template<> struct Tetrahedron<4> : public
  CellTopologyTraits< 3 , 4 , 4 ,
                      MakeTypeList< Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  >::type ,
                      TetrahedronEdgeNodeMap ,
                      MakeTypeList< Triangle<3>  ,
                                    Triangle<3>  ,
                                    Triangle<3>  ,
                                    Triangle<3>  >::type ,
                      TetrahedronSideNodeMap >
{ typedef Tetrahedron<4> base ; };

template<> struct Tetrahedron<10> : public
  CellTopologyTraits< 3 , 4 , 10 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      TetrahedronEdgeNodeMap ,
                      MakeTypeList< Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  >::type ,
                      TetrahedronSideNodeMap >
{ typedef Tetrahedron<4> base ; };

template<> struct Tetrahedron<8> : public
  CellTopologyTraits< 3 , 4 , 8 ,
                      MakeTypeList< Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  >::type ,
                      TetrahedronEdgeNodeMap ,
                      MakeTypeList< Triangle<4>  ,
                                    Triangle<4>  ,
                                    Triangle<4>  ,
                                    Triangle<4>  >::type ,
                      TetrahedronSideNodeMap >
{ typedef Tetrahedron<4> base ; };

//------------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 ,   5 > ,
                IndexList< 1 , 2 ,   6 > ,
                IndexList< 2 , 3 ,   7 > ,
                IndexList< 3 , 0 ,   8 > ,
                IndexList< 0 , 4 ,   9 > ,
                IndexList< 1 , 4 ,  10 > ,
                IndexList< 2 , 4 ,  11 > ,
                IndexList< 3 , 4 ,  12 > >::type PyramidEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0 , 1 , 4 ,       5 , 10 ,  9 > ,
                IndexList< 1 , 2 , 4 ,       6 , 11 , 10 > ,
                IndexList< 2 , 3 , 4 ,       7 , 12 , 11 > ,
                IndexList< 3 , 0 , 4 ,       8 ,  9 , 12 > ,
                IndexList< 0 , 3 , 2 , 1 ,   8 ,  7 ,  6 ,  5 ,  13 >
  >::type PyramidFaceNodeMap ;

template<> struct Pyramid<5> : public
  CellTopologyTraits< 3 , 5 , 5 ,
                      MakeTypeList< Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  >::type ,
                      PyramidEdgeNodeMap ,
                      MakeTypeList< Triangle<3>  ,
                                    Triangle<3>  ,
                                    Triangle<3>  ,
                                    Triangle<3>  ,
                                    Quadrilateral<4>  >::type ,
                      PyramidFaceNodeMap >
{ typedef Pyramid<5> base ; };

template<> struct Pyramid<13> : public
  CellTopologyTraits< 3 , 5 , 13 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      PyramidEdgeNodeMap ,
                      MakeTypeList< Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Quadrilateral<8>  >::type ,
                      PyramidFaceNodeMap >
{ typedef Pyramid<5> base ; };

template<> struct Pyramid<14> : public
  CellTopologyTraits< 3 , 5 , 14 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      PyramidEdgeNodeMap ,
                      MakeTypeList< Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  ,
                                    Quadrilateral<9>  >::type ,
                      PyramidFaceNodeMap >
{ typedef Pyramid<5> base ; };

//------------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 ,   6 > ,
                IndexList< 1 , 2 ,   7 > ,
                IndexList< 2 , 0 ,   8 > ,
                IndexList< 3 , 4 ,  12 > ,
                IndexList< 4 , 5 ,  13 > ,
                IndexList< 5 , 3 ,  14 > ,
                IndexList< 0 , 3 ,   9 > ,
                IndexList< 1 , 4 ,  10 > ,
                IndexList< 2 , 5 ,  11 >
  >::type WedgeEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0 , 1 , 4 , 3 ,   6 , 10 , 12 ,  9 ,  15 > ,
                IndexList< 1 , 2 , 5 , 4 ,   7 , 11 , 13 , 10 ,  16 > ,
                IndexList< 0 , 3 , 5 , 2 ,   9 , 14 , 11 ,  8 ,  17 > ,
                IndexList< 0 , 2 , 1 ,       8 ,  7 ,  6 > ,
                IndexList< 3 , 4 , 5 ,      12 , 13 , 14 >
  >::type WedgeFaceNodeMap ;

template<> struct Wedge<6> : public
  CellTopologyTraits< 3 , 6 , 6 ,
                      MakeTypeList< Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  >::type ,
                      WedgeEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<4>  ,
                                    Quadrilateral<4>  ,
                                    Quadrilateral<4>  ,
                                    Triangle<3>  ,
                                    Triangle<3>  >::type ,
                      WedgeFaceNodeMap >
{ typedef Wedge<6> base ; };

template<> struct Wedge<15> : public
  CellTopologyTraits< 3 , 6 , 15 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      WedgeEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  >::type ,
                      WedgeFaceNodeMap >
{ typedef Wedge<6> base ; };

template<> struct Wedge<18> : public
  CellTopologyTraits< 3 , 6 , 18 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      WedgeEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Triangle<6>  ,
                                    Triangle<6>  >::type ,
                      WedgeFaceNodeMap >
{ typedef Wedge<6> base ; };

//------------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 ,   8 > ,
                IndexList< 1 , 2 ,   9 > ,
                IndexList< 2 , 3 ,  10 > ,
                IndexList< 3 , 0 ,  11 > ,
                IndexList< 4 , 5 ,  16 > ,
                IndexList< 5 , 6 ,  17 > ,
                IndexList< 6 , 7 ,  18 > ,
                IndexList< 7 , 4 ,  19 > ,
                IndexList< 0 , 4 ,  12 > ,
                IndexList< 1 , 5 ,  13 > ,
                IndexList< 2 , 6 ,  14 > ,
                IndexList< 3 , 7 ,  15 > >::type
  HexahedronEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0, 1, 5, 4,   8, 13, 16, 12,   25 > ,
                IndexList< 1, 2, 6, 5,   9, 14, 17, 13,   24 > ,
                IndexList< 2, 3, 7, 6,  10, 15, 18, 14,   26 > ,
                IndexList< 0, 4, 7, 3,  12, 19, 15, 11,   23 > ,
                IndexList< 0, 3, 2, 1,  11, 10,  9,  8,   21 > ,
                IndexList< 4, 5, 6, 7,  16, 17, 18, 19,   22 > >::type
  HexahedronFaceNodeMap ;

//----------------------------------------------------------------------

template<> struct Hexahedron<8> : public
  CellTopologyTraits< 3 , 8 , 8 ,
                      MakeTypeList< Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  ,
                                    Line<2>  >::type ,
                      HexahedronEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<4>  ,
                                    Quadrilateral<4>  ,
                                    Quadrilateral<4>  ,
                                    Quadrilateral<4>  ,
                                    Quadrilateral<4>  ,
                                    Quadrilateral<4>  >::type ,
                      HexahedronFaceNodeMap >
{
  typedef Hexahedron<8> base ;
};

template<> struct Hexahedron<20> : public
  CellTopologyTraits< 3 , 8 , 20 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      HexahedronEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  ,
                                    Quadrilateral<8>  >::type ,
                      HexahedronFaceNodeMap >
{
  typedef Hexahedron<8> base ;
};

template<> struct Hexahedron<27> : public
  CellTopologyTraits< 3 , 8 , 27 ,
                      MakeTypeList< Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  ,
                                    Line<3>  >::type ,
                      HexahedronEdgeNodeMap ,
                      MakeTypeList< Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  ,
                                    Quadrilateral<9>  >::type ,
                      HexahedronFaceNodeMap >
{
  typedef Hexahedron<8> base ;
};

//------------------------------------------------------------------------

template<> inline
const CellTopologyData * getCellTopologyData< Node::Traits >()
{ return getCellTopologyData< Node >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Particle::Traits >()
{ return getCellTopologyData< Particle >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Line<2>::Traits >()
{ return getCellTopologyData< Line<2> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Line<3>::Traits >()
{ return getCellTopologyData< Line<3> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Beam<2>::Traits >()
{ return getCellTopologyData< Beam<2> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Beam<3>::Traits >()
{ return getCellTopologyData< Beam<3> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< ShellLine<2>::Traits >()
{ return getCellTopologyData< ShellLine<2> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< ShellLine<3>::Traits >()
{ return getCellTopologyData< ShellLine<3> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Triangle<3>::Traits >()
{ return getCellTopologyData< Triangle<3> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Triangle<4>::Traits >()
{ return getCellTopologyData< Triangle<4> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Triangle<6>::Traits >()
{ return getCellTopologyData< Triangle<6> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< ShellTriangle<3>::Traits >()
{ return getCellTopologyData< ShellTriangle<3> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< ShellTriangle<6>::Traits >()
{ return getCellTopologyData< ShellTriangle<6> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Quadrilateral<4>::Traits >()
{ return getCellTopologyData< Quadrilateral<4> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Quadrilateral<8>::Traits >()
{ return getCellTopologyData< Quadrilateral<8> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Quadrilateral<9>::Traits >()
{ return getCellTopologyData< Quadrilateral<9> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< ShellQuadrilateral<4>::Traits >()
{ return getCellTopologyData< ShellQuadrilateral<4> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< ShellQuadrilateral<8> ::Traits>()
{ return getCellTopologyData< ShellQuadrilateral<8> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< ShellQuadrilateral<9>::Traits >()
{ return getCellTopologyData< ShellQuadrilateral<9> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Tetrahedron<4>::Traits >()
{ return getCellTopologyData< Tetrahedron<4> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Tetrahedron<8>::Traits >()
{ return getCellTopologyData< Tetrahedron<8> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Tetrahedron<10>::Traits >()
{ return getCellTopologyData< Tetrahedron<10> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Pyramid<5>::Traits >()
{ return getCellTopologyData< Pyramid<5> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Pyramid<13>::Traits >()
{ return getCellTopologyData< Pyramid<13> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Pyramid<14>::Traits >()
{ return getCellTopologyData< Pyramid<14> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Wedge<6>::Traits >()
{ return getCellTopologyData< Wedge<6> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Wedge<15>::Traits >()
{ return getCellTopologyData< Wedge<15> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Wedge<18>::Traits >()
{ return getCellTopologyData< Wedge<18> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Hexahedron<8>::Traits >()
{ return getCellTopologyData< Hexahedron<8> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Hexahedron<20>::Traits >()
{ return getCellTopologyData< Hexahedron<20> >(); }

template<> inline
const CellTopologyData * getCellTopologyData< Hexahedron<27>::Traits >()
{ return getCellTopologyData< Hexahedron<27> >(); }

//------------------------------------------------------------------------

const unsigned * index_identity_array();

const CellTopologyData::Subcell * subcell_nodes_array();

#endif /* DOXYGEN_COMPILE */

/** \} */
} // namespace shards

#endif // Shards_BasicTopologies_hpp

