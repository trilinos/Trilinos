// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file
\brief Definition of Intrepid's canonical cells and dummy cell definitions for user-defined cells.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CELLTEMPLATES_HPP
#define INTREPID_CELLTEMPLATES_HPP

namespace Intrepid {
  
  template<class Scalar>
  const ConnMapTemplate MultiCell<Scalar>::connMapCanonical_[CELL_CANONICAL_MAX][3] ={
  //------------------------------------------------------------------------------------------------       
  { // 0-simplex, i.e. node
    { // node->1cell - the first connMapTemplate struct
      0,             - 
      0,
      {0},
      {},
      {{0}}
     },
     { // node->2cell - the second connMapTemplate struct
       0,
       0,
      {0},
      {},
      {{0}}
     },
     { // node->3cell - the third connMapTemplate struct
       0,
       0,
      {0},
      {},
      {{0}}
     }
  },// end node
  //------------------------------------------------------------------------------------------------        
  { // 1-simplex, i.e. edge
    { // edge->1cell
      1,                                    // topological dimension of the edge is 1
      1,                                    // number of subcells that are 1-cells is 1
      {2},                                  // the only 1-subcell has 2 nodes
      {CELL_EDGE},                          // the type of this subcell is CELL_CELL_EDGE
      {{0,1}}                               // the endpoints are local ndes 0 and 1
    },
    { // edge->2cell
      1,                                    // topological dimension of the edge is 1            
      0,                                    // an edge has no 2-subcells, i.e., 2D subcells
      {0},
      {},
      {{0}}
    },
    { // edge->3cell
      1,                                    // topological dimension of the edge is 1 
      0,                                    // an edge has no 3-subcells, i.e., 3D subcells
      {0},
      {},
      {{0}}
    }
  },  // end edge
  //------------------------------------------------------------------------------------------------        
  { // 2-simplex, i.e. triangle
    { // tri->1cell
      2,                                    // topological dimension of the triangle is 2 
      3,                                    // a triangle has 3 1-subcells (edges), numbered 0,1,2
      {2,2,2},                              // each 1-subcell has 2 vertices
      {CELL_EDGE,CELL_EDGE,CELL_EDGE},      // each 1-subcell is an edge
      {{0,1}, {1,2}, {2,0}}                 // the vertices of edges 0,1,2
    },
    { // tri->2cell
      2,                                    // topological dimension of the triangle is 2 
      1,                                    // a triangle has 1 2-subcells - itself!
      {3},                                  // the 2-subcell has 3 vertices
      {CELL_TRI},                           // the only 2-subcell is a triangle!
      {{0,1,2}}                             // this is the local order of the vertices in the triangle
    },
    { // tri->3cell
      2,                                    // topological dimension of the triangle is 2 
      0,                                    // a triangle has no 3-subcells, i.e., 3D subcells!
      {0},
      {},
      {{0}}
      }
    }, // end tri
    //------------------------------------------------------------------------------------------------        
    { // QUAD
      { // quad->1cell
        2,                                    // topological dimension of the quad is 2
        4,                                    // a quad has 4 1-subcells (edges), numbered 0,1,2,3
        {2,2,2,2},                            // each 1-subcell has 2 vertices
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,       // each 1-subcell is an edge
         CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,0}}          // these are the vertices of edges 0,1,2,3
      },
      { // quad->2cell
        2,                                    // topological dimension of the quad is 2
        1,                                    // a quad has a single 2-cell - the quad itself
        {4},                                  // it has 4 vertices
        {CELL_QUAD},                          // its type is quad
        {{0,1,2,3}}                           // this is the local order of the vertices in the quad
      },
      { // quad->3cell
        2,                                    // topological dimension of the quad is 2
        0,                                    // a quad has no 3-subcells!
        {0},
        {},
        {{0}}
      }
    },  // end QUAD 
    //------------------------------------------------------------------------------------------------            
    {   // 3-simplex, i.e. tetrahedron
      { // tet->1cell
        3,                                    // topological dimension of the tet is 3
        6,                                    // a tet has 6 1-subcells (edges) numbered 0,1,2,3,4,5
        {2,2,2,2,2,2},                        // each 1-subcell has 2 vertices
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,       // each 1-subcell is an edge
         CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,0},                 // these are the vertices of edges 0,1,2,3,4,5
         {0,3}, {1,3}, {2,3}}
      },
      { // tet->2cell
        3,                                    // topological dimension of the tet is 3
        4,                                    // a tet has 4 2-subcells (faces) numbered 0,1,2,3
        {3,3,3,3},                            // all 2-subcells have 3 vertices
        {CELL_TRI,CELL_TRI,CELL_TRI,CELL_TRI},// all 2-subcells are triangles
        {{0,1,3}, {1,2,3}, {0,3,2}, {0,2,1}}  // these are the local vertices of faces 0,1,2,3
      },
      { // tet->3cell
        3,                                    // topological dimension of the tet is 3
        1,                                    // a tet has a single 3-cell - the tet itself!
        {4},
        {CELL_TET},
        {{0,1,2,3}}                           // this is the local vertex order on the tet
      }
    },  // end TET
    //------------------------------------------------------------------------------------------------            
    { // HEXAHEDRON
      { // hex->1cell
        3,                                    // topological dimension of the hex is 3
        12,                                   // a hex has 12 1-subcells (edges) numbered 0,1,...,11
        {2,2,2,2,2,2,2,2,2,2,2,2},            // etc...
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,5}, {2,6}, {3,7}, {4,5}, {5,6}, {6,7}, {7,4}}
      },
      { // hex->2cell
        3,
        6,                                    // a hex has 6 2-subcells (faces) numbered 0,1,2,3,4,5   
        {4,4,4,4,4,4},                        // all 2-subcells have 4 vertices
        {CELL_QUAD,CELL_QUAD,CELL_QUAD,       // all 2-subcells are quads
         CELL_QUAD,CELL_QUAD,CELL_QUAD},
        {{0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7},  // side faces: 0,1,2,3
         {0,3,2,1}, {4,5,6,7}}                        // bottom face: 4; top face: 5
      },                                      // all faces are oriented by outward pointing normal
      { // hex->3cell
        3,
        1,                                    // a hex has 1 3-subcell - the hex itself
        {8},
        {CELL_HEX},
        {{0,1,2,3,4,5,6,7}}                   // this is the local vertex order on the hex
      }
    },  // end HEX
    //------------------------------------------------------------------------------------------------            
    {   // PYRAMID
      { // pyramid->1cell
        3,
        8,
        {2,2,2,2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,4}, {2,4}, {3,4}}
      },
      { // pyramid->2cell
        3,
        5,
        {4,3,3,3,3,3},
        {CELL_QUAD,CELL_TRI,CELL_TRI,CELL_TRI,CELL_TRI},
        {{0,3,2,1}, {0,1,4}, {1,2,4}, {2,3,4}, {3,0,4}}
      },
      { // pyramid->3cell
        3,
        1,
        {5},
        {CELL_PYRAMID},
        {{0,1,2,3,4}}
      }
    },  // end PYRAMID
    //------------------------------------------------------------------------------------------------ 
    {   // PENTAGON
      { // pentagon->1cell
        2,
        5,
        {2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,4}, {4,0}}
      },
      { // pentagon->2cell
        2,
        1,
        {5},
        {CELL_PENTAGON},
        {{0,1,2,3,4}}
      },
      { // pentagon->3cell
        2,
        0,
        {0},
        {},
        {{0}}
      }
    },  // end PENTAGON
    //------------------------------------------------------------------------------------------------ 
    {   // HEXAGON
      { // hexagon->1cell
        2,
        6,
        {2,2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,4}, {4,5},{5,0}}
      },
      { // hexagon->2cell
        2,
        1,
        {6},
        {CELL_HEXAGON},
        {{0,1,2,3,4,5}}
      },
      { // hexagon->3cell
        2,
        0,
        {0},
        {},
        {{0}}
      }
    },  // end HEXAGON
    //------------------------------------------------------------------------------------------------ 
    {   // HEPTAGON
      { // heptagon->1cell
        2,
        7,
        {2,2,2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,4}, {4,5},{5,6},{6,0}}
      },
      { // heptagon->2cell
        2,
        1,
        {7},
        {CELL_HEPTAGON},
        {{0,1,2,3,4,5,6}}
      },
      { // heptagon->3cell
        2,
        0,
        {0},
        {},
        {{0}}
      }
    },  // end HEPTAGON
    //------------------------------------------------------------------------------------------------ 
    {   // OCTAGON
      { // octagon->1cell
        2,
        8,
        {2,2,2,2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,4}, {4,5},{5,6},{6,7},{7,0}}
      },
      { // octagon->2cell
        2,
        1,
        {8},
        {CELL_OCTAGON},
        {{0,1,2,3,4,5,6,7}}
      },
      { // octagon->3cell
        2,
        0,
        {0},
        {},
        {{0}}
      }
    },  // end OCTAGON
    //------------------------------------------------------------------------------------------------ 
    {   // NONAGON
      { // nonagon->1cell
        2,
        9,
        {2,2,2,2,2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,4}, {4,5},{5,6},{6,7},{7,8},{8,0}}
      },
      { // nonagon->2cell
        2,
        1,
        {9},
        {CELL_NONAGON},
        {{0,1,2,3,4,5,6,7,8}}
      },
      { // nonagon->3cell
        2,
        0,
        {0},
        {},
        {{0}}
      }
    },  // end NONAGON
    //------------------------------------------------------------------------------------------------ 
    {   // DECAGON
      { // decagon->1cell
        2,
        10,
        {2,2,2,2,2,2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,4}, {4,5},{5,6},{6,7},{7,8},{8,9},{9,0}}
      },
      { // decagon->2cell
        2,
        1,
        {10},
        {CELL_DECAGON},
        {{0,1,2,3,4,5,6,7,8}}
      },
      { // decagon->3cell
        2,
        0,
        {0},
        {},
        {{0}}
      }
    },  // end DECAGON
    //---------------------------------------------------------------------------------------------// 
    //                                        Prismatic cells                                      //
    //---------------------------------------------------------------------------------------------//   
    {   // TRIPRISM = prism with triangular base (wedge)
      { // triprism->1cell
        3,
        9,
        {2,2,2,2,2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,0}, {0,3}, {1,4}, {2,5}, {3,4}, {4,5}, {5,3}}
      },
      { // triprism->2cell
        3,
        5,
        {4,4,4,3,3},
        {CELL_QUAD,CELL_QUAD,CELL_QUAD,CELL_TRI,CELL_TRI},
        {{0,1,4,3}, {1,2,5,4}, {2,0,3,5}, {0,2,1}, {3,4,5}}
      },
      { // triprism->3cell
        3,
        1,
        {6},
        {CELL_TRIPRISM},
        {{0,1,2,3,4,5}}
      }
    },  // end TRIPRISM
    //------------------------------------------------------------------------------------------------            
    {   // PENTAPRISM = prism with pentagon base
      { // pentaprism->1cell
        3,                                    // topological dimension of pentaprism is 3
        15,                                   // it has 15 1-subcells (edges)
        {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,4}, {4,0},   // edges 0,1,2,3,4 are on the bottom face
         {0,5}, {1,6}, {2,7}, {3,8}, {4,9},   // edges 5,6,7,8,9 are on the vertical faces
         {5,6}, {6,7}, {7,8}, {8,9}, {9,5}}   // edges 10,11,12,13,14 are on the top face
      },
      { // pentaprism->2cell
        3,
        7,                                     // pentaprism has 7 2-cells (faces), numbered 0,..,6
        {4,4,4,4,4,5,5},
        {CELL_QUAD,CELL_QUAD,CELL_QUAD,CELL_QUAD,CELL_QUAD,     // side faces are quads
         CELL_PENTAGON,CELL_PENTAGON},                          // bottom and top faces are pentagons
        {{0,1,6,5}, {1,2,7,6}, {2,3,8,7}, {3,4,9,8}, {4,0,5,9}, // faces 0,1,2,3,4 are side faces
         {0,4,3,2,1}, {5,6,7,8,9}}                              // face 5 is bottom, face 6 is top
      },                                       // all faces are oriented by outward pointing normal
      { // pentaprism->3cell
        3,
        1,
        {10},
        {CELL_PENTAPRISM},
        {{0,1,2,3,4,5,6,7,8,9}}
      }
    },  // end PENTAPRISM
    //------------------------------------------------------------------------------------------------            
    {   // HEXAPRISM = prism with hexagon base
      { // hexaaprism->1cell
        3,                                    // topological dimension of hexaprism is 3
        18,                                   // it has 18 1-subcells (edges)
        {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
        {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
         CELL_EDGE,CELL_EDGE,CELL_EDGE},
        {{0,1}, {1,2}, {2,3}, {3,4},  {4,5},   {5,0},     // edges 0,1,2,3,4,5 are on the bottom face
         {0,6}, {1,7}, {2,8}, {3,9},  {4,10},  {5,11},    // edges 6,7,8,9,10,11 are on the side faces
         {6,7}, {7,8}, {8,9}, {9,10}, {10,11}, {11,6}}    // edges 12,13,14,15,16,17 are on the top face
      },
      { // hexaprism->2cell
        3,
        8,                                    // pentaprism has 8 2-cells (faces), numbered 0,..,7
        {4,4,4,4,4,4,6,6},
        {CELL_QUAD,CELL_QUAD,CELL_QUAD,       // the side faces  are quads
         CELL_QUAD,CELL_QUAD,CELL_QUAD,   
         CELL_PENTAGON,CELL_PENTAGON},        // the bottom and top faces are hexagons
        {{0,1,7,6},  {1,2,8,7},   {2,3,9,8},  // local vertex order on side faces 0,1,2,3,4,5
         {3,4,10,9}, {4,5,11,10}, {5,0,6,11},
         {0,5,4,3,2,1}, {6,7,8,9,10,11}}      // local vertex order on bottom (6) and top (7) faces
      },                                      // all faces are oriented by outward pointing normal
      { // hexaprism->3cell
        3,
        1,
        {10},
        {CELL_HEXAPRISM},
        {{0,1,2,3,4,5,6,7,8,9,10,11}}
      }
    }  // end PENTAPRISM
  }; // end connMapCanonical
  
  //---------------------------------------------------------------------------------------------//            
  //                              For user defined cell templates                                //
  //---------------------------------------------------------------------------------------------//            
  
  template<class Scalar>
  ConnMapTemplate MultiCell<Scalar>::connMapCustom_[CELL_MAX-CELL_CANONICAL_MAX-1][3] = {
    //------------------------------------------------------------------------------------------------            
    {   // POLY0  - empty template for user defined cell
      { // POLY0 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY0 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY0 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    },  // end POLY0
    //------------------------------------------------------------------------------------------------            
    {   // POLY1  - empty template for user defined cell
      { // POLY1 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY1 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY1 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    },  // end POLY1
    //------------------------------------------------------------------------------------------------            
    {   // POLY2  - empty template for user defined cell
      { // POLY2 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY2 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY2 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    },  // end POLY0
    //------------------------------------------------------------------------------------------------            
    {   // POLY3  - empty template for user defined cell
      { // POLY0 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY3 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY3 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    },  // end POLY3
    //------------------------------------------------------------------------------------------------            
    {   // POLY4  - empty template for user defined cell
      { // POLY4 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY4 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY4 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    },  // end POLY4
    //------------------------------------------------------------------------------------------------            
    {   // POLY5  - empty template for user defined cell
      { // POLY5 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY5 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY5 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    },  // end POLY5
    //------------------------------------------------------------------------------------------------            
    {   // POLY6  - empty template for user defined cell
      { // POLY6 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY6 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY6 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    },  // end POLY6
    //------------------------------------------------------------------------------------------------            
    {   // POLY7  - empty template for user defined cell
      { // POLY7 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY7 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY7 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    },  // end POLY7
    //------------------------------------------------------------------------------------------------            
    {   // POLY8  - empty template for user defined cell
      { // POLY8 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY8 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY8 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    },  // end POLY8
    //------------------------------------------------------------------------------------------------            
    {   // POLY9  - empty template for user defined cell
      { // POLY9 -> 1-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY9 -> 2-cell
        0, 0, {0}, {}, {{0}}
      },
      { // POLY9 -> 2cell
        0, 0, {0}, {}, {{0}}
      }
    }  // end POLY9
    //------------------------------------------------------------------------------------------------            
  }; // end connMapCustom
   
   
} // namespace Intrepid
#endif

//==============================================================================
//    D O X Y G E N        D O C U M E N T A T I O N:   P A G E S             //   
//==============================================================================
/*! \page cell_templates Cell templates in Intrepid

\section cell_tempalates_intro Introduction

Admissible cell types and their templates in Intrepid are denoted by an enumerated type ECell. 
There are two kinds of cells and cell tempates 
- canonical cell templates are defined by Intrepid
- custom cell templates are defined by the user. Their names must be one of CELL_POLY0 - CELL_POLY9

\arg CELL_NODE = 0,       // 0-simplex, i.e. node
\arg CELL_EDGE,           // 1-simplex, i.e. edge
\arg CELL_TRI,            // 2-simplex, i.e. triangular cell
\arg CELL_QUAD,           // quadrilateral cell
\arg CELL_TET,            // 3-simplex, i.e. tetrahedral cell
\arg CELL_HEX,            // hexahedral cell
\arg CELL_PYRAMID,        // pyramid cell
\arg CELL_PENTAGON,       // polygon with 5 sides
\arg CELL_HEXAGON,        // polygon with 6 sides
\arg CELL_HEPTAGON,       // polygon with 7 sides
\arg CELL_OCTAGON,        // polygon with 8 sides
\arg CELL_NONAGON,        // polygon with 9 sides
\arg CELL_DECAGON,        // polygon with 10 sides
\arg CELL_TRIPRISM,       // prismatic polyhedron with a triangle base
\arg CELL_PENTAPRISM,     // prismatic polyhedron with a pentagon base
\arg CELL_HEXAPRISM,      // prismatic polyhedron with a hexagon base
\arg CELL_CANONICAL_MAX,  // used as the maximum number of canonical types (current value = 16)
\arg CELL_POLY0,          // user defined cell 
\arg CELL_POLY1,          // user defined cell 
\arg CELL_POLY2,          // user defined cell 
\arg CELL_POLY3,          // user defined cell 
\arg CELL_POLY4,          // user defined cell 
\arg CELL_POLY5,          // user defined cell 
\arg CELL_POLY6,          // user defined cell 
\arg CELL_POLY7,          // user defined cell 
\arg CELL_POLY8,          // user defined cell 
\arg CELL_POLY9,          // user defined cell 
\arg CELL_MAX             // placeholder for looping over all types        (current value = 27)


\section cell_templates_reference Reference cells

  In many cases, a canonical cell can be obtained as an image of a standard "reference" cell
  under a one-to-one polynomial mapping. The reference cell can be thought of as a standard
  instance of the given cell shape. 
  
  If a cell has a standard instance, one can take advantage of this by defining many pieces of data
  such as integration points, only once on the reference cell. Intrepid follows the following
  convention regarding the standard cell shapes, embedded in three dimensional Euclidean space:
  
  CELL_EDGE      -> (-1,0,0),(1,0,0)
  
  CELL_TRI       -> (0,0,0),(1,0,0),(0,1,0)
  
  CELL_QUAD      -> (-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)
  
  TET            -> (0,0,0),(1,0,0), (0,1,0), (0,0,1)
  
  HEX            -> (-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1)  - bottom face 
                    (-1,-1, 1),(1,-1, 1),(1,1, 1),(-1,1, 1)  - top face
  
  TRIPRISM       -> (0,0,0),(1,0,0),(0,1,0)                  - bottom CELL_TRI
                    (0,0,1),(1,0,1),(0,1,1)                  - upper CELL_TRI
  
  PYRAMID        -> (-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)      - CELL_QUAD base and (0,0,1) - top
 

 Reference cells are used in reconstruction methods based on pullbacks, such as FEM. In such methods, 
 basis set on an ambient cell is obtained by pullback of a basis set defined on a standard reference 
 cell. This approach is practical for cells that have polynomial pullbacks. Intrepid supports pullback
 (FEM) reconstructions for CELL_EDGE, CELL_TRI, CELL_QUAD, TET, HEX, PYRAMID and PRISM cells. 
 
 \subsection cell_tempaltes_ref_HEX Reference HEX cell
 
 The vertices, edges and faces of the reference cell3 are numbered as follows using the negative 
 y-axis as standard line of sight:
 @verbatim
          (1,1,1)
     *-------*        *-------*      ^ z
    /|       |       /       /|      |
   / |       |      /       / |      |  ^ y
  /  |       |     /       /  |      | /
 *   |       |    *-------*   |      |/
 |   *-------*    |       |   *  ----/------------>  
 |  /       /     |       |  /      /|            x
 | /       /      |       | /      / |  
 |/       /       |       |/      ^
 *-------*        *-------*    line of sight
 (-1,-1,-1)
 
@endverbatim
The drawing on the left shows: left, bottom and back faces. The front, right
and top faces are shown on the other drawing.

\subsubsection hexvert Vertices

The vertices are numbered first on the bottom face and then by shift on the top face:
@verbatim
    7-------6        7-------6->(1,1,1)
   /|       |       /       /|
  / |       |      /       / |
 /  |       |     /       /  |
4   |       |    4-------5   |
|   3-------2    |       |   2->(1,1,-1)
|  /       /     |       |  /
| /       /      |       | /
|/       /       |       |/
0-------1        0-------1 -> (1,-1,-1)
@endverbatim

@subsubsection hexedg Edges

The local edge numbering is to number first the edges on the bottom, then the vertical edges 
connecting the bottom and top faces, and finally, the top edges.
@verbatim
    *--10---*        *--10---*
   /|       |       /       /|
 11 |       6      11      9 6
 /  7       |     /       /  |
*   |       |    *---8---*   |
|   *---2---*    |       |   *
4  /       /     |       5  /
| 3       1      4       | 1
|/       /       |       |/
*---0---*        *---0---*
@endverbatim
<b> Local orientation of an edge </b> is defined by the order of its endpoints:
@verbatim
    *---<---*        *---<---*
   /|       |       /       /|
  v |       ^      v       ^ ^
 /  ^       |     /       /  |
*   |       |    *--->---*   |
|   *---<---*    |       |   *
|  /       /     |       ^  /
^ v       ^      ^       | ^
|/       /       |       |/
*--->---*        *--->---*
@endverbatim

@subsubsection hexfac Faces

The faces are numbered in this order: front, right, back, left, bottom, top
@verbatim
    *-------*        *-------*
   /|       |       /       /|
  / |   2   |      /   5   / |
 /  |       |     /       /  |
*   |       |    *-------*   |
| 3 *-------*    |       | 1 *
|  /       /     |       |  /
| /   4   /      |   0   | /
|/       /       |       |/
*-------*        *-------*
@endverbatim

<b> Local orientation of a face </b> is defined by the order of its vertices and is
non-unique. For example, face 0 can be ordered as (0,4,5,1) which corresponds to
choosing (0,1,0) as a unit normal, or as (0,1,5,4) which corresponds to choosing (0,-1,0)
as a unit normal.  

The standard orientation is to choose the unit outer normal to each face. 

@verbatim
    ^                ^
    *-------* /      *---|---*
   /|       |       /    |  /|
  / |   2   |      /   5   / |
 /          |     /       / ----> 
*   |       |    *-------*   |
| 3 *-------*    |       | 1 *
|  /       /     |  /    |  /
| /   4   /      | \/ 0  | /
|/       /       |       |/
*-----|-*        *-------*
|
|
v

@endverbatim









*/
