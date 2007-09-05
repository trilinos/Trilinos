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
    \brief Contains all canonical and dummy-custom cell definitions,
           including downward subcell connectivities.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/


#ifndef INTREPID_CELLTEMPLATES_HPP
#define INTREPID_CELLTEMPLATES_HPP

namespace Intrepid {
	
    //-------------------------------------------------------------------------------------//
    //                          Definition of cell templates                               //
    //                                                                                     //
    //  WARNING: Order of connectivity templates must follow the order of the enumerated   //
    //           CellType                                                                  //
    //-------------------------------------------------------------------------------------//
	
  
/* 
 In many cases, a canonical cell can be obtained as an image of a standard "reference" cell
 under a one-to-one polynomial mapping. The reference cell can be thought of as a standard
 instance of the given cell shape. 
 
 If a cell has a standard instance, one can take advantage of this by defining many pieces of data
 such as integration points, only once on the reference cell. Intrepid follows the following
 convention regarding the standard cell shapes, embedded in three dimensional Euclidean space:
 
 EDGE     -> (-1,0,0),(1,0,0)
 
 TRI      -> (0,0,0),(1,0,0),(0,1,0)
 
 QUAD     -> (-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)
 
 TET      -> (0,0,0),(1,0,0), (0,1,0), (0,0,1)
 
 HEX      -> (-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1) - bottom face 
             (-1,-1, 1),(1,-1, 1),(1,1, 1),(-1,1, 1) - top face
 
 PRISM    -> (0,0,0),(1,0,0),(0,1,0) - bottom TRI
             (0,0,1),(1,0,1),(0,1,1) - upper TRI
 
 PYRAMID  -> (-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0) - QUAD base and (0,0,1) - top

 */
  
  
template<class Scalar>
const ConnMapTemplate MultiCell<Scalar>::conn_map_canonical [MAXCANONICAL][3] =
{
  {   // 0-simplex, i.e. node
    { // node->1cell
      0,
      0,
      {0},
      {},
      {{0}}
    },
    { // node->2cell
      0,
      0,
      {0},
      {},
      {{0}}
    },
    { // node->3cell
      0,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end node


  {   // 1-simplex, i.e. edge
    { // edge->1cell
      1,
      1,
      {2},
      {EDGE},
      {{0,1}}
    },
    { // edge->2cell
      1,
      0,
      {0},
      {},
      {{0}}
    },
    { // edge->3cell
      1,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end edge


  {   // 2-simplex, i.e. triangle
    { // tri->1cell
      2,
      3,
      {2,2,2},
      {EDGE,EDGE,EDGE},
      {{0,1}, {1,2}, {2,0}}
    },
    { // tri->2cell
      2,
      1,
      {3},
      {TRI},
      {{0,1,2}}
    },
    { // tri->3cell
      2,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end tri


  {   // quadrilateral
    { // quad->1cell
      2,
      4,
      {2,2,2,2},
      {EDGE,EDGE,EDGE,EDGE},
      {{0,1}, {1,2}, {2,3}, {3,0}}
    },
    { // quad->2cell
      2,
      1,
      {4},
      {QUAD},
      {{0,1,2,3}}
    },
    { // quad->3cell
      2,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end quad 


  {   // 3-simplex, i.e. tetrahedron
    { // tet->1cell
      3,
      6,
      {2,2,2,2,2,2},
      {EDGE,EDGE,EDGE,EDGE,EDGE,EDGE},
      {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}}
    },
    { // tet->2cell
      3,
      4,
      {3,3,3,3},
      {TRI,TRI,TRI,TRI},
      {{0,1,3}, {1,2,3}, {0,3,2}, {0,2,1}}
    },
    { // tet->3cell
      3,
      1,
      {4},
      {TET},
      {{0,1,2,3}}
    }
  },  // end tet


  {   // hexahedron
    { // hex->1cell
      3,
      12,
      {2,2,2,2,2,2,2,2,2,2,2,2},
      {EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE},
      {{0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,5}, {2,6}, {3,7}, {4,5}, {5,6}, {6,7}, {7,4}}
    },
    { // hex->2cell
      3,
      6,
      {4,4,4,4,4,4},
      {QUAD,QUAD,QUAD,QUAD,QUAD,QUAD},
      {{0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7}, {0,3,2,1}, {4,5,6,7}}
    },
    { // hex->3cell
      3,
      1,
      {8},
      {HEX},
      {{0,1,2,3,4,5,6,7}}
    }
  },  // end HEX


  {   // pyramid
    { // pyramid->1cell
      3,
      8,
      {2,2,2,2,2,2,2,2},
      {EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE},
      {{0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,4}, {2,4}, {3,4}}
    },
    { // pyramid->2cell
      3,
      5,
      {4,3,3,3,3,3},
      {TRI,TRI,TRI,TRI,QUAD},
      {{0,1,2,3}, {0,1,4}, {1,2,4}, {2,3,4}, {3,0,4}}
    },
    { // pyramid->3cell
      3,
      1,
      {5},
      {PYRAMID},
      {{0,1,2,3,4}}
    }
  },  // end pyramid


  {   // prism with triangular base (wedge)
    { // prism->1cell
      3,
      9,
      {2,2,2,2,2,2,2,2,2},
      {EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE,EDGE},
      {{0,1}, {1,2}, {2,0}, {0,3}, {1,4}, {2,5}, {3,4}, {4,5}, {5,3}}
    },
    { // prism->2cell
      3,
      5,
      {4,4,4,3,3},
      {QUAD,QUAD,QUAD,TRI,TRI},
      {{0,1,4,3}, {1,2,5,4}, {0,3,5,2}, {0,2,1}, {3,4,5}}
    },
    { // prism->3cell
      3,
      1,
      {6},
      {PRISM},
      {{0,1,2,3,4,5}}
    }
  },  // end prism


  {   // polygon
    { // polygon->1cell
      2,
      0,
      {0},
      {},
      {{0}}
    },
    { // polygon->2cell
      2,
      1,
      {0},
      {POLYGON},
      {{0}}
    },
    { // polygon->3cell
      2,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end polygon 
	
	
  {   // polyhedron
    { // polyhedron->1cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->2cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->3cell
      3,
      1,
      {0},
      {POLYHEDRON},
      {{0}}
    }
  },  // end polyhedron
	
  {   // cellset
    { // cellset->1cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // cellset->2cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // cellset->3cell
      3,
      1,
      {0},
      {CELLSET},
      {{0}}
    }
  }  // end cellset	
};


template<class Scalar>
ConnMapTemplate MultiCell<Scalar>::conn_map_custom[MAXTYPE-MAXCANONICAL-1][3] =
{
  {   // polygon 1
    { // polygon->1cell
      2,
      0,
      {0},
      {},
      {{0}}
    },
    { // polygon->2cell
      2,
      1,
      {0},
      {POLYGON1},
      {{0}}
    },
    { // polygon->3cell
      2,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end polygon 1


  {   // polygon 2
    { // polygon->1cell
      2,
      0,
      {0},
      {},
      {{0}}
    },
    { // polygon->2cell
      2,
      1,
      {0},
      {POLYGON2},
      {{0}}
    },
    { // polygon->3cell
      2,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end polygon 2


  {   // polygon 3
    { // polygon->1cell
      2,
      0,
      {0},
      {},
      {{0}}
    },
    { // polygon->2cell
      2,
      1,
      {0},
      {POLYGON3},
      {{0}}
    },
    { // polygon->3cell
      2,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end polygon 3


  {   // polygon 4
    { // polygon->1cell
      2,
      0,
      {0},
      {},
      {{0}}
    },
    { // polygon->2cell
      2,
      1,
      {0},
      {POLYGON4},
      {{0}}
    },
    { // polygon->3cell
      2,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end polygon 4


  {   // polygon 5
    { // polygon->1cell
      2,
      0,
      {0},
      {},
      {{0}}
    },
    { // polygon->2cell
      2,
      1,
      {0},
      {POLYGON5},
      {{0}}
    },
    { // polygon->3cell
      2,
      0,
      {0},
      {},
      {{0}}
    }
  },  // end polygon 5


  {   // polyhedron 1
    { // polyhedron->1cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->2cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->3cell
      3,
      1,
      {0},
      {POLYHEDRON1},
      {{0}}
    }
  },  // end polyhedron 1


  {   // polyhedron 2
    { // polyhedron->1cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->2cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->3cell
      3,
      1,
      {0},
      {POLYHEDRON2},
      {{0}}
    }
  },  // end polyhedron 2


  {   // polyhedron 3
    { // polyhedron->1cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->2cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->3cell
      3,
      1,
      {0},
      {POLYHEDRON3},
      {{0}}
    }
  },  // end polyhedron 3


  {   // polyhedron 4
    { // polyhedron->1cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->2cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->3cell
      3,
      1,
      {0},
      {POLYHEDRON4},
      {{0}}
    }
  },  // end polyhedron 4


  {   // polyhedron 5
    { // polyhedron->1cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->2cell
      3,
      0,
      {0},
      {},
      {{0}}
    },
    { // polyhedron->3cell
      3,
      1,
      {0},
      {POLYHEDRON5},
      {{0}}
    }
  }   // end polyhedron 5


};

} // namespace Intrepid
#endif

//==============================================================================
//    D O X Y G E N        D O C U M E N T A T I O N:   P A G E S             //   
//==============================================================================
/*!
\page cell_templates Cell templates in Intrepid
 
 \section cell_tempalates_intro Introduction
 
 This section describes the cell types supported by Intrepid. The admissible cell types are
 
 \arg NODE = 0        0-simplex, i.e. node
 \arg EDGE            1-simplex, i.e. edge
 \arg TRI             2-simplex, i.e. triangular cell
 \arg QUAD            quadrilateral cell
 \arg TET             3-simplex, i.e. tetrahedral cell
 \arg HEX             hexahedral cell
 \arg PYRAMID         pyramid cell
 \arg PRISM           prism (wedge) cell
 \arg POLYGON         general (unknown) polygon
 \arg POLYHEDRON      general (unknown) polyhedron
 \arg CELLSET         set of arbitrary cells
 \arg MAXCANONICAL    used as the maximum number of canonical types
 \arg POLYGON1        general polygonal cell (user-defined)
 \arg POLYGON2        general polygonal cell (user-defined)
 \arg POLYGON3        general polygonal cell (user-defined)
 \arg POLYGON4        general polygonal cell (user-defined)
 \arg POLYGON5        general polygonal cell (user-defined)
 \arg POLYHEDRON1     general polyhedral cell (user-defined)
 \arg POLYHEDRON2     general polyhedral cell (user-defined)
 \arg POLYHEDRON3     general polyhedral cell (user-defined)
 \arg POLYHEDRON4     general polyhedral cell (user-defined)
 \arg POLYHEDRON5     general polyhedral cell (user-defined)
 \arg MAXTYPE         placeholder for looping over all types
 

 \section cell_templates_reference Reference cells
 
 Reference cells are used in reconstruction methods based on pullbacks, such as FEM. In such methods, 
 basis set on an ambient cell is obtained by pullback of a basis set defined on a standard reference 
 cell. This approach is practical for cells that have polynomial pullbacks. Intrepid supports pullback
 (FEM) reconstructions for EDGE, TRI, QUAD, TET, HEX, PYRAMID and PRISM cells. 
 
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
          *-------* /      	     *---|---*
         /|       |      	    /    |  /|
        / |   2   |      	   /   5   / |
       /          |     	  /       / ----> 
 <----*   |       |    	     *-------*   |
      | 3 *-------*   		 |       | 1 *
      |  /       /    		 |  /    |  /
      | /   4   /     		 | \/ 0  | /
      |/       /      		 |       |/
      *-----|-*       		 *-------*
            |
            |
            v

@endverbatim
 
 
 
 
 
 
 
 
 
 */
