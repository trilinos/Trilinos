/*------------------------------------------------------------------------*/
/*                shards : Shared Discretization Tools                    */
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
/* Questions? Contact Pavel Bochev      (pbboche@sandia.gov)              */
/*                    H. Carter Edwards (hcedwar@sandia.gov)              */
/*                    Denis Ridzal      (dridzal@sandia.gov).             */
/*------------------------------------------------------------------------*/

#include <Shards_CellTopologyData.h>

#if defined( __cplusplus )
extern "C" {
#endif

int mapCellFaceEdge( const CellTopologyData * cell_topology ,
                     unsigned face_ordinal ,
                     unsigned face_edge_ordinal )
{
  int edge = -1 ;

  if ( cell_topology && face_ordinal < cell_topology->subcell_count[2] ) {
    const CellTopologyData * const face_top =
      cell_topology->subcell[2][face_ordinal].topology ;

    if ( face_edge_ordinal < face_top->subcell_count[1] ) {

      const unsigned face_edge_node_0 =
        face_top->edge[face_edge_ordinal].node[0];
      const unsigned face_edge_node_1 =
        face_top->edge[face_edge_ordinal].node[1];
       
      const unsigned cell_face_edge_node_0 =
        cell_topology->subcell[2][face_ordinal].node[face_edge_node_0];
      const unsigned cell_face_edge_node_1 =
        cell_topology->subcell[2][face_ordinal].node[face_edge_node_1];

      const int edge_count = cell_topology->subcell_count[1] ;

      for ( edge = 0 ; edge < edge_count ; ++edge ) {
        const unsigned cell_edge_node_0 = cell_topology->edge[edge].node[0];
        const unsigned cell_edge_node_1 = cell_topology->edge[edge].node[1];
                         
        if ( ( cell_face_edge_node_0 == cell_edge_node_0 &&
               cell_face_edge_node_1 == cell_edge_node_1 ) ||
             ( cell_face_edge_node_0 == cell_edge_node_1 &&
               cell_face_edge_node_1 == cell_edge_node_0 ) ) {
          break ;
        }
      }
      if ( edge_count == edge ) { edge = -1 ; }
    }
  }
  return edge ;
}

#if defined( __cplusplus )
} /* extern "C" */
#endif

