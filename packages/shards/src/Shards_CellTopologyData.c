// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

