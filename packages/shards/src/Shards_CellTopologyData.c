/*
//@HEADER
// ************************************************************************
//
//                Shards : Shared Discretization Tools
//                 Copyright 2008 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// Export of this program may require a license from the United States
// Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Carter Edwards (hcedwar@sandia.gov),
//                    Pavel Bochev (pbboche@sandia.gov), or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
//@HEADER
*/

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

