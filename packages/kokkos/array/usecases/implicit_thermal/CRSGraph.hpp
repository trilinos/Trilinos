/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <algorithm>
#include <vector>
#include <set>

template< class int_mdarray, class int_multivector>
void init_crsgraph(
  int_mdarray node_elem_offset ,
  int_mdarray node_elem_ids ,
  int_mdarray elem_node_ids ,
  int_mdarray elem_graph_col ,
  int_multivector & graph_col_offset ,
  int_multivector & graph_col_index )
{
  typedef typename int_mdarray::value_type index_type ;
  const int nelem = elem_node_ids.dimension(0);
  const int nnode = node_elem_offset.dimension(0) - 1 ;
  const int nnode_per_elem = elem_node_ids.dimension(1);

  //----------------------------------

  graph_col_offset = KokkosArray::create_multivector< int_multivector >(nnode + 1);

  graph_col_offset(0) = 0 ;

  std::vector< int > col_ids ;

  // Loop over nodes (rows) and elements of each node

  for(int i = 0; i < nnode; i++){
    const int offset_begin = node_elem_offset(i);
    const int offset_end   = node_elem_offset(i+1);

    std::set<int> node_node_ids ;

    for ( int j = offset_begin ; j < offset_end ; ++j ) {

      const int elem_index = node_elem_ids(j,0);

      for ( int k = 0 ; k < nnode_per_elem ; ++k ) {
        const int node_index = elem_node_ids( elem_index , k );

        node_node_ids.insert( node_index );
      }
    }

    graph_col_offset(i+1) = graph_col_offset(i) + node_node_ids.size();

    col_ids.insert( col_ids.end(), node_node_ids.begin(), node_node_ids.end() );
  }

  graph_col_index = KokkosArray::create_multivector< int_multivector >( col_ids.size() );

  for ( size_t i = 0 ; i < col_ids.size() ; ++i ) {
    graph_col_index(i) = col_ids[i] ;
  }

  //----------------------------------
  // Mapping of ( element , row_node , col_node ) -> CRS location

  for ( int elem_id = 0 ; elem_id < nelem ; ++elem_id ) {
    for ( int i = 0 ; i < nnode_per_elem ; ++i ) {
      const index_type node_row = elem_node_ids( elem_id , i );

      const int col_end = graph_col_offset( node_row + 1 );
      const int col_beg = graph_col_offset( node_row );

      for ( int j = 0 ; j < nnode_per_elem ; ++j ) {
        const index_type node_col = elem_node_ids( elem_id , j );

        int col_search = col_beg ;

        for ( int len = col_end - col_beg ; 0 < len ; ) {

          const int half = len >> 1;
          const int middle = col_search + half ;

          if ( graph_col_index(middle) < node_col ){
            col_search = middle + 1 ;
            len -= half + 1 ;
          }
          else {
            len = half ;
          }
        }

        elem_graph_col( elem_id , i , j ) = col_search ;
      }
    }
  }
}

