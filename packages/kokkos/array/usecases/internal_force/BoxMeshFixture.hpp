/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_BOXMESHFIXTURE_HPP
#define KOKKOS_BOXMESHFIXTURE_HPP

#include <stdexcept>
#include <Kokkos_MDArrayView.hpp>

//  construct a structured, rectangular prism mesh of Hex elements,
//  with dimensions given by elems_x, elems_y, elems_z

template < class int_mdarray , class scalar_mdarray >
class BoxMeshFixture {
public:
  int_mdarray    elem_node_ids ;
  int_mdarray    node_elem_offset ;
  int_mdarray    node_elem_ids ;
  scalar_mdarray node_coords ;

  const int elem_count_x ;
  const int elem_count_y ;
  const int elem_count_z ;

  const int nelems;
  const int nnodes;

  //        7 -------------- 6
  //       /|               /|
  //      / |              / |
  //     /  |             /  |
  //    /   |            /   |
  //   4 -------------- 5    |
  //   |    |           |    |
  //   |    |           |    |
  //   |    3 --------- | -- 2
  //   |   /            |   /
  //   |  /             |  /
  //   | /              | /
  //   |/               |/
  //   0 -------------- 1
  //
  //   Coordinate system:
  //     X = grid-i
  //     Y = grid-j
  //     Z = grid-k
  //
  //   Z   Y
  //   |  /
  //   | /
  //   |/
  //   O-----X
  //

  void verify_connectivity_and_coordinates() const
  {
    for ( int i = 0 ; i < nnodes ; ++i ) {
      for ( int j = node_elem_offset(i) ;
                j < node_elem_offset(i+1) ; ++j ) {
        const int elem_index = node_elem_ids(j,0);
        const int node_local = node_elem_ids(j,1);
        const int en_id = elem_node_ids( elem_index , node_local );
        if ( i != en_id ) {
          throw std::runtime_error( std::string("node_elem_ids mapping failure") );
        }
      }
    }

    for(int k = 0; k < elem_count_z; k++){
      for(int j = 0; j < elem_count_y; j++){
        for(int i = 0; i < elem_count_x; i++){

          const int elem_index = k * elem_count_x * elem_count_y +
                                 j * elem_count_x + i;

          int elem_node_coords[3][8];

          elem_node_coords[ 0 ][ 0 ] = i;
          elem_node_coords[ 1 ][ 0 ] = j;
          elem_node_coords[ 2 ][ 0 ] = k ;

          elem_node_coords[ 0 ][ 1 ] = i + 1;
          elem_node_coords[ 1 ][ 1 ] = j;
          elem_node_coords[ 2 ][ 1 ] = k ;

          elem_node_coords[ 0 ][ 2 ] = i + 1;
          elem_node_coords[ 1 ][ 2 ] = j + 1;
          elem_node_coords[ 2 ][ 2 ] = k;

          elem_node_coords[ 0 ][ 3 ] = i;
          elem_node_coords[ 1 ][ 3 ] = j + 1 ;
          elem_node_coords[ 2 ][ 3 ] = k;

          elem_node_coords[ 0 ][ 4 ] = i;
          elem_node_coords[ 1 ][ 4 ] = j;
          elem_node_coords[ 2 ][ 4 ] = k + 1;

          elem_node_coords[ 0 ][ 5 ] = i + 1;
          elem_node_coords[ 1 ][ 5 ] = j;
          elem_node_coords[ 2 ][ 5 ] = k + 1;

          elem_node_coords[ 0 ][ 6 ] = i + 1;
          elem_node_coords[ 1 ][ 6 ] = j + 1;
          elem_node_coords[ 2 ][ 6 ] = k + 1;

          elem_node_coords[ 0 ][ 7 ] = i;
          elem_node_coords[ 1 ][ 7 ] = j + 1;
          elem_node_coords[ 2 ][ 7 ] = k + 1 ;


          for ( int nn = 0 ; nn < 8 ; ++nn ) {
            for ( int nc = 0 ; nc < 3 ; ++nc ) {
              const int en_coord = elem_node_coords[ nc ][ nn ];
              const int en_id    = elem_node_ids( elem_index , nn );
              const int n_coord  = (int) node_coords( en_id , nc );

              if ( en_coord != n_coord ) {
                throw std::runtime_error( std::string("elem_node_coord mapping failure") );
              }
            }
          }
        }
      }
    }
  }

  BoxMeshFixture( const int elems_x ,
                  const int elems_y ,
                  const int elems_z )
  : elem_node_ids()
  , node_elem_offset()
  , node_elem_ids()
  , node_coords()
  , elem_count_x( elems_x )
  , elem_count_y( elems_y )
  , elem_count_z( elems_z )
  , nelems(elems_x*elems_y*elems_z)
  , nnodes( (elems_x+1)*(elems_y+1)*(elems_z+1))
  {
    const int nx = elems_x + 1;
    const int ny = elems_y + 1;
    const int nz = elems_z + 1;

    elem_node_ids    = Kokkos::create_mdarray< int_mdarray    >(nelems, 8);
    node_coords      = Kokkos::create_mdarray< scalar_mdarray >(nnodes, 3);

    // Initialize node coordinates of grid.

    for(int k = 0; k < nz ; k++){
      for(int j = 0; j < ny ; j++){
        for(int i = 0; i < nx ; i++){
          const int node_index  = k * nx * ny + j * nx + i;
          node_coords(node_index,0) = i ;
          node_coords(node_index,1) = j ;
          node_coords(node_index,2) = k ;
        }
      }
    }

    // Initialize element-node connectivity:

    for(int k = 0; k < elems_z; k++){
      for(int j = 0; j < elems_y; j++){
        for(int i = 0; i < elems_x; i++){

          const int elem_index = k * elems_x * elems_y + j * elems_x + i;

          const int node0_index = k * nx * ny + j * nx + i;

          const int node_indices[8] = {
            node0_index ,
            node0_index + 1 ,
            node0_index + 1 + nx ,
            node0_index     + nx ,

            node0_index          + nx * ny ,
            node0_index + 1      + nx * ny ,
            node0_index + 1 + nx + nx * ny ,
            node0_index     + nx + nx * ny };

          for ( int n = 0 ; n < 8 ; ++n ) {
            elem_node_ids(elem_index, n) = node_indices[n] ;
          }
        }
      }
    }

    // node-element connectivity:

    node_elem_offset = Kokkos::create_mdarray< int_mdarray >(nnodes + 1);

    int_mdarray node_elem_count = Kokkos::create_mdarray< int_mdarray >( nnodes );

    for(int i = 0; i < nnodes + 1 ; i++){
      node_elem_count( i ) = 0 ;
    }

    for ( int i = 0 ; i < nelems ; ++i ) {
      for ( int n = 0 ; n < 8 ; ++n ) {
        node_elem_count( elem_node_ids(i,n) )++ ;
      }
    }

    node_elem_offset(0) = 0 ;
    for(int i = 0; i < nnodes; ++i ){
      node_elem_offset( i + 1 ) = node_elem_offset(i) + node_elem_count(i);
      node_elem_count( i ) = 0 ;
    }

    const int count_node_elem = node_elem_offset( nnodes );

    node_elem_ids = Kokkos::create_mdarray< int_mdarray >( count_node_elem , 2 );

    for ( int i = 0 ; i < nelems ; ++i ) {
      for ( int n = 0 ; n < 8 ; ++n ) {
        const int nid = elem_node_ids(i, n);
        const int j = node_elem_offset(nid) + node_elem_count(nid);

        node_elem_ids( j , 0 ) = i ;
        node_elem_ids( j , 1 ) = n ;

        ++node_elem_count( nid );
      }
    }
  }
};

#endif /* #ifndef KOKKOS_BOXMESHFIXTURE_HPP */

