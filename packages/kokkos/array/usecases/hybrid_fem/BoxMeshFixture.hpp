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

#ifndef KOKKOS_BOXMESHFIXTURE_HPP
#define KOKKOS_BOXMESHFIXTURE_HPP

/* #define KOKKOS_ARRAY_BOUNDS_CHECK 1 */

#include <stdexcept>
#include <sstream>

#include <KokkosArray_View.hpp>
#include <KokkosArray_CrsArray.hpp>

#include <BoxMeshPartition.hpp>
#include <FEMesh.hpp>

//----------------------------------------------------------------------------

template< class FEMeshType >
void
box_mesh_fixture_verify(
  const typename FEMeshType::node_coords_type::HostMirror   & node_coords ,
  const typename FEMeshType::elem_node_ids_type::HostMirror & elem_node_ids ,
  const typename FEMeshType::node_elem_ids_type::HostMirror & node_elem_ids )
{
  typedef typename FEMeshType::size_type size_type ;
  typedef typename FEMeshType::node_coords_type node_coords_type ;
  typedef typename node_coords_type::value_type coords_type ;

  enum { element_node_count = 8 };

  const size_type node_count_total = node_coords.dimension(0);
  const size_type elem_count_total = elem_node_ids.dimension(0);

  for ( size_type node_index = 0 ;
                  node_index < node_count_total ; ++node_index ) {

    for ( size_type
            j = node_elem_ids.row_map[ node_index ] ;
            j < node_elem_ids.row_map[ node_index + 1 ] ; ++j ) {

      const size_type elem_index = node_elem_ids.entries(j,0);
      const size_type node_local = node_elem_ids.entries(j,1);
      const size_type en_id      = elem_node_ids(elem_index,node_local);

      if ( node_index != en_id ) {
        std::ostringstream msg ;
        msg << "BoxMeshFixture node_elem_ids error"
            << " : node_index(" << node_index
            << ") entry(" << j 
            << ") elem_index(" << elem_index
            << ") node_local(" << node_local
            << ") elem_node_id(" << en_id
            << ")" ;
        throw std::runtime_error( msg.str() );
      }
    }
  }

  const coords_type elem_node_local_coord[ element_node_count ][3] =
    { { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 1 , 1 , 0 } , { 0 , 1 , 0 } ,
      { 0 , 0 , 1 } , { 1 , 0 , 1 } , { 1 , 1 , 1 } , { 0 , 1 , 1 } };

  for ( size_type elem_index = 0 ;
                  elem_index < elem_count_total; ++elem_index ) {

    coords_type elem_node_coord[ element_node_count ][3] ;

    for ( size_type nn = 0 ; nn < element_node_count ; ++nn ) {
      const size_type node_index = elem_node_ids( elem_index , nn );

      for ( size_type nc = 0 ; nc < 3 ; ++nc ) {
        elem_node_coord[nn][nc] = node_coords( node_index , nc );
      }
    }

    for ( size_type nn = 0 ; nn < element_node_count ; ++nn ) {
      if ( elem_node_coord[nn][0] !=
           elem_node_coord[ 0][0] + elem_node_local_coord[nn][0] ||
           elem_node_coord[nn][1] !=
           elem_node_coord[ 0][1] + elem_node_local_coord[nn][1] ||
           elem_node_coord[nn][2] !=
           elem_node_coord[ 0][2] + elem_node_local_coord[nn][2] ) {
        throw std::runtime_error(
          std::string("elem_node_coord mapping failure") );
      }
    }
  }
}

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

template < typename Scalar , class Device >
HybridFEM::FEMesh< Scalar , 8 , Device >
box_mesh_fixture( const size_t proc_count ,
                  const size_t proc_local ,
                  const size_t nodes_x ,
                  const size_t nodes_y ,
                  const size_t nodes_z )
{
  enum { element_node_count = 8 };

  typedef typename Device::size_type  size_type ;
  typedef Scalar                      scalar_type ;
  typedef HybridFEM::FEMesh< scalar_type , element_node_count , Device > femesh_type ;

  const size_t elem_node_local_coord[ element_node_count ][3] =
    { { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 1 , 1 , 0 } , { 0 , 1 , 0 } ,
      { 0 , 0 , 1 } , { 1 , 0 , 1 } , { 1 , 1 , 1 } , { 0 , 1 , 1 } };

  const BoxBounds use_boxes ;

  BoxType node_box_global ;
  BoxType node_box_local_used ;

  node_box_global[0][0] = 0 ; node_box_global[0][1] = nodes_x ;
  node_box_global[1][0] = 0 ; node_box_global[1][1] = nodes_y ;
  node_box_global[2][0] = 0 ; node_box_global[2][1] = nodes_z ;

  std::vector< BoxType >              node_box_parts( proc_count );
  std::vector<size_t>                 node_used_id_map ;
  std::vector<size_t>                 node_part_counts ;
  std::vector< std::vector<size_t> >  node_send_map ;
  size_t node_count_interior = 0 ;
  size_t node_count_owned    = 0 ;
  size_t node_count_total    = 0 ;

  box_partition_rcb( node_box_global , node_box_parts );

  box_partition_maps( node_box_global , node_box_parts ,
                      use_boxes ,
                      proc_local ,
                      node_box_local_used ,
                      node_used_id_map ,
                      node_count_interior ,
                      node_count_owned ,
                      node_count_total ,
                      node_part_counts ,
                      node_send_map );

  size_t node_count_send = 0 ;
  for ( size_t i = 0 ; i < node_send_map.size() ; ++i ) {
    node_count_send += node_send_map[i].size();
  }

  BoxType node_box_local_owned = node_box_parts[ proc_local ];

  const long local_elems_x =
    ( node_box_local_used[0][1] - node_box_local_used[0][0] ) - 1 ;
  const long local_elems_y =
    ( node_box_local_used[1][1] - node_box_local_used[1][0] ) - 1 ;
  const long local_elems_z =
    ( node_box_local_used[2][1] - node_box_local_used[2][0] ) - 1 ;

  const size_t elem_count_total = std::max( long(0) , local_elems_x ) *
                                  std::max( long(0) , local_elems_y ) *
                                  std::max( long(0) , local_elems_z );

  const long interior_elems_x =
    ( node_box_local_owned[0][1] - node_box_local_owned[0][0] ) - 1 ;
  const long interior_elems_y =
    ( node_box_local_owned[1][1] - node_box_local_owned[1][0] ) - 1 ;
  const long interior_elems_z =
    ( node_box_local_owned[2][1] - node_box_local_owned[2][0] ) - 1 ;

  const size_t elem_count_interior = std::max( long(0) , interior_elems_x ) *
                                     std::max( long(0) , interior_elems_y ) *
                                     std::max( long(0) , interior_elems_z );

  size_t recv_msg_count = 0 ;
  size_t send_msg_count = 0 ;
  size_t send_count = 0 ;

  for ( size_t i = 1 ; i < proc_count ; ++i ) {
    if ( node_part_counts[i] ) ++recv_msg_count ;
    if ( node_send_map[i].size() ) {
      ++send_msg_count ;
      send_count += node_send_map[i].size();
    }
  }

  femesh_type mesh ;

  typedef typename femesh_type::node_coords_type    node_coords_type ;
  typedef typename femesh_type::elem_node_ids_type  elem_node_ids_type ;
  typedef typename femesh_type::node_elem_ids_type  node_elem_ids_type ;

  if ( node_count_total ) {
    mesh.node_coords =
      KokkosArray::create< node_coords_type >( std::string("node_coords"), node_count_total );
  }

  if ( elem_count_total ) {
    mesh.elem_node_ids =
      KokkosArray::create< elem_node_ids_type >( std::string("elem_node_ids"), elem_count_total );
  }

  mesh.parallel_data_map.assign( node_count_interior ,
                                 node_count_owned ,
                                 node_count_total ,
                                 recv_msg_count ,
                                 send_msg_count ,
                                 send_count );

  typename node_coords_type::HostMirror node_coords =
    KokkosArray::create_mirror( mesh.node_coords );

  typename elem_node_ids_type::HostMirror elem_node_ids =
    KokkosArray::create_mirror( mesh.elem_node_ids );

  //------------------------------------
  // node coordinates of grid.

  for ( size_t iz = node_box_local_used[2][0] ;
               iz < node_box_local_used[2][1] ; ++iz ) {

  for ( size_t iy = node_box_local_used[1][0] ;
               iy < node_box_local_used[1][1] ; ++iy ) {

  for ( size_t ix = node_box_local_used[0][0] ;
               ix < node_box_local_used[0][1] ; ++ix ) {
    const size_t node_local_id =
      box_map_id( node_box_local_used , node_used_id_map , ix , iy , iz );
    node_coords( node_local_id , 0 ) = ix ;
    node_coords( node_local_id , 1 ) = iy ;
    node_coords( node_local_id , 2 ) = iz ;
  }}}

  //------------------------------------
  // Initialize element-node connectivity:
  // Order elements that only depend on owned nodes first.
  // These elements could be computed while waiting for
  // received node data.

  size_t elem_index_interior = 0 ;
  size_t elem_index_boundary = elem_count_interior ;

  for ( size_t iz = node_box_local_used[2][0] ;
               iz < node_box_local_used[2][1] - 1 ; ++iz ) {
  for ( size_t iy = node_box_local_used[1][0] ;
               iy < node_box_local_used[1][1] - 1 ; ++iy ) {
  for ( size_t ix = node_box_local_used[0][0] ;
               ix < node_box_local_used[0][1] - 1 ; ++ix ) {

    size_t elem_index ;

    // If lower and upper nodes are owned then element is interior
    if ( contain( node_box_local_owned, ix,   iy,   iz ) &&
         contain( node_box_local_owned, ix+1, iy+1, iz+1 ) ) {
      elem_index = elem_index_interior++ ;
    }
    else {
      elem_index = elem_index_boundary++ ;
    }

    for ( size_t nn = 0 ; nn < element_node_count ; ++nn ) {
      const size_t jx = ix + elem_node_local_coord[nn][0] ;
      const size_t jy = iy + elem_node_local_coord[nn][1] ;
      const size_t jz = iz + elem_node_local_coord[nn][2] ;

      const size_t node_local_id =
        box_map_id( node_box_local_used , node_used_id_map , jx , jy , jz );

      elem_node_ids( elem_index , nn ) = node_local_id ;
    }
  }}}

  //------------------------------------
  // Populate node->element connectivity:

  std::vector<size_t> node_elem_work( node_count_total , (size_t) 0 );

  for ( size_t i = 0 ; i < elem_count_total ; ++i ) {
    for ( size_t n = 0 ; n < element_node_count  ; ++n ) {
      ++node_elem_work[ elem_node_ids(i,n) ];
    }
  }

  mesh.node_elem_ids =
    KokkosArray::create_crsarray< node_elem_ids_type >( node_elem_work );

  typename node_elem_ids_type::HostMirror
    node_elem_ids = KokkosArray::create_mirror( mesh.node_elem_ids );

  for ( size_t i = 0 ; i < node_count_total ; ++i ) {
    node_elem_work[i] = node_elem_ids.row_map[i];
  }

  // Looping in element order insures the list of elements
  // is sorted by element index.

  for ( size_t i = 0 ; i < elem_count_total ; ++i ) {
    for ( size_t n = 0 ; n < element_node_count ; ++n ) {
      const size_type nid = elem_node_ids(i, n);
      const size_type j = node_elem_work[nid] ; ++node_elem_work[nid] ;

      node_elem_ids.entries( j , 0 ) = i ;
      node_elem_ids.entries( j , 1 ) = n ;
    }
  }
  //------------------------------------

  box_mesh_fixture_verify<femesh_type>( node_coords , elem_node_ids , node_elem_ids );

  KokkosArray::deep_copy( mesh.node_coords ,   node_coords );
  KokkosArray::deep_copy( mesh.elem_node_ids , elem_node_ids );
  KokkosArray::deep_copy( mesh.node_elem_ids , node_elem_ids );

  //------------------------------------
  // Communication lists:
  {
    recv_msg_count = 0 ;
    send_msg_count = 0 ;
    send_count = 0 ;

    for ( size_t i = 1 ; i < proc_count ; ++i ) {

      // Order sending starting with the local processor rank 
      // to try to smooth out the amount of messages simultaneously
      // send to a particular processor.

      const int proc = ( proc_local + i ) % proc_count ;
      if ( node_part_counts[i] ) {
        mesh.parallel_data_map.host_recv(recv_msg_count,0) = proc ;
        mesh.parallel_data_map.host_recv(recv_msg_count,1) = node_part_counts[i] ;
        ++recv_msg_count ;
      }
      if ( node_send_map[i].size() ) {
        mesh.parallel_data_map.host_send(send_msg_count,0) = proc ;
        mesh.parallel_data_map.host_send(send_msg_count,1) = node_send_map[i].size() ;
        for ( size_t j = 0 ; j < node_send_map[i].size() ; ++j , ++send_count ) {
          mesh.parallel_data_map.host_send_item(send_count) = node_send_map[i][j] - node_count_interior ;
        }
        ++send_msg_count ;
      }
    }
  }

  return mesh ;
};

#endif /* #ifndef KOKKOS_BOXMESHFIXTURE_HPP */


