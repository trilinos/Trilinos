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

#include <stdexcept>
#include <sstream>

#include <Kokkos_CrsArray.hpp>
#include <Kokkos_MDArray.hpp>
#include <Kokkos_MultiVector.hpp>

#include <BoxMeshPartition.hpp>

//----------------------------------------------------------------------------


template< typename IntType ,
          typename FloatType ,
          unsigned ElemNodeCount ,
          class Device >
struct FEMeshFixture {

  static const IntType element_node_count = ElemNodeCount ;

  // typedef Array< FloatType[3] , Device > node_coords_type ;
  // typedef Array< IntType[ElemNodeCount] , Device > elem_node_ids_type ;

  typedef Kokkos::MDArray<  FloatType ,  Device >  node_coords_type ;
  typedef Kokkos::MDArray<  IntType ,    Device >  elem_node_ids_type ;
  typedef Kokkos::CrsArray< IntType[2] , Device >  node_elem_ids_type ;
  typedef Kokkos::CrsArray< void ,       Device >  node_part_type ;
  typedef Kokkos::CrsArray< unsigned ,   Device >  node_send_type ;

  node_coords_type    node_coords ;
  elem_node_ids_type  elem_node_ids ;
  node_elem_ids_type  node_elem_ids ;
  node_part_type      node_part ;
  node_send_type      node_send ;
};

//----------------------------------------------------------------------------
//  construct a structured, rectangular prism mesh of Hex elements,
//  with dimensions given by elems_x, elems_y, elems_z

template < typename Scalar , class Device >
class BoxMeshFixture :
  public FEMeshFixture< typename Device::size_type /* local indices */ ,
                        Scalar                     /* coordinate value */ ,
                        8                          /* nodes per element */ ,
                        Device >
{
public:
  typedef typename Device::size_type  index_type ;
  typedef Scalar                      scalar_type ;

  static const size_t element_node_count = 8 ;

  typedef FEMeshFixture< index_type, scalar_type, element_node_count, Device >
    fixture_dev_type ;

  typedef typename fixture_dev_type::node_coords_type    node_coords_type ;
  typedef typename fixture_dev_type::elem_node_ids_type  elem_node_ids_type ;
  typedef typename fixture_dev_type::node_elem_ids_type  node_elem_ids_type ;
  typedef typename fixture_dev_type::node_part_type      node_part_type ;
  typedef typename fixture_dev_type::node_send_type      node_send_type ;

  using fixture_dev_type::node_coords ;
  using fixture_dev_type::elem_node_ids ;
  using fixture_dev_type::node_elem_ids ;
  using fixture_dev_type::node_part ;
  using fixture_dev_type::node_send ;

  typedef FEMeshFixture< index_type , scalar_type, element_node_count,
                         typename Kokkos::HostMapped< Device >::type >
    fixture_host_type ;
  
  BoxType node_box_global ;
  BoxType node_box_local_owned ;
  BoxType node_box_local_used ;

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

  static void
  verify_connectivity_and_coordinates( const fixture_host_type & h_mesh )
  {
    const size_t node_count = h_mesh.node_coords.dimension(0);
    const size_t elem_count = h_mesh.elem_node_ids.dimension(0);

    for ( size_t node_index = 0 ; node_index < node_count; ++node_index ) {

      for ( index_type
              j = h_mesh.node_elem_ids.row_entry_begin( node_index ) ;
              j < h_mesh.node_elem_ids.row_entry_end( node_index ) ; ++j ) {

        const index_type elem_index = h_mesh.node_elem_ids(j,0);
        const index_type node_local = h_mesh.node_elem_ids(j,1);
        const index_type en_id = h_mesh.elem_node_ids(elem_index,node_local);

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

    const size_t elem_node_local_coord[ element_node_count ][3] =
      { { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 1 , 1 , 0 } , { 0 , 1 , 0 } ,
        { 0 , 0 , 1 } , { 1 , 0 , 1 } , { 1 , 1 , 1 } , { 0 , 1 , 1 } };

    for ( size_t elem_index = 0 ; elem_index < elem_count; ++elem_index ) {
      size_t elem_node_coord[ element_node_count ][3] ;

      for ( index_type nn = 0 ; nn < element_node_count ; ++nn ) {
        const index_type node_index = h_mesh.elem_node_ids( elem_index , nn );

        for ( index_type nc = 0 ; nc < 3 ; ++nc ) {
          elem_node_coord[nn][nc] = h_mesh.node_coords( node_index , nc );
        }
      }

      for ( index_type nn = 0 ; nn < element_node_count ; ++nn ) {
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

private:

  void populate_node_element( fixture_host_type & h_mesh )
  {
    const size_t node_count = h_mesh.node_coords.dimension(0);
    const size_t elem_count = h_mesh.elem_node_ids.dimension(0);

    std::vector<size_t> node_elem_work( node_count , (size_t) 0 );

    for ( size_t i = 0 ; i < elem_count ; ++i ) {
      for ( index_type n = 0 ; n < element_node_count  ; ++n ) {
        ++node_elem_work[ h_mesh.elem_node_ids(i,n) ];
      }
    }

    node_elem_ids =
      Kokkos::create_crsarray< node_elem_ids_type >( node_elem_work );

    h_mesh.node_elem_ids =
      Kokkos::create_mirror( node_elem_ids , Kokkos::Impl::MirrorUseView() );

    for ( size_t i = 0 ; i < node_count ; ++i ) {
      node_elem_work[i] = h_mesh.node_elem_ids.row_entry_begin(i);
    }

    // Looping in element order insures the list of elements
    // is sorted by element index.

    for ( size_t i = 0 ; i < elem_count ; ++i ) {
      for ( size_t n = 0 ; n < element_node_count ; ++n ) {
        const index_type nid = h_mesh.elem_node_ids(i, n);
        const index_type j = node_elem_work[nid] ; ++node_elem_work[nid] ;

        h_mesh.node_elem_ids( j , 0 ) = i ;
        h_mesh.node_elem_ids( j , 1 ) = n ;
      }
    }

    Kokkos::deep_copy( node_elem_ids , h_mesh.node_elem_ids );
  }

public:

  BoxMeshFixture( const size_t proc_count ,
                  const size_t proc_local ,
                  const size_t nodes_x ,
                  const size_t nodes_y ,
                  const size_t nodes_z )
  : FEMeshFixture< index_type, scalar_type, element_node_count, Device >()
  , node_box_global()
  , node_box_local_owned()
  , node_box_local_used()
  {
    const size_t ghost_layer = 1 ;

    const size_t elem_node_local_coord[ element_node_count ][3] =
      { { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 1 , 1 , 0 } , { 0 , 1 , 0 } ,
        { 0 , 0 , 1 } , { 1 , 0 , 1 } , { 1 , 1 , 1 } , { 0 , 1 , 1 } };


    node_box_global[0][0] = 0 ; node_box_global[0][1] = nodes_x ;
    node_box_global[1][0] = 0 ; node_box_global[1][1] = nodes_y ;
    node_box_global[2][0] = 0 ; node_box_global[2][1] = nodes_z ;

    std::vector< BoxType >              node_box_parts( proc_count );
    std::vector<size_t>                 node_used_id_map ;
    std::vector<size_t>                 node_part_counts ;
    std::vector< std::vector<size_t> >  node_send_map ;
    size_t node_count_interior ;
    size_t node_count_owned ;
    size_t node_count_used ;

    box_partition_rcb( node_box_global , node_box_parts );

    box_partition_maps( node_box_global , node_box_parts , ghost_layer ,
                        proc_local ,
                        node_box_local_used ,
                        node_used_id_map ,
                        node_count_interior ,
                        node_count_owned ,
                        node_count_used ,
                        node_part_counts ,
                        node_send_map );

    node_box_local_owned = node_box_parts[ proc_local ];

    const size_t node_count = count( node_box_local_used );

    const size_t elem_count_owned =
      std::max( long(0),( long( node_box_local_owned[0][1] - node_box_local_owned[0][0] ) - 1 )) *
      std::max( long(0),( long( node_box_local_owned[1][1] - node_box_local_owned[1][0] ) - 1 )) *
      std::max( long(0),( long( node_box_local_owned[2][1] - node_box_local_owned[2][0] ) - 1 ));

    const size_t elem_count =
      std::max( long(0),( long( node_box_local_used[0][1] - node_box_local_used[0][0] ) - 1 )) *
      std::max( long(0),( long( node_box_local_used[1][1] - node_box_local_used[1][0] ) - 1 )) *
      std::max( long(0),( long( node_box_local_used[2][1] - node_box_local_used[2][0] ) - 1 ));

    node_coords   = Kokkos::create_mdarray< node_coords_type >( node_count, 3);
    elem_node_ids = Kokkos::create_mdarray< elem_node_ids_type >( elem_count,  element_node_count );
    
    node_part  = Kokkos::create_crsarray< node_part_type >( node_part_counts );
    node_send  = Kokkos::create_crsarray< node_send_type >( node_send_map );

    fixture_host_type h_mesh ;

    h_mesh.node_coords   = Kokkos::create_mirror( node_coords , Kokkos::Impl::MirrorUseView() );
    h_mesh.elem_node_ids = Kokkos::create_mirror( elem_node_ids , Kokkos::Impl::MirrorUseView() );

    // Initialize node coordinates of grid.

    for ( size_t iz = node_box_local_used[2][0] ;
                 iz < node_box_local_used[2][1] ; ++iz ) {

    for ( size_t iy = node_box_local_used[1][0] ;
                 iy < node_box_local_used[1][1] ; ++iy ) {

    for ( size_t ix = node_box_local_used[0][0] ;
                 ix < node_box_local_used[0][1] ; ++ix ) {
      const size_t node_local_id =
        box_map_id( node_box_local_used , node_used_id_map , ix , iy , iz );
      h_mesh.node_coords( node_local_id , 0 ) = ix ;
      h_mesh.node_coords( node_local_id , 1 ) = iy ;
      h_mesh.node_coords( node_local_id , 2 ) = iz ;
    }}}

    // Initialize element-node connectivity:
    // Order elements that only depend on owned nodes first.
    // These elements could be computed while waiting for
    // received node data.

    size_t elem_index_interior = 0 ;
    size_t elem_index_boundary = elem_count_owned ;

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

      for ( index_type nn = 0 ; nn < element_node_count ; ++nn ) {
        const size_t jx = ix + elem_node_local_coord[nn][0] ;
        const size_t jy = iy + elem_node_local_coord[nn][1] ;
        const size_t jz = iz + elem_node_local_coord[nn][2] ;

        const size_t node_local_id =
          box_map_id( node_box_local_used , node_used_id_map , jx , jy , jz );

        h_mesh.elem_node_ids( elem_index , nn ) = node_local_id ;
      }
    }}}

    populate_node_element( h_mesh );

    verify_connectivity_and_coordinates( h_mesh );

    Kokkos::deep_copy( node_coords ,   h_mesh.node_coords );
    Kokkos::deep_copy( elem_node_ids , h_mesh.elem_node_ids );
    Kokkos::deep_copy( node_elem_ids , h_mesh.node_elem_ids );
  }
};

#endif /* #ifndef KOKKOS_BOXMESHFIXTURE_HPP */

