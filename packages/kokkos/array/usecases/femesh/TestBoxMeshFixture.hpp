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

#ifndef TESTFEMESHBOXFIXTURE_HPP
#define TESTFEMESHBOXFIXTURE_HPP

#include <stdio.h>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>
#include <BoxMeshFixture.hpp>

#include <ParallelDistributedComm.hpp>

//----------------------------------------------------------------------------

namespace TestFEMesh {

template< class FEMeshType >
struct VerifyPack ;

template< class FEMeshType >
struct VerifyUnpack  ;

}

//----------------------------------------------------------------------------

#ifdef HAVE_MPI

namespace TestFEMesh {

template< typename identifier_integer_type ,
          typename coordinate_scalar_type ,
          unsigned ElemNodeCount ,
          class Device >
void verify_parallel(
  comm::Machine machine ,
  const FEMeshFixture< identifier_integer_type ,
                       coordinate_scalar_type ,
                       ElemNodeCount ,
                       Device > & fixture )
{
  typedef FEMeshFixture< identifier_integer_type ,
                         coordinate_scalar_type ,
                         ElemNodeCount ,
                         Device > femesh_type ;

  const size_t chunk_size = 3 ;

  // Communicate node coordinates to verify communication and setup.

  AsyncExchange< coordinate_scalar_type , Device >
    exchange( machine , fixture.node_part , fixture.node_send , chunk_size );

  Kokkos::parallel_for( exchange.send_count , VerifyPack<femesh_type>( fixture , exchange.send_buffer ) );

  exchange.send();

  // Could do something else here ...

  // Wait for incoming data

  exchange.wait_receive();

  // Unpack recv buffer

  unsigned local[2] ;

  local[0] = exchange.recv_count ;
  local[1] =
    Kokkos::parallel_reduce( exchange.recv_count ,
      VerifyUnpack< femesh_type >( fixture , exchange.recv_buffer , exchange.host_recv_part.row_entry_begin(1) ) );

  unsigned global[2] = { 0 , 0 };

  MPI_Allreduce( local , global ,
                 2 , MPI_UNSIGNED , MPI_SUM , machine.mpi_comm );

  if ( 0 == comm::rank( machine ) ) {
    std::cout << "TestFEMesh::verify_parallel "
              << "NP(" << comm::size( machine )
              << ") verified_nodes(" << global[0]
              << ") failed_nodes(" << global[1]
              << ")" << std::endl ;
  }

  if ( global[1] ) {
    throw std::runtime_error( std::string("coordinate exchange failed") );
  }
}

} // namespace TestFEMesh

#else /* ! #ifdef HAVE_MPI */

namespace TestFEMesh {

template< typename identifier_integer_type ,
          typename coordinate_scalar_type ,
          unsigned ElemNodeCount ,
          class Device >
void verify_parallel(
  comm::Machine ,
  const FEMeshFixture< identifier_integer_type ,
                       coordinate_scalar_type ,
                       ElemNodeCount ,
                       Device > & )
{
}

} // namespace TestFEMesh

#endif /* ! #ifdef HAVE_MPI */

//----------------------------------------------------------------------------

template< class Device >
void test_box_fixture( comm::Machine machine ,
                       const size_t nodes_nx ,
                       const size_t nodes_ny ,
                       const size_t nodes_nz )
{
  typedef int coordinate_scalar_type ;
  typedef BoxMeshFixture< coordinate_scalar_type, Device > box_mesh_type ;
  typedef typename box_mesh_type::fixture_dev_type mesh_fixture_type ;

  const size_t proc_count = comm::size( machine );
  const size_t proc_local = comm::rank( machine ) ;

  const box_mesh_type
    box_mesh( proc_count, proc_local, nodes_nx, nodes_ny, nodes_nz );

  const mesh_fixture_type & fixture = box_mesh ;

  TestFEMesh::verify_parallel( machine , fixture );
}

#endif /* #ifndef TESTFEMESHBOXFIXTURE_HPP */

//----------------------------------------------------------------------------

namespace TestFEMesh {

template< typename identifier_integer_type ,
          typename coordinate_scalar_type ,
          unsigned ElemNodeCount >
struct VerifyPack
  < FEMeshFixture< identifier_integer_type ,
                   coordinate_scalar_type ,
                   ElemNodeCount ,
                   KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE              device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef Kokkos::Impl::MemoryView< coordinate_scalar_type , device_type > buffer_type ;

  typedef FEMeshFixture< identifier_integer_type ,
                         coordinate_scalar_type ,
                         ElemNodeCount ,
                         device_type >
    femesh_type ;

  const femesh_type   femesh ;
  const buffer_type   send_buffer ;

  VerifyPack(
    const femesh_type & arg_femesh ,
    const buffer_type & arg_send_buffer )
  : femesh( arg_femesh )
  , send_buffer( arg_send_buffer )
  {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i ) const
  {
    size_type k = i * 3 ;
    const size_type node_id = femesh.node_send(i);
    send_buffer[k]   = femesh.node_coords(node_id,0);
    send_buffer[k+1] = femesh.node_coords(node_id,1);
    send_buffer[k+2] = femesh.node_coords(node_id,2);
  }
};

template< typename identifier_integer_type ,
          typename coordinate_scalar_type ,
          unsigned ElemNodeCount >
struct VerifyUnpack 
  < FEMeshFixture< identifier_integer_type ,
                   coordinate_scalar_type ,
                   ElemNodeCount ,
                   KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE              device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef size_type                        value_type ;
  typedef Kokkos::Impl::MemoryView< coordinate_scalar_type , device_type > buffer_type ;

  typedef FEMeshFixture< identifier_integer_type ,
                         coordinate_scalar_type ,
                         ElemNodeCount ,
                         device_type >
    femesh_type ;

  const femesh_type femesh ;
  const buffer_type recv_buffer ;
  const size_type  node_id_begin ;

  VerifyUnpack(
    const femesh_type & arg_femesh ,
    const buffer_type & arg_recv_buffer ,
    const size_type     arg_node_id_begin )
  : femesh( arg_femesh )
  , recv_buffer( arg_recv_buffer )
  , node_id_begin( arg_node_id_begin )
  {}

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source ; }

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const size_type i , value_type & update ) const
  {
    const size_type node_id = i + node_id_begin ;
    const size_type k = i * 3 ;

    const coordinate_scalar_type x = recv_buffer[k];
    const coordinate_scalar_type y = recv_buffer[k+1];
    const coordinate_scalar_type z = recv_buffer[k+2];

    if ( x != femesh.node_coords(node_id,0) ||
         y != femesh.node_coords(node_id,1) ||
         z != femesh.node_coords(node_id,2) ) {
      printf("TestFEMesh::VerifyUnpack failed at node %d\n",(int)node_id);
      ++update ;
    }
  }
};

}

