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

namespace TestFEMesh {

template< typename coordinate_scalar_type >
struct VerifyUnpack<
  KokkosArray::Array< coordinate_scalar_type[3] , KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE              device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef size_type                        value_type ;
  typedef KokkosArray::View< coordinate_scalar_type[] , device_type > buffer_type ;

  typedef KokkosArray::Array< coordinate_scalar_type[3] , device_type > array_type ;

private:

  array_type  node_coords ;
  buffer_type buffer ;
  size_type   node_begin ;

public:

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
    const size_type node_id = i + node_begin ;
    const size_type k = i * 3 ;

    const coordinate_scalar_type x = buffer[k];
    const coordinate_scalar_type y = buffer[k+1];
    const coordinate_scalar_type z = buffer[k+2];

    if ( x != node_coords(node_id,0) ||
         y != node_coords(node_id,1) ||
         z != node_coords(node_id,2) ) {
      printf("TestFEMesh::VerifyUnpack failed at %d : node %d : { %d %d %d } != { %d %d %d }\n",
             (int)i,(int)node_id, x,y,z,
             node_coords(node_id,0),
             node_coords(node_id,1),
             node_coords(node_id,2));
      ++update ;
    }
  }

  static inline
  size_type unpack( const array_type  & arg_node_coords ,
                    const size_type     arg_node_begin ,
                    const size_type     arg_node_count ,
                    const buffer_type & arg_buffer )
  {
    VerifyUnpack op ;
    op.node_coords = arg_node_coords ;
    op.buffer      = arg_buffer ;
    op.node_begin  = arg_node_begin ;
    return KokkosArray::parallel_reduce( arg_node_count , op );
  }
};

//----------------------------------------------------------------------------
}

