/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <GramSchmidt.hpp>

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_Timer.hpp>

#include <LinAlgBLAS.hpp>

#include <iostream>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Test {


template< typename ScalarQ ,
          typename ScalarR ,
          class DeviceType , class Management >
void modified_gram_schmidt(

  const Kokkos::View< ScalarQ ** ,
                           Kokkos::LayoutLeft ,
                           DeviceType ,
                           Management > & Q ,

  const Kokkos::View< ScalarR ** ,
                           Kokkos::LayoutLeft ,
                           DeviceType ,
                           Management > & R ,

  comm::Machine machine )
{
  const Kokkos::ALL ALL ;

  typedef Kokkos::View< ScalarQ * ,
                             Kokkos::LayoutLeft ,
                             DeviceType ,
                             Kokkos::MemoryUnmanaged >
    vector_view_type ;

  const typename
    Kokkos::View< ScalarR** ,
                       Kokkos::LayoutLeft ,
                       DeviceType >::
    HostMirror hostR = Kokkos::create_mirror_view( R );

  const int length = Q.dimension_0();
  const int count  = Q.dimension_1();

  for ( int j = 0 ; j < count ; ++j ) {

    const vector_view_type  Qj = Kokkos::subview( Q , ALL , j );

    // reads  += length
    // writes += 0
    // flops  += 1 + 2 * length
    const double norm_Qj = Kokkos::norm2( length , Qj , machine );

    hostR(j,j) = norm_Qj ;

    // reads  += length
    // writes += length
    // flops  += 1 + length
    Kokkos::scale( length , 1.0 / norm_Qj , Qj );

    for ( int k = j + 1 ; k < count ; ++k ) {

      const vector_view_type  Qk = Kokkos::subview( Q , ALL , k );

      // reads  += 2 * length
      // writes += 0
      // flops  += 2 * length
      const double Qj_dot_Qk =
        Kokkos::dot( length , Qj , Qk , machine );

      hostR(j,k) = Qj_dot_Qk ;

      // reads  += 2 * length
      // writes += length
      // flops += 2 * length
      Kokkos::axpy( length , - Qj_dot_Qk , Qj , Qk );
    }
  }

  // reads  += 0
  // writes += count * count
  Kokkos::deep_copy( R , hostR );
}

struct ModifedGramSchmidCounts {
  size_t flops ;
  size_t reads ;
  size_t writes ;

  ModifedGramSchmidCounts( const unsigned length ,
                           const unsigned count )
  : flops(0) , reads(0), writes(0)
  {
    for ( unsigned j = 0 ; j < count ; ++j ) {
      // norm2:
      reads  += length ;
      flops  += 1 + 2 * length ;
      // scale:
      reads  += length ;
      writes += length ;
      flops  += 1 + length ;
      for ( unsigned k = j+1 ; k < count ; ++k ) {
        // dot
        reads  += 2 * length ;
        flops  += 2 * length ;
        // axpy
        reads  += 2 * length ;
        writes += length ;
        flops  += 2 * length ;
      }
    }
  }
};

//----------------------------------------------------------------------------

template<>
void driver_modified_gram_schmidt
<
#if defined( __CUDACC__ )
Kokkos::Cuda
#else
Kokkos::Threads
#endif
>
  ( const int length_begin ,
    const int length_end ,
    const int count ,
    const int iter ,
    comm::Machine machine )
{
#if defined( __CUDACC__ )
  typedef Kokkos::Cuda Device ;
#else
  typedef Kokkos::Threads Device ;
#endif

  const int comm_size = comm::size( machine );
  const int comm_rank = comm::rank( machine );

  if ( comm_rank == 0 ) {

    std::cout << ( Kokkos::Impl::is_same<Device,Kokkos::Cuda>::value ?
                   "\"Cuda\"" : "\"Threads\"" )
              << " , \"Double Precision\""
              << std::endl ;

    std::cout << "\"Length\" , \"Count\" , \"millisec\""
              << " , \"Gflops\" , \"Greads\" , \"Gwrites\""
              << " ,   \"Gflops/s\" , \"Read GB/s\" , \"Write GB/s\""
              << std::endl ;
  }

  for ( int length = length_begin ; length < length_end ; length *= 2 ) {

    const ModifedGramSchmidCounts counts( length , count );

    const int local_length_upper = ( length + comm_size - 1 ) / comm_size ;
    const int local_begin  = std::min( length , local_length_upper * comm_rank );
    const int local_next   = std::min( length , local_length_upper * ( comm_rank + 1 ) );
    const int local_length = local_next - local_begin ;

    typedef Kokkos::View< double ** ,
                               Kokkos::LayoutLeft ,
                               Device > matrix_double_type ;

    const matrix_double_type Q( "Q" , local_length , count );
    const matrix_double_type R( "R" , count , count );

    const matrix_double_type::HostMirror hQ =
      Kokkos::create_mirror_view( Q );

    for ( int j = 0 ; j < count ; ++j ) {
      for ( int i = 0 ; i < local_length ; ++i ) {
        hQ(i,j) = ( i + 1 ) * ( j + 1 );
      }
    }

    double dt_min = 0 ;

    for ( int j = 0 ; j < iter ; ++j ) {
      Kokkos::deep_copy( Q , hQ );

      Kokkos::Impl::Timer timer ;

      modified_gram_schmidt( Q , R , machine );

      const double dt = comm::max( machine , timer.seconds() );

      if ( 0 == j || dt < dt_min ) dt_min = dt ;
    }

    if ( 0 == comm_rank ) {

      const double milli_sec   = dt_min * 1.0e3 ;
      const double giga_flops  = ( counts.flops / dt_min ) / 1.0e9 ;
      const double GB_reads  = ( counts.reads * sizeof(double) ) / ( dt_min * 1.0e9 );
      const double GB_writes = ( counts.writes * sizeof(double) ) / ( dt_min * 1.0e9 );

      std::cout << length     << " , "
                << count      << " , "
                << milli_sec  << " , "
                << double(counts.flops) / 1.0e9 << " , "
                << double(counts.reads) / 1.0e9 << " , "
                << double(counts.writes) / 1.0e9 << " ,   "
                << giga_flops << " , "
                << GB_reads << " , "
                << GB_writes
                << std::endl ;
    }
  }
}

}

