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

#ifndef KOKKOS_MDARRAYUTIL_HPP
#define KOKKOS_MDARRAYUTIL_HPP

#include <Kokkos_HostMDArrayView.hpp>
#include <Kokkos_CudaMDArrayView.hpp>

namespace Kokkos {

template < typename Scalar , class DeviceMapDest , class DeviceMapSrc >
struct DeepCopy ;

template < typename Scalar , class DeviceMapDest , class DeviceMapSrc >
void deep_copy( const MDArray<Scalar,DeviceMapDest> & destination ,
                const MDArray<Scalar,DeviceMapSrc> & source )
{
  DeepCopy<Scalar,DeviceA,DeviceB>::run( destination , source );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// Deep copy from Cuda to Host
// A class for partial specialization.
template< typename Scalar >
struct DeepCopy< Scalar , HostMap , CudaMap > {

  typedef size_t size_type ;

  static void run( const MDArrayView<Scalar,HostMap> & dst ,
                   const MDArrayView<Scalar,CudaMap> & src )
  {
    bool ok_shape = dst.rank() == src.rank();
    size_type count = 1 ;
    for ( size_type i = 0 ; ok_shape && i < dst.rank() ; ++i ) {
      ok_shape = dst.dimension(i) == src.dimension(i);
      count *= dst.dimension(i);
    }

    if ( ! ok_shape ) { throw_deep_copy_shape_error( dst , src ); }

    const size_type nbytes = sizeof(Scalar) * count ;

    Scalar * const src_ptr_cuda = src.address_on_device();
    Scalar * const dst_ptr_host = dst.address_on_device();
    Scalar * const tmp_ptr_host = (Scalar *) malloc( nbytes );

    // Copy into temporary buffer and then do permutation between layouts.

    cudaMemcpy( tmp_ptr_host , src_ptr_cuda , nbytes , CudaMemcpyDeviceToHost );

    // Now have to permute between the two maps,
    // requires knowledge of the maps' implementation.
    //
    // HostMapped array is left-to-right strided
    // CudaMapped array has last dimension fasted striding,
    // and then left-to-right for remaining dimensions

    const size_type n_parallel = dst.dimension( dst.rank() - 1 );
    const size_type n_base = count / n_parallel ;

    for ( size_type iP = 0 ; iP < n_parallel ; ++iP ) {
      for ( size_type i = 0 ; i < n_base ; ++i ) {
        dst_ptr_host[ i + iP * n_base ] = tmp_ptr_host[ iP + i * n_parallel ];
      }
    }

    free( tmp_ptr_host );
  }
};

// Deep copy from Host to Cuda
template< typename Scalar >
struct DeepCopy< Scalar , CudaMap , HostMap > {

  static void run( const MDArrayView<Scalar,CudaMap> & dst ,
                   const MDArrayView<Scalar,HostMap> & src )
  {
    bool ok_shape = dst.rank() == src.rank();
    size_type count = 1 ;
    for ( size_type i = 0 ; ok_shape && i < dst.rank() ; ++i ) {
      ok_shape = dst.dimension(i) == src.dimension(i);
      count *= dst.dimension(i);
    }

    if ( ! ok_shape ) { throw_deep_copy_shape_error( dst , src ); }

    const size_type nbytes = sizeof(Scalar) * count ;

    Scalar * const dst_ptr_cuda = dst.address_on_device();
    Scalar * const src_ptr_host = src.address_on_device();
    Scalar * const tmp_ptr_host = (Scalar *) malloc( nbytes );

    // Permute between the two maps,
    // requires knowledge of the maps' implementation.
    //
    // HostMapped array is left-to-right strided
    // CudaMapped array has last dimension fasted striding,
    // and then left-to-right for remaining dimensions

    const size_type n_parallel = dst.dimension( dst.rank() - 1 );
    const size_type n_base = count / n_parallel ;

    for ( size_type i = 0 ; i < n_base ; ++i ) {
      for ( size_type iP = 0 ; iP < n_parallel ; ++iP ) {
        tmp_ptr_host[ iP + i * n_parallel ] = src_ptr_host[ i + iP * n_base ] ;
      }
    }
    // Copy into temporary buffer and then do permutation between layouts.

    cudaMemcpy( dst_ptr_cuda , tmp_ptr_host , nbytes , CudaMemcpyHostToDevice );

    free( tmp_ptr_host );
  }
};


} // namespace Kokkos

#endif /* KOKKOS_MDARRAYUTIL_HPP */


