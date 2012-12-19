/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOSARRAY_CUDA_TEXTURE
#define KOKKOSARRAY_CUDA_TEXTURE

#include <stdlib.h>

#include <KokkosArray_Macros.hpp>
#include <KokkosArray_View.hpp>

namespace KokkosArray {
namespace Impl {

/** \brief
 *
 *  Usage
 *  1) Declare at global file scope:  CudaTexture<MyType>::type  my_texture_obj ;
 *  2) Bind in host code:             CudaTexture<MyType>::bind( my_texture_obj , my_ptr );
 *  3) Access from Cuda code:         CudaTexture<MyType>::fetch( my_texture_obj , i );
 */

template< typename T >
struct CudaTexture ;

} // namespace Impl
} // namespace KokkosArray

#if defined( __CUDACC__ )

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template<>
struct CudaTexture<double>
{
  typedef texture<int2,cudaTextureType1D,cudaReadModeElementType> texture_type ;

  static inline
  void bind( const texture_type & tex_ref , const double * ptr , const size_t count )
  {
    const cudaChannelFormatDesc desc = cudaCreateChannelDesc<int2>();

    cudaBindTexture( NULL , & tex_ref , (void*) ptr , & desc , sizeof(double) * count );
  }

  static inline
  void unbind( const texture_type & tex_ref )
  { cudaUnbindTexture( & tex_ref ); }

  static inline __device__
  double fetch( texture_type tex_ref , int i )
  {
    const int2 v = tex1Dfetch(tex_ref,i);
    return __hiloint2double(v.y,v.x);
  }
};

//----------------------------------------------------------------------------

template< typename T >
struct CudaTexture
{
  typedef texture<T,cudaTextureType1D,cudaReadModeElementType> texture_type ;

  static inline
  void bind( const texture_type & tex_ref , const T * ptr , const size_t count )
  {
    const cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>();

    cudaBindTexture( NULL , & tex_ref , (void*) ptr , & desc , sizeof(T) * count );
  }

  static inline
  void unbind( const texture_type & tex_ref )
  { cudaUnbindTexture( & tex_ref ); }

  static inline __device__
  T fetch( texture_type tex_ref , int i )
  { return tex1Dfetch(tex_ref,i); }
};

} // namespace Impl
} // namespace KokkosArray

#endif /* #if defined( __CUDACC__ ) */

#endif /* #ifndef KOKKOSARRAY_CUDA_TEXTURE */

