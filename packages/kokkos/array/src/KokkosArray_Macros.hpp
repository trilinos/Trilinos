/*
//@HEADER
// ************************************************************************
//
//                             KokkosArray
//         Manycore Performance-Portable Multidimensional Arrays
//
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_MACROS_HPP
#define KOKKOSARRAY_MACROS_HPP

//----------------------------------------------------------------------------

namespace KokkosArray {

class Host;
class Cuda;

} //namespace KokkosArray

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template <typename MemorySpace>
struct Restrict_Function_To_Device
{ static void fail() {} };

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------

#if defined( __CUDACC__ )

/*  Compiling with a CUDA compiler for host and device code.
 *  Declare functions available on host and device, or device only.
 */

#define KOKKOSARRAY_INLINE_FUNCTION        inline __device__ __host__
#define KOKKOSARRAY_INLINE_DEVICE_FUNCTION inline __device__

#if defined( __CUDA_ARCH__ )

/*  Compiling with CUDA compiler for device code.
 *  The value of __CUDA_ARCH__ is an integer value XY0
 *  identifying compilation for 'compute_XY' architecture.
 */

#define KOKKOSARRAY_RESTRICT_CODE_TO_DEVICE

/*  If compiling for device code cannot do expression checking */

#if defined( KOKKOSARRAY_EXPRESSION_CHECK )
#undef  KOKKOSARRAY_EXPRESSION_CHECK
#endif

/*  Compile-time restriction that code is being compiled with Cuda compiler
 *  to run on the Cuda device.
 */
namespace KokkosArray { namespace Impl {

template<>
struct Restrict_Function_To_Device<KokkosArray::Cuda>
{
  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  static void success() {}
};

}} //namespace KokkosArray::Impl

#else

/*  Compiling with CUDA compiler for host code. */

#endif

//----------------------------------------------------------------------------

#elif defined( _OPENMP )

/*  Compiling with an OpenMP compiler.
 *  The value of _OPENMP is an integer value YYYYMM 
 *  where YYYY and MM are the year and month designation
 *  of the supported OpenMP API version.
 */

#define KOKKOSARRAY_INLINE_FUNCTION         inline
#define KOKKOSARRAY_INLINE_DEVICE_FUNCTION  inline

#define KOKKOSARRAY_RESTRICT_CODE_TO_DEVICE

//----------------------------------------------------------------------------

#else //on host

#define KOKKOSARRAY_INLINE_FUNCTION         inline
#define KOKKOSARRAY_INLINE_DEVICE_FUNCTION  inline

#define KOKKOSARRAY_RESTRICT_CODE_TO_DEVICE

namespace KokkosArray { namespace Impl {

template<>
struct Restrict_Function_To_Device<KokkosArray::Host>
{ static void success() {} };

}} //namespace KokkosArray::Impl

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifdef  KOKKOSARRAY_EXPRESSION_CHECK
#define KOKKOSARRAY_CHECK_EXPR(expr) do{expr;}while(false)
#else
#define KOKKOSARRAY_CHECK_EXPR(expr) do{/*expr;*/}while(false)
#endif

#endif // KOKKOSARRAY_MACROS_HPP

