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

#if defined( KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION ) || \
    defined( KOKKOSARRAY_MACRO_DEVICE ) || \
    defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION ) || \
    defined( KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION )

#error "Including <KokkosArray_Cuda_macros.hpp> with macros already defined"

#elif ! defined( __CUDACC__ )

/*  If not compiling with Cuda compiler
 *  then pure-device functions are not available
 */

#define KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION /* */
#define KOKKOSARRAY_MACRO_DEVICE                      KokkosArray::Cuda
/* #define KOKKOSARRAY_MACRO_DEVICE_FUNCTION */
#define KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION      /* */
#define KOKKOSARRAY_MACRO_DEVICE_CAN_THROW( expr ) expr
#define KOKKOSARRAY_MACRO_CHECK( expr )  expr

#else

/*  If compiling with Cuda compiler
 *  then must clarify what functions are only available on the device
 *  versus available on both the device and host.
 */

#define KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION /* */
#define KOKKOSARRAY_MACRO_DEVICE                      KokkosArray::Cuda
#define KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION    __device__ __host__
#define KOKKOSARRAY_MACRO_DEVICE_FUNCTION             __device__

/* Device-only functions are not compiled unless compiling for the device. */
#if defined( __CUDA_ARCH__ )
#define KOKKOSARRAY_MACRO_DEVICE_CAN_THROW( expr )  /* */
#define KOKKOSARRAY_MACRO_CHECK( expr )             /* */
#elif ! defined( KOKKOSARRAY_BOUNDS_CHECK )
#define KOKKOSARRAY_MACRO_DEVICE_CAN_THROW( expr )  expr
#define KOKKOSARRAY_MACRO_CHECK( expr )             /* */
#else
#define KOKKOSARRAY_MACRO_DEVICE_CAN_THROW( expr )  expr
#define KOKKOSARRAY_MACRO_CHECK( expr )             expr
#endif

#endif

