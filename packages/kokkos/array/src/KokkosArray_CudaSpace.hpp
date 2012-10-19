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

#ifndef KOKKOSARRAY_CUDASPACE_HPP
#define KOKKOSARRAY_CUDASPACE_HPP

#include <iosfwd>
#include <typeinfo>
#include <string>

#include <KokkosArray_HostSpace.hpp>

/*--------------------------------------------------------------------------*/


namespace KokkosArray {

/** \brief  Cuda memory management */

class CudaSpace {
public:

  typedef CudaSpace     memory_space ;
  typedef unsigned int  size_type ;

  /** \brief  Allocate a contiguous block of memory on the Cuda device
   *          with size = scalar_size * scalar_count.
   *
   *  The input label is associated with the block of memory.
   *  The block of memory is tracked via reference counting where
   *  allocation gives it a reference count of one.
   *
   *  Allocation may only occur on the master thread of the process.
   */
  static void * allocate( const std::string    & label ,
                          const std::type_info & scalar_type ,
                          const size_t           scalar_size ,
                          const size_t           scalar_count );

#if ! defined( __CUDA_ARCH__ )
  /** \brief  Increment the reference count of the block of memory
   *          in which the input pointer resides.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void increment( const void * );

  /** \brief  Decrement the reference count of the block of memory
   *          in which the input pointer resides.  If the reference
   *          count falls to zero the memory is deallocated.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void decrement( const void * );
#else
  static __device__ void increment( const void * ) {}
  static __device__ void decrement( const void * ) {}
#endif

  /** \brief  Print all tracked memory to the output stream. */
  static void print_memory_view( std::ostream & );

  /** \brief  Retrieve label associated with the input pointer */
  static std::string query_label( const void * );

  /*--------------------------------*/

  /** \brief  Query the preferred value of 'scalar_count' which
   *          would given the best performing alignement for
   *          memory accesses.
   */
  static 
  size_t preferred_alignment( size_t scalar_size , size_t scalar_count );

  /*--------------------------------*/

  static void access_error();
  static void access_error( const void * const );
};

//----------------------------------------------------------------------------

template<>
struct DeepCopy<HostSpace,CudaSpace> {
  DeepCopy( void * dst , const void * src , size_t );
};

template<>
struct DeepCopy<CudaSpace,HostSpace> {
  DeepCopy( void * dst , const void * src , size_t );
};

template<>
struct DeepCopy<CudaSpace,CudaSpace> {
  DeepCopy( void * dst , const void * src , size_t );
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Compiler specifics macros

#if defined( __CUDACC__ )

#if defined( KOKKOSARRAY_INLINE_FUNCTION )
#undef KOKKOSARRAY_INLINE_FUNCTION
#endif
#define KOKKOSARRAY_INLINE_FUNCTION __device__ __host__ inline


#if defined( __CUDA_ARCH__ ) /* Cuda compiling code to run on the Cuda device */

#if ( __CUDA_ARCH__ < 200 )
#error "Cuda device capability >= 2.0 is required"
#endif

#if defined( KOKKOSARRAY_EXECUTION_SPACE )
#undef KOKKOSARRAY_EXECUTION_SPACE
#endif
#define KOKKOSARRAY_EXECUTION_SPACE CudaSpace

#endif  /* #if defined( __CUDA_ARCH__ ) */

/* Force asserts to trap memory space access errors. */
#if defined( NDEBUG )
#undef NDEBUG
#include <assert.h>
#define NDEBUG
#else
#include <assert.h>
#endif

#endif /* #if defined( __CUDACC__ ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

/** \brief  Cuda accessing Cuda is good */
template<>
struct VerifyExecutionSpaceCanAccessDataSpace< CudaSpace , CudaSpace >
{
  KOKKOSARRAY_INLINE_FUNCTION static void verify( void ) {}
  KOKKOSARRAY_INLINE_FUNCTION static void verify( const void * ) {}
};

/** \brief  Cuda accessing non-Cuda is bad */
template< class DataSpace >
struct VerifyExecutionSpaceCanAccessDataSpace< CudaSpace , DataSpace >
{
  KOKKOSARRAY_INLINE_FUNCTION static void verify(void)
  {
    const int Cuda_execution_on_Cuda_memory = 0 ;
    assert(Cuda_execution_on_Cuda_memory);
  }

  KOKKOSARRAY_INLINE_FUNCTION static void verify( const void * )
  {
    const int Cuda_execution_on_Cuda_memory = 0 ;
    assert(Cuda_execution_on_Cuda_memory);
  }
};

/** \brief  Produce error message when trying to access Cuda 
 *          memory on the host.
 */
template< class ExecutionSpace >
struct VerifyExecutionSpaceCanAccessDataSpace< ExecutionSpace , CudaSpace >
{
  inline static void verify( void ) { CudaSpace::access_error(); }
  inline static void verify( const void * p ) { CudaSpace::access_error(p); }
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOSARRAY_CUDASPACE_HPP */


