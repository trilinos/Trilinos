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

#ifndef KOKKOSARRAY_HOSTSPACE_HPP
#define KOKKOSARRAY_HOSTSPACE_HPP

#include <iosfwd>
#include <typeinfo>
#include <string>

#include <KokkosArray_Macros.hpp>

#if defined( __INTEL_COMPILER ) && 0
  /*
   *  Optimization level -O3 generates a seg-fault
   *  with this parameter in the hybrid parallel
   *  nonlinear use case boundary condition
   *  residual enforcement function.
   */

#define KOKKOSARRAY_ASSUME_ALIGNED( M , P ) __assume_aligned( P , M :: MEMORY_ALIGNMENT )

#else

#define KOKKOSARRAY_ASSUME_ALIGNED( M , P ) /* */

#endif

/*--------------------------------------------------------------------------*/

namespace KokkosArray {

/** \brief  Memory management on the host for devices */

class HostSpace {
public:

  enum { MEMORY_ALIGNMENT = 64 };
  enum { WORK_ALIGNMENT   =  8 };

  typedef HostSpace  memory_space ;

#if defined( __INTEL_COMPILER )
  typedef int        size_type ;
#else
  typedef size_t     size_type ;
#endif

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

  /*--------------------------------*/

  /** \brief  Print all tracked memory to the output stream. */
  static void print_memory_view( std::ostream & );

  /** \brief  Retrieve label associated with the input pointer */
  static std::string query_label( const void * );

  /*--------------------------------*/

  static 
  size_t preferred_alignment( size_t scalar_size , size_t scalar_count );

  /*--------------------------------*/
  /* Functions unique to the Host memory space */

  /** \brief  Assert called from the original, master thread.  */
  static void assert_master_thread( const char * const );
};

//----------------------------------------------------------------------------

template< class DstSpace , class SrcSpace >
struct DeepCopy ;

template<>
struct DeepCopy<HostSpace,HostSpace> {
  DeepCopy( void * dst , const void * src , size_t );
};

//----------------------------------------------------------------------------

template< class ExecutionSpace , class DataSpace >
struct VerifyExecutionSpaceCanAccessDataSpace ;

template<>
struct VerifyExecutionSpaceCanAccessDataSpace< HostSpace , HostSpace >
{
  inline static void verify(void) {}
  inline static void verify(const void *) {}
};

//----------------------------------------------------------------------------

} // namespace KokkosArray

#endif /* #define KOKKOSARRAY_HOSTSPACE_HPP */

