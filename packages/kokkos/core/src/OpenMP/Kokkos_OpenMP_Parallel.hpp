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

#ifndef KOKKOS_OPENMP_PARALLEL_HPP
#define KOKKOS_OPENMP_PARALLEL_HPP

#include <omp.h>

#include <Kokkos_Parallel.hpp>
#include <Kokkos_ParallelReduce.hpp>
#include <Host/Kokkos_Host_Thread.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class FunctorType , class WorkSpec >
class ParallelFor< FunctorType , WorkSpec , OpenMP > {
public:

  typedef Kokkos::HostSpace::size_type  size_type ;

  inline
  ParallelFor( const FunctorType & functor ,
               const size_type work_count )
  {
#if defined( __INTEL_COMPILER )
    enum { vectorize = is_same<WorkSpec,VectorParallel>::value };
#else
    enum { vectorize = 0 };
#endif

    OpenMP::assert_ready("Kokkos::OpenMP - parallel_for");

#pragma omp parallel
    {
      HostThread & thread = * OpenMP::get_host_thread();

      const std::pair< size_type , size_type > range =
        thread.work_range( work_count );

      if ( ! vectorize ) {
        for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
          functor( iwork );
        }
      }
      else {
#if defined( __INTEL_COMPILER )
#pragma simd
#pragma ivdep
        for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
          functor( iwork );
        }
#endif
      }
    }
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType >
class ParallelReduce< FunctorType , size_t , OpenMP > {
public:

  typedef Kokkos::HostSpace::size_type  size_type ;

  typedef typename ReduceAdapter< FunctorType >::pointer_type pointer_type ;

  typedef ReduceAdapter< FunctorType > ReduceOp ;

  const ReduceOp reduce ;

  inline
  ParallelReduce( const FunctorType  & functor ,
                  const size_t         work_count ,
                  pointer_type         result )
    : reduce( functor )
  {
#if defined( __INTEL_COMPILER )
    enum { work_align = is_same<WorkSpec,VectorParallel>::value &&
                        power_of_two<HostSpace::WORK_ALIGNMENT>::value 
                      ? HostSpace::WORK_ALIGNMENT : 1 };
    enum { work_mask  = work_align - 1 };
#else
    enum { work_align = 1 };
    enum { work_mask  = 0 };
#endif

    OpenMP::assert_ready("Kokkos::OpenMP - parallel_reduce");

    OpenMP::resize_reduce_scratch( reduce.value_size() * work_align );

#pragma omp parallel
    {
      HostThread & thread = * OpenMP::get_host_thread();

      const std::pair< size_type , size_type > range =
        thread.work_range( work_count );

      if ( ! work_mask ) {
        reduce.init( thread.reduce_data() );

        for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
          functor( iwork , reduce.reference( thread.reduce_data() ) );
        }
      }
      else {
#if defined( __INTEL_COMPILER )

#pragma simd
#pragma ivdep
        for ( size_type j = 0 ; j < work_align ; ++j ) {
          reduce.init( thread.reduce_data() , j );
        }

#pragma simd vectorlength(work_align)
#pragma ivdep
        for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
          functor( iwork , reduce.reference( thread.reduce_data() , iwork & work_mask ) );
        }

        reduce.template join< work_align >( thread.reduce_data() );
#endif
      }
    }

    {
      const unsigned n = omp_get_max_threads();
      HostThread & root = * OpenMP::get_host_thread();
      for ( unsigned i = 1 ; i < n ; ++i ) {
        HostThread & th = * OpenMP::get_host_thread(i);
        reduce.join( root.reduce_data() , th.reduce_data() );
      }
      // reduce.finalize( root.reduce_data() );
    }

    {
      pointer_type ptr = (pointer_type) OpenMP::root_reduce_scratch();

      for ( unsigned i = 0 ; i < reduce.value_count(); ++i ) {
        result[i] = ptr[i] ;
      }
    }
  }

  void wait() {}
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_OPENMP_PARALLEL_HPP */

