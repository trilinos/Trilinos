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
#include <OpenMP/Kokkos_OpenMPexec.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class WorkSpec >
class ParallelFor< FunctorType , WorkSpec , ::Kokkos::OpenMP >
{
public:

  inline
  ParallelFor( const FunctorType & functor , const size_t work_count )
  {
    OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_for");
    OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_for");

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread( omp_get_thread_num() );

      const std::pair< size_t , size_t > range = exec.work_range( work_count );

      for ( size_t iwork = range.first ; iwork < range.second ; ++iwork ) {
        functor( iwork );
      }
    }
/* END #pragma omp parallel */
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class WorkSpec >
class ParallelReduce< FunctorType , WorkSpec , Kokkos::OpenMP >
{
public:
  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;

  inline
  ParallelReduce( const FunctorType & functor ,
                  const size_t        work_count ,
                  pointer_type        result = 0 )
  {
    OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_reduce");
    OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_reduce");

    OpenMPexec::resize_reduce_scratch( Reduce::value_size( functor ) );

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread( omp_get_thread_num() );

      const std::pair<size_t,size_t> range = exec.work_range( work_count );

      typename Reduce::reference_type update = exec.reduce_reference( functor );

      functor.init( update );

      for ( size_t iw = range.first ; iw < range.second ; ++iw ) {
        functor( iw , update );
      }
    }
/* END #pragma omp parallel */

    {
      const int n = omp_get_max_threads();
      OpenMPexec & root = * OpenMPexec::get_thread(0);

      for ( int i = 1 ; i < n ; ++i ) {
        functor.join( root.reduce_reference( functor ) ,
                      OpenMPexec::get_thread(i)->reduce_reference( functor ) );
      }

      Reduce::final( functor , root.reduce_pointer( functor ) );

      if ( result ) {
        const pointer_type ptr = root.reduce_pointer( functor );
        const int n = Reduce::value_count( functor );

        for ( int i = 0 ; i < n ; ++i ) { result[i] = ptr[i] ; }
      }
    }
  }

  void wait() {}
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_USE_PRAGMA_SIMD )

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class ParallelReduce< FunctorType , VectorParallel , ::Kokkos::OpenMP >
{
public:
  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;

  inline
  ParallelReduce( const FunctorType & functor ,
                  const size_t        work_count ,
                  pointer_type        result = 0 )
  {
    typedef integral_constant< size_t , OpenMPexec::VECTOR_LENGTH >     vector_length ;
    typedef integral_constant< size_t , OpenMPexec::VECTOR_LENGTH - 1 > vector_mask ;

    OpenMPexec::verify_is_process("Kokkos::OpenMP parallel_reduce");
    OpenMPexec::verify_initialized("Kokkos::OpenMP parallel_reduce");

    OpenMPexec::resize_reduce_scratch( Reduce::value_size( functor ) * vector_length::value );

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread( omp_get_thread_num() );

      const std::pair<size_t,size_t> range = exec.work_range( work_count );

#pragma simd
#pragma ivdep
      for ( size_t iv = 0 ; iv < vector_length::value ; ++iv ) {
        functor.init( exec.reduce_reference( functor , iv ) );
      }

#pragma simd vectorlength( vector_length::value )
#pragma ivdep
      for ( size_t iw = range.first ; iw < range.second ; ++iw ) {
        functor( iw , exec.reduce_reference( functor , iw & vector_mask::value ) );
      }

      for ( size_t iv = 1 ; iv < vector_length::value ; ++iv ) {
        functor.join( exec.reduce_reference( functor , 0 ) ,
                      exec.reduce_reference( functor , iv ) );
      }
    }
/* END #pragma omp parallel */

    {
      const int n = omp_get_max_threads();
      OpenMPexec & root = * OpenMPexec::get_thread(0);

      for ( int i = 1 ; i < n ; ++i ) {
        functor.join( root.reduce_reference( functor ) ,
                      OpenMPexec::get_thread(i)->reduce_reference( functor ) );
      }

      Reduce::final( functor , root.reduce_pointer( functor ) );

      if ( result ) {
        const pointer_type ptr = root.reduce_pointer( functor );
        const int n = Reduce::value_count( functor );

        for ( int i = 0 ; i < n ; ++i ) { result[i] = ptr[i] ; }
      }
    }
  }

  void wait() {}
};

} // namespace Impl
} // namespace Kokkos

#endif /* #if defined( KOKKOS_USE_PRAGMA_SIMD ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class ParallelReduce< FunctorType , ParallelWorkRequest , ::Kokkos::OpenMP >
{
public:
  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;

  inline
  ParallelReduce( const FunctorType         & functor ,
                  const ParallelWorkRequest & work ,
                  pointer_type                result = 0 )
  {
    OpenMPexec::assert_ready("Kokkos::OpenMP - parallel_reduce");

    OpenMPexec::resize_reduce_scratch( Reduce::value_size( functor ) );

#pragma omp parallel
    {
      OpenMPexec & exec = * OpenMPexec::get_thread( omp_get_thread_num() );

      typename Reduce::reference_type update = exec.reduce_reference( functor );

      functor.init( update );

      for ( exec.team_work_init( work.league_size ) ; exec.team_work_avail() ; exec.team_work_next() ) {
        functor( OpenMP( exec ) , update );
      }
    }
/* END #pragma omp parallel */

    {
      const int n = omp_get_max_threads();
      OpenMPexec & root = * OpenMPexec::get_thread(0);

      for ( int i = 1 ; i < n ; ++i ) {
        functor.join( root.reduce_reference( functor ) ,
                      OpenMPexec::get_thread(i)->reduce_reference( functor ) );
      }

      Reduce::final( functor , root.reduce_pointer( functor ) );

      if ( result ) {
        const pointer_type ptr = root.reduce_pointer( functor );
        const int n = Reduce::value_count( functor );

        for ( int i = 0 ; i < n ; ++i ) { result[i] = ptr[i] ; }
      }
    }
  }

  void wait() {}
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_OPENMP_PARALLEL_HPP */

