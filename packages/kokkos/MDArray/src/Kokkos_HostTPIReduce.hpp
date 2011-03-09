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

#ifndef KOKKOS_HOSTTPIREDUCE_HPP
#define KOKKOS_HOSTTPIREDUCE_HPP

#include <algorithm>
#include <TPI.h>
#include <Kokkos_HostTPI.hpp>
#include <Kokkos_HostTPIValueView.hpp>
#include <Kokkos_ParallelReduce.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------

namespace {

template< class FunctorType >
struct TPI_ParallelReduce {

  typedef HostTPI::size_type               size_type ;
  typedef typename FunctorType::value_type value_type ;

  const FunctorType m_functor ;
  const size_type   m_work_count ;

  TPI_ParallelReduce( const FunctorType & functor ,
                      const size_type     work_count )
    : m_functor( functor ), m_work_count( work_count ) {}

  static void run( TPI_Work * work )
  {
    const TPI_ParallelReduce & self = *((const TPI_ParallelReduce *) work->info );
    value_type & update = *((value_type *) work->reduce );

    // Number of threads:   work->count
    // Rank of this thread: work->rank
    const size_type work_inc   = ( self.m_work_count + work->count - 1 ) / work->count ;
    const size_type work_begin = work_inc * work->rank ;
    const size_type work_end   = std::min( work_begin + work_inc , self.m_work_count );

    for ( size_type iwork = work_begin ; iwork < work_end ; ++iwork ) {
      self.m_functor( iwork , update );
    }
  }

  static void join( TPI_Work * work , const void * reduce )
  {
          value_type & update = *((      value_type *) work->reduce );
    const value_type & source = *((const value_type *) reduce );

    FunctorType::join( update , source );
  }

  static void init( TPI_Work * work )
  {
    value_type & update = *((value_type *) work->reduce );

    FunctorType::init( update );
  }
};

} // namespace

//----------------------------------------------------------------------------

template< class FunctorType >
struct ParallelReduce< FunctorType , void , HostTPI >
{
  typedef HostTPI::size_type                size_type ;
  typedef typename FunctorType::value_type  value_type ;
  typedef TPI_ParallelReduce< FunctorType > tpi_functor ;

  static value_type run( const size_type     work_count ,
                         const FunctorType & functor )
  {
    value_type value ;

    FunctorType::init( value );

    tpi_functor tmp( functor , work_count );

    TPI_Run_threads_reduce( & tpi_functor::run , & tmp ,
                            & tpi_functor::join ,
                            & tpi_functor::init ,
                            sizeof(value_type) , & value );

    return value ;
  }
};

//----------------------------------------------------------------------------

template< class FunctorType >
struct ParallelReduce<
  FunctorType ,
  ValueView< typename FunctorType::value_type , HostTPI > ,
  HostTPI >
{
  typedef HostTPI::size_type                size_type ;
  typedef typename FunctorType::value_type  value_type ;
  typedef ValueView< value_type , HostTPI > view_type ;
  typedef TPI_ParallelReduce< FunctorType > tpi_functor ;

  static void run( const size_type     work_count ,
                   const FunctorType & functor ,
                   const view_type   & view )
  {
    FunctorType::init( *view );

    tpi_functor tmp( functor , work_count );

    TPI_Run_threads_reduce( & tpi_functor::run , & tmp ,
                            & tpi_functor::join ,
                            & tpi_functor::init ,
                            sizeof(value_type) , view.address_on_device  );
  }
};

//----------------------------------------------------------------------------

template< class FunctorType , class FinalizeType >
struct ParallelReduce< FunctorType , FinalizeType , HostTPI > {

  typedef HostTPI::size_type                size_type ;
  typedef typename FunctorType::value_type  value_type ;
  typedef TPI_ParallelReduce< FunctorType > tpi_functor ;

  static void run( const size_type      work_count ,
                   const FunctorType  & functor ,
                   const FinalizeType & finalize )
  {
    value_type result ;

    FunctorType::init( result );

    tpi_functor tmp( functor , work_count );

    TPI_Run_threads_reduce( & tpi_functor::run , & tmp ,
                            & tpi_functor::join ,
                            & tpi_functor::init ,
                            sizeof(value_type) , & result );

    finalize( result );
  }
};

//----------------------------------------------------------------------------

}

#endif /* KOKKOS_HOSTTPIREDUCE_HPP */

