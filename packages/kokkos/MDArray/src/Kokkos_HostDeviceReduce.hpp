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

#ifndef KOKKOS_HOSTDEVICEREDUCE_HPP
#define KOKKOS_HOSTDEVICEREDUCE_HPP

#include <algorithm>
#include <TPI.h>

namespace Kokkos {

class HostDevice ;

template< class FunctorType , class DeviceType > struct ParallelReduce ;

template< class FunctorType >
inline
void parallel_reduce( const FunctorType & functor ,
                      const typename FunctorType::reduce_type & result )
{
  ParallelReduce< FunctorType , typename FunctorType::device_type >
    ::run( functor , result );
}

//----------------------------------------------------------------------------

template< class FunctorType >
struct ParallelReduce< FunctorType , HostDevice > {

  typedef          HostDevice               device_type ;
  typedef typename HostDevice::size_type    size_type ;
  typedef typename FunctorType::reduce_type reduce_type ;
  typedef typename reduce_type::value_type  reduce_value_type ;

  FunctorType m_functor ;
  reduce_type m_result ;

  ParallelReduce( const FunctorType & functor , const reduce_type & result )
    : m_functor( functor ), m_result( result ) {}

  static void run_functor_on_tpi( TPI_Work * work )
  {
    const ParallelReduce & self      = *((const ParallelReduce *) work->info );
    const FunctorType & functor      = self.m_functor ;
    const reduce_type & local_result = *((const reduce_type *) work->reduce );

    const int work_count = functor.work_count();
    const int work_inc   = ( work_count + work->count - 1 ) / work->count ;
    const int work_begin = work_inc * work->rank ;
    const int work_end   = std::max( work_begin + work_inc , work_count );

    for ( int iwork = work_begin ; iwork < work_end ; ++iwork ) {
      functor( iwork , local_result );
    }
  }

  static void run_join_on_tpi( TPI_Work * work , const void * reduce )
  {
    const reduce_type & local_result = *((const reduce_type *) work->reduce );
    const reduce_type & local_source = *((const reduce_type *) reduce );

    FunctorType::join( local_result , local_source );
  }

  static void run_init_on_tpi( TPI_Work * work )
  {
    const ParallelReduce & self      = *((const ParallelReduce *) work->info );
          reduce_type & local_result = *((      reduce_type *) work->reduce );

    reduce_value_type * const value =
      (reduce_value_type*) ( ( & local_result ) + 1 );

    local_result.assign_on_device( self.m_result ); // Copy the shape
    local_result.assign_on_device( value ); // Set the pointer

    FunctorType::init( local_result ); // Initialize the values
  }

  static void run( const FunctorType & functor , const reduce_type & result )
  {
    size_type reduce_count = 1 ;

    for ( size_type i = 0 ; i < result.rank() ; ++i ) {
      reduce_count *= result.dimension(i);
    }

    size_type reduce_total_size = sizeof(reduce_type) +
                                  sizeof(reduce_value_type) * reduce_count ;

    ParallelReduce tmp( functor , result ); 

    TPI_Run_threads_reduce( & run_functor_on_tpi , & tmp ,
                            & run_join_on_tpi , & run_init_on_tpi ,
                            reduce_total_size , & tmp.m_result );
  }
};

}

#endif /* KOKKOS_HOSTDEVICEREDUCE_HPP */

