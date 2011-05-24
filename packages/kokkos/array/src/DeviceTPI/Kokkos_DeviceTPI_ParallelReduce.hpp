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

#ifndef KOKKOS_DEVICETPI_PARALLELREDUCE_HPP
#define KOKKOS_DEVICETPI_PARALLELREDUCE_HPP

#include <algorithm>
#include <TPI.h>

namespace Kokkos {

template< class FunctorType >
class ParallelReduce< FunctorType , void , DeviceTPI > {
public:
  typedef DeviceTPI              device_type ;
  typedef device_type::size_type size_type ;
  typedef typename FunctorType::value_type value_type ;

  const FunctorType  m_work_functor ;
  const size_type    m_work_count ;

  ParallelReduce( const size_type work_count ,
                  const FunctorType & functor )
    : m_work_functor( functor )
    , m_work_count( work_count )
    {}

private:

  // self.m_work_count == total work count
  // work->count       == number of threads

  static void run_work_on_tpi( TPI_Work * work )
  {
    const ParallelReduce & self = *((const ParallelReduce *) work->info );
    value_type           & dst  = *((value_type *) work->reduce );

    const size_type work_inc   = (self.m_work_count + work->count - 1) / work->count ;
    const size_type work_begin = work_inc * work->rank ;
    const size_type work_end   = std::max( work_begin + work_inc , self.m_work_count );

    for ( size_type iwork = work_begin ; iwork < work_end ; ++iwork ) {
      self.m_work_functor( iwork , dst );
    }
  }

  static void run_init_on_tpi( TPI_Work * work )
  {
    value_type & dst = *((value_type *) work->reduce );
    FunctorType::init( dst );
  }

  static void run_join_on_tpi( TPI_Work * work , const void * reduce )
  {
    value_type & dst = *((value_type *) work->reduce );
    const value_type & src = *((const value_type *) reduce );
    FunctorType::join( dst , src );
  }

public:

  static value_type run( const size_type      work_count ,
                         const FunctorType  & work_functor ,
                               value_type   & result )
  {
    FunctorType::init( result );

    // Make a copy just like other devices will have to.

    device_type::set_dispatch_functor();

    const ParallelReduce tmp( work_count , work_functor );

    device_type::clear_dispatch_functor();

    TPI_Run_threads_reduce( & run_work_on_tpi , & tmp ,
                            & run_join_on_tpi ,
                            & run_init_on_tpi ,
                            sizeof(value_type) ,
                            & result );
  }
};

template< class FunctorType , class FinalizeType >
class ParallelReduce< FunctorType , FinalizeType , DeviceTPI > {
public:
  typedef DeviceTPI              device_type ;
  typedef device_type::size_type size_type ;
  typedef typename FunctorType::value_type value_type ;

  const FunctorType  m_work_functor ;
  const FinalizeType m_work_finalize ;
  const size_type    m_work_count ;

  ParallelReduce( const size_type work_count ,
                  const FunctorType & functor ,
                  const FinalizeType & finalize )
    : m_work_functor( functor )
    , m_work_finalize( finalize )
    , m_work_count( work_count )
    {}

private:

  // self.m_work_count == total work count
  // work->count       == number of threads

  static void run_work_on_tpi( TPI_Work * work )
  {
    const ParallelReduce & self = *((const ParallelReduce *) work->info );
    value_type           & dst  = *((value_type *) work->reduce );

    const size_type work_inc   = (self.m_work_count + work->count - 1) / work->count ;
    const size_type work_begin = work_inc * work->rank ;
    const size_type work_end   = std::min( work_begin + work_inc , self.m_work_count );

    for ( size_type iwork = work_begin ; iwork < work_end ; ++iwork ) {
      self.m_work_functor( iwork , dst );
    }
  }

  static void run_init_on_tpi( TPI_Work * work )
  {
    value_type & dst = *((value_type *) work->reduce );
    FunctorType::init( dst );
  }

  static void run_join_on_tpi( TPI_Work * work , const void * reduce )
  {
    value_type & dst = *((value_type *) work->reduce );
    const value_type & src = *((const value_type *) reduce );
    FunctorType::join( dst , src );
  }

public:

  static void run( const size_type      work_count ,
                   const FunctorType  & work_functor ,
                   const FinalizeType & work_finalize )
  {
    value_type result ;

    FunctorType::init( result );

    // Make a copy just like other devices will have to.

    device_type::set_dispatch_functor();

    const ParallelReduce tmp( work_count , work_functor , work_finalize );

    device_type::clear_dispatch_functor();

    TPI_Run_threads_reduce( & run_work_on_tpi , & tmp ,
                            & run_join_on_tpi ,
                            & run_init_on_tpi ,
                            sizeof(value_type) ,
                            & result );

    tmp.m_work_finalize( result );
  }
};

} // namespace Kokkos

#endif /* KOKKOS_DEVICETPI_PARALLELREDUCE_HPP */

