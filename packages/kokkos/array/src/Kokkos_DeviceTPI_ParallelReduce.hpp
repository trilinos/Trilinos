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

#include <Kokkos_ParallelReduce.hpp>

#include <algorithm>
#include <vector>
#include <TPI.h>

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ReduceTraits >
class ParallelReduce< FunctorType , ReduceTraits , typename ReduceTraits::value_type , DeviceTPI > {
public:
  typedef          DeviceTPI   ::size_type   size_type ;
  typedef typename ReduceTraits::value_type  value_type ;

  const FunctorType  m_work_functor ;
  const size_type    m_work_count ;

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
    ReduceTraits::init( dst );
  }

  static void run_join_on_tpi( TPI_Work * work , const void * reduce )
  {
    volatile value_type & dst = *((value_type *) work->reduce );
    const volatile value_type & src = *((const value_type *) reduce );
    ReduceTraits::join( dst , src );
  }

  ParallelReduce( const size_type work_count ,
                  const FunctorType & functor )
    : m_work_functor( functor )
    , m_work_count( work_count )
    {}

public:

  static void execute( const size_type     work_count ,
                       const FunctorType & functor ,
                             value_type  & result )
  {
    DeviceTPI::set_dispatch_functor();

    ParallelReduce driver( work_count , functor );

    DeviceTPI::clear_dispatch_functor();

    ReduceTraits::init( result );

    TPI_Run_threads_reduce( & run_work_on_tpi , & driver ,
                            & run_join_on_tpi ,
                            & run_init_on_tpi ,
                            sizeof(value_type) ,
                            & result );
  }
};


template< class FunctorType , class ReduceTraits , class FinalizeType >
class ParallelReduce< FunctorType , ReduceTraits , FinalizeType , DeviceTPI > {
public:

  typedef          DeviceTPI   ::size_type   size_type ;
  typedef typename ReduceTraits::value_type  value_type ;

  static void execute( const size_type      work_count ,
                       const FunctorType  & functor ,
                       const FinalizeType & finalize )
  {
    value_type result ;

    ParallelReduce< FunctorType , ReduceTraits , value_type , DeviceTPI >
      ::execute( work_count , functor , result );

    finalize( result );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class ReduceTraits , class FinalizeType >
class MultiFunctorParallelReduceMember<void,ReduceTraits,FinalizeType,DeviceTPI> {
public:
  typedef          DeviceTPI   ::size_type   size_type ;
  typedef typename ReduceTraits::value_type  value_type ;

  virtual ~MultiFunctorParallelReduceMember() {}

  virtual void execute( const size_type thread_count ,
                        const size_type thread_rank ,
                        value_type & update ) const = 0 ;
};

template< class FunctorType , class ReduceTraits , class FinalizeType >
class MultiFunctorParallelReduceMember< FunctorType , ReduceTraits , FinalizeType , DeviceTPI >
  : public MultiFunctorParallelReduceMember<void,ReduceTraits,FinalizeType,DeviceTPI> {
public:
  typedef          DeviceTPI   ::size_type   size_type ;
  typedef typename ReduceTraits::value_type  value_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  MultiFunctorParallelReduceMember( const size_type work_count ,
                                    const FunctorType & work_functor )
    : m_work_functor( work_functor )
    , m_work_count(   work_count )
    {}

  virtual void execute( const size_type thread_count ,
                        const size_type thread_rank ,
                        value_type & update ) const
  {
    const size_type work_inc   = (m_work_count + thread_count - 1) / thread_count ;
    const size_type work_begin = work_inc * thread_rank ;
    const size_type work_end   = std::min( work_begin + work_inc , m_work_count );

    for ( size_type iwork = work_begin ; iwork < work_end ; ++iwork ) {
      m_work_functor( iwork , update );
    }
  }
};

} // namespace Impl

template< class ReduceTraits >
class MultiFunctorParallelReduce< ReduceTraits , typename ReduceTraits::value_type , DeviceTPI > {
private:
  typedef          DeviceTPI   ::size_type   size_type ;
  typedef typename ReduceTraits::value_type value_type ;

  typedef MultiFunctorParallelReduce< ReduceTraits , value_type , DeviceTPI > self_type ;
  typedef Impl::MultiFunctorParallelReduceMember<void,ReduceTraits,value_type,DeviceTPI> MemberType ;

  typedef std::vector< MemberType * > MemberContainer ;
  typedef typename MemberContainer::const_iterator MemberIterator ;

  MemberContainer m_member_functors ;

  // self.m_work_count == total work count
  // work->count       == number of threads

  static void run_work_on_tpi( TPI_Work * work )
  {
    const self_type & self = *((const self_type *) work->info );
    value_type      & dst  = *((value_type *) work->reduce );

    for ( MemberIterator m  = self.m_member_functors.begin() ;
                         m != self.m_member_functors.end() ; ++m ) {

      (*m)->execute( work->count , work->rank , dst );
    }
  }

  static void run_init_on_tpi( TPI_Work * work )
  {
    value_type & dst = *((value_type *) work->reduce );
    ReduceTraits::init( dst );
  }

  static void run_join_on_tpi( TPI_Work * work , const void * reduce )
  {
    volatile value_type & dst = *((value_type *) work->reduce );
    const volatile value_type & src = *((const value_type *) reduce );
    ReduceTraits::join( dst , src );
  }

public:

  value_type result ;

  MultiFunctorParallelReduce()
    : m_member_functors()
    , result()
    {} 

  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & functor )
  {
    MemberType * m = new Impl::MultiFunctorParallelReduceMember<FunctorType,ReduceTraits,value_type,DeviceTPI>( work_count , functor );
    m_member_functors.push_back( m );
  }

  void execute()
  {
    ReduceTraits::init( result );

    TPI_Run_threads_reduce( & run_work_on_tpi , this  ,
                            & run_join_on_tpi ,
                            & run_init_on_tpi ,
                            sizeof(value_type) ,
                            & result );
  }
};

template< class ReduceTraits , class FinalizeType >
class MultiFunctorParallelReduce< ReduceTraits , FinalizeType , DeviceTPI > {
public:

  typedef          DeviceTPI   ::size_type   size_type ;
  typedef typename ReduceTraits::value_type value_type ;

private:

  MultiFunctorParallelReduce< ReduceTraits , value_type , DeviceTPI > m_impl ;

public:

  FinalizeType result ;

  MultiFunctorParallelReduce() : m_impl() {} 

  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & functor )
  { m_impl.push_back( work_count , functor ); }

  void execute()
  {
    m_impl.execute();
    result( m_impl.result );
  }
};

} // namespace Kokkos

#endif /* KOKKOS_DEVICETPI_PARALLELREDUCE_HPP */

