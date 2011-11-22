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

#ifndef KOKKOS_DEVICENUMA_PARALLELREDUCE_HPP
#define KOKKOS_DEVICENUMA_PARALLELREDUCE_HPP

#include <Kokkos_DeviceNUMA.hpp>
#include <Kokkos_ParallelReduce.hpp>

#include <algorithm>
#include <vector>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


namespace Kokkos {
namespace Impl {

template< class FunctorType , class ReduceTraits , class FinalizeType >
class ParallelReduce< FunctorType ,
                      ReduceTraits ,
                      FinalizeType ,
                      DeviceNUMA > : public DeviceNUMAWorker {
private:
  typedef          DeviceNUMA::size_type   size_type ;
  typedef typename ReduceTraits ::value_type  value_type ;

  const FunctorType   m_work_functor ;
  const FinalizeType  m_finalize ;
  const size_type     m_work_count ;

  // Virtual method defined in DeviceNUMAWorker
  void execute_on_thread( DeviceNUMAThread & this_thread ) const
  {
    value_type update ; // This thread's reduction value

    ReduceTraits::init( update );

    // Iterate this thread's work

    const std::pair<size_type,size_type> range =
      this_thread.work_range( m_work_count );

    for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
      m_work_functor( iwork , update );
    }

    // Fan-in reduction of other threads' reduction data:
    this_thread.reduce< ReduceTraits >( update );

    if ( 0 == this_thread.rank() ) {
      // Root of the fan-in reduction
      m_finalize( update );
    }
  }

  ParallelReduce( const size_type      work_count ,
                  const FunctorType  & functor ,
                  const FinalizeType & finalize )
    : DeviceNUMAWorker()
    , m_work_functor( functor )
    , m_finalize( finalize )
    , m_work_count( work_count )
    {}

public:

  static void execute( const size_type      work_count ,
                       const FunctorType  & functor ,
                       const FinalizeType & finalize )
  {
    DeviceNUMA::memory_space::set_dispatch_functor();

    ParallelReduce driver( work_count , functor , finalize );

    DeviceNUMA::memory_space::clear_dispatch_functor();

    DeviceNUMAWorker::execute( driver );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , typename ValueType >
class DeviceNUMAMultiFunctorParallelReduceMember ;

template< typename ValueType >
class DeviceNUMAMultiFunctorParallelReduceMember<void,ValueType> {
public:
  DeviceNUMAMultiFunctorParallelReduceMember() {}

  virtual ~DeviceNUMAMultiFunctorParallelReduceMember() {}

  virtual void execute_on_thread( DeviceNUMAThread & this_thread ,
                                  ValueType & update ) const = 0 ;
};
  
template< class FunctorType , typename ValueType >
class DeviceNUMAMultiFunctorParallelReduceMember
  : public DeviceNUMAMultiFunctorParallelReduceMember<void,ValueType> {
public:
  typedef DeviceNUMA::size_type size_type ;
    
  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  ~DeviceNUMAMultiFunctorParallelReduceMember() {}

  DeviceNUMAMultiFunctorParallelReduceMember(
    const FunctorType & work_functor ,
    const size_type work_count )
    : m_work_functor( work_functor )
    , m_work_count(   work_count )
    {}
    
  void execute_on_thread( DeviceNUMAThread & this_thread ,
                          ValueType & update ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( m_work_count );

    for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
      m_work_functor( iwork , update );
    }
  }
};  

} // namespace Impl
  
template< class ReduceTraits , class FinalizeType >
class MultiFunctorParallelReduce< ReduceTraits , FinalizeType , DeviceNUMA > : public Impl::DeviceNUMAWorker {
private:
  typedef          DeviceNUMA::size_type   size_type ;
  typedef typename ReduceTraits ::value_type value_type ;

  typedef Impl::DeviceNUMAMultiFunctorParallelReduceMember<void,value_type> MemberType ;

  typedef std::vector< MemberType * > MemberContainer ;

  typedef typename MemberContainer::const_iterator MemberIterator ;

  MemberContainer m_member_functors ;
  FinalizeType    m_finalize ;

  // Virtual method defined in DeviceNUMAWorker
  void execute_on_thread( Impl::DeviceNUMAThread & this_thread ) const
  {
    value_type update ; // This thread's reduction value

    ReduceTraits::init( update );

    for ( MemberIterator m  = m_member_functors.begin() ;
                         m != m_member_functors.end() ; ++m ) {
      (*m)->execute_on_thread( this_thread , update );
    }

    // Fan-in reduction of other threads' reduction data:
    this_thread.reduce< ReduceTraits >( update );

    if ( 0 == this_thread.rank() ) {
      // Root of the fan-in reduction
      m_finalize( update );
    }
  }

public:

  MultiFunctorParallelReduce( const FinalizeType & finalize )
    : m_member_functors()
    , m_finalize( finalize )
    {}

  ~MultiFunctorParallelReduce()
  {
    while ( ! m_member_functors.empty() ) {
      delete m_member_functors.back();
      m_member_functors.pop_back();
    }
  }

  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & functor )
  {
    typedef Impl::DeviceNUMAMultiFunctorParallelReduceMember<FunctorType,value_type> member_work_type ;

    MemberType * const m = new member_work_type( functor , work_count );

    m_member_functors.push_back( m );
  }

  void execute() const
  {
    DeviceNUMA::memory_space::set_dispatch_functor();
    DeviceNUMAWorker::execute( *this );
    DeviceNUMA::memory_space::clear_dispatch_functor();
  }
};

} // namespace Kokkos

#endif /* KOKKOS_DEVICENUMA_PARALLELREDUCE_HPP */

