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

#ifndef KOKKOSARRAY_HOST_PARALLELREDUCE_HPP
#define KOKKOSARRAY_HOST_PARALLELREDUCE_HPP

#include <KokkosArray_ParallelReduce.hpp>

#include <algorithm>
#include <vector>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename ValueType >
class FunctorAssignment< ValueType , Host >
{
public:

  ValueType & m_result ;

  FunctorAssignment( ValueType & result )
    : m_result( result ) {}

  void operator()( const ValueType & value ) const
    { m_result = value ; }
};

template< typename ValueType , class LayoutType >
class FunctorAssignment< View< ValueType , LayoutType , Host > , Host >
{
public:

  typedef View< ValueType , LayoutType , Host > view_type ;

  view_type m_result ;

  FunctorAssignment( const view_type & view )
    : m_result( view ) {}

  void operator()( const ValueType & value ) const
    { *m_result = value ; }
};

template< class FunctorType , class ReduceTraits , class FinalizeType >
class ParallelReduce< FunctorType , ReduceTraits , FinalizeType , Host >
  : public HostThreadWorker {
private:
  typedef          Host::size_type           size_type ;
  typedef typename ReduceTraits ::value_type value_type ;

  const FunctorType   m_work_functor ;
  const FinalizeType  m_finalize ;
  const size_type     m_work_count ;

  // Virtual method defined in HostThreadWorker
  void execute_on_thread( HostThread & this_thread ) const
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
    this_thread.reduce< ReduceTraits >( update , m_finalize );
    this_thread.return_barrier();
  }

  ParallelReduce( const size_type      work_count ,
                  const FunctorType  & functor ,
                  const FinalizeType & finalize )
    : HostThreadWorker()
    , m_work_functor( functor )
    , m_finalize( finalize )
    , m_work_count( work_count )
    {}

public:

  static void execute( const size_type      work_count ,
                       const FunctorType  & functor ,
                       const FinalizeType & finalize )
  {
    ParallelReduce driver( work_count , functor , finalize );

    HostThreadWorker::execute( driver );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class FunctorType , typename ValueType >
class HostMultiFunctorParallelReduceMember ;

template< typename ValueType >
struct HostMultiFunctorParallelReduceMember<void,ValueType> {

  virtual ~HostMultiFunctorParallelReduceMember() {}

  virtual void apply( HostThread & , ValueType & ) const = 0 ;
};


template< class FunctorType , typename ValueType >
class HostMultiFunctorParallelReduceMember
  : public HostMultiFunctorParallelReduceMember<void,ValueType> {
public:
  typedef Host::size_type size_type ;
    
  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  ~HostMultiFunctorParallelReduceMember() {}

  HostMultiFunctorParallelReduceMember(
    const FunctorType & work_functor ,
    const size_type work_count )
    : m_work_functor( work_functor )
    , m_work_count(   work_count )
    {}
    
  // virtual method
  void apply( HostThread & this_thread ,
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
class MultiFunctorParallelReduce< ReduceTraits , FinalizeType , Host >
  : public Impl::HostThreadWorker {
private:
  typedef          Host::size_type            size_type ;
  typedef typename ReduceTraits ::value_type  value_type ;
  typedef Impl::HostMultiFunctorParallelReduceMember<void,value_type> worker_type ;

  typedef std::vector< worker_type * > MemberContainer ;

  typedef typename MemberContainer::const_iterator MemberIterator ;

  MemberContainer m_member_functors ;
  FinalizeType    m_finalize ;

  // Virtual method defined in HostThreadWorker
  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    value_type update ; // This thread's reduction value

    ReduceTraits::init( update );

    for ( MemberIterator m  = m_member_functors.begin() ;
                         m != m_member_functors.end() ; ++m ) {
      (*m)->apply( this_thread , update );
    }

    // Fan-in reduction of other threads' reduction data:
    this_thread.reduce< ReduceTraits >( update , m_finalize );
    this_thread.return_barrier();
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
    typedef Impl::HostMultiFunctorParallelReduceMember<FunctorType,value_type> member_work_type ;

    worker_type * const m = new member_work_type( functor , work_count );

    m_member_functors.push_back( m );
  }

  void execute() const
  {
    Impl::HostThreadWorker::execute( *this );
  }
};

} // namespace KokkosArray

#endif /* KOKKOSARRAY_HOST_PARALLELREDUCE_HPP */

