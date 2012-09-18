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

#ifndef KOKKOSARRAY_HOST_PARALLELFOR_HPP
#define KOKKOSARRAY_HOST_PARALLELFOR_HPP

#include <KokkosArray_ParallelFor.hpp>

#include <algorithm>
#include <vector>

namespace KokkosArray {
namespace Impl {

template< class FunctorType >
class ParallelFor< FunctorType , Host > : public HostThreadWorker<void> {
private:

  typedef Host::size_type  size_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  // virtual function from HostThreadWorker:
  void execute_on_thread( HostThread & this_thread ) const
  {
    const std::pair< size_type , size_type > range =
      this_thread.work_range( m_work_count );

    for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
      m_work_functor( iwork );
    }

    this_thread.barrier();
  }

  ParallelFor( const size_type work_count , const FunctorType & functor )
    : HostThreadWorker<void>()
    , m_work_functor( functor )
    , m_work_count( work_count )
    {}

public:

  inline
  static void execute( const size_type     work_count ,
                       const FunctorType & functor )
  {
    ParallelFor driver( work_count , functor );

    HostThreadWorker<void>::execute( driver );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class FunctorType >
class HostMultiFunctorParallelForMember : public HostThreadWorker<void> {
public:
  typedef Host::size_type size_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  ~HostMultiFunctorParallelForMember() {}

  HostMultiFunctorParallelForMember(
    const FunctorType & work_functor ,
    const size_type work_count )
    : m_work_functor( work_functor )
    , m_work_count(   work_count )
    {}

  void execute_on_thread( HostThread & this_thread ) const
  {
    const std::pair< size_type , size_type > range =
      this_thread.work_range( m_work_count );

    for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
      m_work_functor( iwork );
    }
  }
};

} // namespace Impl

template<>
class MultiFunctorParallelFor< Host >
  : public Impl::HostThreadWorker<void> {
private:
  typedef Host::size_type               size_type ;
  typedef Impl::HostThreadWorker<void>  worker_type ;
  
  typedef std::vector< worker_type * > MemberContainer ;

  typedef MemberContainer::const_iterator MemberIterator ;

  MemberContainer m_member_functors ;

  // Virtual method defined in HostThreadWorker
  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    for ( MemberIterator m  = m_member_functors.begin() ;
                         m != m_member_functors.end() ; ++m ) {
      (*m)->execute_on_thread( this_thread );
    }

    this_thread.barrier();
  }
  
public: 

  MultiFunctorParallelFor()
    : Impl::HostThreadWorker<void>()
    , m_member_functors() {}
    
  ~MultiFunctorParallelFor()
  { 
    while ( ! m_member_functors.empty() ) {
      delete m_member_functors.back();
      m_member_functors.pop_back();
    }
  }
  
  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & functor )
  {
    typedef Impl::HostMultiFunctorParallelForMember<FunctorType> member_work_type ;

    worker_type * const m = new member_work_type( functor , work_count );

    m_member_functors.push_back( m );
  }

  void execute() const
  {
    Impl::HostThreadWorker<void>::execute( *this );
  }
};

} // namespace KokkosArray

#endif /* KOKKOSARRAY_HOST_PARALLELFOR_HPP */

