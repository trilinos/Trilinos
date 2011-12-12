/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_DEVICEPTHREAD_PARALLELFOR_HPP
#define KOKKOS_DEVICEPTHREAD_PARALLELFOR_HPP

#include <Kokkos_ParallelFor.hpp>

#include <algorithm>
#include <vector>

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class ParallelFor< FunctorType , DevicePthread > : public DevicePthreadWorker {
private:

  typedef DevicePthread::size_type size_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;
  const size_type   m_work_per_thread ;

  virtual
  void execute_on_thread( Impl::DevicePthreadController & this_thread ) const
  {
    // Iterate this thread's work
    size_type iwork = m_work_per_thread * this_thread.rank();

    const size_type work_end =
      std::min( iwork + m_work_per_thread , m_work_count );

    for ( ; iwork < work_end ; ++iwork ) {
      m_work_functor( iwork );
    }

    this_thread.barrier();
  }

  ParallelFor( const size_type work_count ,
               const FunctorType & functor )
    : DevicePthreadWorker()
    , m_work_functor( functor )
    , m_work_count( work_count )
    , m_work_per_thread( DevicePthreadWorker::work_per_thread( work_count ) )
    {}

public:

  static void execute( const size_type     work_count ,
                       const FunctorType & functor )
  {
    DevicePthread::memory_space::set_dispatch_functor();

    ParallelFor driver( work_count , functor );

    DevicePthread::memory_space::clear_dispatch_functor();

    DevicePthread::execute( driver );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class DevicePthreadMultiFunctorParallelForMember ;

template<>
class DevicePthreadMultiFunctorParallelForMember<void> {
public:
  DevicePthreadMultiFunctorParallelForMember() {}

  virtual ~DevicePthreadMultiFunctorParallelForMember() {}

  virtual void execute_on_thread( DevicePthread::size_type thread_rank ) const = 0 ;
};

template< class FunctorType >
class DevicePthreadMultiFunctorParallelForMember
  : public DevicePthreadMultiFunctorParallelForMember<void> {
public:
  typedef DevicePthread::size_type size_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;
  const size_type   m_work_per_thread ;

  ~DevicePthreadMultiFunctorParallelForMember() {}

  DevicePthreadMultiFunctorParallelForMember(
    const FunctorType & work_functor ,
    const size_type work_count ,
    const size_type work_per_thread )
    : m_work_functor( work_functor )
    , m_work_count(   work_count )
    , m_work_per_thread( work_per_thread )
    {}

  void execute_on_thread( size_type thread_rank ) const
  {
    // Iterate this thread's work
    size_type iwork = m_work_per_thread * thread_rank ;

    const size_type work_end = std::min( iwork + m_work_per_thread , m_work_count );

    for ( ; iwork < work_end ; ++iwork ) {
      m_work_functor( iwork );
    }
  }
};

} // namespace Impl

template<>
class MultiFunctorParallelFor< DevicePthread > : public Impl::DevicePthreadWorker {
private:
  typedef DevicePthread::size_type size_type ;
  
  typedef Impl::DevicePthreadMultiFunctorParallelForMember<void> MemberType ;

  typedef std::vector< MemberType * > MemberContainer ;

  typedef MemberContainer::const_iterator MemberIterator ;

  MemberContainer m_member_functors ;

  // Virtual method defined in DevicePthreadWorker
  void execute_on_thread( Impl::DevicePthreadController & this_thread ) const
  {
    const size_type thread_rank = this_thread.rank();

    for ( MemberIterator m  = m_member_functors.begin() ;
                         m != m_member_functors.end() ; ++m ) {
      (*m)->execute_on_thread( thread_rank );
    }

    this_thread.barrier();
  }
  
public: 

  MultiFunctorParallelFor() : m_member_functors() {}
    
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
    typedef Impl::DevicePthreadMultiFunctorParallelForMember<FunctorType> member_work_type ;

    DevicePthread::memory_space::set_dispatch_functor();
    MemberType * const m = new member_work_type( functor , work_count , Impl::DevicePthreadWorker::work_per_thread( work_count ) );
    DevicePthread::memory_space::clear_dispatch_functor();

    m_member_functors.push_back( m );
  }

  void execute() const { DevicePthread::execute( *this ); }
};

} // namespace Kokkos

#endif /* KOKKOS_DEVICEPTHREAD_PARALLELFOR_HPP */

