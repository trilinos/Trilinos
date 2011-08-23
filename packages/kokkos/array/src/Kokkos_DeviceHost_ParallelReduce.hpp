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

#ifndef KOKKOS_DEVICEHOST_PARALLELREDUCE_HPP
#define KOKKOS_DEVICEHOST_PARALLELREDUCE_HPP

#include <Kokkos_ParallelReduce.hpp>
#include <Kokkos_DeviceHost.hpp>
#include <vector>

namespace Kokkos {
namespace Impl {

/** \brief For any serial device on host memory */
template< class FunctorType , class ReduceTraits , class MDArrayMap >
class ParallelReduce< FunctorType , ReduceTraits ,
                      void , Serial< HostMemory , MDArrayMap > > {
public:
  typedef DeviceHost::size_type             size_type ;
  typedef typename ReduceTraits::value_type  value_type ;

  const FunctorType  m_work_functor ;
  const size_type    m_work_count ;

private:

  ParallelReduce( const size_type     work_count ,
                  const FunctorType & functor )
    : m_work_functor( functor )
    , m_work_count( work_count )
    {}

public:

  static
  void execute( const size_type     work_count ,
                const FunctorType & functor ,
                      value_type  & result )
  {
    DeviceHost::set_dispatch_functor();

    const ParallelReduce driver( work_count , functor );

    DeviceHost::clear_dispatch_functor();

    FunctorType::init( result );

    for ( size_type iwork = 0 ; iwork < driver.m_work_count ; ++iwork ) {
      driver.m_work_functor(iwork,result);
    }
  }
};

template< class FunctorType , class ReduceTraits , class FinalizeType , class MDArrayMap >
class ParallelReduce< FunctorType , ReduceTraits , FinalizeType ,
                      Serial< HostMemory , MDArrayMap > > {
public:
  typedef DeviceHost::size_type             size_type ;
  typedef typename ReduceTraits::value_type  value_type ;

  const FunctorType  m_work_functor ;
  const FinalizeType m_finalize ;
  const size_type    m_work_count ;

private:

  ParallelReduce( const size_type      work_count ,
                  const FunctorType  & functor ,
                  const FinalizeType & finalize )
    : m_work_functor( functor )
    , m_finalize( finalize )
    , m_work_count( work_count )
  {}

public:

  static
  void execute( const size_type      work_count ,
                const FunctorType  & functor ,
                const FinalizeType & finalize )
  {
    DeviceHost::set_dispatch_functor();

    ParallelReduce driver( work_count , functor , finalize );

    DeviceHost::clear_dispatch_functor();

    value_type result ;

    FunctorType::init( result );

    for ( size_type iwork = 0 ; iwork < driver.m_work_count ; ++iwork ) {
      driver.m_work_functor(iwork,result);
    }

    driver.m_finalize( result );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template < class ReduceTraits , class FinalizeType >
class MultiFunctorParallelReduceMember< void , ReduceTraits , FinalizeType , DeviceHost > {
public:

  typedef typename ReduceTraits::value_type value_type ;

  virtual ~MultiFunctorParallelReduceMember() {}

  virtual void execute( value_type & update ) const = 0 ;
};

template < class FunctorType , class ReduceTraits , class FinalizeType >
class MultiFunctorParallelReduceMember< FunctorType , ReduceTraits , FinalizeType , DeviceHost >
  : public MultiFunctorParallelReduceMember< void , ReduceTraits , FinalizeType , DeviceHost >
{
public:
  typedef          DeviceHost  ::size_type  size_type ;
  typedef typename ReduceTraits::value_type value_type ;

  FunctorType m_functor ;
  size_type   m_work_count ;

  MultiFunctorParallelReduceMember( const size_type     work_count ,
                                    const FunctorType & functor )
    : m_functor( functor )
    , m_work_count( work_count )
    {}

  void execute( value_type & update ) const
  {
    for ( size_type iwork = 0 ; iwork < m_work_count ; ++iwork ) {
      m_functor( iwork , update );
    }
  }
};

} // namespace Impl

template< class ReduceTraits >
class MultiFunctorParallelReduce< ReduceTraits , typename ReduceTraits::value_type , DeviceHost > {
private:
  typedef          DeviceHost  ::size_type   size_type ;
  typedef typename ReduceTraits::value_type value_type ;

  typedef MultiFunctorParallelReduce< ReduceTraits , value_type , DeviceHost > self_type ;
  typedef Impl::MultiFunctorParallelReduceMember<void,ReduceTraits,value_type,DeviceHost> MemberType ;
  typedef std::vector< MemberType * > MemberContainer ;
  typedef typename MemberContainer::const_iterator MemberIterator ;

  MemberContainer m_member_functors ;

public:

  value_type result ;

  MultiFunctorParallelReduce()
    : m_member_functors()
    , result()
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
    MemberType * m = new Impl::MultiFunctorParallelReduceMember<FunctorType,ReduceTraits,value_type,DeviceHost>( work_count , functor );
    m_member_functors.push_back( m );
  }

  void execute()
  {
    ReduceTraits::init( result );
    for ( MemberIterator m  = m_member_functors.begin();
                         m != m_member_functors.end(); ++m ) {
      (*m)->execute( result );
    }
  }
};

template< class ReduceTraits , class FinalizeType >
class MultiFunctorParallelReduce< ReduceTraits , FinalizeType , DeviceHost > {
private:

  MultiFunctorParallelReduce< ReduceTraits , typename ReduceTraits::value_type , DeviceHost > m_impl ;

public:

  FinalizeType result ;

  MultiFunctorParallelReduce() : m_impl() {}

  template< class FunctorType >
  void push_back( const DeviceHost::size_type work_count ,
                  const FunctorType & functor )
  { m_impl.push_back( work_count , functor ); }

  void execute()
  {
    m_impl.execute();
    result( m_impl.result );
  }
};

} // namespace Kokkos

#endif /* KOKKOS_DEVICEHOST_PARALLELREDUCE_HPP */

