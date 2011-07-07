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

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ReduceTraits >
class ParallelReduce< FunctorType , ReduceTraits , void , DeviceHost > {
public:
  typedef DeviceHost                         device_type ;
  typedef device_type::size_type             size_type ;
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
    device_type::set_dispatch_functor();

    const ParallelReduce driver( work_count , functor );

    device_type::clear_dispatch_functor();

    FunctorType::init( result );

    for ( size_type iwork = 0 ; iwork < driver.m_work_count ; ++iwork ) {
      driver.m_work_functor(iwork,result);
    }
  }
};

template< class FunctorType , class ReduceTraits , class FinalizeType >
class ParallelReduce< FunctorType , ReduceTraits , FinalizeType , DeviceHost > {
public:
  typedef DeviceHost                         device_type ;
  typedef device_type::size_type             size_type ;
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
    device_type::set_dispatch_functor();

    ParallelReduce driver( work_count , functor , finalize );

    device_type::clear_dispatch_functor();

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

#endif /* KOKKOS_DEVICEHOST_PARALLELREDUCE_HPP */

