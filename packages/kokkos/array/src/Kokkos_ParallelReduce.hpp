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

#ifndef KOKKOS_PARALLELREDUCE_HPP
#define KOKKOS_PARALLELREDUCE_HPP

#include <cstddef>

#include <Kokkos_ArrayForwardDeclarations.hpp>
#include <Kokkos_ValueView.hpp>
#include <impl/Kokkos_Timer.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------

namespace Impl {

template< class FunctorType ,
          class FinalizeType ,
          class DeviceType = typename FunctorType::device_type >
class ParallelReduce {
public:
  // non-void FinalizeType post-processes the reduction result
  static void execute( const size_t work_count ,
                       const FunctorType  & functor ,
                       const FinalizeType & finalize );

  // Void FinalizeType returns the reduction result
  static void execute( const size_t work_count ,
                       const FunctorType  & functor ,
                       typename FunctorType::value_type & result );
};

} // namespace Impl

//----------------------------------------------------------------------------

template< class FunctorType >
typename FunctorType::value_type
parallel_reduce( const size_t work_count ,
                 const FunctorType & functor )
{
  typename FunctorType::value_type result ;

  Impl::ParallelReduce< FunctorType , void >
    ::execute( work_count , functor , result );

  return result ;
}

template< class FunctorType >
void parallel_reduce( const size_t work_count ,
                      const FunctorType & functor ,
                      typename FunctorType::value_type & result )
{
  Impl::ParallelReduce< FunctorType , void >
    ::execute( work_count , functor , result );
}

template< class FunctorType , class FinalizeType >
void parallel_reduce( const size_t work_count ,
                      const FunctorType & functor ,
                      const FinalizeType & finalize )
{
  Impl::ParallelReduce< FunctorType , FinalizeType >
    ::execute( work_count , functor , finalize );
}

template< class FunctorType >
void parallel_reduce( const size_t work_count ,
                      const FunctorType & functor ,
                      typename FunctorType::value_type & result ,
                      double & seconds )
{
  typedef typename FunctorType::device_type device_type ;

  const Impl::Timer timer ;

  Impl::ParallelReduce< FunctorType , void >
    ::execute( work_count , functor , result );

  device_type::wait_functor_completion();

  seconds = timer.seconds() ;
}

template< class FunctorType , class FinalizeType >
void parallel_reduce( const size_t work_count ,
                      const FunctorType & functor ,
                      const FinalizeType & finalize ,
                      double & seconds )
{
  typedef typename FunctorType::device_type device_type ;

  const Impl::Timer timer ;

  Impl::ParallelReduce< FunctorType , FinalizeType >
    ::execute( work_count , functor , finalize );

  device_type::wait_functor_completion();

  seconds = timer.seconds() ;
}

} // namespace Kokkos

//----------------------------------------------------------------------------
// Partial specializations for known devices

#if defined( KOKKOS_DEVICE_HOST )
#include <DeviceHost/Kokkos_DeviceHost_ParallelReduce.hpp>
#endif

#if defined( KOKKOS_DEVICE_TPI )
#include <DeviceTPI/Kokkos_DeviceTPI_ParallelReduce.hpp>
#endif

#if defined( KOKKOS_DEVICE_CUDA )
#include <DeviceCuda/Kokkos_DeviceCuda_ParallelReduce.hpp>
#endif

//----------------------------------------------------------------------------

#endif /* KOKKOS_DEVICEHOST_PARALLELREDUCE_HPP */

