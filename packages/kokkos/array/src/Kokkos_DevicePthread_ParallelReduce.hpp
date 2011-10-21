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

#ifndef KOKKOS_DEVICEPTHREAD_PARALLELREDUCE_HPP
#define KOKKOS_DEVICEPTHREAD_PARALLELREDUCE_HPP

#include <Kokkos_DevicePthread.hpp>
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
                      DevicePthread > : public DevicePthreadWorker {
public:
  typedef          DevicePthread::size_type   size_type ;
  typedef typename ReduceTraits ::value_type  value_type ;

  const FunctorType   m_work_functor ;
  const FinalizeType  m_finalize ;

  // Virtual method defined in DevicePthreadWorker
  void execute_on_thread( DevicePthreadController & this_thread ) const
  {
    value_type update ; // This thread's reduction value

    ReduceTraits::init( update );

    // Iterate this thread's work
    size_type iwork = DevicePthreadWorker::m_work_portion * this_thread.rank();

    const size_type work_end =
      std::min( iwork + DevicePthreadWorker::m_work_portion ,
                        DevicePthreadWorker::m_work_count );

    for ( ; iwork < work_end ; ++iwork ) {
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
    : DevicePthreadWorker( work_count )
    , m_work_functor( functor )
    , m_finalize( finalize )
    {}

public:

  static void execute( const size_type      work_count ,
                       const FunctorType  & functor ,
                       const FinalizeType & finalize )
  {
    DevicePthread::memory_space::set_dispatch_functor();

    ParallelReduce driver( work_count , functor , finalize );

    DevicePthread::memory_space::clear_dispatch_functor();

    DevicePthread::execute( driver );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_DEVICEPTHREAD_PARALLELREDUCE_HPP */

