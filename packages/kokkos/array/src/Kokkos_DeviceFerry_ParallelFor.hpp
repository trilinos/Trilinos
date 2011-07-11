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

#ifndef KOKKOS_DEVICEFERRY_PARALLELFOR_HPP
#define KOKKOS_DEVICEFERRY_PARALLELFOR_HPP

#include <Kokkos_ParallelFor.hpp>

#include <Kokkos_DeviceFerry_macros.hpp>

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

#include <stdio.h>
#include <algorithm>

#pragma offload_attribute(push, target(mic))
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task.h>
#pragma offload_attribute(pop)



namespace Kokkos {
namespace Impl {


template< class FunctorType >
class ParallelFor< FunctorType , DeviceFerry > {
public:
  typedef DeviceFerry::size_type size_type ;
  
  const FunctorType m_functor ;

#if 1
	class __declspec(target(mic)) ParallelTBB
	{
	public:
	const FunctorType * m_functor;

	void operator() (const tbb::blocked_range<size_type> & r) const {
	#ifdef __MIC__
		for(size_type i = r.begin() ; i != r.end(); ++i) {
			(*m_functor)(i);
		}
	#endif
	}
	ParallelTBB(const FunctorType *f) : m_functor(f) { }
	};
#endif

private:
  ParallelFor();
  ParallelFor(const FunctorType & f ) : m_functor( f )  { }  

public:
  static void execute( size_type work_count ,
                       const FunctorType & functor )
  {
    DeviceFerry::set_dispatch_functor();

    ParallelFor driver( functor );

    DeviceFerry::clear_dispatch_functor();
    
    ParallelFor *tmp_driver = &driver;
    
	//Perform on device
	
	#pragma offload target(mic) in(work_count) in(tmp_driver : length(1))
	{
#if 0
//		#pragma omp parallel for
		for(size_type i = 0; i <= work_count; ++i) {
			(tmp_functor->m_functor)(i);
		}
#else	
		ParallelTBB tmp_tbb(&(tmp_driver->m_functor));
		tbb::task_scheduler_init init;
		tbb::parallel_for(tbb::blocked_range<size_type>(0,work_count) , tmp_tbb , tbb::auto_partitioner());	
#endif
	}
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* defined (KOKKOS_MACRO_DEVICE_FUNCTION) */

#include <Kokkos_DeviceClear_macros.hpp>

#endif /* KOKKOS_DEVICEFerry_PARALLELREDUCE_HPP */


