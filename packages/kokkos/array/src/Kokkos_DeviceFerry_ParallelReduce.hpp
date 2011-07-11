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

#ifndef KOKKOS_DEVICEFERRY_PARALLELREDUCE_HPP
#define KOKKOS_DEVICEFERRY_PARALLELREDUCE_HPP

#include <Kokkos_ParallelReduce.hpp>

#include <vector>
#include <iostream>

#include <Kokkos_DeviceFerry_macros.hpp>

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

#pragma offload_attribute(push, target(mic))
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task.h>
#pragma offload_attribute(pop)


namespace Kokkos {
namespace Impl {

template< class FunctorType , class ReduceTraits ,class FinalizeType >
class ParallelReduce< FunctorType , ReduceTraits , FinalizeType , DeviceFerry > {
public:
  typedef DeviceFerry::size_type size_type ;
  typedef typename FunctorType::value_type value_type ;

	const FunctorType m_functor;
	const FinalizeType m_finalize;
	
private:
	ParallelReduce();
	ParallelReduce(const FunctorType &func, const FinalizeType &final) : m_functor(func) , m_finalize(final) { }
	
public:	
class __declspec(target(mic)) ParallelTBB
{
public:
	const ReduceTraits * m_functor ;
	value_type m_value;
	
	ParallelTBB(const FunctorType *f) : m_functor(f) 
	{
		FunctorType::init(m_value);
	}

public:


	ParallelTBB() : m_functor(NULL) 
	{
		FunctorType::init(m_value);
	}
	
	ParallelTBB(ParallelTBB & rhs , tbb::split) : m_functor(rhs.m_functor) 
	{
		FunctorType::init(m_value);
	}
		
	void operator() (const tbb::blocked_range<size_type> & range)
	{
		#ifdef __MIC__
		for(size_type i = range.begin() ; i != range.end() ; ++i)
		{
			(*m_functor)(i,m_value);
		}
		#endif
	}
	
	void join ( ParallelTBB & rhs) 
	{
		FunctorType::join(m_value, rhs.m_value);
	}
};	
public:

  static void execute( size_type work_count ,
                       const FunctorType & functor ,
                       const FinalizeType & finalize )
  {
    DeviceFerry::set_dispatch_functor();
    // Make a locally owned copy for consistent behavior
    // across all devices.

    ParallelReduce driver( functor , finalize);

    DeviceFerry::clear_dispatch_functor();
    
    ParallelReduce *tmp_driver = &driver;

//Perform on device
	#pragma offload target(mic) in(work_count) inout(tmp_driver : length(1)) 
	{
		ParallelTBB tmp_tbb(&(tmp_driver->m_functor));
		tbb::task_scheduler_init init;
		tbb::parallel_reduce(tbb::blocked_range<size_type>(0,work_count) , tmp_tbb , tbb::auto_partitioner() );	   

		tmp_driver->m_finalize(tmp_tbb.m_value);
    }
  }
};

//----------------------------------------------------------------------------

template< class FunctorType , class ReduceTraits >
class ParallelReduce< FunctorType , ReduceTraits , void , DeviceFerry > 
{
public:
  typedef DeviceFerry::size_type             size_type ;
  typedef typename FunctorType::value_type value_type ;

  struct AssignValueFunctor {

    value_type & ref ;

    AssignValueFunctor( value_type & arg_ref ) : ref( arg_ref ) {}

    AssignValueFunctor( const AssignValueFunctor & rhs ) : ref( rhs.ref ) {}

    void operator()( const value_type & val ) const { ref = val ; }
  };

  static void execute( const size_type     work_count ,
                       const FunctorType & functor ,
                             value_type  & result )
  {
    ParallelReduce< FunctorType , ReduceTraits AssignValueFunctor , DeviceFerry >
      ::execute( work_count , functor , AssignValueFunctor( result ) );
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* defined (KOKKOS_MACRO_DEVICE_FUNCTION) */

#include <Kokkos_DeviceClear_macros.hpp>

#endif /* KOKKOS_DEVICEFerry_PARALLELREDUCE_HPP */

