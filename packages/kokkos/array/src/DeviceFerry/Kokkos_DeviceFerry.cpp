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
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

//#pragma offload_attribute ( push , target(mic)) 
#include <Kokkos_DeviceFerry.hpp>
#include <impl/Kokkos_MemoryInfo.hpp>
//#pragma offload_attribute(pop)
namespace Kokkos {

namespace {
/*--------------------------------------------------------------------------*/

class DeviceFerry_Impl {
public:
	Impl::MemoryInfoSet m_allocations ;

	DeviceFerry_Impl();
	~DeviceFerry_Impl();

	static DeviceFerry_Impl & singleton();

	void * allocate_memory(
		const std::string    & label ,
		const std::type_info & type ,
		const size_t member_size ,
		const size_t member_count );

	void deallocate_memory( void * ptr );

private:
	DeviceFerry_Impl( const DeviceFerry_Impl & );
	DeviceFerry_Impl & operator = ( const DeviceFerry_Impl & );
};

DeviceFerry_Impl & DeviceFerry_Impl::singleton()
{
  static DeviceFerry_Impl self;
  return self ;
}

DeviceFerry_Impl::DeviceFerry_Impl( )	: m_allocations() { }

DeviceFerry_Impl::~DeviceFerry_Impl()
{
	if ( ! m_allocations.empty() ) {
	std::cerr << "Kokkos::DeviceFerry memory leaks:" << std::endl ;
	m_allocations.print( std::cerr );
	}
}

void * DeviceFerry_Impl::allocate_memory(
	const std::string    & label ,
	const std::type_info & type ,
	const size_t member_size ,
	const size_t member_count )
{
	Impl::MemoryInfo tmp ;

	tmp.m_type  = & type ;
	tmp.m_label = label ;
	tmp.m_size  = member_size ;
	tmp.m_count = member_count ;
	
	long int view;
	size_t view_size = member_size * member_count;
	
	#pragma offload target(mic) in(view_size) out(view)
	{
		long int local_tmp = (long int)malloc(view_size);
		view = local_tmp;
	}
	tmp.m_ptr =(void*) view;


	const bool ok_alloc  = 0 != tmp.m_ptr ;
	const bool ok_insert = ok_alloc && m_allocations.insert( tmp );

	if ( ! ok_alloc || ! ok_insert ) {
	std::ostringstream msg ;
	msg << "Kokkos::DeviceFerry::allocate_memory( " << label
		<< " , " << type.name()
		<< " , " << member_size
		<< " , " << member_count
		<< " ) FAILED " ;
	if ( ok_alloc ) { msg << "memory allocation" ; }
	else            { msg << "with internal error" ; }
	throw std::runtime_error( msg.str() );
	}

	return tmp.m_ptr ;
}

void DeviceFerry_Impl::deallocate_memory( void * ptr )
{
	if ( ! m_allocations.erase( ptr ) ) {
	std::ostringstream msg ;
	msg << "Kokkos::DeviceFerry::deallocate_memory( " << ptr
		<< " ) FAILED memory allocated by this device" ;
	throw std::runtime_error( msg.str() );
	}
  	long int local_tmp = (long int)ptr;
	#pragma offload target(mic) in(local_tmp)
	{
		free((void*)local_tmp);
	}
}

}

/*--------------------------------------------------------------------------*/

void DeviceFerry::initialize( )
{
  DeviceFerry_Impl::singleton( );
}

void * DeviceFerry::allocate_memory(
  const std::string    & label ,
  const std::type_info & type ,
  const size_t member_size ,
  const size_t member_count )
{
	DeviceFerry_Impl & s = DeviceFerry_Impl::singleton();

	return s.allocate_memory( label ,type , member_size , member_count );
}

void DeviceFerry::deallocate_memory( void * ptr )
{
	DeviceFerry_Impl & s = DeviceFerry_Impl::singleton();

	s.deallocate_memory( ptr );
}

void DeviceFerry::print_memory_view( std::ostream & o )
{
	DeviceFerry_Impl & s = DeviceFerry_Impl::singleton();

	s.m_allocations.print( o );
}

/*--------------------------------------------------------------------------*/

unsigned int DeviceFerry::m_launching_kernel = false ;

void DeviceFerry::set_dispatch_functor()
{
	if ( m_launching_kernel ) {
	std::string msg ;
	msg.append( "Kokkos::DeviceFerry::set_dispatch_functor FAILED: " );
	msg.append( "kernel dispatch is already in progress, " );
	msg.append( "a recursive call or forgotten 'clear_dispatch_kernel" );
	throw std::runtime_error( msg );
	}
	m_launching_kernel = true ;
}

void DeviceFerry::clear_dispatch_functor()
{
	if ( ! m_launching_kernel ) {
	std::string msg ;
	msg.append( "Kokkos::DeviceFerry::clear_dispatch_functor FAILED: " );
	msg.append( "no kernel dispatch in progress." );
	throw std::runtime_error( msg );
	}
	m_launching_kernel = false ;
}

void DeviceFerry::wait_functor_completion()
{
  //Ferry synch thread ?
  //cudaThreadSynchronize();
}

} // namespace Kokkos 
