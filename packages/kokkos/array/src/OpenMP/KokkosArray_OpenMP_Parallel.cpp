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

#include <string>
#include <iostream>
#include <stdexcept>
#include <KokkosArray_OpenMP.hpp>
#include <KokkosArray_hwloc.hpp>

namespace KokkosArray {

namespace {

std::pair<unsigned,unsigned> coordinates[ Impl::HostThread::max_thread_count ] ;

}

void OpenMP::assert_not_in_parallel( const char * const function )
{
  if ( omp_in_parallel() ) {
    std::string msg(function);
    msg.append(" ERROR : Cannot be called OMP parallel");
    throw std::runtime_error(msg);
  }
}


void OpenMP::initialize()
{
  assert_not_in_parallel("KokkosArray::OpenMP::initialize");

  const bool ok_inactive = 0 == Impl::HostThread::get_thread_count();

  if ( ok_inactive ) {
#if 1
    const int thread_count  = omp_get_max_threads();
    const int thread_levels = omp_get_max_active_levels();

    std::cout << "KokkosArray::OpenMP::initialize :"
              << " count[" << thread_count << "]"
              << " levels[" << thread_levels << "]"
              << std::endl ;
#endif

#pragma omp parallel
    {
#pragma omp critical
      {
        const int count = omp_get_num_threads();
        const int rank = omp_get_thread_num();
        // Impl::hwloc::bind_this_thread( coordinates[ rank ] );

        Impl::HostThread * const th = new Impl::HostThread();
        Impl::HostThread::set_thread( rank , th );
        th->set_gang_worker( 0 , 1 , rank , count );
      }
    }
// END #pragma omp parallel

    Impl::HostThread::set_thread_relationships();
  }
  else {
    std::ostringstream msg ;

    msg << "KokkosArray::OpenMP::initialize() FAILED" ;

    if ( ! ok_inactive ) {
      msg << " : Device is already active" ;
    }

    KokkosArray::Impl::throw_runtime_exception( msg.str() );
  }
}

void OpenMP::finalize()
{
  assert_not_in_parallel("KokkosArray::OpenMP::finalize");

  resize_reduce_scratch(0);

#pragma omp parallel
  {
    Impl::HostThread * const th =
      Impl::HostThread::clear_thread( omp_get_thread_num() );

    delete th ;
  }
}

//----------------------------------------------------------------------------

void OpenMP::resize_reduce_scratch( unsigned size )
{
  static unsigned m_reduce_size = 0 ;

  const unsigned rem = size % HostSpace::MEMORY_ALIGNMENT ;

  if ( rem ) size += HostSpace::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size ) || ( m_reduce_size < size ) ) {

#pragma omp parallel
    {
#pragma omp critical
      {
        Impl::HostThread::get_thread( omp_get_thread_num() )->resize_reduce(size);
      }
    }

    m_reduce_size = size ;
  }
}

void * OpenMP::root_reduce_scratch()
{
  Impl::HostThread * const th =  Impl::HostThread::get_thread(0);
  return th ? th->reduce_data() : 0 ;
}

} // namespace KokkosArray


