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

#ifndef KOKKOS_DEVICETBB_PARALLELFOR_HPP
#define KOKKOS_DEVICETBB_PARALLELFOR_HPP

#include <Kokkos_ParallelFor.hpp>

#include <algorithm>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class ParallelFor< FunctorType , DeviceTBB > {
public:
  typedef DeviceTBB::size_type size_type ;

  const FunctorType * m_functor ;

  void operator () ( const tbb::blocked_range<size_type> &r ) const {
    for(size_type i = r.begin() ; i != r.end(); ++i) {
	  (*m_functor)(i);
	}
  }
  
private:

  ParallelFor(const FunctorType * f ) : m_functor( f )  { }  

public:

  static void execute( const size_type work_count ,
                       const FunctorType & functor )
  {
    DeviceTBB::set_dispatch_functor();

    ParallelFor driver( & functor );

    DeviceTBB::clear_dispatch_functor();

	tbb::parallel_for(tbb::blocked_range<size_type>(0,work_count) , driver , tbb::auto_partitioner() );
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_DEVICETBB_PARALLELFOR_HPP */

