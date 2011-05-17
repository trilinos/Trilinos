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

#ifndef KOKKOS_DEVICETPI_PARALLELFOR_HPP
#define KOKKOS_DEVICETPI_PARALLELFOR_HPP

#include <algorithm>
#include <TPI.h>

namespace Kokkos {

template< class FunctorType >
class ParallelFor< FunctorType , DeviceTPI > {
public:
  typedef DeviceTPI              device_type ;
  typedef device_type::size_type size_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  ParallelFor( const size_type work_count , const FunctorType & functor )
    : m_work_functor( functor )
    , m_work_count( work_count )
    {}

private:

  // self.m_work_count == total work count
  // work->count       == number of threads

  static void run_on_tpi( TPI_Work * work )
  {
    const ParallelFor & self = *((const ParallelFor *) work->info );

    const size_type work_inc   = (self.m_work_count + work->count - 1) / work->count ;
    const size_type work_begin = work_inc * work->rank ;
    const size_type work_end   = std::min( work_begin + work_inc , self.m_work_count );

    for ( size_type iwork = work_begin ; iwork < work_end ; ++iwork ) {
      self.m_work_functor( iwork );
    }
  }

public:

  static void run( const size_type     work_count ,
                   const FunctorType & work_functor )
  {
    // Make a copy just like other devices will have to.

    device_type::set_dispatch_functor();

    const ParallelFor tmp( work_count , work_functor );

    device_type::clear_dispatch_functor();

    TPI_Run_threads( & run_on_tpi , & tmp , 0 );
  }
};

} // namespace Kokkos

#endif /* KOKKOS_DEVICETPI_PARALLELFOR_HPP */

