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

#ifndef KOKKOS_HOSTDEVICEFOR_HPP
#define KOKKOS_HOSTDEVICEFOR_HPP

#include <algorithm>
#include <TPI.h>

namespace Kokkos {

class HostDevice ;

template< class FunctorType , class DeviceType > class ParallelFor ;

template< class FunctorType >
void run_functor_on_tpi( TPI_Work * work )
{
  const FunctorType & functor = *((const FunctorType *) work->info );

  const int work_count = functor.work_count();
  const int work_inc   = ( work_count + work->count - 1 ) / work->count ;
  const int work_begin = work_inc * work->rank ;
  const int work_end   = std::max( work_begin + work_inc , work_count );

  for ( int iwork = work_begin ; iwork < work_end ; ++iwork ) {
    functor( iwork );
  }
}

template< class FunctorType >
inline
void parallel_for( const FunctorType & functor )
{
  TPI_Run_threads( & run_functor_on_tpi<FunctorType> , & functor , 0 );
};

}

#endif /* KOKKOS_HOSTDEVICEFOR_HPP */

