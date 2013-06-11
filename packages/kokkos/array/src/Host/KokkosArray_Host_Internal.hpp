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

#ifndef KOKKOSARRAY_HOST_INTERNAL_HPP
#define KOKKOSARRAY_HOST_INTERNAL_HPP

#include <utility>

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

/** \brief  The driver subprogram to be run by a spawned thread. */
void host_thread_driver();

/** \name  Internal interface to threading runtime. */
/**@{ */

/** \brief  Span a thread and call 'host_thread_driver()' */
bool host_thread_spawn();

/** \brief  Query if called on the master thread */
bool host_thread_is_master();

/** \brief  Wait for *flag != value */
int  host_thread_wait( volatile int * const flag , const int value );

/** \brief  Hard lock the current thread (e.g., via pthread_mutex */
void host_thread_lock();

/** \brief  Unlock the current thread */
void host_thread_unlock();
/**@} */

//----------------------------------------------------------------------------

void host_thread_mapping( const std::pair<unsigned,unsigned> gang_topo ,
                          const std::pair<unsigned,unsigned> core_use ,
                          const std::pair<unsigned,unsigned> core_topo ,
                          const std::pair<unsigned,unsigned> master_coord ,
                                std::pair<unsigned,unsigned> thread_coord[] );

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_HOST_INTERNAL_HPP */

