/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CORE_HPP
#define KOKKOS_CORE_HPP

//----------------------------------------------------------------------------
// Include the execution space header files for the enabled execution spaces.

#include <Kokkos_Core_fwd.hpp>

#if defined( KOKKOS_HAVE_CUDA )
#include <Kokkos_Cuda.hpp>
#endif

#if defined( KOKKOS_HAVE_OPENMP )
#include <Kokkos_OpenMP.hpp>
#endif

#if defined( KOKKOS_HAVE_SERIAL )
#include <Kokkos_Serial.hpp>
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
#include <Kokkos_Threads.hpp>
#endif

#include <Kokkos_Pair.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Vectorization.hpp>
#include <Kokkos_Atomic.hpp>
#include <Kokkos_hwloc.hpp>

#ifdef KOKKOS_HAVE_CXX11
#include <Kokkos_Complex.hpp>
#endif

//----------------------------------------------------------------------------

namespace Kokkos {

struct InitArguments {
  int num_threads;
  int num_numa;
  int device_id;

  InitArguments() {
    num_threads = -1;
    num_numa = -1;
    device_id = -1;
  }
};

void initialize(int& narg, char* arg[]);

void initialize(const InitArguments& args = InitArguments());

/** \brief  Finalize the spaces that were initialized via Kokkos::initialize */
void finalize();

/** \brief  Finalize all known execution spaces */
void finalize_all();

void fence();

}

#ifdef KOKKOS_HAVE_CXX11
namespace Kokkos {
/* Allocate memory from a memory space.
 * The allocation is tracked in Kokkos memory tracking system, so
 * leaked memory can be identified.
 */
template< class Arg = DefaultExecutionSpace>
void* malloc(const std::string label, size_t count) {
  typedef typename Arg::memory_space MemorySpace;
  return MemorySpace::allocate(label,count);;
}

template< class Arg = DefaultExecutionSpace>
void* malloc(const size_t& count) {
  return malloc<Arg>("DefaultLabel",count);
}

/* Free memory from a memory space.
 */
template< class Arg = DefaultExecutionSpace>
void free(const void* ptr) {
  typedef typename Arg::memory_space MemorySpace;
  MemorySpace::decrement(ptr);
}

template< class Arg = DefaultExecutionSpace>
const void* realloc(const void* old_ptr, size_t size) {
  typedef typename Arg::memory_space MemorySpace;

  //Get information about the old allocation
  const void* start_ptr = MemorySpace::query_start_ptr(old_ptr);
  const size_t old_size = MemorySpace::query_size(old_ptr);
  const std::string label = MemorySpace::query_label(old_ptr);

  if (start_ptr != old_ptr) Impl::throw_runtime_exception("Calling Kokkos::realloc<MemorySpace> is only allowed with pointers which are obtained by a call to Kokkos::malloc<MemorySpace>");

  if (old_size == size) return old_ptr;

  //Do the new allocation
  void* new_ptr = malloc<MemorySpace>(label,size);

  //Copy old data to the new allocation
  Impl::DeepCopy<MemorySpace,MemorySpace>(new_ptr,old_ptr,size>old_size?old_size:size);

  //Elliminate the old allocation
  free<MemorySpace>(old_ptr);

  return new_ptr;
}
}
#endif

#endif
