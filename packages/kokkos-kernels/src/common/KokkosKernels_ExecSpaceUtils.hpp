/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include "Kokkos_Core.hpp"
#include "Kokkos_Atomic.hpp"

#ifndef _KOKKOSKERNELSUTILSEXECSPACEUTILS_HPP
#define _KOKKOSKERNELSUTILSEXECSPACEUTILS_HPP


namespace KokkosKernels{

namespace Impl{

enum ExecSpaceType{Exec_SERIAL, Exec_OMP, Exec_PTHREADS, Exec_QTHREADS, Exec_CUDA};
template <typename ExecutionSpace>
inline ExecSpaceType kk_get_exec_space_type(){
  ExecSpaceType exec_space = Exec_SERIAL;
#if defined( KOKKOS_ENABLE_SERIAL )
  if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
    exec_space = Exec_SERIAL;
  }
#endif

#if defined( KOKKOS_ENABLE_THREADS )
  if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
    exec_space =  Exec_PTHREADS;
  }
#endif

#if defined( KOKKOS_ENABLE_OPENMP )
  if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
    exec_space = Exec_OMP;
  }
#endif

#if defined( KOKKOS_ENABLE_CUDA )
  if (Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value){
    exec_space = Exec_CUDA;
  }
#endif

#if defined( KOKKOS_ENABLE_QTHREAD)
  if (Kokkos::Impl::is_same< Kokkos::Qthread, ExecutionSpace >::value){
    exec_space = Exec_QTHREADS;
  }
#endif
  return exec_space;

}


inline int kk_get_suggested_vector_size(
    const size_t nr, const  size_t nnz, const ExecSpaceType exec_space){
  int suggested_vector_size_ = 1;
  switch (exec_space){
  default:
    break;
  case Exec_SERIAL:
  case Exec_OMP:
  case Exec_PTHREADS:
  case Exec_QTHREADS:
    break;
  case Exec_CUDA:

    if (nr > 0)
      suggested_vector_size_ = nnz / double (nr) + 0.5;
    if (suggested_vector_size_ < 3){
      suggested_vector_size_ = 2;
    }
    else if (suggested_vector_size_ <= 6){
      suggested_vector_size_ = 4;
    }
    else if (suggested_vector_size_ <= 12){
      suggested_vector_size_ = 8;
    }
    else if (suggested_vector_size_ <= 24){
      suggested_vector_size_ = 16;
    }
    else {
      suggested_vector_size_ = 32;
    }
    break;
  }
  return suggested_vector_size_;

}


inline int kk_get_suggested_team_size(const int vector_size, const ExecSpaceType exec_space){
  if (exec_space == Exec_CUDA){
    return 256 / vector_size;
  }
  else {
    return 1;
  }
}

}
}

#endif
