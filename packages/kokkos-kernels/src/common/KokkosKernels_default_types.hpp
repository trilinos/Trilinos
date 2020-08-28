/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef KOKKOSKERNELS_DEFAULT_TYPES_H
#define KOKKOSKERNELS_DEFAULT_TYPES_H

#include "Kokkos_Core.hpp"        //for LayoutLeft/LayoutRight
#include <KokkosKernels_config.h> //for all the ETI #cmakedefine macros

#if defined(KOKKOSKERNELS_INST_ORDINAL_INT)
  typedef int default_lno_t;
#elif defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T)
  typedef int64_t default_lno_t;
#else
  #error "Expect INT and/or INT64_T to be enabled as ORDINAL (lno_t) types"
#endif
  //Prefer int as the default offset type, because cuSPARSE doesn't support size_t for rowptrs.
#if defined(KOKKOSKERNELS_INST_OFFSET_INT)
  typedef int default_size_type;
#elif defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  typedef size_t default_size_type;
#else
  #error "Expect SIZE_T and/or INT to be enabled as OFFSET (size_type) types"
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  typedef Kokkos::LayoutLeft default_layout;
#elif defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  typedef Kokkos::LayoutRight default_layout;
#else
  #error "Expect LAYOUTLEFT and/or LAYOUTRIGHT to be enabled as layout types"
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
  typedef double default_scalar;
#elif defined(KOKKOSKERNELS_INST_FLOAT)
  typedef float default_scalar;
#else
  #error "Expect at least one real-valued scalar type (double or float) to be enabled"
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  typedef Kokkos::Cuda default_device;
#elif defined(KOKKOS_ENABLE_OPENMP)
  typedef Kokkos::OpenMP default_device;
#elif defined(KOKKOS_ENABLE_PTHREAD) || defined(KOKKOS_ENABLE_THREADS)
  typedef Kokkos::Threads default_device;
#else
  typedef Kokkos::Serial default_device;
#endif

#endif // KOKKOSKERNELS_DEFAULT_TYPES_H
