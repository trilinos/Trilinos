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

#include "Kokkos_Core.hpp"         //for LayoutLeft/LayoutRight
#include <KokkosKernels_config.h>  //for all the ETI #cmakedefine macros

#if defined(KOKKOSKERNELS_INST_ORDINAL_INT)
using default_lno_t = int;
#elif defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T)
using default_lno_t     = int64_t;
#else
using default_lno_t     = int;
#endif
// Prefer int as the default offset type, because cuSPARSE doesn't support
// size_t for rowptrs.
#if defined(KOKKOSKERNELS_INST_OFFSET_INT)
using default_size_type = int;
#elif defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)
using default_size_type = size_t;
#else
using default_size_type = int;
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
using default_layout = Kokkos::LayoutLeft;
#elif defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
using default_layout    = Kokkos::LayoutRight;
#else
using default_layout    = Kokkos::LayoutLeft;
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
using default_scalar = double;
#elif defined(KOKKOSKERNELS_INST_FLOAT)
using default_scalar    = float;
#elif defined(KOKKOSKERNELS_INST_HALF)
using default_scalar    = Kokkos::Experimental::half_t;
#else
using default_scalar = double;
#endif

#if defined(KOKKOS_ENABLE_CUDA)
using default_device = Kokkos::Cuda;
#elif defined(KOKKOS_ENABLE_HIP)
using default_device    = Kokkos::Experimental::HIP;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
using default_device    = Kokkos::Experimental::OpenMPTarget;
#elif defined(KOKKOS_ENABLE_OPENMP)
using default_device = Kokkos::OpenMP;
#elif defined(KOKKOS_ENABLE_PTHREAD) || defined(KOKKOS_ENABLE_THREADS)
using default_device = Kokkos::Threads;
#else
using default_device = Kokkos::Serial;
#endif

#endif  // KOKKOSKERNELS_DEFAULT_TYPES_H
