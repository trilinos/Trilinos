/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2015) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef MORKON_KOKKOS_VIEW_TYPES_H
#define MORKON_KOKKOS_VIEW_TYPES_H

/*
   Determine default device type and necessary includes

   Compiled device type is mutually exclusive, can be MORKON_CUDA, MORKON_THREADS, MORKON_OPENMP, or nothing (SERIAL).
   Kokkos library supports KOKKOS_HAVE_CUDA, KOKKOS_HAVE_PTHREAD, KOKKOS_HAVE_OPENMP, or nothing (serial)

*/

#include <KokkosCore_config.h>

// All the below may or may not use HWLOC
#if defined( KOKKOS_HAVE_HWLOC )
#    include <Kokkos_hwloc.hpp>
#endif

#include <Kokkos_Serial.hpp>

#if defined( KOKKOS_HAVE_CUDA ) && defined(MORKON_CUDA)

#ifdef _OPENMP
  #include <Kokkos_OpenMP.hpp>
#else
  #include <Kokkos_Threads.hpp>
#endif

#   include <Kokkos_Cuda.hpp>
#   include <cuda.h>
#   include <cuda_runtime.h>
typedef typename Kokkos::Cuda default_kokkos_device_t;

#elif defined( KOKKOS_HAVE_OPENMP) && defined(MORKON_OPENMP)

#   include <Kokkos_OpenMP.hpp>
typedef typename Kokkos::OpenMP default_kokkos_device_t;

#elif defined( KOKKOS_HAVE_PTHREAD ) && defined(MORKON_THREADS)

#   include <Kokkos_Threads.hpp>
typedef typename Kokkos::Threads default_kokkos_device_t;

#else  // Serial || MPI only

typedef typename Kokkos::Serial default_kokkos_device_t;

#endif

// Different algorithm traits present in MORKON
struct morkon_base_impl {

};

#include <Kokkos_View.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Timer.hpp>

#endif
