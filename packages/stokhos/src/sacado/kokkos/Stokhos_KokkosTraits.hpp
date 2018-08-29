// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_KOKKOS_TRAITS_HPP
#define STOKHOS_KOKKOS_TRAITS_HPP

#include "Sacado_Traits.hpp"

#include "Kokkos_Macros.hpp"
#include "Kokkos_Core_fwd.hpp"

namespace Sacado {

#ifdef KOKKOS_ENABLE_SERIAL
  template <>
  struct StringName< Kokkos::Serial > {
    static std::string eval() { return "Kokkos::Serial"; }
  };
#endif

#ifdef KOKKOS_ENABLE_THREADS
  template <>
  struct StringName< Kokkos::Threads > {
    static std::string eval() { return "Kokkos::Threads"; }
  };
#endif

#ifdef KOKKOS_ENABLE_OPENMP
  template <>
  struct StringName< Kokkos::OpenMP > {
    static std::string eval() { return "Kokkos::OpenMP"; }
  };
#endif

#ifdef KOKKOS_ENABLE_CUDA
  template <>
  struct StringName< Kokkos::Cuda > {
    static std::string eval() { return "Kokkos::Cuda"; }
  };
#endif

}

#endif // STOKHOS_KOKKOS_TRAITS_HPP
