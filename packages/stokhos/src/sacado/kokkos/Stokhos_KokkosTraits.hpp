// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
