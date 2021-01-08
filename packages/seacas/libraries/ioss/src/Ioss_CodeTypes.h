// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_code_types_h
#define IOSS_code_types_h

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#ifdef _WIN64
#define ssize_t __int64
#else
#define ssize_t long
#endif
#endif

namespace Ioss {
  using IntVector   = std::vector<int>;
  using Int64Vector = std::vector<int64_t>;
  using NameList    = std::vector<std::string>;
  using IJK_t       = std::array<int, 3>;
} // namespace Ioss

inline const std::string IOSS_SCALAR() { return std::string("scalar"); }
inline const std::string IOSS_VECTOR_2D() { return std::string("vector_2d"); }
inline const std::string IOSS_VECTOR_3D() { return std::string("vector_3d"); }
inline const std::string IOSS_SYM_TENSOR() { return std::string("sym_tensor_33"); }

#if defined(BUILT_IN_SIERRA)
#define SEACAS_HAVE_MPI
/* #undef IOSS_THREADSAFE */
/* #undef SEACAS_HAVE_KOKKOS */
/* #undef SEACAS_HAVE_DATAWAREHOUSE */
#define SEACAS_HAVE_EXODUS
#define SEACAS_HAVE_CGNS
#define SEACAS_HAVE_PAMGEN
#else
#include <SEACASIoss_config.h>
#endif

#if defined(IOSS_THREADSAFE)
#include <mutex>
#endif

#if defined(SEACAS_HAVE_MPI)
#include <mpi.h>
#define PAR_UNUSED(x)
#else
#define PAR_UNUSED(x)                                                                              \
  do {                                                                                             \
    (void)(x);                                                                                     \
  } while (0)

#ifndef MPI_COMM_SELF
#define MPI_COMM_SELF 0
#endif
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
using MPI_Comm = int;
#endif
#endif

#ifdef SEACAS_HAVE_KOKKOS
#include <Kokkos_Core.hpp> // for Kokkos::complex
#endif

#include <complex>
#if defined(FOUR_BYTE_REAL)
//'FOUR_BYTE_REAL' is a sierra macro which may or may not be defined
// when this header is compiled...
// If FOUR_BYTE_REAL is defined then we know we need float, otherwise
// stick with double.
using Complex = std::complex<float>;
#ifdef SEACAS_HAVE_KOKKOS
using Kokkos_Complex = Kokkos::complex<float>;
#endif
#else
using Complex = std::complex<double>;
#ifdef SEACAS_HAVE_KOKKOS
using Kokkos_Complex = Kokkos::complex<double>;
#endif
#endif
#endif

#if defined(IOSS_THREADSAFE)
#define IOSS_FUNC_ENTER(m) std::lock_guard<std::mutex> guard(m)
#else

#if defined IOSS_TRACE
#include <Ioss_Tracer.h>
#define IOSS_FUNC_ENTER(m) Ioss::Tracer m(__func__)
#else
#define IOSS_FUNC_ENTER(m)
#endif
#endif

#ifndef IOSS_DEBUG_OUTPUT
#define IOSS_DEBUG_OUTPUT 0
#endif
