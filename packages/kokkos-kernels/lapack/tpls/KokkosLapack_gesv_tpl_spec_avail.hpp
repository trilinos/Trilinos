//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_HPP_
#define KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_HPP_

namespace KokkosLapack {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class AMatrix, class BXMV>
struct gesv_tpl_spec_avail {
  enum : bool { value = false };
};

// Generic Host side LAPACK (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK

#define KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_LAPACK(SCALAR, LAYOUT, MEMSPACE) \
  template <class ExecSpace>                                              \
  struct gesv_tpl_spec_avail<                                             \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {          \
    enum : bool { value = true };                                         \
  };

KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_LAPACK(double, Kokkos::LayoutLeft,
                                        Kokkos::HostSpace)
KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_LAPACK(float, Kokkos::LayoutLeft,
                                        Kokkos::HostSpace)
KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<double>,
                                        Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<float>,
                                        Kokkos::LayoutLeft, Kokkos::HostSpace)
/*
#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_LAPACK( double, Kokkos::LayoutRight,
Kokkos::HostSpace) #endif
#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_LAPACK( float, Kokkos::LayoutRight,
Kokkos::HostSpace) #endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_LAPACK( Kokkos::complex<double>,
Kokkos::LayoutRight, Kokkos::HostSpace) #endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_LAPACK( Kokkos::complex<float>,
Kokkos::LayoutRight, Kokkos::HostSpace) #endif
*/
#endif

// MAGMA
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA

#define KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_MAGMA(SCALAR, LAYOUT, MEMSPACE)  \
  template <class ExecSpace>                                              \
  struct gesv_tpl_spec_avail<                                             \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {          \
    enum : bool { value = true };                                         \
  };

KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_MAGMA(double, Kokkos::LayoutLeft,
                                       Kokkos::CudaSpace)
KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_MAGMA(float, Kokkos::LayoutLeft,
                                       Kokkos::CudaSpace)
KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<double>,
                                       Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<float>,
                                       Kokkos::LayoutLeft, Kokkos::CudaSpace)

/*
#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_MAGMA( double, Kokkos::LayoutRight,
Kokkos::CudaSpace) #endif
#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_MAGMA( float, Kokkos::LayoutRight,
Kokkos::CudaSpace) #endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_MAGMA(
Kokkos::complex<double>,Kokkos::LayoutRight, Kokkos::CudaSpace) #endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSLAPACK_GESV_TPL_SPEC_AVAIL_MAGMA( Kokkos::complex<float>,
Kokkos::LayoutRight, Kokkos::CudaSpace) #endif
*/
#endif

}  // namespace Impl
}  // namespace KokkosLapack

#endif
