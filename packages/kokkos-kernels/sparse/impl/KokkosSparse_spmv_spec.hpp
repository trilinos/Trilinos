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
#ifndef KOKKOSSPARSE_IMPL_SPMV_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPMV_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv_handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosSparse_spmv_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct spmv_eti_spec_avail {
  enum : bool { value = false };
};
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector,
          const bool integerScalarType = std::is_integral_v<typename AMatrix::non_const_value_type>>
struct spmv_mv_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPMV_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,        \
                                         MEM_SPACE_TYPE)                                                              \
  template <>                                                                                                         \
  struct spmv_eti_spec_avail<                                                                                         \
      EXEC_SPACE_TYPE,                                                                                                \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>,    \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      Kokkos::View<SCALAR_TYPE const*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                   \
      Kokkos::View<SCALAR_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                        \
    enum : bool { value = true };                                                                                     \
  };

#define KOKKOSSPARSE_SPMV_MV_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,     \
                                            MEM_SPACE_TYPE)                                                           \
  template <>                                                                                                         \
  struct spmv_mv_eti_spec_avail<                                                                                      \
      EXEC_SPACE_TYPE,                                                                                                \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>,    \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      Kokkos::View<SCALAR_TYPE const**, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                   \
      Kokkos::View<SCALAR_TYPE**, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                        \
    enum : bool { value = true };                                                                                     \
  };

// Include the actual specialization declarations
#include <KokkosSparse_spmv_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spmv_eti_spec_avail.hpp>

#include <KokkosSparse_spmv_mv_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spmv_mv_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::spmv (sparse matrix - dense
///   vector multiply) for single vectors (1-D Views).
///
/// For the implementation of KokkosSparse::spmv for multivectors (2-D
/// Views), see the SPMV_MV struct below.
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector,
          bool tpl_spec_avail = spmv_tpl_spec_avail<ExecutionSpace, Handle, AMatrix, XVector, YVector>::value,
          bool eti_spec_avail = spmv_eti_spec_avail<ExecutionSpace, Handle, AMatrix, XVector, YVector>::value>
struct SPMV {
  typedef typename YVector::non_const_value_type coefficient_type;

  static void spmv(const ExecutionSpace& space, Handle* handle, const char mode[], const coefficient_type& alpha,
                   const AMatrix& A, const XVector& x, const coefficient_type& beta, const YVector& y);
};

// Unification layer
/// \brief Implementation of KokkosBlas::spmv (sparse matrix - dense
///   vector multiply) for multiple vectors at a time (multivectors)
///   and possibly multiple coefficients at a time.
///
/// This struct implements the following operations:
///
///   1. Y(:,j) := beta(j) * Y(:,j) + alpha(j) * Op(A) * X(:,j)
///   2. Y(:,j) := beta(j) * Y(:,j) + alpha * Op(A) * X(:,j)
///   3. Y(:,j) := beta * Y(:,j) + alpha(j) * Op(A) * X(:,j)
///   4. Y(:,j) := beta * Y(:,j) + alpha * Op(A) * X(:,j)
///
/// In #1 and #2 above, beta is a 1-D View of coefficients, one for
/// each column of Y.  In #1 and #3 above, alpha is a 1-D View of
/// coefficients, one for each column of X.  Otherwise, alpha
/// resp. beta are each a single coefficient.  In all of these
/// operations, X and Y are 2-D Views ("multivectors").  A is a sparse
/// matrix, and Op(A) is either A itself, its transpose, or its
/// conjugate transpose, depending on the 'mode' argument.
///
/// The last template parameter (integerScalarType) indicates whether the
/// matrix's entries have integer type.  Per Github Issue #700, we
/// don't optimize as heavily for that case, in order to reduce build
/// times and library sizes.
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector,
          const bool integerScalarType = std::is_integral_v<typename AMatrix::non_const_value_type>,
          bool tpl_spec_avail = spmv_mv_tpl_spec_avail<ExecutionSpace, Handle, AMatrix, XVector, YVector>::value,
          bool eti_spec_avail = spmv_mv_eti_spec_avail<ExecutionSpace, Handle, AMatrix, XVector, YVector>::value>
struct SPMV_MV {
  typedef typename YVector::non_const_value_type coefficient_type;

  static void spmv_mv(const ExecutionSpace& space, Handle* handle, const char mode[], const coefficient_type& alpha,
                      const AMatrix& A, const XVector& x, const coefficient_type& beta, const YVector& y);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of spmv for single vectors (1-D Views).
// Unification layer
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct SPMV<ExecutionSpace, Handle, AMatrix, XVector, YVector, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YVector::non_const_value_type coefficient_type;

  static void spmv(const ExecutionSpace& space, Handle* handle, const char mode[], const coefficient_type& alpha,
                   const AMatrix& A, const XVector& x, const coefficient_type& beta, const YVector& y) {
    typedef Kokkos::ArithTraits<coefficient_type> KAT;

    if (beta == KAT::zero()) {
      spmv_beta<ExecutionSpace, Handle, AMatrix, XVector, YVector, 0>(space, handle, mode, alpha, A, x, beta, y);
    } else if (beta == KAT::one()) {
      spmv_beta<ExecutionSpace, Handle, AMatrix, XVector, YVector, 1>(space, handle, mode, alpha, A, x, beta, y);
    } else if (beta == -KAT::one()) {
      spmv_beta<ExecutionSpace, Handle, AMatrix, XVector, YVector, -1>(space, handle, mode, alpha, A, x, beta, y);
    } else {
      spmv_beta<ExecutionSpace, Handle, AMatrix, XVector, YVector, 2>(space, handle, mode, alpha, A, x, beta, y);
    }
  }
};

//! Full specialization of spmv_mv for single vectors (2-D Views).
// Unification layer
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct SPMV_MV<ExecutionSpace, Handle, AMatrix, XVector, YVector, false, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YVector::non_const_value_type coefficient_type;

  // TODO: pass handle through to implementation and use tuning parameters
  static void spmv_mv(const ExecutionSpace& space, Handle* /* handle */, const char mode[],
                      const coefficient_type& alpha, const AMatrix& A, const XVector& x, const coefficient_type& beta,
                      const YVector& y) {
    typedef Kokkos::ArithTraits<coefficient_type> KAT;
    if (alpha == KAT::zero()) {
      spmv_alpha_mv<ExecutionSpace, AMatrix, XVector, YVector, 0>(space, mode, alpha, A, x, beta, y);
    } else if (alpha == KAT::one()) {
      spmv_alpha_mv<ExecutionSpace, AMatrix, XVector, YVector, 1>(space, mode, alpha, A, x, beta, y);
    } else if (alpha == -KAT::one()) {
      spmv_alpha_mv<ExecutionSpace, AMatrix, XVector, YVector, -1>(space, mode, alpha, A, x, beta, y);
    } else {
      spmv_alpha_mv<ExecutionSpace, AMatrix, XVector, YVector, 2>(space, mode, alpha, A, x, beta, y);
    }
  }
};

template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct SPMV_MV<ExecutionSpace, Handle, AMatrix, XVector, YVector, true, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YVector::non_const_value_type coefficient_type;

  static void spmv_mv(const ExecutionSpace& space, Handle* handle, const char mode[], const coefficient_type& alpha,
                      const AMatrix& A, const XVector& x, const coefficient_type& beta, const YVector& y) {
    static_assert(std::is_integral_v<typename AMatrix::non_const_value_type>,
                  "This implementation is only for integer Scalar types.");
    KokkosKernels::Experimental::Controls defaultControls;
    for (size_t j = 0; j < x.extent(1); ++j) {
      auto x_j = Kokkos::subview(x, Kokkos::ALL(), j);
      auto y_j = Kokkos::subview(y, Kokkos::ALL(), j);
      typedef SPMV<ExecutionSpace, Handle, AMatrix, decltype(x_j), decltype(y_j)> impl_type;
      impl_type::spmv(space, handle, mode, alpha, A, x_j, beta, y_j);
    }
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosSparse

//
// Macro for declaration of full specialization of
// KokkosSparse::Impl::SpMV.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,         \
                                        MEM_SPACE_TYPE)                                                               \
  extern template struct SPMV<                                                                                        \
      EXEC_SPACE_TYPE,                                                                                                \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>,    \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      Kokkos::View<SCALAR_TYPE const*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                   \
      Kokkos::View<SCALAR_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      false, true>;

#define KOKKOSSPARSE_SPMV_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,         \
                                        MEM_SPACE_TYPE)                                                               \
  template struct SPMV<                                                                                               \
      EXEC_SPACE_TYPE,                                                                                                \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>,    \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      Kokkos::View<SCALAR_TYPE const*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                   \
      Kokkos::View<SCALAR_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      false, true>;

#define KOKKOSSPARSE_SPMV_MV_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,      \
                                           MEM_SPACE_TYPE)                                                            \
  extern template struct SPMV_MV<                                                                                     \
      EXEC_SPACE_TYPE,                                                                                                \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>,    \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      Kokkos::View<SCALAR_TYPE const**, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                   \
      Kokkos::View<SCALAR_TYPE**, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      std::is_integral_v<SCALAR_TYPE>, false, true>;

#define KOKKOSSPARSE_SPMV_MV_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,      \
                                           MEM_SPACE_TYPE)                                                            \
  template struct SPMV_MV<                                                                                            \
      EXEC_SPACE_TYPE,                                                                                                \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>,    \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      Kokkos::View<SCALAR_TYPE const**, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                   \
      Kokkos::View<SCALAR_TYPE**, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      std::is_integral_v<SCALAR_TYPE>, false, true>;

#include <KokkosSparse_spmv_tpl_spec_decl.hpp>

#include <KokkosSparse_spmv_mv_tpl_spec_decl.hpp>

#endif  // KOKKOSSPARSE_IMPL_SPMV_SPEC_HPP_
