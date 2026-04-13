// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSSPARSE_SPMV_SELLMATRIX_SPEC_HPP_
#define KOKKOSSPARSE_SPMV_SELLMATRIX_SPEC_HPP_

#include "KokkosSparse_SellMatrix.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosSparse_spmv_sellmatrix_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class XVector, class YVector>
struct spmv_sellmatrix_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPMV_SELLMATRIX_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,  \
                                                    EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                      \
  template <>                                                                                             \
  struct spmv_sellmatrix_eti_spec_avail<                                                                  \
      EXEC_SPACE_TYPE,                                                                                    \
      KokkosSparse::Experimental::SellMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                       \
                                             Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                                             Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>, \
      Kokkos::View<SCALAR_TYPE const*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                       \
      Kokkos::View<SCALAR_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                            \
    enum : bool { value = true };                                                                         \
  };

// Include the actual specialization declarations
#include <KokkosSparse_spmv_sellmatrix_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spmv_sellmatrix_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::spmv (sparse matrix - dense
///   vector multiply) for single vectors (1-D Views).
template <
    class ExecutionSpace, class AMatrix, class XVector, class YVector,
    bool tpl_sellmatrix_spec_avail = spmv_sellmatrix_tpl_spec_avail<ExecutionSpace, AMatrix, XVector, YVector>::value,
    bool eti_sellmatrix_spec_avail = spmv_sellmatrix_eti_spec_avail<ExecutionSpace, AMatrix, XVector, YVector>::value>
struct SPMV_SELLMATRIX {
  using coefficient_type = typename YVector::non_const_value_type;

  static void spmv(const ExecutionSpace& space, const char mode[], const coefficient_type& alpha, const AMatrix& A,
                   const XVector& x, const coefficient_type& beta, const YVector& y);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of spmv for single vectors (1-D Views).
// Unification layer
template <class ExecutionSpace, class AMatrix, class XVector, class YVector>
struct SPMV_SELLMATRIX<ExecutionSpace, AMatrix, XVector, YVector, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  using coefficient_type = typename YVector::non_const_value_type;

  static void spmv(const ExecutionSpace& space, const char mode[], const coefficient_type& alpha, const AMatrix& A,
                   const XVector& x, const coefficient_type& beta, const YVector& y) {
    spmv_sellmatrix_impl<ExecutionSpace, AMatrix, XVector, YVector>(space, mode, alpha, A, x, beta, y);
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosSparse

//
// Macro for declaration of full specialization of
// KokkosSparse::Impl::SPMV_SELLMATRIX.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSSPARSE_SPMV_SELLMATRIX_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,     \
                                                   EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                         \
  extern template struct SPMV_SELLMATRIX<                                                                   \
      EXEC_SPACE_TYPE,                                                                                      \
      ::KokkosSparse::Experimental::SellMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                       \
                                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>, \
      Kokkos::View<SCALAR_TYPE const*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                         \
      Kokkos::View<SCALAR_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
      false, true>;

#include <generated_specializations_hpp/KokkosSparse_spmv_sellmatrix_eti_spec_decl.hpp>

#define KOKKOSSPARSE_SPMV_SELLMATRIX_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,     \
                                                   EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                         \
  template struct SPMV_SELLMATRIX<                                                                          \
      EXEC_SPACE_TYPE,                                                                                      \
      ::KokkosSparse::Experimental::SellMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                       \
                                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>, \
      Kokkos::View<SCALAR_TYPE const*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                         \
      Kokkos::View<SCALAR_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
      false, true>;

#include <KokkosSparse_spmv_sellmatrix_tpl_spec_decl.hpp>

#endif  // KOKKOSSPARSE_SPMV_SELLMATRIX_SPEC_HPP_
