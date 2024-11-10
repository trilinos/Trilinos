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
#ifndef KOKKOSSPARSE_IMPL_TRSV_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_TRSV_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_BsrMatrix.hpp"

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosSparse_trsv_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class CrsMatrixType, class DomainMultiVectorType, class RangeMultiVectorType>
struct trsv_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_TRSV_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,        \
                                         MEM_SPACE_TYPE)                                                              \
  template <>                                                                                                         \
  struct trsv_eti_spec_avail<                                                                                         \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                      \
    enum : bool { value = true };                                                                                     \
  };                                                                                                                  \
                                                                                                                      \
  template <>                                                                                                         \
  struct trsv_eti_spec_avail<                                                                                         \
      KokkosSparse::Experimental::BsrMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                                    \
                                            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                          \
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,              \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                      \
    enum : bool { value = true };                                                                                     \
  };

// Include the actual specialization declarations
#include <KokkosSparse_trsv_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_trsv_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::trsv (sparse matrix - dense
///   vector multiply) for single vectors (1-D Views).

template <class CrsMatrixType, class DomainMultiVectorType, class RangeMultiVectorType,
          bool tpl_spec_avail = trsv_tpl_spec_avail<CrsMatrixType, DomainMultiVectorType, RangeMultiVectorType>::value,
          bool eti_spec_avail = trsv_eti_spec_avail<CrsMatrixType, DomainMultiVectorType, RangeMultiVectorType>::value>
struct TRSV {
  static void trsv(const char uplo[], const char trans[], const char diag[], const CrsMatrixType &A,
                   DomainMultiVectorType B, RangeMultiVectorType X);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of trsv for multi vectors.
// Unification layer
template <class CrsMatrixType, class DomainMultiVectorType, class RangeMultiVectorType>
struct TRSV<CrsMatrixType, DomainMultiVectorType, RangeMultiVectorType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void trsv(const char uplo[], const char trans[], const char diag[], const CrsMatrixType &A,
                   DomainMultiVectorType B,
                   RangeMultiVectorType X)  // X is the output MV
  {
    using Wrap = Sequential::TrsvWrap<CrsMatrixType, DomainMultiVectorType, RangeMultiVectorType>;
    if (trans[0] == 'N' || trans[0] == 'n') {    // no transpose
      if (uplo[0] == 'L' || uplo[0] == 'l') {    // lower triangular
        if (diag[0] == 'U' || diag[0] == 'u') {  // unit diagonal
          Wrap::lowerTriSolveCsrUnitDiag(X, A, B);
        } else {  // non unit diagonal
          Wrap::lowerTriSolveCsr(X, A, B);
        }
      } else {                                   // upper triangular
        if (diag[0] == 'U' || diag[0] == 'u') {  // unit diagonal
          Wrap::upperTriSolveCsrUnitDiag(X, A, B);
        } else {  // non unit diagonal
          Wrap::upperTriSolveCsr(X, A, B);
        }
      }
    } else if (trans[0] == 'T' || trans[0] == 't') {  // transpose
      if (uplo[0] == 'L' || uplo[0] == 'l') {         // lower triangular
        // Transposed lower tri CSR => upper tri CSC.
        if (diag[0] == 'U' || diag[0] == 'u') {  // unit diagonal
          Wrap::upperTriSolveCscUnitDiag(X, A, B);
        } else {  // non unit diagonal
          Wrap::upperTriSolveCsc(X, A, B);
        }
      } else {  // upper triangular
        // Transposed upper tri CSR => lower tri CSC.
        if (diag[0] == 'U' || diag[0] == 'u') {  // unit diagonal
          Wrap::lowerTriSolveCscUnitDiag(X, A, B);
        } else {  // non unit diagonal
          Wrap::lowerTriSolveCsc(X, A, B);
        }
      }
    } else if (trans[0] == 'C' || trans[0] == 'c') {  // conj transpose
      if (uplo[0] == 'L' || uplo[0] == 'l') {         // lower triangular
        // Transposed lower tri CSR => upper tri CSC.
        if (diag[0] == 'U' || diag[0] == 'u') {  // unit diagonal
          Wrap::upperTriSolveCscUnitDiagConj(X, A, B);
        } else {  // non unit diagonal
          Wrap::upperTriSolveCscConj(X, A, B);
        }
      } else {  // upper triangular
        // Transposed upper tri CSR => lower tri CSC.
        if (diag[0] == 'U' || diag[0] == 'u') {  // unit diagonal
          Wrap::lowerTriSolveCscUnitDiagConj(X, A, B);
        } else {  // non unit diagonal
          Wrap::lowerTriSolveCscConj(X, A, B);
        }
      }
    }
  }
};

#endif
}  // namespace Impl
}  // namespace KokkosSparse

//
// Macro for declaration of full specialization of
// KokkosSparse::Impl::Dot for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSSPARSE_TRSV_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,         \
                                        MEM_SPACE_TYPE)                                                               \
  extern template struct TRSV<                                                                                        \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      false, true>;                                                                                                   \
                                                                                                                      \
  extern template struct TRSV<                                                                                        \
      KokkosSparse::Experimental::BsrMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                                    \
                                            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                          \
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,              \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      false, true>;

#define KOKKOSSPARSE_TRSV_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,         \
                                        MEM_SPACE_TYPE)                                                               \
  template struct TRSV<                                                                                               \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      false, true>;                                                                                                   \
                                                                                                                      \
  template struct TRSV<                                                                                               \
      KokkosSparse::Experimental::BsrMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                                    \
                                            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                          \
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,              \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                                  \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      false, true>;

#include <KokkosSparse_trsv_tpl_spec_decl.hpp>

#endif  // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
