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
#ifndef KOKKOSSPARSE_IMPL_SPMV_BSRMATRIX_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPMV_BSRMATRIX_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

#include "KokkosSparse_BsrMatrix.hpp"
#include "KokkosSparse_spmv_handle.hpp"
#include "KokkosKernels_Error.hpp"
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosSparse_spmv_bsrmatrix_impl.hpp>
#include "KokkosSparse_spmv_bsrmatrix_impl_v42.hpp"
#endif

namespace KokkosSparse {
namespace Impl {

// default is no eti available
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct spmv_bsrmatrix_eti_spec_avail {
  enum : bool { value = false };
};

template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector,
          const bool integerScalarType = std::is_integral_v<typename AMatrix::non_const_value_type>>
struct spmv_mv_bsrmatrix_eti_spec_avail {
  enum : bool { value = false };
};

#define KOKKOSSPARSE_SPMV_BSRMATRIX_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,            \
                                                   EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                                \
  template <>                                                                                                      \
  struct spmv_bsrmatrix_eti_spec_avail<                                                                            \
      EXEC_SPACE_TYPE,                                                                                             \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>, \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                               \
                                              Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                     \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,         \
      Kokkos::View<SCALAR_TYPE const *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                     \
    enum : bool { value = true };                                                                                  \
  };

#define KOKKOSSPARSE_SPMV_MV_BSRMATRIX_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,         \
                                                      EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                             \
  template <>                                                                                                      \
  struct spmv_mv_bsrmatrix_eti_spec_avail<                                                                         \
      EXEC_SPACE_TYPE,                                                                                             \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>, \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                               \
                                              Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                     \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,         \
      Kokkos::View<SCALAR_TYPE const **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                     \
    enum : bool { value = true };                                                                                  \
  };

}  // namespace Impl
}  // namespace KokkosSparse

// Include which ETIs are available
#include <KokkosSparse_spmv_bsrmatrix_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spmv_bsrmatrix_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spmv_mv_bsrmatrix_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// declaration
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector,
          bool tpl_spec_avail = spmv_bsrmatrix_tpl_spec_avail<ExecutionSpace, Handle, AMatrix, XVector, YVector>::value,
          bool eti_spec_avail = spmv_bsrmatrix_eti_spec_avail<ExecutionSpace, Handle, AMatrix, XVector, YVector>::value>
struct SPMV_BSRMATRIX {
  typedef typename YVector::non_const_value_type YScalar;

  static void spmv_bsrmatrix(const ExecutionSpace &space, Handle *handle, const char mode[], const YScalar &alpha,
                             const AMatrix &A, const XVector &x, const YScalar &beta, const YVector &y);
};

// declaration
template <
    class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector,
    const bool integerScalarType = std::is_integral_v<typename AMatrix::non_const_value_type>,
    bool tpl_spec_avail = spmv_mv_bsrmatrix_tpl_spec_avail<ExecutionSpace, Handle, AMatrix, XVector, YVector>::value,
    bool eti_spec_avail = spmv_mv_bsrmatrix_eti_spec_avail<ExecutionSpace, Handle, AMatrix, XVector, YVector>::value>
struct SPMV_MV_BSRMATRIX {
  typedef typename YVector::non_const_value_type YScalar;

  static void spmv_mv_bsrmatrix(const ExecutionSpace &space, Handle *handle, const char mode[], const YScalar &alpha,
                                const AMatrix &A, const XVector &x, const YScalar &beta, const YVector &y);
};

// actual implementations to be compiled
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct SPMV_BSRMATRIX<ExecutionSpace, Handle, AMatrix, XVector, YVector, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YVector::non_const_value_type YScalar;

  static void spmv_bsrmatrix(const ExecutionSpace &space, Handle *handle, const char mode[], const YScalar &alpha,
                             const AMatrix &A, const XVector &X, const YScalar &beta, const YVector &Y) {
    const bool modeIsNoTrans        = (mode[0] == NoTranspose[0]);
    const bool modeIsConjugate      = (mode[0] == Conjugate[0]);
    const bool modeIsConjugateTrans = (mode[0] == ConjugateTranspose[0]);
    const bool modeIsTrans          = (mode[0] == Transpose[0]);

    // use V41 if requested
    if (handle->algo == SPMV_BSR_V41) {
      if (modeIsNoTrans || modeIsConjugate) {
        return Bsr::spMatVec_no_transpose(space, handle, alpha, A, X, beta, Y, modeIsConjugate);
      } else if (modeIsTrans || modeIsConjugateTrans) {
        return Bsr::spMatVec_transpose(space, handle, alpha, A, X, beta, Y, modeIsConjugateTrans);
      }
    }

    // use V42 if possible
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>() || handle->algo == SPMV_BSR_V42) {
      if (modeIsNoTrans) {
        ::KokkosSparse::Impl::apply_v42(space, alpha, A, X, beta, Y);
        return;
      }
    }

    // fall back to V41 all else fails
    if (modeIsNoTrans || modeIsConjugate) {
      return Bsr::spMatVec_no_transpose(space, handle, alpha, A, X, beta, Y, modeIsConjugate);
    } else if (modeIsTrans || modeIsConjugateTrans) {
      return Bsr::spMatVec_transpose(space, handle, alpha, A, X, beta, Y, modeIsConjugateTrans);
    }

    {
      std::stringstream ss;
      ss << __FILE__ << ":" << __LINE__ << " ";
      ss << "Internal logic error: no applicable BsrMatrix SpMV implementation "
            ". Please report this";
      throw std::runtime_error(ss.str());
    }
  }
};

template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct SPMV_MV_BSRMATRIX<ExecutionSpace, Handle, AMatrix, XVector, YVector, false, false,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YVector::non_const_value_type YScalar;

  enum class Method {
    Fallback,    ///< Don't use tensor cores
    TensorCores  ///< use tensor cores
  };

  static void spmv_mv_bsrmatrix(const ExecutionSpace &space, Handle *handle, const char mode[], const YScalar &alpha,
                                const AMatrix &A, const XVector &X, const YScalar &beta, const YVector &Y) {
#if defined(KOKKOS_ENABLE_CUDA) && (defined(KOKKOS_ARCH_AMPERE) || defined(KOKKOS_ARCH_VOLTA))
    Method method = Method::Fallback;
    {
      // try to use tensor cores if requested
      if (handle->algo == SPMV_BSR_TC) method = Method::TensorCores;
      if (!KokkosSparse::Impl::TensorCoresAvailable<ExecutionSpace, AMatrix, XVector, YVector>::value) {
        method = Method::Fallback;
      }
      // can't use tensor cores unless mode is no-transpose
      if (mode[0] != KokkosSparse::NoTranspose[0]) method = Method::Fallback;
#if KOKKOS_HALF_T_IS_FLOAT
      // disable tensor cores when Kokkos half is actually a float
      method = Method::Fallback;
#endif  // KOKKOS_HALF_T_IS_FLOAT
    }
#endif  // AMPERE || VOLTA

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ARCH_AMPERE)
    {
      typedef Kokkos::Experimental::half_t Half;
      typedef typename AMatrix::non_const_value_type AScalar;
      typedef typename XVector::non_const_value_type XScalar;

      /* Ampere has double += double * double and float += half * half

      use whichever is requested.
      If none requested, used mixed precision if the inputs are mixed, otherwise
      use double
      */
      if (Method::TensorCores == method) {
        auto precision = handle->bsr_tc_precision;
        switch (precision) {
          case KokkosSparse::Experimental::Bsr_TC_Precision::Mixed: {
            BsrMatrixSpMVTensorCoreDispatcher<ExecutionSpace, AMatrix, half, XVector, half, YVector, float, 16, 16,
                                              16>::dispatch(space, alpha, A, X, beta, Y);
            return;
          }
          case KokkosSparse::Experimental::Bsr_TC_Precision::Double: {
            BsrMatrixSpMVTensorCoreDispatcher<ExecutionSpace, AMatrix, double, XVector, double, YVector, double, 8, 8,
                                              4>::dispatch(space, alpha, A, X, beta, Y);
            return;
          }
          case KokkosSparse::Experimental::Bsr_TC_Precision::Automatic:
          default: {
            constexpr bool operandsHalfHalfFloat = std::is_same<AScalar, Half>::value &&
                                                   std::is_same<XScalar, Half>::value &&
                                                   std::is_same<YScalar, float>::value;
            if (operandsHalfHalfFloat) {
              BsrMatrixSpMVTensorCoreDispatcher<ExecutionSpace, AMatrix, half, XVector, half, YVector, float, 16, 16,
                                                16>::dispatch(space, alpha, A, X, beta, Y);
              return;
            } else {
              BsrMatrixSpMVTensorCoreDispatcher<ExecutionSpace, AMatrix, double, XVector, double, YVector, double, 8, 8,
                                                4>::dispatch(space, alpha, A, X, beta, Y);
              return;
            }
          }
        }
      }
    }
#elif defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ARCH_VOLTA)
    {
      /* Volta has float += half * half
         use it for all matrices
      */
      if (Method::TensorCores == method) {
        BsrMatrixSpMVTensorCoreDispatcher<ExecutionSpace, AMatrix, half, XVector, half, YVector, float, 16, 16,
                                          16>::dispatch(space, alpha, A, X, beta, Y);
        return;
      }
    }
#endif  // defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ARCH_AMPERE)

    const bool modeIsNoTrans        = (mode[0] == NoTranspose[0]);
    const bool modeIsConjugate      = (mode[0] == Conjugate[0]);
    const bool modeIsConjugateTrans = (mode[0] == ConjugateTranspose[0]);
    const bool modeIsTrans          = (mode[0] == Transpose[0]);

    // use V41 if requested
    if (handle->algo == SPMV_BSR_V41) {
      if (modeIsNoTrans || modeIsConjugate) {
        return Bsr::spMatMultiVec_no_transpose(space, handle, alpha, A, X, beta, Y, modeIsConjugate);
      } else if (modeIsTrans || modeIsConjugateTrans) {
        return Bsr::spMatMultiVec_transpose(space, handle, alpha, A, X, beta, Y, modeIsConjugateTrans);
      }
    }

    // use V42 if possible
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>() || handle->algo == SPMV_BSR_V42) {
      if (modeIsNoTrans) {
        ::KokkosSparse::Impl::apply_v42(space, alpha, A, X, beta, Y);
        return;
      }
    }

    // use V41 as the ultimate fallback
    if (modeIsNoTrans || modeIsConjugate) {
      return Bsr::spMatMultiVec_no_transpose(space, handle, alpha, A, X, beta, Y, modeIsConjugate);
    } else if (modeIsTrans || modeIsConjugateTrans) {
      return Bsr::spMatMultiVec_transpose(space, handle, alpha, A, X, beta, Y, modeIsConjugateTrans);
    }

    {
      std::stringstream ss;
      ss << __FILE__ << ":" << __LINE__ << " ";
      ss << "Internal logic error: no applicable BsrMatrix SpMV implementation "
            ". Please report this";
      throw std::runtime_error(ss.str());
    }
  }
};

template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct SPMV_MV_BSRMATRIX<ExecutionSpace, Handle, AMatrix, XVector, YVector, true, false,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef typename YVector::non_const_value_type YScalar;

  static void spmv_mv_bsrmatrix(const ExecutionSpace &space, Handle *handle, const char mode[], const YScalar &alpha,
                                const AMatrix &A, const XVector &X, const YScalar &beta, const YVector &Y) {
    static_assert(std::is_integral_v<typename AMatrix::non_const_value_type>,
                  "This implementation is only for integer Scalar types.");
    for (size_t j = 0; j < X.extent(1); ++j) {
      const auto x_j = Kokkos::subview(X, Kokkos::ALL(), j);
      auto y_j       = Kokkos::subview(Y, Kokkos::ALL(), j);
      typedef SPMV_BSRMATRIX<ExecutionSpace, Handle, AMatrix, decltype(x_j), decltype(y_j)> impl_type;
      impl_type::spmv_bsrmatrix(space, handle, mode, alpha, A, x_j, beta, y_j);
    }
  }
};
#endif  // !defined(KOKKOSKERNELS_ETI_ONLY) ||
        // KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
}  // namespace Impl
}  // namespace KokkosSparse

// declare / instantiate the vector version
// Instantiate with A,x,y are all the requested Scalar type (no instantiation of
// mixed-precision operands)
#define KOKKOSSPARSE_SPMV_BSRMATRIX_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,             \
                                                  EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                                 \
  extern template struct SPMV_BSRMATRIX<                                                                           \
      EXEC_SPACE_TYPE,                                                                                             \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>, \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                               \
                                              Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                     \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,         \
      Kokkos::View<SCALAR_TYPE const *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      false, true>;

#define KOKKOSSPARSE_SPMV_BSRMATRIX_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,             \
                                                  EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                                 \
  template struct SPMV_BSRMATRIX<                                                                                  \
      EXEC_SPACE_TYPE,                                                                                             \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>, \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                               \
                                              Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                     \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,         \
      Kokkos::View<SCALAR_TYPE const *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      false, true>;

// declare / instantiate the 2D MV version
// Instantiate with A,x,y are all the requested Scalar type (no instantiation of
// mixed-precision operands)
#define KOKKOSSPARSE_SPMV_MV_BSRMATRIX_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,          \
                                                     EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                              \
  extern template struct SPMV_MV_BSRMATRIX<                                                                        \
      EXEC_SPACE_TYPE,                                                                                             \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>, \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                               \
                                              Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                     \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,         \
      Kokkos::View<SCALAR_TYPE const **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      std::is_integral_v<SCALAR_TYPE>, false, true>;

#define KOKKOSSPARSE_SPMV_MV_BSRMATRIX_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,          \
                                                     EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                              \
  template struct SPMV_MV_BSRMATRIX<                                                                               \
      EXEC_SPACE_TYPE,                                                                                             \
      KokkosSparse::Impl::SPMVHandleImpl<EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SCALAR_TYPE, OFFSET_TYPE, ORDINAL_TYPE>, \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR_TYPE, const ORDINAL_TYPE,                               \
                                              Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                     \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,         \
      Kokkos::View<SCALAR_TYPE const **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      std::is_integral_v<SCALAR_TYPE>, false, true>;

#include <KokkosSparse_spmv_bsrmatrix_tpl_spec_decl.hpp>

#endif  // KOKKOSSPARSE_IMPL_SPMV_BSRMATRIX_SPEC_HPP_
