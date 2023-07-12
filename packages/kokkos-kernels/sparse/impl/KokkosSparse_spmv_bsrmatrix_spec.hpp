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
#include "KokkosKernels_Controls.hpp"
#include "KokkosKernels_Error.hpp"
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosSparse_spmv_bsrmatrix_impl.hpp>
#endif

namespace KokkosSparse {
namespace Experimental {
namespace Impl {

// default is no eti available
template <class AT, class AO, class AD, class AM, class AS, class XT, class XL,
          class XD, class XM, class YT, class YL, class YD, class YM>
struct spmv_bsrmatrix_eti_spec_avail {
  enum : bool { value = false };
};

template <class AT, class AO, class AD, class AM, class AS, class XT, class XL,
          class XD, class XM, class YT, class YL, class YD, class YM,
          const bool integerScalarType =
              std::is_integral<typename std::decay<AT>::type>::value>
struct spmv_mv_bsrmatrix_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPMV_BSRMATRIX_ETI_SPEC_AVAIL(                       \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    MEM_SPACE_TYPE)                                                       \
  template <>                                                             \
  struct spmv_bsrmatrix_eti_spec_avail<                                   \
      const SCALAR_TYPE, const ORDINAL_TYPE,                              \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE,         \
      SCALAR_TYPE const *, LAYOUT_TYPE,                                   \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,     \
      SCALAR_TYPE *, LAYOUT_TYPE,                                         \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > {                         \
    enum : bool { value = true };                                         \
  };

#define KOKKOSSPARSE_SPMV_MV_BSRMATRIX_ETI_SPEC_AVAIL(                    \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    MEM_SPACE_TYPE)                                                       \
  template <>                                                             \
  struct spmv_mv_bsrmatrix_eti_spec_avail<                                \
      const SCALAR_TYPE, const ORDINAL_TYPE,                              \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE,         \
      SCALAR_TYPE const **, LAYOUT_TYPE,                                  \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,     \
      SCALAR_TYPE **, LAYOUT_TYPE,                                        \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > {                         \
    enum : bool { value = true };                                         \
  };

// Include which ETIs are available
#include <KokkosSparse_spmv_bsrmatrix_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spmv_bsrmatrix_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spmv_mv_bsrmatrix_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Experimental {
namespace Impl {

// declaration
template <class AT, class AO, class AD, class AM, class AS, class XT, class XL,
          class XD, class XM, class YT, class YL, class YD, class YM,
          bool tpl_spec_avail = spmv_bsrmatrix_tpl_spec_avail<
              AT, AO, AD, AM, AS, XT, XL, XD, XM, YT, YL, YD, YM>::value,
          bool eti_spec_avail = spmv_bsrmatrix_eti_spec_avail<
              AT, AO, AD, AM, AS, XT, XL, XD, XM, YT, YL, YD, YM>::value>
struct SPMV_BSRMATRIX {
  typedef BsrMatrix<AT, AO, AD, AM, AS> AMatrix;
  typedef Kokkos::View<XT, XL, XD, XM> XVector;
  typedef Kokkos::View<YT, YL, YD, YM> YVector;
  typedef typename YVector::non_const_value_type YScalar;

  static void spmv_bsrmatrix(
      const KokkosKernels::Experimental::Controls &controls, const char mode[],
      const YScalar &alpha, const AMatrix &A, const XVector &x,
      const YScalar &beta, const YVector &y);
};

// declaration
template <class AT, class AO, class AD, class AM, class AS, class XT, class XL,
          class XD, class XM, class YT, class YL, class YD, class YM,
          const bool integerScalarType =
              std::is_integral<typename std::decay<AT>::type>::value,
          bool tpl_spec_avail = spmv_mv_bsrmatrix_tpl_spec_avail<
              AT, AO, AD, AM, AS, XT, XL, XD, XM, YT, YL, YD, YM>::value,
          bool eti_spec_avail = spmv_mv_bsrmatrix_eti_spec_avail<
              AT, AO, AD, AM, AS, XT, XL, XD, XM, YT, YL, YD, YM>::value>
struct SPMV_MV_BSRMATRIX {
  typedef BsrMatrix<AT, AO, AD, AM, AS> AMatrix;
  typedef Kokkos::View<XT, XL, XD, XM> XVector;
  typedef Kokkos::View<YT, YL, YD, YM> YVector;
  typedef typename YVector::non_const_value_type YScalar;

  static void spmv_mv_bsrmatrix(
      const KokkosKernels::Experimental::Controls &controls, const char mode[],
      const YScalar &alpha, const AMatrix &A, const XVector &x,
      const YScalar &beta, const YVector &y);
};

// actual implementations to be compiled
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

template <class AT, class AO, class AD, class AM, class AS, class XT, class XL,
          class XD, class XM, class YT, class YL, class YD, class YM>
struct SPMV_BSRMATRIX<AT, AO, AD, AM, AS, XT, XL, XD, XM, YT, YL, YD, YM, false,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef BsrMatrix<AT, AO, AD, AM, AS> AMatrix;
  typedef Kokkos::View<XT, XL, XD, XM> XVector;
  typedef Kokkos::View<YT, YL, YD, YM> YVector;
  typedef typename YVector::non_const_value_type YScalar;

  static void spmv_bsrmatrix(
      const KokkosKernels::Experimental::Controls &controls, const char mode[],
      const YScalar &alpha, const AMatrix &A, const XVector &X,
      const YScalar &beta, const YVector &Y) {
    //
    if ((mode[0] == KokkosSparse::NoTranspose[0]) ||
        (mode[0] == KokkosSparse::Conjugate[0])) {
      bool useConjugate = (mode[0] == KokkosSparse::Conjugate[0]);
      return Bsr::spMatVec_no_transpose(controls, alpha, A, X, beta, Y,
                                        useConjugate);
    } else if ((mode[0] == KokkosSparse::Transpose[0]) ||
               (mode[0] == KokkosSparse::ConjugateTranspose[0])) {
      bool useConjugate = (mode[0] == KokkosSparse::ConjugateTranspose[0]);
      return Bsr::spMatVec_transpose(controls, alpha, A, X, beta, Y,
                                     useConjugate);
    }
  }
};

template <class AT, class AO, class AD, class AM, class AS, class XT, class XL,
          class XD, class XM, class YT, class YL, class YD, class YM>
struct SPMV_MV_BSRMATRIX<AT, AO, AD, AM, AS, XT, XL, XD, XM, YT, YL, YD, YM,
                         false, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef BsrMatrix<AT, AO, AD, AM, AS> AMatrix;
  typedef Kokkos::View<XT, XL, XD, XM> XVector;
  typedef Kokkos::View<YT, YL, YD, YM> YVector;
  typedef typename YVector::non_const_value_type YScalar;

  enum class Method {
    Fallback,    ///< Don't use tensor cores
    TensorCores  ///< use tensor cores
  };

  /// Precision to use in the tensor core implementation
  enum class Precision {
    Automatic,  ///< Use Double, unless operations match mixed precision
    Double,     ///< fp64 += fp64 * fp64
    Mixed       ///< fp32 += fp16 * fp16
  };

  static void spmv_mv_bsrmatrix(
      const KokkosKernels::Experimental::Controls &controls, const char mode[],
      const YScalar &alpha, const AMatrix &A, const XVector &X,
      const YScalar &beta, const YVector &Y) {
#if defined(KOKKOS_ARCH_AMPERE) || defined(KOKKOS_ARCH_VOLTA)
    Method method = Method::Fallback;
    {
      typedef typename AMatrix::non_const_value_type AScalar;
      typedef typename XVector::non_const_value_type XScalar;
      // try to use tensor cores if requested
      if (controls.getParameter("algorithm") == "experimental_bsr_tc")
        method = Method::TensorCores;
      // can't use tensor cores for complex
      if (Kokkos::Details::ArithTraits<YScalar>::is_complex)
        method = Method::Fallback;
      if (Kokkos::Details::ArithTraits<XScalar>::is_complex)
        method = Method::Fallback;
      if (Kokkos::Details::ArithTraits<AScalar>::is_complex)
        method = Method::Fallback;
      // can't use tensor cores outside GPU
      if (!KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename AMatrix::execution_space>())
        method = Method::Fallback;
      if (!KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename XVector::execution_space>())
        method = Method::Fallback;
      if (!KokkosKernels::Impl::kk_is_gpu_exec_space<
              typename YVector::execution_space>())
        method = Method::Fallback;
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
        Precision precision = Precision::Automatic;
        if (controls.getParameter("tc_precision") == "mixed")
          precision = Precision::Mixed;
        else if (controls.getParameter("tc_precision") == "double")
          precision = Precision::Double;

        switch (precision) {
          case Precision::Mixed: {
            BsrMatrixSpMVTensorCoreDispatcher<AMatrix, half, XVector, half,
                                              YVector, float, 16, 16,
                                              16>::dispatch(alpha, A, X, beta,
                                                            Y);
            return;
          }
          case Precision::Double: {
            BsrMatrixSpMVTensorCoreDispatcher<AMatrix, double, XVector, double,
                                              YVector, double, 8, 8,
                                              4>::dispatch(alpha, A, X, beta,
                                                           Y);
            return;
          }
          case Precision::Automatic:  // fallthrough
          default: {
            constexpr bool operandsHalfHalfFloat =
                std::is_same<AScalar, Half>::value &&
                std::is_same<XScalar, Half>::value &&
                std::is_same<YScalar, float>::value;
            if (operandsHalfHalfFloat) {
              BsrMatrixSpMVTensorCoreDispatcher<AMatrix, half, XVector, half,
                                                YVector, float, 16, 16,
                                                16>::dispatch(alpha, A, X, beta,
                                                              Y);
              return;
            } else {
              BsrMatrixSpMVTensorCoreDispatcher<AMatrix, double, XVector,
                                                double, YVector, double, 8, 8,
                                                4>::dispatch(alpha, A, X, beta,
                                                             Y);
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
        BsrMatrixSpMVTensorCoreDispatcher<AMatrix, half, XVector, half, YVector,
                                          float, 16, 16, 16>::dispatch(alpha, A,
                                                                       X, beta,
                                                                       Y);
        return;
      }
    }
#endif  // KOKKOS_ARCH

    if ((mode[0] == KokkosSparse::NoTranspose[0]) ||
        (mode[0] == KokkosSparse::Conjugate[0])) {
      bool useConjugate = (mode[0] == KokkosSparse::Conjugate[0]);
      return Bsr::spMatMultiVec_no_transpose(controls, alpha, A, X, beta, Y,
                                             useConjugate);
    } else if ((mode[0] == KokkosSparse::Transpose[0]) ||
               (mode[0] == KokkosSparse::ConjugateTranspose[0])) {
      bool useConjugate = (mode[0] == KokkosSparse::ConjugateTranspose[0]);
      return Bsr::spMatMultiVec_transpose(controls, alpha, A, X, beta, Y,
                                          useConjugate);
    }
  }
};

template <class AT, class AO, class AD, class AM, class AS, class XT, class XL,
          class XD, class XM, class YT, class YL, class YD, class YM>
struct SPMV_MV_BSRMATRIX<AT, AO, AD, AM, AS, XT, XL, XD, XM, YT, YL, YD, YM,
                         true, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  typedef BsrMatrix<AT, AO, AD, AM, AS> AMatrix;
  typedef Kokkos::View<XT, XL, XD, XM> XVector;
  typedef Kokkos::View<YT, YL, YD, YM> YVector;
  typedef typename YVector::non_const_value_type YScalar;

  static void spmv_mv_bsrmatrix(
      const KokkosKernels::Experimental::Controls &controls, const char mode[],
      const YScalar &alpha, const AMatrix &A, const XVector &X,
      const YScalar &beta, const YVector &Y) {
    static_assert(std::is_integral<AT>::value,
                  "This implementation is only for integer Scalar types.");
    typedef SPMV_BSRMATRIX<AT, AO, AD, AM, AS, typename XVector::value_type *,
                           XL, XD, XM, typename YVector::value_type *, YL, YD,
                           YM>
        impl_type;
    for (typename AMatrix::non_const_size_type j = 0; j < X.extent(1); ++j) {
      const auto x_j = Kokkos::subview(X, Kokkos::ALL(), j);
      auto y_j       = Kokkos::subview(Y, Kokkos::ALL(), j);
      impl_type::spmv_bsrmatrix(controls, mode, alpha, A, x_j, beta, y_j);
    }
  }
};
#endif  // !defined(KOKKOSKERNELS_ETI_ONLY) ||
        // KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosSparse

// declare / instantiate the vector version
// Instantiate with A,x,y are all the requested Scalar type (no instantiation of
// mixed-precision operands)
#define KOKKOSSPARSE_SPMV_BSRMATRIX_ETI_SPEC_DECL(                        \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    MEM_SPACE_TYPE)                                                       \
  extern template struct SPMV_BSRMATRIX<                                  \
      const SCALAR_TYPE, const ORDINAL_TYPE,                              \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE,         \
      SCALAR_TYPE const *, LAYOUT_TYPE,                                   \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,     \
      SCALAR_TYPE *, LAYOUT_TYPE,                                         \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, false, true>;

#define KOKKOSSPARSE_SPMV_BSRMATRIX_ETI_SPEC_INST(                        \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    MEM_SPACE_TYPE)                                                       \
  template struct SPMV_BSRMATRIX<                                         \
      const SCALAR_TYPE, const ORDINAL_TYPE,                              \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE,         \
      SCALAR_TYPE const *, LAYOUT_TYPE,                                   \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,     \
      SCALAR_TYPE *, LAYOUT_TYPE,                                         \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, false, true>;

// declare / instantiate the 2D MV version
// Instantiate with A,x,y are all the requested Scalar type (no instantiation of
// mixed-precision operands)
#define KOKKOSSPARSE_SPMV_MV_BSRMATRIX_ETI_SPEC_DECL(                         \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,     \
    MEM_SPACE_TYPE)                                                           \
  extern template struct SPMV_MV_BSRMATRIX<                                   \
      const SCALAR_TYPE, const ORDINAL_TYPE,                                  \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE,             \
      SCALAR_TYPE const **, LAYOUT_TYPE,                                      \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,         \
      SCALAR_TYPE **, LAYOUT_TYPE,                                            \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>,                                \
      std::is_integral<typename std::decay<SCALAR_TYPE>::type>::value, false, \
      true>;

#define KOKKOSSPARSE_SPMV_MV_BSRMATRIX_ETI_SPEC_INST(                         \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,     \
    MEM_SPACE_TYPE)                                                           \
  template struct SPMV_MV_BSRMATRIX<                                          \
      const SCALAR_TYPE, const ORDINAL_TYPE,                                  \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE,             \
      SCALAR_TYPE const **, LAYOUT_TYPE,                                      \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,         \
      SCALAR_TYPE **, LAYOUT_TYPE,                                            \
      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>,                                \
      std::is_integral<typename std::decay<SCALAR_TYPE>::type>::value, false, \
      true>;

#include <KokkosSparse_spmv_bsrmatrix_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosSparse_spmv_bsrmatrix_eti_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosSparse_spmv_mv_bsrmatrix_eti_spec_decl.hpp>

#endif  // KOKKOSSPARSE_IMPL_SPMV_BSRMATRIX_SPEC_HPP_
