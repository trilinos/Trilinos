/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
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

  static void spmv_mv_bsrmatrix(
      const KokkosKernels::Experimental::Controls &controls, const char mode[],
      const YScalar &alpha, const AMatrix &A, const XVector &X,
      const YScalar &beta, const YVector &Y) {
#if defined(KOKKOS_ARCH_AMPERE) || defined(KOKKOS_ARCH_VOLTA)
    // user explicitly requests a particular precision
    bool requestMixed  = false;
    bool requestDouble = false;
    if (controls.isParameter("tc_precision")) {
      if (controls.getParameter("tc_precision") == "mixed") {
        requestMixed = true;
      } else if (controls.getParameter("tc_precision") == "double") {
        requestDouble = true;
      }
    }
    //
    bool use_tc = false;
    if ((controls.isParameter("algorithm")) &&
        (controls.getParameter("algorithm") == "experimental_bsr_tc")) {
      if (Kokkos::Details::ArithTraits<YScalar>::is_complex == false)
        use_tc = true;
    }
#endif

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ARCH_AMPERE)
    typedef typename XVector::non_const_value_type XScalar;
    typedef typename AMatrix::non_const_value_type AScalar;
    typedef Kokkos::Experimental::half_t Half;

    /* Ampere has double += double * double and float += half * half

    use whichever is requested.
    If none requested, used mixed precision if the inputs are mixed, otherwise
    use double
    */

    // input precision matches a tensor core fragment type
    constexpr bool operandsHalfHalfFloat = std::is_same<AScalar, Half>::value &&
                                           std::is_same<XScalar, Half>::value &&
                                           std::is_same<YScalar, float>::value;

    if (use_tc) {
      if (requestMixed) {
        BsrMatrixSpMVTensorCoreDispatcher<AMatrix, half, XVector, half, YVector,
                                          float, 16, 16, 16>::dispatch(alpha, A,
                                                                       X, beta,
                                                                       Y);
        return;
      } else if (requestDouble) {
        BsrMatrixSpMVTensorCoreDispatcher<AMatrix, double, XVector, double,
                                          YVector, double, 8, 8,
                                          4>::dispatch(alpha, A, X, beta, Y);
        return;
      } else if (operandsHalfHalfFloat) {
        BsrMatrixSpMVTensorCoreDispatcher<AMatrix, half, XVector, half, YVector,
                                          float, 16, 16, 16>::dispatch(alpha, A,
                                                                       X, beta,
                                                                       Y);
        return;
      } else {
        BsrMatrixSpMVTensorCoreDispatcher<AMatrix, double, XVector, double,
                                          YVector, double, 8, 8,
                                          4>::dispatch(alpha, A, X, beta, Y);
        return;
      }
    }
#elif defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ARCH_VOLTA)
    /* Volta has float += half * half
       use it for all matrices
    */
    if (use_tc) {
      if (requestDouble) {
        KokkosKernels::Impl::throw_runtime_exception(
            "KokkosSparse::spmv[algorithm=experimental_bsr_tc] "
            "tc_precision=double unsupported KOKKOS_ARCH_VOLTA");
      }
      BsrMatrixSpMVTensorCoreDispatcher<AMatrix, half, XVector, half, YVector,
                                        float, 16, 16, 16>::dispatch(alpha, A,
                                                                     X, beta,
                                                                     Y);
      (void)requestMixed;  // unused
      return;
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
