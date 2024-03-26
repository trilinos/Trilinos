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
// ************************************************************************
//@HEADER

/// @file KokkosSparse_LUPrec.hpp

#ifndef KK_LU_PREC_HPP
#define KK_LU_PREC_HPP

#include <KokkosSparse_Preconditioner.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosSparse_sptrsv.hpp>
#include <KokkosSparse_trsv.hpp>

namespace KokkosSparse {
namespace Experimental {

/// \class LUPrec
/// \brief  This class is for applying LU preconditioning.
///         It takes L and U and the apply method returns U^inv L^inv x
/// \tparam CRS the CRS type of L and U
///
/// Preconditioner provides the following methods
///   - initialize() Does nothing; members initialized upon object construction.
///   - isInitialized() returns true
///   - compute() Does nothing; members initialized upon object construction.
///   - isComputed() returns true
///
template <class CRS, class KernelHandle>
class LUPrec : public KokkosSparse::Experimental::Preconditioner<CRS> {
 public:
  using ScalarType = typename std::remove_const<typename CRS::value_type>::type;
  using EXSP       = typename CRS::execution_space;
  using MEMSP      = typename CRS::memory_space;
  using DEVICE     = typename Kokkos::Device<EXSP, MEMSP>;
  using karith     = typename Kokkos::ArithTraits<ScalarType>;
  using View1d     = typename Kokkos::View<ScalarType *, DEVICE>;

 private:
  // trsm takes host views
  CRS _L, _U;
  View1d _tmp, _tmp2;
  mutable KernelHandle _khL;
  mutable KernelHandle _khU;

 public:
  //! Constructor:
  template <class CRSArg>
  LUPrec(const CRSArg &L, const CRSArg &U)
      : _L(L),
        _U(U),
        _tmp("LUPrec::_tmp", L.numPointRows()),
        _tmp2("LUPrec::_tmp", L.numPointRows()),
        _khL(),
        _khU() {
    KK_REQUIRE_MSG(L.numPointRows() == U.numPointRows(),
                   "LUPrec: L.numRows() != U.numRows()");

    _khL.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, L.numRows(),
                              true);
    _khU.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, U.numRows(),
                              false);
  }

  //! Destructor.
  virtual ~LUPrec() {
    _khL.destroy_sptrsv_handle();
    _khU.destroy_sptrsv_handle();
  }

  template <
      typename Matrix,
      typename std::enable_if<is_crs_matrix<Matrix>::value>::type * = nullptr>
  void apply_impl(const Kokkos::View<const ScalarType *, DEVICE> &X,
                  const Kokkos::View<ScalarType *, DEVICE> &Y,
                  const char transM[] = "N", ScalarType alpha = karith::one(),
                  ScalarType beta = karith::zero()) const {
    // tmp = trsv(L, x); //Apply L^inv to x
    // y = trsv(U, tmp); //Apply U^inv to tmp

    KK_REQUIRE_MSG(transM[0] == NoTranspose[0],
                   "LUPrec::apply only supports 'N' for transM");

    sptrsv_symbolic(&_khL, _L.graph.row_map, _L.graph.entries);
    sptrsv_solve(&_khL, _L.graph.row_map, _L.graph.entries, _L.values, X, _tmp);

    sptrsv_symbolic(&_khU, _U.graph.row_map, _U.graph.entries);
    sptrsv_solve(&_khU, _U.graph.row_map, _U.graph.entries, _U.values, _tmp,
                 _tmp2);

    KokkosBlas::axpby(alpha, _tmp2, beta, Y);
  }

  template <
      typename Matrix,
      typename std::enable_if<is_bsr_matrix<Matrix>::value>::type * = nullptr>
  void apply_impl(const Kokkos::View<const ScalarType *, DEVICE> &X,
                  const Kokkos::View<ScalarType *, DEVICE> &Y,
                  const char transM[] = "N", ScalarType alpha = karith::one(),
                  ScalarType beta = karith::zero()) const {
    // tmp = trsv(L, x); //Apply L^inv to x
    // y = trsv(U, tmp); //Apply U^inv to tmp

    KK_REQUIRE_MSG(transM[0] == NoTranspose[0],
                   "LUPrec::apply only supports 'N' for transM");

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
    using Layout = Kokkos::LayoutLeft;
#else
    using Layout = Kokkos::LayoutRight;
#endif

    // trsv is implemented for MV so we need to convert our views
    using UView2d = typename Kokkos::View<
        ScalarType **, Layout, DEVICE,
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;
    using UView2dc = typename Kokkos::View<
        const ScalarType **, Layout, DEVICE,
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;
    UView2dc X2d(X.data(), X.extent(0), 1);
    UView2d Y2d(Y.data(), Y.extent(0), 1),
        tmp2d(_tmp.data(), _tmp.extent(0), 1),
        tmp22d(_tmp2.data(), _tmp2.extent(0), 1);

    KokkosSparse::trsv("L", "N", "N", _L, X2d, tmp2d);
    KokkosSparse::trsv("U", "N", "N", _U, tmp2d, tmp22d);

    KokkosBlas::axpby(alpha, _tmp2, beta, Y);
  }

  ///// \brief Apply the preconditioner to X, putting the result in Y.
  /////
  ///// \tparam XViewType Input vector, as a 1-D Kokkos::View
  ///// \tparam YViewType Output vector, as a nonconst 1-D Kokkos::View
  /////
  ///// \param transM [in] Not used.
  ///// \param alpha [in] Not used
  ///// \param beta [in] Not used.
  /////
  ///// It takes L and U and the stores U^inv L^inv X in Y
  //
  virtual void apply(const Kokkos::View<const ScalarType *, DEVICE> &X,
                     const Kokkos::View<ScalarType *, DEVICE> &Y,
                     const char transM[] = "N",
                     ScalarType alpha    = karith::one(),
                     ScalarType beta     = karith::zero()) const {
    apply_impl<CRS>(X, Y, transM, alpha, beta);
  }
  //@}

  //! Set this preconditioner's parameters.
  void setParameters() {}

  void initialize() {}

  //! True if the preconditioner has been successfully initialized, else false.
  bool isInitialized() const { return true; }

  void compute() {}

  //! True if the preconditioner has been successfully computed, else false.
  bool isComputed() const { return true; }

  //! True if the preconditioner implements a transpose operator apply.
  bool hasTransposeApply() const { return true; }
};

}  // namespace Experimental
}  // End namespace KokkosSparse

#endif
