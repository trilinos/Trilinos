// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

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
  using size_type  = typename CRS::size_type;
  using EXSP       = typename CRS::execution_space;
  using MEMSP      = typename CRS::memory_space;
  using DEVICE     = typename Kokkos::Device<EXSP, MEMSP>;
  using karith     = typename KokkosKernels::ArithTraits<ScalarType>;
  using View1d     = typename Kokkos::View<ScalarType *, DEVICE>;

 private:
  // trsm takes host views
  CRS L_, U_;
  View1d _tmp, _tmp2;
  mutable KernelHandle _khL;
  mutable KernelHandle _khU;

 public:
  //! Constructor:
  template <class CRSArg>
  LUPrec(const CRSArg &L, const CRSArg &U, const size_type block_size = 0)
      : L_(L), U_(U), _tmp("LUPrec::_tmp", L.numPointRows()), _tmp2("LUPrec::_tmp", L.numPointRows()), _khL(), _khU() {
    KK_REQUIRE_MSG(L.numPointRows() == U.numPointRows(), "LUPrec: L.numRows() != U.numRows()");

    _khL.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, L.numRows(), true, block_size);
    _khU.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, U.numRows(), false, block_size);
  }

  //! Destructor.
  virtual ~LUPrec() {
    _khL.destroy_sptrsv_handle();
    _khU.destroy_sptrsv_handle();
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
  virtual void apply(const Kokkos::View<const ScalarType *, DEVICE> &X, const Kokkos::View<ScalarType *, DEVICE> &Y,
                     const char transM[] = "N", ScalarType alpha = karith::one(),
                     ScalarType beta = karith::zero()) const override {
    KK_REQUIRE_MSG(transM[0] == NoTranspose[0], "LUPrec::apply only supports 'N' for transM");

    KokkosSparse::sptrsv_symbolic(&_khL, L_.graph.row_map, L_.graph.entries);
    KokkosSparse::sptrsv_solve(&_khL, L_.graph.row_map, L_.graph.entries, L_.values, X, _tmp);

    KokkosSparse::sptrsv_symbolic(&_khU, U_.graph.row_map, U_.graph.entries);
    KokkosSparse::sptrsv_solve(&_khU, U_.graph.row_map, U_.graph.entries, U_.values, _tmp, _tmp2);

    KokkosBlas::axpby(alpha, _tmp2, beta, Y);
  }
  //@}

  //! Set this preconditioner's parameters.
  void setParameters() override {}

  void initialize() override {}

  //! True if the preconditioner has been successfully initialized, else false.
  bool isInitialized() const override { return true; }

  void compute() override {}

  //! True if the preconditioner has been successfully computed, else false.
  bool isComputed() const override { return true; }

  //! True if the preconditioner implements a transpose operator apply.
  bool hasTransposeApply() const override { return true; }
};

}  // namespace Experimental
}  // End namespace KokkosSparse

#endif
