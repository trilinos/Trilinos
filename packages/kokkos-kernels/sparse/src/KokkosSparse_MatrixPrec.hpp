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

/// @file KokkosSparse_MatrixPrec.hpp

#ifndef KK_MATRIX_PREC_HPP
#define KK_MATRIX_PREC_HPP

#include <KokkosSparse_Preconditioner.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosBlas.hpp>
#include <KokkosSparse_spmv.hpp>

namespace KokkosSparse {

namespace Experimental {

/// @file KokkosSparse_MatrixPrec.hpp
/// \class MatrixPrec
/// \brief  This is a simple class to use if one
///         already has a matrix representation of their
///         preconditioner M.  The class applies an
///         SpMV with M as the preconditioning step.
/// \tparam CRS the type of compressed matrix
///
/// MatrixPrec provides the following methods
///   - initialize() Does nothing; Matrix initialized upon object construction.
///   - isInitialized() returns true
///   - compute() Does nothing; Matrix initialized upon object construction.
///   - isComputed() returns true
///
template <class CRS>
class MatrixPrec : public KokkosSparse::Experimental::Preconditioner<CRS> {
 private:
  CRS _A;

 public:
  using ScalarType = typename std::remove_const<typename CRS::value_type>::type;
  using EXSP       = typename CRS::execution_space;
  using MEMSP      = typename CRS::memory_space;
  using karith     = typename Kokkos::ArithTraits<ScalarType>;

  //! Constructor:
  template <class CRSArg>
  MatrixPrec(const CRSArg &mat) : _A(mat) {}

  //! Destructor.
  virtual ~MatrixPrec() {}

  ///// \brief Apply the preconditioner to X, putting the result in Y.
  /////
  ///// \tparam XViewType Input vector, as a 1-D Kokkos::View
  ///// \tparam YViewType Output vector, as a nonconst 1-D Kokkos::View
  /////
  ///// \param transM [in] "N" for non-transpose, "T" for transpose, "C"
  /////   for conjugate transpose.  All characters after the first are
  /////   ignored.  This works just like the BLAS routines.
  ///// \param alpha [in] Input coefficient of M*x
  ///// \param beta [in] Input coefficient of Y
  /////
  ///// If the result of applying this preconditioner to a vector X is
  ///// \f$M \cdot X\f$, then this method computes \f$Y = \beta Y + \alpha M
  ///\cdot X\f$.
  ///// The typical case is \f$\beta = 0\f$ and \f$\alpha = 1\f$.
  //
  virtual void apply(const Kokkos::View<const ScalarType *, Kokkos::Device<EXSP, MEMSP>> &X,
                     const Kokkos::View<ScalarType *, Kokkos::Device<EXSP, MEMSP>> &Y, const char transM[] = "N",
                     ScalarType alpha = karith::one(), ScalarType beta = karith::zero()) const {
    KokkosSparse::spmv(transM, alpha, _A, X, beta, Y);
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
