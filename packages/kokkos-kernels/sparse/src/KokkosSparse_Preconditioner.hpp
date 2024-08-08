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

/// @file KokkosSparse_Preconditioner.hpp

#ifndef KK_PREC_HPP
#define KK_PREC_HPP

#include <Kokkos_Core.hpp>
#include <KokkosKernels_Controls.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace KokkosSparse {
namespace Experimental {

/// \class Preconditioner
/// \brief Interface for KokkosKernels preconditioners
/// \tparam CRS Type of the compressed matrix
///
/// Preconditioner provides the following methods
///   - initialize() performs all operations based on the graph of the
///     matrix (without considering the numerical values)
///   - isInitialized() returns true if the preconditioner has been
///     successfully initialized
///   - compute() computes everything required to apply the
///     preconditioner, using the matrix's values (and assuming that the
///     graph structure of the matrix has not changed)
///   - isComputed() returns true if the preconditioner has been
///     successfully computed, false otherwise.
///
/// Implementations of compute() must internally call initialize() if
/// isInitialized() returns false. The preconditioner is applied by
/// apply().
/// Every time that initialize() is called, the object destroys all the
/// previously allocated information, and reinitializes the preconditioner.
/// Every time compute() is called, the object recomputes the actual values of
/// the preconditioner.
template <class CRS>
class Preconditioner {
 public:
  using ScalarType = typename std::remove_const<typename CRS::value_type>::type;
  using EXSP       = typename CRS::execution_space;
  using MEMSP      = typename CRS::memory_space;
  using karith     = typename Kokkos::ArithTraits<ScalarType>;

  //! Constructor:
  Preconditioner() {}

  //! Destructor.
  virtual ~Preconditioner() {}

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
                     ScalarType alpha = karith::one(), ScalarType beta = karith::zero()) const = 0;
  //@}

  //! Set this preconditioner's parameters.
  virtual void setParameters() = 0;

  /// @brief Set up the graph structure of this preconditioner.
  ///
  /// If the graph structure of the constructor's input matrix has
  /// changed, or if you have not yet called initialize(), you must
  /// call initialize() before you may call compute() or apply().
  ///
  /// Thus, initialize() corresponds to the "symbolic factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  virtual void initialize() = 0;

  //! True if the preconditioner has been successfully initialized, else false.
  virtual bool isInitialized() const = 0;

  /// @brief Set up the numerical values in this preconditioner.
  ///
  /// If the values of the constructor's input matrix have changed, or
  /// if you have not yet called compute(), you must call compute()
  /// before you may call apply().
  ///
  /// Thus, compute() corresponds to the "numeric factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  virtual void compute() = 0;

  //! True if the preconditioner has been successfully computed, else false.
  virtual bool isComputed() const = 0;

  //! True if the preconditioner implements a transpose operator apply.
  virtual bool hasTransposeApply() const { return false; }
};

}  // namespace Experimental
}  // End namespace KokkosSparse

#endif
