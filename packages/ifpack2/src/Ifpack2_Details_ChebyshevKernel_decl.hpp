// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_Details_ChebyshevKernel_decl.hpp
/// \brief Declaration of Ifpack2::Details::ChebyshevKernel
///
/// \note This file and its contents are implementation details of
///   Ifpack2.

#ifndef IFPACK2_DETAILS_CHEBYSHEVKERNEL_DECL_HPP
#define IFPACK2_DETAILS_CHEBYSHEVKERNEL_DECL_HPP

#include "Ifpack2_config.h"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_MultiVector_fwd.hpp"
#include "Tpetra_Operator_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Teuchos_RCP.hpp"
#include <memory>

namespace Ifpack2 {
namespace Details {

/// \brief Compute scaled damped residual for Chebyshev.
///
/// This is an implementation detail of Ifpack2::Chebyshev.  Given a
/// linear system A*X=B and an "inverse diagonal" matrix (stored as a
/// vector) D_inv, it computes a "scaled damped" residual vector W :=
/// alpha*D_inv*(B-A*X) + beta*W.
///
/// \tparam TpetraOperatorType Specialization of Tpetra::Operator.
///
/// \note To Ifpack2 developers: We can't fuse this with X := X + W,
///   because data dependencies in the input X are not elementwise
///   (unless A is diagonal).
template<class TpetraOperatorType>
class ChebyshevKernel {
private:
  using SC = typename TpetraOperatorType::scalar_type;
  using LO = typename TpetraOperatorType::local_ordinal_type;
  using GO = typename TpetraOperatorType::global_ordinal_type;
  using NT = typename TpetraOperatorType::node_type;

  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using multivector_type = Tpetra::MultiVector<SC, LO, GO, NT>;
  using operator_type = Tpetra::Operator<SC, LO, GO, NT>;
  using vector_type = Tpetra::Vector<SC, LO, GO, NT>;

public:
  ChebyshevKernel (const Teuchos::RCP<const operator_type>& A,
                   const bool useNativeSpMV=false);

  void
  setMatrix (const Teuchos::RCP<const operator_type>& A);

  void
  compute (multivector_type& W,
           const SC& alpha,
           vector_type& D_inv,
           multivector_type& B,
           multivector_type& X,
           const SC& beta);

private:
  using import_type = Tpetra::Import<LO, GO, NT>;
  using export_type = Tpetra::Export<LO, GO, NT>;

  Teuchos::RCP<const operator_type> A_op_;
  Teuchos::RCP<const crs_matrix_type> A_crs_;
  Teuchos::RCP<const import_type> imp_;
  Teuchos::RCP<const export_type> exp_;
  std::unique_ptr<vector_type> X_colMap_;
  std::unique_ptr<multivector_type> V1_;

  Teuchos::RCP<vector_type> W_vec_, B_vec_, X_vec_;

  // External override to not fuse operations into a single kernel
  // And use native blas/SpMV operations
  bool useNativeSpMV_;

  // Do the Import, if needed, and return the column Map version of X.
  vector_type&
  importVector (vector_type& X_domMap);

  bool canFuse (const multivector_type& B) const;

  void
  unfusedCase (multivector_type& W,
               const SC& alpha,
               vector_type& D_inv,
               multivector_type& B,
               const operator_type& A,
               multivector_type& X,
               const SC& beta);

  void
  fusedCase (vector_type& W,
             const SC& alpha,
             vector_type& D_inv,
             vector_type& B,
             const crs_matrix_type& A,
             vector_type& X,
             const SC& beta);
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_CHEBYSHEVKERNEL_DECL_HPP
