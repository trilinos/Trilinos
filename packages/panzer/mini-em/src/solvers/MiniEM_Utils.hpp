// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_Utils_hpp_
#define _MiniEM_Utils_hpp_

#include "PanzerAdaptersSTK_config.hpp"
#include "Teuchos_RCP.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Panzer_NodeType.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Thyra_EpetraLinearOp.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "Epetra_CombineMode.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Thyra_EpetraThyraWrappers.hpp"
#endif
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Teko_Utilities.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"

namespace mini_em {

  /// \brief Writes the given Tpetra-backed operator to Matrix Market files prefixed with s (row/col/domain/range maps and the sparse matrix, or just the diagonal if the operator is not a CRS matrix).
  void writeOut(const std::string & s,const Thyra::LinearOpBase<double> & op);

  /// \brief Describes (prints structural info about) the given Tpetra-backed operator to out, if out is non-null.
  void describeMatrix(const std::string & s,const Thyra::LinearOpBase<double> & op,Teuchos::RCP<Teuchos::FancyOStream> out);

  /// \brief Calls describeMatrix(), then also calls writeOut() (writing to "s.mm") if doWrite is true.
  void describeAndWriteMatrix(const std::string & s, const Thyra::LinearOpBase<double> & op, Teuchos::RCP<Teuchos::FancyOStream> out, const bool doWrite);

  /// \brief Returns the underlying Tpetra::CrsMatrix wrapped by a Tpetra-backed Thyra operator. Asserts if op is not Tpetra-backed or not a CRS matrix.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node=Tpetra::Map<>::node_type>
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > get_Tpetra_CrsMatrix(const Thyra::LinearOpBase<double> & op);

#ifdef PANZER_HAVE_EPETRA_STACK
  Teuchos::RCP<const Epetra_CrsMatrix> get_Epetra_CrsMatrix(const Thyra::LinearOpBase<double> & op);

  Teuchos::RCP<const Epetra_CrsMatrix> get_Epetra_CrsMatrix(const Thyra::DiagonalLinearOpBase<double> & op, const Epetra_Comm& comm);
#endif

  /// \brief Builds a (scaled) identity operator with the same row map and range/domain spaces as the given Tpetra-backed operator.
  Teko::LinearOp getIdentityMatrix(const Teko::LinearOp& op, double scaling);

  /// \brief Returns true if op is a Tpetra-backed operator whose underlying Tpetra operator is not a Tpetra::CrsMatrix (i.e. it is matrix-free).
  bool isMatrixFreeOperator(const Teko::LinearOp& op);

  /// \brief Computes a lumped (row-sum) diagonal approximation of op's inverse, by applying op to the all-ones vector and taking the reciprocal.
  Teko::LinearOp getLumpedInverseDiagonal(const Teko::LinearOp& op);
}
#endif
