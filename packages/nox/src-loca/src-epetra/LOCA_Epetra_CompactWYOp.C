// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Epetra_config.h"
#include "LOCA_Epetra_CompactWYOp.H"

#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Epetra::CompactWYOp::CompactWYOp(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<const Epetra_Operator>& jacOperator,
    const Teuchos::RCP<const Epetra_MultiVector>& A_multiVec,
    const Teuchos::RCP<const Epetra_MultiVector>& Y_x_multiVec,
    const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& Y_p_matrix,
    const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& T_matrix) :
  globalData(global_data),
  label("LOCA::Epetra::CompactWYOp"),
  localMap(Y_x_multiVec->NumVectors(), 0, jacOperator->Comm()),
  J(jacOperator),
  A(A_multiVec),
  Y_x(Y_x_multiVec),
  Y_p(View, localMap, Y_p_matrix->values(), Y_p_matrix->stride(),
      Y_p_matrix->numCols()),
  T(View, localMap, T_matrix->values(), T_matrix->stride(),
    T_matrix->numCols()),
  tmpMat1(NULL),
  tmpMV(NULL),
  dblas()
{
}

LOCA::Epetra::CompactWYOp::~CompactWYOp()
{
  if (tmpMat1 != NULL)
    delete tmpMat1;
  if (tmpMV != NULL)
    delete tmpMV;
}

int
LOCA::Epetra::CompactWYOp::SetUseTranspose(bool UseTranspose)
{
  if (UseTranspose == false)
    return 0;
  else {
    globalData->locaErrorCheck->throwError(
                 "LOCA::Epetra::CompactWYOp::SetUseTranspose",
                 "Operator does not support transpose");
    return -1;
  }
}

int
LOCA::Epetra::CompactWYOp::Apply(const Epetra_MultiVector& Input,
                 Epetra_MultiVector& Result) const
{
  if (tmpMat1 == NULL || tmpMV == NULL) {
    globalData->locaErrorCheck->throwError(
                 "LOCA::Epetra::CompactWYOp::Apply()",
                 "Must call init() before Apply()!");
    return -1;
  }

  // Apply Householder transformation using temporary vector
  applyCompactWY(Input, *tmpMV, *tmpMat1);

  // Compute J*tmpMV
  J->Apply(*tmpMV, Result);

  // Compute J*tmpMV + A*tmpMat1
  if (A.get() != NULL)
    Result.Multiply('N', 'N', 1.0, *A, *tmpMat1, 1.0);

  return 0;
}

int
LOCA::Epetra::CompactWYOp::ApplyInverse(const Epetra_MultiVector& /* cInput */,
                    Epetra_MultiVector& /* Result */) const
{
  globalData->locaErrorCheck->throwError(
      "LOCA::Epetra::CompactWYOp::ApplyInverse",
      "Operator does not support ApplyInverse");
    return -1;
}

double
LOCA::Epetra::CompactWYOp::NormInf() const
{
  double Jn;
  std::vector<double> an(A->NumVectors());

  Jn = J->NormInf();
  A->NormInf(&an[0]);

  for (unsigned int i=0; i<an.size(); i++)
    Jn += an[i];

  return Jn;
}


const char*
LOCA::Epetra::CompactWYOp::Label () const
{
  return const_cast<char*>(label.c_str());
}

bool
LOCA::Epetra::CompactWYOp::UseTranspose() const
{
  return false;
}

bool
LOCA::Epetra::CompactWYOp::HasNormInf() const
{
  return J->HasNormInf();
}

const Epetra_Comm &
LOCA::Epetra::CompactWYOp::Comm() const
{
  return J->Comm();
}
const Epetra_Map&
LOCA::Epetra::CompactWYOp::OperatorDomainMap() const
{
  return J->OperatorDomainMap();
}

const Epetra_Map&
LOCA::Epetra::CompactWYOp::OperatorRangeMap() const
{
  return J->OperatorRangeMap();
}

void
LOCA::Epetra::CompactWYOp::init(const Epetra_MultiVector& x)
{
  if (tmpMat1 != NULL)
    delete tmpMat1;
  if (tmpMV != NULL)
    delete tmpMV;

  tmpMat1 = NULL;
  tmpMV = NULL;

  tmpMat1 = new Epetra_MultiVector(localMap, 1, false);
  tmpMV = new Epetra_MultiVector(x.Map(), 1, false);
}

void
LOCA::Epetra::CompactWYOp::finish()
{
  if (tmpMat1 != NULL)
    delete tmpMat1;
  if (tmpMV != NULL)
    delete tmpMV;

  tmpMat1 = NULL;
  tmpMV = NULL;
}

void
LOCA::Epetra::CompactWYOp::applyCompactWY(const Epetra_MultiVector& x,
                      Epetra_MultiVector& result_x,
                      Epetra_MultiVector& result_p) const
{
  // Compute Y_x^T*x
  result_p.Multiply('T', 'N', 1.0, *Y_x, x, 0.0);

  // Compute T*(Y_x^T*x)
  dblas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
         Teuchos::NON_UNIT_DIAG, result_p.MyLength(),
         result_p.NumVectors(), 1.0, T.Values(), T.MyLength(),
         result_p.Values(), result_p.MyLength());

  // Compute x = x + Y_x*T*(Y_x^T*x)
  result_x = x;
  result_x.Multiply('N', 'N', 1.0, *Y_x, result_p, 1.0);

  // Compute result_p = Y_p*T*(Y_x^T*x)
  dblas.TRMM(Teuchos::LEFT_SIDE, Teuchos::LOWER_TRI, Teuchos::NO_TRANS,
         Teuchos::UNIT_DIAG, result_p.MyLength(),
         result_p.NumVectors(), 1.0, Y_p.Values(), Y_p.MyLength(),
         result_p.Values(), result_p.MyLength());

}
