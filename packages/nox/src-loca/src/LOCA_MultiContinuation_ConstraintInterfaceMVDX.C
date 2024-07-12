// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterfaceMVDX::multiplyDX(
               double alpha,
               const NOX::Abstract::MultiVector& input_x,
               NOX::Abstract::MultiVector::DenseMatrix& result_p) const
{


  if (!isDXZero()) {
    const NOX::Abstract::MultiVector* dgdx = getDX();
    input_x.multiply(alpha, *dgdx, result_p);
  }
  else
    result_p.putScalar(0.0);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterfaceMVDX::addDX(
                      Teuchos::ETransp transb,
                  double alpha,
                      const NOX::Abstract::MultiVector::DenseMatrix& b,
                  double beta,
                  NOX::Abstract::MultiVector& result_x) const
{
  if (!isDXZero()) {
    const NOX::Abstract::MultiVector* dgdx = getDX();
    result_x.update(transb, alpha, *dgdx, b, beta);
  }
  else
    result_x.scale(beta);

  return NOX::Abstract::Group::Ok;
}
