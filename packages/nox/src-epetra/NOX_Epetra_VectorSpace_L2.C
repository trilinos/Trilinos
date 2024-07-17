// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Epetra_VectorSpace_L2.H"
#include "Epetra_Vector.h"

NOX::Epetra::VectorSpaceL2::VectorSpaceL2()
{

}

NOX::Epetra::VectorSpaceL2::~VectorSpaceL2()
{

}

double NOX::Epetra::VectorSpaceL2::
innerProduct(const Epetra_Vector& a, const Epetra_Vector& b) const
{
  double dot;
  a.Dot(b, &dot);
  return dot;
}

double NOX::Epetra::VectorSpaceL2::
norm(const Epetra_Vector& a, NOX::Abstract::Vector::NormType type) const
{
  double n;
  switch (type) {
  case NOX::Abstract::Vector::MaxNorm:
    a.NormInf(&n);
    break;
  case NOX::Abstract::Vector::OneNorm:
    a.Norm1(&n);
    break;
  case NOX::Abstract::Vector::TwoNorm:
  default:
   a.Norm2(&n);
   break;
  }
  return n;
}
