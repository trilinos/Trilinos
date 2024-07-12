// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Epetra_VectorSpace_ScaledL2.H"
#include "Epetra_Vector.h"

NOX::Epetra::VectorSpaceScaledL2::
VectorSpaceScaledL2(const Teuchos::RCP<NOX::Epetra::Scaling>& s,
            NOX::Epetra::Scaling::ScaleType st) :
  scalingPtr(s),
  scaleType(st)
{

}

NOX::Epetra::VectorSpaceScaledL2::~VectorSpaceScaledL2()
{

}

double NOX::Epetra::VectorSpaceScaledL2::
innerProduct(const Epetra_Vector& a, const Epetra_Vector& b) const
{
  if ( Teuchos::is_null(tmpVectorPtr) )
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(a));

  *tmpVectorPtr = a;

  if (scaleType == NOX::Epetra::Scaling::Left) {
    // Do twice on a instead of once on a and once on b.
    scalingPtr->applyLeftScaling(*tmpVectorPtr, *tmpVectorPtr);
    scalingPtr->applyLeftScaling(*tmpVectorPtr, *tmpVectorPtr);
  }
  else {
    // Do twice on a instead of once on a and once on b.
    scalingPtr->applyRightScaling(*tmpVectorPtr, *tmpVectorPtr);
    scalingPtr->applyRightScaling(*tmpVectorPtr, *tmpVectorPtr);
  }

  double dot;
  tmpVectorPtr->Dot(b, &dot);
  return dot;
}

double NOX::Epetra::VectorSpaceScaledL2::
norm(const Epetra_Vector& a, NOX::Abstract::Vector::NormType type) const
{
  if ( Teuchos::is_null(tmpVectorPtr) )
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(a));

  *tmpVectorPtr = a;

  if (scaleType == NOX::Epetra::Scaling::Left) {
    scalingPtr->applyLeftScaling(*tmpVectorPtr, *tmpVectorPtr);
  }
  else {
    scalingPtr->applyRightScaling(*tmpVectorPtr, *tmpVectorPtr);
  }

  double value;
  switch (type) {
  case NOX::Abstract::Vector::MaxNorm:
    tmpVectorPtr->NormInf(&value);
    break;
  case NOX::Abstract::Vector::OneNorm:
    tmpVectorPtr->Norm1(&value);
    break;
  case NOX::Abstract::Vector::TwoNorm:
  default:
   tmpVectorPtr->Norm2(&value);
   break;
  }
  return value;
}
