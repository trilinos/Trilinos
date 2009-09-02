// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
