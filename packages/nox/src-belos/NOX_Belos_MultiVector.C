// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_Belos_MultiVector.H"

NOX::Belos::MultiVector::MultiVector(NOX::Abstract::MultiVector& noxMultiVec)
  : vecPtr(&noxMultiVec),
    ownsVec(false)
{
}

NOX::Belos::MultiVector::MultiVector(const NOX::Belos::MultiVector& source)
  : vecPtr(source.vecPtr->clone(NOX::DeepCopy)),
    ownsVec(true)
{
}

NOX::Belos::MultiVector::MultiVector(NOX::Abstract::MultiVector& noxMultiVec,
			       bool ownsVecFlag)
  : vecPtr(&noxMultiVec),
    ownsVec(ownsVecFlag)
{
}

NOX::Belos::MultiVector::~MultiVector()
{
  if (ownsVec)
    delete vecPtr;
}

NOX::Belos::MultiVector&
NOX::Belos::MultiVector::operator=(const NOX::Belos::MultiVector& source)
{
  if (this != &source) {
    *vecPtr = *source.vecPtr;
  }

  return *this;
}

::Belos::MultiVec<double>*
NOX::Belos::MultiVector::Clone(const int numvecs)
{
  NOX::Abstract::MultiVector *newVec = vecPtr->clone(numvecs);
  return new NOX::Belos::MultiVector(*newVec, true);
}

::Belos::MultiVec<double>*
NOX::Belos::MultiVector::CloneCopy()
{
  NOX::Abstract::MultiVector *newVec = vecPtr->clone(NOX::ShapeCopy);
  return new NOX::Belos::MultiVector(*newVec, true);
}

::Belos::MultiVec<double>*
NOX::Belos::MultiVector::CloneCopy(int index[], int numvecs)
{
  vector<int> idx(index, index+numvecs);
  NOX::Abstract::MultiVector *newVec = vecPtr->subCopy(idx);
  return new NOX::Belos::MultiVector(*newVec, true);
}
 
::Belos::MultiVec<double>*
NOX::Belos::MultiVector::CloneView(int index[], int numvecs)
{
  vector<int> idx(index, index+numvecs);
  NOX::Abstract::MultiVector *newVec = vecPtr->subView(idx);
  return new NOX::Belos::MultiVector(*newVec, true);
}

int
NOX::Belos::MultiVector::GetVecLength() const
{
  return vecPtr->length();
}

int
NOX::Belos::MultiVector::GetNumberVecs() const
{
  return vecPtr->numVectors();
}

void
NOX::Belos::MultiVector::MvTimesMatAddMv(double alpha, 
				   ::Belos::MultiVec<double>& A, 
				   Teuchos::SerialDenseMatrix<int,double>& B,
				   double beta)
{
  // Cast A to a NOX::Belos::MultiVector
  NOX::Belos::MultiVector& nox_belos_A = 
    dynamic_cast<NOX::Belos::MultiVector&>(A);

  vecPtr->update(Teuchos::NO_TRANS, alpha, *(nox_belos_A.vecPtr), B, beta);
}

void
NOX::Belos::MultiVector:: MvAddMv(double alpha, 
			       ::Belos::MultiVec<double>& A, 
			       double beta, 
			       ::Belos::MultiVec<double>& B)
{
  // Cast A to a NOX::Belos::MultiVector
  NOX::Belos::MultiVector& nox_belos_A = 
    dynamic_cast<NOX::Belos::MultiVector&>(A);

  // Cast B to a NOX::Belos::MultiVector
  NOX::Belos::MultiVector& nox_belos_B = 
    dynamic_cast<NOX::Belos::MultiVector&>(B);

  vecPtr->update(alpha, *(nox_belos_A.vecPtr), beta, *(nox_belos_B.vecPtr),
		 0.0);
}

void
NOX::Belos::MultiVector::MvTransMv(double alpha, 
				::Belos::MultiVec<double>& A, 
				Teuchos::SerialDenseMatrix<int,double>& B)
{
  // Cast A to a NOX::Belos::MultiVector
  NOX::Belos::MultiVector& nox_belos_A = 
    dynamic_cast<NOX::Belos::MultiVector&>(A);

  vecPtr->multiply(alpha, *(nox_belos_A.vecPtr), B);
}

::Belos::ReturnType 
NOX::Belos::MultiVector::MvNorm(double *normvec, 
				::Belos::NormType norm_type)
{
  NOX::Abstract::Vector::NormType nox_norm_type;
  if (norm_type == ::Belos::OneNorm)
    nox_norm_type = NOX::Abstract::Vector::OneNorm;
  else if (norm_type == ::Belos::TwoNorm)
    nox_norm_type = NOX::Abstract::Vector::TwoNorm;
  else
    nox_norm_type = NOX::Abstract::Vector::MaxNorm;

  vector<double> res(vecPtr->numVectors());
  vecPtr->norm(res, nox_norm_type);

  for (unsigned int i=0; i<res.size(); i++)
    normvec[i] = res[i];

  return ::Belos::Ok;
}

void
NOX::Belos::MultiVector::SetBlock(::Belos::MultiVec<double>& A, int index[], 
				  int numvecs)
{
  vector<int> idx(index, index+numvecs);

  // Cast A to a NOX::Belos::MultiVector
  NOX::Belos::MultiVector& nox_belos_A = 
    dynamic_cast<NOX::Belos::MultiVector&>(A);

  vecPtr->setBlock(*(nox_belos_A.vecPtr), idx);
}

void
NOX::Belos::MultiVector::MvRandom()
{
  vecPtr->random();
}

void
NOX::Belos::MultiVector::MvInit(double alpha)
{
  vecPtr->init(alpha);
}

void
NOX::Belos::MultiVector::MvPrint(ostream& os)
{
  vecPtr->print();
}

NOX::Abstract::MultiVector&
NOX::Belos::MultiVector::getNoxMultiVector()
{
  return *vecPtr;
}

const NOX::Abstract::MultiVector&
NOX::Belos::MultiVector::getNoxMultiVector() const
{
  return *vecPtr;
}

::Belos::ReturnType
NOX::Belos::MultiVector::noxReturnTypeToBelos(
			    NOX::Abstract::Group::ReturnType noxStatus) const
{
  if (noxStatus == NOX::Abstract::Group::Ok ||
      noxStatus == NOX::Abstract::Group::NotConverged)
    return ::Belos::Ok;
  else if (noxStatus == NOX::Abstract::Group::NotDefined)
    return ::Belos::Undefined;
  else
    return ::Belos::Error;
}

  
