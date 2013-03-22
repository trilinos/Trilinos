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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include "NOX_Belos_MultiVector.H"

NOX::Belos::MultiVector::MultiVector(NOX::Abstract::MultiVector& noxMultiVec)
  : vecPtr(&noxMultiVec),
    ownsVec(false)
{
}

NOX::Belos::MultiVector::MultiVector(const NOX::Belos::MultiVector& source)
  : vecPtr(source.vecPtr->clone(NOX::DeepCopy).release()),
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
  NOX::Abstract::MultiVector *newVec = vecPtr->clone(numvecs).release();
  return new NOX::Belos::MultiVector(*newVec, true);
}

::Belos::MultiVec<double>*
NOX::Belos::MultiVector::CloneCopy()
{
  NOX::Abstract::MultiVector *newVec = vecPtr->clone(NOX::ShapeCopy).release();
  return new NOX::Belos::MultiVector(*newVec, true);
}

::Belos::MultiVec<double>*
NOX::Belos::MultiVector::CloneCopy(int index[], int numvecs)
{
  std::vector<int> idx(index, index+numvecs);
  NOX::Abstract::MultiVector *newVec = vecPtr->subCopy(idx).release();
  return new NOX::Belos::MultiVector(*newVec, true);
}
 
::Belos::MultiVec<double>*
NOX::Belos::MultiVector::CloneView(int index[], int numvecs)
{
  std::vector<int> idx(index, index+numvecs);
  NOX::Abstract::MultiVector *newVec = vecPtr->subView(idx).release();
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

  std::vector<double> res(vecPtr->numVectors());
  vecPtr->norm(res, nox_norm_type);

  for (unsigned int i=0; i<res.size(); i++)
    normvec[i] = res[i];

  return ::Belos::Ok;
}

void
NOX::Belos::MultiVector::SetBlock(::Belos::MultiVec<double>& A, int index[], 
				  int numvecs)
{
  std::vector<int> idx(index, index+numvecs);

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
  vecPtr->print(os);
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

  
