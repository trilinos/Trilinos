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

#include "NOX_Epetra_MultiVector.H"
#include "NOX_Epetra_Vector.H"

#include "Epetra_LocalMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

NOX::Epetra::MultiVector::MultiVector(int numvecs) 
  : noxEpetraVectors(numvecs)
{
  for (unsigned int i=0; i<noxEpetraVectors.size(); i++)
    noxEpetraVectors[i] = NULL;
}

NOX::Epetra::MultiVector::
MultiVector(const Teuchos::RCP<Epetra_MultiVector>& source, 
	    NOX::CopyType type,
	    NOX::Epetra::MultiVector::MemoryType memoryType)
  : noxEpetraVectors(source->NumVectors())
{
  if (memoryType == NOX::Epetra::MultiVector::CreateView) {
    epetraMultiVec = source;
  }
  else {
  
    switch (type) {
      
    case DeepCopy:		// default behavior
      
      epetraMultiVec = Teuchos::rcp(new Epetra_MultiVector(*source)); 
      break;
      
    case ShapeCopy:
      
      epetraMultiVec = 
	Teuchos::rcp(new Epetra_MultiVector(source->Map(), 
					    source->NumVectors())); 
      break;  
    }
  }

  for (unsigned int i=0; i<noxEpetraVectors.size(); i++)
    noxEpetraVectors[i] = NULL;
}

NOX::Epetra::MultiVector::MultiVector(const Epetra_MultiVector& source, 
				      NOX::CopyType type) 
  : noxEpetraVectors(source.NumVectors())
{
  switch (type) {

  case DeepCopy:		// default behavior

    epetraMultiVec = Teuchos::rcp(new Epetra_MultiVector(source)); 
    break;

  case ShapeCopy:

    epetraMultiVec = 
      Teuchos::rcp(new Epetra_MultiVector(source.Map(), source.NumVectors())); 
    break;  

  }

  for (unsigned int i=0; i<noxEpetraVectors.size(); i++)
    noxEpetraVectors[i] = NULL;
}

NOX::Epetra::MultiVector::MultiVector(const NOX::Epetra::MultiVector& source, 
				      NOX::CopyType type) 
  : noxEpetraVectors(source.epetraMultiVec->NumVectors())
{

  switch (type) {

  case DeepCopy:		// default behavior

    epetraMultiVec = Teuchos::rcp(new Epetra_MultiVector(source.getEpetraMultiVector())); 
    break;

  case ShapeCopy:

    epetraMultiVec = 
      Teuchos::rcp(new Epetra_MultiVector(source.getEpetraMultiVector().Map(), 
					  source.numVectors())); 
    break;  

  }
  
  for (unsigned int i=0; i<noxEpetraVectors.size(); i++)
    noxEpetraVectors[i] = NULL;
}

NOX::Epetra::MultiVector::~MultiVector()
{
  for (unsigned int i=0; i<noxEpetraVectors.size(); i++)
    if (noxEpetraVectors[i] != NULL)
      delete noxEpetraVectors[i];
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::operator=(const Epetra_MultiVector& source)
{
  *epetraMultiVec = source;
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::operator=(const NOX::Abstract::MultiVector& source)
{
  return operator=(dynamic_cast<const NOX::Epetra::MultiVector&>(source));
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::operator=(const NOX::Epetra::MultiVector& source)
{
  *epetraMultiVec = *source.epetraMultiVec;
  return *this;
}

Epetra_MultiVector& 
NOX::Epetra::MultiVector::getEpetraMultiVector()
{
  return *epetraMultiVec;
}

const Epetra_MultiVector& 
NOX::Epetra::MultiVector::getEpetraMultiVector() const
{
  return *epetraMultiVec;
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::init(double value)
{
  epetraMultiVec->PutScalar(value);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::random(bool useSeed, int seed)
{
  if (useSeed)
    epetraMultiVec->SetSeed(seed);
  epetraMultiVec->Random();
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::setBlock(const NOX::Abstract::MultiVector& source, 
				   const std::vector<int>& index) {
  return setBlock(dynamic_cast<const NOX::Epetra::MultiVector&>(source),
		  index);
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::setBlock(const NOX::Epetra::MultiVector& source, 
				   const std::vector<int>& index) {
  double* vecPtr;
  double *sourceVecPtr;
  int ind;
  int vecLength = epetraMultiVec->MyLength();
  
  source.checkIndex(index.size()-1);
  for (unsigned int i=0; i<index.size(); i++) {
    ind = index[i];
    checkIndex(ind);
    vecPtr = epetraMultiVec->operator[] (ind);
    sourceVecPtr = source.epetraMultiVec->operator[] (i);
    for (int j = 0; j<vecLength; j++)
      vecPtr[j] = sourceVecPtr[j];
  }
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::augment(const NOX::Abstract::MultiVector& source) {
  return augment(dynamic_cast<const NOX::Epetra::MultiVector&>(source));
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::augment(const NOX::Epetra::MultiVector& source) {
  int numVecs = numVectors();
  int sourceNumVecs = source.numVectors();
  int totalNumVecs = numVecs + sourceNumVecs;
  int vecLength = epetraMultiVec->MyLength();
  double* vecPtr;
  double *tmpVecPtr;
  
  Epetra_MultiVector *tmp = new Epetra_MultiVector(epetraMultiVec->Map(), 
						   totalNumVecs);
  for (int i=0; i<numVecs; i++) {
    vecPtr = epetraMultiVec->operator[] (i);
    tmpVecPtr = tmp->operator[] (i);
    for (int j = 0; j<vecLength; j++)
      tmpVecPtr[j] = vecPtr[j];
  }

  for (int i=0; i<sourceNumVecs; i++) {
    vecPtr = source.epetraMultiVec->operator[] (i);
    tmpVecPtr = tmp->operator[] (numVecs + i);
    for (int j = 0; j<vecLength; j++)
      tmpVecPtr[j] = vecPtr[j];
  }

  epetraMultiVec = Teuchos::rcp(tmp);
					     
  return *this;
}

NOX::Abstract::Vector&
NOX::Epetra::MultiVector::operator [] (int i)
{
  if ( i < 0 || i > (int) noxEpetraVectors.size() ) {
    std::cerr << "NOX::Epetra::MultiVector::operator[]:  Error!  Invalid index " 
	 << i << std::endl;
    throw "NOX::Epetra Error";
  }
  if (noxEpetraVectors[i] == NULL) {
    Teuchos::RCP<Epetra_Vector> epetra_vec = 
      Teuchos::rcp(epetraMultiVec->operator() (i), false);
    noxEpetraVectors[i] = 
      new NOX::Epetra::Vector(epetra_vec, NOX::Epetra::Vector::CreateView);
  }
  return *(noxEpetraVectors[i]);
}

const NOX::Abstract::Vector&
NOX::Epetra::MultiVector::operator [] (int i) const
{
  if ( i < 0 || i > (int) noxEpetraVectors.size() ) {
    std::cerr << "NOX::Epetra::MultiVector::operator[]:  Error!  Invalid index " 
	 << i << std::endl;
    throw "NOX::Epetra Error";
  }
  if (noxEpetraVectors[i] == NULL) {
    Teuchos::RCP<Epetra_Vector> epetra_vec = 
      Teuchos::rcp(epetraMultiVec->operator() (i), false);
    noxEpetraVectors[i] = 
      new NOX::Epetra::Vector(epetra_vec, NOX::Epetra::Vector::CreateView);
  }
  return *(noxEpetraVectors[i]);
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::scale(double gamma)
{
  epetraMultiVec->Scale(gamma);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::update(double alpha, 
				 const NOX::Abstract::MultiVector& a, 
				 double gamma)
{
  return update(alpha, dynamic_cast<const NOX::Epetra::MultiVector&>(a), 
		gamma);
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::update(double alpha, 
				 const NOX::Epetra::MultiVector& a, 
				 double gamma)
{
  epetraMultiVec->Update(alpha, a.getEpetraMultiVector(), gamma);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::update(double alpha, 
				 const NOX::Abstract::MultiVector& a, 
				 double beta, 
				 const NOX::Abstract::MultiVector& b,
				 double gamma)
{
  return update(alpha, dynamic_cast<const NOX::Epetra::MultiVector&>(a), 
		beta, dynamic_cast<const NOX::Epetra::MultiVector&>(b), gamma);
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::update(double alpha, 
				 const NOX::Epetra::MultiVector& a, 
				 double beta, 
				 const NOX::Epetra::MultiVector& b,
				 double gamma)
{
  epetraMultiVec->Update(alpha, a.getEpetraMultiVector(), 
			 beta, b.getEpetraMultiVector(), gamma);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::update(
			     Teuchos::ETransp transb, 
			     double alpha, 
			     const NOX::Abstract::MultiVector& a, 
			     const NOX::Abstract::MultiVector::DenseMatrix& b,
			     double gamma)
{
  return update(transb, alpha, 
		dynamic_cast<const NOX::Epetra::MultiVector&>(a), 
		b, gamma);
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::update(
			     Teuchos::ETransp transb, 
			     double alpha, 
			     const NOX::Epetra::MultiVector& a, 
			     const NOX::Abstract::MultiVector::DenseMatrix& b,
			     double gamma)
{
  // Remove const from b
  NOX::Abstract::MultiVector::DenseMatrix& nc_b = 
    const_cast<NOX::Abstract::MultiVector::DenseMatrix&>(b);

  // Create a replicated-local Epetra_MultiVector using b (view)
  const int izero = 0;
  Epetra_LocalMap localMap(b.numRows(), izero, epetraMultiVec->Map().Comm());
  Epetra_MultiVector B(View, localMap, nc_b.values(), b.stride(), b.numCols());

  char epetra_trans_b;
  if (transb == Teuchos::NO_TRANS)
    epetra_trans_b = 'N';
  else
    epetra_trans_b = 'T';
  
  epetraMultiVec->Multiply('N', epetra_trans_b, alpha, 
			   a.getEpetraMultiVector(), B, gamma);
  return *this;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Epetra::MultiVector::clone(CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::MultiVector> newVec = 
    Teuchos::rcp(new NOX::Epetra::MultiVector(*epetraMultiVec, type));
  return newVec;
}

Teuchos::RCP<NOX::Abstract::MultiVector> 
NOX::Epetra::MultiVector::clone(int numvecs) const
{
  Teuchos::RCP<NOX::Epetra::MultiVector> newVec = 
    Teuchos::rcp(new NOX::Epetra::MultiVector(numvecs));
  newVec->epetraMultiVec = 
    Teuchos::rcp(new Epetra_MultiVector(epetraMultiVec->Map(), numvecs));
  return newVec;
}

Teuchos::RCP<NOX::Abstract::MultiVector> 
NOX::Epetra::MultiVector::subCopy(const std::vector<int>& index) const
{
  int numvecs = index.size();
  Teuchos::RCP<NOX::Epetra::MultiVector> newVec = 
    Teuchos::rcp(new NOX::Epetra::MultiVector(numvecs));
  newVec->epetraMultiVec = 
    Teuchos::rcp(new Epetra_MultiVector(Copy, *epetraMultiVec,
					const_cast<int*>(&index[0]), 
					numvecs));
  return newVec;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Epetra::MultiVector::subView(const std::vector<int>& index) const
{
  int numvecs = index.size();
  Teuchos::RCP<NOX::Epetra::MultiVector> newVec = 
    Teuchos::rcp(new NOX::Epetra::MultiVector(numvecs));
  newVec->epetraMultiVec = 
    Teuchos::rcp(new Epetra_MultiVector(View, *epetraMultiVec,
					const_cast<int*>(&index[0]), 
					numvecs));
  return newVec;
}

void 
NOX::Epetra::MultiVector::norm(std::vector<double>& result,
			       NOX::Abstract::Vector::NormType type) const
{
  switch (type) {
  case NOX::Abstract::Vector::MaxNorm:
    epetraMultiVec->NormInf(&result[0]);
    break;
  case NOX::Abstract::Vector::OneNorm:
    epetraMultiVec->Norm1(&result[0]);
    break;
  case NOX::Abstract::Vector::TwoNorm:
  default:
    epetraMultiVec->Norm2(&result[0]);
    break;
  }
}

void 
NOX::Epetra::MultiVector::multiply(
			     double alpha, 
			     const NOX::Abstract::MultiVector& y,
			     NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  multiply(alpha, dynamic_cast<const NOX::Epetra::MultiVector&>(y), b);
}

void
NOX::Epetra::MultiVector::multiply(
			     double alpha, 
			     const NOX::Epetra::MultiVector& y,
			     NOX::Abstract::MultiVector::DenseMatrix& b) const
{
    // Create a replicated-local Epetra_MultiVector using b (view)
  const int izero = 0;
  Epetra_LocalMap localMap(b.numRows(), izero, epetraMultiVec->Map().Comm());
  Epetra_MultiVector B(View, localMap, b.values(), b.stride(), b.numCols());

  
  B.Multiply('T', 'N', alpha, *(y.epetraMultiVec), *epetraMultiVec, 0.0);
}

int NOX::Epetra::MultiVector::length() const
{
  return epetraMultiVec->GlobalLength();
}

int NOX::Epetra::MultiVector::numVectors() const
{
  return epetraMultiVec->NumVectors();
}

void NOX::Epetra::MultiVector::print(std::ostream& stream) const
{
  epetraMultiVec->Print(stream);
  return;
}

void NOX::Epetra::MultiVector::checkIndex(int idx) const 
{
  if ( idx < 0 || idx >= epetraMultiVec->NumVectors() ) {
    std::cerr << "NOX::Epetra::MultiVector:  Error!  Invalid index " << idx << std::endl;
    throw "NOX::Epetra Error";
  }
}
