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

NOX::Epetra::MultiVector::MultiVector(Epetra_MultiVector& source, 
				      NOX::CopyType type,
				      bool createView)
  : noxEpetraVectors(source.NumVectors())
{
  if (createView) {
    epetraMultiVec = 
      new Epetra_MultiVector(View, source, 0, source.NumVectors());
  }
  else {
  
    switch (type) {
      
    case DeepCopy:		// default behavior
      
      epetraMultiVec = new Epetra_MultiVector(source); 
      break;
      
    case ShapeCopy:
      
      epetraMultiVec = 
	new Epetra_MultiVector(source.Map(), source.NumVectors()); 
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

    epetraMultiVec = new Epetra_MultiVector(source); 
    break;

  case ShapeCopy:

    epetraMultiVec = 
      new Epetra_MultiVector(source.Map(), source.NumVectors()); 
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

    epetraMultiVec = new Epetra_MultiVector(source.getEpetraMultiVector()); 
    break;

  case ShapeCopy:

    epetraMultiVec = 
      new Epetra_MultiVector(source.getEpetraMultiVector().Map(), 
			     source.numVectors()); 
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

  delete epetraMultiVec;
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
				   vector<int>& index) {
  return setBlock(dynamic_cast<const NOX::Epetra::MultiVector&>(source),
		  index);
}

NOX::Abstract::MultiVector& 
NOX::Epetra::MultiVector::setBlock(const NOX::Epetra::MultiVector& source, 
				   vector<int>& index) {
  double* vecPtr;
  double *sourceVecPtr;
  int ind;
  int vecLength = length();
  
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
  int vecLength = length();
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

  delete epetraMultiVec;
  epetraMultiVec = tmp;
					     
  return *this;
}

NOX::Abstract::Vector&
NOX::Epetra::MultiVector::operator [] (int i)
{
  if ( i < 0 || i > noxEpetraVectors.size() ) {
    cerr << "NOX::Epetra::MultiVector::operator[]:  Error!  Invalid index " 
	 << i << endl;
    throw "NOX::Epetra Error";
  }
  if (noxEpetraVectors[i] == NULL) {
    Epetra_Vector* epetra_vec = epetraMultiVec->operator() (i);
    noxEpetraVectors[i] = new NOX::Epetra::Vector(*epetra_vec, NOX::DeepCopy,
						  true);
  }
  return *(noxEpetraVectors[i]);
}

const NOX::Abstract::Vector&
NOX::Epetra::MultiVector::operator [] (int i) const
{
  if ( i < 0 || i > noxEpetraVectors.size() ) {
    cerr << "NOX::Epetra::MultiVector::operator[]:  Error!  Invalid index " 
	 << i << endl;
    throw "NOX::Epetra Error";
  }
  if (noxEpetraVectors[i] == NULL) {
    Epetra_Vector* epetra_vec = epetraMultiVec->operator() (i);
    noxEpetraVectors[i] = new NOX::Epetra::Vector(*epetra_vec, NOX::DeepCopy,
						  true);
  }
  return *(noxEpetraVectors[i]);
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

NOX::Abstract::MultiVector* 
NOX::Epetra::MultiVector::clone(CopyType type) const
{
  NOX::Abstract::MultiVector* newVec = 
    new NOX::Epetra::MultiVector(*epetraMultiVec, type);
  return newVec;
}

NOX::Abstract::MultiVector* 
NOX::Epetra::MultiVector::clone(int numvecs) const
{
  NOX::Epetra::MultiVector* newVec = new NOX::Epetra::MultiVector(numvecs);
  newVec->epetraMultiVec = new Epetra_MultiVector(epetraMultiVec->Map(),
						  numvecs);
  return newVec;
}

NOX::Abstract::MultiVector* 
NOX::Epetra::MultiVector::subCopy(vector<int>& index) const
{
  int numvecs = index.size();
  NOX::Epetra::MultiVector* newVec = new NOX::Epetra::MultiVector(numvecs);
  newVec->epetraMultiVec = new Epetra_MultiVector(Copy, *epetraMultiVec,
						  &index[0], numvecs);
  return newVec;
}

NOX::Abstract::MultiVector* 
NOX::Epetra::MultiVector::subView(vector<int>& index) const
{
  int numvecs = index.size();
  NOX::Epetra::MultiVector* newVec = new NOX::Epetra::MultiVector(numvecs);
  newVec->epetraMultiVec = new Epetra_MultiVector(View, *epetraMultiVec,
						  &index[0], numvecs);
  return newVec;
}

void 
NOX::Epetra::MultiVector::norm(vector<double>& result,
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

void NOX::Epetra::MultiVector::print() const
{
  epetraMultiVec->Print(cout);
  return;
}

void NOX::Epetra::MultiVector::checkIndex(int idx) const 
{
  if ( idx < 0 || idx >= epetraMultiVec->NumVectors() ) {
    cerr << "NOX::Epetra::MultiVector:  Error!  Invalid index " << idx << endl;
    throw "NOX::Epetra Error";
  }
}
