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

#include "NOX_Random.H" //for NOX::Random
#include "LOCA_Extended_Vector.H"  // Class definition
#include "LOCA_Extended_MultiVector.H" // for createMultiVector
#include "LOCA_ErrorCheck.H"
#include "LOCA_Utils.H"

LOCA::Extended::Vector::Vector(int nvecs, int nscalars) :
  vectorPtrs(nvecs), 
  isView(nvecs), 
  numScalars(nscalars), 
  scalarsPtr(NULL)
{
  scalarsPtr = new NOX::Abstract::MultiVector::DenseMatrix(numScalars, 1);
}

LOCA::Extended::Vector::Vector(const LOCA::Extended::Vector& source, 
			       NOX::CopyType type) :
  vectorPtrs(source.vectorPtrs.size()), 
  isView(source.vectorPtrs.size()),
  numScalars(source.numScalars),
  scalarsPtr(NULL)
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++) {
    vectorPtrs[i] = source.vectorPtrs[i]->clone(type);
    isView[i] = false;
  }

  scalarsPtr = new NOX::Abstract::MultiVector::DenseMatrix(*source.scalarsPtr);
  
  if (type != NOX::DeepCopy)
    init(0.0);
}


LOCA::Extended::Vector::~Vector()
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++) {
    if (!isView[i])
      delete vectorPtrs[i];
  }
  delete scalarsPtr;
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::operator=(const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Extended::Vector&>(y));
}

LOCA::Extended::Vector& 
LOCA::Extended::Vector::operator=(const LOCA::Extended::Vector& y)
{ 
  if (this != &y) {

    if (numScalars != y.numScalars) 
      LOCA::ErrorCheck::throwError("LOCA::Extended::Vector::operator=()",
			       "Number of scalars must match in assignment");
    if (vectorPtrs.size() != y.vectorPtrs.size())
      LOCA::ErrorCheck::throwError("LOCA::Extended::Vector::operator=()",
			       "Number of vectors must match in assignment");

    for (unsigned int i=0; i<vectorPtrs.size(); i++)
      *(vectorPtrs[i]) = *(y.vectorPtrs[i]);

    numScalars = y.numScalars;
    //*scalarsPtr = *y.scalarsPtr;
    // copy scalars.  We don't use operator= because it does the wrong thing
    // when either of the scalar arrays are views
    for (int i=0; i<numScalars; i++)
      (*scalarsPtr)(i,0) = (*y.scalarsPtr)(i,0);

  }
  return *this;
}

NOX::Abstract::Vector* 
LOCA::Extended::Vector::clone(NOX::CopyType type) const
{
  return new LOCA::Extended::Vector(*this, type);
}

NOX::Abstract::MultiVector* 
LOCA::Extended::Vector::createMultiVector(
				   const NOX::Abstract::Vector* const* vecs,
				   int numVecs, 
				   NOX::CopyType type) const
{
  // Pointers to sub-vectors of extended vectors
  const NOX::Abstract::Vector** subvecs = 
    new const NOX::Abstract::Vector*[numVecs+1];

  // Pointer to sub-multivector of extended multivector
  NOX::Abstract::MultiVector* submvec;

  // Create empty extended multivector
  LOCA::Extended::MultiVector *mvec = 
    generateMultiVector(numVecs+1, vectorPtrs.size(), numScalars);

  const LOCA::Extended::Vector* evec;

  // Create sub multivectors
  for (unsigned int i=0; i<vectorPtrs.size(); i++) {

    // Get the ith abstract vector from each column
    subvecs[0] = vectorPtrs[i];
    for (int j=0; j<numVecs; j++) {
      evec = dynamic_cast<const LOCA::Extended::Vector*>(vecs[j]);
      subvecs[j+1] = evec->vectorPtrs[i];
    }

    // Create multivector for the ith row
    submvec = vectorPtrs[i]->createMultiVector(subvecs, numVecs+1, type);

    // Set pointer to sub multivec
    mvec->setMultiVectorPtr(i, submvec);

  }

  // Get scalars
  for (int i=0; i<numScalars; i++) 
    mvec->getScalar(i,0) = (*scalarsPtr)(i,0);
  for (int j=0; j<numVecs; j++) {
    evec = dynamic_cast<const LOCA::Extended::Vector*>(vecs[j]);
    for (int i=0; i<numScalars; i++) 
      mvec->getScalar(i,j+1) = (*evec->scalarsPtr)(i,0);
  }

  delete [] subvecs;

  return mvec;
    
}

NOX::Abstract::MultiVector* 
LOCA::Extended::Vector::createMultiVector(int numVecs, 
					  NOX::CopyType type) const
{
  // Pointer to sub-multivector of extended multivector
  NOX::Abstract::MultiVector* submvec;

  // Create empty extended multivector
  LOCA::Extended::MultiVector *mvec = 
    generateMultiVector(numVecs, vectorPtrs.size(), numScalars);

  // Create sub multivectors
  for (unsigned int i=0; i<vectorPtrs.size(); i++) {

    // Create multivector for the ith row
    submvec = vectorPtrs[i]->createMultiVector(numVecs, type);

    // Set pointer to sub multivec
    mvec->setMultiVectorPtr(i, submvec);

  }

  // Get scalars
  if (type == NOX::DeepCopy) {
    for (int j=0; j<numVecs; j++) {
      for (int i=0; i<numScalars; i++)
	mvec->getScalar(i,j) = (*scalarsPtr)(i,0);
    }
  }

  return mvec;
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::init(double gamma)
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->init(gamma);
  scalarsPtr->putScalar(gamma);
  return *this;
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::random(bool useSeed, int seed) {
  if (useSeed)
    NOX::Random::setSeed(seed);
  if (vectorPtrs.size() > 0)
    vectorPtrs[0]->random(useSeed, seed);
  for (unsigned int i=1; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->random();
  scalarsPtr->random();
  return *this;
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::abs(const NOX::Abstract::Vector& y)
{
  const LOCA::Extended::Vector& Y = 
    dynamic_cast<const LOCA::Extended::Vector&>(y);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->abs(*(Y.vectorPtrs[i]));
  for (int i=0; i<numScalars; i++)
    (*scalarsPtr)(i,0) = fabs((*Y.scalarsPtr)(i,0));

  return *this;
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::reciprocal(const NOX::Abstract::Vector& y)
{
  const LOCA::Extended::Vector& Y = 
    dynamic_cast<const LOCA::Extended::Vector&>(y);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->reciprocal(*(Y.vectorPtrs[i]));
  for (int i=0; i<numScalars; i++)
    (*scalarsPtr)(i,0) = 1.0 / (*Y.scalarsPtr)(i,0);
  
  return *this;
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::scale(double gamma)
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->scale(gamma);
  scalarsPtr->scale(gamma);

  return *this;
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::scale(const NOX::Abstract::Vector& y)
{
  const LOCA::Extended::Vector& Y = 
    dynamic_cast<const LOCA::Extended::Vector&>(y);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->scale(*(Y.vectorPtrs[i]));
  for (int i=0; i<numScalars; i++)
    (*scalarsPtr)(i,0) *= (*Y.scalarsPtr)(i,0);

  return *this;
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::update(double alpha, const NOX::Abstract::Vector& y, 
			       double gamma)
{
  const LOCA::Extended::Vector& Y = 
    dynamic_cast<const LOCA::Extended::Vector&>(y);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->update(alpha, *(Y.vectorPtrs[i]), gamma);
  for (int i=0; i<numScalars; i++)
    (*scalarsPtr)(i,0) = alpha*((*Y.scalarsPtr)(i,0)) + gamma*((*scalarsPtr)(i,0));

  return *this;
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::update(double alpha, 
			       const NOX::Abstract::Vector& y, 
			       double beta, const NOX::Abstract::Vector& b,
			       double gamma)
{
  const LOCA::Extended::Vector& Y = 
    dynamic_cast<const LOCA::Extended::Vector&>(y);
  const LOCA::Extended::Vector& B = 
    dynamic_cast<const LOCA::Extended::Vector&>(b);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->update(alpha, *(Y.vectorPtrs[i]), beta, *(B.vectorPtrs[i]),
			  gamma);
  for (int i=0; i<numScalars; i++)
    (*scalarsPtr)(i,0) = alpha*((*Y.scalarsPtr)(i,0)) + beta*((*B.scalarsPtr)(i,0)) + 
      gamma*((*scalarsPtr)(i,0));

  return *this;
}

double 
LOCA::Extended::Vector::norm(NormType type) const
{
  double n = 0.0;
  double nv;

  switch (type) {

  case NOX::Abstract::Vector::MaxNorm:
    for (unsigned int i=0; i<vectorPtrs.size(); i++) 
      n = max(n, vectorPtrs[i]->norm(type));
    n = max(n, scalarsPtr->normInf());
    break;

  case NOX::Abstract::Vector::OneNorm:
    for (unsigned int i=0; i<vectorPtrs.size(); i++) 
      n += vectorPtrs[i]->norm(type);
    n += scalarsPtr->normOne();
    break;

  case NOX::Abstract::Vector::TwoNorm:
  default:
    for (unsigned int i=0; i<vectorPtrs.size(); i++) {
      nv = vectorPtrs[i]->norm(type);
      n += nv*nv;
    }
    nv = scalarsPtr->normFrobenius();
    n += nv*nv;
    n = sqrt(n);
   break;

  }

  return n;
}

double 
LOCA::Extended::Vector::norm(const NOX::Abstract::Vector& weights) const
{
  const LOCA::Extended::Vector& W = 
    dynamic_cast<const LOCA::Extended::Vector&>(weights);
  double n = 0.0;
  double nv;

  for (unsigned int i=0; i<vectorPtrs.size(); i++) {
    nv = vectorPtrs[i]->norm(*(W.vectorPtrs[i]));
    n += nv*nv;
  }
  for (int i=0; i<numScalars; i++) 
    n += ((*W.scalarsPtr)(i,0))*((*scalarsPtr)(i,0))*((*scalarsPtr)(i,0));
  n = sqrt(n);

  return n;
}

double 
LOCA::Extended::Vector::dot(const NOX::Abstract::Vector& y) const
{
  const LOCA::Extended::Vector& Y = 
    dynamic_cast<const LOCA::Extended::Vector&>(y);
  double d = 0.0;
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    d += vectorPtrs[i]->dot(*(Y.vectorPtrs[i]));
  for (int i=0; i<numScalars; i++)
    d += ((*scalarsPtr)(i,0))*((*Y.scalarsPtr)(i,0));

  return d;
}

int 
LOCA::Extended::Vector::length() const
{
  int l = 0;
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    l += vectorPtrs[i]->length();
  l += numScalars;

  return l;
}

void 
LOCA::Extended::Vector::print() const
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->print();
  scalarsPtr->print(cout);
  cout << endl;

}

void 
LOCA::Extended::Vector::setVector(int i, const NOX::Abstract::Vector& v)
{
  if (vectorPtrs[i] != NULL) 
    *(vectorPtrs[i]) = v;
  else
    vectorPtrs[i] = v.clone(NOX::DeepCopy);
  isView[i] = false;
}

void 
LOCA::Extended::Vector::setVectorView(int i, NOX::Abstract::Vector& v)
{
  if (vectorPtrs[i] != NULL) 
    if (!isView[i])
      delete vectorPtrs[i];
  vectorPtrs[i] = &v;
  isView[i] = true;
}

void 
LOCA::Extended::Vector::setScalar(int i, double s)
{
  (*scalarsPtr)(i,0) = s;
}

void
LOCA::Extended::Vector::setScalarArray(double *sv)
{
  delete scalarsPtr;
  scalarsPtr = new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View,
							   sv,
							   numScalars,
							   numScalars,
							   1);
}

const NOX::Abstract::Vector& 
LOCA::Extended::Vector::getVector(int i) const
{
  return *vectorPtrs[i];
}

NOX::Abstract::Vector& 
LOCA::Extended::Vector::getVector(int i)
{
  return *vectorPtrs[i];
}

double 
LOCA::Extended::Vector::getScalar(int i) const
{
  return (*scalarsPtr)(i,0);
}

double& 
LOCA::Extended::Vector::getScalar(int i)
{
  return (*scalarsPtr)(i,0);
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::Extended::Vector::getScalars() const
{
  return *scalarsPtr;
}

NOX::Abstract::MultiVector::DenseMatrix&
LOCA::Extended::Vector::getScalars()
{
  return *scalarsPtr;
}

LOCA::Extended::MultiVector*
LOCA::Extended::Vector::generateMultiVector(int nColumns, int nVectorRows, 
					    int nScalarRows) const
{
  return new LOCA::Extended::MultiVector(nColumns, nVectorRows, nScalarRows);
}
