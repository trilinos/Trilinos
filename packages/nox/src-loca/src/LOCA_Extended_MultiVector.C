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

#include "LOCA_Extended_MultiVector.H"
#include "LOCA_Extended_Vector.H"
#include "LOCA_ErrorCheck.H"

LOCA::Extended::MultiVector::MultiVector(
				   const LOCA::Extended::MultiVector& source,
				   NOX::CopyType type) :
  numColumns(source.numColumns),
  numMultiVecRows(source.numMultiVecRows),
  numScalarRows(source.numScalarRows),
  multiVectorPtrs(numMultiVecRows),
  scalarsPtr(NULL),
  extendedVectorPtrs(numColumns),
  isView(false)
{
  // Copy multi vecs
  for (int i=0; i<numMultiVecRows; i++) 
    multiVectorPtrs[i] = source.multiVectorPtrs[i]->clone(type);

  // Copy scalars
  scalarsPtr = new NOX::Abstract::MultiVector::DenseMatrix(*source.scalarsPtr);

  for (int i=0; i<numColumns; i++) {
    extendedVectorPtrs[i] = NULL;
  }
}

LOCA::Extended::MultiVector::MultiVector(
				   const LOCA::Extended::MultiVector& source,
				   int nColumns) :
  numColumns(nColumns),
  numMultiVecRows(source.numMultiVecRows),
  numScalarRows(source.numScalarRows),
  multiVectorPtrs(numMultiVecRows),
  scalarsPtr(NULL),
  extendedVectorPtrs(numColumns),
  isView(false)
{
  // Clone multivec blocks
  for (int i=0; i<numMultiVecRows; i++) 
    multiVectorPtrs[i] = source.multiVectorPtrs[i]->clone(numColumns);

  for (int i=0; i<numColumns; i++) {
    extendedVectorPtrs[i] = NULL;
  }

  scalarsPtr = new NOX::Abstract::MultiVector::DenseMatrix(numScalarRows,
							   numColumns);
}

LOCA::Extended::MultiVector::MultiVector(
				   const LOCA::Extended::MultiVector& source,
				   vector<int>& index, bool view) :
  numColumns(index.size()),
  numMultiVecRows(source.numMultiVecRows),
  numScalarRows(source.numScalarRows),
  multiVectorPtrs(numMultiVecRows),
  scalarsPtr(NULL),
  extendedVectorPtrs(numColumns),
  isView(view)
{
  // Check indices are valid
  for (unsigned int j=0; j<index.size(); j++) 
    source.checkIndex("LOCA::Extended::MultiVector()", index[j]);

  for (int i=0; i<numColumns; i++) {
    extendedVectorPtrs[i] = NULL;
  }

  // Check if indices are contiguous
  bool isCont = isContiguous(index);
  
  if (view) {
   
    // Copy multivectors
    for (int i=0; i<numMultiVecRows; i++) 
      multiVectorPtrs[i] = source.multiVectorPtrs[i]->subView(index);

    // Copy Scalars
    if (isCont) {
      double *vals = source.scalarsPtr->values() + 
	source.scalarsPtr->numRows()*index[0];
      scalarsPtr = 
	new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View,
						    vals,
						    numScalarRows,
						    numScalarRows,
						    numColumns);
    }
    else {
      LOCA::ErrorCheck::throwError(
		   "LOCA::Extended::MultiVector()",
		   "Sub-view with non-contiguous indices is not supported");
    }

  }
  else {
    
    // Copy multivectors
    for (int i=0; i<numMultiVecRows; i++) 
      multiVectorPtrs[i] = source.multiVectorPtrs[i]->subCopy(index);

    // Copy scalars
    if (isCont) {
      double *vals = source.scalarsPtr->values() + 
	source.scalarsPtr->numRows()*index[0];
      scalarsPtr = 
	new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::Copy,
						    vals,
						    numScalarRows,
						    numScalarRows,
						    numColumns);
    }
    else {
      scalarsPtr = 
	new NOX::Abstract::MultiVector::DenseMatrix(numScalarRows,
						    numColumns);
      for (int j=0; j<numColumns; j++)
	for (int i=0; i<numScalarRows; i++)
	  (*scalarsPtr)(i,j) = (*source.scalarsPtr)(i,index[j]);
    }
  }
}

LOCA::Extended::MultiVector::~MultiVector()
{
  for (int i=0; i<numMultiVecRows; i++)
    delete multiVectorPtrs[i];
  delete scalarsPtr;
  for (int i=0; i<numColumns; i++) {
    if (extendedVectorPtrs[i] != NULL) 
      delete extendedVectorPtrs[i];
  }
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::init(double gamma)
{
  // Initialize multivecs
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->init(gamma);

  // Initialize scalars
  scalarsPtr->putScalar(gamma);

  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::random(bool useSeed, int seed)
{
  // Fill multivecs with random values
  multiVectorPtrs[0]->random(useSeed, seed);
  for (int i=1; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->random();

  // Fill scalars with random values
  scalarsPtr->random();

  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::operator=(
				   const NOX::Abstract::MultiVector& source)
{
  return operator=(dynamic_cast<const LOCA::Extended::MultiVector&>(source));
}

LOCA::Extended::MultiVector& 
LOCA::Extended::MultiVector::operator=(
				   const LOCA::Extended::MultiVector& source)
{
  if (this != &source) {

    // Verify dimensions are consistent
    checkDimensions("LOCA::Extended::MultiVector::operator=()", source);

    // Copy multivecs
    for (int i=0; i<numMultiVecRows; i++)
      *(multiVectorPtrs[i]) = *(source.multiVectorPtrs[i]);

    // Copy scalars (don't use operator= since it destroys views)
    scalarsPtr->assign(*source.scalarsPtr);
  }

  return *this;
}

NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::setBlock(const NOX::Abstract::MultiVector& source,
				      vector<int>& index)
{
  return setBlock(dynamic_cast<const LOCA::Extended::MultiVector&>(source), 
		  index);
}

NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::setBlock(
				   const LOCA::Extended::MultiVector& source, 
				   vector<int>& index)
{
  // Verify dimensions are consistent
  if (source.numMultiVecRows != numMultiVecRows || 
      source.numScalarRows != numScalarRows) 
    LOCA::ErrorCheck::throwError("LOCA::Extended::MultiVector::setBlock()",
	"Size of supplied multivector is incompatible with this multivector");
  if (source.numColumns != index.size()) {
    LOCA::ErrorCheck::throwError("LOCA::Extended::MultiVector::setBlock()",
	"Size of supplied index vector is incompatible with this multivector");
  }

  // Set block in each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->setBlock(*(source.multiVectorPtrs[i]),index);

  // Set scalar vectors
  for (unsigned int j=0; j<index.size(); j++) {
    checkIndex("LOCA::Extended::MultiVector::augment()", index[j]);
    for (int i=0; i<numScalarRows; i++)
      (*scalarsPtr)(i,index[j]) = (*source.scalarsPtr)(i,j);
  }

  return *this;
}

NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::augment(const NOX::Abstract::MultiVector& source)
{
  return augment(dynamic_cast<const LOCA::Extended::MultiVector&>(source));
}

NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::augment(const LOCA::Extended::MultiVector& source)
{
  if (isView) {
    LOCA::ErrorCheck::throwError(
		   "LOCA::Extended::MultiVector::augment()",
		   "Augmenting a multivector view is not supported");
  }

  // Verify dimensions are consistent
  if (source.numMultiVecRows != numMultiVecRows ||
      source.numScalarRows != numScalarRows) 
    LOCA::ErrorCheck::throwError("LOCA::Extended::MultiVector::augment()",
	"Size of supplied multivector is incompatible with this multivector");

  // Augment each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->augment(*(source.multiVectorPtrs[i]));

  // Augment columns to scalar vectors
  scalarsPtr->reshape(numScalarRows, numColumns + source.numColumns);
  for (int j=0; j<source.numColumns; j++)
    for (int i=0; i<numScalarRows; i++)
      (*scalarsPtr)(i,numColumns + j) = (*source.scalarsPtr)(i,j);

  // Expand extended column vector ptrs
  extendedVectorPtrs.resize(numColumns + source.numColumns);
  for (int i=0; i<source.numColumns; i++)
    extendedVectorPtrs[numColumns + i] = NULL;

  // Update number of columns
  numColumns += source.numColumns;

  return *this;
}

NOX::Abstract::Vector&
LOCA::Extended::MultiVector::operator [] (int i)
{
  return getVector(i);
}

const NOX::Abstract::Vector&
LOCA::Extended::MultiVector::operator [] (int i) const
{
  return getVector(i);
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::update(double alpha, 
				    const NOX::Abstract::MultiVector& a, 
				    double gamma)
{
  return update(alpha, dynamic_cast<const LOCA::Extended::MultiVector&>(a), 
		gamma);
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::update(double alpha, 
				    const LOCA::Extended::MultiVector& a, 
				    double gamma)
{
  // Verify dimensions are consistent
  checkDimensions("LOCA::Extended::MultiVector::update()", a);

  // update each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->update(alpha, *(a.multiVectorPtrs[i]), gamma);

  // update scalars
  for (int j=0; j<numColumns; j++)
    for (int i=0; i<numScalarRows; i++)
      (*scalarsPtr)(i,j) = 
	gamma * (*scalarsPtr)(i,j) + alpha * (*a.scalarsPtr)(i,j);

  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::update(double alpha, 
				    const NOX::Abstract::MultiVector& a, 
				    double beta, 
				    const NOX::Abstract::MultiVector& b,
				    double gamma)
{
  return update(alpha, dynamic_cast<const LOCA::Extended::MultiVector&>(a),
		beta, dynamic_cast<const LOCA::Extended::MultiVector&>(b), 
		gamma);
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::update(double alpha, 
				    const LOCA::Extended::MultiVector& a, 
				    double beta, 
				    const LOCA::Extended::MultiVector& b,
				    double gamma)
{
  // Verify dimensions are consistent
  checkDimensions("LOCA::Extended::MultiVector::update()", a);
  checkDimensions("LOCA::Extended::MultiVector::update()", b);

  // update each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->update(alpha, *(a.multiVectorPtrs[i]), 
			       beta, *(b.multiVectorPtrs[i]),
			       gamma);

  // update scalars
  for (int j=0; j<numColumns; j++)
    for (int i=0; i<numScalarRows; i++)
      (*scalarsPtr)(i,j) = 
	gamma * (*scalarsPtr)(i,j) + alpha * (*a.scalarsPtr)(i,j) + 
	beta * (*b.scalarsPtr)(i,j);

  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::update(
			     Teuchos::ETransp transb, 
			     double alpha, 
			     const NOX::Abstract::MultiVector& a, 
			     const NOX::Abstract::MultiVector::DenseMatrix& b, 
			     double gamma)
{
  return update(transb, alpha, 
		dynamic_cast<const LOCA::Extended::MultiVector&>(a), 
		b, gamma);
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::update(
			    Teuchos::ETransp transb, 
			    double alpha, 
			    const LOCA::Extended::MultiVector& a, 
			    const NOX::Abstract::MultiVector::DenseMatrix& b, 
			    double gamma)
{
  // Verify dimensions are consistent
  if (a.numMultiVecRows != numMultiVecRows || 
      a.numScalarRows != numScalarRows) 
    LOCA::ErrorCheck::throwError("LOCA::Extended::MultiVector::update()",
      "Size of supplied multivector is incompatible with this multivector");

  if (transb == Teuchos::NO_TRANS) {
    if (a.numColumns != b.numRows() || numColumns != b.numCols()) 
      LOCA::ErrorCheck::throwError("LOCA::Extended::MultiVector::update()",
	    "Size of supplied matrix is incompatible with this multivector");
  }
  else
    if (a.numColumns != b.numCols() || numColumns != b.numRows()) 
      LOCA::ErrorCheck::throwError("LOCA::Extended::MultiVector::update()",
	    "Size of supplied matrix is incompatible with this multivector");


  // update each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->update(transb, alpha, *(a.multiVectorPtrs[i]), b, 
			       gamma);

  // update scalars
  scalarsPtr->multiply(Teuchos::NO_TRANS, transb, alpha, *(a.scalarsPtr), 
		       b, gamma);

  return *this;
}

NOX::Abstract::MultiVector* 
LOCA::Extended::MultiVector::clone(NOX::CopyType type) const
{
  return new LOCA::Extended::MultiVector(*this, type);
}

NOX::Abstract::MultiVector* 
LOCA::Extended::MultiVector::clone(int numvecs) const
{
  return new LOCA::Extended::MultiVector(*this, numvecs);
}

NOX::Abstract::MultiVector* 
LOCA::Extended::MultiVector::subCopy(vector<int>& index) const
{
  return new LOCA::Extended::MultiVector(*this, index, false);
}

NOX::Abstract::MultiVector* 
LOCA::Extended::MultiVector::subView(vector<int>& index) const
{
  return new LOCA::Extended::MultiVector(*this, index, true);
}

void
LOCA::Extended::MultiVector::norm(vector<double>& result,
		       NOX::Abstract::Vector::NormType type) const
{

  // Make sure result vector is of appropriate size
  if (result.size() != numColumns)
    result.resize(numColumns);

  // Zero out result vector
  for (int i=0; i<numColumns; i++)
    result[i] = 0.0;

  // Intermediate vector to hold norms of a multivector
  vector<double> vecNorm(result);

  switch (type) {

  case NOX::Abstract::Vector::MaxNorm:
    
    // compute max norm of all multivecs
    for (int i=0; i<numMultiVecRows; i++) {
      multiVectorPtrs[i]->norm(vecNorm, type);
      for (int j=0; j<numColumns; j++) {
	if (result[j] < vecNorm[j])
	  result[j] = vecNorm[j];
      }
    }

    // compute max norm of each scalar vector column
    for (int j=0; j<numColumns; j++) {
      for (int i=0; i<numScalarRows; i++)
	if (result[j] < (*scalarsPtr)(i,j))
	  result[j] = (*scalarsPtr)(i,j);
    }
    break;

  case NOX::Abstract::Vector::OneNorm:

    // compute one-norm of all multivecs
    for (int i=0; i<numMultiVecRows; i++) {
      multiVectorPtrs[i]->norm(vecNorm, type);
      for (int j=0; j<numColumns; j++) {
	result[j] += vecNorm[j];
      }
    }

    // compute one-norm of each scalar vector column
    for (int j=0; j<numColumns; j++) {
      for (int i=0; i<numScalarRows; i++)
	result[j] += fabs((*scalarsPtr)(i,j));
    }
    break;

  case NOX::Abstract::Vector::TwoNorm:
  default:

    // compute two-norm of all multivecs
    for (int i=0; i<numMultiVecRows; i++) {
      multiVectorPtrs[i]->norm(vecNorm, type);
      for (int j=0; j<numColumns; j++) {
	result[j] += vecNorm[j]*vecNorm[j];
      }
    }

    // compute two-norm of each scalar vector column
    for (int j=0; j<numColumns; j++) {
      for (int i=0; i<numScalarRows; i++)
	result[j] += (*scalarsPtr)(i,j) * (*scalarsPtr)(i,j);
      result[j] = sqrt(result[j]);
    }
    break;
  }
}

void 
LOCA::Extended::MultiVector::multiply(
			    double alpha, 
			    const NOX::Abstract::MultiVector& y,
			    NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  multiply(alpha, dynamic_cast<const LOCA::Extended::MultiVector&>(y), b);
}

void 
LOCA::Extended::MultiVector::multiply(
			    double alpha, 
			    const LOCA::Extended::MultiVector& y,
			    NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  // Verify dimensions are consistent
  if (y.numMultiVecRows != numMultiVecRows || y.numColumns != b.numRows() ||
      y.numScalarRows != numScalarRows || numColumns != b.numCols()) 
    LOCA::ErrorCheck::throwError("LOCA::Extended::MultiVector::multiply()",
      "Size of supplied multivector/matrix is incompatible with this multivector");

  // Zero out b
  b.putScalar(0.0);

  // Create temporary matrix to hold product for each multivec
  NOX::Abstract::MultiVector::DenseMatrix tmp(b);

  // Compute and sum products for each multivec
  for (int i=0; i<numMultiVecRows; i++) {
    multiVectorPtrs[i]->multiply(alpha, *(y.multiVectorPtrs[i]), tmp);
    b += tmp;
  }

  // Compute and add in product for scalars
  b.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, *y.scalarsPtr,
	     *scalarsPtr, 1.0);

}

int 
LOCA::Extended::MultiVector::length() const
{
  int len = 0;

  // Sum lengths of all multivec block rows
  for (int i=0; i<numMultiVecRows; i++)
    len += multiVectorPtrs[i]->length();

  // Add in number of scalar rows
  len += numScalarRows;

  return len;
}

int 
LOCA::Extended::MultiVector::numVectors() const
{
  return numColumns;
}

void 
LOCA::Extended::MultiVector::print() const
{
  // Print multi-vecs
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->print();

  // Print scalars
  scalarsPtr->print(cout);
}

const NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::getMultiVector(int i) const
{
  checkIndex("LOCA::Extended::MultiVector::getMultiVector()", i);

  return *(multiVectorPtrs[i]);
}

NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::getMultiVector(int i)
{
  checkIndex("LOCA::Extended::MultiVector::getMultiVector()", i);

  return *(multiVectorPtrs[i]);
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::Extended::MultiVector::getScalars() const 
{
  return *scalarsPtr;
}

NOX::Abstract::MultiVector::DenseMatrix&
LOCA::Extended::MultiVector::getScalars()
{
  return *scalarsPtr;
}

const double&
LOCA::Extended::MultiVector::getScalar(int i, int j) const
{
  checkIndex("LOCA::Extended::MultiVector::getScalar()", i, j);

  return (*scalarsPtr)(i,j);
}

double&
LOCA::Extended::MultiVector::getScalar(int i, int j)
{
  checkIndex("LOCA::Extended::MultiVector::getScalar()", i, j);

  return (*scalarsPtr)(i,j);
}

LOCA::Extended::Vector&
LOCA::Extended::MultiVector::getVector(int i) 
{
  // Verify index is valid
  checkIndex("LOCA::Extended::MultiVector::vector()", i);

  // Create extended vector view if necessary
  if (extendedVectorPtrs[i] == NULL) {
    extendedVectorPtrs[i] = generateVector(numMultiVecRows, numScalarRows);
    for (int k=0; k<numMultiVecRows; k++) 
      extendedVectorPtrs[i]->setVectorView(k, (*multiVectorPtrs[k])[i]);
    extendedVectorPtrs[i]->setScalarArray((*scalarsPtr)[i]);
  }

  return *(extendedVectorPtrs[i]);
}

const LOCA::Extended::Vector&
LOCA::Extended::MultiVector::getVector(int i) const
{
  // Verify index is valid
  checkIndex("LOCA::Extended::MultiVector::vector()", i);

  // Create extended vector view if necessary
  if (extendedVectorPtrs[i] == NULL) {
    extendedVectorPtrs[i] = generateVector(numMultiVecRows, numScalarRows);
    for (int k=0; k<numMultiVecRows; k++) 
      extendedVectorPtrs[i]->setVectorView(k, (*multiVectorPtrs[k])[i]);
    extendedVectorPtrs[i]->setScalarArray((*scalarsPtr)[i]);
  }

  return *(extendedVectorPtrs[i]);
}

int
LOCA::Extended::MultiVector::getNumScalarRows() const
{
  return numScalarRows;
}

int
LOCA::Extended::MultiVector::getNumMultiVectors() const
{
  return numMultiVecRows;
}

LOCA::Extended::MultiVector::MultiVector(int nColumns, int nVectorRows,
					 int nScalarRows) :
  numColumns(nColumns),
  numMultiVecRows(nVectorRows),
  numScalarRows(nScalarRows),
  multiVectorPtrs(numMultiVecRows),
  scalarsPtr(NULL),
  extendedVectorPtrs(numColumns),
  isView(false)
{
  for (int i=0; i<numColumns; i++) {
    extendedVectorPtrs[i] = NULL;
  }

  scalarsPtr = new NOX::Abstract::MultiVector::DenseMatrix(numScalarRows, 
							   numColumns);
}

LOCA::Extended::Vector* 
LOCA::Extended::MultiVector::generateVector(int nVecs, int nScalarRows) const
{
  return new LOCA::Extended::Vector(nVecs, nScalarRows);
}

void
LOCA::Extended::MultiVector::setMultiVectorPtr(
					 int i, 
					 NOX::Abstract::MultiVector* v)
{
  checkIndex("LOCA::Extended::MultiVector::setMultiVectorPtr()", i);

  multiVectorPtrs[i] = v;
}

void
LOCA::Extended::MultiVector::checkDimensions(
				  const string& callingFunction,
				  const LOCA::Extended::MultiVector& a) const
{
  if (a.numMultiVecRows != numMultiVecRows || a.numColumns != numColumns ||
      a.numScalarRows != numScalarRows)
    LOCA::ErrorCheck::throwError(callingFunction,
      "Size of supplied multivector is incompatible with this multivector");
}

void 
LOCA::Extended::MultiVector::checkIndex(const string& callingFunction,
					int i) const 
{
  if ( i < 0 || i >= numColumns ) 
    LOCA::ErrorCheck::throwError(callingFunction, "Invalid column index");
}

void 
LOCA::Extended::MultiVector::checkIndex(const string& callingFunction,
					int i, int j) const 
{
  if ( i < 0 || i >= numScalarRows ) 
    LOCA::ErrorCheck::throwError(callingFunction, "Invalid row index");
  if ( j < 0 || j >= numColumns ) 
    LOCA::ErrorCheck::throwError(callingFunction, "Invalid column index");
}

bool
LOCA::Extended::MultiVector::isContiguous(const vector<int>& index) const 
{
  for (unsigned int i=0; i<index.size(); i++) {
    if (index[i] != index[0] + i)
      return false;
  }
  return true;
}
