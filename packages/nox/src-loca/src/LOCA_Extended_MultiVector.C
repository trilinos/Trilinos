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

LOCA::Extended::MultiVector::MultiVector(int nColumns, int nVecs, 
					 int nScalarRows, bool view) :
  numColumns(nColumns),
  numMultiVecRows(nVecs),
  numScalarRows(nScalarRows),
  multiVectorPtrs(nVecs),
  scalarVectorPtrs(nColumns),
  extendedVectorPtrs(nColumns),
  isView(view)
{
  for (int i=0; i<numColumns; i++) {
    if (!isView) {
      scalarVectorPtrs[i] = 
	new NOX::Abstract::MultiVector::DenseMatrix(numScalarRows,1);
    }
    extendedVectorPtrs[i] = NULL;
  }
}

LOCA::Extended::MultiVector::MultiVector(
				   const LOCA::Extended::MultiVector& source,
				   NOX::CopyType type) :
  numColumns(source.numColumns),
  numMultiVecRows(source.numMultiVecRows),
  numScalarRows(source.numScalarRows),
  multiVectorPtrs(numMultiVecRows),
  scalarVectorPtrs(numColumns),
  extendedVectorPtrs(numColumns),
  isView(false)
{
  // Copy multi vecs
  for (int i=0; i<numMultiVecRows; i++) 
    multiVectorPtrs[i] = source.multiVectorPtrs[i]->clone(type);

  // Copy scalars
  for (int i=0; i<numColumns; i++) {
    scalarVectorPtrs[i] = 
      new NOX::Abstract::MultiVector::DenseMatrix(*(source.scalarVectorPtrs[i]));
    extendedVectorPtrs[i] = NULL;
  }
}

LOCA::Extended::MultiVector::~MultiVector()
{
  for (int i=0; i<numMultiVecRows; i++)
    delete multiVectorPtrs[i];
  for (int i=0; i<numColumns; i++) {
    if (!isView) 
      delete scalarVectorPtrs[i];
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
  for (int i=0; i<numColumns; i++)
    scalarVectorPtrs[i]->putScalar(gamma);

  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::random(bool useSeed, int seed)
{
  // Fill multivecs with random values
  multiVectorPtrs[0]->random(useSeed, seed);
  for (int i=1; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->random();

  // Fill scalar vectors with random values
  for (int i=0; i<numColumns; i++)
    scalarVectorPtrs[i]->random();

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
    checkDimensions(source);

    // Copy multivecs
    for (int i=0; i<numMultiVecRows; i++)
      *(multiVectorPtrs[i]) = *(source.multiVectorPtrs[i]);

    // Copy scalar vectors
    for (int i=0; i<numColumns; i++)
      *(scalarVectorPtrs[i]) = *(source.scalarVectorPtrs[i]);
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
  checkDimensions(source, index);

  // Set block in each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->setBlock(*(source.multiVectorPtrs[i]),index);

  // Set scalar vectors
  for (unsigned int j=0; j<index.size(); j++) {
    *(scalarVectorPtrs[index[j]]) = *(source.scalarVectorPtrs[j]);
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
  // Verify dimensions are consistent
  checkAugmentDimensions(source);

  // Augment each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->augment(*(source.multiVectorPtrs[i]));

  // Augment columns to scalar vectors
  scalarVectorPtrs.resize(numColumns + source.numColumns);
  for (int i=0; i<source.numColumns; i++)
    scalarVectorPtrs[numColumns + i] = 
      new NOX::Abstract::MultiVector::DenseMatrix(*(source.scalarVectorPtrs[i]));

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
  // Verify index is valid
  checkIndex(i);

  // Create extended vector view if necessary
  if (extendedVectorPtrs[i] == NULL) {
    extendedVectorPtrs[i] = 
      new LOCA::Extended::Vector(numMultiVecRows, numScalarRows);
    for (int k=0; k<numMultiVecRows; k++) 
      extendedVectorPtrs[i]->setVectorView(k, multiVectorPtrs[k]->operator[](i));
    for (int k=0; k<numScalarRows; k++)
      extendedVectorPtrs[i]->setScalar(k,scalarVectorPtrs[i]->operator()(k,1));
  }

  return *(extendedVectorPtrs[i]);
}

const NOX::Abstract::Vector&
LOCA::Extended::MultiVector::operator [] (int i) const
{
  // Verify index is valid
  checkIndex(i);

  // Create extended vector view if necessary
  if (extendedVectorPtrs[i] == NULL) {
    extendedVectorPtrs[i] = 
      new LOCA::Extended::Vector(numMultiVecRows, numScalarRows);
    for (int k=0; k<numMultiVecRows; k++) 
      extendedVectorPtrs[i]->setVectorView(k, multiVectorPtrs[k]->operator[](i));
    for (int k=0; k<numScalarRows; k++)
      extendedVectorPtrs[i]->setScalar(k,scalarVectorPtrs[i]->operator()(k,1));
  }

  return *(extendedVectorPtrs[i]);
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
  checkDimensions(a);

  // update each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->update(alpha, *(a.multiVectorPtrs[i]), gamma);

  // update each scalar vector
  for (int j=0; j<numColumns; j++) {
    for (int i=0; i<numScalarRows; i++)
      scalarVectorPtrs[j]->operator()(i,1) = 
	gamma *   scalarVectorPtrs[j]->operator()(i,1) + 
	alpha * a.scalarVectorPtrs[j]->operator()(i,1);
  }

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
  checkDimensions(a);
  checkDimensions(b);

  // update each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->update(alpha, *(a.multiVectorPtrs[i]), 
			       beta, *(b.multiVectorPtrs[i]),
			       gamma);

  // update each scalar vector
  for (int j=0; j<numColumns; j++) {
    for (int i=0; i<numScalarRows; i++)
      scalarVectorPtrs[j]->operator()(i,1) = 
	gamma *   scalarVectorPtrs[j]->operator()(i,1) + 
	alpha * a.scalarVectorPtrs[j]->operator()(i,1) +
	beta  * b.scalarVectorPtrs[j]->operator()(i,1);
  }

  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::update(
			     double alpha, 
			     const NOX::Abstract::MultiVector& a, 
			     const NOX::Abstract::MultiVector::DenseMatrix& b, 
			     double gamma)
{
  return update(alpha, dynamic_cast<const LOCA::Extended::MultiVector&>(a), 
		b, gamma);
}

NOX::Abstract::MultiVector& 
LOCA::Extended::MultiVector::update(
			    double alpha, 
			    const LOCA::Extended::MultiVector& a, 
			    const NOX::Abstract::MultiVector::DenseMatrix& b, 
			    double gamma)
{
  // Verify dimensions are consistent
  checkDimensions(a, b);

  // update each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->update(alpha, *(a.multiVectorPtrs[i]), b, gamma);

  // update scalar vectors
  double tmp;
  for (int j=0; j<numColumns; j++) {
    for (int i=0; i<numScalarRows; i++) {
      tmp = 0.0;
      for (int k=0; k<a.numColumns; k++)
	tmp += a.scalarVectorPtrs[k]->operator()(i,1) * b(k,j);
      scalarVectorPtrs[j]->operator()(i,1) = 
	gamma * scalarVectorPtrs[j]->operator()(i,1) + alpha * tmp;
    }
  }

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
  // Create empty extended multi-vec of appropriate size
  LOCA::Extended::MultiVector* tmp = 
    new LOCA::Extended::MultiVector(numvecs, numMultiVecRows, numScalarRows); 
  
  // Clone multivec blocks
  for (int i=0; i<numMultiVecRows; i++) 
    tmp->multiVectorPtrs[i] = multiVectorPtrs[i]->clone(numvecs);
  
  return tmp;
}

NOX::Abstract::MultiVector* 
LOCA::Extended::MultiVector::subCopy(vector<int>& index) const
{
  int numvecs = index.size();

  // Verify dimensions are consistent
  checkDimensions(*this, index);

  // Create extended multi-vec of appropriate size
  LOCA::Extended::MultiVector* tmp = 
    new LOCA::Extended::MultiVector(numvecs, numMultiVecRows, numScalarRows); 
  
  // Clone multivec blocks
  for (int i=0; i<numMultiVecRows; i++) 
    tmp->multiVectorPtrs[i] = multiVectorPtrs[i]->subCopy(index);

  // Copy scalar vectors
  for (int i=0; i<numvecs; i++)
    *(tmp->scalarVectorPtrs[i]) = *(scalarVectorPtrs[index[i]]);

  return tmp;
}

NOX::Abstract::MultiVector* 
LOCA::Extended::MultiVector::subView(vector<int>& index) const
{
  int numvecs = index.size();

  // Verify dimensions are consistent
  checkDimensions(*this, index);

  // Create extended multi-vec of appropriate size
  LOCA::Extended::MultiVector* tmp = 
    new LOCA::Extended::MultiVector(numvecs, numMultiVecRows, numScalarRows,
				    true); 
  
  // Clone multivec blocks
  for (int i=0; i<numMultiVecRows; i++) 
    tmp->multiVectorPtrs[i] = multiVectorPtrs[i]->subView(index);

  // Copy scalar vector views
  for (int i=0; i<numvecs; i++)
    tmp->scalarVectorPtrs[i] = scalarVectorPtrs[index[i]];

  return tmp;
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
  double tmp;

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
    for (int i=0; i<numColumns; i++) {
      tmp = scalarVectorPtrs[i]->normInf();
      if (result[i] < tmp)
	result[i] = tmp;
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
    for (int i=0; i<numColumns; i++) {
      result[i] += scalarVectorPtrs[i]->normOne();
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
    for (int i=0; i<numColumns; i++) {
      tmp = scalarVectorPtrs[i]->normFrobenius();
      result[i] += tmp*tmp;
      result[i] = sqrt(result[i]);
    }
    break;
  }
}

void 
LOCA::Extended::MultiVector::dot(
			    double alpha, 
			    const NOX::Abstract::MultiVector& y,
			    NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  dot(alpha, dynamic_cast<const LOCA::Extended::MultiVector&>(y), b);
}

void 
LOCA::Extended::MultiVector::dot(
			    double alpha, 
			    const LOCA::Extended::MultiVector& y,
			    NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  // Zero out b
  b.putScalar(0.0);

  // Create temporary matrix to hold dot product for each multivec
  NOX::Abstract::MultiVector::DenseMatrix tmp(b);

  // Compute and sum dot products for each multivec
  for (int i=0; i<numMultiVecRows; i++) {
    multiVectorPtrs[i]->dot(alpha, *(y.multiVectorPtrs[i]), tmp);
    b += tmp;
  }

  // Compute and add in dot product for scalars
  for (int i=0; i<b.numRows(); i++)
    for (int j=0; j<b.numCols(); j++)
      for (int k=0; k<numScalarRows; k++)
	b(i,j) += scalarVectorPtrs[i]->operator()(k,1) * y.scalarVectorPtrs[j]->operator()(k,1);

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
  for (int i=0; i<numScalarRows; i++) {
    for (int j = 0; j<numColumns; j++)
      cout << scalarVectorPtrs[j]->operator()(i,1) << " ";
    cout << endl;
  }
}

void
LOCA::Extended::MultiVector::setMultiVector(
					 int i, 
					 const NOX::Abstract::MultiVector& v)
{
  checkIndex(i);

  multiVectorPtrs[i] = v.clone(NOX::DeepCopy);
}

void
LOCA::Extended::MultiVector::setMultiVectorPtr(
					 int i, 
					 NOX::Abstract::MultiVector* v)
{
  checkIndex(i);

  multiVectorPtrs[i] = v;
}

void
LOCA::Extended::MultiVector::setScalar(int i, int j, double s) 
{
  checkIndex(i,j);

  scalarVectorPtrs[j]->operator()(i,1) = s;
}

const NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::getMultiVector(int i) const
{
  checkIndex(i);

  return *(multiVectorPtrs[i]);
}

NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::getMultiVector(int i)
{
  checkIndex(i);

  return *(multiVectorPtrs[i]);
}

const double&
LOCA::Extended::MultiVector::getScalar(int i, int j) const
{
  checkIndex(i, j);

  return scalarVectorPtrs[j]->operator()(i,1);
}

double&
LOCA::Extended::MultiVector::getScalar(int i, int j)
{
  checkIndex(i, j);

  return scalarVectorPtrs[j]->operator()(i,1);
}

void
LOCA::Extended::MultiVector::checkDimensions(
				  const LOCA::Extended::MultiVector& a) const
{
  if (a.numMultiVecRows != numMultiVecRows || a.numColumns != numColumns ||
      a.numScalarRows != numScalarRows) {
    cerr << "LOCA::Extended::MultiVector:  Error!  Size of supplied "
	 << " multivector is incompatible with this multivector" << endl;
    throw "NOX Error";
  }
}

void
LOCA::Extended::MultiVector::checkDimensions(
					const LOCA::Extended::MultiVector& a,
					vector<int>& index) const
{
  if (a.numMultiVecRows != numMultiVecRows || 
      a.numScalarRows != numScalarRows) {
    cerr << "LOCA::Extended::MultiVector:  Error!  Size of supplied "
	 << " multivector is incompatible with this multivector" << endl;
    throw "NOX Error";
  }
  if (a.numColumns != index.size()) {
    cerr << "LOCA::Extended::MultiVector:  Error!  Size of supplied "
	 << " multivector is incompatible supplied index vector" << endl;
    throw "NOX Error";
  }
  for (unsigned int i=0; i<index.size(); i++) {
    if (index[i] <= 0 || index[i] >= numColumns) {
      cerr << "LOCA::Extended::MultiVector:  Error!  Supplied index" 
	   << " vector is invalid" << endl;
      throw "NOX Error";
    }
  }
}

void
LOCA::Extended::MultiVector::checkDimensions(
		       const LOCA::Extended::MultiVector& a,
		       const NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  if (a.numMultiVecRows != numMultiVecRows || a.numColumns != b.numRows() ||
      a.numScalarRows != numScalarRows || numColumns != b.numCols()) {
    cerr << "LOCA::Extended::MultiVector:  Error!  Size of supplied "
	 << " multivector/matrix is incompatible with this multivector" 
	 << endl;
    throw "NOX Error";
  }
}

void
LOCA::Extended::MultiVector::checkAugmentDimensions(
				  const LOCA::Extended::MultiVector& a) const
{
  if (a.numMultiVecRows != numMultiVecRows ||
      a.numScalarRows != numScalarRows) {
    cerr << "LOCA::Extended::MultiVector:  Error!  Size of supplied "
	 << " multivector is incompatible with this multivector" << endl;
    throw "NOX Error";
  }
}

void 
LOCA::Extended::MultiVector::checkIndex(int i) const 
{
  if ( i < 0 || i >= numColumns ) {
    cerr << "LOCA::Extended::MultiVector:  Error!  Invalid column index " 
	 << i << endl;
    throw "NOX Error";
  }
}

void 
LOCA::Extended::MultiVector::checkIndex(int i, int j) const 
{
  if ( i < 0 || i >= numColumns ) {
    cerr << "LOCA::Extended::MultiVector:  Error!  Invalid column index " 
	 << i << endl;
    throw "NOX Error";
  }
  if ( j < 0 || j >= numScalarRows ) {
    cerr << "LOCA::Extended::MultiVector:  Error!  Invalid scalar row  index " 
	 << j << endl;
    throw "NOX Error";
  }
}
