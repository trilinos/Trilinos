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

#include "NOX_LAPACK_MultiVector.H"
#include "NOX_LAPACK_Vector.H"
#include "NOX_BLAS_Wrappers.H"
#include "NOX_Random.H"

NOX::LAPACK::MultiVector::MultiVector() {}

NOX::LAPACK::MultiVector::MultiVector(int n, int nVecs) :
  x(n*nVecs), 
  vecPtrs(nVecs),
  noxVectors(nVecs),
  len(n),
  numCols(nVecs),
  ownsCols(true)
{
  checkDimensions();
  for (int i=0; i<nVecs; i++) {
    vecPtrs[i] = &(x[n*i]);
    noxVectors[i] = NULL;
  }
}

NOX::LAPACK::MultiVector::MultiVector(int n, int nVecs, double *v) :
  x(n*nVecs), 
  vecPtrs(nVecs),
  noxVectors(nVecs),
  len(n),
  numCols(nVecs),
  ownsCols(true)
{
  checkDimensions();
  for (int i=0; i<n*nVecs; i++)
    x[i] = v[i];
  for (int i=0; i<nVecs; i++) {
    vecPtrs[i] = &(x[n*i]);
    noxVectors[i] = NULL;
  }
}

NOX::LAPACK::MultiVector::MultiVector(int n, int nVecs, double **v) :
  x(n*nVecs), 
  vecPtrs(nVecs),
  noxVectors(nVecs),
  len(n),
  numCols(nVecs),
  ownsCols(true)
{
  checkDimensions();
  for (int j=0; j<nVecs; j++)
    for (int i=0; i<n; i++)
      x[n*j+i] = v[j][i];
  for (int i=0; i<nVecs; i++) {
    vecPtrs[i] = &(x[n*i]);
    noxVectors[i] = NULL;
  }
}

NOX::LAPACK::MultiVector::MultiVector(const NOX::LAPACK::MultiVector& source,
				      NOX::CopyType type) :
  x(source.numCols * source.len), 
  vecPtrs(source.numCols),
  noxVectors(source.numCols),
  len(source.len),
  numCols(source.numCols),
  ownsCols(true)
{
  for (int j=0; j<numCols; j++) { 
    for (int i=0; i<len; i++) 
      x[j*len+i] = source.vecPtrs[j][i];
    vecPtrs[j] = &(x[len*j]);
    noxVectors[j] = NULL;
  }
}

NOX::LAPACK::MultiVector::~MultiVector()
{
  for (unsigned int i=0; i<noxVectors.size(); i++)
    if (noxVectors[i] != NULL)
      delete noxVectors[i];
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::init(double value)
{
  for (int j=0; j<numCols; j++)
    for (int i=0; i<len; i++)
      vecPtrs[j][i] = value;

  return *this;
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::random(bool useSeed, int seed)
{
  if (useSeed)
    NOX::Random::setSeed(seed);

  for (int j=0; j<numCols; j++)
    for (int i=0; i<len; i++)
      vecPtrs[j][i] = NOX::Random::number();

  return *this;
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::operator=(const NOX::LAPACK::MultiVector& source)
{
  if (this != &source) {
    checkSize(source.numCols);

    for (int j=0; j<numCols; j++)
      for (int i=0; i<len; i++)
	vecPtrs[j][i] = source.vecPtrs[j][i];
  }

  return *this;
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::operator=(const NOX::Abstract::MultiVector& source)
{
  return operator=(dynamic_cast<const NOX::LAPACK::MultiVector&>(source));
}

NOX::Abstract::MultiVector&
NOX::LAPACK::MultiVector::setBlock(const NOX::Abstract::MultiVector& source, 
				   vector<int>& index)
{
  return setBlock(dynamic_cast<const NOX::LAPACK::MultiVector&>(source),
		  index);
}

NOX::Abstract::MultiVector&
NOX::LAPACK::MultiVector::setBlock(const NOX::LAPACK::MultiVector& source, 
				   vector<int>& index)
{
  int ind;

  source.checkSize(index.size()-1);
  for (unsigned int j=0; j<index.size(); j++) {
    ind = index[j];
    checkIndex(ind);
    for (int i=0; i<len; i++)
      vecPtrs[ind][i] = source.vecPtrs[j][i];
  }

  return *this;
}

NOX::Abstract::MultiVector&
NOX::LAPACK::MultiVector::augment(const NOX::Abstract::MultiVector& source)
{
  return augment(dynamic_cast<const NOX::LAPACK::MultiVector&>(source));
}

NOX::Abstract::MultiVector&
NOX::LAPACK::MultiVector::augment(const NOX::LAPACK::MultiVector& source)
{
  if (!ownsCols) {
    cerr << "NOX::LAPACK::MultiVector:  Error!  Cannot call augment()" << endl 
	 << "on a MultiVector created as a View" << endl;
    throw "NOX::LAPACK Error";
  }

  int newSize = numCols + source.numCols;
  x.resize(newSize*len);
  vecPtrs.resize(newSize);
  for (int i=0; i<source.numCols*len; i++)
    x[numCols*len+i] = source.x[i];
  for (int i=0; i<newSize; i++)
    vecPtrs[i] = &(x[len*i]);
  numCols = newSize;

  return *this;
}

NOX::Abstract::Vector&
NOX::LAPACK::MultiVector::operator [] (int i)
{
  if ( i < 0 || i > numCols ) {
    cerr << "NOX::LAPACK::MultiVector::operator[]:  Error!  Invalid index " 
	 << i << endl;
    throw "NOX::LAPACK Error";
  }
  if (noxVectors[i] == NULL) {
    noxVectors[i] = new NOX::LAPACK::Vector(len, vecPtrs[i], true);
  }
  return *(noxVectors[i]);
}

const NOX::Abstract::Vector&
NOX::LAPACK::MultiVector::operator [] (int i) const
{
  if ( i < 0 || i > numCols ) {
    cerr << "NOX::LAPACK::MultiVector::operator[]:  Error!  Invalid index " 
	 << i << endl;
    throw "NOX::LAPACK Error";
  }
  if (noxVectors[i] == NULL) {
    noxVectors[i] = new NOX::LAPACK::Vector(len, vecPtrs[i], true);
  }
  return *(noxVectors[i]);
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::update(double alpha, 
				 const NOX::Abstract::MultiVector& a, 
				 double gamma)
{
  return update(alpha, dynamic_cast<const NOX::LAPACK::MultiVector&>(a),
		gamma);
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::update(double alpha, 
				 const NOX::LAPACK::MultiVector& a, 
				 double gamma)
{
  checkSize(a.numCols);
  for (int j=0; j<numCols; j++)
    for (int i=0; i<len; i++)
      vecPtrs[j][i] = gamma*vecPtrs[j][i] + alpha * a.vecPtrs[j][i];

  return *this;
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::update(double alpha, 
				 const NOX::Abstract::MultiVector& a, 
				 double beta, 
				 const NOX::Abstract::MultiVector& b,
				 double gamma)
{
  return update(alpha, dynamic_cast<const NOX::LAPACK::MultiVector&>(a),
		beta, dynamic_cast<const NOX::LAPACK::MultiVector&>(b), gamma);
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::update(double alpha, 
				 const NOX::LAPACK::MultiVector& a, 
				 double beta, 
				 const NOX::LAPACK::MultiVector& b,
				 double gamma)
{
  checkSize(a.numCols);
  checkSize(b.numCols);
  for (int j=0; j<numCols; j++)
    for (int i=0; i<len; i++)
      vecPtrs[j][i] = gamma*vecPtrs[j][i] + alpha * a.vecPtrs[j][i] + 
	beta * b.vecPtrs[j][i];

  return *this;
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::update(
			    double alpha, 
			    const NOX::Abstract::MultiVector& a, 
			    const NOX::Abstract::MultiVector::DenseMatrix& b, 
			    double gamma)
{
  return update(alpha, dynamic_cast<const NOX::LAPACK::MultiVector&>(a), b,
		gamma);
}

NOX::Abstract::MultiVector& 
NOX::LAPACK::MultiVector::update(
			    double alpha, 
			    const NOX::LAPACK::MultiVector& a, 
			    const NOX::Abstract::MultiVector::DenseMatrix& b, 
			    double gamma)
{
  int b_rows = b.numRows();
  int b_cols = b.numCols();
  a.checkSize(b_rows);
  checkSize(b_cols);

  // Create strided versions of a and c (this) if they are not already strided
  // Note:  This is not the best test.  A multivector could not own its columns
  // and still be strided, i.e., if the indices in the view were contiguous
  NOX::LAPACK::MultiVector *sa;
  NOX::LAPACK::MultiVector *sc;
  if (a.ownsCols)
    sa = const_cast<NOX::LAPACK::MultiVector*>(&a);
  else
    sa = new NOX::LAPACK::MultiVector(a);

  if (ownsCols)
    sc = this;
  else
    sc = new NOX::LAPACK::MultiVector(*this);
 
  // Get underlying arrays
  double *av = &(sa->x[0]);
  double *bv = b.values();
  double *cv = &(sc->x[0]);

  // Call BLAS routine
  DGEMM_F77("N", "N", &len, &numCols, &b_cols, &alpha, av, &len, bv, &b_rows,
	    &gamma, cv, &len);

  // Copy solution into this and delete temporaries
  if (!ownsCols) {
    *this = *sc;
    delete sc;
  }
  if (!a.ownsCols)
    delete sa;

  return *this;
}

NOX::Abstract::MultiVector* 
NOX::LAPACK::MultiVector::clone(CopyType type) const
{
  return new NOX::LAPACK::MultiVector(*this, type);
}

NOX::Abstract::MultiVector* 
NOX::LAPACK::MultiVector::clone(int numvecs) const
{
  return new NOX::LAPACK::MultiVector(len, numvecs);
}

NOX::Abstract::MultiVector* 
NOX::LAPACK::MultiVector::subCopy(vector<int>& index) const
{
  NOX::LAPACK::MultiVector* tmp = 
    new NOX::LAPACK::MultiVector(len, index.size());
  int ind;
  for (unsigned int j=0; j<index.size(); j++) {
    ind = index[j];
    checkIndex(ind);
    for (int i=0; i<len; i++)
      tmp->vecPtrs[j][i] = vecPtrs[ind][i];
  }

  return tmp;
}

NOX::Abstract::MultiVector* 
NOX::LAPACK::MultiVector::subView(vector<int>& index) const
{
  checkSize(index.size());
  NOX::LAPACK::MultiVector* tmp = new NOX::LAPACK::MultiVector;
  tmp->vecPtrs.resize(index.size());
  tmp->len = len;
  tmp->numCols = index.size();
  tmp->ownsCols = false;
  int ind;
  for (unsigned int j=0; j<index.size(); j++) {
    ind = index[j];
    checkIndex(ind);
    tmp->vecPtrs[j] = vecPtrs[ind];
  }

  return tmp;
}

void
NOX::LAPACK::MultiVector::norm(vector<double>& result, 
			       NOX::Abstract::Vector::NormType type) const
{
  if (result.size() != numCols)
    result.resize(numCols);

  switch (type) {

  case NOX::Abstract::Vector::MaxNorm:
    double value;
    for (int j=0; j<numCols; j++) {
      value = fabs(vecPtrs[j][0]);
      for (int i=1; i<len; i++)
	if (value < fabs(vecPtrs[j][i]))
	  value = fabs(vecPtrs[j][i]);
      result[j] = value;
    }
    break;

  case NOX::Abstract::Vector::OneNorm:
    for (int j=0; j<numCols; j++) 
      result[j] = DASUM_F77(&len, vecPtrs[j], &NOX::LAPACK::i_one);
    break;

  case NOX::Abstract::Vector::TwoNorm:
  default:
    for (int j=0; j<numCols; j++) 
      result[j] = DNRM2_F77(&len, vecPtrs[j], &NOX::LAPACK::i_one);
   break;

  }
}

void NOX::LAPACK::MultiVector::dot(
			   double alpha, 
			   const NOX::Abstract::MultiVector& y,
			   NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  dot(alpha, dynamic_cast<const NOX::LAPACK::MultiVector&>(y), b);
}

void NOX::LAPACK::MultiVector::dot(
			   double alpha, 
			   const NOX::LAPACK::MultiVector& y,
			   NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  checkSize(y.numCols);

  // Create strided versions of y and c (this) if they are not already strided
  // Note:  This is not the best test.  A multivector could not own its columns
  // and still be strided, i.e., if the indices in the view were contiguous
  const NOX::LAPACK::MultiVector *sy, *sc;
  if (y.ownsCols)
    sy = &y;
    //sy = const_cast<NOX::LAPACK::MultiVector *>(&y);
  else
    sy = new NOX::LAPACK::MultiVector(y);

  if (ownsCols)
    sc = this;
  else
    sc = new NOX::LAPACK::MultiVector(*this);
 
  // Get underlying arrays
  const double *yv = &(sy->x[0]);
  double *bv = b.values();
  const double *cv = &(sc->x[0]);

  // Call BLAS routine
  DGEMM_F77("T", "N", &len, &numCols, &numCols, &alpha, cv, &len, yv, &len,
	    &NOX::LAPACK::d_zero, bv, &numCols);

  // Delete temporaries
  if (!ownsCols) {
    delete sc;
  }
  if (!y.ownsCols)
    delete sy;
}

double& 
NOX::LAPACK::MultiVector::operator() (int i, int j)
{
  checkIndex(i,j);
  return vecPtrs[j][i];
}

const double& 
NOX::LAPACK::MultiVector::operator () (int i, int j) const
{
  checkIndex(i,j);
  return vecPtrs[j][i];
}

int 
NOX::LAPACK::MultiVector::length() const
{
  return len;
}

int 
NOX::LAPACK::MultiVector::numVectors() const
{
  return numCols;
}

void 
NOX::LAPACK::MultiVector::print() const
{
  for (int j=0; j<numCols; j++) {
    cout << "[ ";
    for (int i=0; i<len; i++)
      cout << vecPtrs[j][i] << " ";
    cout << "]" << endl;
  }
}

void 
NOX::LAPACK::MultiVector::checkIndex(int j) const 
{
  if ( j < 0 || j >= numCols ) {
    cerr << "NOX::LAPACK::MultiVector:  Error!  Invalid column index " 
	 << j << endl;
    throw "NOX::LAPACK Error";
  }
}

void 
NOX::LAPACK::MultiVector::checkIndex(int i, int j) const 
{
  if ( i < 0 || i >= len || j < 0 || j >= numCols ) {
    cerr << "NOX::LAPACK::MultiVector:  Error!  Invalid index (" 
	 << i << ", " << j << ")" << endl;
    throw "NOX::LAPACK Error";
  }
}

void 
NOX::LAPACK::MultiVector::checkDimensions() const
{
  if (len <= 0 || numCols <= 0) {
    cerr << "NOX::LAPACK::MultiVector:  Error!  Invalid dimensions "
	 << len << ", " << numCols << endl;
    throw "NOX::LAPACK Error";
  }
}

void
NOX::LAPACK::MultiVector::checkSize(int sz) const
{
  if (numCols != sz) {
    cerr << "NOX::LAPACK::MultiVector:  Error!  Size of supplied multivector "
	 << "is incompatible with this multivector" << endl;
    throw "NOX::LAPACK Error";
  }
}
