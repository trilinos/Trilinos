// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Extended_MultiVector.H"
#include "LOCA_Extended_Vector.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Extended::MultiVector::MultiVector(
                   const LOCA::Extended::MultiVector& source,
                   NOX::CopyType type) :
  globalData(source.globalData),
  numColumns(source.numColumns),
  numMultiVecRows(source.numMultiVecRows),
  numScalarRows(source.numScalarRows),
  multiVectorPtrs(numMultiVecRows),
  scalarsPtr(),
  extendedVectorPtrs(numColumns),
  isView(false)
{
  // Copy multi vecs
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i] = source.multiVectorPtrs[i]->clone(type);

  // Copy scalars

  scalarsPtr =
    Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::Copy,
                                                             *source.scalarsPtr));

  for (int i=0; i<numColumns; i++) {
    extendedVectorPtrs[i] = Teuchos::null;
  }
}

LOCA::Extended::MultiVector::MultiVector(
                   const LOCA::Extended::MultiVector& source,
                   int nColumns) :
  globalData(source.globalData),
  numColumns(nColumns),
  numMultiVecRows(source.numMultiVecRows),
  numScalarRows(source.numScalarRows),
  multiVectorPtrs(numMultiVecRows),
  scalarsPtr(),
  extendedVectorPtrs(numColumns),
  isView(false)
{
  // Clone multivec blocks
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i] = source.multiVectorPtrs[i]->clone(numColumns);

  for (int i=0; i<numColumns; i++) {
    extendedVectorPtrs[i] = Teuchos::null;
  }

  scalarsPtr =
    Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(numScalarRows,
                                 numColumns));
}

LOCA::Extended::MultiVector::MultiVector(
                   const LOCA::Extended::MultiVector& source,
                   const std::vector<int>& index, bool view) :
  globalData(source.globalData),
  numColumns(index.size()),
  numMultiVecRows(source.numMultiVecRows),
  numScalarRows(source.numScalarRows),
  multiVectorPtrs(numMultiVecRows),
  scalarsPtr(),
  extendedVectorPtrs(numColumns),
  isView(view)
{
  // Check indices are valid
  for (unsigned int j=0; j<index.size(); j++)
    source.checkIndex("LOCA::Extended::MultiVector()", index[j]);

  for (int i=0; i<numColumns; i++) {
    extendedVectorPtrs[i] = Teuchos::null;
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
    Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View,
                                 vals,
                                 numScalarRows,
                                 numScalarRows,
                                 numColumns));
    }
    else {
      globalData->locaErrorCheck->throwError(
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
    Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::Copy,
                                 vals,
                                 numScalarRows,
                                 numScalarRows,
                                 numColumns));
    }
    else {
      scalarsPtr =
    Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(numScalarRows,
                                 numColumns));
      for (int j=0; j<numColumns; j++)
    for (int i=0; i<numScalarRows; i++)
      (*scalarsPtr)(i,j) = (*source.scalarsPtr)(i,index[j]);
    }
  }
}

LOCA::Extended::MultiVector::~MultiVector()
{
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
  operator=(dynamic_cast<const LOCA::Extended::MultiVector&>(source));
  return *this;
}

LOCA::Extended::MultiVector&
LOCA::Extended::MultiVector::operator=(
                   const LOCA::Extended::MultiVector& source)
{
  if (this != &source) {

    // Verify dimensions are consistent
    checkDimensions("LOCA::Extended::MultiVector::operator=()", source);

    globalData = source.globalData;

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
                      const std::vector<int>& index)
{
  return setBlock(dynamic_cast<const LOCA::Extended::MultiVector&>(source),
          index);
}

NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::setBlock(
                   const LOCA::Extended::MultiVector& source,
                   const std::vector<int>& index)
{
  // Verify dimensions are consistent
  if (source.numMultiVecRows != numMultiVecRows ||
      source.numScalarRows != numScalarRows)
    globalData->locaErrorCheck->throwError(
    "LOCA::Extended::MultiVector::setBlock()",
    "Size of supplied multivector is incompatible with this multivector");
  if (static_cast<unsigned int>(source.numColumns) != index.size()) {
   globalData->locaErrorCheck->throwError(
    "LOCA::Extended::MultiVector::setBlock()",
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
    globalData->locaErrorCheck->throwError(
           "LOCA::Extended::MultiVector::augment()",
           "Augmenting a multivector view is not supported");
  }

  // Verify dimensions are consistent
  if (source.numMultiVecRows != numMultiVecRows ||
      source.numScalarRows != numScalarRows)
    globalData->locaErrorCheck->throwError(
    "LOCA::Extended::MultiVector::augment()",
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
    extendedVectorPtrs[numColumns + i] = Teuchos::null;

  // Update number of columns
  numColumns += source.numColumns;

  return *this;
}

NOX::Abstract::Vector&
LOCA::Extended::MultiVector::operator [] (int i)
{
  return *(getVector(i));
}

const NOX::Abstract::Vector&
LOCA::Extended::MultiVector::operator [] (int i) const
{
  return *(getVector(i));
}

NOX::Abstract::MultiVector&
LOCA::Extended::MultiVector::scale(double gamma)
{
  // scale each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->scale(gamma);

  // scale scalars
  scalarsPtr->scale(gamma);

  return *this;
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
    globalData->locaErrorCheck->throwError(
      "LOCA::Extended::MultiVector::update()",
      "Size of supplied multivector is incompatible with this multivector");

  if (transb == Teuchos::NO_TRANS) {
    if (a.numColumns != b.numRows() || numColumns != b.numCols())
      globalData->locaErrorCheck->throwError(
        "LOCA::Extended::MultiVector::update()",
        "Size of supplied matrix is incompatible with this multivector");
  }
  else
    if (a.numColumns != b.numCols() || numColumns != b.numRows())
      globalData->locaErrorCheck->throwError(
        "LOCA::Extended::MultiVector::update()",
        "Size of supplied matrix is incompatible with this multivector");


  // update each multivec
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->update(transb, alpha, *(a.multiVectorPtrs[i]), b,
                   gamma);

  // update scalars
  if (numScalarRows > 0)
    scalarsPtr->multiply(Teuchos::NO_TRANS, transb, alpha, *(a.scalarsPtr),
             b, gamma);

  return *this;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Extended::MultiVector::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new LOCA::Extended::MultiVector(*this, type));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Extended::MultiVector::clone(int numvecs) const
{
  return Teuchos::rcp(new LOCA::Extended::MultiVector(*this, numvecs));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Extended::MultiVector::subCopy(const std::vector<int>& index) const
{
  return Teuchos::rcp(new LOCA::Extended::MultiVector(*this, index, false));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Extended::MultiVector::subView(const std::vector<int>& index) const
{
  return Teuchos::rcp(new LOCA::Extended::MultiVector(*this, index, true));
}

void
LOCA::Extended::MultiVector::norm(std::vector<double>& result,
               NOX::Abstract::Vector::NormType type) const
{

  // Make sure result vector is of appropriate size
  if (result.size() != static_cast<unsigned int>(numColumns))
    result.resize(numColumns);

  // Zero out result vector
  for (int i=0; i<numColumns; i++)
    result[i] = 0.0;

  // Intermediate vector to hold norms of a multivector
  std::vector<double> vecNorm(result);

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
    globalData->locaErrorCheck->throwError(
  "LOCA::Extended::MultiVector::multiply()",
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
  if (numScalarRows > 0)
    b.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, *y.scalarsPtr,
           *scalarsPtr, 1.0);

}

NOX::size_type
LOCA::Extended::MultiVector::length() const
{
  NOX::size_type len = 0;

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
LOCA::Extended::MultiVector::print(std::ostream& stream) const
{
  // Print multi-vecs
  for (int i=0; i<numMultiVecRows; i++)
    multiVectorPtrs[i]->print(stream);

  // Print scalars
  scalarsPtr->print(stream);
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::Extended::MultiVector::getMultiVector(int i) const
{
  checkVectorRowIndex("LOCA::Extended::MultiVector::getMultiVector()", i);

  return multiVectorPtrs[i];
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Extended::MultiVector::getMultiVector(int i)
{
  checkVectorRowIndex("LOCA::Extended::MultiVector::getMultiVector()", i);

  return multiVectorPtrs[i];
}

Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Extended::MultiVector::getScalars() const
{
  return scalarsPtr;
}

Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Extended::MultiVector::getScalars()
{
  return scalarsPtr;
}

Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Extended::MultiVector::getScalarRows(int num_rows, int row) const
{
  return
    Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View,
                                 *scalarsPtr,
                                 num_rows,
                                 numColumns,
                                 row,
                                 0));
}

Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Extended::MultiVector::getScalarRows(int num_rows, int row)
{
  return
    Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View,
                                 *scalarsPtr,
                                 num_rows,
                                 numColumns,
                                 row,
                                 0));
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

Teuchos::RCP<LOCA::Extended::Vector>
LOCA::Extended::MultiVector::getVector(int i)
{
  // Verify index is valid
  checkIndex("LOCA::Extended::MultiVector::vector()", i);

  // Create extended vector view if necessary
  if (extendedVectorPtrs[i] == Teuchos::null) {
    extendedVectorPtrs[i] = generateVector(numMultiVecRows, numScalarRows);
    for (int k=0; k<numMultiVecRows; k++)
      extendedVectorPtrs[i]->setVectorView(
                    k,
                    Teuchos::rcp(&(*multiVectorPtrs[k])[i],
                             false));
    if (numScalarRows > 0)
      extendedVectorPtrs[i]->setScalarArray((*scalarsPtr)[i]);
  }

  return extendedVectorPtrs[i];
}

Teuchos::RCP<const LOCA::Extended::Vector>
LOCA::Extended::MultiVector::getVector(int i) const
{
  // Verify index is valid
  checkIndex("LOCA::Extended::MultiVector::vector()", i);

  // Create extended vector view if necessary
  if (extendedVectorPtrs[i] == Teuchos::null) {
    extendedVectorPtrs[i] = generateVector(numMultiVecRows, numScalarRows);
    for (int k=0; k<numMultiVecRows; k++)
      extendedVectorPtrs[i]->setVectorView(
                      k,
                      Teuchos::rcp(&(*multiVectorPtrs[k])[i],
                           false));
    if (numScalarRows > 0)
      extendedVectorPtrs[i]->setScalarArray((*scalarsPtr)[i]);
  }

  return extendedVectorPtrs[i];
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

LOCA::Extended::MultiVector::MultiVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            int nColumns, int nVectorRows,
            int nScalarRows) :
  globalData(global_data),
  numColumns(nColumns),
  numMultiVecRows(nVectorRows),
  numScalarRows(nScalarRows),
  multiVectorPtrs(numMultiVecRows),
  scalarsPtr(),
  extendedVectorPtrs(numColumns),
  isView(false)
{
  for (int i=0; i<numColumns; i++) {
    extendedVectorPtrs[i] = Teuchos::null;
  }

  scalarsPtr =
    Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(numScalarRows,
                                 numColumns));
}

Teuchos::RCP<LOCA::Extended::Vector>
LOCA::Extended::MultiVector::generateVector(int nVecs, int nScalarRows) const
{
  return Teuchos::rcp(new LOCA::Extended::Vector(globalData, nVecs,
                         nScalarRows));
}

void
LOCA::Extended::MultiVector::setMultiVectorPtr(
               int i,
               Teuchos::RCP<NOX::Abstract::MultiVector> v)
{
  checkVectorRowIndex("LOCA::Extended::MultiVector::setMultiVectorPtr()",i);

  multiVectorPtrs[i] = v;
}

void
LOCA::Extended::MultiVector::checkDimensions(
                  const std::string& callingFunction,
                  const LOCA::Extended::MultiVector& a) const
{
  if (a.numMultiVecRows != numMultiVecRows || a.numColumns != numColumns ||
      a.numScalarRows != numScalarRows)
    globalData->locaErrorCheck->throwError(callingFunction,
      "Size of supplied multivector is incompatible with this multivector");
}

void
LOCA::Extended::MultiVector::checkIndex(const std::string& callingFunction,
                    int i) const
{
  if ( i < 0 || i >= numColumns )
    globalData->locaErrorCheck->throwError(callingFunction,
                        "Invalid column index");
}

void
LOCA::Extended::MultiVector::checkVectorRowIndex(const std::string& callingFunction,
                         int i) const
{
  if ( i < 0 || i >= numMultiVecRows)
    globalData->locaErrorCheck->throwError(callingFunction,
                        "Invalid vector row index");
}
void
LOCA::Extended::MultiVector::checkIndex(const std::string& callingFunction,
                    int i, int j) const
{
  if ( i < 0 || i >= numScalarRows )
    globalData->locaErrorCheck->throwError(callingFunction,
                        "Invalid row index");
  if ( j < 0 || j >= numColumns )
    globalData->locaErrorCheck->throwError(callingFunction,
                        "Invalid column index");
}

bool
LOCA::Extended::MultiVector::isContiguous(const std::vector<int>& index) const
{
  for (unsigned int i=0; i<index.size(); i++) {
    if (static_cast<unsigned int>(index[i]) != index[0] + i)
      return false;
  }
  return true;
}
