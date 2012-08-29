/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <fei_ParameterSet.hpp>
#include "fei_Matrix_Local.hpp"
#include "fei_Matrix_core.hpp"
#include "fei_sstream.hpp"
#include "fei_fstream.hpp"

namespace fei {

Matrix_Local::Matrix_Local(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
               fei::SharedPtr<fei::SparseRowGraph> sparseRowGraph)
 : matrixGraph_(matrixGraph),
   sparseRowGraph_(sparseRowGraph),
   coefs_(sparseRowGraph->packedColumnIndices.size(), 0.0),
   stateChanged_(false),
   work_data1D_(),
   work_data2D_()
{
}

Matrix_Local::~Matrix_Local()
{
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Matrix>
Matrix_Local::create_Matrix_Local(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                                bool blockEntry)
{
  fei::SharedPtr<fei::SparseRowGraph> srg =
    matrixGraph->createGraph(blockEntry, true);
  fei::SharedPtr<fei::Matrix> mat(new fei::Matrix_Local(matrixGraph, srg));
  return(mat);
}

//----------------------------------------------------------------------------
const char*
Matrix_Local::typeName()
{ return( "fei::Matrix_Local" ); }

//----------------------------------------------------------------------------
int
Matrix_Local::parameters(const fei::ParameterSet& /*paramset*/)
{
  return(0);
}

//----------------------------------------------------------------------------
int
Matrix_Local::parameters(int /*numParams*/, const char* const* /*paramStrings*/)
{
  return(0);
}

fei::SharedPtr<fei::MatrixGraph>
Matrix_Local::getMatrixGraph() const
{ return( matrixGraph_ ); }

void
Matrix_Local::setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
{ matrixGraph_ = matrixGraph; }

int
Matrix_Local::getGlobalNumRows() const
{ return( sparseRowGraph_->rowNumbers.size() ); }

int
Matrix_Local::getLocalNumRows() const
{ return( getGlobalNumRows() ); }

int
Matrix_Local::getRowIndex(int rowNumber) const
{
  int* rows = &(sparseRowGraph_->rowNumbers[0]);
  int numRows = getLocalNumRows();
  return( fei::binarySearch(rowNumber, rows, numRows) );
}

int
Matrix_Local::getRowLength(int row, int& length) const
{
  int idx = getRowIndex(row);
  if (idx < 0) return(idx);

  length = sparseRowGraph_->rowOffsets[idx+1] -
           sparseRowGraph_->rowOffsets[idx];
  return(0);
}

int
Matrix_Local::putScalar(double scalar)
{
  for(unsigned i=0; i<coefs_.size(); ++i) coefs_[i] = scalar;
  stateChanged_ = true;
  return(0);
}

int
Matrix_Local::copyOutRow(int row, int len, double* coefs, int* indices) const
{
  int idx = getRowIndex(row);
  if (idx < 0) return(idx);

  int offset = sparseRowGraph_->rowOffsets[idx];
  int length = sparseRowGraph_->rowOffsets[idx+1]-offset;
  if (length > len) length = len;

  for(int i=0; i<length; ++i) {
    indices[i] = sparseRowGraph_->packedColumnIndices[offset+i];
    coefs[i] = coefs_[offset+i];
  }

  return(0);
}

int
Matrix_Local::giveToMatrix(int numRows, const int* rows,
                      int numCols, const int* cols,
                      const double* const* values,
                      bool sumInto,
                      int format)
{
  if (numRows == 0 || numCols == 0) {
    return(0);
  }

  if (format != FEI_DENSE_ROW && format != FEI_DENSE_COL) {
    return(-1);
  }

  const double** myvalues = const_cast<const double**>(values);
  if (format != FEI_DENSE_ROW) {
    fei::Matrix_core::copyTransposeToWorkArrays(numRows, numCols, values,
                              work_data1D_, work_data2D_);
    myvalues = &work_data2D_[0];
  }

  for(int i=0; i<numRows; ++i) {
    int idx = getRowIndex(rows[i]);
    if (idx < 0) {
      throw std::runtime_error("fei::Matrix_Local::sumIn ERROR, row not found.");
    }

    int offset = sparseRowGraph_->rowOffsets[idx];
    int len = sparseRowGraph_->rowOffsets[idx+1] - offset;

    int* colInds = &(sparseRowGraph_->packedColumnIndices[offset]);
    double* coefs   = &(coefs_[offset]);

    for(int j=0; j<numCols; ++j) {
      int idx2 = fei::binarySearch(cols[j], colInds, len);
      if (idx2 < 0) {
        throw std::runtime_error("fei::Matrix_Local::sumIn ERROR, col not found.");
      }

      if (sumInto) {
        coefs[idx2] += myvalues[i][j];
      }
      else {
        coefs[idx2] = myvalues[i][j];
      }
    }
  }

  stateChanged_ = true;
  return(0);
}

int
Matrix_Local::sumIn(int numRows, const int* rows,
                    int numCols, const int* cols,
                    const double* const* values,
                    int format)
{
  return( giveToMatrix(numRows, rows, numCols, cols, values,
                       true, format) );
}

int
Matrix_Local::copyIn(int numRows, const int* rows,
                       int numCols, const int* cols,
                       const double* const* values,
                      int format)
{
  return( giveToMatrix(numRows, rows, numCols, cols, values,
                       false, format) );
}

int
Matrix_Local::sumInFieldData(int fieldID,
                               int idType,
                               int rowID,
                               int colID,
                               const double* const* data,
                               int format)
{
  fei::SharedPtr<fei::VectorSpace> rspace = matrixGraph_->getRowSpace();
  fei::SharedPtr<fei::VectorSpace> cspace = matrixGraph_->getColSpace();

  int fieldSize = (int)rspace->getFieldSize(fieldID);
  std::vector<int> indices(2*fieldSize);

  rspace->getGlobalIndex(idType, rowID, fieldID, indices[0]);
  for(int i=1; i<fieldSize; ++i) {
    indices[i] = indices[0]+i;
  }

  cspace->getGlobalIndex(idType, colID, fieldID, indices[fieldSize]);
  for(int i=1; i<fieldSize; ++i) {
    indices[fieldSize+i] = indices[fieldSize]+i;
  }

  return( giveToMatrix(fieldSize, &indices[0], fieldSize, &indices[fieldSize],
                       data, true, format) );
}

int
Matrix_Local::sumInFieldData(int fieldID,
                               int idType,
                               int rowID,
                               int colID,
                               const double* data,
                               int format)
{
  fei::SharedPtr<fei::VectorSpace> rspace = matrixGraph_->getRowSpace();

  int fieldSize = (int)rspace->getFieldSize(fieldID);
  std::vector<const double*> data2D(fieldSize);

  int offset = 0;
  for(int i=0; i<fieldSize; ++i) {
    data2D[i] = &data[offset];
    offset += fieldSize;
  }

  return( sumInFieldData(fieldID, idType, rowID, colID,
                         &data2D[0], format) );
}

int
Matrix_Local::sumIn(int blockID, int connectivityID,
                    const double* const* values,
                    int format)
{
  int numIndices = matrixGraph_->getConnectivityNumIndices(blockID);
  std::vector<int> indices(numIndices);

  matrixGraph_->getConnectivityIndices(blockID, connectivityID,
                                       numIndices, &indices[0], numIndices);

  return( giveToMatrix(numIndices, &indices[0], numIndices, &indices[0],
                       values, true, format) );
}

int
Matrix_Local::globalAssemble()
{ return(0); }

int
Matrix_Local::multiply(fei::Vector* x,
                       fei::Vector* y)
{
  FEI_COUT << "fei::Matrix_Local::multiply NOT IMPLEMENTED."<<FEI_ENDL;
  return(-1);
}

void
Matrix_Local::setCommSizes()
{
}

int
Matrix_Local::gatherFromOverlap(bool accumulate)
{
  (void)accumulate;
  return(0);
}

int
Matrix_Local::writeToFile(const char* filename,
                          bool matrixMarketFormat)
{
  fei::SharedPtr<fei::VectorSpace> vspace = matrixGraph_->getRowSpace();

  MPI_Comm comm = vspace->getCommunicator();

  int localProc = fei::localProc(comm);

  FEI_OSTRINGSTREAM osstr;
  osstr << filename << "." << localProc << ".mtx";
  std::string fullname = osstr.str();

  FEI_OFSTREAM ofstr(fullname.c_str(), IOS_OUT);

  return( writeToStream(ofstr, matrixMarketFormat) );
}

int
Matrix_Local::writeToStream(FEI_OSTREAM& ostrm,
                            bool matrixMarketFormat)
{
  static char mmbanner[] = "%%MatrixMarket matrix coordinate real general";

  fei::SharedPtr<fei::VectorSpace> rspace = matrixGraph_->getRowSpace();
  fei::SharedPtr<fei::VectorSpace> cspace = matrixGraph_->getColSpace();

  int numRows = getLocalNumRows();
  int numCols = cspace->getEqnNumbers().size();
  int nnz = coefs_.size();

  if (matrixMarketFormat) {
    ostrm << mmbanner << FEI_ENDL;
    ostrm << numRows << " " << numCols << " " << nnz << FEI_ENDL;
  }
  else {
    ostrm << numRows << " " << numCols << " "<< FEI_ENDL;
  }

  std::vector<int>& rowNumbers = sparseRowGraph_->rowNumbers;
  std::vector<int>& rowOffsets = sparseRowGraph_->rowOffsets;
  std::vector<int>& colIndices = sparseRowGraph_->packedColumnIndices;

  ostrm.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
  ostrm.precision(13);

  int offset = 0;
  for(unsigned i=0; i<rowNumbers.size(); ++i) {
    int rowlen = rowOffsets[i+1]-rowOffsets[i];

    for(int j=0; j<rowlen; ++j) {
      if (matrixMarketFormat) {
        ostrm << rowNumbers[i]+1 << " " << colIndices[offset]+1
           << " " << coefs_[offset] << FEI_ENDL;
      }
      else {
        ostrm << rowNumbers[i] << " " << colIndices[offset]
           << " " << coefs_[offset] << FEI_ENDL;
      }
      ++offset;
    }
  }

  return(0);
}

bool
Matrix_Local::usingBlockEntryStorage()
{ return( false ); }

void
Matrix_Local::markState()
{
  stateChanged_ = false;
}

bool
Matrix_Local::changedSinceMark()
{ return(stateChanged_); }

const std::vector<int>&
Matrix_Local::getRowNumbers() const
{ return( sparseRowGraph_->rowNumbers ); }

const std::vector<int>&
Matrix_Local::getRowOffsets() const
{ return( sparseRowGraph_->rowOffsets ); }

const std::vector<int>&
Matrix_Local::getColumnIndices() const
{ return( sparseRowGraph_->packedColumnIndices ); }

const std::vector<double>&
Matrix_Local::getCoefs() const
{ return( coefs_ ); }

}//namespace fei

