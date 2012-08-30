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
#include "fei_MatrixReducer.hpp"
#include "fei_EqnComm.hpp"
#include "fei_Matrix_core.hpp"
#include "fei_sstream.hpp"
#include "fei_fstream.hpp"
#include <fei_CommUtils.hpp>

namespace fei {

MatrixReducer::MatrixReducer(fei::SharedPtr<fei::Reducer> reducer,
                             fei::SharedPtr<fei::Matrix> target)
 : reducer_(reducer),
   target_(target),
   globalAssembleCalled_(false),
   changedSinceMark_(false)
{
  fei::SharedPtr<fei::VectorSpace> vspace =
    target->getMatrixGraph()->getRowSpace();
  MPI_Comm comm = vspace->getCommunicator();
  int numLocal = reducer_->getLocalReducedEqns().size();
  fei::SharedPtr<fei::EqnComm> eqnComm(new fei::EqnComm(comm, numLocal));
  fei::Matrix_core* target_core =
    dynamic_cast<fei::Matrix_core*>(target_.get());
    if (target_core == NULL) {
    throw std::runtime_error("fei::MatrixReducer ERROR, target matrix not dynamic_cast-able to fei::Matrix_core.");
  }

  target_core->setEqnComm(eqnComm);
}

MatrixReducer::~MatrixReducer()
{
}

int
MatrixReducer::parameters(const fei::ParameterSet& paramset)
{
  return(target_->parameters(paramset));
}

void
MatrixReducer::setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
{
  target_->setMatrixGraph(matrixGraph);
}

int
MatrixReducer::getGlobalNumRows() const
{
  return(target_->getMatrixGraph()->getRowSpace()->getGlobalNumIndices());
}

int
MatrixReducer::getLocalNumRows() const
{
  return(target_->getMatrixGraph()->getRowSpace()->getNumIndices_Owned());
}

int MatrixReducer::putScalar(double scalar)
{ return(target_->putScalar(scalar)); }

int
MatrixReducer::getRowLength(int row, int& length) const
{
  if (reducer_->isSlaveEqn(row)) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixReducer::getRowLength ERROR, row="<<row<<" is a slave eqn. You can't get a slave row from the reduced matrix.";
    throw std::runtime_error(osstr.str());
  }

  int reducedrow = reducer_->translateToReducedEqn(row);
  return(target_->getRowLength(reducedrow, length));
}

int
MatrixReducer::copyOutRow(int row, int len, double* coefs, int* indices) const
{
  if (reducer_->isSlaveEqn(row)) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixReducer::copyOutRow ERROR, requested row ("<<row
      <<") is a slave eqn. You can't get a slave row from the reduced matrix.";
    throw std::runtime_error(osstr.str());
  }

  int reducedrow = reducer_->translateToReducedEqn(row);
  int err = target_->copyOutRow(reducedrow, len, coefs, indices);
  for(int i=0; i<len; ++i) {
    indices[i] = reducer_->translateFromReducedEqn(indices[i]);
  }
  return(err);
}

int
MatrixReducer::sumIn(int numRows, const int* rows,
                     int numCols, const int* cols,
                     const double* const* values,
                     int format)
{
  int err = reducer_->addMatrixValues(numRows, rows, numCols, cols,
                                      values, true, *target_, format);
  return(err);
}

int
MatrixReducer::copyIn(int numRows, const int* rows,
                      int numCols, const int* cols,
                      const double* const* values,
                      int format)
{
  int err = reducer_->addMatrixValues(numRows, rows, numCols, cols,
                                      values, false, *target_, format);
  return(err);
}

int
MatrixReducer::sumInFieldData(int fieldID,
                              int idType,
                              int rowID,
                              int colID,
                              const double* const* data,
                              int format)
{
  fei::SharedPtr<fei::VectorSpace> rowSpace =
    target_->getMatrixGraph()->getRowSpace();
  fei::SharedPtr<fei::VectorSpace> colSpace =
    target_->getMatrixGraph()->getColSpace();

  unsigned fieldSize = rowSpace->getFieldSize(fieldID);
  std::vector<int> indices(fieldSize*2);
  int* rowIndices = &indices[0];
  int* colIndices = rowIndices+fieldSize;

  rowSpace->getGlobalIndices(1, &rowID, idType, fieldID, rowIndices);
  colSpace->getGlobalIndices(1, &colID, idType, fieldID, colIndices);

  if (format != FEI_DENSE_ROW) {
    throw std::runtime_error("MatrixReducer: bad format");
  }

  int err = reducer_->addMatrixValues(fieldSize, rowIndices,
                                      fieldSize, colIndices,
                                      data, true, *target_, format);
  return(err);
}

int
MatrixReducer::sumInFieldData(int fieldID,
			       int idType,
			       int rowID,
			       int colID,
			       const double* data,
			       int format)
{
  fei::SharedPtr<fei::VectorSpace> rowSpace =
    target_->getMatrixGraph()->getRowSpace();

  unsigned fieldSize = rowSpace->getFieldSize(fieldID);

  std::vector<const double*> data_2d(fieldSize);
  for(unsigned i=0; i<fieldSize; ++i) {
    data_2d[i] = &data[i*fieldSize];
  }

  return(sumInFieldData(fieldID, idType, rowID, colID, &data_2d[0], format));
}

int
MatrixReducer::sumIn(int blockID, int connectivityID,
		      const double* const* values,
		      int format)
{
  fei::SharedPtr<fei::MatrixGraph> matGraph = getMatrixGraph();
  int numRowIndices, numColIndices, dummy;
  matGraph->getConnectivityNumIndices(blockID, numRowIndices, numColIndices);

  std::vector<int> indices(numRowIndices+numColIndices);
  int* rowIndices = &indices[0];
  int* colIndices = rowIndices+numRowIndices;

  matGraph->getConnectivityIndices(blockID, connectivityID,
                                   numRowIndices, rowIndices, dummy,
                                   numColIndices, colIndices, dummy);

  return(sumIn(numRowIndices, rowIndices, numColIndices, colIndices,
               values, format));
}

int
MatrixReducer::globalAssemble()
{
  reducer_->assembleReducedMatrix(*target_);
  return(target_->globalAssemble());
}

int MatrixReducer::multiply(fei::Vector* x, fei::Vector* y)
{ return(target_->multiply(x, y)); }

int MatrixReducer::gatherFromOverlap(bool accumulate)
{
  reducer_->assembleReducedMatrix(*target_);
  target_->setCommSizes();
  return(target_->gatherFromOverlap(accumulate));
}

int MatrixReducer::writeToFile(const char* filename,
			    bool matrixMarketFormat)
{
  static char mmbanner[] = "%%MatrixMarket matrix coordinate real general";
  std::vector<int>& localrows = reducer_->getLocalReducedEqns();
  int localNumRows = localrows.size();

  int globalNNZ = 0;
  int localNNZ = 0;

  for(int i=0; i<localNumRows; ++i) {
    int len;
    CHK_ERR( target_->getRowLength(localrows[i], len) );
    localNNZ += len;
  }

  MPI_Comm comm = getMatrixGraph()->getRowSpace()->getCommunicator();

  CHK_MPI( fei::GlobalSum(comm, localNNZ, globalNNZ) );
  int globalNumRows = 0;
  CHK_MPI( fei::GlobalSum(comm, localNumRows, globalNumRows) );

  int globalNumCols = globalNumRows;

  for(int p=0; p<fei::numProcs(comm); ++p) {
    fei::Barrier(comm);
    if (p != fei::localProc(comm)) continue;

    FEI_OFSTREAM* outFile = NULL;
    if (p==0) {
      outFile = new FEI_OFSTREAM(filename, IOS_OUT);
      FEI_OFSTREAM& ofs = *outFile;
      if (matrixMarketFormat) {
        ofs << mmbanner << FEI_ENDL;
        ofs <<globalNumRows<< " " <<globalNumCols<< " " <<globalNNZ<<FEI_ENDL;
      }
      else {
        ofs <<globalNumRows<< " " <<globalNumCols<<FEI_ENDL;
      }
    }
    else outFile = new FEI_OFSTREAM(filename, IOS_APP);

    outFile->setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
    outFile->precision(13);
    FEI_OFSTREAM& ofs = *outFile;

    int rowLength;
    std::vector<int> work_indices;
    std::vector<double> work_data1D;

    for(int i=0; i<localNumRows; ++i) {
      int row = localrows[i];
      CHK_ERR( target_->getRowLength(row, rowLength) );

      work_indices.resize(rowLength);
      work_data1D.resize(rowLength);

      int* indPtr = &work_indices[0];
      double* coefPtr = &work_data1D[0];

      CHK_ERR( target_->copyOutRow(row, rowLength, coefPtr, indPtr) );

      for(int j=0; j<rowLength; ++j) {
        if (matrixMarketFormat) {
          ofs << row+1 <<" "<<indPtr[j]+1<<" "<<coefPtr[j]<<FEI_ENDL;
        }
        else {
          ofs << row <<" "<<indPtr[j]<<" "<<coefPtr[j]<<FEI_ENDL;
        }
      }
    }

    delete outFile;
  }

  return(0);
}

int MatrixReducer::writeToStream(FEI_OSTREAM& ostrm,
			      bool matrixMarketFormat)
{
  return(target_->writeToStream(ostrm, matrixMarketFormat));
}

void MatrixReducer::markState()
{ target_->markState(); }

bool MatrixReducer::changedSinceMark()
{ return(target_->changedSinceMark()); }

}//namespace fei

