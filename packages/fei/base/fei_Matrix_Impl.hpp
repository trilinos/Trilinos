/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Matrix_Impl_hpp_
#define _fei_Matrix_Impl_hpp_

#include <fei_fwd.hpp>
#include <fei_mpi.h>
#include <fei_defs.h>
#include "fei_iostream.hpp"
#include "fei_fstream.hpp"
#include "fei_sstream.hpp"
#include <fei_ArrayUtils.hpp>
#include <fei_MatrixTraits.hpp>
#include <fei_MatrixTraits_LinProbMgr.hpp>
#include <fei_MatrixTraits_LinSysCore.hpp>
#include <fei_MatrixTraits_FEData.hpp>
#include <fei_MatrixTraits_FillableMat.hpp>
#include <fei_FillableVec.hpp>
#include <fei_FillableMat.hpp>

#include <snl_fei_FEMatrixTraits.hpp>
#include <snl_fei_FEMatrixTraits_FED.hpp>
#include <snl_fei_BlockMatrixTraits.hpp>
#include <fei_Matrix.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix_core.hpp>
#include <snl_fei_Utils.hpp>

#undef fei_file
#define fei_file "fei_Matrix_Impl.hpp"
#include <fei_ErrMacros.hpp>

namespace fei {

  /** To be used for local assembly of shared data. Provides operations for
      gathering the overlapped data (locally-stored shared data) to a non-
      overlapped data distribution (e.g., send shared data to owning processor)
      and vice-versa for scattering non-overlapped data to the overlapped
      distribution.

      When shared data that is not locally-owned is assembled into this
      object, it will be held locally until the above-mentioned gather
      operation is performed. When data that is locally-owned is assembled
      into this object, it will be passed directly to the underlying algebraic
      (non-overlapping) matrix.
  */
  template<typename T>
  class Matrix_Impl : public fei::Matrix, public fei::Matrix_core {
  public:
    /** Constructor */
    Matrix_Impl(fei::SharedPtr<T> matrix,
               fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                int numLocalEqns);

    /** Destructor */
    virtual ~Matrix_Impl();

    /** Return a name describing the run-time type
        of this object.
    */
    const char* typeName()
      {
        if (haveBlockMatrix()) {
          return(snl_fei::BlockMatrixTraits<T>::typeName());
        }
        else if (haveFEMatrix()) {
          return(snl_fei::FEMatrixTraits<T>::typeName());
        }
        else {
          return(fei::MatrixTraits<T>::typeName());
        }
      }

    /** Parameters method
     */
    int parameters(const fei::ParameterSet& paramset);

    /** Obtain the underlying matrix object. Note that this will generally be
     only the locally-owned portion of the matrix, not including any data
    corresponding to shared-but-not-owned nodes, etc. */
    fei::SharedPtr<T> getMatrix() { return( matrix_ ); }
    const fei::SharedPtr<T> getMatrix() const { return( matrix_ ); }

    fei::SharedPtr<fei::MatrixGraph> getMatrixGraph() const
      {return( Matrix_core::getMatrixGraph() ); }

    /** Set the fei::MatrixGraph associated with this matrix */
    void setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph);

    /** Get the global number of rows in the matrix.
     */
    int getGlobalNumRows() const
      {
        if ((int)(globalOffsets().size()) < numProcs()+1) return(-1);
        return(globalOffsets()[numProcs()]);
      }

    /** Get the local number of rows in the matrix.
     */
    int getLocalNumRows() const
      {
        return(lastLocalOffset() - firstLocalOffset() + 1);
      }

    /** Set a specified scalar throughout the matrix. */
    int putScalar(double scalar);

   /** Get the length of a row of the matrix.
       @param row Global 0-based equation number
       @param length Output. Length of the row.
       @return error-code non-zero if any error occurs.
   */
    int getRowLength(int row, int& length) const;

   /** Obtain a copy of the coefficients and indices for a row of the matrix.
       @param row Global 0-based equation number
       @param len Length of the caller-allocated coefs and indices arrays
       @param coefs Caller-allocated array, length 'len', to be filled with
       coefficients
       @param indices Caller-allocated array, length 'len', to be filled with
       indices. (These indices will be global 0-based equation numbers.)
       @return error-code non-zero if any error occurs.
   */
    int copyOutRow(int row, int len, double* coefs, int* indices) const;

    /** Sum coefficients into the matrix, adding them to any coefficients that
        may already exist at the specified row/column locations.

        @param numRows
        @param rows
        @param numCols
        @param cols
        @param values
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major. Others not recognized
     */
    int sumIn(int numRows, const int* rows,
              int numCols, const int* cols,
              const double* const* values,
              int format=0);

    /** Copy coefficients into the matrix, overwriting any coefficients that
        may already exist at the specified row/column locations.

        @param numRows
        @param rows
        @param numCols
        @param cols
        @param values
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major. Others not recognized
    */
    int copyIn(int numRows, const int* rows,
               int numCols, const int* cols,
               const double* const* values,
               int format=0);

    /** Sum coefficients into the matrix, specifying row/column locations by
        identifier/fieldID pairs.
        @param fieldID Input. field-identifier for which data is being input.
        @param idType Input. The identifier-type of the identifiers.
        @param rowID Input. Identifier in row-space, for which data is being
        input.
        @param colID Input. Identifier in column-space, for which data is being
        input.
        @param data Input. C-style table of data. num-rows is the field-size
        (i.e., number of scalar components that make up the field) of 'fieldID',
        as is num-columns.
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major. Others not recognized
        @return error-code 0 if successful
    */
    int sumInFieldData(int fieldID,
                       int idType,
                       int rowID,
                       int colID,
                       const double* const* data,
                       int format=0);

    /** Sum coefficients into the matrix, specifying row/column locations by
        identifier/fieldID pairs.
        @param fieldID Input. field-identifier for which data is being input.
        @param idType Input. The identifier-type of the identifiers.
        @param rowID Input. Identifier in row-space, for which data is being
        input.
        @param colID Input. Identifier in column-space, for which data is being
        input.
        @param data Input. 1-D list representing a packed table of data. Data may
        be backed in row-major or column-major order and this may be specified with
        the 'format' argument. The "table" of data is of size num-rows X num-columns
        and num-rows is the field-size (i.e., number of scalar components that
        make up the field) of 'fieldID', as is num-columns.
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major. Others not recognized
        @return error-code 0 if successful
    */
    int sumInFieldData(int fieldID,
                       int idType,
                       int rowID,
                       int colID,
                       const double* data,
                       int format=0);

    /** Sum coefficients, associated with a connectivity-block that was
        initialized on the MatrixGraph object, into this matrix.

        @param blockID
        @param connectivityID
        @param values
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major. Others not recognized
     */
    int sumIn(int blockID, int connectivityID,
              const double* const* values,
              int format=0);

    /** Perform any necessary internal communications/synchronizations or other
        operations appropriate at end of data input. For some implementations this
        will be a no-op.
    */
    int globalAssemble();

    /** Form a matrix-vector product y = 'this' * x
     */
    int multiply(fei::Vector* x,
                 fei::Vector* y);

    /** After local overlapping data has been input, (e.g., element-data for a
        finite-element application) call this method to have data that 
        corresponds to shared identifiers be communicated from sharing-but-not-
        owning processors, to owning processors.
    */
    int gatherFromOverlap(bool accumulate = true);

    /** Implementation of fei::Matrix::writeToFile */
    int writeToFile(const char* filename,
                    bool matrixMarketFormat=true);

    /** Implementation of fei::Matrix::writeToStream */

    int writeToStream(FEI_OSTREAM& ostrm,
                      bool matrixMarketFormat=true);

    bool usingBlockEntryStorage()
      {
        return(haveBlockMatrix());
      }

    /** for experts only */
    int giveToUnderlyingMatrix(int numRows, const int* rows,
                               int numCols, const int* cols,
                               const double* const* values,
                               bool sumInto,
                               int format);

    /** for experts only */
    int giveToUnderlyingBlockMatrix(int row,
                                    int rowDim,
                                    int numCols,
                                    const int* cols,
                                    const int* LDAs,
                                    const int* colDims,
                                    const double* const* values,
                                    bool sumInto);

    void markState();

    bool changedSinceMark();

  private:
    int giveToMatrix(int numRows, const int* rows,
                     int numCols, const int* cols,
                     const double* const* values,
                     bool sumInto,
                     int format);
 
    int giveToBlockMatrix(int numRows, const int* rows,
                          int numCols, const int* cols,
                          const double* const* values,
                          bool sumInto);

    fei::SharedPtr<T> matrix_;
    bool globalAssembleCalled_;
    bool changedSinceMark_;
    std::string dbgprefix_;
  };//class Matrix_Impl
}//namespace fei

//----------------------------------------------------------------------------
template<typename T>
inline int fei::Matrix_Impl<T>::globalAssemble()
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"globalAssemble"<<FEI_ENDL;
  }

  globalAssembleCalled_ = true;

  if (haveBlockMatrix()) {
    return( snl_fei::BlockMatrixTraits<T>::globalAssemble(matrix_.get()) );
  }
  return( fei::MatrixTraits<T>::globalAssemble(matrix_.get()) );
}

//----------------------------------------------------------------------------
template<typename T>
inline void fei::Matrix_Impl<T>::setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
{
  Matrix_core::setMatrixGraph(matrixGraph);
}

//----------------------------------------------------------------------------
template<typename T>
inline void fei::Matrix_Impl<T>::markState()
{
  changedSinceMark_ = false;
}

//----------------------------------------------------------------------------
template<typename T>
inline bool fei::Matrix_Impl<T>::changedSinceMark()
{
  return(changedSinceMark_);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::giveToUnderlyingMatrix(int numRows, const int* rows,
                                               int numCols, const int* cols,
                                               const double* const* values,
                                               bool sumInto,
                                               int format)
{
  if (format != FEI_DENSE_ROW) {
    ERReturn(-1);
  }

  if (output_level_ > fei::BRIEF_LOGS && output_stream_ != 0) {
    FEI_OSTREAM& os = *output_stream_;
    for(int i=0; i<numRows; ++i) {
      os << dbgprefix_<<"giveToUnderlyingMatrix ";
      for(int j=0; j<numCols; ++j) {
        os << "("<<rows[i]<<","<<cols[j]<<","<<values[i][j]<<") ";
      }
      os << FEI_ENDL;
    }
  }

  int err = fei::MatrixTraits<T>::putValuesIn(matrix_.get(), numRows, rows,
                                        numCols, cols, values, sumInto);
  if (err != 0) {
    return(err);
  }

  changedSinceMark_ = true;
  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::giveToUnderlyingBlockMatrix(int row,
                                                    int rowDim,
                                                    int numCols,
                                                    const int* cols,
                                                    const int* LDAs,
                                                    const int* colDims,
                                                    const double* const* values,
                                                    bool sumInto)
{
  if (sumInto) {
    if ( snl_fei::BlockMatrixTraits<T>::sumIn(matrix_.get(),
                                         row, rowDim, numCols, cols,
                                         LDAs, colDims, values) != 0) {
      ERReturn(-1);
    }
  }
  else {
    if ( snl_fei::BlockMatrixTraits<T>::copyIn(matrix_.get(),
                                          row, rowDim, numCols, cols,
                                          LDAs, colDims, values) != 0) {
      ERReturn(-1);
    }
  }

  changedSinceMark_ = true;

  return(0);
}

#include <fei_macros.hpp>
#include <fei_chk_mpi.hpp>

#include <fei_TemplateUtils.hpp>

#include <fei_Pattern.hpp>
#include <fei_VectorSpace.hpp>

#include <fei_ConnectivityBlock.hpp>

#include <snl_fei_PointBlockMap.hpp>

#include <fei_Record.hpp>

#include <snl_fei_Utils.hpp>

//----------------------------------------------------------------------------
template<typename T>
fei::Matrix_Impl<T>::Matrix_Impl(fei::SharedPtr<T> matrix,
                           fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                                int numLocalEqns)
  : Matrix_core(matrixGraph, numLocalEqns),
    matrix_(matrix),
    globalAssembleCalled_(false),
    changedSinceMark_(true),
    dbgprefix_("MatImpl: ")
{
  if (strcmp(snl_fei::FEMatrixTraits<T>::typeName(), "unsupported")) {
    setFEMatrix(true);
  }
  else {
    setFEMatrix(false);
  }

  if (strcmp(snl_fei::BlockMatrixTraits<T>::typeName(), "unsupported")) {
    setBlockMatrix(true);
  }
  else {
    setBlockMatrix(false);
  }
}

//----------------------------------------------------------------------------
template<typename T>
fei::Matrix_Impl<T>::~Matrix_Impl()
{
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::parameters(const fei::ParameterSet& paramset)
{
  Matrix_core::parameters(paramset);
  return 0;
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::getRowLength(int row, int& length) const
{
  if (haveBlockMatrix()) {
    return( snl_fei::BlockMatrixTraits<T>::getPointRowLength(matrix_.get(), row, length) );
  }
  else {
    int code = fei::MatrixTraits<T>::getRowLength(matrix_.get(), row, length);
    if (code < 0) {
      //matrix row probably not local. maybe it's shared.
      int proc = getOwnerProc(row);
      const FillableMat* remote_mat = getRemotelyOwnedMatrix(proc);
      if (remote_mat->hasRow(row)) {
        const FillableVec* row_entries = remote_mat->getRow(row);
        length = row_entries->size();
      }
      else {
        length = 0;
      }
    }
  }
  return 0;
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::putScalar(double scalar)
{
  if (output_level_ > fei::BRIEF_LOGS && output_stream_ != 0) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"putScalar("<<scalar<<")"<<FEI_ENDL;
  }

  if (haveFEMatrix()) {
    if (scalar != 0.0) return(-1);
    CHK_ERR( snl_fei::FEMatrixTraits<T>::reset(matrix_.get()) );
  }
  else if (haveBlockMatrix()) {
    if (globalAssembleCalled_ == true) {
      CHK_ERR( snl_fei::BlockMatrixTraits<T>::putScalar(matrix_.get(), scalar) );
    }
  }
  else {
    CHK_ERR( fei::MatrixTraits<T>::setValues(matrix_.get(), scalar) );
  }

  putScalar_remotelyOwned(scalar);

  changedSinceMark_ = true;

  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::copyOutRow(int row, int len,
                                   double* coefs, int* indices) const
{
  if (len < 1) {
    return 0;
  }

  if (haveBlockMatrix()) {
    int dummy;
    return( snl_fei::BlockMatrixTraits<T>::copyOutPointRow(matrix_.get(), firstLocalOffset(),
                                                  row, len,
                                                  coefs, indices, dummy));
  }
  else {
    int code = fei::MatrixTraits<T>::copyOutRow(matrix_.get(), row, len,
                                        coefs, indices);
    if (code < 0) {
      int proc = getOwnerProc(row);
      const FillableMat* remote_mat = getRemotelyOwnedMatrix(proc);
      if (remote_mat->hasRow(row)) {
        const FillableVec* row_entries = remote_mat->getRow(row);
        FillableVec::const_iterator it=row_entries->begin(); 
        for(int i=0; it!=row_entries->end(); ++it, ++i) {
          indices[i] = it->first;
          coefs[i] = it->second;
        }
      }
    }
  }
  return 0;
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::sumIn(int numRows, const int* rows,
                              int numCols, const int* cols,
                              const double* const* values,
                              int format)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"sumIn"<<FEI_ENDL;
    if (output_level_ >= fei::FULL_LOGS) {
      for(int i=0; i<numRows; ++i) {
        for(int j=0; j<numCols; ++j) {
          os << "("<<rows[i]<<","<<cols[j]<<","<<values[i][j]<<") ";
        }
        os << FEI_ENDL;
      }
    }
  }

  return( giveToMatrix( numRows, rows, numCols, cols, values, true, format) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::copyIn(int numRows, const int* rows,
                               int numCols, const int* cols,
                               const double* const* values,
                               int format)
{
  if (output_level_ > fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << "copyIn"<<FEI_ENDL;
    if (output_level_ >= fei::FULL_LOGS) {
      for(int i=0; i<numRows; ++i) {
        for(int j=0; j<numCols; ++j) {
          os << "("<<rows[i]<<","<<cols[j]<<","<<values[i][j]<<") ";
        }
        os << FEI_ENDL;
      }
    }
  }
  return( giveToMatrix( numRows, rows, numCols, cols, values, false, format) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::sumInFieldData(int fieldID,
                                       int idType,
                                       int rowID,
                                       int colID,
                                       const double* const* data,
                                       int format)
{
  fei::SharedPtr<fei::VectorSpace> vspace = vecSpace();

  int fieldSize = vspace->getFieldSize(fieldID);

  if (fieldSize <= 0) ERReturn(-1);

  work_indices_.resize(fieldSize*2);
  int* indicesPtr = &work_indices_[0];
  int i;

  CHK_ERR( vspace->getGlobalIndices(1, &rowID, idType, fieldID, indicesPtr));
  for(i=1; i<fieldSize; ++i) {
    indicesPtr[i] = indicesPtr[i-1]+1;
  }

  CHK_ERR( vspace->getGlobalIndices(1, &colID, idType, fieldID,
                                        &(indicesPtr[fieldSize])) );
  for(i=fieldSize+1; i<2*fieldSize; ++i) {
    indicesPtr[i] = indicesPtr[i-1]+1;
  }

  CHK_ERR( sumIn(fieldSize, indicesPtr, fieldSize, &(indicesPtr[fieldSize]),
                 data, format) );

  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::sumInFieldData(int fieldID,
                                       int idType,
                                       int rowID,
                                       int colID,
                                       const double* data,
                                       int format)
{
  fei::SharedPtr<fei::VectorSpace> vspace = vecSpace();

  int fieldSize = vspace->getFieldSize(fieldID);

  if (fieldSize <= 0) ERReturn(-1);

  work_data2D_.resize(fieldSize);

  const double** data2DPtr = &work_data2D_[0];
  for(int i=0; i<fieldSize; ++i) {
    data2DPtr[i] = &(data[i*fieldSize]);
  }

  CHK_ERR( sumInFieldData(fieldID, idType, rowID, colID, data2DPtr, format) );

  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::sumIn(int blockID, int connectivityID,
                              const double* const* values, int format)
{
  if (output_level_ > fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"sumIn blkID=" << blockID
       << ", connID=" << connectivityID << FEI_ENDL;
  }

  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();

  if (mgraph.get() == NULL) ERReturn(-1);

  const fei::ConnectivityBlock* cblock = mgraph->getConnectivityBlock(blockID);
  if (cblock==NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::Matrix_Impl::sumIn ERROR, unable to "
          << "look up connectivity-block with ID "<<blockID;
    throw std::runtime_error(osstr.str());
  }

  bool symmetric = cblock->isSymmetric();
  const fei::Pattern* pattern = cblock->getRowPattern();
  const fei::Pattern* colpattern = symmetric ? NULL : cblock->getColPattern();
  const fei::Record<int>*const* rowConn = cblock->getRowConnectivity(connectivityID);

  fei::SharedPtr<fei::VectorSpace> rspace = mgraph->getRowSpace();
  rspace->getGlobalIndices(pattern, rowConn, work_indices2_);

  int numRowIndices = work_indices2_.size();
  int* rowIndices = &work_indices2_[0];

  if (haveFEMatrix() || haveBlockMatrix()) {
    FieldDofMap<int>& fdofmap = rspace->getFieldDofMap();
    const std::map<int,int>& connIDs = cblock->getConnectivityIDs();
    std::map<int,int>::const_iterator
      iter = connIDs.find(connectivityID);
    if (iter == connIDs.end()) ERReturn(-1);
    int connOffset = iter->second;
    int numIDs = pattern->getNumIDs();
    const int* indicesPerID = pattern->getNumIndicesPerID();
    const int* fieldsPerID = pattern->getNumFieldsPerID();
    const int* fieldIDs = pattern->getFieldIDs();

    int numDofs = 0;
    for(int ii=0; ii<numIDs; ++ii) {
      numDofs += indicesPerID[ii];
    }

    const int* numIndicesPerID = pattern->getNumIndicesPerID();
    work_indices_.resize(numIDs+numDofs);
    int i, *nodeNumbers = &work_indices_[0];
    int* dof_ids = nodeNumbers+numIDs;

    int foffset = 0;
    int doffset = 0;
    for(i=0; i<numIDs; ++i) {
      nodeNumbers[i] = rowConn[i]->getNumber();
      for(int ii=0; ii<fieldsPerID[i]; ++ii) {
        int fieldSize = rspace->getFieldSize(fieldIDs[foffset]);
        int dof_id = fdofmap.get_dof_id(fieldIDs[foffset++], 0);
        for(int j=0; j<fieldSize; ++j) {
          dof_ids[doffset++] = dof_id + j;
        }
      }
    }

    if (haveFEMatrix()) {
      CHK_ERR( snl_fei::FEMatrixTraits<T>::sumInElemMatrix(matrix_.get(), blockID, connOffset,
                                                  numIDs, nodeNumbers,
                                                  numIndicesPerID, dof_ids, values) );
      changedSinceMark_ = true;
    }

    if (haveBlockMatrix()) {
      if (format != FEI_DENSE_ROW && format != FEI_DENSE_COL &&
          format != FEI_BLOCK_DIAGONAL_ROW) {
        fei::console_out() << "fei::Matrix_Impl::sumIn ERROR, for block-matrix, format must"
                 << " be FEI_DENSE_ROW or FEI_DENSE_COL."<<FEI_ENDL;
        ERReturn(-1);
      }

      int numPtIndices = pattern->getNumIndices();
      int* ptIndices = &work_indices2_[0];

      int numPtColIndices = symmetric ? numPtIndices : colpattern->getNumIndices();

      int len = numPtIndices*numPtColIndices;

      if (format == FEI_BLOCK_DIAGONAL_ROW) {
        len = 0;
        for(i=0; i<numIDs; ++i) {
          len += numIndicesPerID[i]*numIndicesPerID[i];
        }
      }

      double* ccvalues = new double[len];
      //Ouch! Data copy! Remember to optimize this later...
      if (format == FEI_BLOCK_DIAGONAL_ROW) {
        snl_fei::copy2DBlockDiagToColumnContig(numIDs, numIndicesPerID,
                                               values, format, ccvalues);
      }
      else {
        snl_fei::copy2DToColumnContig(numPtIndices, numPtColIndices, values, format,
                                      ccvalues);
      }

      int ptRowOffset = 0;
      for(i=0; i<numIDs; ++i) {

        if (rowConn[i]->getOwnerProc() == localProc()) {

          int numColIDs = numIDs;
          int* colNodeNums = nodeNumbers;
          const int* colDims = numIndicesPerID;
          int LDA = numPtColIndices;
          if (format == FEI_BLOCK_DIAGONAL_ROW) {
            numColIDs = 1;
            colNodeNums = &(nodeNumbers[i]);
            colDims = &(numIndicesPerID[i]);
            LDA = numIndicesPerID[i];
          }

          CHK_ERR( snl_fei::BlockMatrixTraits<T>::sumIn(matrix_.get(),
                                               nodeNumbers[i],
                                               numIndicesPerID[i],
                                               numColIDs, colNodeNums,
                                               colDims, LDA,
                                               &(ccvalues[ptRowOffset])) );
          changedSinceMark_ = true;
        }
        else {
          if (output_level_ >= fei::FULL_LOGS && output_stream_ != NULL) {
            FEI_OSTREAM& os = *output_stream_;
            for(int ii=0; ii<numIndicesPerID[i]; ++ii) {
              os << "#   remote pt-row " << ptIndices[ptRowOffset]+ii <<" ";
              for(int jj=0; jj<numPtIndices; ++jj) {
                os << "("<<ptIndices[jj]<<","<<values[ptRowOffset+ii][jj]<<") ";
              }
              os << FEI_ENDL;
            }
          }

          for(int ii=0; ii<numIndicesPerID[i]; ++ii) {
            int p=eqnComm_->getOwnerProc(ptIndices[ptRowOffset]+ii);
            FillableMat* remote_mat = getRemotelyOwnedMatrix(p);
            remote_mat->sumInRow(ptIndices[ptRowOffset]+ii,
                                      ptIndices,
                                      values[ptRowOffset+ii],
                                      numPtIndices);
            changedSinceMark_ = true;
          }
        }

        ptRowOffset += numIndicesPerID[i];
      }

      delete [] ccvalues;
    }

    return(0);
  }

  int numColIndices = symmetric ? numRowIndices : colpattern->getNumIndices();
  int* colIndices = rowIndices;
  const fei::Record<int>*const* colConn = NULL;

  if (!symmetric) {
    colConn = cblock->getColConnectivity(connectivityID);
    mgraph->getColSpace()->getGlobalIndices(colpattern,
                                            colConn, work_indices_);
    colIndices = &work_indices_[0];
  }

  if (symmetric) {
    if (format == FEI_DENSE_ROW || format == FEI_DENSE_COL) {
      CHK_ERR( sumIn(numRowIndices, rowIndices, numRowIndices, rowIndices,
                     values, format) );
    }
    else if (format == FEI_BLOCK_DIAGONAL_ROW) {
      int totalnumfields = pattern->getTotalNumFields();
      const int* fieldIDs = pattern->getFieldIDs();
      fei::SharedPtr<fei::VectorSpace> rowspace = mgraph->getRowSpace();
      fei::VectorSpace* rowspaceptr = rowspace.get();
      int ioffset = 0;
      for(int i=0; i<totalnumfields; ++i) {
        int fieldsize = rowspaceptr->getFieldSize(fieldIDs[i]);
        if (ioffset+fieldsize > numRowIndices) {
          FEI_OSTRINGSTREAM osstr;
          osstr<<"snl_fei::sumIn, format=FEI_BLOCK_DIAGONAL_ROW, block-sizes"
               << " not consistent with total num-indices."<<FEI_ENDL;
          throw std::runtime_error(osstr.str());
        }

        CHK_ERR( sumIn(fieldsize, &(rowIndices[ioffset]),
                       fieldsize, &(rowIndices[ioffset]),
                       &(values[ioffset]), FEI_DENSE_ROW) );
        ioffset += fieldsize;
      }
    }
    else {
      FEI_OSTRINGSTREAM osstr;
      osstr << "fei::Matrix_Impl::sumIn, format="<<format<<" not supported."
            << FEI_ENDL;
      throw std::runtime_error(osstr.str());
    }
  }
  else {
    if (format != FEI_DENSE_ROW && format != FEI_DENSE_COL) {
      FEI_OSTRINGSTREAM osstr;
      osstr << "fei::Matrix_Impl::sumIn, format="<<format<<" not valid with"
            << " un-symmetric matrix contributions."<<FEI_ENDL;
      throw std::runtime_error(osstr.str());
    }

    CHK_ERR( sumIn(numRowIndices, rowIndices, numColIndices, colIndices,
                   values, format) );
  }

  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::multiply(fei::Vector* x,
                                 fei::Vector* y)
{
  return( fei::MatrixTraits<T>::matvec(matrix_.get(), x, y) );
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::gatherFromOverlap(bool accumulate)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    (*output_stream_) << dbgprefix_<<"gatherFromOverlap"<<FEI_ENDL;
  }

  CHK_ERR( Matrix_core::gatherFromOverlap(accumulate) );

  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::giveToMatrix(int numRows, const int* rows,
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
    copyTransposeToWorkArrays(numRows, numCols, values,
                              work_data1D_, work_data2D_);
    myvalues = &work_data2D_[0];
  }

  if ((int)work_ints_.size() < numRows) {
    work_ints_.resize(numRows);
  }

  if (haveBlockMatrix()) {
    return( giveToBlockMatrix(numRows, rows, numCols, cols,
                              myvalues, sumInto) );
  }

  int i; 
  int numRemote = 0;
  int* workIntPtr = &work_ints_[0];
  for(i=0; i<numRows; ++i) {
    int row = rows[i];
    if (row < firstLocalOffset() || row > lastLocalOffset()) {
      ++numRemote;
      workIntPtr[i] = 1;
    }
    else {
      workIntPtr[i] = 0;
    }
  }

  if (numRemote < 1) {
    int err = giveToUnderlyingMatrix(numRows, rows, numCols, cols, myvalues,
                                    sumInto, FEI_DENSE_ROW);
    if (err != 0) {
      FEI_OSTRINGSTREAM osstr;
      osstr << "fei::Matrix_Impl::giveToMatrix ERROR: err="<<err
        << " returned from giveToUnderlyingMatrix.";
      throw std::runtime_error(osstr.str());
    }
    return(0);
  }

  for(i=0; i<numRows; ++i) {
    int row = rows[i];
    const double*const rowvalues = myvalues[i];

    if (workIntPtr[i] > 0) {
      int proc = eqnComm_->getOwnerProc(row);
      FillableMat* remote_mat = getRemotelyOwnedMatrix(proc);

      if (output_level_ >= fei::FULL_LOGS && output_stream_ != NULL) {
        FEI_OSTREAM& os = *output_stream_;
        os << dbgprefix_<<" remote_mat["<<proc<<"]: ";
        for(int jj=0; jj<numCols; ++jj) {
          os << "("<<row<<","<<cols[jj]<<","<<rowvalues[jj]<<") ";
        }
        os << FEI_ENDL;
      }

      if (sumInto) {
        remote_mat->sumInRow(row, cols, rowvalues, numCols);
        changedSinceMark_ = true;
      }
      else {
        remote_mat->putRow(row, cols, rowvalues, numCols);
        changedSinceMark_ = true;
      }

    }
    else {
      CHK_ERR( giveToUnderlyingMatrix(1, &row, numCols, cols, &rowvalues,
                                      sumInto, 0) );
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::giveToBlockMatrix(int numRows, const int* rows,
                                          int numCols, const int* cols,
                                          const double* const* values,
                                          bool sumInto)
{
  if (output_level_ >= fei::FULL_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << "#  giveToBlockMatrix numRows: " << numRows
       << ", numCols: " << numCols << FEI_ENDL;
  }

  if (numRows == 0 || numCols == 0) {
    return(0);
  }

  snl_fei::PointBlockMap* pointBlockMap = vecSpace()->getPointBlockMap();

  if (sumInto && numProcs() == 1) {
    std::vector<int> temp((numRows+numCols)*2);
    int* tempPtr = &temp[0];
    int* blkRows = tempPtr;
    int* blkRowOffsets = tempPtr+numRows;
    int* blkCols = blkRowOffsets+numRows;
    int* blkColOffsets = blkCols+numCols;

    CHK_ERR( convertPtToBlk(numRows, rows, numCols, cols,
                            blkRows, blkRowOffsets,
                            blkCols, blkColOffsets) );

    std::vector<int> blockRows, blockRowSizes;
    std::vector<int> blockCols, blockColSizes;
    for(int i=0; i<numRows; ++i) fei::sortedListInsert(blkRows[i], blockRows);
    for(int i=0; i<numCols; ++i) fei::sortedListInsert(blkCols[i], blockCols);

    int rowSizeTotal = 0, colSizeTotal = 0;

    for(size_t i=0; i<blockRows.size(); ++i) {
      int size = pointBlockMap->getBlkEqnSize(blockRows[i]);
      blockRowSizes.push_back(size);
      rowSizeTotal += size;
    }
    for(size_t i=0; i<blockCols.size(); ++i) {
      int size = pointBlockMap->getBlkEqnSize(blockCols[i]);
      blockColSizes.push_back(size);
      colSizeTotal += size;
    }
    std::vector<double> coefs_1d(rowSizeTotal*colSizeTotal, 0.0);
    double* coefs1dPtr = &coefs_1d[0];
    std::vector<double*> coefs_2d(blockRows.size()*blockCols.size());
    double** coefs2dPtr = &coefs_2d[0];

    int blkCounter = 0;
    int offset = 0;
    for(size_t i=0; i<blockRows.size(); ++i) {
      for(size_t j=0; j<blockCols.size(); ++j) {
        coefs2dPtr[blkCounter++] = &(coefs1dPtr[offset]);
        offset += blockRowSizes[i]*blockColSizes[j];
      }
    }

    for(int i=0; i<numRows; ++i) {
      int rowind = fei::binarySearch(blkRows[i], blockRows);
      int rowsize = blockRowSizes[rowind];

      for(int j=0; j<numCols; ++j) {
        int colind = fei::binarySearch(blkCols[j], blockCols);
        int pos = blkColOffsets[j]*rowsize + blkRowOffsets[i];

        coefs2dPtr[rowind*blockCols.size()+colind][pos] += values[i][j];
      }
    }

    for(size_t i=0; i<blockRows.size(); ++i) {
      CHK_ERR( giveToUnderlyingBlockMatrix(blockRows[i],
                                           blockRowSizes[i],
                                           blockCols.size(),
                                           &blockCols[0],
                                           &blockColSizes[0],
                                           &blockColSizes[0],
                                           &coefs2dPtr[i*blockCols.size()],
                                           true) );
    }

    return(0);
  }

  int maxBlkEqnSize = pointBlockMap->getMaxBlkEqnSize();
  int coefBlkLen = maxBlkEqnSize*maxBlkEqnSize*2;

  for(int i=0; i<numRows; ++i) {
    int row = rows[i];

    if (row < firstLocalOffset() || row > lastLocalOffset()) {
      int proc = eqnComm_->getOwnerProc(row);
      FillableMat* remote_mat = getRemotelyOwnedMatrix(proc);
      if (sumInto) {
        remote_mat->sumInRow(row, cols, values[i], numCols);
        changedSinceMark_ = true;
      }
      else {
        remote_mat->putRow(row, cols, values[i], numCols);
        changedSinceMark_ = true;
      }
      continue;
    }

    int blockRow = pointBlockMap->eqnToBlkEqn(row);
    int blockRowSize = pointBlockMap->getBlkEqnSize(blockRow);
    int blockRowOffset = pointBlockMap->getBlkEqnOffset(blockRow, row);

    int blockRowLength = 0;
    CHK_ERR( snl_fei::BlockMatrixTraits<T>::getRowLength(matrix_.get(), blockRow,
                                                blockRowLength) );

    std::vector<int> blkCols(blockRowLength);
    int* blkCols_ptr = &blkCols[0];
    std::vector<int> blkColDims(blockRowLength);
    int* blkColDims_ptr = &blkColDims[0];
    std::vector<double> coefs_1D(blockRowLength*coefBlkLen);
    double* coefs_1D_ptr = &coefs_1D[0];
    int coefsLen = coefs_1D.size();
    std::vector<double*> coefs_2D(blockRowLength);
    double** coefs_2D_ptr = &coefs_2D[0];

    std::vector<int> LDAs(blockRowLength);
    std::fill(LDAs.begin(), LDAs.end(), blockRowSize);
    std::fill(coefs_1D.begin(), coefs_1D.end(), 0.0);

    int checkRowLen = 0;
    CHK_ERR( snl_fei::BlockMatrixTraits<T>::copyOutRow(matrix_.get(),
                                              blockRow, blockRowLength,
                                              blockRowSize,
                                              blkCols_ptr,
                                              blkColDims_ptr,
                                              coefs_1D_ptr,
                                              coefsLen,
                                              checkRowLen) );
    int coefs_1D_offset = 0;
    for(int j=0; j<checkRowLen; ++j) {
      coefs_2D_ptr[j] = &(coefs_1D_ptr[coefs_1D_offset]);
      coefs_1D_offset += blockRowSize*blkColDims_ptr[j];
    }

    for(int j=0; j<numCols; ++j) {
      int blockCol = pointBlockMap->eqnToBlkEqn(cols[j]);
      int blkOffset= pointBlockMap->getBlkEqnOffset(blockCol, cols[j]);

      for(int jj=0; jj<blockRowLength; ++jj) {

        if (blockCol == blkCols_ptr[jj]) {
          if (sumInto) {
            coefs_2D_ptr[jj][blkOffset*blockRowSize+blockRowOffset] += values[i][j];
          }
          else {
            coefs_2D_ptr[jj][blkOffset*blockRowSize+blockRowOffset] = values[i][j];
          }

          break;
        }
      }
    }

    //Now put the block-row back into the matrix
    CHK_ERR( giveToUnderlyingBlockMatrix(blockRow, blockRowSize,
                                         blockRowLength, blkCols_ptr,
                                         &LDAs[0],
                                         blkColDims_ptr,
                                         coefs_2D_ptr,
                                         false) );
  }

  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::writeToFile(const char* filename,
                                    bool matrixMarketFormat)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();
  if (mgraph.get() == NULL) {
    return(-1);
  }

  if (haveFEMatrix()) {
    return(0);//temporary no-op...
  }

  if (haveBlockMatrix()) {
    CHK_ERR( snl_fei::BlockMatrixTraits<T>::globalAssemble(matrix_.get()) );
  }

  static char mmbanner[] = "%%MatrixMarket matrix coordinate real general";

  fei::SharedPtr<fei::VectorSpace> vspace = mgraph->getRowSpace();

  int globalNNZ = 0;
  int globalNumRows = vspace->getGlobalNumIndices();
  int localNumRows = vspace->getNumIndices_Owned();

  fei::SharedPtr<fei::VectorSpace> cspace = mgraph->getColSpace();
  int globalNumCols = globalNumRows;
  if (cspace.get() != NULL) {
    globalNumCols = cspace->getGlobalNumIndices();
  }

  std::vector<int> indices_owned;
  int localNNZ = 0;
  CHK_ERR( vspace->getIndices_Owned(indices_owned) );
  int* rowsPtr = &indices_owned[0];
  for(int i=0; i<localNumRows; ++i) {
    int len;
    CHK_ERR( getRowLength(rowsPtr[i], len) );
    localNNZ += len;
  }

  CHK_MPI( GlobalSum(getCommunicator(), localNNZ, globalNNZ) );

  for(int p=0; p<numProcs(); ++p) {
    fei::Barrier(getCommunicator());
    if (p != localProc()) continue;

    FEI_OFSTREAM* outFile = NULL;
    if (p==0) {
      outFile = new FEI_OFSTREAM(filename, IOS_OUT);
      FEI_OFSTREAM& ofs = *outFile;
      if (matrixMarketFormat) {
        ofs << mmbanner << FEI_ENDL;
        ofs <<globalNumRows << " " <<globalNumCols << " " <<globalNNZ <<FEI_ENDL;
      }
      else {
        ofs <<globalNumRows << " " <<globalNumCols <<FEI_ENDL;
      }
    }
    else outFile = new FEI_OFSTREAM(filename, IOS_APP);

    outFile->setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
    outFile->precision(13);
    FEI_OFSTREAM& ofs = *outFile;

    int rowLength;

    for(int i=firstLocalOffset(); i<=lastLocalOffset(); ++i) {
      CHK_ERR( getRowLength(i, rowLength) );

      work_indices_.resize(rowLength);
      work_data1D_.resize(rowLength);

      int* indPtr = &work_indices_[0];
      double* coefPtr = &work_data1D_[0];

      CHK_ERR( copyOutRow(i, rowLength, coefPtr, indPtr) );

      fei::insertion_sort_with_companions<double>(rowLength,
                                            indPtr, coefPtr);

      for(int j=0; j<rowLength; ++j) {
        if (matrixMarketFormat) {
          ofs << i+1 <<" "<<indPtr[j]+1<<" "<<coefPtr[j]<<FEI_ENDL;
        }
        else {
          ofs << i <<" "<<indPtr[j]<<" "<<coefPtr[j]<<FEI_ENDL;
        }
      }
    }

    delete outFile;
  }

  return(0);
}

//----------------------------------------------------------------------------
template<typename T>
int fei::Matrix_Impl<T>::writeToStream(FEI_OSTREAM& ostrm,
                                      bool matrixMarketFormat)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();
  if (mgraph.get() == NULL) {
    return(-1);
  }

  if (haveFEMatrix()) {
    return(0);//temporary no-op...
  }

  if (haveBlockMatrix()) {
    CHK_ERR( snl_fei::BlockMatrixTraits<T>::globalAssemble(matrix_.get()) );
  }

  static char mmbanner[] = "%%MatrixMarket matrix coordinate real general";

  fei::SharedPtr<fei::VectorSpace> vspace = mgraph->getRowSpace();

  int globalNNZ = 0;
  int localNumRows = vspace->getNumIndices_Owned();
  std::vector<int> indices_owned;
  int localNNZ = 0;
  CHK_ERR( vspace->getIndices_Owned(indices_owned));
  int* rowsPtr = &indices_owned[0];
  for(int i=0; i<localNumRows; ++i) {
    int len;
    CHK_ERR( getRowLength(rowsPtr[i], len) );
    localNNZ += len;
  }

  CHK_MPI( fei::GlobalSum(getCommunicator(), localNNZ, globalNNZ) );

  IOS_FMTFLAGS oldf = ostrm.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);

  for(int p=0; p<numProcs(); ++p) {
    fei::Barrier(getCommunicator());
    if (p != localProc()) continue;

    if (p==0) {
      int globalSize = globalOffsets()[numProcs()]-1;
      if (matrixMarketFormat) {
        ostrm << mmbanner << FEI_ENDL;
        ostrm << globalSize << " " << globalSize << " " << globalNNZ << FEI_ENDL;
      }
      else {
        ostrm << globalSize << " " << globalSize << FEI_ENDL;
      }
    }

    int rowLength;

    for(int i=firstLocalOffset(); i<=lastLocalOffset(); ++i) {
      CHK_ERR( getRowLength(i, rowLength) );

      work_indices_.resize(rowLength);
      work_data1D_.resize(rowLength);

      int* indPtr = &work_indices_[0];
      double* coefPtr = &work_data1D_[0];

      CHK_ERR( copyOutRow(i, rowLength, coefPtr, indPtr) );

      for(int j=0; j<rowLength; ++j) {
        ostrm << i << " " << indPtr[j] << " " << coefPtr[j] << FEI_ENDL;
      }
    }
  }

  ostrm.setf(oldf, IOS_FLOATFIELD);

  return(0);
}

#endif // _fei_Matrix_Impl_hpp_

