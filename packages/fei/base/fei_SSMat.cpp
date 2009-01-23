/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <limits>
#include <cmath>

#include <feiArray.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_EqnBuffer.hpp>

#include <fei_SSVec.hpp>
#include <fei_SSMat.hpp>

#undef fei_file
#define fei_file "fei_SSMat.cpp"
#include <fei_ErrMacros.hpp>

//==============================================================================
SSMat::SSMat(int alloc_increment, int row_alloc_increment)
  : structurallySymmetric(false),
    whichConstructor_(SS_Constr_Default),
    rowNumbers_(NULL),
    rows_(NULL),
    highWaterMark_(0),
    row_alloc_incr_(row_alloc_increment)
{
  (void)alloc_increment;
  rowNumbers_ = new feiArray<int>(0, alloc_increment);

  //We want to be able to clear and re-fill a matrix without
  //destroying and re-allocating memory, and without causing memory leaks. When
  //clearing the rows, we will set the length of each row to 0 without
  //de-allocating the row, and we will set the table's first dimension (num-rows)
  //to 0, also without de-allocating. When we then come back to re-fill the
  //matrix, we need to know that the previous allocated-length contains not just 
  //pointers-to-rows, but pointers-to-allocated-rows. If the table's first 
  //dimension was over-allocated initially, we wouldn't know whether the
  //pointers-to-row were just un-initialized pointers, or actually allocated.

  rows_ = new feiArray<SSVec*>(0, alloc_increment);
}

//==============================================================================
SSMat::SSMat(const SSMat& src)
  : structurallySymmetric(src.structurallySymmetric),
    whichConstructor_(SS_Constr_Default),
    rowNumbers_(NULL),
    rows_(NULL),
    highWaterMark_(0),
    row_alloc_incr_(0)
{
  *this = src;
}

//==============================================================================
SSMat& SSMat::operator=(const SSMat& src)
{
  structurallySymmetric = src.structurallySymmetric;
  whichConstructor_ = SS_Constr_Default;
  delete rowNumbers_;
  SSVec** tmp_rows_ptr = rows_->dataPtr();
  for(int i=0; i<rows_->length(); ++i) {
    delete tmp_rows_ptr[i];
  }
  delete rows_;
  rowNumbers_ = new feiArray<int>(*(src.rowNumbers_));
  rows_ = new feiArray<SSVec*>(src.rowNumbers_->length(), 64);
  row_alloc_incr_ = src.row_alloc_incr_;
  feiArray<SSVec*>& srcrows = *(src.rows_);
  SSVec** rowsPtr = rows_->dataPtr();
  for(int i=0; i<src.rowNumbers_->length(); ++i) {
    rowsPtr[i] = new SSVec(*(srcrows[i]));
  }

  int alloclen = rows_->allocatedLength();
  if (alloclen > rows_->length()) {
    rowsPtr = rows_->dataPtr();
    for(int j=rows_->length(); j<alloclen; ++j) {
      rowsPtr[j] = NULL;
    }
  }

  highWaterMark_ = rows_->length();

  return(*this);
}

//==============================================================================
SSMat& SSMat::operator+=(const SSMat& src)
{
  whichConstructor_ = SS_Constr_Default;
  feiArray<int>& rowNumbers = *(src.rowNumbers_);
  feiArray<SSVec*>& srcrows = *(src.rows_);

  for(int i=0; i<rowNumbers.length(); ++i) {
    SSVec* row = srcrows[i];
    sumInRow(rowNumbers[i], row->indices().dataPtr(),
             row->coefs().dataPtr(), row->length());
  }

  highWaterMark_ = rows_->length();

  return(*this);
}

//==============================================================================
SSMat::SSMat(EqnBuffer& eqnBuf)
  : structurallySymmetric(false),
    whichConstructor_(SS_Constr_EqnBuf),
    rowNumbers_(NULL),
    rows_(NULL),
    highWaterMark_(0),
    row_alloc_incr_(0)
{
  rowNumbers_ = &(eqnBuf.eqnNumbersPtr());
  rows_ = &(eqnBuf.eqns());
}

//==============================================================================
SSMat::~SSMat()
{
  if (whichConstructor_ == SS_Constr_EqnBuf) return;

  if (rows_->length() > highWaterMark_) {
    highWaterMark_ = rows_->length();
  }

  int len = highWaterMark_;
  rows_->resize(len);

  if (whichConstructor_ != SS_Constr_EqnBuf) {
    for(int i=0; i<len; i++) {
      delete (*rows_)[i]; (*rows_)[i] = NULL;
    }
  }

  delete rows_; rows_ = NULL;
  delete rowNumbers_; rowNumbers_ = NULL;
}


//==============================================================================
void SSMat::logicalClear()
{
  if (whichConstructor_ != SS_Constr_Default) {
    throw std::runtime_error("fei SSMat::logicalClear, not constructed with SS_Constr_Default.");
  }

  if (rows_->length() > highWaterMark_) { 
    highWaterMark_ = rows_->length();
  }

  rowNumbers_->resize(0);
  SSVec** rowsPtr = rows_->dataPtr();

  for(int i=0; i<rows_->length(); i++) {
    rowsPtr[i]->indices().resize(0);
    rowsPtr[i]->coefs().resize(0);
  }

  rows_->resize(0);
}

//------------------------------------------------------------------------------
int SSMat::numNonzeros()
{
  int nnz = 0;
  SSVec** rows = rows_->dataPtr();
  for(int i=0; i<rows_->length(); ++i) {
    nnz += rows[i]->length();
  }

  return(nnz);
}

//==============================================================================
int SSMat::getMinCol()
{
  int numrows = rows_->length();
  if (numrows < 1) return(-1);

  SSVec** rows = rows_->dataPtr();
  int mincol = std::numeric_limits<int>::max();

  for(int i=0; i<numrows; ++i) {
    int rowlen = rows[i]->length();
    if (rowlen < 1) continue;

    int thiscol = rows[i]->indices().dataPtr()[0];
    if (thiscol < mincol) mincol = thiscol;
  }

  return( mincol );
}

//==============================================================================
int SSMat::getMaxCol()
{
  int numrows = rows_->length();
  if (numrows < 1) return(-1);

  SSVec** rows = rows_->dataPtr();
  int maxcol = std::numeric_limits<int>::min();

  for(int i=0; i<numrows; ++i) {
    int rowlen = rows[i]->length();
    if (rowlen < 1) continue;

    int thiscol = rows[i]->indices().dataPtr()[rowlen-1];
    if (thiscol > maxcol) maxcol = thiscol;
  }

  return( maxcol );
}

//==============================================================================
int SSMat::matMat(SSMat& inMat, SSMat& result, bool storeResultZeros)
{
  feiArray<int>& inRowNumbers = inMat.getRowNumbers();
  int* inRowNumbersPtr = inRowNumbers.dataPtr();
  int inRowNumbersLen = inRowNumbers.length();
  SSVec** inRows = inMat.getRows().dataPtr();

  bool inMatSorted = false;
  if (inMat.whichConstructor_==SS_Constr_Default ||
      inMat.whichConstructor_==SS_Constr_EqnBuf) inMatSorted = true;

  //clear the result matrix without deleting any of its memory.
  result.logicalClear();

  SSVec** rowsPtr = rows_->dataPtr();
  int* rowNumbersPtr = rowNumbers_->dataPtr();

  //double fei_eps = std::numeric_limits<double>::epsilon();
  double fei_eps = 1.e-49;

  //loop down the rows of 'this' matrix
  for(int i=0; i<rowNumbers_->length(); i++) {
    int row = rowNumbersPtr[i];

    feiArray<int>& indicesRow = rowsPtr[i]->indices();
    int* indPtr = indicesRow.dataPtr();
    feiArray<double>& coefRow = rowsPtr[i]->coefs();
    double* coefPtr = coefRow.dataPtr();

    SSVec* resultRow = result.getRow(row, true);

    for(int j=0; j<indicesRow.length(); j++) {
      int col = indPtr[j];
      double coef = coefPtr[j];

      int rindex;
      if (inMatSorted) {
        rindex = snl_fei::binarySearch(col, inRowNumbersPtr,
                                       inRowNumbersLen);
      }
      else rindex = inRowNumbers.find(col);

      if (rindex < 0) continue;

      feiArray<int>& inIndicesRow = inRows[rindex]->indices();
      int* inIndPtr = inIndicesRow.dataPtr();
      feiArray<double>& inCoefsRow = inRows[rindex]->coefs();
      double* inCoefPtr = inCoefsRow.dataPtr();

      for(int k=0; k<inIndicesRow.length(); k++) {
	double resultCoef = coef * inCoefPtr[k];
	int resultCol = inIndPtr[k];

	//result(row, resultCol) += this(row, col)*inMat(col, resultCol)

        if (std::abs(resultCoef) > fei_eps || storeResultZeros) {
          resultRow->addEntry(resultCol, resultCoef);
        }
      }
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
SSVec* SSMat::rowContainingCol(int col, int& offsetInRow)
{
  SSVec** rows = rows_->dataPtr();
  for(int i=0; i<rows_->length(); ++i) {
    SSVec* row = rows[i];
    int* indices = row->indices().dataPtr();
    int len = row->length();
    int offset = snl_fei::binarySearch(col, indices, len);
    if (offset > -1) {
      offsetInRow = offset;
      return(row);
    }
  }

  return(NULL);
}

//==============================================================================
int SSMat::matTransMat(SSMat& inMat, SSMat& result, bool storeResultZeros)
{
  //clear the result matrix without deleting any of its memory.
  result.logicalClear();

  feiArray<int>& inRowNumbers = inMat.getRowNumbers();
  if (inRowNumbers.length() < 1) {
    return(0);
  }

  int* inRowNumbersPtr = inRowNumbers.dataPtr();
  int inRowNumbersLen = inRowNumbers.length();
  SSVec** inRows = inMat.getRows().dataPtr();

  bool inMatSorted = false;
  if (inMat.whichConstructor_==SS_Constr_Default ||
      inMat.whichConstructor_==SS_Constr_EqnBuf) inMatSorted = true;

  SSVec** rowsPtr = rows_->dataPtr();
  int* rowNumbersPtr = rowNumbers_->dataPtr();
  int numrows = rowNumbers_->length();

  int startsearch = 0, dummy;
  int endsearch = inRowNumbersLen-1;
  int rindex = 0;
  feiArray<double> resultCoefs;

  //loop down the rows of 'this' matrix
  for(int i=0; i<numrows; i++) {
    int row = rowNumbersPtr[i];

    int inRow = rindex < 0 ? -1 : inRowNumbersPtr[rindex];

    if (row < inRow) continue;

    if (row != inRow) {
      if (inMatSorted) {
	rindex = snl_fei::binarySearch(row, inRowNumbersPtr, inRowNumbersLen,
				       startsearch, endsearch, dummy);
      }
      else rindex = inRowNumbers.find(row);
      if (rindex < 0) continue;
    }

    startsearch = rindex+1;

    int rowlen = rowsPtr[i]->length();
    int* indPtr = rowsPtr[i]->indices().dataPtr();
    double* coefPtr = rowsPtr[i]->coefs().dataPtr();

    feiArray<int>& inIndicesRow = inRows[rindex]->indices();
    int inIndicesRowLen = inIndicesRow.length();
    int* inIndPtr = inIndicesRow.dataPtr();
    double* inCoefPtr = inRows[rindex]->coefs().dataPtr();

    resultCoefs.resize(inIndicesRowLen);
    double* resultCoefsPtr = resultCoefs.dataPtr();

    for(int j=0; j<rowlen; j++) {
      int col = indPtr[j];
      double coef = coefPtr[j];

      SSVec* resultrow = result.getRow(col, true);

      for(int k=0; k<inIndicesRowLen; k++) {
	resultCoefsPtr[k] = coef * inCoefPtr[k];
      }

      resultrow->addEntries_sortedInput(inIndicesRowLen, resultCoefsPtr,
					inIndPtr, storeResultZeros);
    }
  }

  return(0);
}

//==============================================================================
int SSMat::matVec(SSVec& inVec, SSVec& result)
{
  feiArray<int>& resultIndices = result.indices();
  feiArray<double>& resultCoefs = result.coefs();

  resultIndices.resize(0);
  resultCoefs.resize(0);

  feiArray<int>& inVecIndices = inVec.indices();
  int* inVecIndicesPtr = inVecIndices.dataPtr();
  int inVecIndicesLen = inVecIndices.length();

  feiArray<double>& inVecCoefs = inVec.coefs();
  double* inVecCoefPtr = inVecCoefs.dataPtr();
  bool inVecSorted = inVec.whichConstructor_==SS_Constr_Default ? true : false;
  SSVec** rowsPtr = rows_->dataPtr();
  int* rowNumbersPtr = rowNumbers_->dataPtr();

  for(int i=0; i<rowNumbers_->length(); i++) {
    feiArray<int>& indicesRow = rowsPtr[i]->indices();
    int* indPtr = indicesRow.dataPtr();
    feiArray<double>& coefsRow = rowsPtr[i]->coefs();
    double* coefPtr = coefsRow.dataPtr();

    double sum = 0.0;
    bool entry = false;

    for(int j=0; j<indicesRow.length(); j++) {
      int col = indPtr[j];
      double coef = coefPtr[j];

      int index;
      if (inVecSorted) {
        index = snl_fei::binarySearch(col, inVecIndicesPtr,
                                      inVecIndicesLen);
      }
      else index = inVecIndices.find(col);
      if (index >= 0) {
	sum += coef*inVecCoefPtr[index];
	entry = true;
      }
    }

    if (entry) {
      resultIndices.append(rowNumbersPtr[i]);
      resultCoefs.append(sum);
    }
  }

  return(0);
}

//==============================================================================
int SSMat::matTransVec(SSVec& inVec, SSVec& result)
{
  feiArray<int>& resultIndices = result.indices();
  resultIndices.resize(0);
  int* resultIndicesPtr = resultIndices.dataPtr();
  int resultIndLen = resultIndices.length();
  feiArray<double>& resultCoefs = result.coefs();
  resultCoefs.resize(0);
  double* resultCoefsPtr = resultCoefs.dataPtr();

  feiArray<int>& inVecIndices = inVec.indices();
  int* inVecIndicesPtr = inVecIndices.dataPtr();
  int inVecIndicesLen = inVecIndices.length();

  feiArray<double>& inVecCoefs = inVec.coefs();
  double* inVecCoefPtr = inVecCoefs.dataPtr();
  bool inVecSorted = inVec.whichConstructor_==SS_Constr_Default ? true : false;
  int inVecLen = inVecIndices.length();

  int startsearch = 0, dummy;
  int endsearch = inVecLen-1;
  SSVec** rowsPtr = rows_->dataPtr();
  int* rowNumbersPtr = rowNumbers_->dataPtr();

  int index = 0, insertPoint;

  for(int i=0; i<rowNumbers_->length(); i++) {
    int row = rowNumbersPtr[i];

    int inVecInd = index < 0 ? -1 : inVecIndicesPtr[index];

    if (row < inVecInd) continue;

    if (row != inVecInd) {
      if (inVecSorted) index = snl_fei::binarySearch(row, inVecIndicesPtr,
						     inVecIndicesLen,
						     startsearch, endsearch, dummy);
      else index = inVecIndices.find(row);
      if (index < 0) continue;
    }

    startsearch = index+1;

    feiArray<int>& indicesRow = rowsPtr[i]->indices();
    int* indPtr = indicesRow.dataPtr();
    feiArray<double>& coefsRow = rowsPtr[i]->coefs();
    double* coefPtr = coefsRow.dataPtr();

    double inVecCoef = inVecCoefPtr[index];

    for(int j=0; j<indicesRow.length(); j++) {
      int col = indPtr[j];
      double coef = coefPtr[j] * inVecCoef;

      int thisindex = snl_fei::binarySearch(col, resultIndicesPtr,
                                        resultIndLen, insertPoint);
      if (thisindex >= 0) resultCoefsPtr[thisindex] += coef;
      else {
	resultIndices.insert(col, insertPoint);
	resultCoefs.insert(coef, insertPoint);
	resultIndicesPtr = resultIndices.dataPtr();
	resultIndLen = resultIndices.length();
	resultCoefsPtr = resultCoefs.dataPtr();
      }
    }
  }

  return(0);
}

//==============================================================================
int SSMat::coefficientPointer(int row, int col, double*& coefPtr)
{
  int index = rowNumbers_->find(row);

  if (index < 0) return(1);

  feiArray<int>& indicesRow = (*rows_)[index]->indices();

  int colIndex = indicesRow.find(col);

  if (colIndex < 0) return(1);

  feiArray<double>& coefRow = (*rows_)[index]->coefs();

  coefPtr = &( coefRow[colIndex] );
  return(0);
}

//==============================================================================
int SSMat::sumInCoef(int row, int col, double coef)
{
  if (whichConstructor_ != SS_Constr_Default) return(-1);

  int rowIndex, colIndex;

  createPosition(row, col, rowIndex, colIndex);

  (*rows_)[rowIndex]->coefs()[colIndex] += coef;

  return(0);
}

//==============================================================================
int SSMat::putCoef(int row, int col, double coef)
{
  if (whichConstructor_ != SS_Constr_Default) return(-1);

  int rowIndex, colIndex;

  createPosition(row, col, rowIndex, colIndex);

  (*rows_)[rowIndex]->coefs()[colIndex] = coef;

  return(0);
}

//==============================================================================
int SSMat::sumInRow(int row, const int* cols, const double* coefs, int len)
{
  SSVec* ssrow = getRow(row, true);
  ssrow->addEntries(len, coefs, cols);
  return(0);
}

//==============================================================================
int SSMat::putRow(int row, const int* cols, const double* coefs, int len)
{
  SSVec* ssrow = getRow(row, true);
  ssrow->putEntries(len, coefs, cols);
  return(0);
}

//==============================================================================
void SSMat::createPosition(int row, int col)
{
  int dummy1, dummy2;

  createPosition(row, col, dummy1, dummy2);
}

//==============================================================================
SSVec* SSMat::getRow(int row, bool create_if_necessary)
{
  bool matIsSorted = whichConstructor_==SS_Constr_Default ? true : false;

  int index = -1, insertPoint = -1;
  if (matIsSorted) {
    index = snl_fei::binarySearch(row, rowNumbers_->dataPtr(),
                                  rowNumbers_->length(), insertPoint);
  }
  else {
    index = rowNumbers_->find(row);
  }

  if (index < 0) {
    if (create_if_necessary && insertPoint > -1) {
      SSVec* newrow = insertRow(row, insertPoint);
      return(newrow);
    }
    return(NULL);
  }
  else {
    return(rows_->dataPtr()[index]);
  }
}

//==============================================================================
void SSMat::createPosition(int row, int col, int& rowIndex, int& colIndex)
{
  if (whichConstructor_ != SS_Constr_Default) {
    throw std::runtime_error("fei SSMat::createPosition: not constructed with SS_Constr_Default.");
  }

  int insertPoint = -1;
  rowIndex = snl_fei::binarySearch(row, rowNumbers_->dataPtr(),
                                   rowNumbers_->length(), insertPoint);

  if (rowIndex < 0) {
    insertRow(row, insertPoint);
    rowIndex = insertPoint;
  }

  SSVec* matrixrow = rows_->dataPtr()[rowIndex];
  feiArray<int>& indices = matrixrow->indices();
  feiArray<double>& coefs= matrixrow->coefs();

  colIndex = snl_fei::binarySearch(col, indices.dataPtr(),
                                   indices.length(), insertPoint);

  if (colIndex < 0) {
    indices.insert(col, insertPoint);
    coefs.insert(0.0, insertPoint);
    colIndex = insertPoint;
  }
}

//==============================================================================
void SSMat::appendRow(int row)
{
  //If the default constructor wasn't used to create this SSMat instance, then
  //we can't alter the data.
  if (whichConstructor_ != SS_Constr_Default) {
    throw std::runtime_error("fei SSMat::appendRow: not constructed with SS_Constr_Default.");
  }

  //append the new row-number in the rowNumbers_ array.
  rowNumbers_->append(row);

  int allocLen = rows_->allocatedLength();
  int numRows = rows_->length();

  //allocLen is the number of rows that have ever been allocated for this
  //matrix; it is the "high-water-mark" of num-rows for this matrix. (If
  //indices_->length() < indices_->allocatedLength(), it means that the
  //'logicalClear()' function has been called. In that case, we simply 
  //"re-activate" a previously allocated row rather than allocating a new one.

  if (numRows < allocLen) {
    rows_->resize(numRows+1);
  }
  else {
    //append a new SSVec* to the rows_ table.
    rows_->append(new SSVec(row_alloc_incr_));
  }
}

//==============================================================================
SSVec* SSMat::insertRow(int row, int index)
{
  //If the default constructor wasn't used to create this SSMat instance, then
  //we can't alter the data.
  if (whichConstructor_ != SS_Constr_Default) {
    throw std::runtime_error("fei SSMat::insertRow: not constructed with SS_Constr_Default.");
  }

  //insert the new row-number in the rowNumbers_ array.
  rowNumbers_->insert(row, index);

  int numRows = rows_->length();

  //If rows_->length() < highWaterMark_, it means that the
  //'logicalClear()' function has been called. In that case, we simply 
  //"re-activate" a previously allocated row rather than allocating a new one.

  SSVec* returnValue = NULL;

  if (numRows < highWaterMark_) {
    rows_->resize(numRows+1);
    SSVec** rowsPtr = rows_->dataPtr();
    returnValue = rowsPtr[numRows];

    for(int i=numRows; i>index; i--) {
      rowsPtr[i] = rowsPtr[i-1];
    }

    rowsPtr[index] = returnValue;
  }
  else {
    //append a new SSVec* to the rows_ table.
    returnValue = new SSVec(row_alloc_incr_);
    rows_->insert(returnValue, index);
  }

  return(returnValue);
}

//----------------------------------------------------------------------------
bool SSMat::operator==(const SSMat& lhs)
{
  if (*rowNumbers_ != *(lhs.rowNumbers_)) return(false);

  SSVec** rows = rows_->dataPtr();
  SSVec** lhs_rows = lhs.rows_->dataPtr();

  for(int i=0; i<rows_->length(); ++i) {
    feiArray<int>& indices = rows[i]->indices();
    const feiArray<int>& lhsindices = lhs_rows[i]->indices();
    feiArray<double>& coefs = rows[i]->coefs();
    const feiArray<double>& lhscoefs = lhs_rows[i]->coefs();

    if (indices != lhsindices || coefs != lhscoefs) {
      return(false);
    }
  }

  return(true);
}

//----------------------------------------------------------------------------
bool SSMat::operator!=(const SSMat& lhs)
{
  return( !(*this == lhs) );
}

//----------------------------------------------------------------------------
void SSMat::writeToStream(FEI_OSTREAM& os)
{
  int nnz = 0;
  int i, j;
  feiArray<int> cols;
  for(i=0; i<rowNumbers_->length(); ++i) {
    SSVec& row = *((*rows_)[i]);
    nnz += row.length();
    int* colindices = row.indices().dataPtr();
    for(j=0; j<row.length(); ++j) {
      snl_fei::sortedListInsert(colindices[j], cols);
    }
  }

  int numCols = cols.length();
  int numRows = rowNumbers_->length();
  int* rowNumPtr = rowNumbers_->dataPtr();

  os << "SSMat sizes: "<<numRows << " " << numCols << " " << nnz << FEI_ENDL;
  for(i=0; i<numRows; ++i) {
    SSVec& row = *((*rows_)[i]);
    os << rowNumPtr[i] << ": ";
    for(j=0; j<row.length(); ++j) {
      os << "("<< row.indices()[j]<<","<<row.coefs()[j]<<") ";
    }
    os << FEI_ENDL;
  }
}
