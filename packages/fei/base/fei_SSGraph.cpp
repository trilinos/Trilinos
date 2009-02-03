/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <feiArray.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_SSGraph.hpp>

#undef fei_file
#define fei_file "fei_SSGraph.cpp"
#include <fei_ErrMacros.hpp>

//------------------------------------------------------------------------------
SSGraph::SSGraph()
{
  rows_ = new feiArray<int>;
  rowLengths_ = new feiArray<int>;

  //The second constructor argument to feiArray is the 'allocIncrement', or the
  //amount by which the internal array grows when a new entry is added. We want
  //it to be 1 in this case, because each entry is a pointer-to-array which
  //represents a row of the table and we don't want to over-allocate the table.
  //The reason
  //for this is that we want to be able to clear and re-fill a matrix without
  //destroying and re-allocating memory, and without causing memory leaks. When
  //clearing the table, we set the length of each row to 0 without changing the
  //de-allocating the row, and we set the table's first dimension (num-rows) to
  //0, also without de-allocating. When we then come back to re-fill the matrix,
  //we need to know that the previous allocated-length contains not just 
  //pointers-to-rows, but pointers-to-allocated-rows. If the table's first 
  //dimension was over-allocated initially, we wouldn't know whether the
  //pointers-to-row were just un-initialized pointers, or actually allocated.

  indices_ = new feiArray<feiArray<int>*>(0,1);

  whichConstructor_ = SS_Constr_Default;
}

//------------------------------------------------------------------------------
SSGraph::SSGraph(int numRows, const int* rowNumbers,
	     int numCols, const int* colIndices)
{
  rows_ = new feiArray<int>(numRows, numRows, (int*)rowNumbers);
  rowLengths_ = new feiArray<int>(numRows, numRows);

  indices_ = new feiArray<feiArray<int>*>(numRows, numRows);

  (*(indices_))[0] = new feiArray<int>(numCols, numCols, (int*)colIndices);

  for(int i=0; i<numRows; i++) {
    (*rowLengths_)[i] = numCols;

    if (i > 0) (*(indices_))[i] = (*(indices_))[0];
  }

  whichConstructor_ = SS_Constr_RawArrays;
}

//------------------------------------------------------------------------------
SSGraph::SSGraph(int numRows, const int* rowNumbers,
	     int numColsPerRow, const int* rowColOffsets,
	     const int* colIndices)
{
  rows_ = new feiArray<int>(numRows, numRows, (int*)rowNumbers);
  rowLengths_ = new feiArray<int>(numRows, numRows);
  *rowLengths_ = numColsPerRow;

  indices_ = new feiArray<feiArray<int>*>(numRows, numRows);

  for(int i=0; i<numRows; i++) {
    int offset = rowColOffsets[i];

    (*indices_)[i] = new feiArray<int>(numColsPerRow, numColsPerRow,
				       (int*)(&(colIndices[offset])));
  }

  whichConstructor_ = SS_Constr_RawArrays2;
}

//------------------------------------------------------------------------------
SSGraph::~SSGraph()
{
  if (whichConstructor_ == SS_Constr_EqnBuf) return;

  int len = indices_->allocatedLength();
  indices_->resize(len);

  for(int i=0; i<len; i++) {
    if (whichConstructor_ == SS_Constr_Default || 
	whichConstructor_ == SS_Constr_RawArrays2 ||
	i == 0) delete (*indices_)[i];
  }

  delete indices_;
  delete rows_;
  delete rowLengths_;
}


//------------------------------------------------------------------------------
void SSGraph::logicalClear()
{
  if (whichConstructor_ != SS_Constr_Default) {
    throw std::runtime_error("fei SSGraph::logicalClear: not constructed with SS_Constr_Default");
  }

  rows_->resize(0);
  rowLengths_->resize(0);

  for(int i=0; i<indices_->length(); i++) {
    (*(indices_))[i]->resize(0);
  }

  indices_->resize(0);
}

//------------------------------------------------------------------------------
void SSGraph::createPosition(int row, int col)
{
  int dummy1, dummy2;

  createPosition(row, col, dummy1, dummy2);
}

//------------------------------------------------------------------------------
void SSGraph::createPosition(int row, int col, int& rowIndex, int& colIndex)
{
  if (whichConstructor_ != SS_Constr_Default) {
    throw std::runtime_error("fei SSGraph::createPosition: not constructed with SS_Constr_Default");
  }

  int insertPoint = -1;
  rowIndex = snl_fei::binarySearch(row, *rows_, insertPoint);

  if (rowIndex < 0) {
    insertRow(row, insertPoint);
    rowIndex = insertPoint;
  }

  colIndex = snl_fei::binarySearch(col, *((*indices_)[rowIndex]), insertPoint);

  if (colIndex < 0) {
    (*indices_)[rowIndex]->insert(col, insertPoint);
    (*(rowLengths_))[rowIndex]++;
    colIndex = insertPoint;
  }
}

//------------------------------------------------------------------------------
void SSGraph::appendRow(int row)
{
  //If the default constructor wasn't used to create this SSGraph instance, then
  //we can't alter the data.
  if (whichConstructor_ != SS_Constr_Default) {
    throw std::runtime_error("fei SSGraph::createPosition: not constructed with SS_Constr_Default");
  }

  //append the new row-number in the rows_ array.
  rows_->append(row);
  rowLengths_->append(0);

  int allocLen = indices_->allocatedLength();
  int numRows = indices_->length();

  //allocLen is the number of rows that have ever been allocated for this
  //matrix; it is the "high-water-mark" of num-rows for this matrix. (If
  //indices_->length() < indices_->allocatedLength(), it means that the
  //'logicalClear()' function has been called. In that case, we simply 
  //"re-activate" a previously allocated row rather than allocating a new one.

  if (numRows < allocLen) {
    indices_->resize(numRows+1);
  }
  else {
    //append a new feiArray<int>* to the indices_ table.
    indices_->append(new feiArray<int>);
  }
}

//------------------------------------------------------------------------------
void SSGraph::insertRow(int row, int index)
{
  //If the default constructor wasn't used to create this SSGraph instance, then
  //we can't alter the data.
  if (whichConstructor_ != SS_Constr_Default) {
    throw std::runtime_error("fei SSGraph::createPosition: not constructed with SS_Constr_Default");
  }

  //insert the new row-number in the rows_ array.
  rows_->insert(row, index);
  rowLengths_->insert(0, index);

  int allocLen = indices_->allocatedLength();
  int numRows = indices_->length();

  //allocLen is the number of rows that have ever been allocated for this
  //matrix; it is the "high-water-mark" of num-rows for this matrix. (If
  //indices_->length() < indices_->allocatedLength(), it means that the
  //'logicalClear()' function has been called. In that case, we simply 
  //"re-activate" a previously allocated row rather than allocating a new one.

  if (numRows < allocLen) {
    indices_->resize(numRows+1);
    feiArray<feiArray<int>*>& indRef = *indices_;
    feiArray<int>* tmpI = indRef[numRows];

    for(int i=numRows; i>index; i--) {
      indRef[i] = indRef[i-1];
    }

    indRef[index] = tmpI;
  }
  else {
    //append a new feiArray<int>* to the indices_ table.
    indices_->insert(new feiArray<int>, index);
  }
}



