/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_macros.hpp>

#include <fei_mpi.h>

#include <test_utils/test_AztecWrappers.hpp>

#include <fei_ArrayUtils.hpp>
#include <test_utils/LibraryFactory.hpp>

#ifdef HAVE_FEI_AZTECOO
#include <fei_Aztec_Map.hpp>
#include <fei_Aztec_LSVector.hpp>
#include <fei_AztecDMSR_Matrix.hpp>
#endif

#undef fei_file
#define fei_file "test_AztecWrappers.cpp"

#include <fei_ErrMacros.hpp>

#ifdef HAVE_FEI_AZTECOO
int compare_DMSR_contents(fei_trilinos::AztecDMSR_Matrix& matrix, int localOffset,
			  std::vector<std::vector<int> >& colIndices,
			  std::vector<std::vector<double> >& values);

int fill_DMSR(fei_trilinos::AztecDMSR_Matrix& matrix, int localOffset,
	      std::vector<std::vector<int> >& colIndices,
	      std::vector<std::vector<double> >& values, bool sumInto);
#endif

test_AztecWrappers::test_AztecWrappers(MPI_Comm comm)
 : tester(comm)
{
}

test_AztecWrappers::~test_AztecWrappers()
{
}

int test_AztecWrappers::runtests()
{
  if (numProcs_ < 2) {
    CHK_ERR( serialtest1() );
  }

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_AztecWrappers::serialtest1()
{
  return(0);
}

int test_AztecWrappers::test1()
{
#ifdef HAVE_FEI_AZTECOO
  int localSize = 3, globalSize = localSize*numProcs_;
  int localOffset = localSize*localProc_;
  int i;

  std::vector<int> update(localSize);
  for(i=0; i<localSize; i++) update[i] = localOffset+i;

  fei::SharedPtr<fei_trilinos::Aztec_Map> map(
    new fei_trilinos::Aztec_Map(globalSize, localSize, &update[0], localOffset, comm_));

  fei_trilinos::AztecDMSR_Matrix* matrix = new fei_trilinos::AztecDMSR_Matrix(map);

  std::vector<int> elemrows(localSize);
  std::vector<int> elemcols(globalSize);
  double** elemcoefs = new double*[localSize];
  for(int j=0; j<globalSize; ++j) elemcols[j] = j;
  for(i=0; i<localSize; ++i) {
    elemrows[i] = localOffset+i;
    elemcoefs[i] = new double[globalSize];
    for(int j=0; j<globalSize; ++j) {
      elemcoefs[i][j] = (double)(localOffset+i+j);
    }
  }
  
  std::vector<std::vector<int> > colIndices(localSize);
  std::vector<std::vector<double> > values(localSize);
  std::vector<int> rowLengths(localSize);
  std::vector<int*> colPtrs(localSize);
  int nnzeros = 0;

  for(i=0; i<localSize; i++) {
    int diagEntry = 0;
    int row = i+localOffset;
    for(int j=0; j<globalSize; j++) {
      int col = j;
      if (col == row) diagEntry = 1;
      colIndices[i].push_back(col);
      values[i].push_back((double)(row+col));
    }
    rowLengths[i] = colIndices[i].size() - diagEntry;
    nnzeros += rowLengths[i] + 1;
    colPtrs[i] = &(colIndices[i][0]);
  }

  matrix->allocate( &rowLengths[0] );

  if (!(matrix->isAllocated())) {
    ERReturn(-1);
  }

  if (matrix->getNumNonZeros() != nnzeros) {
    ERReturn(-1);
  }

  CHK_ERR( fill_DMSR(*matrix, localOffset, colIndices, values, true) );

  int* rowinds = &elemrows[0];
  int* colinds = &elemcols[0];

  CHK_ERR( matrix->sumIntoRow(localSize, rowinds, globalSize, colinds, elemcoefs) );

  for(i=0; i<localSize; ++i) {
    for(int j=0; j<globalSize; ++j) values[i][j] *= 2.0;
  }

  CHK_ERR( compare_DMSR_contents(*matrix, localOffset, colIndices, values) );

  for(i=0; i<localSize; ++i) {
    for(int j=0; j<globalSize; ++j) values[i][j] /= 2.0;
  }

  CHK_ERR( fill_DMSR(*matrix, localOffset, colIndices, values, false) );

  CHK_ERR( compare_DMSR_contents(*matrix, localOffset, colIndices, values) );

  if (matrix->writeToFile("A_Az_notFilled.mtx") != true) {
    ERReturn(-1);
  }

  if (matrix->readFromFile("A_Az_notFilled.mtx") != true) {
    ERReturn(-1);
  }

  CHK_ERR( compare_DMSR_contents(*matrix, localOffset, colIndices, values) );

  matrix->fillComplete();

  if (!(matrix->isFilled())) {
    ERReturn(-1);
  }

  if (matrix->writeToFile("A_Az_filled.mtx") != true) {
    ERReturn(-1);
  }

  CHK_ERR( compare_DMSR_contents(*matrix, localOffset, colIndices, values) );

  CHK_ERR( fill_DMSR(*matrix, localOffset, colIndices, values, false) );

  CHK_ERR( compare_DMSR_contents(*matrix, localOffset, colIndices, values) );

  matrix->put(0.0);

  CHK_ERR( fill_DMSR(*matrix, localOffset, colIndices, values, true) );

  CHK_ERR( matrix->sumIntoRow(localSize, rowinds, globalSize, colinds, elemcoefs) );

  for(i=0; i<localSize; ++i) {
    for(int j=0; j<globalSize; ++j) values[i][j] *= 2.0;
    delete [] elemcoefs[i];
  }
  delete [] elemcoefs;

  CHK_ERR( compare_DMSR_contents(*matrix, localOffset, colIndices, values) );

  if (matrix->writeToFile("A_Az_filled2.mtx") != true) {
    ERReturn(-1);
  }

  if (matrix->readFromFile("A_Az_filled2.mtx") != true) {
    ERReturn(-1);
  }

  CHK_ERR( compare_DMSR_contents(*matrix, localOffset, colIndices, values) );

  delete matrix;
#endif
  return(0);
}

int test_AztecWrappers::test2()
{
  return(0);
}

int test_AztecWrappers::test3()
{
  return(0);
}

int test_AztecWrappers::test4()
{
  return(0);
}

#ifdef HAVE_FEI_AZTECOO
//==============================================================================
int compare_DMSR_contents(fei_trilinos::AztecDMSR_Matrix& matrix, int localOffset,
			  std::vector<std::vector<int> >& colIndices,
			  std::vector<std::vector<double> >& values)
{
  int localSize = colIndices.size();

  for(int i=0; i<localSize; i++) {
    int row = i+localOffset;
    int rowLen = matrix.rowLength(row);
    if (rowLen == 0) ERReturn(-1);
    int* tempInd = new int[rowLen];
    double* tempVal = new double[rowLen];
    std::vector<int> sortedInd;
    std::vector<double> sortedVal;

    int tmpLen = rowLen;
    matrix.getRow(row, tmpLen, tempVal, tempInd);
    if (tmpLen != rowLen) ERReturn(-1);

    for(int j=0; j<tmpLen; j++) {
      int offset = fei::sortedListInsert(tempInd[j], sortedInd);
      if (offset < 0) ERReturn(-1);
      sortedVal.insert(sortedVal.begin()+offset, tempVal[j]);
    }

    delete [] tempInd;
    delete [] tempVal;

    std::vector<int>& colInds = colIndices[i];
    if (sortedInd != colInds) {
      ERReturn(-1);
    }
    std::vector<double>& vals = values[i];
    if (sortedVal != vals) {
      ERReturn(-1);
    }
  }

  return(0);
}

//==============================================================================
int fill_DMSR(fei_trilinos::AztecDMSR_Matrix& matrix, int localOffset,
	      std::vector<std::vector<int> >& colIndices,
	      std::vector<std::vector<double> >& values, bool sumInto)
{
  int localSize = colIndices.size();

  for(int i=0; i<localSize; i++) {
    int row = localOffset + i;
    int rowLen = values[i].size();

    if (sumInto) {
      CHK_ERR( matrix.sumIntoRow(row, rowLen,
				 &(values[i][0]), &(colIndices[i][0])) );
    }
    else {
      CHK_ERR( matrix.putRow(row, rowLen,
			     &(values[i][0]), &(colIndices[i][0])) );
    }
  }

  return(0);
}
#endif
