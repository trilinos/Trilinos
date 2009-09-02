/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_macros.hpp>

#include <test_utils/test_Algebraic.hpp>

#include <fei_VectorSpace.hpp>
#include <fei_MatrixGraph_Impl2.hpp>
#include <fei_SparseRowGraph.hpp>

#undef fei_file
#define fei_file "test_Algebraic.cpp"
#include <fei_ErrMacros.hpp>

test_Algebraic::test_Algebraic(MPI_Comm comm)
 : tester(comm)
{
}

test_Algebraic::~test_Algebraic()
{
}

int test_Algebraic::runtests()
{
  if (numProcs_ < 2) {
    CHK_ERR( serialtest1() );
    CHK_ERR( serialtest2() );
  }

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_Algebraic::serialtest1()
{
  int i, numRows = 10;
  fei::SharedPtr<fei::VectorSpace> vspace(new fei::VectorSpace(comm_));

  int idType = 0;

  vspace->defineIDTypes(1, &idType);

  std::vector<int> rowNumbers(numRows);
  for(i=0; i<numRows; ++i) {
    rowNumbers[i] = i;
  }

  CHK_ERR( vspace->addDOFs(idType, numRows, &rowNumbers[0]) );

  CHK_ERR( vspace->initComplete() );

  int index = -1;
  CHK_ERR( vspace->getGlobalIndex(idType, rowNumbers[3], index) );

  if (index != 3) {
    ERReturn(-1);
  }

  int numDOF = vspace->getNumDegreesOfFreedom(idType, rowNumbers[3]);

  if (numDOF != 1) {
    ERReturn(-1);
  }

  std::vector<int> globalOffsets;

  vspace->getGlobalIndexOffsets(globalOffsets);

  if (globalOffsets[0] != 0) {
    ERReturn(-1);
  }

  if (globalOffsets[1] != numRows) {
    ERReturn(-1);
  }

  fei::SharedPtr<fei::VectorSpace> dummy;
  fei::MatrixGraph_Impl2 mgraph(vspace, dummy);

  std::vector<int> rowOffsets(numRows+1);
  std::vector<int> packedColumnIDs(numRows);
  for(i=0; i<numRows; ++i) {
    rowOffsets[i] = i;
    packedColumnIDs[i] = i;
  }
  rowOffsets[numRows] = numRows;

  CHK_ERR( mgraph.initConnectivity(idType, numRows,
				   &rowNumbers[0],
				   &rowOffsets[0],
				   &packedColumnIDs[0]) );

  CHK_ERR( mgraph.initComplete() );

  fei::SharedPtr<fei::SparseRowGraph> localgraph = mgraph.createGraph(false);

  int mnumRows = localgraph->rowNumbers.size();
  int* mrowOffsets = &(localgraph->rowOffsets[0]);
  int mnumNonzeros = localgraph->packedColumnIndices.size();
  int* mpackedColumnIndices = &(localgraph->packedColumnIndices[0]);

  if (mnumRows != numRows) {
    ERReturn(-1);
  }

  if (mnumNonzeros != numRows) {
    ERReturn(-1);
  }

  for(i=0; i<numRows; ++i) {
    if ((mrowOffsets[i+1]-mrowOffsets[i]) != 1) {
      ERReturn(-1);
    }
    if (mpackedColumnIndices[i] != packedColumnIDs[i]) {
      ERReturn(-1);
    }
  }

  return(0);
}

int test_Algebraic::serialtest2()
{
  int i, numRows = 10;
  fei::SharedPtr<fei::VectorSpace> vspace(new fei::VectorSpace(comm_));

  int idType = 0;

  vspace->defineIDTypes(1, &idType);

  std::vector<int> rowNumbers(numRows);
  for(i=0; i<numRows; ++i) {
    rowNumbers[i] = i;
  }

  CHK_ERR( vspace->addDOFs(idType, numRows, &rowNumbers[0]) );

  CHK_ERR( vspace->initComplete() );

  int index = -1;
  CHK_ERR( vspace->getGlobalIndex(idType, rowNumbers[3], index) );

  if (index != 3) {
    ERReturn(-1);
  }

  int numDOF = vspace->getNumDegreesOfFreedom(idType, rowNumbers[3]);

  if (numDOF != 1) {
    ERReturn(-1);
  }

  std::vector<int> globalOffsets;

  vspace->getGlobalIndexOffsets(globalOffsets);

  if (globalOffsets[0] != 0) {
    ERReturn(-1);
  }

  if (globalOffsets[1] != numRows) {
    ERReturn(-1);
  }

  fei::SharedPtr<fei::VectorSpace> dummy;
  fei::MatrixGraph_Impl2 mgraph(vspace, dummy);

  std::vector<int> rowLengths(numRows);
  std::vector<int> packedColumnIDs(numRows);
  std::vector<int*> columnIDs(numRows);

  for(i=0; i<numRows; ++i) {
    rowLengths[i] = 1;
    packedColumnIDs[i] = i;
    columnIDs[i] = &(packedColumnIDs[i]);
  }

  CHK_ERR( mgraph.initConnectivity(idType, numRows,
				   &rowNumbers[0],
				   &rowLengths[0],
				   &columnIDs[0]) );

  CHK_ERR( mgraph.initComplete() );

  fei::SharedPtr<fei::SparseRowGraph> localgraph = mgraph.createGraph(false);

  int mnumRows = localgraph->rowNumbers.size();
  int* mrowOffsets = &(localgraph->rowOffsets[0]);
  int mnumNonzeros = localgraph->packedColumnIndices.size();
  int* mpackedColumnIndices = &(localgraph->packedColumnIndices[0]);

  if (mnumRows != numRows) {
    ERReturn(-1);
  }

  if (mnumNonzeros != numRows) {
    ERReturn(-1);
  }

  for(i=0; i<numRows; ++i) {
    if ((mrowOffsets[i+1]-mrowOffsets[i]) != 1) {
      ERReturn(-1);
    }
    if (mpackedColumnIndices[i] != packedColumnIDs[i]) {
      ERReturn(-1);
    }
  }

  return(0);
}

int test_Algebraic::test1()
{
  return(0);
}

int test_Algebraic::test2()
{
  return(0);
}

int test_Algebraic::test3()
{
  return(0);
}

int test_Algebraic::test4()
{
  return(0);
}
