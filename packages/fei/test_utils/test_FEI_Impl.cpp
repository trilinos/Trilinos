/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_macros.hpp>
#include <fei_sstream.hpp>

#include <test_utils/test_FEI_Impl.hpp>

#include <feiArray.hpp>
#include <fei_FEI_Impl.hpp>
#include <fei_Aztec_LinSysCore.hpp>
#include <test_utils/FEData.hpp>
#include <fei_LibraryWrapper.hpp>

#include <test_utils/testData.hpp>

#undef fei_file
#define fei_file "test_FEI_Impl.cpp"
#include <fei_ErrMacros.hpp>

test_FEI_Impl::test_FEI_Impl(MPI_Comm comm)
 : tester(comm)
{
}

test_FEI_Impl::~test_FEI_Impl()
{
}

int test_FEI_Impl::runtests()
{
  if (numProcs_<2) {
    CHK_ERR( serialtest1() );
  }

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_FEI_Impl::serialtest1()
{
  return(0);
}

int test_FEI_Impl::compareCoefs(int n,
			     const double*const* coefs1,
			     const double*const* coefs2)
{
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      double diff = coefs1[i][j] - coefs2[i][j];
      if (diff > 1.e-20 || diff < -1.e-20) {
	return(-1);
      }
    }
  }

  return(0);
}

int test_FEI_Impl::test1()
{
  if (numProcs_ > 1) {
    return(0);
  }

  testData* testdata = new testData(localProc_, numProcs_);

  fei::SharedPtr<LinearSystemCore> linSys(new Aztec_LinSysCore(comm_));
  fei::SharedPtr<LibraryWrapper> wrapper(new LibraryWrapper(linSys));
  fei::SharedPtr<fei::FEI_Impl> fei(new fei::FEI_Impl(wrapper, comm_, 0));

  std::string param0("name test1");
  FEI_OSTRINGSTREAM osstr;
  osstr << "debugOutput ";
  if (path_.empty()) osstr << ".";
  else osstr << path_;

  std::string param1 = osstr.str();

  int numParams = 2;
  char** params = new char*[numParams];
  params[0] = const_cast<char*>(param0.c_str());
  params[1] = const_cast<char*>(param1.c_str());

  //call the parameters function a couple of times to test the fei's internal
  //method for merging string lists when parameters is called more than once.
  CHK_ERR( fei->parameters(1, &params[0]) );
  CHK_ERR( fei->parameters(1, &params[1]) );
  CHK_ERR( fei->parameters(2, params) );

  delete [] params;

  CHK_ERR( fei->setIDLists(1, &(testdata->ids[0]),
			   1, &(testdata->ids[0])) );

  CHK_ERR( fei->initFields(testdata->fieldIDs.size(),
				&(testdata->fieldSizes[0]),
				&(testdata->fieldIDs[0])) );

  int numNodesPerElem = testdata->ids.size();
  feiArray<int> numFieldsPerNode(numNodesPerElem);
  numFieldsPerNode = 1;
  feiArray<int*>nodalFieldIDs(numNodesPerElem);
  nodalFieldIDs = &(testdata->fieldIDs[0]);

  CHK_ERR( fei->initElemBlock(0, //blockID
				   1, //numElements
				   numNodesPerElem,
				   numFieldsPerNode.dataPtr(),
				   nodalFieldIDs.dataPtr(),
				   0, //numElemDofFieldsPerElement
				   NULL, //elemDofFieldIDs
				   0)); //interleaveStrategy

  CHK_ERR( fei->initElem(0, //blockID
			      0, //elemID
			      &(testdata->ids[0])) );

  feiArray<int*> sharingProcs2D(testdata->sharedIDs.size());
  int offset = 0;
  for(int i=0; i<(int)testdata->numSharingProcsPerID.size(); ++i) {
    sharingProcs2D[i] = &(testdata->sharingProcs[offset]);
    offset += testdata->numSharingProcsPerID[i];
  }

  if (testdata->sharedIDs.size() > 0) {
    CHK_ERR( fei->initSharedNodes(testdata->sharedIDs.size(),
				     &(testdata->sharedIDs[0]),
				     &(testdata->numSharingProcsPerID[0]),
				     sharingProcs2D.dataPtr()) );
  }

  int patternID = 0;
  int numRowIDs = 1;
  int* numFieldsPerRow = new int[1];
  numFieldsPerRow[0] = 1;
  int* rowFieldIDs = new int[1];
  rowFieldIDs[0] = testdata->fieldIDs[0];

  int rowIDType = FEI_NODE;

  CHK_ERR( fei->initCoefAccessPattern(patternID, numRowIDs, numFieldsPerRow,
				      &rowFieldIDs, numRowIDs, numFieldsPerRow,
				      &rowFieldIDs,
				      0) );//interleaveStrategy

  delete [] numFieldsPerRow;
  delete [] rowFieldIDs;

  CHK_ERR( fei->initCoefAccess(patternID,
			       &rowIDType, &(testdata->ids[0]),
			       &rowIDType, &(testdata->ids[0])) );

  CHK_ERR( fei->initComplete() );

  feiArray<double> rhsData(testdata->ids.size());
  rhsData = 1.0;

  double one = 1.0;
  CHK_ERR( fei->setMatScalars(1, &(testdata->ids[0]), &one) );
  CHK_ERR( fei->setRHSScalars(1, &(testdata->ids[0]), &one) );

  CHK_ERR( fei->setCurrentMatrix(testdata->ids[0]) );
  CHK_ERR( fei->setCurrentRHS(testdata->ids[0]) );

  CHK_ERR( fei->putIntoRHS(FEI_NODE, testdata->fieldIDs[0],
			   testdata->ids.size(),
			   &(testdata->ids[0]),
			   rhsData.dataPtr()) );

  double* matrixEntries = new double[16];
  double** matrixEntriesPtr = new double*[4];
  for(int ii=0; ii<4; ++ii) matrixEntriesPtr[ii] = matrixEntries+ii*4;
  for(int jj=0; jj<16; ++jj) matrixEntries[jj] = 1.0;

  
  CHK_ERR( fei->putIntoMatrix(patternID, &rowIDType, &(testdata->ids[0]),
			      &rowIDType, &(testdata->ids[0]),
			      matrixEntriesPtr) );

  CHK_ERR( fei->putIntoRHS(patternID, &rowIDType, &(testdata->ids[0]),
			   matrixEntries) );

  CHK_ERR( fei->sumIntoRHS(patternID, &rowIDType, &(testdata->ids[0]),
			   matrixEntries) );

  CHK_ERR( fei->sumIntoMatrix(patternID, &rowIDType, &(testdata->ids[0]),
			      &rowIDType, &(testdata->ids[0]),
			      matrixEntriesPtr) );

  CHK_ERR( fei->getFromMatrix(patternID, &rowIDType, &(testdata->ids[0]),
			      &rowIDType, &(testdata->ids[0]),
			      matrixEntriesPtr) );

  delete [] matrixEntries;
  delete [] matrixEntriesPtr;

  int numBCNodes = 2;
  GlobalID* BCNodeIDs = &(testdata->ids[0]);
  int BCFieldID = testdata->fieldIDs[0];
  double* values = new double[numBCNodes];
  int* offsetsIntoField = new int[numBCNodes];
  for(int ii=0; ii<numBCNodes; ++ii) {
    values[ii] = 1.0;
    offsetsIntoField[ii] = 0;
  }

  CHK_ERR( fei->loadNodeBCs(numBCNodes, BCNodeIDs, BCFieldID,
			    offsetsIntoField, values) );

  delete [] offsetsIntoField;
  delete [] values;

  CHK_ERR( fei->loadComplete() );

  int numActiveNodes = 0;
  CHK_ERR( fei->getNumLocalNodes(numActiveNodes) );
  if (numActiveNodes != (int)testdata->ids.size()) {
    ERReturn(-1);
  }

  GlobalID* localNodes = new GlobalID[numActiveNodes];
  CHK_ERR( fei->getLocalNodeIDList(numActiveNodes, localNodes, numActiveNodes) );

  int totalFieldSize = 0;
  for(int ii=0; ii<(int)testdata->fieldSizes.size(); ++ii) {
    totalFieldSize += testdata->fieldSizes[ii];
  }

  double* soln = new double[numActiveNodes*totalFieldSize];
  int* offsets = new int[numActiveNodes+1];

  CHK_ERR( fei->getNodalSolution(numActiveNodes, localNodes,
				 offsets, soln) );
  delete [] offsets;
  delete [] soln;
  delete [] localNodes;

  CHK_ERR( fei->resetInitialGuess() );

  int fieldSize = 0;
  CHK_ERR( fei->getFieldSize(testdata->fieldIDs[0], fieldSize) );

  double initTime, loadTime, solveTime, solnReturnTime;
  CHK_ERR( fei->cumulative_cpu_times(initTime, loadTime, solveTime,
				      solnReturnTime) );

  delete testdata;

  return(0);
}

int test_FEI_Impl::test2()
{
  fei::SharedPtr<testData> testdata(new testData(localProc_, numProcs_));
  fei::SharedPtr<FiniteElementData> fedata(new FEData(comm_));
  fei::SharedPtr<LibraryWrapper> wrapper(new LibraryWrapper(fedata));
  fei::SharedPtr<fei::FEI_Impl> fei(new fei::FEI_Impl(wrapper, comm_, 0));

  std::string param0("name test1");
  FEI_OSTRINGSTREAM osstr;
  osstr << "debugOutput ";
  if (path_.empty()) osstr << ".";
  else osstr << path_;

  std::string param1 = osstr.str();

  int numParams = 2;
  char** params = new char*[numParams];
  params[0] = const_cast<char*>(param0.c_str());
  params[1] = const_cast<char*>(param1.c_str());

  //call the parameters function a couple of times to test the fei's internal
  //method for merging string lists when parameters is called more than once.
  CHK_ERR( fei->parameters(1, &params[0]) );
  CHK_ERR( fei->parameters(1, &params[1]) );
  CHK_ERR( fei->parameters(2, params) );

  delete [] params;

  CHK_ERR( fei->setIDLists(1, &(testdata->ids[0]),
			   1, &(testdata->ids[0])) );
  CHK_ERR( fei->initFields(testdata->fieldIDs.size(),
				&(testdata->fieldSizes[0]),
				&(testdata->fieldIDs[0])) );

  int numNodesPerElem = testdata->ids.size();
  feiArray<int> numFieldsPerNode(numNodesPerElem);
  numFieldsPerNode = 1;
  feiArray<int*>nodalFieldIDs(numNodesPerElem);
  nodalFieldIDs = &(testdata->fieldIDs[0]);

  CHK_ERR( fei->initElemBlock(0, //blockID
				   1, //numElements
				   numNodesPerElem,
				   numFieldsPerNode.dataPtr(),
				   nodalFieldIDs.dataPtr(),
				   0, //numElemDofFieldsPerElement
				   NULL, //elemDofFieldIDs
				   0)); //interleaveStrategy

  CHK_ERR( fei->initElem(0, //blockID
			      0, //elemID
			      &(testdata->ids[0])) );

  feiArray<int*> sharingProcs2D(testdata->sharedIDs.size());
  int i, offset = 0;
  for(i=0; i<(int)testdata->numSharingProcsPerID.size(); ++i) {
    sharingProcs2D[i] = &(testdata->sharingProcs[offset]);
    offset += testdata->numSharingProcsPerID[i];
  }

  if (testdata->sharedIDs.size() > 0) {
    CHK_ERR( fei->initSharedNodes(testdata->sharedIDs.size(),
				     &(testdata->sharedIDs[0]),
				     &(testdata->numSharingProcsPerID[0]),
				     sharingProcs2D.dataPtr()) );
  }

  int patternID = 0;
  int numRowIDs = 1;
  int* numFieldsPerRow = new int[1];
  numFieldsPerRow[0] = 1;
  int* rowFieldIDs = new int[1];
  rowFieldIDs[0] = testdata->fieldIDs[0];

  int rowIDType = FEI_NODE;

  CHK_ERR( fei->initCoefAccessPattern(patternID, numRowIDs, numFieldsPerRow,
				      &rowFieldIDs, numRowIDs, numFieldsPerRow,
				      &rowFieldIDs,
				      0) );//interleaveStrategy

  delete [] numFieldsPerRow;
  delete [] rowFieldIDs;

  CHK_ERR( fei->initCoefAccess(patternID,
			       &rowIDType, &(testdata->ids[0]),
			       &rowIDType, &(testdata->ids[0])) );

  CHK_ERR( fei->initComplete() );

  feiArray<double> rhsData(testdata->ids.size());
  rhsData = 1.0;

  double one = 1.0;
  CHK_ERR( fei->setMatScalars(1, &(testdata->ids[0]), &one) );
  CHK_ERR( fei->setRHSScalars(1, &(testdata->ids[0]), &one) );

  CHK_ERR( fei->setCurrentMatrix(testdata->ids[0]) );
  CHK_ERR( fei->setCurrentRHS(testdata->ids[0]) );

  int ii;

  int numBCNodes = 2;
  GlobalID* BCNodeIDs = &(testdata->ids[0]);
  int BCFieldID = testdata->fieldIDs[0];
  double* values = new double[numBCNodes];
  int* offsetsIntoField = new int[numBCNodes];
  for(ii=0; ii<numBCNodes; ++ii) {
    values[ii] = 1.0;
    offsetsIntoField[ii] = 0;
  }

  CHK_ERR( fei->loadNodeBCs(numBCNodes, BCNodeIDs, BCFieldID,
			    offsetsIntoField, values) );

  delete [] values;
  delete [] offsetsIntoField;

  CHK_ERR( fei->loadComplete() );

  int numActiveNodes = 0;
  CHK_ERR( fei->getNumLocalNodes(numActiveNodes) );
  if (numActiveNodes != (int)testdata->ids.size()) {
    ERReturn(-1);
  }

  GlobalID* localNodes = new GlobalID[numActiveNodes];
  CHK_ERR( fei->getLocalNodeIDList(numActiveNodes, localNodes, numActiveNodes) );

  int totalFieldSize = 0;
  for(ii=0; ii<(int)testdata->fieldSizes.size(); ++ii) {
    totalFieldSize += testdata->fieldSizes[ii];
  }

  double* soln = new double[numActiveNodes*totalFieldSize];
  int* offsets = new int[numActiveNodes+1];

  CHK_ERR( fei->getNodalSolution(numActiveNodes, localNodes,
				 offsets, soln) );
  delete [] offsets;
  delete [] soln;
  delete [] localNodes;

  CHK_ERR( fei->resetInitialGuess() );

  int fieldSize = 0;
  CHK_ERR( fei->getFieldSize(testdata->fieldIDs[0], fieldSize) );

  double initTime, loadTime, solveTime, solnReturnTime;
  CHK_ERR( fei->cumulative_cpu_times(initTime, loadTime, solveTime,
				      solnReturnTime) );

  return(0);
}

int test_FEI_Impl::test3()
{
  return(0);
}

int test_FEI_Impl::test4()
{
  return(0);
}
