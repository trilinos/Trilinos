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

#include <fei_macros.hpp>
#include <fei_mpi.h>

#include <test_utils/test_SNL_FEI_Structure.hpp>

#include <SNL_FEI_Structure.hpp>

#include <test_utils/testData.hpp>

#undef fei_file
#define fei_file "test_SNL_FEI_Structure.cpp"
#include <fei_ErrMacros.hpp>

test_SNL_FEI_Structure::test_SNL_FEI_Structure(MPI_Comm comm)
 : tester(comm)
{
}

test_SNL_FEI_Structure::~test_SNL_FEI_Structure()
{
}

int test_SNL_FEI_Structure::runtests()
{
  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_SNL_FEI_Structure::test1()
{
  testData* testdata = new testData(localProc_, numProcs_);

  SNL_FEI_Structure structure(comm_);

  CHK_ERR( structure.initFields(testdata->fieldIDs.size(),
				&(testdata->fieldSizes[0]),
				&(testdata->fieldIDs[0])) );

  int numNodesPerElem = testdata->ids.size();
  std::vector<int> numFieldsPerNode(numNodesPerElem, 1);
  std::vector<int*>nodalFieldIDs(numNodesPerElem, &(testdata->fieldIDs[0]));

  CHK_ERR( structure.initElemBlock(0, //blockID
				   1, //numElements
				   numNodesPerElem,
				   &numFieldsPerNode[0],
				   &nodalFieldIDs[0],
				   0, //numElemDofFieldsPerElement
				   NULL, //elemDofFieldIDs
				   0)); //interleaveStrategy

  CHK_ERR( structure.initElem(0, //blockID
			      0, //elemID
			      &(testdata->ids[0])) );

  std::vector<int*> sharingProcs2D(testdata->sharedIDs.size());
  int i, offset = 0;
  for(i=0; i<(int)testdata->numSharingProcsPerID.size(); ++i) {
    sharingProcs2D[i] = &(testdata->sharingProcs[offset]);
    offset += testdata->numSharingProcsPerID[i];
  }

  if (testdata->sharedIDs.size() > 0) {
    CHK_ERR( structure.initSharedNodes(testdata->sharedIDs.size(),
      testdata->sharedIDs.size() ? &(testdata->sharedIDs[0]) : 0,
      testdata->numSharingProcsPerID.size() ? &(testdata->numSharingProcsPerID[0]) : 0,
      &sharingProcs2D[0]) );
  }

  CHK_ERR( structure.initComplete() );

  int numActiveNodes = structure.getNumActiveNodes();
  if (numActiveNodes != (int)testdata->ids.size()) {
    ERReturn(-1);
  }

  int fieldSize = structure.getFieldSize(testdata->fieldIDs[0]);

  int numLocalEqns = fieldSize*2;
  if (localProc_ == 0) numLocalEqns += 2;
  int checkNumLocalEqns = structure.getNumLocalEqns();
  if (numLocalEqns != checkNumLocalEqns) {
    ERReturn(-1);
  }

  int numGlobalEqns = fieldSize*(numProcs_*2 + 2);
  int checkNumGlobalEqns = structure.getNumGlobalEqns();
  if (checkNumGlobalEqns != numGlobalEqns) {
    ERReturn(-1);
  }

  std::vector<int> rowLengths;
  CHK_ERR( structure.getMatrixRowLengths(rowLengths) );

  int numNonzeros = 0;
  for(size_t j=0; j<rowLengths.size(); ++j) {
    numNonzeros += rowLengths[j];
  }

  std::vector<int> colIndices_1d(numNonzeros);
  std::vector<int*> colIndPtrs(rowLengths.size());

  offset = 0;
  for(size_t j=0; j<rowLengths.size(); ++j) {
    colIndPtrs[j] = &(colIndices_1d[offset]);
    offset += rowLengths[j];
  }

  CHK_ERR( structure.getMatrixStructure(&colIndPtrs[0],
					rowLengths) );

  delete testdata;

  return(0);
}

int test_SNL_FEI_Structure::test2()
{
  testData* testdata = new testData(localProc_, numProcs_);

  SNL_FEI_Structure structure(comm_);

  CHK_ERR( structure.initFields(testdata->fieldIDs.size(),
				&(testdata->fieldSizes[0]),
				&(testdata->fieldIDs[0])) );

  int numNodesPerElem = testdata->ids.size();
  std::vector<int> numFieldsPerNode(numNodesPerElem, testdata->fieldIDs.size());
  std::vector<int*>nodalFieldIDs(numNodesPerElem, &(testdata->fieldIDs[0]));
  std::vector<int> elemDofFieldIDs = testdata->fieldIDs;

  CHK_ERR( structure.initElemBlock(0, //blockID
				   1, //numElements
				   numNodesPerElem,
				   &numFieldsPerNode[0],
				   &nodalFieldIDs[0],
				   elemDofFieldIDs.size(),
				   &elemDofFieldIDs[0],
				   0)); //interleaveStrategy

  CHK_ERR( structure.initElem(0, //blockID
			      0, //elemID
			      &(testdata->ids[0])) );

  std::vector<int*> sharingProcs2D(testdata->sharedIDs.size());
  int i, offset = 0;
  for(i=0; i<(int)testdata->numSharingProcsPerID.size(); ++i) {
    sharingProcs2D[i] = &(testdata->sharingProcs[offset]);
    offset += testdata->numSharingProcsPerID[i];
  }

  if (testdata->sharedIDs.size() > 0) {
    CHK_ERR( structure.initSharedNodes(testdata->sharedIDs.size(),
      testdata->sharedIDs.size() ? &(testdata->sharedIDs[0]) : 0,
      testdata->numSharingProcsPerID.size() ? &(testdata->numSharingProcsPerID[0]) : 0,
      &sharingProcs2D[0]) );
  }

  CHK_ERR( structure.initComplete() );

  int numActiveNodes = structure.getNumActiveNodes();
  if (numActiveNodes != (int)testdata->ids.size()) {
    ERReturn(-1);
  }

  int numEqnsPerNode = 0;
  for(i=0; i<(int)testdata->fieldSizes.size(); ++i) {
    numEqnsPerNode += testdata->fieldSizes[i];
  }

  int numLocalEqns = 3*numEqnsPerNode;//2 nodes + elem-dofs

  if (localProc_ == 0) numLocalEqns += 2*numEqnsPerNode;
  int checkNumLocalEqns = structure.getNumLocalEqns();
  if (numLocalEqns != checkNumLocalEqns) {
    ERReturn(-1);
  }

  int numGlobalEqns = (numProcs_*3+2)*numEqnsPerNode;
  int checkNumGlobalEqns = structure.getNumGlobalEqns();
  if (checkNumGlobalEqns != numGlobalEqns) {
    ERReturn(-1);
  }

  std::vector<int> rowLengths;
  CHK_ERR( structure.getMatrixRowLengths(rowLengths) );

  int numNonzeros = 0;
  for(size_t j=0; j<rowLengths.size(); ++j) {
    numNonzeros += rowLengths[j];
  }

  std::vector<int> colIndices_1d(numNonzeros);
  std::vector<int*> colIndPtrs(rowLengths.size());

  offset = 0;
  for(size_t j=0; j<rowLengths.size(); ++j) {
    colIndPtrs[j] = &(colIndices_1d[offset]);
    offset += rowLengths[j];
  }

  CHK_ERR( structure.getMatrixStructure(&colIndPtrs[0],
					rowLengths) );

  delete testdata;

  return(0);
}

int test_SNL_FEI_Structure::test3()
{
  return(0);
}

int test_SNL_FEI_Structure::test4()
{
  return(0);
}
