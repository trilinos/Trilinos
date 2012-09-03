// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Testing of GraphModel built from Xpetra matrix input adapters.
//

/*! \brief Test of GraphModel interface.
 *
 *  \todo test all methods of GraphModel
 *  \todo test with GraphInput
 */

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <bitset>

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayView.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::ArrayView;

typedef Zoltan2::StridedData<lno_t, scalar_t> input_t;

using std::string;
using std::vector;

void printGraph(lno_t nrows, const gno_t *v, const lno_t *elid, 
    const gno_t *egid,
    const int *owner, const lno_t *idx, const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  if (rank == 0){
    if (owner)
      std::cout << "Global graph:" << std::endl;
    else
      std::cout << "Local graph:" << std::endl;
  }

  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << "Rank " << p << std::endl;
      if (owner){
        for (lno_t i=0; i < nrows; i++){
          std::cout << "  Vtx " << *v++ << ": ";
          if (elid)
            for (lno_t j=idx[i]; j < idx[i+1]; j++)
              std::cout << *elid++ << " (" << *owner++ << ") ";
          else
            for (lno_t j=idx[i]; j < idx[i+1]; j++)
              std::cout << *egid++ << " (" << *owner++ << ") ";
          std::cout << std::endl;
        }
        std::cout.flush();
      }
      else{
        for (lno_t i=0; i < nrows; i++){
          std::cout << "  Vtx " << i << ": ";
          if (elid)
            for (lno_t j=idx[i]; j < idx[i+1]; j++)
              std::cout << *elid++ << " ";
          else
            for (lno_t j=idx[i]; j < idx[i+1]; j++)
              std::cout << *egid++ << " ";
          std::cout << std::endl;
        }
        std::cout.flush();
      }
    }
    comm->barrier();
  }
  comm->barrier();
}

void checkModel(RCP<const Tpetra::CrsMatrix<scalar_t, lno_t, gno_t> > &M,
    const RCP<const Comm<int> > &comm,
    bool consecutiveIdsRequested, bool noSelfEdges)
{
  int fail=0;
  int rank = comm->getRank();
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags;
  if (consecutiveIdsRequested)
    modelFlags.set(Zoltan2::IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  if (noSelfEdges)
    modelFlags.set(Zoltan2::SELF_EDGES_MUST_BE_REMOVED);

  if (rank==0){
    std::cout << "     Request consecutive IDs " << consecutiveIdsRequested;
    std::cout << ", remove self edges " << noSelfEdges << std::endl;;
  }

  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;
  typedef Zoltan2::MatrixInput<tcrsMatrix_t> base_adapter_t;
  typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> adapter_t;

  adapter_t tmi(M);

  Zoltan2::GraphModel<base_adapter_t> *model = NULL;
  const base_adapter_t *baseTmi = &tmi;

  try{
    model = new Zoltan2::GraphModel<base_adapter_t>(baseTmi, env, 
      comm, modelFlags);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "Creating xpetra graph model", 1)

  // Test the GraphModel interface

  lno_t nLocalRows = M->getNodeNumRows();
  lno_t nLocalNonZeros = M->getNodeNumEntries();
  gno_t nGlobalRows =  M->getGlobalNumRows();
  gno_t nGlobalNonZeros = M->getGlobalNumEntries();

  if (model->getLocalNumVertices() != size_t(nLocalRows))
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumVertices", 1)

  if (model->getGlobalNumVertices() != size_t(nGlobalRows))
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getGlobalNumVertices", 1)

  if (noSelfEdges){
    if (model->getGlobalNumEdges() >  size_t(nGlobalNonZeros))
      fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "getGlobalNumEdges", 1)

    if (model->getLocalNumGlobalEdges() > size_t(nLocalNonZeros))
      fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumGlobalEdges", 1)
  }
  else{
    if (model->getGlobalNumEdges() !=  size_t(nGlobalNonZeros))
      fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "getGlobalNumEdges", 1)

    if (model->getLocalNumGlobalEdges() != size_t(nLocalNonZeros))
      fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumGlobalEdges", 1)
  }

  ArrayView<const gno_t> vertexGids;
  ArrayView<input_t> coords;  // not implemented yet
  ArrayView<input_t> wgts;    // not implemented yet

  try{
    model->getVertexList(vertexGids, coords, wgts);
  }
  catch (std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexList", 1)

  if (vertexGids.size() != nLocalRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexList", 1)
  
  ArrayView<const gno_t> edgeGids;
  ArrayView<const int> procIds;
  ArrayView<const lno_t> offsets;
  size_t numEdges=0;

  try{
    numEdges = model->getEdgeList(edgeGids, procIds, offsets, wgts);
  }
  catch(std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getEdgeList", 1)

  lno_t numLocalEdges = 0;
  for (size_t i=0; i < numEdges; i++){
    if (procIds[i] == rank)
      numLocalEdges++;
  }

  if (numEdges != model->getLocalNumGlobalEdges())
    fail = 1;

  TEST_FAIL_AND_EXIT(*comm, !fail, "getEdgeList size", 1)

  if (nGlobalRows < 200){
    printGraph(nLocalRows, vertexGids.getRawPtr(), NULL,
      edgeGids.getRawPtr(), procIds.getRawPtr(), offsets.getRawPtr(), comm);
  }
  else{
    if (rank==0) 
      std::cout << "    " << nGlobalRows << " total rows" << std::endl;
  }

  // Get graph restricted to this process

  ArrayView<const lno_t> localEdges;
  ArrayView<const lno_t> localOffsets;
  size_t numLocalNeighbors=0;

  try{
    numLocalNeighbors= model->getLocalEdgeList(localEdges, localOffsets, wgts);
  }
  catch(std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalEdgeList", 1)

  if (size_t(numLocalEdges) != numLocalNeighbors)
    fail = 1;

  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalEdgeList size", 1)

  if (nGlobalRows < 200){
    if (numLocalEdges == 0){
      if (rank == 0)
        std::cout << "  Graph of local edges is empty" << std::endl; 
    }
    else{
      printGraph(nLocalRows, vertexGids.getRawPtr(), 
        localEdges.getRawPtr(), NULL, NULL, localOffsets.getRawPtr(), comm);
    }
  }

  // Get graph restricted to this process

  delete model;

  if (rank==0) std::cout << "    OK" << std::endl;
}

void testGraphModel(string fname, gno_t xdim, gno_t ydim, gno_t zdim,
    const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;

  // Input generator
  UserInputForTests *input;

  if (fname.size() > 0)
    input = new UserInputForTests(testDataFilePath, fname, comm, true);
  else
    input = new UserInputForTests(xdim,ydim,zdim,string(""), comm, true);

  RCP<tcrsMatrix_t> M = input->getTpetraCrsMatrix();

  bool consecutiveIds=true;
  bool noSelfEdges=true;

  // Row Ids of test input are already consecutive

  RCP<const tcrsMatrix_t> Mconsec = rcp_const_cast<const tcrsMatrix_t>(M);

  RCP<const Tpetra::CrsGraph<lno_t, gno_t> > graph = Mconsec->getCrsGraph();

  printTpetraGraph<lno_t, gno_t>(comm, *graph, cout, 100, 
    "Graph with consecutive IDs");

//  checkModel(Mconsec, comm, consecutiveIds, noSelfEdges);
//  checkModel(Mconsec, comm, !consecutiveIds, !noSelfEdges);
//  checkModel(Mconsec, comm, !consecutiveIds, noSelfEdges);

  // Do a round robin migration so that global IDs are not consecutive.

  Array<gno_t> myNewRows;
  for (size_t i=rank; i < Mconsec->getGlobalNumRows(); i+=nprocs)
    myNewRows.push_back(i);

  RCP<const tcrsMatrix_t> Mnonconsec = 
    Zoltan2::XpetraTraits<tcrsMatrix_t>::doMigration(
      Mconsec, myNewRows.size(), myNewRows.getRawPtr());

  graph = Mnonconsec->getCrsGraph();

  printTpetraGraph<lno_t, gno_t>(comm, *graph, cout, 100, 
    "Graph with non-consecutive IDs");

  checkModel(Mnonconsec, comm, consecutiveIds, noSelfEdges);
  checkModel(Mnonconsec, comm, !consecutiveIds, !noSelfEdges);
  checkModel(Mnonconsec, comm, !consecutiveIds, noSelfEdges);

  delete input;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  string fname("simple");

  if (rank==0)
    std::cout << fname << std::endl;

  testGraphModel(fname, 0, 0, 0, comm);

  if (rank==0)
    std::cout << "PASS" << std::endl;

  return 0;
}

