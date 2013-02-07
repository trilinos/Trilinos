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
 *  \todo test with GraphAdapter: add testGraphAdapter which is 
           like testMatrixAdapter except is uses GraphAdapter
           queries and it may have edges weights.
 *  \todo Address the TODOs in the code below.
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

typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;
typedef Tpetra::Map<lno_t, gno_t, node_t> tmap_t;
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

void testMatrixAdapter(RCP<const Tpetra::CrsMatrix<scalar_t, lno_t, gno_t> > &M,
    const RCP<const Comm<int> > &comm,
    bool idsAreConsecutive,
    int rowWeightDim, int nnzDim, int coordDim,
    bool consecutiveIdsRequested, bool removeSelfEdges)
{
  int fail=0;
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  lno_t nLocalRows = M->getNodeNumRows();
  lno_t nLocalNZ = M->getNodeNumEntries();
  gno_t nGlobalRows =  M->getGlobalNumRows();
  gno_t nGlobalNZ = M->getGlobalNumEntries();

  // Create a matrix input adapter.

  typedef Zoltan2::MatrixAdapter<tcrsMatrix_t> base_adapter_t;
  typedef Zoltan2::XpetraCrsMatrixAdapter<tcrsMatrix_t> adapter_t;

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags;
  if (consecutiveIdsRequested)
    modelFlags.set(Zoltan2::IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  if (removeSelfEdges)
    modelFlags.set(Zoltan2::SELF_EDGES_MUST_BE_REMOVED);

  scalar_t **coords=NULL;
  scalar_t **rowWeights=NULL;

  if (coordDim > 0){
    coords = new scalar_t * [coordDim];
    for (int i=0; i < coordDim; i++){
      coords[i] = new scalar_t [nLocalRows];
      for (lno_t j=0; j < nLocalRows; j++){
        coords[i][j] = 100000*i + j;
      }
    }
  }

  if (rowWeightDim > 0){
    rowWeights = new scalar_t * [rowWeightDim];
    for (int i=0; i < rowWeightDim; i++){
      if (nnzDim == i)
        rowWeights[i] = NULL;
      else{
        rowWeights[i] = new scalar_t [nLocalRows];
        for (lno_t j=0; j < nLocalRows; j++){
          rowWeights[i][j] = 200000*i + j;
        }
      }
    }
  }

  adapter_t tmi(M, coordDim, rowWeightDim);

  for (int i=0; i < rowWeightDim; i++){
    if (nnzDim == i)
      tmi.setRowWeightIsNumberOfNonZeros(i);
    else
      tmi.setRowWeights(i, rowWeights[i], 1);
  }

  for (int i=0; i < coordDim; i++){
    tmi.setRowCoordinates(i, coords[i], 1);
  }

  int numLocalDiags = M->getNodeNumDiags();
  int numGlobalDiags = M->getGlobalNumDiags();

  const RCP<const tmap_t> rowMap = M->getRowMap();
  const RCP<const tmap_t> colMap = M->getColMap();

  // How many neighbor vertices are not on my process.

  // we know GIDs are consecutive
  int minLocalGID = rowMap->getMinGlobalIndex();
  int maxLocalGID = rowMap->getMaxGlobalIndex();

  int *numNbors = new int [nLocalRows];
  int *numLocalNbors = new int [nLocalRows];
  bool *haveDiag = new bool [nLocalRows];
  gno_t totalLocalNbors = 0;

  for (lno_t i=0; i < nLocalRows; i++){
    numLocalNbors[i] = 0;
    haveDiag[i] = false;
    ArrayView<const lno_t> idx;
    ArrayView<const scalar_t> val;
    M->getLocalRowView(i, idx, val);
    numNbors[i] = idx.size();

    for (lno_t j=0; j < idx.size(); j++){
      if (idx[j] == i){
        haveDiag[i] = true;
      }
      else{
        gno_t gidVal = colMap->getGlobalElement(idx[j]);
        if (gidVal >= minLocalGID && gidVal <= maxLocalGID){
          numLocalNbors[i]++;
          totalLocalNbors++;
        }
      }
    }
  }

  // Create a GraphModel based on this input data.

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

  if (model->getLocalNumVertices() != size_t(nLocalRows))
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumVertices", 1)

  if (model->getGlobalNumVertices() != size_t(nGlobalRows))
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getGlobalNumVertices", 1)

  size_t num = (removeSelfEdges ? (nLocalNZ-numLocalDiags) : nLocalNZ);

  if (model->getLocalNumGlobalEdges() != num)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumGlobalEdges", 1)

  num = totalLocalNbors;
  if (!removeSelfEdges){
    for (lno_t i=0; i < nLocalRows; i++)
      if (haveDiag[i])
        num++;
  }

  if (model->getLocalNumLocalEdges() != num)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumLocalEdges", 1)

  num = (removeSelfEdges ? (nGlobalNZ-numGlobalDiags) : nGlobalNZ);

  if (model->getGlobalNumEdges() != num)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getGlobalNumEdges", 1)

  if (model->getVertexWeightDim() != rowWeightDim)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexWeightDim", 1)

  if (model->getEdgeWeightDim() != 0)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getEdgeWeightDim", 1)

  if (model->getCoordinateDim() != coordDim)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getCoordinateDim", 1)

  ArrayView<const gno_t> vertexGids;
  ArrayView<input_t> crds;
  ArrayView<input_t> wgts;

  try{
    model->getVertexList(vertexGids, crds, wgts);
  }
  catch (std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexList", 1)

  if (vertexGids.size() != nLocalRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexList", 1)

  // We know model stores things in same order we gave it.

  if (idsAreConsecutive){
    for (lno_t i=0; i < nLocalRows; i++){
      if (vertexGids[i] != minLocalGID + i) {
        fail = 1;
        break;
      }
    }
  }
  else{  // round robin ids
    gno_t myGid = rank;
    for (lno_t i=0; i < nLocalRows; i++, myGid += nprocs){
      if (vertexGids[i] != myGid){
        fail = 1;
        break;
      }
    }
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "vertex gids", 1)

  if ((crds.size() != coordDim) || (wgts.size() != rowWeightDim))
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "coord or weight array size", 1)

  for (int i=0; !fail && i < coordDim; i++){
    for (lno_t j=0; j < nLocalRows; j++){
      if (crds[i][j] != 200000*i + j){
        fail = 1;
        break;
      }
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "coord values", 1)

  for (int i=0; !fail && i < rowWeightDim; i++){
    if (nnzDim == i){
      for (lno_t j=0; j < nLocalRows; j++){
        scalar_t val = numNbors[j];
        if (removeSelfEdges && haveDiag[j])
          val -= 1;
        if (wgts[i][j] != val){
          fail = 1;
          break;
        }
      }
    }
    else{
      for (lno_t j=0; j < nLocalRows; j++){
        if (wgts[i][j] != 100000*i + j){
          fail = 1;
          break;
        }
      }
    }
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "row weight values", 1)
  
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

  TEST_FAIL_AND_EXIT(*comm, wgts.size() == 0, "edge weights present", 1)

  num = 0;
  for (ArrayView<const lno_t>::size_type i=0; i < offsets.size()-1; i++){
    size_t edgeListSize = offsets[i+1] - offsets[i];
    num += edgeListSize;
    size_t val = numNbors[i];
    if (removeSelfEdges && haveDiag[i])
      val--;
    if (edgeListSize != val){
      fail = 1;
      break;
    }
  }

  TEST_FAIL_AND_EXIT(*comm, numEdges==num, "getEdgeList size", 1)

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

#ifdef THIS_CODE_DOESNT_FAIL
  num = 0;
  for (size_t i=0; i < localOffsets.size()-1; i++){
    size_t edgeListSize = localOffsets[i+1] - localOffsets[i];
    num += edgeListSize;
    size_t val = numLocalNbors[i];
    if (removeSelfEdges && haveDiag[i])
      val--;
    if (edgeListSize != val){
      fail = 1;
      break;
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalEdgeList list sizes", 1)

  TEST_FAIL_AND_EXIT(*comm, numLocalNeighbors==num, "getLocalEdgeList size", 1)

  if (size_t(totalLocalNbors) != numLocalNeighbors)
    fail = 1;

  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalEdgeList size", 1)

  if (nGlobalRows < 200){
    if (totalLocalNbors == 0){
      if (rank == 0)
        std::cout << "  Graph of local edges is empty" << std::endl; 
    }
    else{
      printGraph(nLocalRows, vertexGids.getRawPtr(), 
        localEdges.getRawPtr(), NULL, NULL, localOffsets.getRawPtr(), comm);
    }
  }
#endif

  delete model;

  if (nLocalRows){
    delete [] numNbors;
    delete [] numLocalNbors;
    delete [] haveDiag;

    if (rowWeightDim > 0){
      for (int i=0; i < rowWeightDim; i++){
        if (rowWeights[i])
          delete [] rowWeights[i];
      }
      delete [] rowWeights;
    }

    if (coordDim > 0){
      for (int i=0; i < coordDim; i++){
        if (coords[i])
          delete [] coords[i];
      }
      delete [] coords;
    }
  }

  if (rank==0) std::cout << "    OK" << std::endl;
}

void testGraphModel(string fname, gno_t xdim, gno_t ydim, gno_t zdim,
    const RCP<const Comm<int> > &comm,
    int rowWeightDim, int nnzDim, int coordDim,
    bool consecutiveIdsRequested, bool removeSelfEdges)
{
  int rank = comm->getRank();

  if (rank==0){
    cout << endl << "=======================" << endl;
    if (fname.size() > 0)
      cout << endl << "Test parameters: file name " << fname << endl;
    else{
      cout << endl << "Test parameters: dimension ";
      cout  << xdim << "x" << ydim << "x" << zdim << endl;
    }

    cout << "Vertex weight dim: " << rowWeightDim << endl;
    if (nnzDim >= 0)
     cout << "  Dimension " << nnzDim << " is number of neighbors" << endl;

    cout << "Coordinate dim: " << coordDim << endl;
    cout << "Request consecutive vertex gids: ";
    cout << (consecutiveIdsRequested ? "yes" : "no") << endl;
    cout << "Request to remove self edges: ";
    cout << (removeSelfEdges ? "yes" : "no") << endl;
  }

  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;

  // Input generator
  UserInputForTests *input;

  if (fname.size() > 0)
    input = new UserInputForTests(testDataFilePath, fname, comm, true);
  else
    input = new UserInputForTests(xdim,ydim,zdim,string(""), comm, true);

  RCP<tcrsMatrix_t> M = input->getTpetraCrsMatrix();

  // Row Ids of test input are already consecutive

  RCP<const tcrsMatrix_t> Mconsec = rcp_const_cast<const tcrsMatrix_t>(M);

  RCP<const Tpetra::CrsGraph<lno_t, gno_t> > graph = Mconsec->getCrsGraph();

  printTpetraGraph<lno_t, gno_t>(comm, *graph, cout, 100, 
    "Graph with consecutive IDs");

  bool idsAreConsecutive = true;

  testMatrixAdapter(Mconsec, comm,  idsAreConsecutive,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

#if 0
  testGraphAdapter(Mconsec, comm, idsAreConsecutive,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);
#endif

  // Do a round robin migration so that global IDs are not consecutive.

#ifdef TODO_THESE_HAVE_BEEN_TESTED
  Array<gno_t> myNewRows;
  for (size_t i=rank; i < Mconsec->getGlobalNumRows(); i+=nprocs)
    myNewRows.push_back(i);

  RCP<const tcrsMatrix_t> Mnonconsec = 
    Zoltan2::XpetraTraits<tcrsMatrix_t>::doMigration(
      Mconsec, myNewRows.size(), myNewRows.getRawPtr());

  graph = Mnonconsec->getCrsGraph();

  printTpetraGraph<lno_t, gno_t>(comm, *graph, cout, 100, 
    "Graph with non-consecutive IDs");

  idsAreConsecutive = false;

  testMatrixAdapter(Mnonconsec, comm, idsAreConsecutive,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

#if 0
  testGraphAdapter(Mnonconsec, comm, idsAreConsecutive,
     rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);
#endif

#endif

  delete input;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  int rowWeightDim=0;
  int nnzDim = -1; 
  int coordDim=0;
  bool consecutiveIdsRequested=false, removeSelfEdges=false;
  string fname("simple");

  testGraphModel(fname, 0, 0, 0, comm,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

#ifdef TODO_THESE_HAVE_BEEN_TESTED
  rowWeightDim = 1;

  testGraphModel(fname, 0, 0, 0, comm,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  nnzDim = 1;

  testGraphModel(fname, 0, 0, 0, comm,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  rowWeightDim = 2;
  coordDim = 3;

  testGraphModel(fname, 0, 0, 0, comm,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  consecutiveIdsRequested = true;

  testGraphModel(fname, 0, 0, 0, comm,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  removeSelfEdges = true;

  testGraphModel(fname, 0, 0, 0, comm,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  consecutiveIdsRequested = false;

  testGraphModel(fname, 0, 0, 0, comm,
    rowWeightDim, nnzDim, coordDim,
    consecutiveIdsRequested, removeSelfEdges);
#endif

  if (rank==0)
    cout << "PASS" << endl;

  return 0;
}

