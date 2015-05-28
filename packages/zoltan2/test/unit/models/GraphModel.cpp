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
           like testAdapter except is uses GraphAdapter
           queries and it may have edges weights.
 *  \todo Address the TODOs in the code below.
 */

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_InputTraits.hpp>

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

typedef Zoltan2::BasicUserTypes<zscalar_t, zzgid_t, zlno_t, zgno_t> simpleUser_t;

typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t>      tcrsMatrix_t;
typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t>                 tcrsGraph_t;
typedef Tpetra::Map<zlno_t, zgno_t, znode_t>                      tmap_t;

typedef Zoltan2::BasicVectorAdapter<simpleUser_t>              simpleVAdapter_t;

typedef Zoltan2::MatrixAdapter<tcrsMatrix_t,simpleUser_t>      baseMAdapter_t;
typedef Zoltan2::GraphAdapter<tcrsGraph_t,simpleUser_t>        baseGAdapter_t;

typedef Zoltan2::XpetraCrsMatrixAdapter<tcrsMatrix_t,simpleUser_t> xMAdapter_t;
typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t,simpleUser_t>   xGAdapter_t;

using std::string;
using std::vector;

/////////////////////////////////////////////////////////////////////////////
void printGraph(zlno_t nrows, const zgno_t *v,
    const zlno_t *elid, const zgno_t *egid,
    const int *owner, const zlno_t *idx,
    const RCP<const Comm<int> > &comm)
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
        for (zlno_t i=0; i < nrows; i++){
          std::cout << "  Vtx " << *v++ << ": ";
          if (elid)
            for (zlno_t j=idx[i]; j < idx[i+1]; j++)
              std::cout << *elid++ << " (" << *owner++ << ") ";
          else
            for (zlno_t j=idx[i]; j < idx[i+1]; j++)
              std::cout << *egid++ << " (" << *owner++ << ") ";
          std::cout << std::endl;
        }
        std::cout.flush();
      }
      else{
        for (zlno_t i=0; i < nrows; i++){
          std::cout << "  Vtx " << i << ": ";
          if (elid)
            for (zlno_t j=idx[i]; j < idx[i+1]; j++)
              std::cout << *elid++ << " ";
          else
            for (zlno_t j=idx[i]; j < idx[i+1]; j++)
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

/////////////////////////////////////////////////////////////////////////////
template <typename BaseAdapter, typename Adapter, typename MatrixOrGraph>
void testAdapter(
    RCP<const MatrixOrGraph> &M,
    RCP<const Tpetra::CrsGraph<zlno_t, zgno_t> > &Mgraph,
    const RCP<const Comm<int> > &comm,
    bool idsAreConsecutive,
    int nVtxWeights, int nEdgeWeights, int nnzWgtIdx, int coordDim,
    bool consecutiveIdsRequested, bool removeSelfEdges)
{
  typedef Zoltan2::StridedData<zlno_t, zscalar_t> input_t;

  int fail=0;
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  zlno_t nLocalRows = M->getNodeNumRows();
  zlno_t nLocalNZ = M->getNodeNumEntries();
  zgno_t nGlobalRows =  M->getGlobalNumRows();
  zgno_t nGlobalNZ = M->getGlobalNumEntries();

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags;
  if (consecutiveIdsRequested)
    modelFlags.set(Zoltan2::IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  if (removeSelfEdges)
    modelFlags.set(Zoltan2::SELF_EDGES_MUST_BE_REMOVED);

  // Set up some fake input
  zscalar_t **coords=NULL;
  zscalar_t **rowWeights=NULL;

  if (coordDim > 0){
    coords = new zscalar_t * [coordDim];
    for (int i=0; i < coordDim; i++){
      coords[i] = new zscalar_t [nLocalRows];
      for (zlno_t j=0; j < nLocalRows; j++){
        coords[i][j] = 100000*i + j;
      }
    }
  }

  if (nVtxWeights > 0){
    rowWeights = new zscalar_t * [nVtxWeights];
    for (int i=0; i < nVtxWeights; i++){
      if (nnzWgtIdx == i)
        rowWeights[i] = NULL;
      else{
        rowWeights[i] = new zscalar_t [nLocalRows];
        for (zlno_t j=0; j < nLocalRows; j++){
          rowWeights[i][j] = 200000*i + j;
        }
      }
    }
  }

  if (nEdgeWeights > 0){
    printf("TODO:  STILL NEED TO TEST EDGE WEIGHTS!\n");
  }

  // Create a matrix or graph input adapter.

  Adapter tmi(M, nVtxWeights);

  for (int i=0; i < nVtxWeights; i++){
    if (nnzWgtIdx == i)
      tmi.setWeightIsDegree(i);
    else
      tmi.setWeights(rowWeights[i], 1, i);
  }

  zzgid_t *gids = NULL;

  simpleVAdapter_t *via = NULL;

  if (coordDim > 0) {
    gids = new zzgid_t[nLocalRows];
    for (zlno_t i = 0; i < nLocalRows; i++)
      gids[i] = M->getRowMap()->getGlobalElement(i);
    via = new simpleVAdapter_t(nLocalRows, gids, coords[0],
                                           (coordDim > 1 ? coords[1] : NULL), 
                                           (coordDim > 2 ? coords[2] : NULL),
                                            1,1,1);
    tmi.setCoordinateInput(via);
  }

  int numLocalDiags = M->getNodeNumDiags();
  int numGlobalDiags = M->getGlobalNumDiags();

  const RCP<const tmap_t> rowMap = M->getRowMap();
  const RCP<const tmap_t> colMap = M->getColMap();

  // How many neighbor vertices are not on my process.

  int *numNbors = new int [nLocalRows];
  int *numLocalNbors = new int [nLocalRows];
  bool *haveDiag = new bool [nLocalRows];
  zgno_t totalLocalNbors = 0;

  for (zlno_t i=0; i < nLocalRows; i++){
    numLocalNbors[i] = 0;
    haveDiag[i] = false;
    ArrayView<const zlno_t> idx;
    Mgraph->getLocalRowView(i, idx);
    numNbors[i] = idx.size();

    for (zlno_t j=0; j < idx.size(); j++){
      if (idx[j] == i){
        haveDiag[i] = true;
        numLocalNbors[i]++;
        totalLocalNbors++;
      }
      else{
        zgno_t gidVal = colMap->getGlobalElement(idx[j]);
        if (rowMap->getLocalElement(gidVal) !=
            Teuchos::OrdinalTraits<zlno_t>::invalid()) {
          // nbor is local to this process
          numLocalNbors[i]++;
          totalLocalNbors++;
        }
      }
    }
  }

  // Create a GraphModel based on this input data.

  if (rank == 0) std::cout << "        Creating GraphModel" << std::endl;
  Zoltan2::GraphModel<BaseAdapter> *model = NULL;
  const BaseAdapter *baseTmi = dynamic_cast<BaseAdapter *>(&tmi);

  try{
    model = new Zoltan2::GraphModel<BaseAdapter>(baseTmi, env, 
                                                 comm, modelFlags);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "Creating graph model", 1)

  // Test the GraphModel interface

  if (rank == 0) std::cout << "        Checking counts" << std::endl;
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

  num = (removeSelfEdges ? (totalLocalNbors - numLocalDiags) : totalLocalNbors);

  if (model->getLocalNumLocalEdges() != num)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumLocalEdges", 1)

  num = (removeSelfEdges ? (nGlobalNZ-numGlobalDiags) : nGlobalNZ);

  if (model->getGlobalNumEdges() != num)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getGlobalNumEdges", 1)

  if (model->getNumWeightsPerVertex() != nVtxWeights)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getNumWeightsPerVertex", 1)

  if (model->getNumWeightsPerEdge() != 0)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getNumWeightsPerEdge", 1)

  if (model->getCoordinateDim() != coordDim)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getCoordinateDim", 1)

  if (rank == 0) std::cout << "        Checking vertices" << std::endl;
  ArrayView<const zgno_t> vertexGids;
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
  TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexList size", 1)

  // We know model stores things in same order we gave it.
  if (idsAreConsecutive){
    zgno_t minLocalGID = rowMap->getMinGlobalIndex();
    for (zlno_t i=0; i < nLocalRows; i++){
      if (vertexGids[i] != minLocalGID + i) {
        fail = 1;
        break;
      }
    }
  }
  else{  // round robin ids
    if (consecutiveIdsRequested) {
      zgno_t myFirstRow;
      zgno_t tnLocalRows = nLocalRows;
      scan(*comm, Teuchos::REDUCE_SUM, 1, &tnLocalRows, &myFirstRow);
      myFirstRow -= nLocalRows;
      for (zlno_t i=0; i < nLocalRows; i++){
        if (vertexGids[i] != myFirstRow+i){
          std::cout << rank << " Row " << i << " of " << nLocalRows
                    << " myFirstRow+i " << myFirstRow+i
                    << " vertexGids " << vertexGids[i] 
                    << std::endl;
          fail = 1;
          break;
        }
      }
    }
    else {
      zgno_t myGid = rank;
      for (zlno_t i=0; i < nLocalRows; i++, myGid += nprocs){
        if (vertexGids[i] != myGid){
          std::cout << rank << " Row " << i << " of " << nLocalRows
                    << " myGid " << myGid << " vertexGids " << vertexGids[i] 
                    << std::endl;
          fail = 1;
          break;
        }
      }
    }
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "vertex gids", 1)

  if ((crds.size() != coordDim) || (wgts.size() != nVtxWeights))
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "coord or weight array size", 1)

  for (int i=0; !fail && i < coordDim; i++){
    for (zlno_t j=0; j < nLocalRows; j++){
      if (crds[i][j] != 100000*i + j){
        fail = 1;
        break;
      }
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "coord values", 1)

  for (int i=0; !fail && i < nVtxWeights; i++){
    if (nnzWgtIdx == i){
      for (zlno_t j=0; j < nLocalRows; j++){
        zscalar_t val = numNbors[j];
        if (removeSelfEdges && haveDiag[j])
          val -= 1;
        if (wgts[i][j] != val){
          fail = 1;
          break;
        }
      }
    }
    else{
      for (zlno_t j=0; j < nLocalRows; j++){
        if (wgts[i][j] != 200000*i + j){
          fail = 1;
          break;
        }
      }
    }
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "row weight values", 1)
  
  if (rank == 0) std::cout << "        Checking edges" << std::endl;
  ArrayView<const zgno_t> edgeGids;
  ArrayView<const int> procIds;
  ArrayView<const zlno_t> offsets;
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
  for (ArrayView<const zlno_t>::size_type i=0; i < offsets.size()-1; i++){
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
    if (rank == 0)
      std::cout << "Printing global graph now " << nGlobalRows << std::endl;
    printGraph(nLocalRows, vertexGids.getRawPtr(), NULL,
      edgeGids.getRawPtr(), procIds.getRawPtr(), offsets.getRawPtr(), comm);
  }
  else{
    if (rank==0) 
      std::cout << "    " << nGlobalRows << " total rows" << std::endl;
  }

  // Get graph restricted to this process

  if (rank == 0) std::cout << "        Checking local edges" << std::endl;
  ArrayView<const zlno_t> localEdges;
  ArrayView<const zlno_t> localOffsets;
  size_t numLocalNeighbors=0;

  try{
    numLocalNeighbors= model->getLocalEdgeList(localEdges, localOffsets, wgts);
  }
  catch(std::exception &e){
    std::cout << rank << ") Error " << e.what() << std::endl;
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalEdgeList", 1)
  TEST_FAIL_AND_EXIT(*comm, ((localOffsets.size()-1) == nLocalRows),
                     "getLocalEdgeList localOffsets.size", 1)

  num = 0;
  for (zlno_t i=0; i < nLocalRows; i++){
    size_t edgeListSize = localOffsets[i+1] - localOffsets[i];
    num += edgeListSize;
    size_t val = numLocalNbors[i];
    if (removeSelfEdges && haveDiag[i])
      val--;
    if (edgeListSize != val){
        std::cout << rank << "vtx " << i << " of " << localOffsets.size()
                  << " Number of local edges in model " << edgeListSize
                  << " not equal to expected number of edges " << val 
                  << std::endl;
      fail = 1;
      break;
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalEdgeList list sizes", 1)

  TEST_FAIL_AND_EXIT(*comm, numLocalNeighbors==num,
                     "getLocalEdgeList sum size", 1)

  fail = ((removeSelfEdges ? size_t(totalLocalNbors-numLocalDiags)
                           : size_t(totalLocalNbors))
          != numLocalNeighbors);
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalEdgeList total size", 1)

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

  if (rank == 0) std::cout << "        Cleaning up" << std::endl;
  delete model;

  if (nLocalRows){
    delete [] numNbors;
    delete [] numLocalNbors;
    delete [] haveDiag;

    if (nVtxWeights > 0){
      for (int i=0; i < nVtxWeights; i++){
        if (rowWeights[i])
          delete [] rowWeights[i];
      }
      delete [] rowWeights;
    }

    if (coordDim > 0){
      delete via;
      delete [] gids;
      for (int i=0; i < coordDim; i++){
        if (coords[i])
          delete [] coords[i];
      }
      delete [] coords;
    }
  }

  if (rank==0) std::cout << "    OK" << std::endl;
}

/////////////////////////////////////////////////////////////////////////////
void testGraphModel(string fname, zgno_t xdim, zgno_t ydim, zgno_t zdim,
    const RCP<const Comm<int> > &comm,
    int nVtxWeights, int nnzWgtIdx, int coordDim,
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

    cout << "Num Vertex Weights: " << nVtxWeights << endl;
    if (nnzWgtIdx >= 0)
     cout << "  Dimension " << nnzWgtIdx << " is number of neighbors" << endl;

    cout << "Coordinate dim: " << coordDim << endl;
    cout << "Request consecutive vertex gids: ";
    cout << (consecutiveIdsRequested ? "yes" : "no") << endl;
    cout << "Request to remove self edges: ";
    cout << (removeSelfEdges ? "yes" : "no") << endl;
  }

  // Input generator
  UserInputForTests *uinput;

  if (fname.size() > 0)
    uinput = new UserInputForTests(testDataFilePath, fname, comm, true);
  else
    uinput = new UserInputForTests(xdim,ydim,zdim,string(""), comm, true, true);

  RCP<tcrsMatrix_t> M = uinput->getUITpetraCrsMatrix();

  // Row Ids of test input are already consecutive

  RCP<const tcrsMatrix_t> Mconsec = rcp_const_cast<const tcrsMatrix_t>(M);

  RCP<const Tpetra::CrsGraph<zlno_t, zgno_t> > graph = Mconsec->getCrsGraph();

  printTpetraGraph<zlno_t, zgno_t>(comm, *graph, cout, 100, 
    "Graph with consecutive IDs");

  if (rank == 0) 
    std::cout << "   TEST MatrixAdapter with Consecutive IDs" << std::endl;
  bool idsAreConsecutive = true;

  testAdapter<baseMAdapter_t,xMAdapter_t,tcrsMatrix_t>(Mconsec, graph, comm,
                                                       idsAreConsecutive,
                                                       nVtxWeights, 0, 
                                                       nnzWgtIdx, coordDim,
                                                       consecutiveIdsRequested,
                                                       removeSelfEdges);
  if (rank == 0) 
    std::cout << "   TEST GraphAdapter with Consecutive IDs" << std::endl;
  testAdapter<baseGAdapter_t,xGAdapter_t,tcrsGraph_t>(graph, graph, comm,
                                                      idsAreConsecutive,
                                                      nVtxWeights, 1, 
                                                      nnzWgtIdx, coordDim,
                                                      consecutiveIdsRequested,
                                                      removeSelfEdges);

  // Do a round robin migration so that global IDs are not consecutive.

  Array<zgno_t> myNewRows;
  int nprocs = comm->getSize();
  for (size_t i=rank; i < Mconsec->getGlobalNumRows(); i+=nprocs)
    myNewRows.push_back(i);

  RCP<const tcrsMatrix_t> Mnonconsec = rcp_const_cast<const tcrsMatrix_t>(
    Zoltan2::XpetraTraits<tcrsMatrix_t>::doMigration(
      *Mconsec, myNewRows.size(), myNewRows.getRawPtr()));

  graph = Mnonconsec->getCrsGraph();

  printTpetraGraph<zlno_t, zgno_t>(comm, *graph, cout, 100, 
    "Graph with non-consecutive IDs");

  if (rank == 0)
    std::cout << "   TEST MatrixAdapter with Round-Robin IDs" << std::endl;
  idsAreConsecutive = false;

  testAdapter<baseMAdapter_t,xMAdapter_t,tcrsMatrix_t>(Mnonconsec, graph, comm,
                                                       idsAreConsecutive,
                                                       nVtxWeights, 0, 
                                                       nnzWgtIdx, coordDim,
                                                       consecutiveIdsRequested,
                                                       removeSelfEdges);

  if (rank == 0)
    std::cout << "   TEST GraphAdapter with Round-Robin IDs" << std::endl;
  testAdapter<baseGAdapter_t,xGAdapter_t,tcrsGraph_t>(graph, graph, comm,
                                                      idsAreConsecutive,
                                                      nVtxWeights, 0, 
                                                      nnzWgtIdx, coordDim,
                                                      consecutiveIdsRequested,
                                                      removeSelfEdges);

  delete uinput;
}

/////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  int nVtxWeights=0;
  int nnzWgtIdx = -1; 
  int coordDim=0;
  bool consecutiveIdsRequested=false, removeSelfEdges=false;
  string fname("simple");

  if (rank == 0)
    std::cout << "TESTING base case" << std::endl;
  testGraphModel(fname, 0, 0, 0, comm,
    nVtxWeights, nnzWgtIdx, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  if (rank == 0)
    std::cout << "TESTING with row weights" << std::endl;
  nVtxWeights = 1;

  testGraphModel(fname, 0, 0, 0, comm,
    nVtxWeights, nnzWgtIdx, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  if (rank == 0)
    std::cout << "TESTING with weights = nnz" << std::endl;
  nnzWgtIdx = 1;

  testGraphModel(fname, 0, 0, 0, comm,
    nVtxWeights, nnzWgtIdx, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  if (rank == 0)
    std::cout << "TESTING with multiple row weights and coords" << std::endl;
  nVtxWeights = 2;
  coordDim = 3;

  testGraphModel(fname, 0, 0, 0, comm,
    nVtxWeights, nnzWgtIdx, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  if (rank == 0)
    std::cout << "TESTING with consecutiveIdsRequested" << std::endl;
  consecutiveIdsRequested = true;

  testGraphModel(fname, 0, 0, 0, comm,
    nVtxWeights, nnzWgtIdx, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  if (rank == 0)
    std::cout << "TESTING with removeSelfEdges" << std::endl;
  removeSelfEdges = true;

  testGraphModel(fname, 0, 0, 0, comm,
    nVtxWeights, nnzWgtIdx, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  if (rank == 0)
    std::cout << "TESTING with consecutiveIdsRequested=false" << std::endl;
  consecutiveIdsRequested = false;

  testGraphModel(fname, 0, 0, 0, comm,
    nVtxWeights, nnzWgtIdx, coordDim,
    consecutiveIdsRequested, removeSelfEdges);

  if (rank==0)
    cout << "PASS" << endl;

  return 0;
}

