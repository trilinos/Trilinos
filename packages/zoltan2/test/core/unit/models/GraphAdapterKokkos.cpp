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
// Testing of GraphAdapter built from Xpetra matrix input adapters.
//

/*! \brief Test of GraphAdapter interface.
 *
 */

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


const int SMALL_NUMBER_OF_ROWS = 5;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;

typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> simpleUser_t;

typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t>     tcrsMatrix_t;
typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t>                 tcrsGraph_t;
typedef Tpetra::Map<zlno_t, zgno_t, znode_t>                      tmap_t;

typedef Zoltan2::BasicVectorAdapter<simpleUser_t>              simpleVAdapter_t;

typedef Zoltan2::MatrixAdapter<tcrsMatrix_t,simpleUser_t>      baseMAdapter_t;
typedef Zoltan2::GraphAdapter<tcrsGraph_t,simpleUser_t>        baseGAdapter_t;

typedef Zoltan2::XpetraCrsMatrixAdapter<tcrsMatrix_t,simpleUser_t> xMAdapter_t;
typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t,simpleUser_t>   xGAdapter_t;
typedef typename xGAdapter_t::offset_t    zoffset_t;

using std::string;
using std::vector;

/////////////////////////////////////////////////////////////////////////////
int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();

  int nVtxWeights=5;
  int nnzWgtIdx = -1;
  string fname("simple");

  if (rank == 0)
    std::cout << "TESTING base case (global)" << std::endl;

  // Input generator
  UserInputForTests *uinput;

  uinput = new UserInputForTests(testDataFilePath, fname, comm, true);

  RCP<tcrsMatrix_t> M = uinput->getUITpetraCrsMatrix();
  zlno_t nLocalRows = M->getLocalNumRows();

  // Weights:
  zscalar_t **rowWeights=NULL;
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

  RCP<const tcrsMatrix_t> Mconsec = rcp_const_cast<const tcrsMatrix_t>(M);
  RCP<const Tpetra::CrsGraph<zlno_t, zgno_t> > graph = Mconsec->getCrsGraph();
  xGAdapter_t tmi(graph, nVtxWeights);
  for (int i=0; i < nVtxWeights; i++){
    if (nnzWgtIdx == i)
      tmi.setWeightIsDegree(i);
    else
      tmi.setWeights(rowWeights[i], 1, i);
  }

  simpleVAdapter_t *via = NULL;

  // Set up some fake input
  zscalar_t **coords=NULL;
  int coordDim= 3;

  if (coordDim > 0){
    coords = new zscalar_t * [coordDim];
    for (int i=0; i < coordDim; i++){
      coords[i] = new zscalar_t [nLocalRows];
      for (zlno_t j=0; j < nLocalRows; j++){
        coords[i][j] = 100000*i + j;
      }
    }
  }


  zgno_t *gids = NULL;
  if (coordDim > 0) {
    gids = new zgno_t[nLocalRows];
    for (zlno_t i = 0; i < nLocalRows; i++)
      gids[i] = M->getRowMap()->getGlobalElement(i);
    via = new simpleVAdapter_t(nLocalRows, gids, coords[0],
                                           (coordDim > 1 ? coords[1] : NULL),
                                           (coordDim > 2 ? coords[2] : NULL),
                                            1,1,1);
    tmi.setCoordinateInput(via);
  }


  // TEST of getIDsView, getIDsKokkosView (compatible with all adapters) and getVertexIDsView (graphAdapter)

  const zgno_t *vertexIds;
  const zgno_t *ids;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kIds;

  int fail=0;
  if (tmi.getLocalNumVertices() != tmi.getLocalNumIDs())
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumVertices != tmi.getLocalNumIDs() for graphAdapter", 1)


  tmi.getVertexIDsView(vertexIds);
  tmi.getIDsView(ids);
  tmi.getIDsKokkosView(kIds);
  for(size_t i = 0; i < tmi.getLocalNumIDs(); ++i) {
   if (ids[i] != vertexIds[i])
     fail = 1;

   if (kIds(i) != vertexIds[i])
     fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "ids != vertexIds != kIds", 1)

// TEST of getWeightsView, getWeightsKokkosView and getVertexWeightsView (GRAPH_VERTEX)
  for (int w=0; !fail && w < nVtxWeights; w++){
    const zscalar_t *wgts;
    const zscalar_t *vwgts;
    Kokkos::View<zscalar_t **, typename znode_t::device_type> wkgts;
    int stride;
    tmi.getWeightsView(wgts, stride, w);
    tmi.getWeightsKokkosView(wkgts);
    tmi.getVertexWeightsView(vwgts, stride, w);
    for(size_t i = 0; i < tmi.getLocalNumIDs(); ++i) {
      if (wgts[stride*i] != vwgts[stride*i])
        fail = 1;

      if (wkgts(i, w) != vwgts[stride*i])
        fail = 1;
    }

  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "wgts != vwgts != wkgts", 1)

// TEST of getEdgesKokkosView, getEdgesView (GRAPH_VERTEX)
  const zoffset_t *offsets;
  const zgno_t *adjids;
  Kokkos::View<const zoffset_t *, typename znode_t::device_type> kOffsets;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kAdjids;

  tmi.getEdgesView(offsets, adjids);
  tmi.getEdgesKokkosView(kOffsets, kAdjids);
  for(size_t i = 0; i < tmi.getLocalNumVertices() + 1; ++i) {
    if (offsets[i] != kOffsets(i))
      fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "offsets != kOffsets", 1)

  for(size_t i = 0; i < tmi.getLocalNumEdges(); ++i) {
    if (adjids[i] != kAdjids(i))
      fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "adjids != kAdjids", 1)

  // TEST of getCoordinatesView, getVertexCoords,
  // getCoordinatesKokkosView (AdapterWithCoordsWrapper) and getVertexCoordsKokkos

  typedef Zoltan2::StridedData<zlno_t, zscalar_t> input_t;
  input_t *coordArray = new input_t [coordDim];

  for (int dim=0; dim < coordDim; dim++){
    int stride;
    const zscalar_t *coords=NULL;
    try{
      tmi.getCoordinateInput()->getCoordinatesView(coords, stride, dim);
    }
    Z2_FORWARD_EXCEPTIONS;

    ArrayRCP<const zscalar_t> cArray(coords, 0, tmi.getLocalNumIDs()*stride, false);
    coordArray[dim] = input_t(cArray, stride);
  }
  ArrayRCP<input_t> crds = arcp<input_t>(coordArray, 0, coordDim, true);
  Kokkos::View<zscalar_t **, Kokkos::LayoutLeft, typename znode_t::device_type> kCrds;

  tmi.getCoordinateInput()->getCoordinatesKokkosView(kCrds);

  for(size_t i = 0; i < tmi.getLocalNumIDs(); ++i) {
    for (int j = 0; j < coordDim; j++) {
      if (crds[j][i] != kCrds(i, j))
        fail = 1;
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "crds != kCrds", 1)

  // clean
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

  delete uinput;

  if (rank==0)
  std::cout << "PASS" << std::endl;

  return 0;
}

