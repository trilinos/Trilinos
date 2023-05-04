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
#include <Zoltan2_PamgenMeshAdapter.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <bitset>


// Teuchos includes
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayView.hpp>


const int SMALL_NUMBER_OF_ROWS = 5;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;

typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> simpleUser_t;
typedef Zoltan2::BasicUserTypes<> basic_user_t;

typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t>     tcrsMatrix_t;
typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t>                 tcrsGraph_t;
typedef Tpetra::Map<zlno_t, zgno_t, znode_t>                      tmap_t;

typedef Zoltan2::BasicVectorAdapter<simpleUser_t>              simpleVAdapter_t;

typedef Zoltan2::MatrixAdapter<tcrsMatrix_t,simpleUser_t>      baseMAdapter_t;
typedef Zoltan2::GraphAdapter<tcrsGraph_t,simpleUser_t>        baseGAdapter_t;
typedef Zoltan2::MeshAdapter<simpleUser_t>                  baseMeshAdapter_t;

typedef Zoltan2::XpetraCrsMatrixAdapter<tcrsMatrix_t,simpleUser_t> xMAdapter_t;
typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t,simpleUser_t>   xGAdapter_t;

typedef Zoltan2::PamgenMeshAdapter<basic_user_t> inputAdapter_t;
typedef inputAdapter_t::base_adapter_t mesh_base_adapter_t;

typedef typename xGAdapter_t::offset_t    zoffset_t;

using std::string;
using std::vector;

void testGraphModelWithMatrixAdapter(Zoltan2::ModelFlags modelFlag, string fname,  const RCP<const Comm<int> > &comm, int nVtxWeights
                                    ,int nnzWgtIdx) {
  int rank = comm->getRank();
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
  xMAdapter_t tmi(Mconsec, nVtxWeights);
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

  if (rank == 0) std::cout << "        Creating GraphModel" << std::endl;
  Zoltan2::GraphModel<baseMAdapter_t> *model = NULL;
  RCP<const baseMAdapter_t> baseTmi = rcp(dynamic_cast<baseMAdapter_t *>(&tmi),false);
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags;
  modelFlags.set(modelFlag);

  int fail=0;
  try{
    model = new Zoltan2::GraphModel<baseMAdapter_t>(baseTmi, env,
                                                 comm, modelFlags);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "Creating graph model", 1)

  typedef Zoltan2::StridedData<zlno_t, zscalar_t> input_t;

  // TEST of getVertexList and getVertexListKokkos
  ArrayView<const zgno_t> vertexIds;
  ArrayView<input_t> wgts;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kVertexIds;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> kWgts;
  auto nbVertices = model->getVertexList(vertexIds, wgts);
  auto kNbVertices = model->getVertexListKokkos(kVertexIds, kWgts);
  auto nbWeightPerVertex = model->getNumWeightsPerVertex();
  if(nbVertices != kNbVertices)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getVertexList != getVertexListKokkos", 1)
  }

  for(size_t i = 0; i < nbVertices; ++i) {
     if (vertexIds[i] != kVertexIds(i))
         fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "vertexIds != kVertexIds", 1)

  for(size_t w = 0; w < nbWeightPerVertex; ++w) {
      for(size_t i = 0; i < nbVertices; ++i) {
          if (wgts[w][i] != kWgts(i, w))
              fail = 1;
      }

  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "weight outputs of getVertexList and getVertexListKokkos are different", 1)


  // TEST of getVertexCoords and getVertexCoordsKokkos
  ArrayView<input_t> cxyz;
  Kokkos::View<zscalar_t **,Kokkos::LayoutLeft, typename znode_t::device_type> kcxyz;
  auto nbVerticesFromCoords = model->getVertexCoords(cxyz);
  auto kNbVerticesFromKCoords = model->getVertexCoordsKokkos(kcxyz);
  if(nbVerticesFromCoords != kNbVerticesFromKCoords)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getVertexCoords != getVertexCoordsKokkos", 1)
  }

  for(size_t d = 0; d < model->getCoordinateDim(); ++d) {
    for(size_t i = 0; i < nbVerticesFromCoords; ++i) {
      if (cxyz[d][i] != kcxyz(i, d))
        fail = 1;
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "coordinates ouputs of getVertexCoords  and getVertexCoordsKokkos are different", 1)

  // TEST of getEdgeList and getEdgeListKokkos
  ArrayView<const zgno_t> edgesIds;
  ArrayView<const zoffset_t> offsets;
  ArrayView<input_t> eWgts;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kEdgeIds;
  Kokkos::View<const zoffset_t *, typename znode_t::device_type> kOffsets;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> kEWgts;
  auto nbEdges =  model->getEdgeList(edgesIds, offsets, eWgts);
  auto kNbEdges =  model->getEdgeListKokkos(kEdgeIds, kOffsets, kEWgts);
  if(nbEdges != kNbEdges)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getEdgeList != getEdgeListKokkos", 1)
  }
  for(size_t i = 0; i < nbEdges; ++i) {
    if (edgesIds[i] != kEdgeIds(i))
      fail = 1;
  }

  for(size_t i = 0; i <offsets.size(); ++i) {
    if (offsets[i] != kOffsets(i))
      fail = 1;
  }
  if (model->getNumWeightsPerEdge() != 0)
      fail = 1;

  TEST_FAIL_AND_EXIT(*comm, !fail, "getEdgeList != getEdgeListKokkos", 1)

  // TEST of getVertexDist and getVertexDistKokkos: only available for consecutiveIds
  bool consecutiveIdsRequired = modelFlags.test(Zoltan2::GENERATE_CONSECUTIVE_IDS);
  if(consecutiveIdsRequired)  {
    ArrayView<size_t> vtxDist;
    Kokkos::View<size_t *, typename znode_t::device_type> kVtxDist;
    model->getVertexDist(vtxDist);
    model->getVertexDistKokkos(kVtxDist);
    for(size_t i = 0; i < comm->getSize() + 1; ++i) {
      if (vtxDist[i] != kVtxDist(i))
        fail = 1;
    }
    TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexDist != getVertexDistKokkos", 1)
  }

      delete model;
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
}

void testGraphModelWithGraphAdapter(Zoltan2::ModelFlags modelFlag, string fname,  const RCP<const Comm<int> > &comm, int nVtxWeights
                                    ,int nnzWgtIdx) {
  int rank = comm->getRank();
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
  // We don't use the edges weight parameter because the getWeightsView of the GraphAdapter is not able to call
  // the appropriate getEdgeWeightsView method.
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

  if (rank == 0) std::cout << "        Creating GraphModel" << std::endl;
  Zoltan2::GraphModel<baseGAdapter_t> *model = NULL;
  RCP<const baseGAdapter_t> baseTmi = rcp(dynamic_cast<baseGAdapter_t *>(&tmi),false);
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags;
  modelFlags.set(modelFlag);

  int fail=0;
  try{
    model = new Zoltan2::GraphModel<baseGAdapter_t>(baseTmi, env,
                                                 comm, modelFlags);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "Creating graph model", 1)

  typedef Zoltan2::StridedData<zlno_t, zscalar_t> input_t;

  // TEST of getVertexList and getVertexListKokkos
  ArrayView<const zgno_t> vertexIds;
  ArrayView<input_t> wgts;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kVertexIds;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> kWgts;
  auto nbVertices = model->getVertexList(vertexIds, wgts);
  auto kNbVertices = model->getVertexListKokkos(kVertexIds, kWgts);
  auto nbWeightPerVertex = model->getNumWeightsPerVertex();
  if(nbVertices != kNbVertices)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getVertexList != getVertexListKokkos", 1)
  }

  for(size_t i = 0; i < nbVertices; ++i) {
    if (vertexIds[i] != kVertexIds(i))
      fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "vertexIds != kVertexIds", 1)

  for(size_t w = 0; w < nbWeightPerVertex; ++w) {
    for(size_t i = 0; i < nbVertices; ++i) {
      if (wgts[w][i] != kWgts(i, w))
        fail = 1;
    }

  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "weight outputs of getVertexWeightsView and getVertexCoordsKokkos are different", 1)

  // TEST of getVertexCoords and getVertexCoordsKokkos
  ArrayView<input_t> cxyz;
  Kokkos::View<zscalar_t **,Kokkos::LayoutLeft, typename znode_t::device_type> kcxyz;
  auto nbVerticesFromCoords = model->getVertexCoords(cxyz);
  auto kNbVerticesFromKCoords = model->getVertexCoordsKokkos(kcxyz);
  if(nbVerticesFromCoords != kNbVerticesFromKCoords)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getVertexCoords != getVertexCoordsKokkos", 1)
  }

  for(size_t d = 0; d < model->getCoordinateDim(); ++d) {
    for(size_t i = 0; i < nbVerticesFromCoords; ++i) {
      if (cxyz[d][i] != kcxyz(i, d))
        fail = 1;
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "coordinates ouputs of getVertexCoords  and getVertexCoordsKokkos are different", 1)

  // TEST of getEdgeList and getEdgeListKokkos
  ArrayView<const zgno_t> edgesIds;
  ArrayView<const zoffset_t> offsets;
  ArrayView<input_t> eWgts;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kEdgeIds;
  Kokkos::View<const zoffset_t *, typename znode_t::device_type> kOffsets;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> kEWgts;
  auto nbEdges =  model->getEdgeList(edgesIds, offsets, eWgts);
  auto kNbEdges =  model->getEdgeListKokkos(kEdgeIds, kOffsets, kEWgts);
  if(nbEdges != kNbEdges)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getEdgeList != getEdgeListKokkos", 1)
  }
  for(size_t i = 0; i < nbEdges; ++i) {
    if (edgesIds[i] != kEdgeIds(i))
      fail = 1;
  }

  for(size_t i = 0; i < model->getLocalNumObjects() +1; ++i) {
    if (offsets[i] != kOffsets(i))
      fail = 1;
  }

  // Test of eWgts and kEWgts is currently not possible because getWeightsView method is not implemented
  // for GRAPH_EDGE GraphEntityType

  // TEST of getVertexDist and getVertexDistKokkos: only available for consecutiveIds
  bool consecutiveIdsRequired = modelFlags.test(Zoltan2::GENERATE_CONSECUTIVE_IDS);
  if(consecutiveIdsRequired)  {
    ArrayView<size_t> vtxDist;
    Kokkos::View<size_t *, typename znode_t::device_type> kVtxDist;
    model->getVertexDist(vtxDist);
    model->getVertexDistKokkos(kVtxDist);
    for(size_t i = 0; i < comm->getSize() + 1; ++i) {
      if (vtxDist[i] != kVtxDist(i))
        fail = 1;
    }
    TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexDist != getVertexDistKokkos", 1)
  }

  delete model;
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
}

void testGraphModelWithMeshAdapter(Zoltan2::ModelFlags modelFlag, string fname,  const RCP<const Comm<int> > &comm, int nVtxWeights
                                    ,int nnzWgtIdx) {

  int rank = comm->getRank();
  int numProcs = comm->getSize();
  /***************************************************************************/
  /*************************** GET XML INPUTS ********************************/
  /***************************************************************************/

  // default values for command-line arguments
  std::string xmlMeshInFileName("Poisson.xml");

  // Read xml file into parameter list
  Teuchos::ParameterList inputMeshList;

  if(xmlMeshInFileName.length()) {
    if (rank == 0) {
      std::cout << "\nReading parameter list from the XML file \""
                <<xmlMeshInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile(xmlMeshInFileName,
                                         Teuchos::inoutArg(inputMeshList));
    if (rank == 0) {
      inputMeshList.print(std::cout,2,true,true);
      std::cout << "\n";
    }
  }
  else {
    std::cout << "Cannot read input file: " << xmlMeshInFileName << "\n";
  }

  // Get pamgen mesh definition
  std::string meshInput = Teuchos::getParameter<std::string>(inputMeshList,
                                                             "meshInput");
  // Get dimensions
  int dim = 3;
  // Input generator
  // Generate mesh with Pamgen
  long long maxInt = std::numeric_limits<long long>::max();
  Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);

  // Creating mesh adapter
  if (rank == 0) std::cout << "Creating mesh adapter ... \n\n";

  typedef Zoltan2::PamgenMeshAdapter<basic_user_t> inputAdapter_t;

  inputAdapter_t ia(*comm, "region");

  if (rank == 0) std::cout << "REGION-BASED TEST" << std::endl;

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags;
  modelFlags.set(modelFlag);

  int fail=0;

  if (rank == 0) std::cout << "REGION-BASED TEST" << std::endl;
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));


  if (rank == 0) std::cout << "        Creating GraphModel" << std::endl;
  Zoltan2::GraphModel<mesh_base_adapter_t> *model = NULL;
  RCP<const mesh_base_adapter_t> baseInputAdapter = rcp(dynamic_cast<mesh_base_adapter_t *>(&ia),false);

  try{
    model = new Zoltan2::GraphModel<mesh_base_adapter_t>(baseInputAdapter, env, comm, modelFlags);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "Creating graph model", 1)

  typedef Zoltan2::StridedData<zlno_t, zscalar_t> input_t;

  // TEST of getVertexList and getVertexListKokkos
  ArrayView<const zgno_t> vertexIds;
  ArrayView<input_t> wgts;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kVertexIds;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> kWgts;
  auto nbVertices = model->getVertexList(vertexIds, wgts);
  auto kNbVertices = model->getVertexListKokkos(kVertexIds, kWgts);
  auto nbWeightPerVertex = model->getNumWeightsPerVertex();
  if(nbVertices != kNbVertices)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getVertexList != getVertexListKokkos", 1)
  }

  for(size_t i = 0; i < nbVertices; ++i) {
  if (vertexIds[i] != kVertexIds(i))
      fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "vertexIds != kVertexIds", 1)

  for(size_t w = 0; w < nbWeightPerVertex; ++w) {
    for(size_t i = 0; i < nbVertices; ++i) {
      if (wgts[w][i] != kWgts(i, w))
        fail = 1;
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "weight outputs of getVertexWeightsView and getVertexCoordsKokkos are different", 1)

  // TEST of getVertexCoords and getVertexCoordsKokkos
  ArrayView<input_t> cxyz;
  Kokkos::View<zscalar_t **,Kokkos::LayoutLeft, typename znode_t::device_type> kcxyz;
  auto nbVerticesFromCoords = model->getVertexCoords(cxyz);
  auto kNbVerticesFromKCoords = model->getVertexCoordsKokkos(kcxyz);
  if(nbVerticesFromCoords != kNbVerticesFromKCoords)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getVertexCoords != getVertexCoordsKokkos", 1)
  }

  for(size_t d = 0; d < model->getCoordinateDim(); ++d) {
    for(size_t i = 0; i < nbVerticesFromCoords; ++i) {
      if (cxyz[d][i] != kcxyz(i, d))
        fail = 1;
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "coordinates ouputs of getVertexCoords  and getVertexCoordsKokkos are different", 1)

  // TEST of getEdgeList and getEdgeListKokkos
  ArrayView<const zgno_t> edgesIds;
  ArrayView<const inputAdapter_t::offset_t> offsets;
  ArrayView<input_t> eWgts;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kEdgeIds;
  Kokkos::View<const inputAdapter_t::offset_t *, typename znode_t::device_type> kOffsets;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> kEWgts;
  auto nbEdges =  model->getEdgeList(edgesIds, offsets, eWgts);
  auto kNbEdges =  model->getEdgeListKokkos(kEdgeIds, kOffsets, kEWgts);
  if(nbEdges != kNbEdges)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getEdgeList != getEdgeListKokkos", 1)
  }
  for(size_t i = 0; i < nbEdges; ++i) {
    if (edgesIds[i] != kEdgeIds(i))
        fail = 1;
  }

  for(size_t i = 0; i < model->getLocalNumObjects() +1; ++i) {
           if (offsets[i] != kOffsets(i))
               fail = 1;
  }

  // Test of eWgts and kEWgts is currently not possible because getWeightsView method is not implemented
  // for GRAPH_EDGE GraphEntityType

  // TEST of getVertexDist and getVertexDistKokkos: only available for consecutiveIds
  bool consecutiveIdsRequired = modelFlags.test(Zoltan2::GENERATE_CONSECUTIVE_IDS);
  if(consecutiveIdsRequired)  {
    ArrayView<size_t> vtxDist;
    Kokkos::View<size_t *, typename znode_t::device_type> kVtxDist;
    model->getVertexDist(vtxDist);
    model->getVertexDistKokkos(kVtxDist);
    for(size_t i = 0; i < comm->getSize() + 1; ++i) {
      if (vtxDist[i] != kVtxDist(i))
        fail = 1;
    }
    TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexDist != getVertexDistKokkos", 1)
  }

  delete model;
  // delete mesh
  if (rank == 0) std::cout << "Deleting the mesh ... \n\n";
  Delete_Pamgen_Mesh();
}


void testUsingModelFlags(Zoltan2::ModelFlags modelFlag) {

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int nVtxWeights=5;
  int nnzWgtIdx = -1;
  string fname("simple");

  testGraphModelWithGraphAdapter(modelFlag, fname, comm, nVtxWeights, nnzWgtIdx);
  testGraphModelWithMatrixAdapter(modelFlag, fname, comm, nVtxWeights, nnzWgtIdx);
  testGraphModelWithMeshAdapter(modelFlag, fname, comm, nVtxWeights, nnzWgtIdx);

}
/////////////////////////////////////////////////////////////////////////////
int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();

  if (rank == 0)
    std::cout << "TESTING base case (global)" << std::endl;

  testUsingModelFlags(Zoltan2::GENERATE_CONSECUTIVE_IDS);
  testUsingModelFlags(Zoltan2::REMOVE_SELF_EDGES);
  testUsingModelFlags(Zoltan2::BUILD_LOCAL_GRAPH);

  if (rank==0)
    std::cout << "PASS" << std::endl;
  return 0;
}

