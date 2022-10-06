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
// Testing of HyperGraphModel built from APF mesh adapters.
//

/*! \brief Test of HyperGraphModel interface.
 *
 *  \todo test all methods of HyperGraphModel
 *  \todo add functionality to choose action from command line
 */

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_HyperGraphModel.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_PamgenMeshAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_InputTraits.hpp>
// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayView.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;
typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> simpleUser_t;
typedef Zoltan2::BasicUserTypes<> basic_user_t;
typedef Zoltan2::MeshAdapter<simpleUser_t>                  baseMeshAdapter_t;

typedef Zoltan2::PamgenMeshAdapter<basic_user_t> inputAdapter_t;
typedef inputAdapter_t::base_adapter_t mesh_base_adapter_t;

// TO change
typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t,simpleUser_t>   xGAdapter_t;
typedef typename inputAdapter_t::offset_t    zoffset_t;


void testUsingModelFlags(Zoltan2::ModelFlags modelFlag) {

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();

  // Creating mesh adapter
  if (rank == 0) std::cout << "Creating mesh adapter ... \n\n";

  typedef Zoltan2::PamgenMeshAdapter<basic_user_t> inputAdapter_t;

  inputAdapter_t ia(*comm, "region");

  if (rank == 0) std::cout << "REGION-BASED TEST" << std::endl;

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags;
  modelFlags.set(modelFlag);

  int fail=0;

  RCP<Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));
  Teuchos::ParameterList &pl = env->getParametersNonConst();
  pl.set("hypergraph_model_type","ghosting");

  const baseMeshAdapter_t *base_ia = dynamic_cast<const baseMeshAdapter_t*>(&ia);
  Zoltan2::modelFlag_t graphFlags_;
  RCP<const baseMeshAdapter_t> baseInputAdapter_(base_ia,false);

  RCP<const Zoltan2::Environment> envConst = Teuchos::rcp_const_cast<const Zoltan2::Environment>(env);
  Zoltan2::HyperGraphModel<inputAdapter_t> model(baseInputAdapter_,envConst,comm,
                                                 graphFlags_,Zoltan2::HYPEREDGE_CENTRIC);

  typedef Zoltan2::StridedData<zlno_t, zscalar_t> input_t;

  // TEST of getVertexList and getVertexListKokkos
  ArrayView<const zgno_t> vertexIds;
  ArrayView<input_t> wgts;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kVertexIds;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> kWgts;
  auto nbVertices = model.getVertexList(vertexIds, wgts);
  auto kNbVertices = model.getVertexListKokkos(kVertexIds, kWgts);
  auto nbWeightPerVertex = model.getNumWeightsPerVertex();
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
  auto nbVerticesFromCoords = model.getVertexCoords(cxyz);
  auto kNbVerticesFromKCoords = model.getVertexCoordsKokkos(kcxyz);
  if(nbVerticesFromCoords != kNbVerticesFromKCoords)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getVertexCoords != getVertexCoordsKokkos", 1)
  }

  for(size_t d = 0; d < model.getCoordinateDim(); ++d) {
    for(size_t i = 0; i < nbVerticesFromCoords; ++i) {
      if (cxyz[d][i] != kcxyz(i, d))
        fail = 1;
    }
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "coordinates ouputs of getVertexCoords  and getVertexCoordsKokkos are different", 1)

  // TEST of getEdgeList and getEdgeListKokkos
  ArrayView<const zgno_t> edgesIds;
  ArrayView<input_t> eWgts;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kEdgeIds;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> kEWgts;
  auto nbEdges =  model.getEdgeList(edgesIds, eWgts);
  auto kNbEdges =  model.getEdgeListKokkos(kEdgeIds, kEWgts);
  if(nbEdges != kNbEdges)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getEdgeList != getEdgeListKokkos", 1)
  }
  for(size_t i = 0; i < nbEdges; ++i) {
    if (edgesIds[i] != kEdgeIds(i))
      fail = 1;
  }

  // TEST of getPinList and getPinListKokkos
  ArrayView<const zgno_t> pIds;
  ArrayView<const zoffset_t> pOffsets;
  ArrayView<input_t> pWgts;
  Kokkos::View<const zgno_t *, typename znode_t::device_type> kPIds;
  Kokkos::View<const zoffset_t *, typename znode_t::device_type> kPOffsets;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> kPWgts;
  auto nbPins =  model.getPinList(pIds, pOffsets, pWgts);
  auto kNbPins =  model.getPinListKokkos(kPIds, kPOffsets, kPWgts);
  if(nbPins != kNbPins)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getPinList != getPinListKokkos", 1)
  }
  for(size_t i = 0; i < nbPins; ++i) {
    if (pIds[i] != kPIds(i))
      fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getPinList and  getPinListKokkos are different for pIds", 1)

  for(size_t i = 0; i < kPOffsets.size(); ++i) {
    if (pOffsets[i] != kPOffsets(i))
      fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getPinList and  getPinListKokkos are different for pOffsets", 1)
  if(kPWgts.size() !=0 ||  pWgts.size() != 0)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Weight for edges detected", 1)
  }

  // TEST of getPinList and getOwnedListKokkos
  ArrayView<bool> isOwner;
  Kokkos::View<bool *, typename znode_t::device_type> kIsOwner;
  auto nbIsOwner =  model.getOwnedList(isOwner);
  auto kNbIsOwner =  model.getOwnedListKokkos(kIsOwner);
  if(nbIsOwner != kNbIsOwner)
  {
    fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "Return of getOwnedList != getOwnedListKokkos", 1)
  }

  for(size_t i = 0; i < nbIsOwner; ++i) {
    if (isOwner[i] != kIsOwner(i))
      fail = 1;
  }

}

int main(int narg, char *arg[]) {

  Tpetra::ScopeGuard tscope(&narg, &arg);

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

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

  testUsingModelFlags(Zoltan2::GENERATE_CONSECUTIVE_IDS);
  testUsingModelFlags(Zoltan2::REMOVE_SELF_EDGES);
  testUsingModelFlags(Zoltan2::BUILD_LOCAL_GRAPH);

  if (rank==0)
    std::cout<<"PASS\n";

  return 0;
}
