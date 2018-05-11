// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

// Panzer STK
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"

// Panzer
#include "Panzer_Workset.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_IntegrationDescriptor.hpp"
#include "Panzer_BasisDescriptor.hpp"
#include "Panzer_Evaluator_DomainInterface.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_L2Projection.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_BlockedDOFManager.hpp"

#include "Tpetra_CrsMatrix.hpp"

// Teuchos
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include <memory>

TEUCHOS_UNIT_TEST(l2_projection, dof_manager)
{
  using namespace Teuchos;
  using namespace panzer;
  using namespace panzer_stk;

  RCP<MpiComm<int>> comm = rcp(new MpiComm<int>(MPI_COMM_WORLD));

  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();
  TEUCHOS_ASSERT(numProcs == 2);

  RCP<panzer_stk::STK_Interface> mesh;
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("Z Blocks",1);
    pl->set("X Elements",2);
    pl->set("Y Elements",2);
    pl->set("Z Elements",1);
    pl->set("X Procs",2);
    pl->set("Y Procs",1);
    pl->set("Z Procs",1);
    panzer_stk::CubeHexMeshFactory factory;
    factory.setParameterList(pl);
    mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  }

  // Build Worksets

  const int basisOrder = 1;
  BasisDescriptor hgradBD(basisOrder,"HGrad");
  BasisDescriptor hcurlBD(basisOrder,"HCurl");
  BasisDescriptor hdivBD(basisOrder,"HDiv");

  const int intOrder = 2;
  IntegrationDescriptor integrationDescriptor(intOrder,IntegrationDescriptor::VOLUME);

  WorksetNeeds worksetNeeds;
  worksetNeeds.addBasis(hgradBD);
  worksetNeeds.addBasis(hcurlBD);
  worksetNeeds.addBasis(hcurlBD);
  worksetNeeds.addIntegrator(integrationDescriptor);

  RCP<WorksetFactory> worksetFactory(new WorksetFactory(mesh));
  std::vector<std::string> eBlockNames;
  mesh->getElementBlockNames(eBlockNames);
  std::map<std::string,WorksetNeeds> eblockNeeds;
  for (const auto& block : eBlockNames)
    eblockNeeds[block] = worksetNeeds;
  RCP<WorksetContainer> worksetContainer(new WorksetContainer(worksetFactory,eblockNeeds));

  // Build Connection Manager
  using LO = int;
  using GO = panzer::Ordinal64;
  const RCP<panzer::ConnManager<LO,GO> > connManager = rcp(new panzer_stk::STKConnManager<GO>(mesh));

  // Set up bases for projections
  auto cellTopology = mesh->getCellTopology(eBlockNames[0]);

  auto hgradBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hgradBD.getType(),hgradBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hgradFP(new panzer::Intrepid2FieldPattern(hgradBasis));

  auto curlBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hcurlBD.getType(),hcurlBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hcurlFP(new panzer::Intrepid2FieldPattern(curlBasis));

  auto divBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hdivBD.getType(),hdivBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hdivFP(new panzer::Intrepid2FieldPattern(divBasis));


  // ****************
  // Testing ghosting
  // ****************
  out.setShowProcRank(true);
  out.setOutputToRootOnly(myRank);
  connManager->buildConnectivity(*hgradFP);
  for (const auto& block : eBlockNames) {
    const auto& localElements = connManager->getElementBlock(block);
    const auto& ghostElements = connManager->getNeighborElementBlock(block);
    out << "block=" << block << ", numOwnedElements = " << localElements.size() 
        << ", numGhostedElements = " << ghostElements.size() << std::endl;
    for (const auto& e : localElements) {
      auto cellSize = connManager->getConnectivitySize(e);
      auto conn = connManager->getConnectivity(e);
      std::stringstream os;
      os << "owned lid=" << e << ", cellSize=" << cellSize << ", conn=";
      for (int i=0; i < cellSize; ++i)
        os << conn[i] << ",";
      out << os.str() << std::endl;
    }
    for (const auto& e : ghostElements) {
      auto cellSize = connManager->getConnectivitySize(e);
      auto conn = connManager->getConnectivity(e);
      std::stringstream os;
      os << "ghosted lid=" << e << ", cellSize=" << cellSize << ", conn=";
      for (int i=0; i < cellSize; ++i)
        os << conn[i] << ",";
      out << os.str() << std::endl;
    }
  }
  out.setOutputToRootOnly(0);
  // ****************
  // ****************





  // Build source DOF Manager that mimics multi-fluid plasma dof manager
  RCP<panzer::DOFManager<LO,GO>> sourceGlobalIndexer = rcp(new panzer::DOFManager<LO,GO>(connManager,*comm->getRawMpiComm()));
  sourceGlobalIndexer->addField("Chaff0",hgradFP);
  sourceGlobalIndexer->addField("Chaff1",hgradFP);
  sourceGlobalIndexer->addField("Chaff2",hgradFP);
  sourceGlobalIndexer->addField("E_Field",hcurlFP);
  sourceGlobalIndexer->addField("B_Field",hdivFP);
  sourceGlobalIndexer->buildGlobalUnknowns();

  

  // Build Target DOF Manager (Scalar fields on hgrad)
  RCP<panzer::DOFManager<LO,GO>> targetGlobalIndexer = rcp(new panzer::DOFManager<LO,GO>(connManager,*comm->getRawMpiComm()));
  targetGlobalIndexer->addField("Projection Nodal Scalar",hgradFP);
  targetGlobalIndexer->addField("E1",hgradFP);
  targetGlobalIndexer->addField("E2",hgradFP);
  targetGlobalIndexer->addField("B0",hgradFP);
  targetGlobalIndexer->addField("B1",hgradFP);
  targetGlobalIndexer->addField("B2",hgradFP);
  targetGlobalIndexer->buildGlobalUnknowns();

  // Build projection objects
  panzer::L2Projection<LO,GO> projectionFactory;
  projectionFactory.setup(hgradBD,integrationDescriptor,comm,connManager,eBlockNames,worksetContainer);

  // Build mass matrix
  auto massMatrix = projectionFactory.buildMassMatrix();
  massMatrix->print(out);
  massMatrix->getRowMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);
  massMatrix->getColMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);

  // Test a projection using the mass matrix
  

  // Build rhs matrix
  auto rhsMatrix_Chaff1 = projectionFactory.buildRHSMatrix(sourceGlobalIndexer,Teuchos::null,"Chaff1",hgradBD); // Project from scalar basis
  auto rhsMatrix_E0 = projectionFactory.buildRHSMatrix(sourceGlobalIndexer,Teuchos::null,"E_Field",hcurlBD,true,0); // Project from vector basis
  auto rhsMatrix_E1 = projectionFactory.buildRHSMatrix(sourceGlobalIndexer,Teuchos::null,"E_Field",hcurlBD,true,1); // Project from vector basis
  auto rhsMatrix_E2 = projectionFactory.buildRHSMatrix(sourceGlobalIndexer,Teuchos::null,"E_Field",hcurlBD,true,2); // Project from vector basis

  // Build RHS Multivector target
  const auto targetRangeMap = massMatrix->getRangeMap();
  using NODE = Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>;
  const int numVectors = 6; // 3 E, 3 B
  const bool zeroOut = true;
  Tpetra::MultiVector<double,LO,GO,NODE> rhsMV(targetRangeMap,numVectors,zeroOut);
}
