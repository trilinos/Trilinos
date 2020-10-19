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
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
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
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"

//Intrepid2
#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

// Teuchos
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"

// Solver
#include "BelosSolverFactory.hpp"
#include "BelosMultiVecTraits_Tpetra.hpp"
#include "BelosOperatorTraits_Tpetra.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosLinearProblem.hpp"

#include <memory>
#include <vector>
#include <string>
#include <cmath>

TEUCHOS_UNIT_TEST(L2Projection, ToNodal)
{
  using namespace Teuchos;
  using namespace panzer;
  using namespace panzer_stk;

  RCP<MpiComm<int>> comm = rcp(new MpiComm<int>(MPI_COMM_WORLD));

  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("Total Time");

  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();
  const int numXElements = 3;
  const int numYElements = numXElements;
  const int numZElements = numXElements;
  const double boxLength = 1.0;
  TEUCHOS_ASSERT(numXElements >= numProcs);

  RCP<panzer_stk::STK_Interface> mesh;
  {
    PANZER_FUNC_TIME_MONITOR("L2Projection: mesh construction");
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("Z Blocks",1);
    pl->set("X Elements",numXElements);
    pl->set("Y Elements",numYElements);
    pl->set("Z Elements",numZElements);
    pl->set("X Procs",numProcs);
    pl->set("Y Procs",1);
    pl->set("Z Procs",1);
    pl->set("X0",0.0);
    pl->set("Y0",0.0);
    pl->set("Z0",0.0);
    pl->set("Xf",boxLength);
    pl->set("Yf",boxLength);
    pl->set("Zf",boxLength);
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
  worksetNeeds.addBasis(hdivBD);
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
  using GO = panzer::GlobalOrdinal;
  timer->start("ConnManager ctor");
  const RCP<panzer::ConnManager> connManager = rcp(new panzer_stk::STKConnManager(mesh));
  timer->stop("ConnManager ctor");

  // Set up bases for projections
  auto cellTopology = mesh->getCellTopology(eBlockNames[0]);
  auto dim = cellTopology->getDimension();

  auto hgradBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hgradBD.getType(),hgradBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hgradFP(new panzer::Intrepid2FieldPattern(hgradBasis));

  auto curlBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hcurlBD.getType(),hcurlBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hcurlFP(new panzer::Intrepid2FieldPattern(curlBasis));

  auto divBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hdivBD.getType(),hdivBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hdivFP(new panzer::Intrepid2FieldPattern(divBasis));

  // ****************
  // Testing ghosting
  // ****************
  // This shows that the new workset construction path does not
  // provide a full one-ring but only supplies ghosted elements on
  // faces of owned elements.
  const bool printGhostingInfo = true;
  if (printGhostingInfo) {
    PANZER_FUNC_TIME_MONITOR("L2Projection: test ghosting (build connectivity)");
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
  }
  // ****************
  // ****************

  // Build source DOF Manager that mimics multi-fluid plasma dof manager
  timer->start("Build sourceGlobalIndexer");
  RCP<panzer::DOFManager> sourceGlobalIndexer = rcp(new panzer::DOFManager(connManager,*comm->getRawMpiComm()));
  sourceGlobalIndexer->addField("PHI",hgradFP); // Electric Potential for ES
  sourceGlobalIndexer->addField("Chaff0",hgradFP); // Dummy
  sourceGlobalIndexer->addField("Chaff1",hgradFP); // Dummy
  sourceGlobalIndexer->addField("Chaff2",hgradFP); // Dummy
  sourceGlobalIndexer->addField("E_Field",hcurlFP); // Electric Field for EM
  sourceGlobalIndexer->addField("B_Field",hdivFP); // Magnetic Field for EM
  sourceGlobalIndexer->buildGlobalUnknowns();
  timer->stop("Build sourceGlobalIndexer");

  // Build Target DOF Manager (Separate scalar fields on hgrad)
  timer->start("Build targetGlobalIndexer");
  RCP<panzer::DOFManager> targetGlobalIndexer = rcp(new panzer::DOFManager(connManager,*comm->getRawMpiComm()));
  targetGlobalIndexer->addField("Projection to Mesh Vertices",hgradFP);
  targetGlobalIndexer->buildGlobalUnknowns();
  timer->stop("Build targetGlobalIndexer");

  // Build projection factory
  timer->start("projectionFactory.setup()");
  panzer::L2Projection projectionFactory;
  worksetContainer->setGlobalIndexer(sourceGlobalIndexer);
  projectionFactory.setup(hgradBD,integrationDescriptor,comm,connManager,eBlockNames,worksetContainer);
  timer->stop("projectionFactory.setup()");

  TEST_ASSERT(nonnull(projectionFactory.getTargetGlobalIndexer()));

  // Build mass matrix
  timer->start("projectionFactory.buildMassMatrix()");
  auto massMatrix = projectionFactory.buildMassMatrix();
  timer->stop("projectionFactory.buildMassMatrix()");
  massMatrix->print(out);
  massMatrix->getRowMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);
  massMatrix->getColMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);

  // Build rhs matrix
  timer->start("projectionFactory.buildRHSMatrix()");
  const int xDir = 0;
  const int yDir = 1;
  const int zDir = 2;
  using NodeType = Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>;
  auto rhsMatrix_PHI = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"PHI",hgradBD);          // Project value from scalar basis
  auto rhsMatrix_DPHI_DX = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"PHI",hgradBD,xDir); // Project gradient from scalar basis
  auto rhsMatrix_DPHI_DY = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"PHI",hgradBD,yDir); // Project gradient from scalar basis
  auto rhsMatrix_DPHI_DZ = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"PHI",hgradBD,zDir); // Project gradient from scalar basis
  auto rhsMatrix_E0 = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"E_Field",hcurlBD,xDir);  // Project value from vector basis
  auto rhsMatrix_E1 = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"E_Field",hcurlBD,yDir);  // Project value from vector basis
  auto rhsMatrix_E2 = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"E_Field",hcurlBD,zDir);  // Project value from vector basis
  auto rhsMatrix_B0 = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"B_Field",hdivBD,xDir);   // Project value from vector basis
  auto rhsMatrix_B1 = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"B_Field",hdivBD,yDir);   // Project value from vector basis
  auto rhsMatrix_B2 = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"B_Field",hdivBD,zDir);   // Project value from vector basis
  timer->stop("projectionFactory.buildRHSMatrix()");

  // Store in vector for convenience
  std::vector<RCP<Tpetra::CrsMatrix<double,LO,GO,NodeType>>> rhsMatrices;
  rhsMatrices.push_back(rhsMatrix_PHI);
  rhsMatrices.push_back(rhsMatrix_DPHI_DX);
  rhsMatrices.push_back(rhsMatrix_DPHI_DY);
  rhsMatrices.push_back(rhsMatrix_DPHI_DZ);
  rhsMatrices.push_back(rhsMatrix_E0);
  rhsMatrices.push_back(rhsMatrix_E1);
  rhsMatrices.push_back(rhsMatrix_E2);
  rhsMatrices.push_back(rhsMatrix_B0);
  rhsMatrices.push_back(rhsMatrix_B1);
  rhsMatrices.push_back(rhsMatrix_B2);

  // Create a Names vector for output
  std::vector<std::string> names;
  names.push_back("PHI");
  names.push_back("DPHI_DX");
  names.push_back("DPHI_DY");
  names.push_back("DPHI_DZ");
  names.push_back("E0");
  names.push_back("E1");
  names.push_back("E2");
  names.push_back("B0");
  names.push_back("B1");
  names.push_back("B2");

  // Allocate the source vector
  timer->start("Allocate Source Vector");
  using VectorType = Tpetra::Vector<double,LO,GO,NodeType>;
  auto sourceValues = rcp(new VectorType(rhsMatrix_PHI->getDomainMap(),true));
  timer->stop("Allocate Source Vector");

  // Fill the source vector.
  timer->start("Fill Source Vector");

  using DynRankView = Kokkos::DynRankView<double,PHX::Device>;
  using DynRankViewIntHost = Kokkos::DynRankView<int,Kokkos::HostSpace>;
  auto hostSourceValues = Kokkos::create_mirror_view(sourceValues->getLocalView<PHX::Device>());
  {
    const int PHI_Index = sourceGlobalIndexer->getFieldNum("PHI");
    const int E_Index = sourceGlobalIndexer->getFieldNum("E_Field");
    const int B_Index = sourceGlobalIndexer->getFieldNum("B_Field");

    std::vector<std::string> elementBlockNames;
    sourceGlobalIndexer->getElementBlockIds(elementBlockNames);
    for (const auto& block : elementBlockNames) {

      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
      const auto worksets = worksetContainer->getWorksets(wd);
      for (const auto& workset : *worksets) {


        //Initialize orientation tools
        DynRankViewIntHost ownedNodesGID("ownedNodesGID", workset.numOwnedCells(), cellTopology->getNodeCount());
        for (int e = 0; e < (int)workset.numOwnedCells(); ++e) {
          auto conn = connManager->getConnectivity(e);
          for(int j=0; j< (int)cellTopology->getNodeCount(); ++j)
            ownedNodesGID(e,j) = conn[j];
        }
        Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device> elemOrts("elemOrts", workset.numOwnedCells());
        Intrepid2::OrientationTools<PHX::Device>::getOrientation(elemOrts, ownedNodesGID, *cellTopology);

        Kokkos::View<LO**,PHX::Device> localIds("projection unit test: LocalIds", workset.numOwnedCells()+workset.numGhostCells()+workset.numVirtualCells(),
                                                sourceGlobalIndexer->getElementBlockGIDCount(block));
        // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
        const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
        sourceGlobalIndexer->getElementLIDs(cellLocalIdsNoGhost,localIds);

        // Create vector to store if LID is owned
        timer->start("Create isOwned view");
        Kokkos::View<bool**,PHX::Device> isOwned("projection unit test: isOwned", workset.numOwnedCells(),localIds.extent(1));
        {
          std::vector<GO> cellGIDs(localIds.extent(1));
          std::vector<bool> cellOwnedIds(localIds.extent(1));
          for (std::size_t cell=0; cell < cellLocalIdsNoGhost.extent(0); ++cell) {
            // Assumes UVM
            sourceGlobalIndexer->getElementGIDs(cellLocalIdsNoGhost(cell),cellGIDs);
            sourceGlobalIndexer->ownedIndices(cellGIDs,cellOwnedIds);
            for (std::size_t i=0; i < cellOwnedIds.size(); ++i)
              isOwned(cell,i) = cellOwnedIds[i];
          }
        }
        timer->stop("Create isOwned view");

        const auto offsetsPHI = sourceGlobalIndexer->getGIDFieldOffsetsKokkos(block,PHI_Index);
        const auto offsetsE = sourceGlobalIndexer->getGIDFieldOffsetsKokkos(block,E_Index);
        const auto offsetsB = sourceGlobalIndexer->getGIDFieldOffsetsKokkos(block,B_Index);
        const auto& basisValues = workset.getBasisValues(hgradBD,integrationDescriptor);
        const auto& coords = basisValues.basis_coordinates;
        const auto& x = sourceValues->getLocalView<PHX::Device>();
        const int numBasisPHI = static_cast<int>(offsetsPHI.extent(0));
        const int numBasisE = static_cast<int>(offsetsE.extent(0));
        const int numBasisB = static_cast<int>(offsetsB.extent(0));

        using li = Intrepid2::Experimental::LagrangianInterpolation<PHX::Device>;

        //Computing HGRAD coefficients for PHI to interpolate function f(x,y,z) = 1+x+2y+3z
        DynRankView basisCoeffsPHI("basisCoeffsPHI", workset.numOwnedCells(), numBasisPHI);
        {
          DynRankView dofCoordsPHI("dofCoordsPHI", workset.numOwnedCells(), numBasisPHI, dim);
          DynRankView dofCoeffsPHI("dofCoeffsPHI", workset.numOwnedCells(), numBasisPHI);

          li::getDofCoordsAndCoeffs(dofCoordsPHI,  dofCoeffsPHI, hgradBasis.getRawPtr(), elemOrts);

          //map the reference Dof coordinates into physical frame
          DynRankView physCoordsPHI("physCoordsPHI", workset.numOwnedCells(), numBasisPHI, dim);
          auto wsCoords = Kokkos::subview(coords.get_view(), std::pair<int,int>(0, workset.numOwnedCells()), Kokkos::ALL(), Kokkos::ALL());
          Intrepid2::CellTools<PHX::Device>::mapToPhysicalFrame(physCoordsPHI,dofCoordsPHI,wsCoords,*cellTopology);

          //evaluate the function f at the coordinates physCoordsPHI
          DynRankView functValuesAtDofCoordsPHI("funPHI", workset.numOwnedCells(), numBasisPHI);
          Kokkos::parallel_for(workset.numOwnedCells(),KOKKOS_LAMBDA (const int& cell) {
            for (int basis=0; basis < numBasisPHI; ++basis)
              functValuesAtDofCoordsPHI(cell, basis) = 1.0 + physCoordsPHI(cell,basis,0) + 2.0 * physCoordsPHI(cell,basis,1)
            + 3.0 * physCoordsPHI(cell,basis,2);
          });

          //compute basis coefficients
          li::getBasisCoeffs(basisCoeffsPHI, functValuesAtDofCoordsPHI, dofCoeffsPHI);
        }


        //Computing HCURL coefficients for E to interpolate the constant vector [1,1,1]
        DynRankView basisCoeffsE("basisCoeffsE", workset.numOwnedCells(), numBasisE);
        {
          DynRankView dofCoordsE("dofCoordsE", workset.numOwnedCells(), numBasisE, dim);
          DynRankView dofCoeffsE("dofCoeffsE", workset.numOwnedCells(), numBasisE, dim);
          li::getDofCoordsAndCoeffs(dofCoordsE,  dofCoeffsE, curlBasis.getRawPtr(), elemOrts);

          //Because the function is constant, we do not need to map the coordinates to physical frame.

          // Evaluate the function (in the physical frame) and map it back to the reference frame
          // In order to map an HCurl function back to the reference frame we need to multiply it
          // by J^T  (J being the Jacobian of the map from reference to physical frame)
          // In our case J is diagonal with entries equal to boxLength/2.0/numXElements,
          // where 2 is the length of the reference cell edge.
          DynRankView functValuesAtDofCoordsE("funE", workset.numOwnedCells(), numBasisE,dim);
          double curlScaling = boxLength/2.0/numXElements;
          Kokkos::parallel_for(workset.numOwnedCells(),KOKKOS_LAMBDA (const int& cell) {
            for (int basis=0; basis < numBasisE; ++basis)
              for (int d = 0; d < (int)dim; d++)
                functValuesAtDofCoordsE(cell,basis,d) = curlScaling*1.0;
          });

          //compute basis coefficients
          li::getBasisCoeffs(basisCoeffsE, functValuesAtDofCoordsE, dofCoeffsE);
        }

   /*
        // Alternative way of computing HCurl coefficients with L2 projection
        #include "Intrepid2_ProjectionTools.hpp"
        using pts = Intrepid2::Experimental::ProjectionTools<PHX::Device>;
        DynRankView basisCoeffsE("basisCoeffsE", workset.numOwnedCells(), numBasisE);
        {
          int targetCubDegree(0);
          Intrepid2::Experimental::ProjectionStruct<PHX::Device,double> projStruct;
          projStruct.createL2DGProjectionStruct(curlBasis, targetCubDegree);
          int numPoints = projStruct.getNumTargetEvalPoints();
          //DynRankView evalPoints("evalPoints", elemOrts, workset.numOwnedCells(), numPoints, dim);
          //pts::getL2EvaluationPoints(evalPoints, curlBasis, &projStruct);

          DynRankView functValuesAtEvalPoints("funE", workset.numOwnedCells(), numPoints ,dim);
          double curlScaling = boxLength/2.0/numXElements;
          Kokkos::parallel_for(workset.numOwnedCells(),KOKKOS_LAMBDA (const int& cell) {
            for (int pt=0; pt < numPoints; ++pt)
              for (int d = 0; d < (int)dim; d++)
                functValuesAtEvalPoints(cell,pt,d) = curlScaling*1.0;
          });
          pts::getL2DGBasisCoeffs(basisCoeffsE,
              functValuesAtEvalPoints,
              elemOrts,
              curlBasis.getRawPtr(),
              &projStruct);
        }
        */


        //Computing HDIV coefficients for B to interpolate the constant vector [1,1,1]
        DynRankView basisCoeffsB("basisCoeffsB", workset.numOwnedCells(), numBasisB);
        {
          DynRankView dofCoordsB("dofCoordsB", workset.numOwnedCells(), numBasisB, dim);
          DynRankView dofCoeffsB("dofCoeffsB", workset.numOwnedCells(), numBasisB, dim);
          li::getDofCoordsAndCoeffs(dofCoordsB,  dofCoeffsB, divBasis.getRawPtr(), elemOrts);


          // Evaluate the function (in the physical frame) and map it back to the reference frame
          // In order to map an HDiv function back to the reference frame we need to multiply it
          // by det(J) J^(-1)   (J being the Jacobian of the map from reference to physical frame)
          // In our case J is diagonal with entries equal to boxLength/2.0/numXElements,
          // where 2 is the length of the reference cell edge, so
          // det(J) J^(-1) = (boxLength/2.0/numXElements)^2 I
          DynRankView functValuesAtDofCoordsB("funB", workset.numOwnedCells(), numBasisB, dim);
          double divScaling = std::pow(boxLength/2.0/numXElements,2);
          Kokkos::parallel_for(workset.numOwnedCells(),KOKKOS_LAMBDA (const int& cell) {
            for (int basis=0; basis < numBasisB; ++basis)
              for (int d = 0; d < (int)dim; ++d)
                functValuesAtDofCoordsB(cell, basis,d) = divScaling * 1.0;
          });

          //compute basis coefficients
          li::getBasisCoeffs(basisCoeffsB, functValuesAtDofCoordsB, dofCoeffsB);
        }



        // fill the vector of basis coefficients x
        Kokkos::parallel_for(workset.numOwnedCells(),KOKKOS_LAMBDA (const int& cell) {
          for (int basis=0; basis < numBasisPHI; ++basis) {
            const int lid = localIds(cell,offsetsPHI(basis));
            if (isOwned(cell,offsetsPHI(basis))) {
              x(lid,0) = basisCoeffsPHI(cell,basis);
            }
          }
          for (int basis=0; basis < numBasisE; ++basis) {
            const int lid = localIds(cell,offsetsE(basis));
            if (isOwned(cell,offsetsE(basis)))
              x(lid,0) = basisCoeffsE(cell,basis);;
          }
          for (int basis=0; basis < numBasisB; ++basis) {
            const int lid = localIds(cell,offsetsB(basis));
            if (isOwned(cell,offsetsB(basis)))
              x(lid,0) = basisCoeffsB(cell,basis);
          }
        });
      }
    }
  }

  timer->stop("Fill Source Vector");

  // Assemble RHS vectors using matvec
  timer->start("Allocate RHS Vectors");
  auto rhs_PHI = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_DPHI_DX = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_DPHI_DY = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_DPHI_DZ = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_E0 = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_E1 = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_E2 = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_B0 = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_B1 = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_B2 = rcp(new VectorType(massMatrix->getDomainMap(),true));
  timer->stop("Allocate RHS Vectors");

  // Store in vector for convenience
  std::vector<RCP<VectorType>> rhs;
  rhs.push_back(rhs_PHI);
  rhs.push_back(rhs_DPHI_DX);
  rhs.push_back(rhs_DPHI_DY);
  rhs.push_back(rhs_DPHI_DZ);
  rhs.push_back(rhs_E0);
  rhs.push_back(rhs_E1);
  rhs.push_back(rhs_E2);
  rhs.push_back(rhs_B0);
  rhs.push_back(rhs_B1);
  rhs.push_back(rhs_B2);

  timer->start("Assemble RHS Vectors");
  for (size_t i=0; i < rhs.size(); ++i)
    rhsMatrices[i]->apply(*sourceValues,*rhs[i]);
  timer->stop("Assemble RHS Vectors");

  // Build RHS Multivector target for efficient consistent matrix solve
  timer->start("Copy RHS Values into MV");
  const auto targetRangeMap = massMatrix->getRangeMap();
  const int numVectors = 10; // PHI, DPHI_DX, DPHI_XY, DPHI_DZ, E, B
  const auto rhsMV = rcp(new Tpetra::MultiVector<double,LO,GO,NodeType>(targetRangeMap,numVectors,false));
  const auto mvView = rhsMV->getLocalView<NodeType>();
  TEST_EQUALITY(mvView.extent(1),static_cast<size_t>(numVectors));
  for (int col=0; col < numVectors; ++col) {
    const auto source = rhs[col]->getLocalView<NodeType>();
    const int numEntries = source.extent(0);
    Kokkos::parallel_for(numEntries, KOKKOS_LAMBDA (const int& i) { mvView(i,col) = source(i,0); });
    typename PHX::Device().fence();
  }
  timer->stop("Copy RHS Values into MV");

  // Solve the multiple matrices
  using MV = Tpetra::MultiVector<double,LO,GO,NodeType>;
  using OP = Tpetra::Operator<double,LO,GO,NodeType>;
  const auto solutionMV = rcp(new Tpetra::MultiVector<double,LO,GO,NodeType>(targetRangeMap,numVectors,true));
  const auto problem = rcp(new Belos::LinearProblem<double,MV,OP>(massMatrix, solutionMV, rhsMV));
  problem->setProblem();

  TEST_ASSERT(nonnull(massMatrix));
  TEST_ASSERT(nonnull(solutionMV));
  TEST_ASSERT(nonnull(rhsMV));
  TEST_ASSERT(nonnull(problem));

  timer->start("Create Belos Solver");
  RCP<ParameterList> pl = parameterList("L2 Consistent Projection Linear Solver");
  pl->set("Num Blocks", 100); // Max Krylov vectors to store
  pl->set("Maximum Iterations", 100);
  pl->set("Convergence Tolerance", 1e-8);
  pl->set("output Frequency",1); // all iterations
  pl->set("Verbosity", Belos::Errors | Belos::Warnings | Belos::IterationDetails| Belos::FinalSummary);
  Belos::SolverFactory<double,MV,OP> belosFactory;
  const auto solver = belosFactory.create("GMRES",pl);
  solver->setProblem(problem);
  timer->stop("Create Belos Solver");

  timer->start("Belos Solve");
  const auto result = solver->solve();
  timer->stop("Belos Solve");

  TEST_EQUALITY(result,Belos::Converged);

  // Check the final values on host
  timer->start("Check CONSISTENT Projected Values on Host");
  {
    const auto hostValues = Kokkos::create_mirror_view(solutionMV->getLocalView<PHX::Device>());
    Kokkos::deep_copy(hostValues,solutionMV->getLocalView<PHX::Device>());
    typename PHX::Device().fence();

    const int phiIndex = 0;
    const int dphiDxIndex = 1;
    const int dphiDyIndex = 2;
    const int dphiDzIndex = 3;
    const int exIndex = 4;
    const int eyIndex = 5;
    const int ezIndex = 6;
    const int bxIndex = 7;
    const int byIndex = 8;
    const int bzIndex = 9;
    double tol = 1.0e-4;

    std::vector<std::string> elementBlockNames;
    sourceGlobalIndexer->getElementBlockIds(elementBlockNames);
    for (const auto& block : elementBlockNames) {

      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
      const auto worksets = worksetContainer->getWorksets(wd);
      for (const auto& workset : *worksets) {
        Kokkos::View<LO**,PHX::Device> localIds("projection unit test: LocalIds", workset.numOwnedCells()+workset.numGhostCells()+workset.numVirtualCells(),
                                                targetGlobalIndexer->getElementBlockGIDCount(block));
        const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
        targetGlobalIndexer->getElementLIDs(cellLocalIdsNoGhost,localIds);

        // Create vector to store if LID is owned
        timer->start("Create isOwned view");
        Kokkos::View<bool**,PHX::Device> isOwned("projection unit test: isOwned", workset.numOwnedCells(),localIds.extent(1));
        {
          std::vector<GO> cellGIDs(localIds.extent(1));
          std::vector<bool> cellOwnedIds(localIds.extent(1));
          for (std::size_t cell=0; cell < cellLocalIdsNoGhost.extent(0); ++cell) {
            // Assumes UVM
            targetGlobalIndexer->getElementGIDs(cellLocalIdsNoGhost(cell),cellGIDs);
            targetGlobalIndexer->ownedIndices(cellGIDs,cellOwnedIds);
            for (std::size_t i=0; i < cellOwnedIds.size(); ++i)
              isOwned(cell,i) = cellOwnedIds[i];
          }
        }
        timer->stop("Create isOwned view");

        const auto offsets = targetGlobalIndexer->getGIDFieldOffsets(block,0);
        const auto& basisValues = workset.getBasisValues(hgradBD,integrationDescriptor);
        const auto& coords = basisValues.basis_coordinates;
        const int numBasis = static_cast<int>(offsets.size());

        for (int cell=0; cell < workset.numOwnedCells(); ++cell) {
          for (int basis=0; basis < numBasis; ++basis) {
            if (isOwned(cell,offsets[basis])) {
              const int lid = localIds(cell,offsets[basis]);
              const double phiGold = 1.0 + coords(cell,basis,0) + 2.0 * coords(cell,basis,1)
                + 3.0 * coords(cell,basis,2);

              //Note: here we are relying on the fact that, for low order HGRAD basis,
              //the basis coefficients are the values of the function at the nodes
              TEST_FLOATING_EQUALITY(hostValues(lid,phiIndex), phiGold, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDxIndex), 1.0, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDyIndex), 2.0, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDzIndex), 3.0, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,exIndex), 1.0, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,eyIndex), 1.0, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,ezIndex), 1.0, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,bxIndex), 1.0, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,byIndex), 1.0, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,bzIndex), 1.0, tol);
            }
          }
        }
      }
    }
  }
  timer->stop("Check CONSISTENT Projected Values on Host");

  // ***************************
  // Now test lumped mass matrix
  // ***************************
  timer->start("projectionFactory.buildInverseLumpedMassMatrix()");
  auto invLumpedMassMatrix = projectionFactory.buildInverseLumpedMassMatrix();
  timer->stop("projectionFactory.buildInverseLumpedMassMatrix()");

  timer->start("apply lumped mass matrix");
  {
    auto x = solutionMV->getLocalView<PHX::Device>();
    auto rhsMV_k = rhsMV->getLocalView<PHX::Device>();
    const auto ilmm = invLumpedMassMatrix->getLocalView<PHX::Device>();
    const int numEntries = static_cast<int>(x.extent(0));
    Kokkos::parallel_for(numEntries,KOKKOS_LAMBDA (const int i)
      {
        for (int field=0; field < numVectors; ++field)
          x(i,field) = ilmm(i,0) * rhsMV_k(i,field);
      });
    typename PHX::Device().fence();
    solutionMV->template modify<PHX::Device>();
  }
  timer->stop("apply lumped mass matrix");

  // Check the final values on host
  timer->start("Check LUMPED Projected Values on Host");
  {
    const auto hostValues = Kokkos::create_mirror_view(solutionMV->getLocalView<PHX::Device>());
    Kokkos::deep_copy(hostValues,solutionMV->getLocalView<PHX::Device>());
    typename PHX::Device().fence();

    const int phiIndex = 0;
    const int dphiDxIndex = 1;
    const int dphiDyIndex = 2;
    const int dphiDzIndex = 3;
    const int exIndex = 4;
    const int eyIndex = 5;
    const int ezIndex = 6;
    const int bxIndex = 7;
    const int byIndex = 8;
    const int bzIndex = 9;
    // NOTE: Lumping generates significant errors when the meshes are
    // very coarse (such as in this unit test). We see
    // superconvergence on the interior points of this axis aligned
    // mesh and for the constant fields, thus can can tighten the
    // tolerances for some fields. Refining the mesh results in more
    // accurate projected values.
    double looseTol = 0.5;
    double superconvergedTol = 1.0e-4;

    std::vector<std::string> elementBlockNames;
    sourceGlobalIndexer->getElementBlockIds(elementBlockNames);
    for (const auto& block : elementBlockNames) {

      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
      const auto worksets = worksetContainer->getWorksets(wd);
      for (const auto& workset : *worksets) {
        Kokkos::View<LO**,PHX::Device> localIds("projection unit test: LocalIds", workset.numOwnedCells()+workset.numGhostCells()+workset.numVirtualCells(),
                                                targetGlobalIndexer->getElementBlockGIDCount(block));
        const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
        targetGlobalIndexer->getElementLIDs(cellLocalIdsNoGhost,localIds);

        // Create vector to store if LID is owned
        timer->start("Create isOwned view");
        Kokkos::View<bool**,PHX::Device> isOwned("projection unit test: isOwned", workset.numOwnedCells(),localIds.extent(1));
        {
          std::vector<GO> cellGIDs(localIds.extent(1));
          std::vector<bool> cellOwnedIds(localIds.extent(1));
          for (std::size_t cell=0; cell < cellLocalIdsNoGhost.extent(0); ++cell) {
            // Assumes UVM
            targetGlobalIndexer->getElementGIDs(cellLocalIdsNoGhost(cell),cellGIDs);
            targetGlobalIndexer->ownedIndices(cellGIDs,cellOwnedIds);
            for (std::size_t i=0; i < cellOwnedIds.size(); ++i)
              isOwned(cell,i) = cellOwnedIds[i];
          }
        }
        timer->stop("Create isOwned view");

        const auto offsets = targetGlobalIndexer->getGIDFieldOffsets(block,0);
        const auto& basisValues = workset.getBasisValues(hgradBD,integrationDescriptor);
        const auto& coords = basisValues.basis_coordinates;
        const int numBasis = static_cast<int>(offsets.size());

        for (int cell=0; cell < workset.numOwnedCells(); ++cell) {
          for (int basis=0; basis < numBasis; ++basis) {
            if (isOwned(cell,offsets[basis])) {
              const int lid = localIds(cell,offsets[basis]);
              const double phiGold = 1.0 + coords(cell,basis,0) + 2.0 * coords(cell,basis,1)
                + 3.0 * coords(cell,basis,2);

              //Note: here we are relying on the fact that, for low order HGRAD basis,
              //the basis coefficients are the values of the function at the nodes
              TEST_FLOATING_EQUALITY(hostValues(lid,phiIndex), phiGold, looseTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDxIndex), 1.0, superconvergedTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDyIndex), 2.0, superconvergedTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDzIndex), 3.0, superconvergedTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,exIndex), 1.0, superconvergedTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,eyIndex), 1.0, superconvergedTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,ezIndex), 1.0, superconvergedTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,bxIndex), 1.0, superconvergedTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,byIndex), 1.0, superconvergedTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,bzIndex), 1.0, superconvergedTol);
            }
          }
        }
      }
    }
  }
  timer->stop("Check LUMPED Projected Values on Host");

  timer->stop("Total Time");

  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction = true;
  options.output_minmax = true;
  options.output_histogram = false;
  options.num_histogram = 5;
  timer->report(out,comm,options);
}

TEUCHOS_UNIT_TEST(L2Projection, CurlMassMatrix)
{
  using namespace Teuchos;
  using namespace panzer;
  using namespace panzer_stk;

  RCP<MpiComm<int>> comm = rcp(new MpiComm<int>(MPI_COMM_WORLD));

  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("Total Time");

  const int numProcs = comm->getSize();
  const int numXElements = 4;
  const int numYElements = 2;
  const double boxLength = 1.0;
  const double boxHeight = 0.5;
  TEUCHOS_ASSERT(numXElements >= numProcs);

  RCP<panzer_stk::STK_Interface> mesh;
  {
    PANZER_FUNC_TIME_MONITOR("L2Projection: mesh construction");
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",numXElements);
    pl->set("Y Elements",numYElements);
    pl->set("X Procs",numProcs);
    pl->set("Y Procs",1);
    pl->set("X0",0.0);
    pl->set("Y0",0.0);
    pl->set("Xf",boxLength);
    pl->set("Yf",boxHeight);
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  }

  // Build Worksets

  const int basisOrder = 1;
  BasisDescriptor hcurlBD(basisOrder,"HCurl");
  BasisDescriptor hdivBD(basisOrder,"HDiv");

  const int intOrder = 2;
  IntegrationDescriptor integrationDescriptor(intOrder,IntegrationDescriptor::VOLUME);

  WorksetNeeds worksetNeeds;
  worksetNeeds.addBasis(hcurlBD);
  worksetNeeds.addBasis(hdivBD);
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
  using GO = panzer::GlobalOrdinal;
  timer->start("ConnManager ctor");
  const RCP<panzer::ConnManager> connManager = rcp(new panzer_stk::STKConnManager(mesh));
  timer->stop("ConnManager ctor");

  // Set up bases for projections
  auto cellTopology = mesh->getCellTopology(eBlockNames[0]);

  auto curlBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hcurlBD.getType(),hcurlBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hcurlFP(new panzer::Intrepid2FieldPattern(curlBasis));

  auto divBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hdivBD.getType(),hdivBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hdivFP(new panzer::Intrepid2FieldPattern(divBasis));

  // Build source DOF Manager for edge and face DOFs
  timer->start("Build sourceGlobalIndexer");
  RCP<panzer::BlockedDOFManager> sourceGlobalIndexer = rcp(new panzer::BlockedDOFManager(connManager,*comm->getRawMpiComm()));
  sourceGlobalIndexer->addField("E_Field",hcurlFP); // Electric Field for EM
  sourceGlobalIndexer->addField("B_Field",hdivFP); // Magnetic Field for EM
  sourceGlobalIndexer->buildGlobalUnknowns();
  timer->stop("Build sourceGlobalIndexer");

  // Build projection factory
  timer->start("projectionFactory.setup()");
  panzer::L2Projection projectionFactory;
  worksetContainer->setGlobalIndexer(sourceGlobalIndexer);
  projectionFactory.setup(hcurlBD,integrationDescriptor,comm,connManager,eBlockNames,worksetContainer);
  timer->stop("projectionFactory.setup()");

  TEST_ASSERT(nonnull(projectionFactory.getTargetGlobalIndexer()));

  // Build mass matrix
  timer->start("projectionFactory.buildMassMatrix()");
  auto curlMassMatrix = projectionFactory.buildMassMatrix();
  timer->stop("projectionFactory.buildMassMatrix()");
  curlMassMatrix->print(out);
  curlMassMatrix->getRowMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);
  curlMassMatrix->getColMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);

  // Build the mass matrix from connectivity knowing that this is a regular quad mesh
  timer->start("build mass matrix from connectivity");

  // get local ids
  auto e_ugi = sourceGlobalIndexer->getFieldDOFManagers()[sourceGlobalIndexer->getFieldBlock(sourceGlobalIndexer->getFieldNum("E_Field"))];
  auto lids = e_ugi->getLIDs();

  // set up a global and ghosted mass matrix
  std::vector<Teuchos::RCP<const panzer::GlobalIndexer>> indexers;
  indexers.push_back(e_ugi);
  panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,LO,GO,panzer::TpetraNodeType> factory(comm,indexers);
  auto connMassMatrix = factory.getTpetraMatrix(0,0);
  auto ghostedMatrix = factory.getGhostedTpetraMatrix(0,0);
  connMassMatrix->resumeFill();
  connMassMatrix->setAllToScalar(0.0);
  ghostedMatrix->resumeFill();
  ghostedMatrix->setAllToScalar(0.0);
  typename PHX::Device().fence();

  // fill in the mass matrix
  // the integral of the edge basis squared over one cell is 4/3
  // the integral of the edge basis times the basis function across from it in the element is 2/3
  const auto localMass = ghostedMatrix->getLocalMatrix();
  const int numElems = lids.extent(0);
  Kokkos::parallel_for(numElems, KOKKOS_LAMBDA (const int& i) {
    double row_values[2]={4.0/3.0,2.0/3.0};
    LO cols[2];
    for(int r = 0; r < 4; r++){
      cols[0] = lids(i,r);
      cols[1] = lids(i,(r+2)%4);
      localMass.sumIntoValues(lids(i,r),cols,2,row_values,false,true);
    }
  });
  typename PHX::Device().fence();
  ghostedMatrix->fillComplete();
  const auto exporter = factory.getGhostedExport(0);
  connMassMatrix->doExport(*ghostedMatrix, *exporter, Tpetra::ADD);
  connMassMatrix->fillComplete();
  timer->stop("build mass matrix from connectivity");

  connMassMatrix->print(out);
  connMassMatrix->getRowMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);
  connMassMatrix->getColMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);

  // compute difference between the two versions of the mass matrix
  using NodeType = Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>;
  auto difference = Tpetra::MatrixMatrix::add<double,LO,GO,NodeType>(1.0,false,*curlMassMatrix,-1.0,false,*connMassMatrix);
  double error = difference->getFrobeniusNorm();
  double norm = connMassMatrix->getFrobeniusNorm();
  double tol = 1.0e-14;
  TEST_COMPARE(error,<,tol*norm);
}

// This is to demonstrate the issue with row sum lumping as opposed to
// proportional lumping.
TEUCHOS_UNIT_TEST(L2Projection, HighOrderTri)
{
  using namespace Teuchos;
  using namespace panzer;
  using namespace panzer_stk;

  RCP<MpiComm<int>> comm = rcp(new MpiComm<int>(MPI_COMM_WORLD));

  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("Total Time");

  const int numProcs = comm->getSize();
  const int numXElements = 3;
  const int numYElements = numXElements;
  const double boxLength = 1.0;
  TEUCHOS_ASSERT(numXElements >= numProcs);

  RCP<panzer_stk::STK_Interface> mesh;
  {
    PANZER_FUNC_TIME_MONITOR("L2Projection: mesh construction");
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",numXElements);
    pl->set("Y Elements",numYElements);
    pl->set("X Procs",numProcs);
    pl->set("Y Procs",1);
    pl->set("X0",0.0);
    pl->set("Y0",0.0);
    pl->set("Xf",boxLength);
    pl->set("Yf",boxLength);
    panzer_stk::SquareTriMeshFactory factory;
    factory.setParameterList(pl);
    mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  }

  // Build Worksets

  const int basisOrder = 1;
  BasisDescriptor hgradBD(basisOrder,"HGrad");

  const int intOrder = 2;
  IntegrationDescriptor integrationDescriptor(intOrder,IntegrationDescriptor::VOLUME);

  WorksetNeeds worksetNeeds;
  worksetNeeds.addBasis(hgradBD);
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
  using GO = panzer::GlobalOrdinal;
  timer->start("ConnManager ctor");
  const RCP<panzer::ConnManager> connManager = rcp(new panzer_stk::STKConnManager(mesh));
  timer->stop("ConnManager ctor");

  // Set up bases for projections
  auto cellTopology = mesh->getCellTopology(eBlockNames[0]);

  auto hgradBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hgradBD.getType(),hgradBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hgradFP(new panzer::Intrepid2FieldPattern(hgradBasis));

  // Build source DOF Manager that mimics multi-fluid plasma dof manager
  timer->start("Build sourceGlobalIndexer");
  RCP<panzer::DOFManager> sourceGlobalIndexer = rcp(new panzer::DOFManager(connManager,*comm->getRawMpiComm()));
  sourceGlobalIndexer->addField("PHI",hgradFP); // Electric Potential for ES
  sourceGlobalIndexer->addField("Chaff0",hgradFP); // Dummy
  sourceGlobalIndexer->addField("Chaff1",hgradFP); // Dummy
  sourceGlobalIndexer->addField("Chaff2",hgradFP); // Dummy
  sourceGlobalIndexer->buildGlobalUnknowns();
  timer->stop("Build sourceGlobalIndexer");

  // Build Target DOF Manager (Separate scalar fields on hgrad)
  timer->start("Build targetGlobalIndexer");
  RCP<panzer::DOFManager> targetGlobalIndexer = rcp(new panzer::DOFManager(connManager,*comm->getRawMpiComm()));
  targetGlobalIndexer->addField("Projection to Mesh Vertices",hgradFP);
  targetGlobalIndexer->buildGlobalUnknowns();
  timer->stop("Build targetGlobalIndexer");

  // Build projection factory
  timer->start("projectionFactory.setup()");
  panzer::L2Projection projectionFactory;
  worksetContainer->setGlobalIndexer(sourceGlobalIndexer);
  projectionFactory.setup(hgradBD,integrationDescriptor,comm,connManager,eBlockNames,worksetContainer);
  timer->stop("projectionFactory.setup()");

  TEST_ASSERT(nonnull(projectionFactory.getTargetGlobalIndexer()));

  // Build mass matrix
  timer->start("projectionFactory.buildMassMatrix()");
  auto massMatrix = projectionFactory.buildMassMatrix();
  timer->stop("projectionFactory.buildMassMatrix()");
  massMatrix->print(out);
  massMatrix->getRowMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);
  massMatrix->getColMap()->describe(out,Teuchos::EVerbosityLevel::VERB_EXTREME);

  // Build rhs matrix
  timer->start("projectionFactory.buildRHSMatrix()");
  const int xDir = 0;
  const int yDir = 1;
  using NodeType = Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>;
  auto rhsMatrix_PHI = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"PHI",hgradBD);          // Project value from scalar basis
  auto rhsMatrix_DPHI_DX = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"PHI",hgradBD,xDir); // Project gradient from scalar basis
  auto rhsMatrix_DPHI_DY = projectionFactory.buildRHSMatrix(*sourceGlobalIndexer,Teuchos::null,"PHI",hgradBD,yDir); // Project gradient from scalar basis
  timer->stop("projectionFactory.buildRHSMatrix()");

  // Store in vector for convenience
  std::vector<RCP<Tpetra::CrsMatrix<double,LO,GO,NodeType>>> rhsMatrices;
  rhsMatrices.push_back(rhsMatrix_PHI);
  rhsMatrices.push_back(rhsMatrix_DPHI_DX);
  rhsMatrices.push_back(rhsMatrix_DPHI_DY);

  // Create a Names vector for output
  std::vector<std::string> names;
  names.push_back("PHI");
  names.push_back("DPHI_DX");
  names.push_back("DPHI_DY");

  // Allocate the source vector
  timer->start("Allocate Source Vector");
  using VectorType = Tpetra::Vector<double,LO,GO,NodeType>;
  auto sourceValues = rcp(new VectorType(rhsMatrix_PHI->getDomainMap(),true));
  timer->stop("Allocate Source Vector");

  // Fill the source vector.
  timer->start("Fill Source Vector");
  auto hostSourceValues = Kokkos::create_mirror_view(sourceValues->getLocalView<PHX::Device>());
  {
    const int PHI_Index = sourceGlobalIndexer->getFieldNum("PHI");

    std::vector<std::string> elementBlockNames;
    sourceGlobalIndexer->getElementBlockIds(elementBlockNames);
    for (const auto& block : elementBlockNames) {

      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
      const auto worksets = worksetContainer->getWorksets(wd);
      for (const auto& workset : *worksets) {

        Kokkos::View<LO**,PHX::Device> localIds("projection unit test: LocalIds", workset.numOwnedCells()+workset.numGhostCells()+workset.numVirtualCells(),
                                                sourceGlobalIndexer->getElementBlockGIDCount(block));
        // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
        const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
        sourceGlobalIndexer->getElementLIDs(cellLocalIdsNoGhost,localIds);

        // Create vector to store if LID is owned
        timer->start("Create isOwned view");
        Kokkos::View<bool**,PHX::Device> isOwned("projection unit test: isOwned", workset.numOwnedCells(),localIds.extent(1));
        {
          std::vector<GO> cellGIDs(localIds.extent(1));
          std::vector<bool> cellOwnedIds(localIds.extent(1));
          for (std::size_t cell=0; cell < cellLocalIdsNoGhost.extent(0); ++cell) {
            // Assumes UVM
            sourceGlobalIndexer->getElementGIDs(cellLocalIdsNoGhost(cell),cellGIDs);
            sourceGlobalIndexer->ownedIndices(cellGIDs,cellOwnedIds);
            for (std::size_t i=0; i < cellOwnedIds.size(); ++i)
              isOwned(cell,i) = cellOwnedIds[i];
          }
        }
        timer->stop("Create isOwned view");

        const auto offsetsPHI = sourceGlobalIndexer->getGIDFieldOffsetsKokkos(block,PHI_Index);
        const auto& basisValues = workset.getBasisValues(hgradBD,integrationDescriptor);
        const auto& coords = basisValues.basis_coordinates;
        const auto& x = sourceValues->getLocalView<PHX::Device>();
        const int numBasisPHI = static_cast<int>(offsetsPHI.extent(0));

        Kokkos::parallel_for(workset.numOwnedCells(),KOKKOS_LAMBDA (const int& cell) {
          for (int basis=0; basis < numBasisPHI; ++basis) {
            const int lid = localIds(cell,offsetsPHI(basis));
            if (isOwned(cell,offsetsPHI(basis)))
              x(lid,0) = 1.0 + coords(cell,basis,0) + 2.0 * coords(cell,basis,1);
          }
        });
      }
    }
  }

  timer->stop("Fill Source Vector");

  // Assemble RHS vectors using matvec
  timer->start("Allocate RHS Vectors");
  auto rhs_PHI = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_DPHI_DX = rcp(new VectorType(massMatrix->getDomainMap(),true));
  auto rhs_DPHI_DY = rcp(new VectorType(massMatrix->getDomainMap(),true));
  timer->stop("Allocate RHS Vectors");

  // Store in vector for convenience
  std::vector<RCP<VectorType>> rhs;
  rhs.push_back(rhs_PHI);
  rhs.push_back(rhs_DPHI_DX);
  rhs.push_back(rhs_DPHI_DY);

  timer->start("Assemble RHS Vectors");
  for (size_t i=0; i < rhs.size(); ++i)
    rhsMatrices[i]->apply(*sourceValues,*rhs[i]);
  timer->stop("Assemble RHS Vectors");

  // Build RHS Multivector target for efficient consistent matrix solve
  timer->start("Copy RHS Values into MV");
  const auto targetRangeMap = massMatrix->getRangeMap();
  const int numVectors = 3; // PHI, DPHI_DX, DPHI_XY
  const auto rhsMV = rcp(new Tpetra::MultiVector<double,LO,GO,NodeType>(targetRangeMap,numVectors,false));
  const auto mvView = rhsMV->getLocalView<NodeType>();
  TEST_EQUALITY(mvView.extent(1),static_cast<size_t>(numVectors));
  for (int col=0; col < numVectors; ++col) {
    const auto source = rhs[col]->getLocalView<NodeType>();
    const int numEntries = source.extent(0);
    Kokkos::parallel_for(numEntries, KOKKOS_LAMBDA (const int& i) { mvView(i,col) = source(i,0); });
    typename PHX::Device().fence();
  }
  timer->stop("Copy RHS Values into MV");

  // Solve the multiple matrices
  using MV = Tpetra::MultiVector<double,LO,GO,NodeType>;
  using OP = Tpetra::Operator<double,LO,GO,NodeType>;
  const auto solutionMV = rcp(new Tpetra::MultiVector<double,LO,GO,NodeType>(targetRangeMap,numVectors,true));
  const auto problem = rcp(new Belos::LinearProblem<double,MV,OP>(massMatrix, solutionMV, rhsMV));
  problem->setProblem();

  TEST_ASSERT(nonnull(massMatrix));
  TEST_ASSERT(nonnull(solutionMV));
  TEST_ASSERT(nonnull(rhsMV));
  TEST_ASSERT(nonnull(problem));

  timer->start("Create Belos Solver");
  RCP<ParameterList> pl = parameterList("L2 Consistent Projection Linear Solver");
  pl->set("Num Blocks", 100); // Max Krylov vectors to store
  pl->set("Maximum Iterations", 100);
  pl->set("Convergence Tolerance", 1e-8);
  pl->set("output Frequency",1); // all iterations
  pl->set("Verbosity", Belos::Errors | Belos::Warnings | Belos::IterationDetails| Belos::FinalSummary);
  Belos::SolverFactory<double,MV,OP> belosFactory;
  const auto solver = belosFactory.create("GMRES",pl);
  solver->setProblem(problem);
  timer->stop("Create Belos Solver");

  timer->start("Belos Solve");
  const auto result = solver->solve();
  timer->stop("Belos Solve");

  TEST_EQUALITY(result,Belos::Converged);

  // Check the final values on host
  out << "Checking CONSISTENT values!" << std::endl;
  timer->start("Check CONSISTENT Projected Values on Host");
  {
    const auto hostValues = Kokkos::create_mirror_view(solutionMV->getLocalView<PHX::Device>());
    Kokkos::deep_copy(hostValues,solutionMV->getLocalView<PHX::Device>());
    typename PHX::Device().fence();

    const int phiIndex = 0;
    const int dphiDxIndex = 1;
    const int dphiDyIndex = 2;
    double tol = 1.0e-4;

    std::vector<std::string> elementBlockNames;
    sourceGlobalIndexer->getElementBlockIds(elementBlockNames);
    for (const auto& block : elementBlockNames) {

      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
      const auto worksets = worksetContainer->getWorksets(wd);
      for (const auto& workset : *worksets) {
        Kokkos::View<LO**,PHX::Device> localIds("projection unit test: LocalIds", workset.numOwnedCells()+workset.numGhostCells()+workset.numVirtualCells(),
                                                targetGlobalIndexer->getElementBlockGIDCount(block));
        const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
        targetGlobalIndexer->getElementLIDs(cellLocalIdsNoGhost,localIds);

        // Create vector to store if LID is owned
        timer->start("Create isOwned view");
        Kokkos::View<bool**,PHX::Device> isOwned("projection unit test: isOwned", workset.numOwnedCells(),localIds.extent(1));
        {
          std::vector<GO> cellGIDs(localIds.extent(1));
          std::vector<bool> cellOwnedIds(localIds.extent(1));
          for (std::size_t cell=0; cell < cellLocalIdsNoGhost.extent(0); ++cell) {
            // Assumes UVM
            targetGlobalIndexer->getElementGIDs(cellLocalIdsNoGhost(cell),cellGIDs);
            targetGlobalIndexer->ownedIndices(cellGIDs,cellOwnedIds);
            for (std::size_t i=0; i < cellOwnedIds.size(); ++i)
              isOwned(cell,i) = cellOwnedIds[i];
          }
        }
        timer->stop("Create isOwned view");

        const auto offsets = targetGlobalIndexer->getGIDFieldOffsets(block,0);
        const auto& basisValues = workset.getBasisValues(hgradBD,integrationDescriptor);
        const auto& coords = basisValues.basis_coordinates;
        const int numBasis = static_cast<int>(offsets.size());

        for (int cell=0; cell < workset.numOwnedCells(); ++cell) {
          for (int basis=0; basis < numBasis; ++basis) {
            if (isOwned(cell,offsets[basis])) {
              const int lid = localIds(cell,offsets[basis]);
              const double phiGold = 1.0 + coords(cell,basis,0) + 2.0 * coords(cell,basis,1);
              TEST_FLOATING_EQUALITY(hostValues(lid,phiIndex), phiGold, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDxIndex), 1.0, tol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDyIndex), 2.0, tol);
            }
          }
        }
      }
    }
  }
  timer->stop("Check CONSISTENT Projected Values on Host");

  // ***************************
  // Now test lumped mass matrix
  // ***************************
  timer->start("projectionFactory.buildInverseLumpedMassMatrix()");
  auto invLumpedMassMatrix = projectionFactory.buildInverseLumpedMassMatrix();
  timer->stop("projectionFactory.buildInverseLumpedMassMatrix()");

  timer->start("apply lumped mass matrix");
  {
    auto x = solutionMV->getLocalView<PHX::Device>();
    auto rhsMV_k = rhsMV->getLocalView<PHX::Device>();
    const auto ilmm = invLumpedMassMatrix->getLocalView<PHX::Device>();
    const int numEntries = static_cast<int>(x.extent(0));
    Kokkos::parallel_for(numEntries,KOKKOS_LAMBDA (const int i)
      {
        for (int field=0; field < numVectors; ++field)
          x(i,field) = ilmm(i,0) * rhsMV_k(i,field);
      });
    typename PHX::Device().fence();
    solutionMV->template modify<PHX::Device>();
  }
  timer->stop("apply lumped mass matrix");

  // Check the final values on host
  out << "Checking LUMPED values!" << std::endl;
  timer->start("Check LUMPED Projected Values on Host");
  {
    const auto hostValues = Kokkos::create_mirror_view(solutionMV->getLocalView<PHX::Device>());
    Kokkos::deep_copy(hostValues,solutionMV->getLocalView<PHX::Device>());
    typename PHX::Device().fence();

    const int phiIndex = 0;
    const int dphiDxIndex = 1;
    const int dphiDyIndex = 2;
    // NOTE: Lumping generates significant errors when the meshes are
    // very coarse (such as in this unit test). We see
    // superconvergence on the interior points of this axis aligned
    // mesh and for the constant fields, thus can can tighten the
    // tolerances for some fields. Refining the mesh results in more
    // accurate projected values.
    double looseTol = 0.5;
    double superconvergedTol = 1.0e-4;

    std::vector<std::string> elementBlockNames;
    sourceGlobalIndexer->getElementBlockIds(elementBlockNames);
    for (const auto& block : elementBlockNames) {

      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
      const auto worksets = worksetContainer->getWorksets(wd);
      for (const auto& workset : *worksets) {
        Kokkos::View<LO**,PHX::Device> localIds("projection unit test: LocalIds", workset.numOwnedCells()+workset.numGhostCells()+workset.numVirtualCells(),
                                                targetGlobalIndexer->getElementBlockGIDCount(block));
        const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
        targetGlobalIndexer->getElementLIDs(cellLocalIdsNoGhost,localIds);

        // Create vector to store if LID is owned
        timer->start("Create isOwned view");
        Kokkos::View<bool**,PHX::Device> isOwned("projection unit test: isOwned", workset.numOwnedCells(),localIds.extent(1));
        {
          std::vector<GO> cellGIDs(localIds.extent(1));
          std::vector<bool> cellOwnedIds(localIds.extent(1));
          for (std::size_t cell=0; cell < cellLocalIdsNoGhost.extent(0); ++cell) {
            // Assumes UVM
            targetGlobalIndexer->getElementGIDs(cellLocalIdsNoGhost(cell),cellGIDs);
            targetGlobalIndexer->ownedIndices(cellGIDs,cellOwnedIds);
            for (std::size_t i=0; i < cellOwnedIds.size(); ++i)
              isOwned(cell,i) = cellOwnedIds[i];
          }
        }
        timer->stop("Create isOwned view");

        const auto offsets = targetGlobalIndexer->getGIDFieldOffsets(block,0);
        const auto& basisValues = workset.getBasisValues(hgradBD,integrationDescriptor);
        const auto& coords = basisValues.basis_coordinates;
        const int numBasis = static_cast<int>(offsets.size());

        for (int cell=0; cell < workset.numOwnedCells(); ++cell) {
          for (int basis=0; basis < numBasis; ++basis) {
            if (isOwned(cell,offsets[basis])) {
              const int lid = localIds(cell,offsets[basis]);
              const double phiGold = 1.0 + coords(cell,basis,0) + 2.0 * coords(cell,basis,1);
              TEST_FLOATING_EQUALITY(hostValues(lid,phiIndex), phiGold, looseTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDxIndex), 1.0, superconvergedTol);
              TEST_FLOATING_EQUALITY(hostValues(lid,dphiDyIndex), 2.0, superconvergedTol);
            }
          }
        }
      }
    }
  }
  timer->stop("Check LUMPED Projected Values on Host");

  timer->stop("Total Time");

  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction = true;
  options.output_minmax = true;
  options.output_histogram = false;
  options.num_histogram = 5;
  timer->report(out,comm,options);
}


// Tests the element block multiplier
TEUCHOS_UNIT_TEST(L2Projection, ElementBlockMultiplier)
{
  using namespace Teuchos;
  using namespace panzer;
  using namespace panzer_stk;

  RCP<MpiComm<int>> comm = rcp(new MpiComm<int>(MPI_COMM_WORLD));

  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("Total Time");

  const int numProcs = comm->getSize();
  const int numXElements = 3;
  const int numYElements = numXElements;
  const double boxLength = 1.0;
  TEUCHOS_ASSERT(numXElements >= numProcs);

  RCP<panzer_stk::STK_Interface> mesh;
  {
    PANZER_FUNC_TIME_MONITOR("L2Projection: mesh construction");
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",numXElements);
    pl->set("Y Elements",numYElements);
    pl->set("X Procs",numProcs);
    pl->set("Y Procs",1);
    pl->set("X0",0.0);
    pl->set("Y0",0.0);
    pl->set("Xf",boxLength);
    pl->set("Yf",boxLength);
    panzer_stk::SquareTriMeshFactory factory;
    factory.setParameterList(pl);
    mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  }

  // Build Worksets

  const int basisOrder = 1;
  BasisDescriptor hgradBD(basisOrder,"HGrad");

  const int intOrder = 2;
  IntegrationDescriptor integrationDescriptor(intOrder,IntegrationDescriptor::VOLUME);

  WorksetNeeds worksetNeeds;
  worksetNeeds.addBasis(hgradBD);
  worksetNeeds.addIntegrator(integrationDescriptor);

  RCP<WorksetFactory> worksetFactory(new WorksetFactory(mesh));
  std::vector<std::string> eBlockNames;
  mesh->getElementBlockNames(eBlockNames);
  std::map<std::string,WorksetNeeds> eblockNeeds;
  for (const auto& block : eBlockNames)
    eblockNeeds[block] = worksetNeeds;
  RCP<WorksetContainer> worksetContainer(new WorksetContainer(worksetFactory,eblockNeeds));

  // Build Connection Manager
  timer->start("ConnManager ctor");
  const RCP<panzer::ConnManager> connManager = rcp(new panzer_stk::STKConnManager(mesh));
  timer->stop("ConnManager ctor");

  // Set up bases for projections
  auto cellTopology = mesh->getCellTopology(eBlockNames[0]);

  auto hgradBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(hgradBD.getType(),hgradBD.getOrder(),*cellTopology);
  RCP<const panzer::FieldPattern> hgradFP(new panzer::Intrepid2FieldPattern(hgradBasis));

  // Build source DOF Manager that mimics multi-fluid plasma dof manager
  timer->start("Build sourceGlobalIndexer");
  RCP<panzer::DOFManager> sourceGlobalIndexer = rcp(new panzer::DOFManager(connManager,*comm->getRawMpiComm()));
  sourceGlobalIndexer->addField("PHI",hgradFP); // Electric Potential for ES
  sourceGlobalIndexer->addField("Chaff0",hgradFP); // Dummy
  sourceGlobalIndexer->addField("Chaff1",hgradFP); // Dummy
  sourceGlobalIndexer->addField("Chaff2",hgradFP); // Dummy
  sourceGlobalIndexer->buildGlobalUnknowns();
  timer->stop("Build sourceGlobalIndexer");

  // Build Target DOF Manager (Separate scalar fields on hgrad)
  timer->start("Build targetGlobalIndexer");
  RCP<panzer::DOFManager> targetGlobalIndexer = rcp(new panzer::DOFManager(connManager,*comm->getRawMpiComm()));
  targetGlobalIndexer->addField("Projection to Mesh Vertices",hgradFP);
  targetGlobalIndexer->buildGlobalUnknowns();
  timer->stop("Build targetGlobalIndexer");

  // Build projection factory
  timer->start("projectionFactory.setup()");
  panzer::L2Projection projectionFactory;
  worksetContainer->setGlobalIndexer(sourceGlobalIndexer);
  projectionFactory.setup(hgradBD,integrationDescriptor,comm,connManager,eBlockNames,worksetContainer);
  timer->stop("projectionFactory.setup()");

  TEST_ASSERT(nonnull(projectionFactory.getTargetGlobalIndexer()));

  // Build mass matrix
  std::unordered_map<std::string,double> ebMultipliers = {{"eblock-0_0",1.0},{"eblock-1_0",1.0}};

  timer->start("projectionFactory.buildMassMatrix() 1.0");
  auto massMatrix_1 = projectionFactory.buildMassMatrix(false,&ebMultipliers);
  timer->stop("projectionFactory.buildMassMatrix() 1.0");

  const double multiplier = 4.0;
  ebMultipliers["eblock-0_0"] = multiplier;
  ebMultipliers["eblock-1_0"] = multiplier;

  timer->start("projectionFactory.buildMassMatrix() 4.0");
  auto massMatrix_4 = projectionFactory.buildMassMatrix(false,&ebMultipliers);
  timer->stop("projectionFactory.buildMassMatrix() 4.0");

  const auto view_1 = massMatrix_1->getLocalValuesView();
  const auto view_4 = massMatrix_4->getLocalValuesView();

  const double tol = 100.0 * std::numeric_limits<double>::epsilon();
  int valuesCheck = 0; // if check > 0 then multipliers broke

  Kokkos::parallel_reduce(view_1.size(),KOKKOS_LAMBDA (const int i, int& update) {
      if (std::fabs(multiplier * view_1(i) - view_4(i)) > tol)
        update += 1;
  }, valuesCheck);
  Kokkos::fence();

  TEST_EQUALITY(valuesCheck, 0);
}
