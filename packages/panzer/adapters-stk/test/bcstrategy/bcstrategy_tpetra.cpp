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

#include <Teuchos_Assert.hpp>
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "Panzer_BC.hpp"
#include "Panzer_BCStrategy.hpp"
#include "Panzer_Traits.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "user_app_BCStrategy_Factory.hpp"
#include <iostream>

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_GlobalData.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_EquationSetFactory.hpp"

using ST = double;
using LO = panzer::LocalOrdinal;
using GO = panzer::GlobalOrdinal;
using NT = panzer::TpetraNodeType;
using TpetraVector = Tpetra::Vector<ST, LO, GO, NT>;
using TpetraLOC = panzer::TpetraLinearObjContainer<ST, LO, GO, NT>;
using TpetraRowMatrix = Tpetra::RowMatrix<ST, LO, GO, NT>;

namespace panzer {

void testInitialzationTpetra(const Teuchos::RCP<Teuchos::ParameterList> &ipb,
                             std::vector<panzer::BC> &bcs);

TEUCHOS_UNIT_TEST(bcstrategy, basic_construction_tpetra) {

  std::size_t bc_id = 0;
  panzer::BCType dirichlet = BCT_Dirichlet;
  std::string sideset_id = "4";
  std::string element_block_id = "fluid";
  std::string dof_name = "UX";
  std::string strategy = "Constant";
  double value = 5.0;
  Teuchos::ParameterList p;
  p.set("Value", value);
  panzer::BC bc(bc_id, dirichlet, sideset_id, element_block_id, dof_name,
                strategy, p);

  Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits>> bcs;

  Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

  user_app::BCFactory my_factory;
  bcs = my_factory.buildBCStrategy(bc, gd);
}

TEUCHOS_UNIT_TEST(bcstrategy, constant_bc_strategy_tpetra) {

  using std::cout;
  using std::endl;
  using Teuchos::RCP;

  auto pl = rcp(new Teuchos::ParameterList);
  pl->set("X Blocks", 1);
  pl->set("Y Blocks", 1);
  pl->set("X Elements", 1);
  pl->set("Y Elements", 1);

  panzer_stk::SquareQuadMeshFactory factory;
  factory.setParameterList(pl);
  auto mesh = factory.buildMesh(MPI_COMM_WORLD);

  auto ipb = Teuchos::parameterList("Physics Blocks");
  std::vector<panzer::BC> bcs;
  testInitialzationTpetra(ipb, bcs);

  Teuchos::RCP<panzer::FieldManagerBuilder> fmb =
      Teuchos::rcp(new panzer::FieldManagerBuilder);

  // build physics blocks
  //////////////////////////////////////////////////////////////
  const std::size_t workset_size = 1;
  auto eqset_factory = Teuchos::rcp(new user_app::MyFactory);
  user_app::BCFactory bc_factory;
  std::vector<Teuchos::RCP<panzer::PhysicsBlock>> physicsBlocks;
  {
    std::map<std::string, std::string> block_ids_to_physics_ids;
    block_ids_to_physics_ids["eblock-0_0"] = "test physics";

    std::map<std::string, Teuchos::RCP<const shards::CellTopology>>
        block_ids_to_cell_topo;
    block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");

    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

    const int default_integration_order = 1;

    panzer::buildPhysicsBlocks(block_ids_to_physics_ids, block_ids_to_cell_topo,
                               ipb, default_integration_order, workset_size,
                               eqset_factory, gd, false, physicsBlocks);
  }

  // build worksets
  //////////////////////////////////////////////////////////////
  auto wkstFactory = Teuchos::rcp(
      new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
  Teuchos::RCP<panzer::WorksetContainer>
      wkstContainer // attach it to a workset container (uses lazy evaluation)
      = Teuchos::rcp(new panzer::WorksetContainer);
  wkstContainer->setFactory(wkstFactory);
  for (size_t i = 0; i < physicsBlocks.size(); i++)
    wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),
                            physicsBlocks[i]->getWorksetNeeds());
  wkstContainer->setWorksetSize(workset_size);

  // build DOF Manager
  /////////////////////////////////////////////////////////////

  // build the connection manager
  const Teuchos::RCP<panzer::ConnManager> conn_manager =
      Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  auto tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

  panzer::DOFManagerFactory globalIndexerFactory;
  RCP<panzer::GlobalIndexer> dofManager =
      globalIndexerFactory.buildGlobalIndexer(
          Teuchos::opaqueWrapper(MPI_COMM_WORLD), physicsBlocks, conn_manager);

  //    auto dofManagerCasted =
  //    Teuchos::rcp_dynamic_cast<panzer::BlockedDOFManager>(dofManager);
  //    Teuchos::RCP<panzer::BlockedTpetraLinearObjFactory<panzer::Traits,ST,LO,GO>
  //    > eLinObjFactory
  //          = Teuchos::rcp(new
  //          panzer::BlockedTpetraLinearObjFactory<panzer::Traits,ST,LO,GO>(tComm.getConst(),
  //          dofManagerCasted));

  // and linear object factory
  auto eLinObjFactory = Teuchos::rcp(
      new panzer::TpetraLinearObjFactory<
          panzer::Traits, double, panzer::LocalOrdinal, panzer::GlobalOrdinal>(
          tComm.getConst(), dofManager));

  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits>> linObjFactory =
      eLinObjFactory;

  // setup field manager build
  /////////////////////////////////////////////////////////////

  // Add in the application specific closure model factory
  panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
  user_app::MyModelFactory_TemplateBuilder cm_builder;
  cm_factory.buildObjects(cm_builder);

  Teuchos::ParameterList closure_models("Closure Models");
  closure_models.sublist("solid")
      .sublist("SOURCE_TEMPERATURE")
      .set<double>("Value", 1.0);
  closure_models.sublist("ion solid")
      .sublist("SOURCE_ION_TEMPERATURE")
      .set<double>("Value", 1.0);

  Teuchos::ParameterList user_data("User Data");

  fmb->setWorksetContainer(wkstContainer);
  fmb->setupVolumeFieldManagers(physicsBlocks, cm_factory, closure_models,
                                *linObjFactory, user_data);
  fmb->setupBCFieldManagers(bcs, physicsBlocks, *eqset_factory, cm_factory,
                            bc_factory, closure_models, *linObjFactory,
                            user_data);

  panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
  panzer::AssemblyEngine_TemplateBuilder builder(fmb, linObjFactory);
  ae_tm.buildObjects(builder);

  auto eGhosted = Teuchos::rcp_dynamic_cast<TpetraLOC>(
      linObjFactory->buildGhostedLinearObjContainer());
  auto eGlobal = Teuchos::rcp_dynamic_cast<TpetraLOC>(
      linObjFactory->buildLinearObjContainer());
  eLinObjFactory->initializeGhostedContainer(TpetraLOC::X | TpetraLOC::DxDt |
                                                 TpetraLOC::F | TpetraLOC::Mat,
                                             *eGhosted);
  eLinObjFactory->initializeContainer(
      TpetraLOC::X | TpetraLOC::DxDt | TpetraLOC::F | TpetraLOC::Mat, *eGlobal);
  panzer::AssemblyEngineInArgs input(eGhosted, eGlobal);

  // MPL TORM input.container_: LinearObjContainer
  const auto x = Teuchos::rcp_dynamic_cast<TpetraLOC>(input.container_);

  x->get_x()->putScalar(1.0);
  input.beta = 1.0;

  ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
  ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

  // Check residual values.  Evaluation should have put (x - 5.0)
  // into each residual.  With initial guess of 1.0, check to make
  // sure each entry in residual has -4.0.  Note that we are using
  // one element with same dirichlet bc on each side, so all nodes
  // have same dirichlet bc applied to it.

  auto f = Teuchos::rcp_dynamic_cast<TpetraLOC>(input.container_)->get_f();
  ST tol = 10.0 * std::numeric_limits<ST>::epsilon();
  auto fData = f->getData(0);
  for (int i = 0; i < fData.size(); ++i) {
    TEST_FLOATING_EQUALITY(fData[i], -4.0, tol);
  }

  // Check Jacobian values.  Should have one on diagonal and zero
  // elsewhere.
  auto jac = Teuchos::rcp_dynamic_cast<TpetraLOC>(input.container_)->get_A();
  for (size_t i = 0; i < jac->getLocalNumRows(); ++i) {
    TpetraRowMatrix::local_inds_host_view_type indices;
    TpetraRowMatrix::values_host_view_type values;
    jac->getLocalRowView(i, indices, values);
    GO gI = jac->getRowMap()->getGlobalElement(i);
    for (size_t j = 0; j < indices.extent(0); j++) {
      GO gJ = jac->getColMap()->getGlobalElement(indices[j]);
      std::cout << "J(" << gI << "," << gJ << ") = " << values[j] << std::endl;
      if (gI == gJ) {
        // Diagonal entry; value should be row GO
        TEST_FLOATING_EQUALITY(values[j], 1.0, tol);
      } else {
        // Nondiagonal entry; value should be 0.0
        TEST_FLOATING_EQUALITY(values[j], 0.0, tol);
      }
    }
  }

  jac->print(std::cout);
}

void testInitialzationTpetra(const Teuchos::RCP<Teuchos::ParameterList> &ipb,
                             std::vector<panzer::BC> &bcs) {
  // Physics block
  Teuchos::ParameterList &physics_block = ipb->sublist("test physics");
  {
    Teuchos::ParameterList &p = physics_block.sublist("a");
    p.set("Type", "Energy");
    p.set("Prefix", "");
    p.set("Model ID", "solid");
    p.set("Basis Type", "HGrad");
    p.set("Basis Order", 1);
  }

  {
    std::size_t bc_id = 0;
    panzer::BCType dirichlet = BCT_Dirichlet;
    std::string sideset_id = "left";
    std::string element_block_id = "eblock-0_0";
    std::string dof_name = "TEMPERATURE";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value", value);
    panzer::BC bc(bc_id, dirichlet, sideset_id, element_block_id, dof_name,
                  strategy, p);
    bcs.push_back(bc);
  }
  {
    std::size_t bc_id = 1;
    panzer::BCType dirichlet = BCT_Dirichlet;
    std::string sideset_id = "right";
    std::string element_block_id = "eblock-0_0";
    std::string dof_name = "TEMPERATURE";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value", value);
    panzer::BC bc(bc_id, dirichlet, sideset_id, element_block_id, dof_name,
                  strategy, p);
    bcs.push_back(bc);
  }
  {
    std::size_t bc_id = 2;
    panzer::BCType dirichlet = BCT_Dirichlet;
    std::string sideset_id = "top";
    std::string element_block_id = "eblock-0_0";
    std::string dof_name = "TEMPERATURE";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value", value);
    panzer::BC bc(bc_id, dirichlet, sideset_id, element_block_id, dof_name,
                  strategy, p);
    bcs.push_back(bc);
  }
  {
    std::size_t bc_id = 3;
    panzer::BCType dirichlet = BCT_Dirichlet;
    std::string sideset_id = "bottom";
    std::string element_block_id = "eblock-0_0";
    std::string dof_name = "TEMPERATURE";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value", value);
    panzer::BC bc(bc_id, dirichlet, sideset_id, element_block_id, dof_name,
                  strategy, p);
    bcs.push_back(bc);
  }
}

} // namespace panzer
