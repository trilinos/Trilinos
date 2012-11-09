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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ParameterList_ObjectBuilders.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Epetra_Comm.h"
#ifdef HAVE_MPI
   #include "Teuchos_DefaultMpiComm.hpp"
   #include "Epetra_MpiComm.h"
#else
   #include "NO_SERIAL_BUILD.h"
#endif

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraVector.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

#include "Thyra_LinearOpTester.hpp"

#include <cstdio> // for get char

namespace panzer {

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs);

  Teuchos::RCP<const Thyra::LinearOpBase<double> >  eLinearOp;
  Teuchos::RCP<const Thyra::LinearOpBase<double> >  tLinearOp;

  Teuchos::RCP<const Thyra::VectorBase<double> >  eVector;
  Teuchos::RCP<const Thyra::VectorBase<double> >  tVector;

  TEUCHOS_UNIT_TEST(assembly_engine, basic_epetra)
  {
    using Teuchos::RCP;
  
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",6);
    pl->set("Y Elements",4);
    
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb = 
      Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    user_app::MyFactory eqset_factory;
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");
      
      std::map<std::string,panzer::InputPhysicsBlock> 
        physics_id_to_input_physics_blocks;
      physics_id_to_input_physics_blocks["test physics"] = ipb;

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                              block_ids_to_cell_topo,
                              physics_id_to_input_physics_blocks,
                              Teuchos::as<int>(mesh->getDimension()), workset_size,
                              eqset_factory,
			      gd,
			      false,
                              physicsBlocks);
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,physicsBlocks,workset_size));

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager 
    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory<int,int> globalIndexerFactory;
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
         = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);
 
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > eLinObjFactory
          = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(Comm.getConst(),dofManager));
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory = eLinObjFactory;

    // setup field manager build
    /////////////////////////////////////////////////////////////
 
    // Add in the application specific closure model factory
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);

    Teuchos::ParameterList user_data("User Data");

    fmb->setWorksetContainer(wkstContainer);
    fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
    fmb->setupBCFieldManagers(bcs,physicsBlocks,eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

    panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
    panzer::AssemblyEngine_TemplateBuilder builder(fmb,linObjFactory);
    ae_tm.buildObjects(builder);

    RCP<panzer::EpetraLinearObjContainer> eGhosted 
       = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(linObjFactory->buildGhostedLinearObjContainer());
    RCP<panzer::EpetraLinearObjContainer> eGlobal
       = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(linObjFactory->buildLinearObjContainer());
    eLinObjFactory->initializeGhostedContainer(panzer::EpetraLinearObjContainer::X |
                                               panzer::EpetraLinearObjContainer::DxDt |
                                               panzer::EpetraLinearObjContainer::F |
                                               panzer::EpetraLinearObjContainer::Mat,*eGhosted);
    eLinObjFactory->initializeContainer(panzer::EpetraLinearObjContainer::X |
                                        panzer::EpetraLinearObjContainer::DxDt |
                                        panzer::EpetraLinearObjContainer::F |
                                        panzer::EpetraLinearObjContainer::Mat,*eGlobal);
    eGhosted->initialize();
    eGlobal->initialize();
    panzer::AssemblyEngineInArgs input(eGhosted,eGlobal);
    input.alpha = 0.0;
    input.beta = 1.0;

    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

    //input.f->Print(std::cout);
    //input.j->Print(std::cout);

    eLinearOp = Thyra::epetraLinearOp(eGlobal->get_A());
    eVector = Thyra::create_Vector(eGlobal->get_f(),eLinearOp->range());
  }

  TEUCHOS_UNIT_TEST(assembly_engine, basic_tpetra)
  {
    using Teuchos::RCP;

    // build global communicator
    Teuchos::RCP<Teuchos::Comm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
  
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",6);
    pl->set("Y Elements",4);
    
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb = 
      Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    user_app::MyFactory eqset_factory;
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");
      
      std::map<std::string,panzer::InputPhysicsBlock> 
        physics_id_to_input_physics_blocks;
      physics_id_to_input_physics_blocks["test physics"] = ipb;

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                              block_ids_to_cell_topo,
                              physics_id_to_input_physics_blocks,
                              Teuchos::as<int>(mesh->getDimension()), workset_size,
                              eqset_factory,
			      gd,
			      false,
                              physicsBlocks);
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,physicsBlocks,workset_size));

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager 
    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory<int,int> globalIndexerFactory;
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
         = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);
 
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory 
          = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,int>(comm,dofManager));

    // setup field manager build
    /////////////////////////////////////////////////////////////
 
    // Add in the application specific closure model factory
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);

    Teuchos::ParameterList user_data("User Data");

    fmb->setWorksetContainer(wkstContainer);
    fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
    fmb->setupBCFieldManagers(bcs,physicsBlocks,eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

    panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
    panzer::AssemblyEngine_TemplateBuilder builder(fmb,linObjFactory);
    ae_tm.buildObjects(builder);

    RCP<panzer::LinearObjContainer> tGhosted = linObjFactory->buildGhostedLinearObjContainer();
    RCP<panzer::LinearObjContainer> tGlobal = linObjFactory->buildLinearObjContainer();
    linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                              panzer::LinearObjContainer::DxDt |
                                              panzer::LinearObjContainer::F |
                                              panzer::LinearObjContainer::Mat,*tGhosted);
    linObjFactory->initializeContainer(panzer::EpetraLinearObjContainer::X |
                                       panzer::EpetraLinearObjContainer::DxDt |
                                       panzer::EpetraLinearObjContainer::F |
                                       panzer::EpetraLinearObjContainer::Mat,*tGlobal);

    // panzer::pauseToAttach();
    tGhosted->initialize();
    tGlobal->initialize();

    panzer::AssemblyEngineInArgs input(tGhosted,tGlobal);
    input.alpha = 0.0;
    input.beta = 1.0;

    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

    RCP<panzer::TpetraLinearObjContainer<double,int,int> > globalCont 
       = Teuchos::rcp_dynamic_cast<panzer::TpetraLinearObjContainer<double,int,int> >(tGlobal);

    globalCont->get_A()->fillComplete();
    Teuchos::RCP<const Tpetra::Operator<double,int,int> > baseOp = globalCont->get_A();
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > rangeSpace = Thyra::createVectorSpace<double>(baseOp->getRangeMap());
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domainSpace = Thyra::createVectorSpace<double>(baseOp->getDomainMap());

    tLinearOp = Thyra::constTpetraLinearOp<double,int,int>(rangeSpace, domainSpace, baseOp);
    tVector = Thyra::constTpetraVector<double,int,int>(Thyra::tpetraVectorSpace<double,int,int>(baseOp->getRangeMap()).getConst(),
                                                       globalCont->get_f().getConst());
  }

  TEUCHOS_UNIT_TEST(assembly_engine, z_basic_epetra_vtpetra)
  {
     TEUCHOS_ASSERT(tLinearOp!=Teuchos::null);
     TEUCHOS_ASSERT(eLinearOp!=Teuchos::null);

     TEUCHOS_ASSERT(tVector!=Teuchos::null);
     TEUCHOS_ASSERT(eVector!=Teuchos::null);

     Thyra::LinearOpTester<double> tester;
     tester.set_all_error_tol(1e-14);
     tester.show_all_tests(true);
     tester.dump_all(true);
     tester.num_random_vectors(200);

     {
        const bool result = tester.compare( *tLinearOp, *eLinearOp, Teuchos::ptrFromRef(out) );
        TEST_ASSERT(result);
     }

     {
        const bool result = Thyra::testRelNormDiffErr(
           "Tpetra",*tVector,
           "Epetra",*eVector,
           "linear_properties_error_tol()", 1e-14,
           "linear_properties_warning_tol()", 1e-14,
           &out);
        TEST_ASSERT(result);
     }
  }

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    panzer::InputEquationSet ies_1;
    ies_1.name = "Energy";
    ies_1.basis = "Q2";
    ies_1.integration_order = 1;
    ies_1.model_id = "solid";
    ies_1.prefix = "";

    panzer::InputEquationSet ies_2;
    ies_2.name = "Energy";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = "ion solid";
    ies_2.prefix = "ION_";

    ipb.physics_block_id = "4";
    ipb.eq_sets.push_back(ies_1);
    ipb.eq_sets.push_back(ies_2);


    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }    
    {
      std::size_t bc_id = 1;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }   
    {
      std::size_t bc_id = 2;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }
  }

}
