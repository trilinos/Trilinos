// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
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
#include "Thyra_TestingTools.hpp"

#include <cstdio> // for get char

namespace panzer {

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
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

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb =
      Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      int default_integration_order = 1;

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
				 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
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
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physicsBlocks.size();i++)
      wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory globalIndexerFactory;
    RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

    Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > eLinObjFactory
          = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));
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
    fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

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

    eLinearOp = Thyra::epetraLinearOp(eGlobal->get_A());
    eVector = Thyra::create_Vector(eGlobal->get_f(),eLinearOp->range());

//     std::cout << *eVector << std::endl;
//     std::cout << *eLinearOp << std::endl;

  }

  TEUCHOS_UNIT_TEST(assembly_engine, basic_epetra_in_stages)
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

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb =
      Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      int default_integration_order = 1;

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
				 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
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
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physicsBlocks.size();i++)
      wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory globalIndexerFactory;
    RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

    Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > eLinObjFactory
          = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));
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
    fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

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

    int flags1 = panzer::AssemblyEngine<panzer::Traits::Residual>::EvaluationFlags::Initialize +panzer::AssemblyEngine<panzer::Traits::Residual>::EvaluationFlags::VolumetricFill;
    int flags2 = panzer::AssemblyEngine<panzer::Traits::Residual>::EvaluationFlags::All - flags1;
    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input, panzer::AssemblyEngine<panzer::Traits::Residual>::EvaluationFlags(flags1));
    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input, panzer::AssemblyEngine<panzer::Traits::Residual>::EvaluationFlags(flags2));
    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

    eLinearOp = Thyra::epetraLinearOp(eGlobal->get_A());
    eVector = Thyra::create_Vector(eGlobal->get_f(),eLinearOp->range());

//     std::cout << *eVector << std::endl;
//     std::cout << *eLinearOp << std::endl;

  }

  TEUCHOS_UNIT_TEST(assembly_engine, dirichlet_only)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;


    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",6);
    pl->set("Y Elements",4);

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb =
      Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      int default_integration_order = 1;

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
				 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
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
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physicsBlocks.size();i++)
      wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory globalIndexerFactory;
    RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

    Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > eLinObjFactory
          = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));
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
    fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

    panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
    panzer::AssemblyEngine_TemplateBuilder builder(fmb,linObjFactory);
    ae_tm.buildObjects(builder);

    RCP<panzer::LinearObjContainer> ghosted = linObjFactory->buildGhostedLinearObjContainer();
    RCP<panzer::LinearObjContainer> global  = linObjFactory->buildLinearObjContainer();
    eLinObjFactory->initializeGhostedContainer(panzer::EpetraLinearObjContainer::X |
                                               panzer::EpetraLinearObjContainer::DxDt |
                                               panzer::EpetraLinearObjContainer::F,*ghosted);
    eLinObjFactory->initializeContainer(panzer::EpetraLinearObjContainer::X |
                                        panzer::EpetraLinearObjContainer::DxDt |
                                        panzer::EpetraLinearObjContainer::F,*global);
    ghosted->initialize();
    global->initialize();

    Thyra::assign(Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(global,true)->get_x_th().ptr(),1.0);
    Thyra::assign(Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(ghosted,true)->get_x_th().ptr(),-33.0);
    panzer::AssemblyEngineInArgs input(ghosted,global);
    input.alpha = 0.0;
    input.beta = 1.0;

    Teuchos::RCP<panzer::LinearObjContainer> counter = ae_tm.getAsObject<panzer::Traits::Residual>()->evaluateOnlyDirichletBCs(input);
    TEST_ASSERT(counter!=Teuchos::null);

    out << "evaluated" << std::endl;
    RCP<Thyra::VectorBase<double> > count = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(counter,true)->get_f_th();
    RCP<Thyra::VectorBase<double> > values = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(global,true)->get_f_th();

    TEST_ASSERT(count!=Teuchos::null);
    TEST_ASSERT(values!=Teuchos::null);

    Teuchos::ArrayRCP<double> count_array, values_array;
    rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(count,true)->getNonconstLocalData(Teuchos::ptrFromRef(count_array));
    rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(values,true)->getNonconstLocalData(Teuchos::ptrFromRef(values_array));

    bool passed = true;
    TEST_EQUALITY(count_array.size(),values_array.size());
    for(Teuchos::ArrayRCP<double>::size_type i=0;i<count_array.size();i++) {
      if(count_array[i]==0.0)
        passed &= (values_array[i]==0.0);
      else {
        passed &= (std::fabs((values_array[i]-(-4.0))/4.0)<1e-14);
        out << values_array[i] << " " << i <<std::endl;
      }
    }
    TEST_ASSERT(passed);

  }

  TEUCHOS_UNIT_TEST(assembly_engine, basic_tpetra)
  {
    using Teuchos::RCP;

    {
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

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb =
      Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      const int default_integration_order = 1;

      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
				 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
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
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physicsBlocks.size();i++)
      wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory globalIndexerFactory;
    RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
          = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal>(comm,dofManager));

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
    fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

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

    RCP<panzer::TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> > globalCont
       = Teuchos::rcp_dynamic_cast<panzer::TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> >(tGlobal);

    Teuchos::RCP<const Tpetra::Operator<double,int,panzer::GlobalOrdinal> > baseOp = globalCont->get_A();
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > rangeSpace = Thyra::createVectorSpace<double>(baseOp->getRangeMap());
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domainSpace = Thyra::createVectorSpace<double>(baseOp->getDomainMap());

    tLinearOp = Thyra::constTpetraLinearOp<double,int,panzer::GlobalOrdinal>(rangeSpace, domainSpace, baseOp);
    tVector = Thyra::constTpetraVector<double,int,panzer::GlobalOrdinal>(Thyra::tpetraVectorSpace<double,int,panzer::GlobalOrdinal>(baseOp->getRangeMap()).getConst(),
                                                       globalCont->get_f().getConst());
    }
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


     // Need to kill global objects so that memory leaks on kokkos ref
     // count pointing doesn't trigger test failure.
     eLinearOp = Teuchos::null;
     tLinearOp = Teuchos::null;
     eVector = Teuchos::null;
     tVector = Teuchos::null;

  }

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    // Physics block
    Teuchos::ParameterList& physics_block = ipb->sublist("test physics");
    {
      Teuchos::ParameterList& p = physics_block.sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",2);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
    }

    // BCs
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
