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
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_GlobalData.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_WorksetContainer.hpp"

#include "Panzer_ResponseEvaluatorFactory_IPCoordinates.hpp"

#include "TestEvaluators.hpp"

#include <vector>
#include <map>
#include <string>

namespace panzer {

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  struct RespFactoryIPCoords_Builder {
    int cubatureDegree;

    template <typename T>
    Teuchos::RCP<ResponseEvaluatorFactoryBase> build() const
    { return Teuchos::rcp(new ResponseEvaluatorFactory_IPCoordinates<T>(cubatureDegree)); }
  };

  TEUCHOS_UNIT_TEST(response_library_stk2, test)
  {
    using Teuchos::RCP;


  #ifdef HAVE_MPI
     Teuchos::RCP<const Teuchos::MpiComm<int> > tcomm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
  #else
     Teuchos::RCP<Teuchos::Comm<int> > tcomm = FAIL
  #endif

    panzer_stk::SquareQuadMeshFactory mesh_factory;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    const std::size_t workset_size = 1;

    panzer::FieldManagerBuilder fmb;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",2);
       pl->set("X0",0.0);
       pl->set("Y0",0.0);
       pl->set("Xf",8.0);
       pl->set("Yf",4.0);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // setup physic blocks
    /////////////////////////////////////////////
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physics_blocks;
    {
       testInitialzation(ipb, bcs);

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
                                  physics_blocks);
    }

    // setup worksets
    /////////////////////////////////////////////

     std::vector<std::string> validEBlocks;
     mesh->getElementBlockNames(validEBlocks);

    // build WorksetContainer
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physics_blocks.size();i++)
      wkstContainer->setNeeds(physics_blocks[i]->elementBlockID(),physics_blocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);

    // setup DOF manager
    /////////////////////////////////////////////
    const Teuchos::RCP<panzer::ConnManager> conn_manager
           = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    Teuchos::RCP<const panzer::GlobalIndexerFactory > indexerFactory
          = Teuchos::rcp(new panzer::DOFManagerFactory);
    const Teuchos::RCP<panzer::GlobalIndexer> dofManager
          = indexerFactory->buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physics_blocks,conn_manager);

    // and linear object factory
    Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > elof
          = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tcomm.getConst(),dofManager));

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lof = elof;

    // setup field manager builder
    /////////////////////////////////////////////

    // Add in the application specific closure model factory
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("Response Model");

    Teuchos::ParameterList user_data("User Data");
    user_data.sublist("Panzer Data").set("Mesh", mesh);
    user_data.sublist("Panzer Data").set("DOF Manager", dofManager);
    user_data.sublist("Panzer Data").set("Linear Object Factory", lof);
    user_data.set<int>("Workset Size",workset_size);

    // IP is in center of element
    user_data.sublist("IP Coordinates").set<int>("Integration Order",1);

    // setup and evaluate ResponseLibrary
    ///////////////////////////////////////////////////

    out << "Adding responses" << std::endl;

    RCP<ResponseLibrary<Traits> > rLibrary
          = Teuchos::rcp(new ResponseLibrary<Traits>(wkstContainer,dofManager,lof));


    ResponseEvaluatorFactory_IPCoordinates_Builder builder;
    builder.cubatureDegree = 1;

    std::vector<std::string> blocks(1);
    blocks[0] = "eblock-0_0";
    rLibrary->addResponse("IPCoordinates-0_0",blocks,builder);
    blocks[0] = "eblock-1_0";
    rLibrary->addResponse("IPCoordinates-1_0",blocks,builder);

    Teuchos::RCP<ResponseBase> resp00 = rLibrary->getResponse<panzer::Traits::Residual>("IPCoordinates-0_0");
    Teuchos::RCP<ResponseBase> resp10 = rLibrary->getResponse<panzer::Traits::Residual>("IPCoordinates-1_0");

    TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_IPCoordinates<panzer::Traits::Residual> >(resp00,true));
    TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_IPCoordinates<panzer::Traits::Residual> >(resp10,true));

    rLibrary->buildResponseEvaluators(physics_blocks,
  				      cm_factory,
                                      closure_models,
  				      user_data,true);

    Teuchos::RCP<panzer::LinearObjContainer> loc = lof->buildLinearObjContainer();
    lof->initializeContainer(panzer::LinearObjContainer::X,*loc);
    Teuchos::RCP<panzer::LinearObjContainer> gloc = lof->buildGhostedLinearObjContainer();
    lof->initializeGhostedContainer(panzer::LinearObjContainer::X,*gloc);

    out << "evaluating VFM" << std::endl;
    panzer::AssemblyEngineInArgs ae_inargs(gloc,loc);
    rLibrary->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
    rLibrary->evaluate<panzer::Traits::Residual>(ae_inargs);

    std::map<std::string,Teuchos::RCP<const std::vector<panzer::Traits::Residual::ScalarT> > > coords;
    coords["eblock-0_0"] = Teuchos::rcp_dynamic_cast<Response_IPCoordinates<panzer::Traits::Residual> >(resp00,true)->getCoords();;
    coords["eblock-1_0"] = Teuchos::rcp_dynamic_cast<Response_IPCoordinates<panzer::Traits::Residual> >(resp10,true)->getCoords();;

    // Debugging
    if (true) {
      Teuchos::RCP<Teuchos::FancyOStream> out2 = Teuchos::getFancyOStream(Teuchos::rcp(&out,false));
      out2->setOutputToRootOnly(-1);
      *out2 << "\nPrinting IP coordinates for block: eblock-0_0" << std::endl;
      for (std::vector<panzer::Traits::Residual::ScalarT>::const_iterator i = (coords["eblock-0_0"])->begin(); i != (coords["eblock-0_0"])->end(); ++i)
	*out2 << "pid = " << tcomm->getRank() << ", val = " << *i << std::endl;
      *out2 << "\nPrinting IP coordinates for block: eblock-1_0" << std::endl;
      for (std::vector<panzer::Traits::Residual::ScalarT>::const_iterator i = (coords["eblock-1_0"])->begin(); i != (coords["eblock-1_0"])->end(); ++i)
	*out2 << "pid = " << tcomm->getSize() << ", val = " << *i << std::endl;
    }


    const double double_tol = 10.0*std::numeric_limits<double>::epsilon();

    // NOTE: if the ordering of elements in STK changes or the
    // ordering of integration points in Intrepid2 changes, this test
    // will break!  It assumes a fixed deterministic ordering.
    if (tcomm->getSize() == 1) {
      // eblock 1
      {
	const std::vector<panzer::Traits::Residual::ScalarT>& values = *(coords["eblock-0_0"]);
	TEST_FLOATING_EQUALITY(values[0], 1.0, double_tol);  // x
	TEST_FLOATING_EQUALITY(values[1], 3.0, double_tol);  // x
	TEST_FLOATING_EQUALITY(values[2], 1.0, double_tol);  // x
	TEST_FLOATING_EQUALITY(values[3], 3.0, double_tol);  // x
	TEST_FLOATING_EQUALITY(values[4], 1.0, double_tol);  // y
	TEST_FLOATING_EQUALITY(values[5], 1.0, double_tol);  // y
	TEST_FLOATING_EQUALITY(values[6], 3.0, double_tol);  // y
	TEST_FLOATING_EQUALITY(values[7], 3.0, double_tol);  // y
      }
      // eblock 2
      {
	const std::vector<panzer::Traits::Residual::ScalarT>& values = *(coords["eblock-1_0"]);
	TEST_FLOATING_EQUALITY(values[0], 5.0, double_tol);  // x
	TEST_FLOATING_EQUALITY(values[1], 7.0, double_tol);  // x
	TEST_FLOATING_EQUALITY(values[2], 5.0, double_tol);  // x
	TEST_FLOATING_EQUALITY(values[3], 7.0, double_tol);  // x
	TEST_FLOATING_EQUALITY(values[4], 1.0, double_tol);  // y
	TEST_FLOATING_EQUALITY(values[5], 1.0, double_tol);  // y
	TEST_FLOATING_EQUALITY(values[6], 3.0, double_tol);  // y
	TEST_FLOATING_EQUALITY(values[7], 3.0, double_tol);  // y
      }
    }
    else if (tcomm->getSize() == 2) {

      if (tcomm->getRank() == 0) {
	// eblock 1
	{
	  const std::vector<panzer::Traits::Residual::ScalarT>& values = *(coords["eblock-0_0"]);
	  TEST_FLOATING_EQUALITY(values[0], 1.0, double_tol);  // x
	  TEST_FLOATING_EQUALITY(values[1], 1.0, double_tol);  // x
	  TEST_FLOATING_EQUALITY(values[2], 1.0, double_tol);  // y
	  TEST_FLOATING_EQUALITY(values[3], 3.0, double_tol);  // y
	}
	// eblock 2
	{
	  const std::vector<panzer::Traits::Residual::ScalarT>& values = *(coords["eblock-1_0"]);
	  TEST_FLOATING_EQUALITY(values[0], 5.0, double_tol);  // x
	  TEST_FLOATING_EQUALITY(values[1], 5.0, double_tol);  // x
	  TEST_FLOATING_EQUALITY(values[2], 1.0, double_tol);  // y
	  TEST_FLOATING_EQUALITY(values[3], 3.0, double_tol);  // y
	}
      }
      else if (tcomm->getRank() == 1) {
	// eblock 1
	{
	  const std::vector<panzer::Traits::Residual::ScalarT>& values = *(coords["eblock-0_0"]);
	  TEST_FLOATING_EQUALITY(values[0], 3.0, double_tol);  // x
	  TEST_FLOATING_EQUALITY(values[1], 3.0, double_tol);  // x
	  TEST_FLOATING_EQUALITY(values[2], 1.0, double_tol);  // y
	  TEST_FLOATING_EQUALITY(values[3], 3.0, double_tol);  // y
	}
	// eblock 2
	{
	  const std::vector<panzer::Traits::Residual::ScalarT>& values = *(coords["eblock-1_0"]);
	  TEST_FLOATING_EQUALITY(values[0], 7.0, double_tol);  // x
	  TEST_FLOATING_EQUALITY(values[1], 7.0, double_tol);  // x
	  TEST_FLOATING_EQUALITY(values[2], 1.0, double_tol);  // y
	  TEST_FLOATING_EQUALITY(values[3], 3.0, double_tol);  // y
	}
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Error - this test can only be run with 1 or 2 processes!");
    }

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
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }


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
