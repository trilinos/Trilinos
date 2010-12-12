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
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManager.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ModelFactory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include <cstdio> // for get char

namespace panzer {

   void pause_to_attach()
   {
      MPI_Comm mpicomm = MPI_COMM_WORLD;
      Teuchos::RCP<Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(
            Teuchos::rcp(new Teuchos::OpaqueWrapper<MPI_Comm>(mpicomm)));
      Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
      out.setShowProcRank(true);
      out.setOutputToRootOnly(-1);

      out << "PID = " << getpid();

      if (comm->getRank() == 0)
         getchar();
      comm->barrier();
   }

  void testInitialzation(const panzer_stk::STK_Interface& mesh,
			 panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs);



  TEUCHOS_UNIT_TEST(field_manager_builder, basic)
  {
    using Teuchos::RCP;
  
    // pause_to_attach();

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",2);
    pl->set("Y Elements",2);
    
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(*mesh, ipb, bcs);

    const std::size_t workset_size = 20;
    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh,ipb, workset_size);

    const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets 
       = panzer_stk::buildBCWorksets(*mesh,ipb,bcs);

    std::map<std::string,std::string> block_ids_to_physics_ids;
    block_ids_to_physics_ids["eblock-0_0"] = "test physics";
    
    std::map<std::string,panzer::InputPhysicsBlock> 
      physics_id_to_input_physics_blocks;
    physics_id_to_input_physics_blocks["test physics"] = ipb;

    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    
    user_app::MyFactory eqset_factory;
				  
    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
      Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);

    user_app::BCFactory bc_factory;

    fmb->setup(conn_manager,
	      MPI_COMM_WORLD,
	      block_ids_to_physics_ids,
	      physics_id_to_input_physics_blocks,
	      volume_worksets,
	      bc_worksets,
	      Teuchos::as<int>(mesh->getDimension()),
	      eqset_factory,
	      bc_factory,
	      workset_size);

    panzer::AssemblyEngine_TemplateManager<panzer::Traits,int,int> ae_tm;
    panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb);
    ae_tm.buildObjects(builder);


    RCP<DOFManager<int,int> > dofManager = 
      ae_tm.getAsObject<panzer::Traits::Residual>()->getManagerBuilder()->getDOFManager();
    RCP<Epetra_Map> ghosted_map = dofManager->getOverlapMap();
    RCP<Epetra_CrsGraph> ghosted_graph = dofManager->getOverlapGraph();
    
    panzer::AssemblyEngineInArgs input;
    input.x = rcp(new Epetra_Vector(*ghosted_map),true);
    input.dxdt = rcp(new Epetra_Vector(*ghosted_map));
    input.f = rcp(new Epetra_Vector(*ghosted_map));
    input.j = rcp(new Epetra_CrsMatrix(Copy, *ghosted_graph));

    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

    input.f->Print(std::cout);
    input.j->Print(std::cout);
  }

  void testInitialzation(const panzer_stk::STK_Interface& mesh,
			 panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    panzer::InputEquationSet ies_1;
    ies_1.name = "Energy";
    ies_1.basis = "Q1";
    ies_1.integration_order = 1;
    ies_1.model_id = 6;
    ies_1.model_factory = "rf";
    ies_1.prefix = "";
    ies_1.params.set<int>("junk", 1);

    ipb.physics_block_id = "test physics";
    ipb.eq_sets.push_back(ies_1);

    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 3.0;
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
      std::size_t bc_id = 2;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 6.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }
    {
      std::size_t bc_id = 3;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "bottom";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 4.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }
  }

}
