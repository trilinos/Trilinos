#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Kokkos_DefaultNode.hpp"
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
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
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

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs);

  TEUCHOS_UNIT_TEST(model_evaluator, basic)
  {
    using Teuchos::RCP;
  
    // pause_to_attach();

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

    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
      Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    user_app::MyFactory eqset_factory;
    user_app::BCFactory bc_factory;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    std::map<std::string,panzer::InputPhysicsBlock> eb_id_to_ipb;

    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";
      
      std::map<std::string,panzer::InputPhysicsBlock> 
        physics_id_to_input_physics_blocks;
      physics_id_to_input_physics_blocks["test physics"] = ipb;
  
      for (std::map<std::string,std::string>::iterator block = block_ids_to_physics_ids.begin();
	   block != block_ids_to_physics_ids.end(); ++block)
	eb_id_to_ipb[block->first] = physics_id_to_input_physics_blocks[block->second];

      bool build_transient_support = true;
      fmb->buildPhysicsBlocks(block_ids_to_physics_ids,
                              physics_id_to_input_physics_blocks,
                              Teuchos::as<int>(mesh->getDimension()), workset_size,
                              eqset_factory,
			      build_transient_support,
                              physicsBlocks);
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh,eb_id_to_ipb, workset_size);

    const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets 
       = panzer_stk::buildBCWorksets(*mesh,eb_id_to_ipb,bcs);

    // build DOF Manager
    /////////////////////////////////////////////////////////////
 
    // build the connection manager 
    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory<int,int> globalIndexerFactory;
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
         = globalIndexerFactory.buildUniqueGlobalIndexer(MPI_COMM_WORLD,physicsBlocks,conn_manager);

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
          = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(Comm,dofManager));

    // setup field manager build
    /////////////////////////////////////////////////////////////
 
    // Add in the application specific closure model factory
    Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > cm_factory = 
      Teuchos::rcp(new panzer::ClosureModelFactory_TemplateManager<panzer::Traits>);
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    (Teuchos::rcp_const_cast<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> >(cm_factory))->buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);

    fmb->setupVolumeFieldManagers(volume_worksets,physicsBlocks,*cm_factory,closure_models,dofManager,*linObjFactory);
    fmb->setupBCFieldManagers(bc_worksets,physicsBlocks,eqset_factory,*cm_factory,bc_factory,closure_models,*linObjFactory);

    typedef double ScalarT;
    typedef std::size_t LO;
    typedef std::size_t GO;
    typedef Kokkos::DefaultNode::DefaultNodeType NODE;
    
    panzer::ModelEvaluator<ScalarT,LO,GO,NODE> model;

    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof =
      Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjFactory<panzer::Traits,int> >(linObjFactory); 

    // Test a transient me
    {
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = true;
      RCP<panzer::ModelEvaluator_Epetra> me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(fmb,ep_lof,p_names,build_transient_support));

      EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();
      
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x));
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot));
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_alpha));
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_beta));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_f));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_W));

      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_sg));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot_sg));
      
      
      RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      RCP<Epetra_Vector> x_dot = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      x->Update(1.0, *(me->get_x_init()), 0.0);
      x_dot->PutScalar(0.0);
      in_args.set_x(x);
      in_args.set_x_dot(x_dot);
      in_args.set_alpha(0.0);
      in_args.set_beta(1.0);
      
      RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
      RCP<Epetra_Operator> J_tmp = me->create_W();
      RCP<Epetra_CrsMatrix> J = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(J_tmp);
      TEST_ASSERT(!Teuchos::is_null(J));
      out_args.set_f(f);
      out_args.set_W(J);

      me->evalModel(in_args, out_args);

    }

    // Test a steady-state me
    {
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = false;
      RCP<panzer::ModelEvaluator_Epetra> me = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(fmb,ep_lof,p_names,build_transient_support));

      EpetraExt::ModelEvaluator::InArgs in_args = me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = me->createOutArgs();
      
      TEST_ASSERT(in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_alpha));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_beta));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_f));
      TEST_ASSERT(out_args.supports(EpetraExt::ModelEvaluator::OUT_ARG_W));

      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_sg));
      TEST_ASSERT(!in_args.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot_sg));
      
      
      RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*me->get_x_map()));
      x->Update(1.0, *(me->get_x_init()), 0.0);
      in_args.set_x(x);
      
      RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(*me->get_f_map()));
      RCP<Epetra_Operator> J_tmp = me->create_W();
      RCP<Epetra_CrsMatrix> J = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(J_tmp);
      TEST_ASSERT(!Teuchos::is_null(J));
      out_args.set_f(f);
      out_args.set_W(J);

      me->evalModel(in_args, out_args);
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
