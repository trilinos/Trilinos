#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Kokkos_DefaultNode.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
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
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_STK_RythmosObserver_Epetra.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ModelFactory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include <cstdio> // for get char

// Piro solver objects
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Piro_ConfigDefs.hpp"
#include "Piro_NOXSolver.hpp"
#include "Piro_RythmosSolver.hpp"

#include "write_solution_data.hpp"

namespace panzer {

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs);

  // *****************************************************************
  // NOX Test
  // *****************************************************************
  TEUCHOS_UNIT_TEST(solver, NOX_steady_state)
  {
    using Teuchos::RCP;
  
    // pauseToAttach();

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",20);
    pl->set("Y Elements",20);
    
    panzer_stk::SquareQuadMeshFactory mesh_factory;
    mesh_factory.setParameterList(pl);
    //RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<panzer_stk::STK_Interface> mesh = 
      mesh_factory.buildUncommitedMesh(MPI_COMM_WORLD);

    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
      Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    user_app::MyFactory eqset_factory;
    user_app::BCFactory bc_factory;
    const std::size_t workset_size = 20;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;

    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      //block_ids_to_physics_ids["eblock-1_0"] = "test physics";
      
      std::map<std::string,panzer::InputPhysicsBlock> 
        physics_id_to_input_physics_blocks;
      physics_id_to_input_physics_blocks["test physics"] = ipb;
  
      fmb->buildPhysicsBlocks(block_ids_to_physics_ids,
                              physics_id_to_input_physics_blocks,
                              Teuchos::as<int>(mesh->getDimension()),
			      workset_size,
                              eqset_factory,
                              physicsBlocks);
    }

   // finish building mesh, set required field variables and mesh bulk data
   ////////////////////////////////////////////////////////////////////////

   {
      std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator physIter;
      for(physIter=physicsBlocks.begin();physIter!=physicsBlocks.end();++physIter) {
         Teuchos::RCP<const panzer::PhysicsBlock> pb = *physIter;
         const std::vector<StrBasisPair> & blockFields = pb->getProvidedDOFs();

         // insert all fields into a set
         std::set<StrBasisPair,StrBasisComp> fieldNames;
         fieldNames.insert(blockFields.begin(),blockFields.end());

         // add basis to DOF manager: block specific
         std::set<StrBasisPair>::const_iterator fieldItr;
         for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {

            mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
         }
      }

      mesh_factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
   }

    // build worksets
    //////////////////////////////////////////////////////////////
    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh,ipb, workset_size);

    const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets 
       = panzer_stk::buildBCWorksets(*mesh,ipb,bcs);

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
 
    fmb->setupVolumeFieldManagers(volume_worksets,physicsBlocks,dofManager,*linObjFactory);
    fmb->setupBCFieldManagers(bc_worksets,physicsBlocks,eqset_factory,bc_factory,*linObjFactory);

    // Print Phalanx DAGs
    fmb->writeVolumeGraphvizDependencyFiles("Panzer_", physicsBlocks);

    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof =
      Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjFactory<panzer::Traits,int> >(linObjFactory); 
    
    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    {
      Teuchos::RCP<Teuchos::Array<std::string> > p_0 = Teuchos::rcp(new Teuchos::Array<std::string>);
      p_0->push_back("viscosity");
      p_names.push_back(p_0);
    }
    bool is_transient = false;
    RCP<panzer::ModelEvaluator_Epetra> ep_me = 
      Teuchos::rcp(new panzer::ModelEvaluator_Epetra(fmb,ep_lof, p_names,
						     is_transient));

    // Get solver params from input file
    RCP<Teuchos::ParameterList> piro_params = rcp(new Teuchos::ParameterList("Piro Parameters"));
    Teuchos::updateParametersFromXmlFile("solver_nox.xml", piro_params.get());
    
    // Build stratimikos solver
    std::string& solver = piro_params->get("Piro Solver","NOX");
    RCP<Teuchos::ParameterList> stratParams;
    
    if (solver=="NOX" || solver=="LOCA") {
      stratParams = Teuchos::rcp(&(piro_params->sublist("NOX").sublist("Direction").
				   sublist("Newton").sublist("Stratimikos Linear Solver").sublist("Stratimikos")),false);
    }
    else if (solver=="Rythmos") {
      stratParams = Teuchos::rcp(&(piro_params->sublist("Rythmos").sublist("Stratimikos")),false);
    }
    
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(stratParams);
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
      createLinearSolveStrategy(linearSolverBuilder);

    // Build Thyra Model Evluator
    RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me = 
      Thyra::epetraModelEvaluator(ep_me,lowsFactory);
    
    RCP<Thyra::ModelEvaluatorDefaultBase<double> > piro;
    if (solver=="NOX") {
      piro = rcp(new Piro::NOXSolver<double>(piro_params, thyra_me));
    }
    else if (solver=="Rythmos") {
      piro = rcp(new Piro::RythmosSolver<double>(piro_params, thyra_me));
    }
    else {
      TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error: Unknown Piro Solver : " << solver);
    }
    
    // Now the (somewhat cumbersome) setting of inputs and outputs
    Thyra::ModelEvaluatorBase::InArgs<double> inArgs = piro->createInArgs();
    const Thyra::ModelEvaluatorBase::InArgs<double> inArgsNominal = piro->getNominalValues();
    int num_p = inArgs.Np();     // Number of *vectors* of parameters
    assert (num_p == 1);  // Logic needs to be generalized -- hardwire to 1 p vector in model
    RCP<Thyra::VectorBase<double> > p1 = Thyra::createMember(*piro->get_p_space(0));
    Thyra::copy(*inArgsNominal.get_p(0), p1.ptr());
    //int numParams = p1->space()->dim(); // Number of parameters in p1 vector
      
    inArgs.set_p(0,p1);

    // Set output arguments to evalModel call
    Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = piro->createOutArgs();
    int num_g = outArgs.Ng(); // Number of *vectors* of responses
    assert (num_g == 1);

    //RCP<Thyra::MultiVectorBase<double> > dgdp =
    //Thyra::createMembers(*thyra_me->get_x_space(),numParams);
    //if (computeSens) outArgs.set_DgDp(0, 0, dgdp);
    
    // Solution vector is returned as extra respons vector
    RCP<Thyra::VectorBase<double> > gx = Thyra::createMember(*thyra_me->get_x_space());
    outArgs.set_g(0,gx);

    // Now, solve the problem and return the responses
    piro->evalModel(inArgs, outArgs);

    //std::cout << *gx << std::endl;

    // Export solution to ghosted vector for exodus output
    RCP<Epetra_Vector> solution = Thyra::get_Epetra_Vector(*(ep_lof->getMap()), gx);
    Epetra_Vector ghosted_solution(*(ep_lof->getGhostedMap()));
    RCP<Epetra_Import> importer = ep_lof->getGhostedImport();
    ghosted_solution.PutScalar(0.0);
    ghosted_solution.Import(*solution,*importer,Insert);

    write_solution_data(*Teuchos::rcp_dynamic_cast<panzer::DOFManager<int,int> >(dofManager),*mesh,
			ghosted_solution);
    mesh->writeToExodus("output.exo");

    // Test solution values on left, middle, and right side of mesh.
    // Note that this is based on the exact 20x20 test mesh on 4
    // processes.  It will fail on more or less processes due to
    // global node re-numbering.
    {
      Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm = DefaultComm<Ordinal>::getComm();

      if (Teuchos::size(*comm) == 4) {
	
	const int gid_left = 0;
	const int gid_middle = 225;
	const int gid_right = 340;
	const int lid_left = solution->Map().LID(gid_left);
	const int lid_middle = solution->Map().LID(gid_middle);
	const int lid_right = solution->Map().LID(gid_right);
	
	//solution->Print(std::cout);
	
	if (lid_left != -1) {
	  double left_value = (*solution)[lid_left];
	  TEST_FLOATING_EQUALITY(left_value, 1.0, 1.0e-9);
	}
	
	if (lid_middle != -1) {
	  double middle_value = (*solution)[lid_middle];
	  TEST_FLOATING_EQUALITY(middle_value, 1.625, 1.0e-9);
	}
	
	if (lid_right != -1) {
	  double right_value = (*solution)[lid_right];
	  TEST_FLOATING_EQUALITY(right_value, 2.0, 1.0e-9);
	}
	
	int globalSuccess_int = -1;
	Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
	TEST_EQUALITY_CONST( globalSuccess_int, 0 );
      }
    }

    // This is for debugging and to test the evaluation of the
    // residual and JAcobian at the same time.  Currently NOX
    // evaluates them sparately, so this method is not tested in the
    // model evaluator.
    {
      using Teuchos::rcp;
      RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*ep_me->get_x_map()));
      RCP<Epetra_Vector> f = rcp(new Epetra_Vector(*ep_me->get_f_map()));
      RCP<Epetra_Operator> J = ep_me->create_W();
      x->PutScalar(0.0);
      EpetraExt::ModelEvaluator::InArgs in_args = ep_me->createInArgs();
      EpetraExt::ModelEvaluator::OutArgs out_args = ep_me->createOutArgs();
      in_args.set_x(x);
      out_args.set_f(f);
      out_args.set_W(J);
      ep_me->evalModel(in_args, out_args);
      //Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(J)->Print(std::cout);
      ep_me->evalModel(in_args, out_args);
      //Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(J)->Print(std::cout);
    }
    
  }

  // *****************************************************************
  // Rythmos Test
  // *****************************************************************
  TEUCHOS_UNIT_TEST(solver, Rythmos_transient)
  {
    using Teuchos::RCP;
    
    // pause_to_attach();

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",20);
    pl->set("Y Elements",20);
    
    panzer_stk::SquareQuadMeshFactory mesh_factory;
    mesh_factory.setParameterList(pl);
    //RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<panzer_stk::STK_Interface> mesh = 
      mesh_factory.buildUncommitedMesh(MPI_COMM_WORLD);

    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
      Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    user_app::MyFactory eqset_factory;
    user_app::BCFactory bc_factory;
    const std::size_t workset_size = 20;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;

    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      //block_ids_to_physics_ids["eblock-1_0"] = "test physics";
      
      std::map<std::string,panzer::InputPhysicsBlock> 
        physics_id_to_input_physics_blocks;
      physics_id_to_input_physics_blocks["test physics"] = ipb;
  
      fmb->buildPhysicsBlocks(block_ids_to_physics_ids,
                              physics_id_to_input_physics_blocks,
                              Teuchos::as<int>(mesh->getDimension()),
			      workset_size,
                              eqset_factory,
                              physicsBlocks);
    }

   // finish building mesh, set required field variables and mesh bulk data
   ////////////////////////////////////////////////////////////////////////

   {
      std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator physIter;
      for(physIter=physicsBlocks.begin();physIter!=physicsBlocks.end();++physIter) {
         Teuchos::RCP<const panzer::PhysicsBlock> pb = *physIter;
         const std::vector<StrBasisPair> & blockFields = pb->getProvidedDOFs();

         // insert all fields into a set
         std::set<StrBasisPair,StrBasisComp> fieldNames;
         fieldNames.insert(blockFields.begin(),blockFields.end());

         // add basis to DOF manager: block specific
         std::set<StrBasisPair>::const_iterator fieldItr;
         for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {

            mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
         }
      }

      mesh_factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
   }

   mesh->setupTransientExodusFile("transient.exo"); 
   
    // build worksets
    //////////////////////////////////////////////////////////////
    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh,ipb, workset_size);

    const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets 
       = panzer_stk::buildBCWorksets(*mesh,ipb,bcs);

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
 
    fmb->setupVolumeFieldManagers(volume_worksets,physicsBlocks,dofManager,*linObjFactory);
    fmb->setupBCFieldManagers(bc_worksets,physicsBlocks,eqset_factory,bc_factory,*linObjFactory);

    // Print Phalanx DAGs
    fmb->writeVolumeGraphvizDependencyFiles("Panzer_", physicsBlocks);

    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof =
      Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjFactory<panzer::Traits,int> >(linObjFactory); 
    
    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    {
      Teuchos::RCP<Teuchos::Array<std::string> > p_0 = Teuchos::rcp(new Teuchos::Array<std::string>);
      p_0->push_back("viscosity");
      p_names.push_back(p_0);
    }
    bool is_transient = true;
    RCP<panzer::ModelEvaluator_Epetra> ep_me = 
      Teuchos::rcp(new panzer::ModelEvaluator_Epetra(fmb, ep_lof, p_names, is_transient));

    // Get solver params from input file
    RCP<Teuchos::ParameterList> piro_params = rcp(new Teuchos::ParameterList("Piro Parameters"));
    Teuchos::updateParametersFromXmlFile("solver_rythmos.xml", piro_params.get());
    
    // Build stratimikos solver
    std::string& solver = piro_params->get("Piro Solver","NOX");
    RCP<Teuchos::ParameterList> stratParams;
    
    if (solver=="NOX" || solver=="LOCA") {
      stratParams = Teuchos::rcp(&(piro_params->sublist("NOX").sublist("Direction").
				   sublist("Newton").sublist("Stratimikos Linear Solver").sublist("Stratimikos")),false);
    }
    else if (solver=="Rythmos") {
      stratParams = Teuchos::rcp(&(piro_params->sublist("Rythmos").sublist("Stratimikos")),false);
    }
    
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(stratParams);
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
      createLinearSolveStrategy(linearSolverBuilder);

    // Build Thyra Model Evluator
    RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me = 
      Thyra::epetraModelEvaluator(ep_me,lowsFactory);
    
    RCP<Thyra::ModelEvaluatorDefaultBase<double> > piro;
    if (solver=="NOX") {
      piro = rcp(new Piro::NOXSolver<double>(piro_params, thyra_me));
    }
    else if (solver=="Rythmos") {
      RCP<panzer_stk::RythmosObserver_Epetra> observer = 
	Teuchos::rcp(new panzer_stk::RythmosObserver_Epetra(mesh, dofManager, ep_lof));
      piro = rcp(new Piro::RythmosSolver<double>(piro_params, thyra_me, observer));
    }
    else {
      TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error: Unknown Piro Solver : " << solver);
    }
    
    // Now the (somewhat cumbersome) setting of inputs and outputs
    Thyra::ModelEvaluatorBase::InArgs<double> inArgs = piro->createInArgs();
    const Thyra::ModelEvaluatorBase::InArgs<double> inArgsNominal = piro->getNominalValues();
    int num_p = inArgs.Np();     // Number of *vectors* of parameters
    assert (num_p == 1);  // Logic needs to be generalized -- hardwire to 1 p vector in model
    RCP<Thyra::VectorBase<double> > p1 = Thyra::createMember(*piro->get_p_space(0));
    Thyra::copy(*inArgsNominal.get_p(0), p1.ptr());
    //int numParams = p1->space()->dim(); // Number of parameters in p1 vector
      
    inArgs.set_p(0,p1);

    // Set output arguments to evalModel call
    Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = piro->createOutArgs();
    int num_g = outArgs.Ng(); // Number of *vectors* of responses
    assert (num_g == 1);

    //RCP<Thyra::MultiVectorBase<double> > dgdp =
    //Thyra::createMembers(*thyra_me->get_x_space(),numParams);
    //if (computeSens) outArgs.set_DgDp(0, 0, dgdp);
    
    // Solution vector is returned as extra respons vector
    RCP<Thyra::VectorBase<double> > gx = Thyra::createMember(*thyra_me->get_x_space());
    outArgs.set_g(0,gx);

    // Now, solve the problem and return the responses
    piro->evalModel(inArgs, outArgs);

    //std::cout << *gx << std::endl;

    // Export solution to ghosted vector for exodus output
    RCP<Epetra_Vector> solution = Thyra::get_Epetra_Vector(*(ep_lof->getMap()), gx);
    Epetra_Vector ghosted_solution(*(ep_lof->getGhostedMap()));
    RCP<Epetra_Import> importer = ep_lof->getGhostedImport();
    ghosted_solution.PutScalar(0.0);
    ghosted_solution.Import(*solution,*importer,Insert);

    write_solution_data(*Teuchos::rcp_dynamic_cast<panzer::DOFManager<int,int> >(dofManager),*mesh,
			ghosted_solution);
    
  }

  // *****************************************************************
  // Support functions
  // *****************************************************************
  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    panzer::InputEquationSet ies_1;
    ies_1.name = "Energy";
    ies_1.basis = "Q1";
    ies_1.integration_order = 2;
    ies_1.model_id = 6;
    ies_1.model_factory = "rf";
    ies_1.prefix = "";
    ies_1.params.set<int>("junk", 1);

    panzer::InputEquationSet ies_2;
    ies_2.name = "Energy";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = 6;
    ies_2.model_factory = "rf";
    ies_2.prefix = "ION_";
    ies_2.params.set<int>("junk", 1);

    ipb.physics_block_id = "4";
    ipb.eq_sets.push_back(ies_1);
    //ipb.eq_sets.push_back(ies_2);


    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 1.0;
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
      double value = 2.0;
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
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      //bcs.push_back(bc);
    }
  }
  
}
