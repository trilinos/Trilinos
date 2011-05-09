#ifndef PANZER_STK_MODEL_EVALUATOR_FACTORY_T_HPP
#define PANZER_STK_MODEL_EVALUATOR_FACTORY_T_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Panzer_config.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_MultiBlockMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_ParameterList_ObjectBuilders.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_Basis.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ClosureModel_Factory.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_InitialCondition_Builder.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Panzer_STK_NOXObserverFactory_Epetra.hpp"
#include "Panzer_STK_RythmosObserverFactory_Epetra.hpp"
#include <vector>

// Piro solver objects
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Piro_ConfigDefs.hpp"
#include "Piro_NOXSolver.hpp"
#include "Piro_RythmosSolver.hpp"

#ifdef HAVE_TEKO
#include "Teko_StratimikosFactory.hpp"
#endif

namespace panzer_stk {
  
  template<typename ScalarT>
  void ModelEvaluatorFactory_Epetra<ScalarT>::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
  {
    paramList->validateParametersAndSetDefaults(*this->getValidParameters(), 0);
    this->setMyParamList(paramList);
  }
  
  template<typename ScalarT>
  Teuchos::RCP<const Teuchos::ParameterList> ModelEvaluatorFactory_Epetra<ScalarT>::getValidParameters() const
  {
    static Teuchos::RCP<const Teuchos::ParameterList> validPL;
    if (is_null(validPL)) {
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList());

      pl->sublist("Physics Blocks");
      pl->sublist("Closure Models");
      pl->sublist("Boundary Conditions");
      pl->sublist("Solution Control");
      pl->sublist("Solver Factories");
      pl->sublist("Mesh");
      pl->sublist("Initial Conditions");
      pl->sublist("Output");
      pl->sublist("Output").set("File Name","panzer.exo");  
      pl->sublist("Assembly");
      pl->sublist("Block ID to Physics ID Mapping");
      pl->sublist("Options");
      pl->sublist("User Data");
      pl->sublist("User Data").sublist("Panzer Data");
     
      validPL = pl;
    }
    return validPL;
  }
  
  template<typename ScalarT>
  void  ModelEvaluatorFactory_Epetra<ScalarT>::buildObjects(const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    TEST_FOR_EXCEPTION(Teuchos::is_null(this->getParameterList()), std::runtime_error,
		       "ParameterList must be set before objects can be built!");

    Teuchos::FancyOStream fout(Teuchos::rcpFromRef(std::cout));
    fout.setOutputToRootOnly(0); 
   
    // this function will need to be broken up eventually and probably
    // have parts moved back into panzer.  Just need to get something
    // running.
 
    Teuchos::ParameterList& p = *this->getNonconstParameterList();

    // Build mesh
    Teuchos::ParameterList& mesh_params = p.sublist("Mesh");
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
    Teuchos::RCP<panzer_stk::STK_Interface> mesh;
    if (mesh_params.get<std::string>("Source") ==  "Exodus File") {
      TEUCHOS_ASSERT(true);
      // mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory(mesh_params.sublist("Exodus File").get<std::string>("File Name")));
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory());
      mesh_factory->setParameterList(Teuchos::rcp(new Teuchos::ParameterList(mesh_params.sublist("Exodus File"))));
    }
    else if (mesh_params.get<std::string>("Source") ==  "Inline Mesh") {

      int dimension = mesh_params.sublist("Inline Mesh").get<int>("Mesh Dimension");
      
      if (dimension == 2) {
	mesh_factory = Teuchos::rcp(new panzer_stk::SquareQuadMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }
      else if(dimension == 3) {
	mesh_factory = Teuchos::rcp(new panzer_stk::CubeHexMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }
      else if(dimension==4) {
	mesh_factory = Teuchos::rcp(new panzer_stk::MultiBlockMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }

    }
   
    const Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = 
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    mesh = mesh_factory->buildUncommitedMesh(*(mpi_comm->getRawMpiComm()));
    
    Teuchos::RCP<std::map<std::string,std::string> > block_ids_to_physics_ids = 
      Teuchos::rcp(new std::map<std::string,std::string>);
    panzer::buildBlockIdToPhysicsIdMap(*block_ids_to_physics_ids, p.sublist("Block ID to Physics ID Mapping"));
    
    Teuchos::RCP<std::map<std::string,panzer::InputPhysicsBlock> > physics_id_to_input_physics_blocks = 
      Teuchos::rcp(new std::map<std::string,panzer::InputPhysicsBlock>);
    panzer::buildInputPhysicsBlocks(*physics_id_to_input_physics_blocks, p.sublist("Physics Blocks"));

    Teuchos::RCP<std::vector<panzer::BC> > bcs = 
      Teuchos::rcp(new std::vector<panzer::BC>);
    panzer::buildBCs(*bcs, p.sublist("Boundary Conditions"));
    
    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = 
      Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);

    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;

    std::size_t workset_size = p.sublist("Assembly").get<std::size_t>("Workset Size");
    std::string field_order = p.sublist("Assembly").get<std::string>("Field Order",""); // control nodal ordering of unknown
                                                                                        // global IDs in linear system

    Teuchos::RCP<const panzer::EquationSetFactory> eqset_factory = 
      p.sublist("Assembly").get<Teuchos::RCP<const panzer::EquationSetFactory> >("Equation Set Factory");

    Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > cm_factory = 
      p.sublist("Assembly").get<Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >("Closure Model Factory");
      
    bool is_transient  = p.sublist("Solution Control").get<std::string>("Piro Solver") == "Rythmos" ? true : false;
    
    fmb->buildPhysicsBlocks(*block_ids_to_physics_ids,
			    *physics_id_to_input_physics_blocks,
			    Teuchos::as<int>(mesh->getDimension()),
			    workset_size,
			    *eqset_factory,
			    is_transient,
			    physicsBlocks);

    // finish building mesh, set required field variables and mesh bulk data
    ////////////////////////////////////////////////////////////////////////
    {
      std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator physIter;
      for(physIter=physicsBlocks.begin();physIter!=physicsBlocks.end();++physIter) {
	Teuchos::RCP<const panzer::PhysicsBlock> pb = *physIter;
	const std::vector<panzer::StrBasisPair> & blockFields = pb->getProvidedDOFs();
	
	// insert all fields into a set
	std::set<panzer::StrBasisPair,panzer::StrBasisComp> fieldNames;
	fieldNames.insert(blockFields.begin(),blockFields.end());
	
	// add basis to DOF manager: block specific
	std::set<panzer::StrBasisPair>::const_iterator fieldItr;
	for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr)
	  mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
      }
   
      mesh_factory->completeMeshConstruction(*mesh,*(mpi_comm->getRawMpiComm()));
    }
      
    mesh->print(fout);
    mesh->setupTransientExodusFile(p.sublist("Output").get<std::string>("File Name")); 

    // build worksets
    //////////////////////////////////////////////////////////////
    std::map<std::string,panzer::InputPhysicsBlock> eb_id_to_ipb;
    for (std::map<std::string,std::string>::iterator block = block_ids_to_physics_ids->begin();
	 block != block_ids_to_physics_ids->end(); ++block)
      eb_id_to_ipb[block->first] = (*physics_id_to_input_physics_blocks)[block->second];

    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh, eb_id_to_ipb, workset_size);

    const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets  = 
      panzer_stk::buildBCWorksets(*mesh,eb_id_to_ipb,*bcs);

    // build DOF Manager
    /////////////////////////////////////////////////////////////
 
    // build the connection manager 
    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory<int,int> globalIndexerFactory;
    Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
      = globalIndexerFactory.buildUniqueGlobalIndexer(*(mpi_comm->getRawMpiComm()),physicsBlocks,conn_manager,field_order);
    
    Teuchos::RCP<const Epetra_Comm> ep_comm = Teuchos::rcp(new Epetra_MpiComm(*(mpi_comm->getRawMpiComm())));

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
      = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(ep_comm,dofManager));

    // Add mesh objects to user data to make available to user ctors
    /////////////////////////////////////////////////////////////
    p.sublist("User Data").sublist("Panzer Data").set("STK Mesh", mesh);
    p.sublist("User Data").sublist("Panzer Data").set("DOF Manager", dofManager);
    p.sublist("User Data").sublist("Panzer Data").set("Linear Object Factory", linObjFactory);

    // setup field manager build
    /////////////////////////////////////////////////////////////
 
    Teuchos::RCP<const panzer::BCStrategyFactory> bc_factory = 
      p.sublist("Assembly").get<Teuchos::RCP<const panzer::BCStrategyFactory> >("BC Factory");

    fmb->setupVolumeFieldManagers(volume_worksets,physicsBlocks,*cm_factory,p.sublist("Closure Models"),dofManager,*linObjFactory,p.sublist("User Data"));
    fmb->setupBCFieldManagers(bc_worksets,physicsBlocks,*eqset_factory,*cm_factory,*bc_factory,p.sublist("Closure Models"),*linObjFactory, p.sublist("User Data"));

    // Print Phalanx DAGs
    {
      bool write_dot_files = false;
      write_dot_files = p.sublist("Options").get("Write Volume Assembly Graphs",write_dot_files);
      if (write_dot_files) {
	std::string prefix = "Panzer_AssemblyGraph_";
	prefix = p.sublist("Options").get("Volume Assembly Graph Prefix",prefix);
	fmb->writeVolumeGraphvizDependencyFiles(prefix, physicsBlocks);
      }
    }

    // build solvers
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof =
      Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjFactory<panzer::Traits,int> >(linObjFactory); 
    
    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    {
      Teuchos::RCP<Teuchos::Array<std::string> > p_0 = Teuchos::rcp(new Teuchos::Array<std::string>);
      p_0->push_back("viscosity");
      p_names.push_back(p_0);
    }
    RCP<panzer::ModelEvaluator_Epetra> ep_me = 
      Teuchos::rcp(new panzer::ModelEvaluator_Epetra(fmb,ep_lof, p_names, is_transient));

    // Setup initial conditions
    /////////////////////////////////////////////////////////////
    {
      bool write_dot_files = false;
      std::string prefix = "Panzer_AssemblyGraph_";
      write_dot_files = p.sublist("Options").get("Write Volume Assembly Graphs",write_dot_files);
      prefix = p.sublist("Options").get("Volume Assembly Graph Prefix",prefix);

      std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > > phx_ic_field_managers;
      panzer::setupInitialConditionFieldManagers(volume_worksets,
						 physicsBlocks,
						 *cm_factory,
						 p.sublist("Initial Conditions"),
						 dofManager,
						 *linObjFactory,
						 p.sublist("User Data"),
						 write_dot_files,
						 prefix,
						 phx_ic_field_managers);

      Teuchos::RCP<panzer::LinearObjContainer> loc = linObjFactory->buildLinearObjContainer();
      Teuchos::RCP<panzer::EpetraLinearObjContainer> eloc = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(loc);
      eloc->x = Teuchos::rcp_const_cast<Epetra_Vector>(ep_me->get_x_init());
      
      panzer::evaluateInitialCondition(fmb->getWorksets(), phx_ic_field_managers, loc, 0.0);

      if (is_transient)
	mesh->writeToExodus(0.0);

    }
   
    // Build stratimikos solver (note that this is a hard coded path to linear solver options in nox list!)
    RCP<Teuchos::ParameterList> strat_params = Teuchos::rcp(new Teuchos::ParameterList);
    std::string solver = p.sublist("Solution Control").get<std::string>("Piro Solver");
    {
      *strat_params = p.sublist("Solution Control").sublist("NOX").sublist("Direction").
	sublist("Newton").sublist("Stratimikos Linear Solver").sublist("Stratimikos");
    }

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    #ifdef HAVE_TEKO
       Teko::addTekoToStratimikosBuilder(linearSolverBuilder);
    #endif
    linearSolverBuilder.setParameterList(strat_params);
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);

    // Build Thyra Model Evluator
    RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me = 
      Thyra::epetraModelEvaluator(ep_me,lowsFactory);
    
    m_physics_me = thyra_me;

    RCP<Teuchos::ParameterList> piro_params = Teuchos::rcp(new Teuchos::ParameterList);
    *piro_params = p.sublist("Solution Control");
    RCP<Thyra::ModelEvaluatorDefaultBase<double> > piro;
    if (solver=="NOX") {
      Teuchos::RCP<const panzer_stk::NOXObserverFactory_Epetra> observer_factory = 
	p.sublist("Solver Factories").get<Teuchos::RCP<const panzer_stk::NOXObserverFactory_Epetra> >("NOX Observer Factory");
      Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo = observer_factory->buildNOXObserver(mesh,dofManager,ep_lof);
      piro_params->sublist("NOX").sublist("Solver Options").set("User Defined Pre/Post Operator", ppo);
      piro = Teuchos::rcp(new Piro::NOXSolver<double>(piro_params, thyra_me));
    }
    else if (solver=="Rythmos") {
      Teuchos::RCP<const panzer_stk::RythmosObserverFactory_Epetra> observer_factory = 
	p.sublist("Solver Factories").get<Teuchos::RCP<const panzer_stk::RythmosObserverFactory_Epetra> >("Rythmos Observer Factory");
      piro = Teuchos::rcp(new Piro::RythmosSolver<double>(piro_params, thyra_me, observer_factory->buildRythmosObserver(mesh,dofManager,ep_lof)));
    } 
    else {
      TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error: Unknown Piro Solver : " << solver);
    }

    m_rome_me = piro;

  }

  template<typename ScalarT>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > ModelEvaluatorFactory_Epetra<ScalarT>::getPhysicsModelEvaluator()
  {
    TEST_FOR_EXCEPTION(Teuchos::is_null(m_physics_me), std::runtime_error,
		       "Objects are not built yet!  Please call buildObjects() member function.");
    return  m_physics_me;
  }
  
  template<typename ScalarT>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > ModelEvaluatorFactory_Epetra<ScalarT>::getResponseOnlyModelEvaluator()
  {
    TEST_FOR_EXCEPTION(Teuchos::is_null(m_rome_me), std::runtime_error,
		       "Objects are not built yet!  Please call buildObjects() member function.");
    return m_rome_me;
  }

}

#endif
