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
#include "Panzer_STK_Utilities.hpp"
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
#include "Panzer_STK_ParameterListCallback.hpp"
#include <vector>

// Piro solver objects
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Piro_ConfigDefs.hpp"
#include "Piro_NOXSolver.hpp"
#include "Piro_RythmosSolver.hpp"

#include "EpetraExt_VectorOut.h"

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

    // for convience cast to an MPI comm
    const Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = 
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
   
    // this function will need to be broken up eventually and probably
    // have parts moved back into panzer.  Just need to get something
    // running.
 
    Teuchos::ParameterList& p = *this->getNonconstParameterList();

    // "parse" parameter list
    Teuchos::ParameterList & mesh_params     = p.sublist("Mesh");
    Teuchos::ParameterList & assembly_params = p.sublist("Assembly");
    Teuchos::ParameterList & solncntl_params = p.sublist("Solution Control");

    Teuchos::ParameterList & user_data_params = p.sublist("User Data");
    Teuchos::ParameterList & panzer_data_params = user_data_params.sublist("Panzer Data");

    // Build mesh factory and uncommitted mesh
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory = this->buildSTKMeshFactory(mesh_params);
    Teuchos::RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(*(mpi_comm->getRawMpiComm()));
    
    // setup physical mappings and boundary conditions
    std::map<std::string,std::string> block_ids_to_physics_ids;
    panzer::buildBlockIdToPhysicsIdMap(block_ids_to_physics_ids, p.sublist("Block ID to Physics ID Mapping"));
    
    std::map<std::string,panzer::InputPhysicsBlock> physics_id_to_input_physics_blocks; 
    panzer::buildInputPhysicsBlocks(physics_id_to_input_physics_blocks, p.sublist("Physics Blocks"));

    std::vector<panzer::BC> bcs;
    panzer::buildBCs(bcs, p.sublist("Boundary Conditions"));
    
    // extract assembly information
    std::size_t workset_size = assembly_params.get<std::size_t>("Workset Size");
    std::string field_order  = assembly_params.get<std::string>("Field Order",""); // control nodal ordering of unknown
                                                                                   // global IDs in linear system

    Teuchos::RCP<const panzer::EquationSetFactory> eqset_factory = 
      assembly_params.get<Teuchos::RCP<const panzer::EquationSetFactory> >("Equation Set Factory");

    Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > cm_factory = 
      assembly_params.get<Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >("Closure Model Factory");
      
    // this is weird...we are accessing the solution control to determine if things are transient
    // it is backwards!
    bool is_transient  = solncntl_params.get<std::string>("Piro Solver") == "Rythmos" ? true : false;
 
    // build physics blocks

    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb = Teuchos::rcp(new panzer::FieldManagerBuilder<int,int>);
    fmb->buildPhysicsBlocks(block_ids_to_physics_ids,
			    physics_id_to_input_physics_blocks,
			    Teuchos::as<int>(mesh->getDimension()),
			    workset_size,
			    *eqset_factory,
			    is_transient,
			    physicsBlocks);

    // finish building mesh, set required field variables and mesh bulk data
    ////////////////////////////////////////////////////////////////////////
    try {
       // this throws some exceptions, catch them as neccessary
       this->finalizeMeshConstruction(*mesh_factory,physicsBlocks,*mpi_comm,*mesh);
    } catch(const panzer_stk::STK_Interface::ElementBlockException & ebexp) {
       fout << "*****************************************\n\n";
       fout << "Element block exception, could not finalize the mesh, printing block and sideset information:\n";
       fout.pushTab(3);
       mesh->printMetaData(fout);
       fout.popTab();
       fout << std::endl;

       throw ebexp;
    } catch(const panzer_stk::STK_Interface::SidesetException & ssexp) {
       fout << "*****************************************\n\n";
       fout << "Sideset exception, could not finalize the mesh, printing block and sideset information:\n";
       fout.pushTab(3);
       mesh->printMetaData(fout);
       fout.popTab();
       fout << std::endl;

       throw ssexp;
    }

    mesh->print(fout);
    mesh->setupTransientExodusFile(p.sublist("Output").get<std::string>("File Name")); 

    // build worksets
    //////////////////////////////////////////////////////////////
    std::map<std::string,panzer::InputPhysicsBlock> eb_id_to_ipb;
    for (std::map<std::string,std::string>::iterator block = block_ids_to_physics_ids.begin();
	 block != block_ids_to_physics_ids.end(); ++block)
      eb_id_to_ipb[block->first] = physics_id_to_input_physics_blocks[block->second];

    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh, eb_id_to_ipb, workset_size);

    const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets  = 
      panzer_stk::buildBCWorksets(*mesh,eb_id_to_ipb,bcs);

    // build DOF Manager
    /////////////////////////////////////////////////////////////
 
    // build the connection manager 
    const Teuchos::RCP<panzer_stk::STKConnManager> stkConn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager = stkConn_manager;

    panzer::DOFManagerFactory<int,int> globalIndexerFactory;
    Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
      = globalIndexerFactory.buildUniqueGlobalIndexer(*(mpi_comm->getRawMpiComm()),physicsBlocks,conn_manager,field_order);
    
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
      = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(mpi_comm,dofManager));

    // Add mesh objects to user data to make available to user ctors
    /////////////////////////////////////////////////////////////
    panzer_data_params.set("STK Mesh", mesh);
    panzer_data_params.set("DOF Manager", dofManager);
    panzer_data_params.set("Linear Object Factory", linObjFactory);

    // setup field manager build
    /////////////////////////////////////////////////////////////
 
    Teuchos::RCP<const panzer::BCStrategyFactory> bc_factory = 
      p.sublist("Assembly").get<Teuchos::RCP<const panzer::BCStrategyFactory> >("BC Factory");

    fmb->setupVolumeFieldManagers(volume_worksets,physicsBlocks,*cm_factory,p.sublist("Closure Models"),dofManager,*linObjFactory,user_data_params);
    fmb->setupBCFieldManagers(bc_worksets,physicsBlocks,*eqset_factory,*cm_factory,*bc_factory,p.sublist("Closure Models"),*linObjFactory,user_data_params);

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

      // Push solution into ghosted epetra vector and then into stk for outputting
      { 
	Epetra_Vector ghosted_solution(*(ep_lof->getGhostedMap()));
	Teuchos::RCP<Epetra_Import> importer = ep_lof->getGhostedImport();
	ghosted_solution.PutScalar(0.0);
	ghosted_solution.Import(*(eloc->x),*importer,Insert);
	
	write_solution_data(*Teuchos::rcp_dynamic_cast<panzer::DOFManager<int,int> >(dofManager),*mesh,
			    ghosted_solution);
      }

      if (is_transient)
	mesh->writeToExodus(0.0);

    }
   
    // Build stratimikos solver (note that this is a hard coded path to linear solver options in nox list!)
    RCP<Teuchos::ParameterList> strat_params = Teuchos::rcp(new Teuchos::ParameterList);
    std::string solver = solncntl_params.get<std::string>("Piro Solver");
    {
      *strat_params = solncntl_params.sublist("NOX").sublist("Direction").
	sublist("Newton").sublist("Stratimikos Linear Solver").sublist("Stratimikos");
    }

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    #ifdef HAVE_TEKO 
    {
       std::string fieldName;

       Teuchos::RCP<Teko::RequestHandler> reqHandler = Teuchos::rcp(new Teko::RequestHandler);
       Teuchos::RCP<const panzer::DOFManager<int,int> > dofs =
          Teuchos::rcp_dynamic_cast<const panzer::DOFManager<int,int> >(dofManager);

       // add in the coordinate parameter list callback handler
       if(determineCoordinateField(*dofs,fieldName)) {
          std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > fieldPatterns;
          fillFieldPatternMap(*dofs,fieldName,fieldPatterns);
          reqHandler->addRequestCallback(Teuchos::rcp(new 
                panzer_stk::ParameterListCallback<int,int>(fieldName,fieldPatterns,stkConn_manager,dofManager)));

          Teuchos::RCP<panzer_stk::ParameterListCallback<int,int> > callback = Teuchos::rcp(new 
                panzer_stk::ParameterListCallback<int,int>(fieldName,fieldPatterns,stkConn_manager,dofManager));
          reqHandler->addRequestCallback(callback);

          bool writeCoordinates = p.sublist("Options").get("Write Coordinates",false);
          if(writeCoordinates) {
             // force parameterlistcallback to build coordinates
             callback->preRequest(Teko::RequestMesg(Teuchos::rcp(new Teuchos::ParameterList())));
             
             // extract coordinate vectors
             const std::vector<double> & xcoords = callback->getXCoordsVector();
             const std::vector<double> & ycoords = callback->getYCoordsVector();
             const std::vector<double> & zcoords = callback->getZCoordsVector();

             // use epetra to write coordinates to matrix market files
             Epetra_MpiComm ep_comm(*mpi_comm->getRawMpiComm());
             Epetra_Map map(-1,xcoords.size(),0,ep_comm);

             Teuchos::RCP<Epetra_Vector> vec;
             switch(mesh->getDimension()) {
             case 3:
                vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&zcoords[0])));
                EpetraExt::VectorToMatrixMarketFile("zcoords.mm",*vec);
             case 2:
                vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&ycoords[0])));
                EpetraExt::VectorToMatrixMarketFile("ycoords.mm",*vec);
             case 1:
                vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&xcoords[0])));
                EpetraExt::VectorToMatrixMarketFile("xcoords.mm",*vec);
                break;
             default:
                TEUCHOS_ASSERT(false);
             }
          }
       }
       // else write_out_the_mesg("Warning: No unique field determines the coordinates, coordinates unavailable!")   

       Teko::addTekoToStratimikosBuilder(linearSolverBuilder,reqHandler);
    }
    #endif
    linearSolverBuilder.setParameterList(strat_params);
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);

    // Build Thyra Model Evluator
    RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me = 
      Thyra::epetraModelEvaluator(ep_me,lowsFactory);
    
    m_physics_me = thyra_me;

    RCP<Teuchos::ParameterList> piro_params = Teuchos::rcp(new Teuchos::ParameterList(solncntl_params));
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

  //! build STK mesh from a mesh parameter list
  template<typename ScalarT>
  Teuchos::RCP<panzer_stk::STK_MeshFactory> ModelEvaluatorFactory_Epetra<ScalarT>::buildSTKMeshFactory(const Teuchos::ParameterList & mesh_params) const
  {
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;

    // first contruct the mesh factory
    if (mesh_params.get<std::string>("Source") ==  "Exodus File") {
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
      else if(dimension==4) { // not really "dimension==4" simply a flag to try this other mesh for testing
	mesh_factory = Teuchos::rcp(new panzer_stk::MultiBlockMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }
    }

    return mesh_factory;
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory_Epetra<ScalarT>::finalizeMeshConstruction(const STK_MeshFactory & mesh_factory,
                                                                       const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                                                       const Teuchos::MpiComm<int> mpi_comm, 
                                                                       STK_Interface & mesh) const
  {
    // finish building mesh, set required field variables and mesh bulk data
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
	  mesh.addSolutionField(fieldItr->first,pb->elementBlockID());
      }
   
      mesh_factory.completeMeshConstruction(mesh,*(mpi_comm.getRawMpiComm()));
    }
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

  template<typename ScalarT>
  bool ModelEvaluatorFactory_Epetra<ScalarT>::determineCoordinateField(
                                   const panzer::DOFManager<int,int> & dofManager,std::string & fieldName) const
  {
     std::vector<string> elementBlocks;
     dofManager.getElementBlockIds(elementBlocks);
 
     // grab fields for first block
     std::set<int> runningFields(dofManager.getFields(elementBlocks[0]));

     // loop over all element blocks intersecting the fields 
     for(std::size_t b=1;b<elementBlocks.size();b++) {
        std::string blockId = elementBlocks[b];

        const std::set<int> & fields = dofManager.getFields(blockId);

        std::set<int> currentFields(runningFields);
        runningFields.clear();
        std::set_intersection(fields.begin(),fields.end(),
                              currentFields.begin(),currentFields.end(),
                              std::inserter(runningFields,runningFields.begin()));
     }

     if(runningFields.size()<1) 
        return false;

     fieldName = dofManager.getFieldString(*runningFields.begin());
     return true;
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory_Epetra<ScalarT>::fillFieldPatternMap(const panzer::DOFManager<int,int> & dofManager,
                                                                  const std::string & fieldName, 
                                                                  std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fieldPatterns) const
  {
     std::vector<string> elementBlocks;
     dofManager.getElementBlockIds(elementBlocks);

     for(std::size_t e=0;e<elementBlocks.size();e++) {
        std::string blockId = elementBlocks[e];
        
        if(dofManager.fieldInBlock(fieldName,blockId))
           fieldPatterns[blockId] =
              Teuchos::rcp_dynamic_cast<const panzer::IntrepidFieldPattern>(dofManager.getFieldPattern(blockId,fieldName),true);
     }
  }

}

#endif
