// Standard headers
#include <cstdio>
#include <iomanip>
#include <unistd.h>

// TrilinosCouplings headers
#include "TrilinosCouplings_config.h"

// Teuchos headers
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Belos headers
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"

// MueLu headers
#include "MueLu.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MutuallyExclusiveTime.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_Utilities.hpp"

// Region MG headers
#include "SetupRegionUtilities.hpp"
#include "SetupRegionVector_def.hpp"
#include "SetupRegionMatrix_def.hpp"
#include "SetupRegionHierarchy_def.hpp"

// Shards headers
#include "Shards_CellTopology.hpp"

// Panzer headers
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_CheckBCConsistency.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_Response_Functional.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_MeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"

// Percept headers
#include "percept/PerceptMesh.hpp"

// Tpetra headers
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_FECrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Include factories for boundary conditions and other Panzer setup
// Most of which is taken from PoissonExample in Panzer_STK
#include "muelu_region_poisson.hpp"


int main(int argc, char *argv[]) {

  using ST = double;
  using LO = panzer::LocalOrdinal;
  using GO = panzer::GlobalOrdinal;
  using NT = panzer::TpetraNodeType;
  using OP = Tpetra::Operator<ST,LO,GO,NT>;
  using MV = Tpetra::MultiVector<ST,LO,GO,NT>; 

  Kokkos::initialize(argc,argv);
  { // Kokkos scope

  /**********************************************************************************/
  /************************************** SETUP *************************************/
  /**********************************************************************************/

  // Setup output stream, MPI, and grab info
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  Teuchos::RCP<const Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

  const int numRanks = comm->getSize();
  const int myRank = comm->getRank();

  if (myRank == 0) {
    *out << "Running MueLu region driver on " << numRanks << " ranks... \n";
  }

  // Parse command line arguments
  Teuchos::CommandLineProcessor clp(false);
  std::string exodusFileName        = "";                  clp.setOption("exodus-mesh",           &exodusFileName,          "Exodus hex mesh filename");
  std::string pamgenFileName        = "cylinder.rtp";      clp.setOption("pamgen-mesh",           &pamgenFileName,          "Pamgen hex mesh filename");
  std::string xmlFileName           = "";                  clp.setOption("xml",                   &xmlFileName,             "MueLu parameters file");
  int mesh_refinements              = 1;                   clp.setOption("mesh-refinements",      &mesh_refinements,        "Uniform mesh refinements");

  bool print_percept_mesh           = false;                clp.setOption("print-percept-mesh", "no-print-percept-mesh", &print_percept_mesh, "Calls perceptMesh's print_info routine");
  bool delete_parent_elements       = false;                clp.setOption("delete-parent-elements", "keep-parent-elements", &delete_parent_elements,"Save the parent elements in the perceptMesh");
  bool useStackedTimer              = false;                clp.setOption("stacked-timer","no-stacked-timer", &useStackedTimer, "use stacked timer");
  bool showTimerSummary             = false;                clp.setOption("show-timer-summary", "no-show-timer-summary", &showTimerSummary, "Switch on/off the timer summary at the end of the run.");

  /* GH: come back to these later
  bool optPrintLocalStats           = false;               clp.setOption("localstats", "nolocalstats", &optPrintLocalStats, "print per-process statistics");
  std::string convergenceLog        = "residual_norm.txt"; clp.setOption("convergence-log",       &convergenceLog,        "file in which the convergence history of the linear solver is stored");
  int         maxIts                = 200;                 clp.setOption("its",                   &maxIts,                "maximum number of solver iterations");
  std::string smootherType          = "Jacobi";            clp.setOption("smootherType",          &smootherType,          "smoother to be used: (None | Jacobi | Gauss | Chebyshev)");
  int         smootherIts           = 2;                   clp.setOption("smootherIts",           &smootherIts,           "number of smoother iterations");
  double      smootherDamp          = 0.67;                clp.setOption("smootherDamp",          &smootherDamp,          "damping parameter for the level smoother");
  double      smootherChebyEigRatio = 2.0;                 clp.setOption("smootherChebyEigRatio", &smootherChebyEigRatio, "eigenvalue ratio max/min used to approximate the smallest eigenvalue for Chebyshev relaxation");
  double      smootherChebyBoostFactor = 1.1;              clp.setOption("smootherChebyBoostFactor", &smootherChebyBoostFactor, "boost factor for Chebyshev smoother");
  double      tol                   = 1e-12;               clp.setOption("tol",                   &tol,                   "solver convergence tolerance");
  bool        scaleResidualHist     = true;                clp.setOption("scale", "noscale",      &scaleResidualHist,     "scaled Krylov residual history");
  bool        serialRandom          = false;               clp.setOption("use-serial-random", "no-use-serial-random", &serialRandom, "generate the random vector serially and then broadcast it");
  bool        keepCoarseCoords      = false;               clp.setOption("keep-coarse-coords", "no-keep-coarse-coords", &keepCoarseCoords, "keep coordinates on coarsest level of region hierarchy");
  bool        coarseSolverRebalance = false;               clp.setOption("rebalance-coarse", "no-rebalance-coarse", &coarseSolverRebalance, "rebalance before AMG coarse grid solve");
  int         rebalanceNumPartitions = 1;                  clp.setOption("numPartitions",         &rebalanceNumPartitions, "number of partitions for rebalancing the coarse grid AMG solve");
  std::string coarseSolverType      = "direct";            clp.setOption("coarseSolverType",      &coarseSolverType,      "Type of solver for (composite) coarse level operator (smoother | direct | amg)");
  std::string unstructured          = "{}";                clp.setOption("unstructured",          &unstructured,          "List of ranks to be treated as unstructured, e.g. {0, 2, 5}");
  std::string coarseAmgXmlFile      = "";                  clp.setOption("coarseAmgXml",          &coarseAmgXmlFile,      "Read parameters for AMG as coarse level solve from this xml file.");
  std::string coarseSmootherXMLFile = "";                  clp.setOption("coarseSmootherXML",     &coarseSmootherXMLFile, "File containing the parameters to use with the coarse level smoother.");
  std::string equilibrate           = "no" ;               clp.setOption("equilibrate",           &equilibrate,           "equilibrate the system (no | diag | 1-norm)");
  int  cacheSize                    = 0;                   clp.setOption("cachesize",               &cacheSize,           "cache size (in KB)");
  std::string cycleType             = "V";                 clp.setOption("cycleType", &cycleType, "{Multigrid cycle type. Possible values: V, W.");
  */
  
  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }
  


  // get xml file from command line if provided, otherwise use default
  std::string  xmlSolverInFileName(xmlFileName);

  // Read xml file into parameter list
  Teuchos::ParameterList inputSolverList;

  if(xmlSolverInFileName.length()) {
    if (myRank == 0)
      *out << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n" << std::endl;
    Teuchos::updateParametersFromXmlFile (xmlSolverInFileName, Teuchos::ptr(&inputSolverList));
  }
  else if (myRank == 0)
    *out << "Using default solver values ..." << std::endl;


  /**********************************************************************************/
  /************************************* MESH ***************************************/
  /**********************************************************************************/
  
  //Teuchos::RCP<Teuchos::Time> meshTimer = Teuchos::TimeMonitor::getNewCounter("Step 1: Mesh generation");
  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
  if(useStackedTimer)
    stacked_timer = rcp(new Teuchos::StackedTimer("MueLu_Driver"));
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
  RCP<Teuchos::TimeMonitor> globalTimeMonitor = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: S - Global Time")));
  RCP<Teuchos::TimeMonitor> tm                = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 1 - Build Mesh and Assign Physics")));  

  Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
  Teuchos::RCP<Teuchos::ParameterList> mesh_pl = Teuchos::rcp(new Teuchos::ParameterList);
  std::string celltype = "Hex";  
  
  if(exodusFileName.length())
  {
    // set the filename and type
    mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory);
    mesh_pl->set("File Name",exodusFileName);
    mesh_pl->set("File Type","Exodus");
    mesh_pl->set("Levels of Uniform Refinement",mesh_refinements); // this may impose 2^(level) coarsening
    mesh_pl->set("Keep Percept Data",true); // this is necessary to gather mesh hierarchy information
    mesh_pl->set("Keep Percept Parent Elements",!delete_parent_elements); // this is necessary to gather mesh hierarchy information

    // check if the user provided conflicting mesh arguments
    if(pamgenFileName.length())
      *out << "Warning: using --exodus-mesh despite --pamgen-mesh being provided as well" << std::endl;
  }
  else if(pamgenFileName.length())
  {
    // set the filename and type
    mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory);
    mesh_pl->set("File Name",pamgenFileName);
    mesh_pl->set("File Type","Pamgen");
    mesh_pl->set("Levels of Uniform Refinement",mesh_refinements); // this may impose 2^(level) coarsening
    mesh_pl->set("Keep Percept Data",true); // this is necessary to gather mesh hierarchy information
    mesh_pl->set("Keep Percept Parent Elements",!delete_parent_elements); // this is necessary to gather mesh hierarchy information
  }
  else
    throw std::runtime_error("no mesh file name found!");

  // set the parameters
  mesh_factory->setParameterList(mesh_pl);
  
  // build the mesh
  Teuchos::RCP<panzer_stk::STK_Interface> mesh;
  mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);

  // setup the physics block
  Teuchos::RCP<Example::EquationSetFactory> eqset_factory = Teuchos::rcp(new Example::EquationSetFactory);
  Example::BCStrategyFactory bc_factory;
  const std::size_t workset_size = 10;
  const int discretization_order = 1;
  
  Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
  std::vector<panzer::BC> bcs;
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;

  std::vector<std::string> eBlocks;
  mesh->getElementBlockNames(eBlocks);

  std::vector<std::string> sidesets;
  mesh->getSidesetNames(sidesets);

  // set physics and boundary conditions
  {
    bool build_transient_support = false;

    const int integration_order = 10;
    Teuchos::ParameterList& p = ipb->sublist("Poisson Physics");
    p.set("Type","Poisson");
    p.set("Model ID","solid");
    p.set("Basis Type","HGrad");
    p.set("Basis Order",discretization_order);
    p.set("Integration Order",integration_order);

    // GH: double-check. this assumes we impose Dirichlet BCs on all boundaries of all physics blocks
    // which may potentially assign Dirichlet BCs to internal block boundaries, which is undesirable  
    for(size_t i=0; i<eBlocks.size(); ++i)
    {
      for(size_t j=0; j<sidesets.size(); ++j)
      {
        std::size_t bc_id = i;    
        panzer::BCType bctype = panzer::BCT_Dirichlet;
        std::string sideset_id = sidesets[j];
        std::string element_block_id = eBlocks[i];
        std::string dof_name = "TEMPERATURE";
        std::string strategy = "Constant";
        double value = 0.0;
        Teuchos::ParameterList p;
        p.set("Value",value);
        panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
                      strategy, p);
        bcs.push_back(bc);
      }    
    const panzer::CellData volume_cell_data(workset_size, mesh->getCellTopology(eBlocks[i]));

    // GobalData sets ostream and parameter interface to physics
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

    // Can be overridden by the equation set
    int default_integration_order = 1;
    
    // the physics block nows how to build and register evaluator with the field manager
    Teuchos::RCP<panzer::PhysicsBlock> pb 
      = Teuchos::rcp(new panzer::PhysicsBlock(ipb,
                                              eBlocks[i], 
                                              default_integration_order,
                                              volume_cell_data,
                                              eqset_factory,
                                              gd,
                                              build_transient_support));

    // we can have more than one physics block, one per element block
    physicsBlocks.push_back(pb);
    }
    panzer::checkBCConsistency(eBlocks,sidesets,bcs);
  }


  // finish building mesh, set required field variables and mesh bulk data
  ////////////////////////////////////////////////////////////////////////

  for(size_t i=0; i<physicsBlocks.size(); ++i)
  {
    Teuchos::RCP<panzer::PhysicsBlock> pb = physicsBlocks[i]; // we are assuming only one physics block

    const std::vector<panzer::StrPureBasisPair> & blockFields = pb->getProvidedDOFs();

    // insert all fields into a set
    std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp> fieldNames;
    fieldNames.insert(blockFields.begin(),blockFields.end());

    // add basis to DOF manager: block specific
    std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp>::const_iterator fieldItr;
    for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr)
       mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
  }
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD); // this is where the mesh refinements are applied

  // GH: query the Pamgen mesh info here; this is where we will construct a region hierarchy
  Teuchos::RCP<percept::PerceptMesh> refinedMesh = mesh->getRefinedMesh();
  if(print_percept_mesh)
    refinedMesh->print_info(*out,"",1,true);

  // check if the parent information is stored
  std::vector<stk::mesh::EntityRank> ranks_to_be_deleted;
  //ranks_to_be_deleted.push_back(stk::topology::ELEMENT_RANK); // GH: we only care about elements for now
  ranks_to_be_deleted.push_back(refinedMesh->side_rank());
  if (refinedMesh->get_spatial_dim() == 3)
    ranks_to_be_deleted.push_back(refinedMesh->edge_rank());

  for (unsigned irank=0; irank < ranks_to_be_deleted.size()-1; irank++) // don't look for children of nodes
    {
      const stk::mesh::BucketVector & buckets = refinedMesh->get_bulk_data()->buckets( ranks_to_be_deleted[irank] );
      int npar=0;
      int nchild=0;
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
	{
	  stk::mesh::Bucket & bucket = **k ;

	  //std::cout << "New bucket" << std::endl;

	  const unsigned num_elements_in_bucket = bucket.size();
	  for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
	    {
	      stk::mesh::Entity element = bucket[iElement];
	      //	      std::cout << element << std::endl;
	      if (!refinedMesh->isParentElement(element, false))
		{
		  ++nchild;
		  refinedMesh->printParent(element,true);
		  //		  stk::mesh::Entity parent = refinedMesh->printParent(element,true);
		  //		  stk::mesh::Entity parent = 

    		  //std::cout << "Child index " << iElement << std::endl;
		}
	      else
		{
		  ++npar;
		  //std::cout << "Parent index " << iElement << std::endl;
		  //parents.insert(element);
		}
	    }
	}
      //std::cout << "Children on topology rank " << irank << ": " << nchild << std::endl;
      //std::cout << "Parents on topology rank " << irank << ": " << npar << std::endl;
    }
  

  // build DOF Manager and linear object factory
  /////////////////////////////////////////////////////////////

  tm = Teuchos::null;
  tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 2 - Build DOF Manager and Worksets")));
  // build the connection manager 
  const Teuchos::RCP<panzer::ConnManager> conn_manager 
    = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  panzer::DOFManagerFactory globalIndexerFactory;
  Teuchos::RCP<panzer::GlobalIndexer> dofManager 
    = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

  // construct some linear algebra object, build object to pass to evaluators
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
    = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,ST,LO,GO>(comm.getConst(),dofManager));

  /*
  ArrayRCP<GO> dofs_to_region;
  // Construct the DofToRegion map based on the Percept mesh and the DOF manager
  {
    std::vector<GO> gids;
    dofManager->getElementGIDs((LO) 0, gids);
    std::cout << "Element 0 GIDs: " << gids((int) 0); 
    for(GO i=1; i<gids.size(); ++i)
    {
      std::cout << ", " << gids(static_cast<LO>(i));
    }
    std::cout << std::endl;
  }
  */

  // build worksets
  ////////////////////////////////////////////////////////

  Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
    = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
  Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
    = Teuchos::rcp(new panzer::WorksetContainer);
  wkstContainer->setFactory(wkstFactory);
  for(size_t i=0;i<physicsBlocks.size();i++) 
    wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
  wkstContainer->setWorksetSize(workset_size);
  wkstContainer->setGlobalIndexer(dofManager);


  // Setup response library for checking the error in this manufactered solution
  ////////////////////////////////////////////////////////////////////////

  Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > errorResponseLibrary
    = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,dofManager,linObjFactory));

  {
    const int integration_order = 10;
    
    panzer::FunctionalResponse_Builder<int,int> builder;
    builder.comm = MPI_COMM_WORLD;
    builder.cubatureDegree = integration_order;
    builder.requiresCellIntegral = true;
    builder.quadPointField = "TEMPERATURE_L2_ERROR";

    errorResponseLibrary->addResponse("L2 Error",eBlocks,builder);

    /*    builder.comm = MPI_COMM_WORLD;
    builder.cubatureDegree = integration_order;
    builder.requiresCellIntegral = true;
    builder.quadPointField = "TEMPERATURE_H1_ERROR";

    errorResponseLibrary->addResponse("H1 Error",eBlocks,builder);
    */
  }


  // setup closure model
  /////////////////////////////////////////////////////////////
 
  // Add in the application specific closure model factory
  panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory; 
  Example::ClosureModelFactory_TemplateBuilder cm_builder;
  cm_factory.buildObjects(cm_builder);

  Teuchos::ParameterList closure_models("Closure Models");
  {
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<std::string>("Type","SIMPLE SOURCE"); // a constant source
    // SOURCE_TEMPERATURE field is required by the PoissonEquationSet
    // required for error calculation
    closure_models.sublist("solid").sublist("TEMPERATURE_L2_ERROR").set<std::string>("Type","L2 ERROR_CALC");
    closure_models.sublist("solid").sublist("TEMPERATURE_L2_ERROR").set<std::string>("Field A","TEMPERATURE");
    closure_models.sublist("solid").sublist("TEMPERATURE_L2_ERROR").set<std::string>("Field B","TEMPERATURE_EXACT");

    /*
    closure_models.sublist("solid").sublist("TEMPERATURE_H1_ERROR").set<std::string>("Type","H1 ERROR_CALC");
    closure_models.sublist("solid").sublist("TEMPERATURE_H1_ERROR").set<std::string>("Field A","TEMPERATURE");
    closure_models.sublist("solid").sublist("TEMPERATURE_H1_ERROR").set<std::string>("Field B","TEMPERATURE_EXACT");
    */
    closure_models.sublist("solid").sublist("TEMPERATURE_EXACT").set<std::string>("Type","TEMPERATURE_EXACT");
  }

  Teuchos::ParameterList user_data("User Data"); // user data can be empty here


  // setup field manager builder
  /////////////////////////////////////////////////////////////

  Teuchos::RCP<panzer::FieldManagerBuilder> fmb 
    = Teuchos::rcp(new panzer::FieldManagerBuilder);
  fmb->setWorksetContainer(wkstContainer);
  fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
  fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,
                            *linObjFactory,user_data);

  fmb->writeVolumeGraphvizDependencyFiles("Poisson", physicsBlocks);


  // setup assembly engine
  /////////////////////////////////////////////////////////////

  panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
  panzer::AssemblyEngine_TemplateBuilder builder(fmb,linObjFactory);
  ae_tm.buildObjects(builder);


  // Finalize construcition of STK writer response library
  /////////////////////////////////////////////////////////////
  {
    user_data.set<int>("Workset Size",workset_size);
    errorResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                                  cm_factory,
                                                  closure_models,
                                                  user_data);
  }


  // assemble linear system
  /////////////////////////////////////////////////////////////

  Teuchos::RCP<panzer::LinearObjContainer> ghostCont = linObjFactory->buildGhostedLinearObjContainer();
  Teuchos::RCP<panzer::LinearObjContainer> container = linObjFactory->buildLinearObjContainer();
  linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                            panzer::LinearObjContainer::F |
                                            panzer::LinearObjContainer::Mat,*ghostCont);
  linObjFactory->initializeContainer(panzer::LinearObjContainer::X |
                                     panzer::LinearObjContainer::F |
                                     panzer::LinearObjContainer::Mat,*container);
  ghostCont->initialize();
  container->initialize();

  panzer::AssemblyEngineInArgs input(ghostCont,container);
  input.alpha = 0;
  input.beta = 1;

  // evaluate physics: This does both the Jacobian and residual at once
  ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);


  // solve linear system
  /////////////////////////////////////////////////////////////


  // GH: convert this to our MueLu region solver
  // convert generic linear object container to tpetra container
  Teuchos::RCP<panzer::TpetraLinearObjContainer<ST,LO,GO> > tp_container
    = Teuchos::rcp_dynamic_cast<panzer::TpetraLinearObjContainer<ST,LO,GO> >(container);

  // Setup the linear solve: notice A is used directly 
  Belos::LinearProblem<ST,MV,OP> problem(tp_container->get_A(), tp_container->get_x(), tp_container->get_f());

  problem.setProblem();

  Teuchos::RCP<Teuchos::ParameterList> pl_belos = Teuchos::rcp(new Teuchos::ParameterList());
  pl_belos->set("Maximum Iterations", 1000);
  pl_belos->set("Convergence Tolerance", 1e-9);

  // build the solver  
  Belos::PseudoBlockGmresSolMgr<ST,MV,OP> solver(Teuchos::rcpFromRef(problem), pl_belos);

  // solve the linear system
  solver.solve();

  // scale by -1 since we solved a residual correction
  tp_container->get_x()->scale(-1.0);
  std::cout << "Solution local length: " << tp_container->get_x()->getLocalLength() << std::endl;
  std::cout << "Solution norm: " << tp_container->get_x()->norm2() << std::endl;
  
  // output data (optional)
  /////////////////////////////////////////////////////////////

  // write out solution to matrix
  {
    // redistribute solution vector to ghosted vector
    linObjFactory->globalToGhostContainer(*container,*ghostCont, panzer::TpetraLinearObjContainer<ST,LO,GO>::X 
					                       | panzer::TpetraLinearObjContainer<ST,LO,GO>::DxDt); 

    // get X Tpetra_Vector from ghosted container
    //Teuchos::RCP<panzer::TpetraLinearObjContainer<ST,LO,GO> > tp_ghostCont = Teuchos::rcp_dynamic_cast<panzer::TpetraLinearObjContainer<ST,LO,GO> >(ghostCont);
    //panzer_stk::write_solution_data(*dofManager,*mesh,*tp_ghostCont->get_x());

    // Due to multiple instances of this test being run at the same
    // time (one for each celltype and each order), we need to
    // differentiate output to prevent race conditions on output
    // file. Multiple runs for different mesh refinement levels for
    // the same celltype/order are ok as they are staged one after
    // another in the ADD_ADVANCED_TEST cmake macro.
    std::ostringstream filename;
    filename << "output_" << celltype << "_p" << discretization_order << ".exo";
    mesh->writeToExodus(filename.str());
  }

  // compute error norm
  /////////////////////////////////////////////////////////////

  if(true)
  {
    panzer::AssemblyEngineInArgs respInput(ghostCont,container);
    respInput.alpha = 0;
    respInput.beta = 1;

    Teuchos::RCP<panzer::ResponseBase> l2_resp = errorResponseLibrary->getResponse<panzer::Traits::Residual>("L2 Error");
    Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > l2_resp_func 
      = Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(l2_resp);
    Teuchos::RCP<Thyra::VectorBase<double> > l2_respVec = Thyra::createMember(l2_resp_func->getVectorSpace());
    l2_resp_func->setVector(l2_respVec);

    /*
    Teuchos::RCP<panzer::ResponseBase> h1_resp = errorResponseLibrary->getResponse<panzer::Traits::Residual>("H1 Error");
    Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > h1_resp_func 
      = Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(h1_resp);
    Teuchos::RCP<Thyra::VectorBase<double> > h1_respVec = Thyra::createMember(h1_resp_func->getVectorSpace());
    h1_resp_func->setVector(h1_respVec);
    */

    errorResponseLibrary->addResponsesToInArgs<panzer::Traits::Residual>(respInput);
    errorResponseLibrary->evaluate<panzer::Traits::Residual>(respInput);

    if(myRank == 0)
    {
      *out << "This is the Basis Order" << std::endl;
      *out << "Basis Order = " << discretization_order << std::endl;
      *out << "This is the L2 Error" << std::endl;
      *out << "L2 Error = " << sqrt(l2_resp_func->value) << std::endl;
      //*out << "This is the H1 Error" << std::endl;
      //*out << "H1 Error = " << sqrt(h1_resp_func->value) << std::endl;
    }
  }

  tm = Teuchos::null;
  globalTimeMonitor = Teuchos::null;

  if (showTimerSummary)
  {
    RCP<ParameterList> reportParams = rcp(new ParameterList);
    const std::string filter = "";
    if (useStackedTimer) {
      Teuchos::StackedTimer::OutputOptions options;
      options.output_fraction = options.output_histogram = options.output_minmax = true;
      stacked_timer->report(*out, comm, options);
    } else {
      std::ios_base::fmtflags ff(out->flags());
      Teuchos::TimeMonitor::report(comm.ptr(), *out, filter, reportParams);
      *out << std::setiosflags(ff);
    }
  }

  } // Kokkos scope
  Kokkos::finalize();

  return 0;
} // main
