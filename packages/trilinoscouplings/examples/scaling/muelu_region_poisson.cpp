// Standard headers
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unistd.h>

// TrilinosCouplings headers
#include "TrilinosCouplings_config.h"

// Teuchos headers
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_YamlParameterListHelpers.hpp"

// Belos headers
#include "BelosBiCGStabSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMueLuAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"
#ifdef HAVE_MUELU_TPETRA
#include <BelosTpetraAdapter.hpp>
#endif

// MueLu headers
#include "MueLu.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MutuallyExclusiveTime.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif

#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif

// MueLu and Xpetra Tpetra stack
#ifdef HAVE_MUELU_TPETRA
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <KokkosBlas1_abs.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix.hpp>
#include <Tpetra_computeRowAndColumnOneNorms.hpp>
#endif

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
#include <Amesos2_config.h>
#include <Amesos2.hpp>
#endif

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
#include "Panzer_TpetraLinearObjFactory.hpp"

// Percept headers
#include <percept/PerceptMesh.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>

// Tpetra headers
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_FECrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Include factories for boundary conditions and other Panzer setup
// Most of which is taken from PoissonExample in Panzer_STK
#include "muelu_region_poisson.hpp"


// need this for debugging purposes... technically any class that has a size()
// method and [] operator that returns a type that can be streamed works,
// but the stream needs to be able to use std::endl; (could replace with \n)
#define DUMPSTDVECTOR(vector, stream) \
  stream << #vector << ".size() = " << vector.size() << std::endl; \
  stream << #vector << "        = (" << vector[0]; \
  for(unsigned int i=1; i < vector.size(); ++i) \
    stream << ", " << vector[i]; \
  stream << ")" << std::endl;



int main(int argc, char *argv[]) {

  // The following typedefs are used so that the two codes will work together.
  // TODO: change everything to use the same types to avoid the following silly code...
  // i.e. check what Drekar does first before making changes
  // Panzer types
  using ST = double;
  using LO = panzer::LocalOrdinal;
  using GO = panzer::GlobalOrdinal;
  using NT = panzer::TpetraNodeType;
  using OP = Tpetra::Operator<ST,LO,GO,NT>;
  using MV = Tpetra::MultiVector<ST,LO,GO,NT>; 

  // MueLu types
  using Scalar = ST; 
  using LocalOrdinal = LO; 
  using GlobalOrdinal = GO;
  using Node = NT;

#include <MueLu_UseShortNames.hpp>

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize (argc, argv);
  Tpetra::initialize (&argc, &argv);
  { // Parallel initialization scope


    ////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// SETUP //////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////


    // Setup output stream, MPI, and grab info
    Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
    out.setOutputToRootOnly(0);

    //Teuchos::RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    //Teuchos::FancyOStream& out = *fancy;
    //out.setOutputToRootOnly(0); // use out on rank 0

    Teuchos::RCP<Teuchos::FancyOStream> fancydebug = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& debug = *fancydebug; // use on all ranks

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    const int numRanks = comm->getSize();
    const int myRank = comm->getRank();

    out << "Running TrilinosCouplings region multigrid driver on " << numRanks << " ranks... \n";

    // Parse command line arguments
    Teuchos::CommandLineProcessor clp(false);
    std::string exodusFileName        = "";                  clp.setOption("exodus-mesh",           &exodusFileName,          "Exodus hex mesh filename (overrides a pamgen-mesh if both specified)");
    std::string pamgenFileName        = "cylinder.rtp";      clp.setOption("pamgen-mesh",           &pamgenFileName,          "Pamgen hex mesh filename");
    std::string xmlFileName           = "";                  clp.setOption("xml",                   &xmlFileName,             "MueLu parameters from an xml file");
    std::string yamlFileName          = "";                  clp.setOption("yaml",                  &yamlFileName,            "MueLu parameters from a yaml file");
    int mesh_refinements              = 1;                   clp.setOption("mesh-refinements",      &mesh_refinements,        "Uniform mesh refinements");
    bool delete_parent_elements       = false;               clp.setOption("delete-parent-elements", "keep-parent-elements", &delete_parent_elements,"Save the parent elements in the perceptMesh");

    // Multigrid options
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
    int         rebalanceNumPartitions = -1;                 clp.setOption("numPartitions",         &rebalanceNumPartitions, "number of partitions for rebalancing the coarse grid AMG solve");
    std::string coarseSolverType      = "direct";            clp.setOption("coarseSolverType",      &coarseSolverType,      "Type of solver for (composite) coarse level operator (smoother | direct | amg)");
    std::string unstructured          = "{}";                clp.setOption("unstructured",          &unstructured,          "List of ranks to be treated as unstructured, e.g. {0, 2, 5}");
    std::string coarseAmgXmlFile      = "";                  clp.setOption("coarseAmgXml",          &coarseAmgXmlFile,      "Read parameters for AMG as coarse level solve from this xml file.");
    std::string coarseSmootherXMLFile = "";                  clp.setOption("coarseSmootherXML",     &coarseSmootherXMLFile, "File containing the parameters to use with the coarse level smoother.");
    int  cacheSize = 0;                                      clp.setOption("cachesize",               &cacheSize,           "cache size (in KB)"); // what does this do?
    std::string cycleType = "V";                             clp.setOption("cycleType", &cycleType, "Multigrid cycle type. Possible values: V, W.");
#ifdef HAVE_MUELU_TPETRA
    std::string equilibrate = "no" ;                         clp.setOption("equilibrate",           &equilibrate,           "equilibrate the system (no | diag | 1-norm)");
#endif
#ifdef HAVE_MUELU_CUDA
    bool profileSetup = false;                               clp.setOption("cuda-profile-setup", "no-cuda-profile-setup", &profileSetup, "enable CUDA profiling for setup");
    bool profileSolve = false;                               clp.setOption("cuda-profile-solve", "no-cuda-profile-solve", &profileSolve, "enable CUDA profiling for solve");
#endif

    // debug options
    bool print_percept_mesh           = false;               clp.setOption("print-percept-mesh", "no-print-percept-mesh", &print_percept_mesh, "Calls perceptMesh's print_info routine");
    bool print_debug_info             = false;               clp.setOption("print-debug-info", "no-print-debug-info", &print_debug_info, "Print more debugging information");
    bool dump_element_vertices        = false;               clp.setOption("dump-element-vertices", "no-dump-element-vertices", &dump_element_vertices, "Dump the panzer_stk mesh vertices, element-by-element");

    // timer options
    bool useStackedTimer              = false;               clp.setOption("stacked-timer","no-stacked-timer", &useStackedTimer, "use stacked timer");
    bool showTimerSummary             = false;               clp.setOption("show-timer-summary", "no-show-timer-summary", &showTimerSummary, "Switch on/off the timer summary at the end of the run.");

    TEUCHOS_ASSERT(mesh_refinements >= 0); // temporarily do this instead of typing as unsigned int to get around the expected 7 arguments error for clp.setOption(...unsigned int...)

    clp.recogniseAllOptions(true);
    switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName != "" && yamlFileName != "", std::runtime_error,
                               "Cannot provide both xml and yaml input files");

    // get xml file from command line if provided, otherwise use default
    std::string  xmlSolverInFileName(xmlFileName);

    // Read xml file into parameter list
    Teuchos::ParameterList inputSolverList;

    if(xmlSolverInFileName.length()) {
        out << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n" << std::endl;
      Teuchos::updateParametersFromXmlFile (xmlSolverInFileName, Teuchos::ptr(&inputSolverList));
    }
    else {
      out << "Using default solver values ..." << std::endl;
    }


    ////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////// MESH AND WORKSETS /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////


    Teuchos::RCP<Teuchos::Time> meshTimer = Teuchos::TimeMonitor::getNewCounter("Step 1: Mesh generation");
    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
    if(useStackedTimer)
      stacked_timer = rcp(new Teuchos::StackedTimer("MueLu_Driver"));
    Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
    RCP<Teuchos::TimeMonitor> globalTimeMonitor = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: S - Global Time")));
    RCP<Teuchos::TimeMonitor> tm                = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 1 - Build Mesh and Assign Physics")));

    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
    Teuchos::RCP<Teuchos::ParameterList> mesh_pl = Teuchos::rcp(new Teuchos::ParameterList);

    // TODO: due to #8475, we have this. we will clean up Panzer and this code once we verify
    // the functionality of the full HHG scheme obtained by using Percept
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory_percept;
    Teuchos::RCP<Teuchos::ParameterList> mesh_pl_percept = Teuchos::rcp(new Teuchos::ParameterList);

    bool use_exodus_mesh = exodusFileName.length() > 0;
    bool use_pamgen_mesh = !use_exodus_mesh;
    bool use_percept = mesh_refinements > 0;

    if(use_exodus_mesh)
    {
      // set the filename and type
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory);
      mesh_pl->set("File Name", exodusFileName);
      mesh_pl->set("File Type", "Exodus");
    }
    else if(use_pamgen_mesh)
    {
      // set the filename and type
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory);
      mesh_pl->set("File Name", pamgenFileName);
      mesh_pl->set("File Type", "Pamgen");

      if(use_percept)
      {
        mesh_pl->set("Levels of Uniform Refinement", mesh_refinements); // this multiplies the number of elements by 2^(dimension*level)

        mesh_factory_percept = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory);
        mesh_pl_percept->set("File Name", pamgenFileName);
        mesh_pl_percept->set("File Type", "Pamgen");
        mesh_pl_percept->set("Levels of Uniform Refinement", mesh_refinements); // this multiplies the number of elements by 2^(dimension*level)
        mesh_pl_percept->set("Keep Percept Data", true); // this is necessary to gather mesh hierarchy information
        mesh_pl_percept->set("Keep Percept Parent Elements", !delete_parent_elements); // this is necessary to gather mesh hierarchy information
      }
    }
    else
      throw std::runtime_error("no mesh file name found!");
    

    mesh_factory_percept->setParameterList(mesh_pl_percept);
    Teuchos::RCP<panzer_stk::STK_Interface> parent_data_mesh;
    parent_data_mesh = mesh_factory_percept->buildMesh(MPI_COMM_WORLD);

    // set the parameters
    mesh_factory->setParameterList(mesh_pl);

    // build the mesh
    Teuchos::RCP<panzer_stk::STK_Interface> mesh;
    mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);

    // setup the physics block
    Teuchos::RCP<Example::EquationSetFactory> eqset_factory = Teuchos::rcp(new Example::EquationSetFactory);
    Example::BCStrategyFactory bc_factory;
    const std::size_t workset_size = 10; // TODO: this may be much larger in practice. experiment with it.
    const int discretization_order = 1;

    // grab the number and names of mesh blocks
    std::vector<std::string> eBlocks;
    mesh->getElementBlockNames(eBlocks);
    std::vector<bool> unstructured_eBlocks(eBlocks.size(), false);
    // TODO: set unstructured blocks based on some sort of input information; for example, using the Exodus ex_get_var* functions

    // grab the number and name of nodesets
    std::vector<std::string> nodesets;
    mesh->getNodesetNames(nodesets);

    // TODO: Dirichlet boundaries on top and bottom of cylinder

    // grab the number and names of sidesets
    std::vector<std::string> sidesets;
    mesh->getSidesetNames(sidesets);

    // create a physics blocks parameter list
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;

    // set physics and boundary conditions on each block
    {
      bool build_transient_support = false;

      const int integration_order = 10;
      Teuchos::ParameterList& p = ipb->sublist("Poisson Physics");
      p.set("Type","Poisson");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",discretization_order);
      p.set("Integration Order",integration_order);

      // TODO: double-check. this assumes we impose Dirichlet BCs on all boundaries of all physics blocks
      // It may potentially assign Dirichlet BCs to internal block boundaries, which is undesirable
      for(size_t i=0; i<eBlocks.size(); ++i)
      {
        for(size_t j=0; j<sidesets.size(); ++j)
        {
          std::size_t bc_id = j;
          panzer::BCType bctype = panzer::BCT_Dirichlet;
          std::string sideset_id = sidesets[j];
          std::string element_block_id = eBlocks[i];
          std::string dof_name = "TEMPERATURE";
          std::string strategy = "Constant";
          double value = 0.0;
	  Teuchos::ParameterList pbc; // this is how the official example does it, so I'll leave it alone for now
          pbc.set("Value",value);
          panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name,
                        strategy, pbc);
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
    }
    panzer::checkBCConsistency(eBlocks,sidesets,bcs);


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

    unsigned int num_dimensions = mesh->getDimension();

    // build DOF Manager and linear object factory
    /////////////////////////////////////////////////////////////

    tm = Teuchos::null;
    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 2 - Build DOF Manager and Worksets")));
    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory globalIndexerFactory;
    Teuchos::RCP<panzer::GlobalIndexer> dofManager = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

    // construct some linear algebra object, build object to pass to evaluators
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,ST,LO,GO>(comm.getConst(),dofManager));

    // build worksets
    ////////////////////////////////////////////////////////

    // build STK workset factory and attach it to a workset container (uses lazy evaluation)
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh));
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physicsBlocks.size();i++)
      wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);
    wkstContainer->setGlobalIndexer(dofManager);


    ////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// CONSTRUCT REGIONS //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////


    // The code in this section assumes that a region hierarchy can be established with Percept. 
    // In the case where a region hierarchy is constructed from an exodus data input, for example, 
    // this implementation will need to be updated.
    // TODO: Assign MPI rank p to region p and collect element IDs. If this region is assigned via
    // percept, reorder the element IDs lexicographically using the utility in the header. 
    // Then collect mesh->identifier(node) for the nodes in lexicographic order for region p, 
    // and put those in quasiRegionGIDs. Coordinates should be able to be extracted from the 
    // stk::mesh::entity node as well.

    tm = Teuchos::null;
    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 3 - Setup Region Information")));

    Teuchos::RCP<stk::mesh::BulkData> bulk_data = mesh->getBulkData();
    Teuchos::RCP<stk::mesh::MetaData> meta_data = mesh->getMetaData();

    std::vector<stk::mesh::Entity> nodes;
    bulk_data->get_entities(mesh->getNodeRank(),meta_data->locally_owned_part(),nodes);
    std::vector<stk::mesh::Entity> cells;
    bulk_data->get_entities(mesh->getElementRank(),meta_data->locally_owned_part(),cells);

    const unsigned int region_cells_per_dim = (1 << mesh_refinements);
    const unsigned int region_nodes_per_dim = region_cells_per_dim + 1;
    const unsigned int region_nodes_per_proc = region_nodes_per_dim*region_nodes_per_dim*region_nodes_per_dim;

    // Do a little bit of output here, regardless of verbosity
    out << "Driver parameters: " << std::endl;
    out << "          Dimension = " << num_dimensions << std::endl;
    out << "  Refinement Levels = " << mesh_refinements << std::endl;
    out << "   Region cells/dim = " << region_cells_per_dim << std::endl;
    out << "   Region nodes/dim = " << region_nodes_per_dim << std::endl;
    debug << "       Rank " << myRank << " cells = " << cells.size() << std::endl;
    debug << "       Rank " << myRank << " nodes = " << nodes.size() << std::endl;
    out << std::endl;

    comm->barrier();

    const unsigned int children_per_element = 1 << (num_dimensions*mesh_refinements);
    if(print_debug_info)
      out << "Number of mesh children = " << children_per_element << std::endl;

    // initialize data here that we will use for the region MG solver
    std::vector<GO> local_child_element_gids; // these don't always start at 0, and changes I'm making to Panzer keep changing this, so I'll store them for now
    std::vector<GO> local_child_element_region_gids;
    std::vector<GO> local_element_gids;
    Array<GO> x_coords(region_nodes_per_proc, 0);
    Array<GO> y_coords(region_nodes_per_proc, 0);
    Array<GO> z_coords(region_nodes_per_proc, 0);
    Array<GO>  quasiRegionGIDs;
    Array<GO>  quasiRegionCoordGIDs;

    quasiRegionGIDs.resize(region_nodes_per_proc);

    // if we use Percept, this is the strategy to follow.
    // this may be turned into its own function later
    if(use_percept)
    {
      // get the Percept mesh from Panzer
      Teuchos::RCP<percept::PerceptMesh> percept_mesh = parent_data_mesh->getRefinedMesh();
      if(print_percept_mesh)
        percept_mesh->print_info(out,"",1,true);

      // ids are linear within stk, but we need an offset because the original mesh info comes first
//      size_t node_id_start = 0;
//      {
//        const stk::mesh::BucketVector & local_buckets = percept_mesh->get_bulk_data()->get_buckets(stk::topology::ELEM_RANK,percept_mesh->get_fem_meta_data()->locally_owned_part());
//        //const stk::mesh::BucketVector & buckets = percept_mesh->get_bulk_data()->buckets(percept_mesh->node_rank());
//        stk::mesh::Bucket & bucket = **local_buckets.begin() ;
//        node_id_start = percept_mesh->id(bucket[0]);
//        if(print_debug_info)
//          debug << "Starting node id = " << node_id_start << std::endl;
//      }
//
//      size_t elem_id_start = 0;
//      {
//        const stk::mesh::BucketVector & local_buckets = percept_mesh->get_bulk_data()->get_buckets(stk::topology::ELEM_RANK,percept_mesh->get_fem_meta_data()->locally_owned_part());
//        //const stk::mesh::BucketVector & buckets = percept_mesh->get_bulk_data()->buckets(percept_mesh->element_rank());
//        stk::mesh::Bucket & bucket = **local_buckets.begin() ;
//        elem_id_start = percept_mesh->id(bucket[0]);
//        if(print_debug_info)
//          debug << "Starting element id = " << elem_id_start << std::endl;
//      }
      //panzer_stk::workset_utils::getIdsAndVertices


      // grab the region information from the mesh that keeps parent elements
      {
        // count parents and children
        int npar=0;
        int nchild=0;

        const stk::mesh::BucketVector & local_buckets = percept_mesh->get_bulk_data()->get_buckets(stk::topology::ELEM_RANK,percept_mesh->get_fem_meta_data()->locally_owned_part());
        for (stk::mesh::BucketVector::const_iterator k = local_buckets.begin(); k != local_buckets.end(); ++k)
        {
          const stk::mesh::Bucket & bucket = **k ;
          if(print_debug_info)
            debug << "New bucket" << std::endl;

          const unsigned int num_elements_in_bucket = bucket.size();
          for (unsigned int iElement = 0; iElement < num_elements_in_bucket; iElement++)
          {
            const stk::mesh::Entity element = bucket[iElement];
            if (!percept_mesh->isParentElement(element, false))
            {
              ++nchild;

              // this is the important part here. take the id of the element and the id of the element's root
              local_child_element_gids.push_back(percept_mesh->id(element));
              local_child_element_region_gids.push_back(percept_mesh->id(percept_mesh->rootOfTree(element)));

              if(print_debug_info)
                debug << "Percept Element = " << percept_mesh->id(element) << std::endl;

              const percept::MyPairIterRelation elem_nodes ( *percept_mesh, element,  stk::topology::NODE_RANK);

              for (unsigned int i_node = 0; i_node < elem_nodes.size(); i_node++)
              {
                const stk::mesh::Entity node = elem_nodes[i_node].entity();
                //local_node_gids.push_back(percept_mesh->id(node));

                if(print_debug_info)
                  debug << "Stk Node = " << percept_mesh->id(node) << std::endl;
              }
            }
            else
            {
              if(print_debug_info)
                debug << "p = " << myRank << ", parent = " << percept_mesh->id(element) << std::endl;
            }
          }
        }

        if(print_debug_info)
        {
          for(unsigned int i=0; i<local_child_element_gids.size(); ++i)
          {
            out << "child= " << local_child_element_gids[i] << " parent= " << local_child_element_region_gids[i] << std::endl;
          }
        }
      }

      // however, the mesh that keeps parent elements does not assemble finite
      // element data properly; therefore, we must now convert from element IDs
      // on the mesh with parent elements to element IDs on the mesh without
      // parent elements.
      {
        const stk::mesh::BucketVector & local_buckets = mesh->getBulkData()->get_buckets(stk::topology::ELEM_RANK,mesh->getMetaData()->locally_owned_part());
        for (stk::mesh::BucketVector::const_iterator k = local_buckets.begin(); k != local_buckets.end(); ++k)
        {
          const stk::mesh::Bucket & bucket = **k;

          const unsigned int num_elements_in_bucket = bucket.size();
          for (unsigned int iElement = 0; iElement < num_elements_in_bucket; iElement++)
          {
            const stk::mesh::Entity element = bucket[iElement];
            local_element_gids.push_back(mesh->getBulkData()->identifier(element));
          }
        }
      }

      // grab the Percept element renumbering
      sleep(myRank);
      const std::vector<unsigned int> percept_lexicographic_elements = renumberPerceptCellsToLexicographic(mesh_refinements);
      DUMPSTDVECTOR(percept_lexicographic_elements, debug);
      DUMPSTDVECTOR(local_child_element_gids, debug);
      DUMPSTDVECTOR(local_element_gids, debug);
      comm->barrier();

      // I think this is no longer necessary. Since elements are locally ordered,
      // percept_lexicographic_elements is exactly the local indexing we need.
      // // apply the renumbering
      //std::vector<unsigned int> renumbered_local_element_gids;
      //for(unsigned int i=0; i<percept_lexicographic_elements.size(); ++i)
      //  renumbered_local_element_gids.push_back(local_element_gids[percept_lexicographic_elements[i]]);
      // // make sure we have the correct number of elements on this rank
      //TEUCHOS_ASSERT(renumbered_local_element_gids.size() == children_per_element);

      // once we have the elements in lexicographic order, we have to loop through 
      // one last time and grab the vertices in lexicographic order too
      // warning: this is not very performative as far as access patterns, etc go
      // I'll revisit performance once it's working. in principle, this isn't a
      // huge computational sink since we're likely not making regions too huge...
      // we'll see what the timers say
      {
        const stk::mesh::FieldBase *coordinatesField = mesh->getMetaData()->get_field(stk::topology::NODE_RANK, "coordinates");
        const stk::mesh::BucketVector & local_buckets = mesh->getBulkData()->get_buckets(stk::topology::ELEM_RANK,mesh->getMetaData()->locally_owned_part());
        for (stk::mesh::BucketVector::const_iterator iBucket = local_buckets.begin(); iBucket != local_buckets.end(); ++iBucket)
        {
          stk::mesh::Bucket & bucket = **iBucket;

          const unsigned int num_elements_in_bucket = bucket.size();
          debug << "bucket size = " << num_elements_in_bucket << std::endl;

          for (unsigned int k = 0; k < region_cells_per_dim; k++) // z direction
          {
            for (unsigned int j = 0; j < region_cells_per_dim; j++) // y direction
            {
              for (unsigned int i = 0; i < region_cells_per_dim; i++) // x direction
              {
                const unsigned int iElement = percept_lexicographic_elements[region_cells_per_dim*region_cells_per_dim*k + region_cells_per_dim*j + i];
                const stk::mesh::Entity element = bucket[iElement];

                // GH: commenting for now... I don't think I need this
                // const unsigned int size = mesh->getBulkData()->num_connectivity(element, stk::topology::NODE_RANK);

                // these are fine... this never has over 8 nodes per element
                // this is what Percept::MyPairIterRelation does under the hood. I'm stealing it since this mesh doesn't have access to Percept
                const stk::mesh::Entity *nodes = mesh->getBulkData()->begin(element, stk::topology::NODE_RANK);
                const stk::mesh::ConnectivityOrdinal *ordinals = mesh->getBulkData()->begin_ordinals(element, stk::topology::NODE_RANK);

                double *nodeCoords[8];
                for (unsigned int i=0; i<8; ++i)
                  nodeCoords[i] = static_cast<double *>(stk::mesh::field_data(*coordinatesField, nodes[i]));

                // subtracting 1 everywhere since STK numbering starts at 1
                quasiRegionGIDs[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = mesh->getBulkData()->identifier(nodes[0]) - 1;
                x_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[0][0];
                y_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[0][1];
                z_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[0][2];

                // if we're on the +z side of a region
                if(k==region_cells_per_dim-1)
                {
                  quasiRegionGIDs[region_nodes_per_dim*region_nodes_per_dim*(k+1) + region_nodes_per_dim*j + i] = mesh->getBulkData()->identifier(nodes[4]) - 1;
                  x_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[4][0];
                  y_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[4][1];
                  z_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[4][2];
                  if(j==region_cells_per_dim-1)
                  {
                    quasiRegionGIDs[region_nodes_per_dim*region_nodes_per_dim*(k+1) + region_nodes_per_dim*(j+1) + i] = mesh->getBulkData()->identifier(nodes[6]) - 1;
                    x_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[6][0];
                    y_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[6][1];
                    z_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[6][2];
                    if(i==region_cells_per_dim-1)
                    {
                      quasiRegionGIDs[region_nodes_per_dim*region_nodes_per_dim*(k+1) + region_nodes_per_dim*(j+1) + i + 1] = mesh->getBulkData()->identifier(nodes[7]) - 1;
                      x_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[7][0];
                      y_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[7][1];
                      z_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[7][2];
                    }
                  }
                  if(i==region_cells_per_dim-1)
                  {
                    quasiRegionGIDs[region_nodes_per_dim*region_nodes_per_dim*(k+1) + region_nodes_per_dim*j + i + 1] = mesh->getBulkData()->identifier(nodes[5]) - 1;
                    x_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[5][0];
                    y_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[5][1];
                    z_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[5][2];
                  }
                }
                // if we're on the +y side of a region
                if(j==region_cells_per_dim-1)
                {
                  quasiRegionGIDs[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*(j+1) + i] = mesh->getBulkData()->identifier(nodes[3]) - 1;
                  x_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[3][0];
                  y_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[3][1];
                  z_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[3][2];
                  if(i==region_cells_per_dim-1)
                  {
                    quasiRegionGIDs[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*(j+1) + i + 1] = mesh->getBulkData()->identifier(nodes[2]) - 1;
                    x_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[2][0];
                    y_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[2][1];
                    z_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[2][2];
                  }
                }
                // if we're on the +x side of a region
                if(i==region_cells_per_dim-1)
                {
                  quasiRegionGIDs[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i + 1] = mesh->getBulkData()->identifier(nodes[1]) - 1;
                  x_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[1][0];
                  y_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[1][1];
                  z_coords[region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i] = nodeCoords[1][2];
                }
              }
            }
          }
        }
      }

      debug << "quasiRegionGIDs = " << quasiRegionGIDs << std::endl;


    } // if use percept

    {
      std::vector<stk::mesh::Entity> elements;
      Kokkos::DynRankView<double,PHX::Device> vertices;
      std::vector<std::size_t> localIds;

      panzer_stk::workset_utils::getIdsAndVertices(*mesh,eBlocks[0],localIds,vertices); // TODO: in Matthias' case, this is likely eBlock[myRank]

      if(dump_element_vertices)
      {
        sleep(myRank);
        for(unsigned int ielem=0; ielem<vertices.extent(0); ++ielem)
          for(unsigned int ivert=0; ivert<vertices.extent(1); ++ivert)
          {
            out << "rank p=" << myRank << " element=" << ielem << " vertex=" << ivert << "   (" << vertices(ielem,ivert,0);
            for(unsigned int idim=1; idim<vertices.extent(2); ++idim) // fenceposting the output
              out << ", " << vertices(ielem,ivert,idim);
            out << ")" << std::endl;
          }
        return 0;
      }
    }

    // Setup response library for checking the error in this manufactured solution
    ////////////////////////////////////////////////////////////////////////

    tm = Teuchos::null;
    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 4 - Other Panzer Setup")));
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > errorResponseLibrary = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,dofManager,linObjFactory));

    {
      // must set a higher integration order to avoid false superconvergence
      const int integration_order = 5;

      panzer::FunctionalResponse_Builder<int,int> builder;
      builder.comm = MPI_COMM_WORLD;
      builder.cubatureDegree = integration_order;
      builder.requiresCellIntegral = true;
      builder.quadPointField = "TEMPERATURE_L2_ERROR";

      errorResponseLibrary->addResponse("L2 Error",eBlocks,builder);

      // TODO: uncomment the H1 errors once things look correct in the L2 norm
//      builder.comm = MPI_COMM_WORLD;
//      builder.cubatureDegree = integration_order;
//      builder.requiresCellIntegral = true;
//      builder.quadPointField = "TEMPERATURE_H1_ERROR";
//
//      errorResponseLibrary->addResponse("H1 Error",eBlocks,builder);
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

      // TODO: uncomment the H1 errors once things look correct in the L2 norm
//      closure_models.sublist("solid").sublist("TEMPERATURE_H1_ERROR").set<std::string>("Type","H1 ERROR_CALC");
//      closure_models.sublist("solid").sublist("TEMPERATURE_H1_ERROR").set<std::string>("Field A","TEMPERATURE");
//      closure_models.sublist("solid").sublist("TEMPERATURE_H1_ERROR").set<std::string>("Field B","TEMPERATURE_EXACT");
      closure_models.sublist("solid").sublist("TEMPERATURE_EXACT").set<std::string>("Type","TEMPERATURE_EXACT");
    }

    Teuchos::ParameterList user_data("User Data"); // user data can be empty here


    // setup field manager builder
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);
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


    // Finalize construction of STK writer response library
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


    ////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// LINEAR SOLVER ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////


    // TODO: this goes away once we finish getting the runtime errors in the region driver section sorted
    tm = Teuchos::null;
    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 5 - Linear Solver")));

    // convert generic linear object container to tpetra container
    Teuchos::RCP<panzer::TpetraLinearObjContainer<ST,LO,GO> > tp_container = Teuchos::rcp_dynamic_cast<panzer::TpetraLinearObjContainer<ST,LO,GO> >(container);
    Teuchos::RCP<MueLu::TpetraOperator<ST,LO,GO,NT> > mueLuPreconditioner;

    if(xmlFileName.size())
    {
      mueLuPreconditioner = MueLu::CreateTpetraPreconditioner(Teuchos::rcp_dynamic_cast<Tpetra::Operator<ST,LO,GO,NT> >(tp_container->get_A()), xmlFileName);
    }
    else
    {
      Teuchos::ParameterList mueLuParamList;
      if(print_debug_info)
      {
        mueLuParamList.set("verbosity", "high");
      }
      else
      {
        mueLuParamList.set("verbosity", "low");
      }
      mueLuParamList.set("max levels", 3);
      mueLuParamList.set("coarse: max size", 10);
      mueLuParamList.set("multigrid algorithm", "sa");
      mueLuPreconditioner = MueLu::CreateTpetraPreconditioner(Teuchos::rcp_dynamic_cast<Tpetra::Operator<ST,LO,GO,NT> >(tp_container->get_A()), mueLuParamList);
    }

    // Setup the linear solve
    Belos::LinearProblem<ST,MV,OP> problem(tp_container->get_A(), tp_container->get_x(), tp_container->get_f());
    problem.setLeftPrec(mueLuPreconditioner);
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
    if(print_debug_info)
    {
      debug << "Solution local length: " << tp_container->get_x()->getLocalLength() << std::endl;
      out << "Solution norm: " << tp_container->get_x()->norm2() << std::endl;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// REGION DRIVER ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::ArrayRCP;
      using Teuchos::TimeMonitor;
      using Teuchos::ParameterList;

      using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;

      // =========================================================================
      // Convenient definitions
      // =========================================================================
      using STS = Teuchos::ScalarTraits<SC>;
      SC zero = STS::zero(), one = STS::one();
      using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
      using real_type = typename STS::coordinateType;


      ParameterList paramList;
      //auto inst = xpetraParameters.GetInstantiation();

      if (yamlFileName != "") {
        Teuchos::updateParametersFromYamlFileAndBroadcast(yamlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
      } else {
        //if (inst == Xpetra::COMPLEX_INT_INT)
        //  xmlFileName = (xmlFileName != "" ? xmlFileName : "muelu_region_poisson_input-complex.xml");
        //else
          xmlFileName = (xmlFileName != "" ? xmlFileName : "muelu_region_poisson_input.xml");
        Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
      }

      Array<RCP<Teuchos::ParameterList> > smootherParams(1); //TODO: this is good, resized to numlevel
      smootherParams[0] = rcp(new Teuchos::ParameterList());
      smootherParams[0]->set("smoother: type",    smootherType);
      smootherParams[0]->set("smoother: sweeps",  smootherIts);
      smootherParams[0]->set("smoother: damping", smootherDamp);
      smootherParams[0]->set("smoother: Chebyshev eigRatio", smootherChebyEigRatio);
      smootherParams[0]->set("smoother: Chebyshev boost factor", smootherChebyBoostFactor);

      bool useUnstructured = false;
      Array<LO> unstructuredRanks = Teuchos::fromStringToArray<LO>(unstructured);
      for(int idx = 0; idx < unstructuredRanks.size(); ++idx) {
        if(unstructuredRanks[idx] == myRank) {useUnstructured = true;}
      }
      

      // =========================================================================
      // Problem construction
      // =========================================================================
#ifdef HAVE_MUELU_OPENMP
      //std::string node_name = Node::name();
      //if(!comm->getRank() && !node_name.compare("OpenMP/Wrapper"))
      //  galeriStream<<"OpenMP Max Threads = "<<omp_get_max_threads()<<std::endl;
#endif


      comm->barrier();
      Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
      if(useStackedTimer)
        stacked_timer = rcp(new Teuchos::StackedTimer("MueLu_Driver"));
      Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
      RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: S - Global Time")));
      RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Build Composite Matrix")));


      RCP<Matrix> A;
      RCP<Map>    nodeMap, dofMap;
      RCP<Vector> X, B;
      RCP<MultiVector>           nullspace;
      RCP<RealValuedMultiVector> coordinates;

      const int numDofsPerNode = 1;
      Teuchos::Array<LO> lNodesPerDim(3);
      if(use_percept)
      {
        lNodesPerDim[0] = region_nodes_per_dim;
        lNodesPerDim[1] = region_nodes_per_dim;
        lNodesPerDim[2] = region_nodes_per_dim;
      }

      // Create map and coordinates from Panzer's dofManager
      // TODO: get coordinates from Panzer
      RCP<panzer::TpetraLinearObjFactory<panzer::Traits,ST,LO,GO> > tp_object_factory = rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,ST,LO,GO>(comm, dofManager));
      nodeMap = Teuchos::rcp_dynamic_cast<Map>(tp_object_factory->getMap());
      //  coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double,LO,GO,Map,RealValuedMultiVector>("3D", nodeMap, galeriList);
      // }
      
      dofMap = Xpetra::MapFactory<LO,GO,Node>::Build(nodeMap, numDofsPerNode);
      A = Teuchos::rcp_dynamic_cast<Matrix>(tp_container->get_A());
      //nullspace = Pr->BuildNullspace(); // TODO: get a nullspace

      X = Teuchos::rcp_dynamic_cast<Vector>(tp_container->get_x());
      B = Teuchos::rcp_dynamic_cast<Vector>(tp_container->get_f());

      if(serialRandom) {
        //Build the seed on rank zero and broadcast it.
        size_t localNumElements = 0;
        if(comm->getRank() == 0) {
          localNumElements = static_cast<size_t>(dofMap->getGlobalNumElements());
        }
        RCP<Map> serialMap = MapFactory::Build(dofMap->lib(),
                                               dofMap->getGlobalNumElements(),
                                               localNumElements,
                                               0,
                                               comm);
        RCP<Vector> Xserial = VectorFactory::Build(serialMap);
        Xserial->setSeed(251743369);
        Xserial->randomize();
        RCP<Import> randomnessImporter = ImportFactory::Build(serialMap, dofMap);
        X->doImport(*Xserial, *randomnessImporter, Xpetra::INSERT);
      } else {
        // we set seed for reproducibility
        Utilities::SetRandomSeed(*comm);
        X->randomize();
      }

      A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

      Teuchos::Array<typename STS::magnitudeType> norms(1);
      B->norm2(norms);
      B->scale(one/norms[0]);

#ifdef MATLAB_COMPARE
      Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("Ax.mm",*B);
      Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("A.mm",*A);
      B->putScalar(zero);
      Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("rhs.mm",*B);
      Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("x.mm",*X);
#endif

      comm->barrier();
      tm = Teuchos::null;

      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - Compute region data")));

      // Set aggregation type for each region
      std::string aggregationRegionType;
      RCP<ParameterList> interfaceParams = rcp(new ParameterList());
      if(useUnstructured) {
        aggregationRegionType = "uncoupled";
      } else {
        aggregationRegionType = "structured";
      }

      const LO numLocalCompositeNodes = lNodesPerDim[0]*lNodesPerDim[1]*lNodesPerDim[2];

      // Rule for boundary duplication
      // For any two ranks that share an interface:
      // the lowest rank owns the interface and the highest rank gets extra nodes

      // 1D example of the relation between Composite, Quasi Region, and Region formats
      //
      // Composite:
      // Rank 0   Rank 1
      // [0 1 2]  [3 4]
      //
      // Quasi Region:
      // Rank 0   Rank 1
      // [0 1 2]  [2 3 4]
      //
      // Region:
      // Rank 0   Rank 1
      // [0 1 2]  [5 3 4]

      // First we count how many nodes the region needs to send and receive
      // and allocate arrays accordingly
      Array<int> boundaryConditions;
      int maxRegPerGID = 0;
      int numInterfaces = 0;
      LO numLocalRegionNodes = 0;
      Array<GO>  sendGIDs;
      Array<int> sendPIDs;
      Array<LO>  rNodesPerDim(3);
      Array<LO>  compositeToRegionLIDs(nodeMap->getNodeNumElements()*numDofsPerNode);
      Array<GO>  interfaceGIDs;
      Array<LO>  interfaceLIDsData;

      // TODO: finish generating the appropriate data that this function typically generates
//      createRegionData(num_dimensions, useUnstructured, numDofsPerNode,
//                       gNodesPerDim(), lNodesPerDim(), procsPerDim(), nodeMap, dofMap,
//                       maxRegPerGID, numLocalRegionNodes, boundaryConditions,
//                       sendGIDs, sendPIDs, numInterfaces, rNodesPerDim,
//                       quasiRegionGIDs, quasiRegionCoordGIDs, compositeToRegionLIDs,
//                       interfaceGIDs, interfaceLIDsData);

      const LO numSend = static_cast<LO>(sendGIDs.size());

      // std::cout << "p=" << myRank << " | numSend=" << numSend << std::endl;
      // << ", numReceive=" << numReceive << std::endl;
      // std::cout << "p=" << myRank << " | receiveGIDs: " << receiveGIDs << std::endl;
      // std::cout << "p=" << myRank << " | receivePIDs: " << receivePIDs << std::endl;
      // std::cout << "p=" << myRank << " | sendGIDs: " << sendGIDs << std::endl;
      // std::cout << "p=" << myRank << " | sendPIDs: " << sendPIDs << std::endl;

      // Second we actually fill the send and receive arrays with appropriate data
      // which will allow us to compute the region and composite maps.
      // Now we can construct a list of GIDs that corresponds to rowMap
      Array<LO>  interfacesDimensions, interfacesLIDs;
      if(useUnstructured) {
        findInterface(num_dimensions, rNodesPerDim, boundaryConditions,
                      interfacesDimensions, interfacesLIDs);

        // std::cout << "p=" << myRank << " | numLocalRegionNodes=" << numLocalRegionNodes
        //           << ", rNodesPerDim: " << rNodesPerDim << std::endl;
        // std::cout << "p=" << myRank << " | boundaryConditions: " << boundaryConditions << std::endl
        //           << "p=" << myRank << " | rNodesPerDim: " << rNodesPerDim << std::endl
        //           << "p=" << myRank << " | interfacesDimensions: " << interfacesDimensions << std::endl
        //           << "p=" << myRank << " | interfacesLIDs: " << interfacesLIDs << std::endl;
      }

      interfaceParams->set<Array<LO> >("interfaces: nodes per dimensions", interfacesDimensions); // nodesPerDimensions);
      interfaceParams->set<Array<LO> >("interfaces: interface nodes",      interfacesLIDs); // interfaceLIDs);

      // std::cout << "p=" << myRank << " | compositeToRegionLIDs: " << compositeToRegionLIDs << std::endl;
      // std::cout << "p=" << myRank << " | quasiRegionGIDs: " << quasiRegionGIDs << std::endl;
      // std::cout << "p=" << myRank << " | interfaceGIDs: " << interfaceGIDs << std::endl;
      // std::cout << "p=" << myRank << " | interfaceLIDsData: " << interfaceLIDsData << std::endl;
      // std::cout << "p=" << myRank << " | interfaceLIDs: " << interfaceLIDs << std::endl;
      // std::cout << "p=" << myRank << " | quasiRegionCoordGIDs: " << quasiRegionCoordGIDs() << std::endl;

      // In our very particular case we know that a node is at most shared by 4 (8) regions in 2D (3D) problems.
      // Other geometries will certainly have different constrains and a parallel reduction using MAX
      // would be appropriate.

      comm->barrier();
      tm = Teuchos::null;

      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3 - Build Region Matrix")));

      RCP<TimeMonitor> tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.1 - Build Region Maps")));

      Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap, colMap;
      Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > revisedRowMap, revisedColMap;
      rowMap = Xpetra::MapFactory<LO,GO,Node>::Build(dofMap->lib(),
                                                     Teuchos::OrdinalTraits<GO>::invalid(),
                                                     quasiRegionGIDs(),
                                                     dofMap->getIndexBase(),
                                                     dofMap->getComm());
      colMap = rowMap;
      revisedRowMap = Xpetra::MapFactory<LO,GO,Node>::Build(dofMap->lib(),
                                                            Teuchos::OrdinalTraits<GO>::invalid(),
                                                            numLocalRegionNodes*numDofsPerNode,
                                                            dofMap->getIndexBase(),
                                                            dofMap->getComm());
      revisedColMap = revisedRowMap;

      // Build objects needed to construct the region coordinates
      Teuchos::RCP<Xpetra::Map<LO,GO,NO> > quasiRegCoordMap = Xpetra::MapFactory<LO,GO,Node>::
          Build(nodeMap->lib(),
                Teuchos::OrdinalTraits<GO>::invalid(),
                quasiRegionCoordGIDs(),
                nodeMap->getIndexBase(),
                nodeMap->getComm());
      Teuchos::RCP<Xpetra::Map<LO,GO,NO> > regCoordMap = Xpetra::MapFactory<LO,GO,Node>::
          Build(nodeMap->lib(),
                Teuchos::OrdinalTraits<GO>::invalid(),
                numLocalRegionNodes,
                nodeMap->getIndexBase(),
                nodeMap->getComm());

      comm->barrier();
      tmLocal = Teuchos::null;
      tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.2 - Build Region Importers")));

      // Setup importers
      RCP<Import> rowImport;
      RCP<Import> colImport;
      rowImport = ImportFactory::Build(dofMap, rowMap);
      colImport = ImportFactory::Build(dofMap, colMap);
      RCP<Import> coordImporter = ImportFactory::Build(nodeMap, quasiRegCoordMap);

      comm->barrier();
      tmLocal = Teuchos::null;
      tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.3 - Import ghost GIDs")));

      Array<GO>  interfaceCompositeGIDs, interfaceRegionGIDs;
      ExtractListOfInterfaceRegionGIDs(revisedRowMap, interfaceLIDsData, interfaceRegionGIDs);

      RCP<Xpetra::MultiVector<LO, LO, GO, NO> > regionsPerGIDWithGhosts;
      RCP<Xpetra::MultiVector<GO, LO, GO, NO> > interfaceGIDsMV;
      MakeRegionPerGIDWithGhosts(nodeMap, revisedRowMap, rowImport,
                                 maxRegPerGID, numDofsPerNode,
                                 lNodesPerDim, sendGIDs, sendPIDs, interfaceLIDsData,
                                 regionsPerGIDWithGhosts, interfaceGIDsMV);

      Teuchos::ArrayRCP<LO> regionMatVecLIDs;
      RCP<Import> regionInterfaceImporter;
      SetupMatVec(interfaceGIDsMV, regionsPerGIDWithGhosts, revisedRowMap, rowImport,
                  regionMatVecLIDs, regionInterfaceImporter);

      comm->barrier();
      tmLocal = Teuchos::null;
      tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.4 - Build QuasiRegion Matrix")));

      std::cout << "About to create quasi region matrix" << std::endl;
      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > quasiRegionMats;
      MakeQuasiregionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A),
                              regionsPerGIDWithGhosts, rowMap, colMap, rowImport,
                              quasiRegionMats, regionMatVecLIDs);
      std::cout << "Done creating quasi region matrix" << std::endl;

      comm->barrier();
      tmLocal = Teuchos::null;
      tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.5 - Build Region Matrix")));

      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats;
      MakeRegionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A), A->getRowMap(), rowMap,
                         revisedRowMap, revisedColMap,
                         rowImport, quasiRegionMats, regionMats);

      // We don't need the composite operator on the fine level anymore. Free it!
      A = Teuchos::null;

      comm->barrier();
      tmLocal = Teuchos::null;

      comm->barrier();
      tm = Teuchos::null;

      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 4 - Build Region Hierarchy")));

      // Setting up parameters before hierarchy construction
      // These need to stay in the driver as they would be provide by an app
      Array<int> regionNodesPerDim;
      RCP<MultiVector> regionNullspace;
      RCP<RealValuedMultiVector> regionCoordinates;

      // Set mesh structure data
      regionNodesPerDim = rNodesPerDim;

      // create nullspace vector
      regionNullspace = MultiVectorFactory::Build(rowMap, nullspace->getNumVectors());
      regionNullspace->doImport(*nullspace, *rowImport, Xpetra::INSERT);
      regionNullspace->replaceMap(revisedRowMap);

      // create region coordinates vector
      regionCoordinates = Xpetra::MultiVectorFactory<real_type,LO,GO,NO>::Build(quasiRegCoordMap,
                                                                                coordinates->getNumVectors());
      regionCoordinates->doImport(*coordinates, *coordImporter, Xpetra::INSERT);
      regionCoordinates->replaceMap(regCoordMap);

      using Tpetra_CrsMatrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      using Tpetra_MultiVector = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

      // Stuff for multi-level algorithm
      //
      // To allow for multi-level schemes with more than two levels, we need to store
      // maps, matrices, vectors, and stuff like that on each level. Since we call the
      // multi-level scheme recursively, this should be reflected in the design of
      // variables.
      //
      // We use MueLu::Hierarchy and MueLu:Level to store each quantity on each level.
      //
      RCP<ParameterList> coarseSolverData = rcp(new ParameterList());
      coarseSolverData->set<std::string>("coarse solver type", coarseSolverType);
      coarseSolverData->set<bool>("coarse solver rebalance", coarseSolverRebalance);
      coarseSolverData->set<int>("coarse rebalance num partitions", rebalanceNumPartitions);
      coarseSolverData->set<std::string>("amg xml file", coarseAmgXmlFile);
      coarseSolverData->set<std::string>("smoother xml file", coarseSmootherXMLFile);
      RCP<ParameterList> hierarchyData = rcp(new ParameterList());


      // Create MueLu Hierarchy Initially...
      // Read MueLu parameter list form xml file
      RCP<ParameterList> mueluParams = Teuchos::rcp(new ParameterList());
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, mueluParams.ptr(), *dofMap->getComm());

      // Insert region-specific data into parameter list
      const std::string userName = "user data";
      Teuchos::ParameterList& userParamList = mueluParams->sublist(userName);
      userParamList.set<int>        ("int num_dimensions", num_dimensions);
      userParamList.set<Array<LO> > ("Array<LO> lNodesPerDim", regionNodesPerDim);
      userParamList.set<std::string>("string aggregationRegionType", aggregationRegionType);
      userParamList.set<Array<LO> > ("Array<LO> nodeOnInterface", interfaceParams->get<Array<LO> >("interfaces: interface nodes"));
      userParamList.set<Array<LO> > ("Array<LO> interfacesDimensions", interfaceParams->get<Array<LO> >("interfaces: nodes per dimensions"));
      if(Teuchos::nonnull(regionCoordinates)) {
        userParamList.set("Coordinates", regionCoordinates);
      }
      if(Teuchos::nonnull(regionNullspace)) {
        userParamList.set("Nullspace", regionNullspace);
      }

      tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("CreateXpetraPreconditioner: Hierarchy")));

      // Create multigrid hierarchy part 1
      RCP<Hierarchy> regHierarchy  = MueLu::CreateXpetraPreconditioner(regionMats, *mueluParams);

      {
        RCP<MueLu::Level> level = regHierarchy->GetLevel(0);
        level->Set<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("rowImport",rowImport);
        level->Set<ArrayView<LocalOrdinal> > ("compositeToRegionLIDs", compositeToRegionLIDs() );
        level->Set<RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >("interfaceGIDs", interfaceGIDsMV);
        level->Set<RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >("regionsPerGIDWithGhosts", regionsPerGIDWithGhosts);
        level->Set<Teuchos::ArrayRCP<LocalOrdinal> >("regionMatVecLIDs", regionMatVecLIDs);
        level->Set<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("regionInterfaceImporter", regionInterfaceImporter);
        level->print( std::cout, MueLu::Extreme );
      }

      tmLocal = Teuchos::null;


      // Create multigrid hierarchy part 2
      createRegionHierarchy(num_dimensions,
                            regionNodesPerDim,
                            aggregationRegionType,
                            interfaceParams,
                            maxRegPerGID,
                            coarseSolverData,
                            smootherParams,
                            hierarchyData,
                            regHierarchy,
                            keepCoarseCoords);

      hierarchyData->print();



      comm->barrier();
      tm = Teuchos::null;

      // Extract the number of levels from the prolongator data structure
      const int numLevels = regHierarchy->GetNumLevels();

      // Set data for fast MatVec
      for(LO levelIdx = 0; levelIdx < numLevels; ++levelIdx) {
        RCP<MueLu::Level> level = regHierarchy->GetLevel(levelIdx);
        RCP<Xpetra::Import<LO, GO, NO> > regionInterfaceImport = level->Get<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("regionInterfaceImporter");
        Teuchos::ArrayRCP<LO>            regionMatVecLIDs1     = level->Get<Teuchos::ArrayRCP<LO> >("regionMatVecLIDs");
        smootherParams[levelIdx]->set("Fast MatVec: interface LIDs",
                                      regionMatVecLIDs1);
        smootherParams[levelIdx]->set("Fast MatVec: interface importer",
                                      regionInterfaceImport);
      }

      // RCP<Teuchos::FancyOStream> fancy2 = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      // Teuchos::FancyOStream& out2 = *fancy2;
      // for(LO levelIdx = 0; levelIdx < numLevels; ++levelIdx) {
      //   out2 << "p=" << myRank << " | regionMatVecLIDs on level " << levelIdx << std::endl;
      //   regionMatVecLIDsPerLevel[levelIdx]->describe(out2, Teuchos::VERB_EXTREME);
      // }

      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 5 - Solve with V-cycle")));

      {
        //    std::cout << myRank << " | Running V-cycle ..." << std::endl;

        TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

        // We first use the non-level container variables to setup the fine grid problem.
        // This is ok since the initial setup just mimics the application and the outer
        // Krylov method.
        //
        // We switch to using the level container variables as soon as we enter the
        // recursive part of the algorithm.
        //

        // Composite residual vector
        RCP<Vector> compRes = VectorFactory::Build(dofMap, true);

        // transform composite vectors to regional layout
        Teuchos::RCP<Vector> quasiRegX;
        Teuchos::RCP<Vector> regX;
        compositeToRegional(X, quasiRegX, regX,
                            revisedRowMap, rowImport);

        RCP<Vector> quasiRegB;
        RCP<Vector> regB;
        compositeToRegional(B, quasiRegB, regB,
                            revisedRowMap, rowImport);
#ifdef DUMP_LOCALX_AND_A
        FILE *fp;
        char str[80];
        sprintf(str,"theMatrix.%d",myRank);
        fp = fopen(str,"w");
        fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
        LO numNzs = 0;
        for (size_t kkk = 0; kkk < regionMats->getNodeNumRows(); kkk++) {
          ArrayView<const LO> AAcols;
          ArrayView<const SC> AAvals;
          regionMats->getLocalRowView(kkk, AAcols, AAvals);
          const int *Acols    = AAcols.getRawPtr();
          const SC  *Avals = AAvals.getRawPtr();
          numNzs += AAvals.size();
        }
        fprintf(fp, "%d %d %d\n",regionMats->getNodeNumRows(),regionMats->getNodeNumRows(),numNzs);

        for (size_t kkk = 0; kkk < regionMats->getNodeNumRows(); kkk++) {
          ArrayView<const LO> AAcols;
          ArrayView<const SC> AAvals;
          regionMats->getLocalRowView(kkk, AAcols, AAvals);
          const int *Acols    = AAcols.getRawPtr();
          const SC  *Avals = AAvals.getRawPtr();
          LO RowLeng = AAvals.size();
          for (LO kk = 0; kk < RowLeng; kk++) {
            fprintf(fp, "%d %d %22.16e\n",kkk+1,Acols[kk]+1,Avals[kk]);
          }
        }
        fclose(fp);
        sprintf(str,"theX.%d",myRank);
        fp = fopen(str,"w");
        ArrayRCP<SC> lX= regX->getDataNonConst(0);
        for (size_t kkk = 0; kkk < regionMats->getNodeNumRows(); kkk++) fprintf(fp, "%22.16e\n",lX[kkk]);
        fclose(fp);
#endif

        RCP<Vector> regRes;
        regRes = VectorFactory::Build(revisedRowMap, true);

        /////////////////////////////////////////////////////////////////////////
        // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
        /////////////////////////////////////////////////////////////////////////

        // Prepare output of residual norm to file
        RCP<std::ofstream> log;
        if (myRank == 0)
        {
          log = rcp(new std::ofstream(convergenceLog.c_str()));
          (*log) << "# num procs = " << dofMap->getComm()->getSize() << "\n"
              << "# iteration | res-norm (scaled=" << scaleResidualHist << ")\n"
              << "#\n";
          *log << std::setprecision(16) << std::scientific;
        }

        // Print type of residual norm to the screen
        if (scaleResidualHist)
          out << "Using scaled residual norm." << std::endl;
        else
          out << "Using unscaled residual norm." << std::endl;


        // Richardson iterations
        magnitude_type normResIni = Teuchos::ScalarTraits<magnitude_type>::zero();
        const int old_precision = std::cout.precision();
        std::cout << std::setprecision(8) << std::scientific;
        int cycle = 0;

        Teuchos::RCP<Vector> regCorrect;
        regCorrect = VectorFactory::Build(revisedRowMap, true);
        for (cycle = 0; cycle < maxIts; ++cycle)
        {
          const Scalar SC_ZERO = Teuchos::ScalarTraits<SC>::zero();
          regCorrect->putScalar(SC_ZERO);
          // Get Stuff out of Hierarchy
          RCP<MueLu::Level> level = regHierarchy->GetLevel(0);
          RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal> > regInterfaceScalings = level->Get<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal> > >("regInterfaceScalings");
          // check for convergence
          {
            ////////////////////////////////////////////////////////////////////////
            // SWITCH BACK TO NON-LEVEL VARIABLES
            ////////////////////////////////////////////////////////////////////////
            computeResidual(regRes, regX, regB, regionMats, *smootherParams[0]);
            scaleInterfaceDOFs(regRes, regInterfaceScalings, true);

            compRes = VectorFactory::Build(dofMap, true);
            regionalToComposite(regRes, compRes, rowImport);

            typename Teuchos::ScalarTraits<Scalar>::magnitudeType normRes = compRes->norm2();
            if(cycle == 0) { normResIni = normRes; }

            if (scaleResidualHist)
              normRes /= normResIni;

            // Output current residual norm to screen (on proc 0 only)
            out << cycle << "\t" << normRes << std::endl;
            if (myRank == 0)
              (*log) << cycle << "\t" << normRes << "\n";

            if (normRes < tol)
              break;
          }

          /////////////////////////////////////////////////////////////////////////
          // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
          /////////////////////////////////////////////////////////////////////////

          bool zeroInitGuess = true;
          scaleInterfaceDOFs(regRes, regInterfaceScalings, false);
          vCycle(0, numLevels, cycleType, regHierarchy,
                 regCorrect, regRes,
                 smootherParams, zeroInitGuess, coarseSolverData, hierarchyData);

          regX->update(one, *regCorrect, one);
        }
        out << "Number of iterations performed for this solve: " << cycle << std::endl;

        std::cout << std::setprecision(old_precision);
        std::cout.unsetf(std::ios::fixed | std::ios::scientific);
      }

      comm->barrier();
      tm = Teuchos::null;
      globalTimeMonitor = Teuchos::null;

      if (showTimerSummary)
      {
        RCP<ParameterList> reportParams = rcp(new ParameterList);
        const std::string filter = "";
        if (useStackedTimer) {
          Teuchos::StackedTimer::OutputOptions options;
          options.output_fraction = options.output_histogram = options.output_minmax = true;
          stacked_timer->report(out, comm, options);
        } else {
          std::ios_base::fmtflags ff(out.flags());
          TimeMonitor::report(comm.ptr(), out, filter, reportParams);
          out << std::setiosflags(ff);
        }
      }

      TimeMonitor::clearCounters();

      return EXIT_SUCCESS;

    }
    
    ////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// OUTPUT RESULTS ///////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////

    tm = Teuchos::null;
    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 6 - Output Data")));

    // TODO: MueLu ordering and Panzer ordering will likely not match here... we'll need to run a conversion
    // write the solution to matrix
    {
      // redistribute solution vector to ghosted vector
      linObjFactory->globalToGhostContainer(*container,*ghostCont, panzer::TpetraLinearObjContainer<ST,LO,GO>::X
                                            | panzer::TpetraLinearObjContainer<ST,LO,GO>::DxDt);

      // get X Tpetra_Vector from ghosted container
      // TODO: there is some magic here with Tpetra objects that needs to be fixed
      //Teuchos::RCP<panzer::TpetraLinearObjContainer<ST,LO,GO> > tp_ghostCont = Teuchos::rcp_dynamic_cast<panzer::TpetraLinearObjContainer<ST,LO,GO> >(ghostCont);
      //panzer_stk::write_solution_data(*dofManager,*mesh,*tp_ghostCont->get_x());

      std::ostringstream filename;
      filename << "regionMG_output" << discretization_order << ".exo";
      mesh->writeToExodus(filename.str());
    }

    // compute the error of the finite element solution
    /////////////////////////////////////////////////////////////

    {
      panzer::AssemblyEngineInArgs respInput(ghostCont,container);
      respInput.alpha = 0;
      respInput.beta = 1;

      Teuchos::RCP<panzer::ResponseBase> l2_resp = errorResponseLibrary->getResponse<panzer::Traits::Residual>("L2 Error");
      Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > l2_resp_func = Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(l2_resp);
      Teuchos::RCP<Thyra::VectorBase<double> > l2_respVec = Thyra::createMember(l2_resp_func->getVectorSpace());
      l2_resp_func->setVector(l2_respVec);


//      Teuchos::RCP<panzer::ResponseBase> h1_resp = errorResponseLibrary->getResponse<panzer::Traits::Residual>("H1 Error");
//      Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > h1_resp_func = Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(h1_resp);
//      Teuchos::RCP<Thyra::VectorBase<double> > h1_respVec = Thyra::createMember(h1_resp_func->getVectorSpace());
//      h1_resp_func->setVector(h1_respVec);


      errorResponseLibrary->addResponsesToInArgs<panzer::Traits::Residual>(respInput);
      errorResponseLibrary->evaluate<panzer::Traits::Residual>(respInput);

      out << "This is the Basis Order" << std::endl;
      out << "Basis Order = " << discretization_order << std::endl;
      out << "This is the L2 Error" << std::endl;
      out << "L2 Error = " << sqrt(l2_resp_func->value) << std::endl;
      //out << "This is the H1 Error" << std::endl;
      //out << "H1 Error = " << sqrt(h1_resp_func->value) << std::endl;
    }

    tm = Teuchos::null;
    globalTimeMonitor = Teuchos::null;

    if (showTimerSummary)
    {
      RCP<ParameterList> reportParams = rcp(new ParameterList);
      const std::string filter = "";
      if (useStackedTimer)
      {
        Teuchos::StackedTimer::OutputOptions options;
        options.output_fraction = options.output_histogram = options.output_minmax = true;
        stacked_timer->report(out, comm, options);
      }
      else
      {
        std::ios_base::fmtflags ff(out.flags());
        Teuchos::TimeMonitor::report(comm.ptr(), out, filter, reportParams);
        out << std::setiosflags(ff);
      }
    }

  } // Parallel initialization scope
  Tpetra::finalize();
  Kokkos::finalize();

  // TODO: Does this need to be scoped again for MPI? Double-check based on review

  return 0;
} // main
