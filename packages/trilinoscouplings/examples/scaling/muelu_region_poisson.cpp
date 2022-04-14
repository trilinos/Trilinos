// Standard headers
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unistd.h>

// Kokkos headers have to go first
#include "Kokkos_View.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"

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
#include "SolveRegionHierarchy_def.hpp"

// Shards headers
#include "Shards_CellTopology.hpp"

// Tpetra headers
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_FECrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Include factories for boundary conditions and other Panzer setup
// Most of which is taken from PoissonExample in Panzer_STK
#include "muelu_region_poisson.hpp"

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


// need this for debugging purposes... technically any class that has a size()
// method and [] operator that returns a type that can be streamed works,
// but the stream needs to be able to use std::endl;
#define DUMPSTDVECTOR(vector, stream) \
  stream << #vector << ".size() = " << vector.size() << std::endl; \
  stream << #vector << "        = (" << vector[0]; \
  for(unsigned int i=1; i < vector.size(); ++i) \
    stream << ", " << vector[i]; \
  stream << ")" << std::endl;

// print Teuchos::ArrayRCPs
template <typename T>
inline void printArrayRCP(Teuchos::ArrayRCP<T> arrayrcp, std::string name="ArrayRCP") {
  std::cout << name << " = [";
  for(unsigned int i=0; i<arrayrcp.size(); ++i)
    std::cout << arrayrcp[i] << " ";
  std::cout << "]" << std::endl;
}

// print std::vectors
template <typename T>
inline void printStdVector(std::vector<T> vector, std::string name="vector") {
  std::cout << name << " = [";
  for(unsigned int i=0; i<vector.size(); ++i)
    std::cout << vector[i] << " ";
  std::cout << "]" << std::endl;
}

// print Kokkos::Views (defined as a macro to avoid template parameter hell) ((use responsibly))
#define PRINT_VIEW1(view)                                \
std::cout << #view << " (" << view.extent(0) << ") = ["; \
for(unsigned int i=0; i<view.extent(0)-1; ++i)           \
  std::cout << view(i) << ", ";                          \
std::cout << view(view.extent(0)-1) << "]" << std::endl;

#define PRINT_VIEW2(view)                                                         \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << ") = " << std::endl;  \
for(unsigned int i=0; i<view.extent(0); ++i) {                                    \
  for(unsigned int j=0; j<view.extent(1); ++j)                                    \
    std::cout << view(i,j) << " ";                                                \
  std::cout << std::endl << " ";                                                  \
}

// some things like multivectors are "2D views" but only appear as 1D in practice, so the above print isn't as pretty
#define PRINT_VIEW2_LINEAR(view)                                                               \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << ") = [" << std::endl;  \
for(unsigned int i=0; i<view.extent(0); ++i)                                      \
  for(unsigned int j=0; j<view.extent(1); ++j)                                    \
    std::cout << view(i,j) << " ";                                                \
std::cout << "]" << std::endl;

#define PRINT_VIEW3(view)                                                                                  \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << "," << view.extent(2) << ") = " << std::endl;  \
for(unsigned int i=0; i<view.extent(0); ++i) {                                                             \
  if(i==0) std::cout << "    [";                                                                           \
  for(unsigned int j=0; j<view.extent(1); ++j) {                                                           \
    for(unsigned int k=0; k<view.extent(2); ++k)                                                           \
      std::cout << view(i,j,k) << " ";                                                                     \
  std::cout << "]" << std::endl << "    [";                                                                \
  }                                                                                                        \
}                                                                                                          \
std::cout << "]" << std::endl;

#define PRINT_DRV(view)                                                                                    \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << "," << view.extent(2) << ") = ["; \
for(unsigned int i=0; i<view.extent(0); ++i)                                                               \
  for(unsigned int j=0; j<view.extent(1); ++j)                                                             \
    for(unsigned int k=0; k<view.extent(2); ++k)                                                           \
      std::cout << view(i,j,k) << " ";                                                                     \
std::cout << "]" << std::endl;

#define PRINT_VAR(var)                         \
std::cout << #var << "=" << var << std::endl;

#define PRINT_VIEW2_MAX(view)                  \
double max = -100000;                          \
for(unsigned int i=0; i<view.extent(0); ++i)   \
  for(unsigned int j=0; j<view.extent(1); ++j) \
    if(view(i,j) > max)                        \
      max = view(i,j);                         \
std::cout << #view << " max=" << max << std::endl;




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


    // ---------------------------------------------------------------------------------
    // ------------------------------------- SETUP -------------------------------------
    // ---------------------------------------------------------------------------------


    // Setup output stream, MPI, and grab info
    Teuchos::FancyOStream out_root(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream out_all(Teuchos::rcpFromRef(std::cout));
    out_root.setOutputToRootOnly(0);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    const int num_ranks = comm->getSize();
    const int my_rank = comm->getRank();

    out_root << "Running TrilinosCouplings region multigrid driver on " << num_ranks << " ranks... \n";

    // Parse command line arguments
    Teuchos::CommandLineProcessor clp(false);
    std::string exodus_name           = "";                  clp.setOption("exodus-mesh",           &exodus_name,             "Exodus hex mesh filename (overrides a pamgen-mesh if both specified)");
    std::string pamgen_name  = "muelu_region_poisson_1.gen"; clp.setOption("pamgen-mesh",           &pamgen_name,             "Pamgen hex mesh filename");
    std::string xmlFileName           = "";                  clp.setOption("xml",                   &xmlFileName,             "MueLu parameters from an xml file");
    std::string yamlFileName          = "";                  clp.setOption("yaml",                  &yamlFileName,            "MueLu parameters from a yaml file");
    int mesh_refinements              = 1;                   clp.setOption("mesh-refinements",      &mesh_refinements,        "Uniform mesh refinements");
    int discretization_order          = 1;                   clp.setOption("discretization-order",  &discretization_order,    "Finite element discretization order");
    int workset_size                  = 10;                  clp.setOption("workset-size",          &workset_size,            "Size of underlying Panzer worksets");
    bool delete_parent_elements       = false;               clp.setOption("delete-parent-elements", "keep-parent-elements", &delete_parent_elements,"Save the parent elements in the perceptMesh");


    // Multigrid options (revisit these later)
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

    TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName != "" && yamlFileName != "", std::runtime_error, "Cannot provide both xml and yaml input files");

    // get xml file from command line if provided, otherwise use default
    std::string  xml_solver_file_name(xmlFileName);

    // Read xml file into parameter list
    Teuchos::ParameterList inputSolverList;

    if(xml_solver_file_name.length()) {
      out_root << "\nReading parameter list from the XML file \"" << xml_solver_file_name << "\" ...\n" << std::endl;
      Teuchos::updateParametersFromXmlFile(xml_solver_file_name, Teuchos::ptr(&inputSolverList));
    }
    else
      out_root << "Using default solver values ..." << std::endl;



    // ---------------------------------------------------------------------------------
    // ----------------------------- MESH AND WORKSETS ---------------------------------
    // ---------------------------------------------------------------------------------
    
    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
    if(useStackedTimer)
      stacked_timer = rcp(new Teuchos::StackedTimer("MueLu_Region_Poisson_Driver"));
    Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
    RCP<Teuchos::TimeMonitor> globalTimeMonitor = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: S - Global Time")));
    RCP<Teuchos::TimeMonitor> tm                = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 1 - Build Mesh and Assign Physics")));

    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
    Teuchos::RCP<Teuchos::ParameterList> mesh_pl = Teuchos::rcp(new Teuchos::ParameterList);

    // TODO: due to #8475, we have this. we will clean up Panzer and this code once we verify
    // the functionality of the full HHG scheme obtained by using Percept
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory_percept;
    Teuchos::RCP<Teuchos::ParameterList> mesh_pl_percept = Teuchos::rcp(new Teuchos::ParameterList);

    bool use_exodus_mesh = exodus_name.length() > 0;
    bool use_pamgen_mesh = !use_exodus_mesh;
    bool use_percept = mesh_refinements > 0;

    if(use_exodus_mesh)
    {
      // set the filename and type
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory);
      mesh_pl->set("File Name", exodus_name);
      mesh_pl->set("File Type", "Exodus");
    }
    else if(use_pamgen_mesh)
    {
      // set the filename and type
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory);
      mesh_pl->set("File Name", pamgen_name);
      mesh_pl->set("File Type", "Pamgen");

      if(use_percept)
      {
        mesh_pl->set("Levels of Uniform Refinement", mesh_refinements); // this multiplies the number of elements by 2^(dimension*level)

        mesh_factory_percept = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory);
        mesh_pl_percept->set("File Name", pamgen_name);
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

    // grab the number and names of mesh blocks
    std::vector<std::string> element_blocks;
    mesh->getElementBlockNames(element_blocks);
    std::vector<bool> unstructured_eBlocks(element_blocks.size(), false);
    // TODO: set unstructured blocks based on some sort of input information; for example, using the Exodus ex_get_var* functions
    // TODO: add options for setting/choosing boundaries on top and bottom of cylinder

    // grab the number and name of nodesets
    std::vector<std::string> nodesets;
    mesh->getNodesetNames(nodesets);

    // grab the number and names of sidesets
    std::vector<std::string> sidesets;
    mesh->getSidesetNames(sidesets);

    // create a physics blocks parameter list
    Teuchos::RCP<Teuchos::ParameterList> physics_block_settings = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physics_blocks;

    // set physics and boundary conditions on each block
    {
      bool build_transient_support = false;

      const int integration_order = 10;
      Teuchos::ParameterList& p = physics_block_settings->sublist("Poisson Physics");
      p.set("Type","Poisson");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",discretization_order);
      p.set("Integration Order",integration_order);

      // TODO: double-check. this assumes we impose Dirichlet BCs on all boundaries of all element blocks
      // It may potentially assign Dirichlet BCs to internal block boundaries, which is undesirable
      for(size_t i=0; i<element_blocks.size(); ++i)
      {
        for(size_t j=0; j<sidesets.size(); ++j)
        {
          std::size_t bc_id = j;
          panzer::BCType bctype = panzer::BCT_Dirichlet;
          std::string sideset_id = sidesets[j];
          std::string element_block_id = element_blocks[i];
          std::string dof_name = "TEMPERATURE";
          std::string strategy = "Constant";
          double value = 0.0;
	        Teuchos::ParameterList pbc; // this is how the official example does it, so I'll leave it alone for now
          pbc.set("Value",value);
          panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name,
                        strategy, pbc);
          bcs.push_back(bc);
        }
        const panzer::CellData volume_cell_data(workset_size, mesh->getCellTopology(element_blocks[i]));

        // GobalData sets ostream and parameter interface to physics
        Teuchos::RCP<panzer::GlobalData> global_data = panzer::createGlobalData();

        // Can be overridden by the equation set
        int default_integration_order = 1;

        // the physics block nows how to build and register evaluator with the field manager
        Teuchos::RCP<panzer::PhysicsBlock> physics_block
        = Teuchos::rcp(new panzer::PhysicsBlock(physics_block_settings,
                                                element_blocks[i],
                                                default_integration_order,
                                                volume_cell_data,
                                                eqset_factory,
                                                global_data,
                                                build_transient_support));

        // we can have more than one physics block, one per element block
        physics_blocks.push_back(physics_block);
      }
    }
    panzer::checkBCConsistency(element_blocks,sidesets,bcs);


    // finish building mesh, set required field variables and mesh bulk data
    for(size_t i=0; i<physics_blocks.size(); ++i)
    {
      Teuchos::RCP<panzer::PhysicsBlock> physics_block = physics_blocks[i]; // we are assuming only one physics block

      const std::vector<panzer::StrPureBasisPair> & blockFields = physics_block->getProvidedDOFs();

      // insert all fields into a set
      std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp> fieldNames;
      fieldNames.insert(blockFields.begin(),blockFields.end());

      // add basis to DOF manager: block specific
      std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp>::const_iterator fieldItr;
      for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr)
        mesh->addSolutionField(fieldItr->first,physics_block->elementBlockID());
    }
    mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD); // this is where the mesh refinements are applied

    unsigned int num_dimensions = mesh->getDimension();
    TEUCHOS_ASSERT(num_dimensions == 2 || num_dimensions == 3); // no other mesh would make sense

    // ---------------------------------------------------------------------------------
    // -------------------------------- BUILD DOF MANAGER ------------------------------
    // ---------------------------------------------------------------------------------
    tm = Teuchos::null;
    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 2 - Build DOF Manager and Worksets")));

    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory global_indexer_factory;
    Teuchos::RCP<panzer::GlobalIndexer> dof_manager = global_indexer_factory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physics_blocks,conn_manager);

    // construct some linear algebra object, build object to pass to evaluators
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linear_object_factory = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,ST,LO,GO>(comm.getConst(),dof_manager));

    // build STK workset factory and attach it to a workset container (uses lazy evaluation)
    Teuchos::RCP<panzer_stk::WorksetFactory> workset_factory = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh));
    Teuchos::RCP<panzer::WorksetContainer> workset_container = Teuchos::rcp(new panzer::WorksetContainer);
    workset_container->setFactory(workset_factory);
    for(size_t i=0;i<physics_blocks.size();i++)
      workset_container->setNeeds(physics_blocks[i]->elementBlockID(),physics_blocks[i]->getWorksetNeeds());
    workset_container->setWorksetSize(workset_size);
    workset_container->setGlobalIndexer(dof_manager);



    // ---------------------------------------------------------------------------------
    // -------------------------------- CONSTRUCT REGIONS ------------------------------
    // ---------------------------------------------------------------------------------
    tm = Teuchos::null;
    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 3 - Setup Region Information")));

    // The code in this section assumes that a region hierarchy can be established with Percept. 
    // In the case where a region hierarchy is constructed from an exodus data input, for example, 
    // this implementation will need to be updated.
    // TODO: Assign MPI rank p to region p and collect element IDs. If this region is assigned via
    // percept, reorder the element IDs lexicographically using the utility in the header. 
    // Then collect mesh->identifier(node) for the nodes in lexicographic order for region p, 
    // and put those in quasiRegionGIDs. Coordinates should be able to be extracted from the 
    // stk::mesh::entity node as well.

    const unsigned int region_cells_per_dim = (1 << mesh_refinements);
    const unsigned int region_cells_per_proc = std::pow(region_cells_per_dim, num_dimensions);
    const unsigned int region_nodes_per_dim = region_cells_per_dim + 1;
    const unsigned int region_nodes_per_proc = std::pow(region_nodes_per_dim, num_dimensions);

    // Do a little bit of output here, regardless of verbosity
    out_root << "Driver parameters: " << std::endl;
    out_root << "          Dimension = " << num_dimensions << std::endl;
    out_root << "  Refinement Levels = " << mesh_refinements << std::endl;
    out_root << "   Region cells/dim = " << region_cells_per_dim << std::endl;
    out_root << "   Region nodes/dim = " << region_nodes_per_dim << std::endl;
    out_root << std::endl;

    // initialize data here that we will use for the region MG solver
    std::vector<GO> local_child_element_gids; // these don't always start at 0, and changes I'm making to Panzer keep changing this, so I'll store them for now
    std::vector<GO> local_child_element_region_gids;
    std::vector<GO> local_element_gids;
    Array<ST> x_coords(region_nodes_per_proc, 0);
    Array<ST> y_coords(region_nodes_per_proc, 0);
    Array<ST> z_coords(region_nodes_per_proc, 0);
    Array<GO>  quasiRegionGIDs(region_nodes_per_proc, 0);
    Array<GO>  quasiRegionCoordGIDs;

    // I'm bothered by how much digging into stk meshes is required in this part because it's not feasible for solvers,
    // so I'm going to assume 1 region per processor and just grab the element reorder using ONLY the DOFManager

    std::vector<size_t> region_cells_lexicographic_order = renumberPerceptCellsToLexicographic(mesh_refinements, num_dimensions);
    Kokkos::View<const LO**, Kokkos::Serial> LIDs = dof_manager->getLIDs();

    TEUCHOS_ASSERT(mesh->getEntityCounts(stk::topology::ELEM_RANK) == num_ranks*region_cells_lexicographic_order.size()); // this checks 1 region per processor
    TEUCHOS_ASSERT(LIDs.extent(0) == region_cells_lexicographic_order.size()); // this checks we have the mesh and reordering have the same number of cells

    Kokkos::DynRankView<ST, Kokkos::Serial> vertices;
    mesh->getElementVertices(region_cells_lexicographic_order,vertices);
    PRINT_VIEW3(vertices)

    // for loops, I need these offsets sometimes to make it more readable
    const unsigned int plane_offset = region_nodes_per_dim*region_nodes_per_dim;
    const unsigned int line_offset = region_nodes_per_dim;

    // TODO: this is hard-coded for 3D...
    for (unsigned int k = 0; k < region_cells_per_dim; k++) { // z direction
      for (unsigned int j = 0; j < region_cells_per_dim; j++) { // y direction
        for (unsigned int i = 0; i < region_cells_per_dim; i++) { // x direction

          const unsigned int i_node = region_nodes_per_dim*region_nodes_per_dim*k + region_nodes_per_dim*j + i;
          const unsigned int i_elem = region_cells_per_dim*region_cells_per_dim*k + region_cells_per_dim*j + i;

          // the bottom-left corner of a cell is always inside the region when following lexicographic ordering
          quasiRegionGIDs[i_node] = LIDs(region_cells_lexicographic_order[i_elem],0);
          x_coords[i_node] = vertices(i_elem,0,0);
          y_coords[i_node] = vertices(i_elem,0,1);
          z_coords[i_node] = vertices(i_elem,0,2);

          // if we're on the +z side of a region
          if(k==region_cells_per_dim-1)
          {
            quasiRegionGIDs[i_node+plane_offset] = LIDs(region_cells_lexicographic_order[i_elem],4);
            x_coords[i_node+plane_offset] = vertices(i_elem,4,0);
            y_coords[i_node+plane_offset] = vertices(i_elem,4,1);
            z_coords[i_node+plane_offset] = vertices(i_elem,4,2);
            // if we're on a +yz edge
            if(j==region_cells_per_dim-1)
            {
              quasiRegionGIDs[i_node+plane_offset+line_offset] = LIDs(region_cells_lexicographic_order[i_elem],7);
              x_coords[i_node+plane_offset+line_offset] = vertices(i_elem,7,0);
              y_coords[i_node+plane_offset+line_offset] = vertices(i_elem,7,1);
              z_coords[i_node+plane_offset+line_offset] = vertices(i_elem,7,2);
              // if we're on an +xyz corner
              if(i==region_cells_per_dim-1)
              {
                quasiRegionGIDs[i_node+plane_offset+line_offset+1] = LIDs(region_cells_lexicographic_order[i_elem],6);
                x_coords[i_node+plane_offset+line_offset+1] = vertices(i_elem,6,0);
                y_coords[i_node+plane_offset+line_offset+1] = vertices(i_elem,6,1);
                z_coords[i_node+plane_offset+line_offset+1] = vertices(i_elem,6,2);
              }
            }
            // if we're on the +xz edge
            if(i==region_cells_per_dim-1)
            {
              quasiRegionGIDs[i_node+plane_offset+1] = LIDs(region_cells_lexicographic_order[i_elem],5);
              x_coords[i_node+plane_offset+1] = vertices(i_elem,5,0);
              y_coords[i_node+plane_offset+1] = vertices(i_elem,5,1);
              z_coords[i_node+plane_offset+1] = vertices(i_elem,5,2);
            }
          }
          // if we're on the +y side of a region
          if(j==region_cells_per_dim-1)
          {
            quasiRegionGIDs[i_node+line_offset] = LIDs(region_cells_lexicographic_order[i_elem],3);
            x_coords[i_node+line_offset] = vertices(i_elem,3,0);
            y_coords[i_node+line_offset] = vertices(i_elem,3,1);
            z_coords[i_node+line_offset] = vertices(i_elem,3,2);
            // if we're on a +xy edge
            if(i==region_cells_per_dim-1)
            {
              quasiRegionGIDs[i_node+line_offset+1] = LIDs(region_cells_lexicographic_order[i_elem],2);
              x_coords[i_node+line_offset+1] = vertices(i_elem,2,0);
              y_coords[i_node+line_offset+1] = vertices(i_elem,2,1);
              z_coords[i_node+line_offset+1] = vertices(i_elem,2,2);
            }
          }
          // if we're on the +x side of a region
          if(i==region_cells_per_dim-1)
          {
            quasiRegionGIDs[i_node+1] = LIDs(region_cells_lexicographic_order[i_elem],1);
            x_coords[i_node+1] = vertices(i_elem,1,0);
            y_coords[i_node+1] = vertices(i_elem,1,1);
            z_coords[i_node+1] = vertices(i_elem,1,2);
          }
        }
      }
    }

    DUMPSTDVECTOR(x_coords, std::cout)
    DUMPSTDVECTOR(y_coords, std::cout)
    DUMPSTDVECTOR(z_coords, std::cout)
    PRINT_VIEW2(LIDs)
    DUMPSTDVECTOR(quasiRegionGIDs, std::cout)
    

//     // if we use Percept, this is the strategy to follow.
//     // this may be turned into its own function later
//     if(use_percept)
//     {
//       // get the Percept mesh from Panzer
//       Teuchos::RCP<percept::PerceptMesh> percept_mesh = parent_data_mesh->getRefinedMesh();
//       if(print_percept_mesh)
//         percept_mesh->print_info(out_root,"",1,true);

//       // ids are linear within stk, but we need an offset because the original mesh info comes first
// //      size_t node_id_start = 0;
// //      {
// //        const stk::mesh::BucketVector & local_buckets = percept_mesh->get_bulk_data()->get_buckets(stk::topology::ELEM_RANK,percept_mesh->get_fem_meta_data()->locally_owned_part());
// //        //const stk::mesh::BucketVector & buckets = percept_mesh->get_bulk_data()->buckets(percept_mesh->node_rank());
// //        stk::mesh::Bucket & bucket = **local_buckets.begin() ;
// //        node_id_start = percept_mesh->id(bucket[0]);
// //        if(print_debug_info)
// //          debug << "Starting node id = " << node_id_start << std::endl;
// //      }
// //
// //      size_t elem_id_start = 0;
// //      {
// //        const stk::mesh::BucketVector & local_buckets = percept_mesh->get_bulk_data()->get_buckets(stk::topology::ELEM_RANK,percept_mesh->get_fem_meta_data()->locally_owned_part());
// //        //const stk::mesh::BucketVector & buckets = percept_mesh->get_bulk_data()->buckets(percept_mesh->element_rank());
// //        stk::mesh::Bucket & bucket = **local_buckets.begin() ;
// //        elem_id_start = percept_mesh->id(bucket[0]);
// //        if(print_debug_info)
// //          debug << "Starting element id = " << elem_id_start << std::endl;
// //      }
//       //panzer_stk::workset_utils::getIdsAndVertices


//       // grab the region information from the mesh that keeps parent elements
//       {
//         // count parents and children
//         int npar=0;
//         int nchild=0;

//         const stk::mesh::BucketVector & local_buckets = percept_mesh->get_bulk_data()->get_buckets(stk::topology::ELEM_RANK,percept_mesh->get_fem_meta_data()->locally_owned_part());
//         for (stk::mesh::BucketVector::const_iterator k = local_buckets.begin(); k != local_buckets.end(); ++k)
//         {
//           const stk::mesh::Bucket & bucket = **k ;
//           if(print_debug_info)
//             debug << "New bucket" << std::endl;

//           const unsigned int num_elements_in_bucket = bucket.size();
//           for (unsigned int iElement = 0; iElement < num_elements_in_bucket; iElement++)
//           {
//             const stk::mesh::Entity element = bucket[iElement];
//             if (!percept_mesh->isParentElement(element, false))
//             {
//               ++nchild;

//               // this is the important part here. take the id of the element and the id of the element's root
//               local_child_element_gids.push_back(percept_mesh->id(element));
//               local_child_element_region_gids.push_back(percept_mesh->id(percept_mesh->rootOfTree(element)));

//               if(print_debug_info)
//                 debug << "Percept Element = " << percept_mesh->id(element) << std::endl;

//               const percept::MyPairIterRelation elem_nodes ( *percept_mesh, element,  stk::topology::NODE_RANK);

//               for (unsigned int i_node = 0; i_node < elem_nodes.size(); i_node++)
//               {
//                 const stk::mesh::Entity node = elem_nodes[i_node].entity();
//                 //local_node_gids.push_back(percept_mesh->id(node));

//                 if(print_debug_info)
//                   debug << "Stk Node = " << percept_mesh->id(node) << std::endl;
//               }
//             }
//             else
//             {
//               if(print_debug_info)
//                 debug << "p = " << my_rank << ", parent = " << percept_mesh->id(element) << std::endl;
//             }
//           }
//         }

//         if(print_debug_info)
//         {
//           for(unsigned int i=0; i<local_child_element_gids.size(); ++i)
//           {
//             out_root << "child= " << local_child_element_gids[i] << " parent= " << local_child_element_region_gids[i] << std::endl;
//           }
//         }
//       }

//       // however, the mesh that keeps parent elements does not assemble finite
//       // element data properly; therefore, we must now convert from element IDs
//       // on the mesh with parent elements to element IDs on the mesh without
//       // parent elements.
//       {
//         const stk::mesh::BucketVector & local_buckets = mesh->getBulkData()->get_buckets(stk::topology::ELEM_RANK,mesh->getMetaData()->locally_owned_part());
//         for (stk::mesh::BucketVector::const_iterator k = local_buckets.begin(); k != local_buckets.end(); ++k)
//         {
//           const stk::mesh::Bucket & bucket = **k;

//           const unsigned int num_elements_in_bucket = bucket.size();
//           for (unsigned int iElement = 0; iElement < num_elements_in_bucket; iElement++)
//           {
//             const stk::mesh::Entity element = bucket[iElement];
//             local_element_gids.push_back(mesh->getBulkData()->identifier(element));
//           }
//         }
//       }

//       // grab the Percept element renumbering
//       sleep(my_rank);
//       const std::vector<unsigned int> percept_lexicographic_elements = renumberPerceptCellsToLexicographic(mesh_refinements);
//       DUMPSTDVECTOR(percept_lexicographic_elements, debug);
//       DUMPSTDVECTOR(local_child_element_gids, debug);
//       DUMPSTDVECTOR(local_element_gids, debug);
//       comm->barrier();

//       // I think this is no longer necessary. Since elements are locally ordered,
//       // percept_lexicographic_elements is exactly the local indexing we need.
//       // // apply the renumbering
//       //std::vector<unsigned int> renumbered_local_element_gids;
//       //for(unsigned int i=0; i<percept_lexicographic_elements.size(); ++i)
//       //  renumbered_local_element_gids.push_back(local_element_gids[percept_lexicographic_elements[i]]);
//       // // make sure we have the correct number of elements on this rank
//       //TEUCHOS_ASSERT(renumbered_local_element_gids.size() == children_per_element);

//       // once we have the elements in lexicographic order, we have to loop through 
//       // one last time and grab the vertices in lexicographic order too
//       // warning: this is not very performative as far as access patterns, etc go
//       // I'll revisit performance once it's working. in principle, this isn't a
//       // huge computational sink since we're likely not making regions too huge...
//       // we'll see what the timers say
//       {
//         const stk::mesh::FieldBase *coordinatesField = mesh->getMetaData()->get_field(stk::topology::NODE_RANK, "coordinates");
//         const stk::mesh::BucketVector & local_buckets = mesh->getBulkData()->get_buckets(stk::topology::ELEM_RANK,mesh->getMetaData()->locally_owned_part());
//         for (stk::mesh::BucketVector::const_iterator iBucket = local_buckets.begin(); iBucket != local_buckets.end(); ++iBucket)
//         {
//           stk::mesh::Bucket & bucket = **iBucket;

//           const unsigned int num_elements_in_bucket = bucket.size();
//           debug << "bucket size = " << num_elements_in_bucket << std::endl;

//         }
//       }

//       debug << "quasiRegionGIDs = " << quasiRegionGIDs << std::endl;


//     } // if use percept

    // {
    //   std::vector<stk::mesh::Entity> elements;
    //   Kokkos::DynRankView<double,PHX::Device> vertices;
    //   std::vector<std::size_t> localIds;

    //   panzer_stk::workset_utils::getIdsAndVertices(*mesh,element_blocks[0],localIds,vertices); // TODO: in Matthias' case, this is likely eBlock[my_rank]

    //   if(dump_element_vertices)
    //   {
    //     sleep(my_rank);
    //     for(unsigned int ielem=0; ielem<vertices.extent(0); ++ielem)
    //       for(unsigned int ivert=0; ivert<vertices.extent(1); ++ivert)
    //       {
    //         out_root << "rank p=" << my_rank << " element=" << ielem << " vertex=" << ivert << "   (" << vertices(ielem,ivert,0);
    //         for(unsigned int idim=1; idim<vertices.extent(2); ++idim) // fenceposting the output
    //           out_root << ", " << vertices(ielem,ivert,idim);
    //         out_root << ")" << std::endl;
    //       }
    //     return 0;
    //   }
    // }



    // ---------------------------------------------------------------------------------
    // --------------------------------- SETUP RESPONSES -------------------------------
    // ---------------------------------------------------------------------------------
    tm = Teuchos::null;
    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 4 - Other Panzer Setup")));
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > error_response_library = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(workset_container,dof_manager,linear_object_factory));

    {
      // must set a higher integration order for error computations to avoid false superconvergence
      const int integration_order = 5;

      panzer::FunctionalResponse_Builder<int,int> builder;
      builder.comm = MPI_COMM_WORLD;
      builder.cubatureDegree = integration_order;
      builder.requiresCellIntegral = true;
      builder.quadPointField = "TEMPERATURE_L2_ERROR";

      error_response_library->addResponse("L2 Error",element_blocks,builder);

      builder.comm = MPI_COMM_WORLD;
      builder.cubatureDegree = integration_order;
      builder.requiresCellIntegral = true;
      builder.quadPointField = "TEMPERATURE_H1_ERROR";

      error_response_library->addResponse("H1 Error",element_blocks,builder);
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

      closure_models.sublist("solid").sublist("TEMPERATURE_H1_ERROR").set<std::string>("Type","H1 ERROR_CALC");
      closure_models.sublist("solid").sublist("TEMPERATURE_H1_ERROR").set<std::string>("Field A","TEMPERATURE");
      closure_models.sublist("solid").sublist("TEMPERATURE_H1_ERROR").set<std::string>("Field B","TEMPERATURE_EXACT");
      closure_models.sublist("solid").sublist("TEMPERATURE_EXACT").set<std::string>("Type","TEMPERATURE_EXACT");
    }

    Teuchos::ParameterList user_data("User Data"); // user data can be empty here


    // setup field manager builder
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);
    fmb->setWorksetContainer(workset_container);
    fmb->setupVolumeFieldManagers(physics_blocks,cm_factory,closure_models,*linear_object_factory,user_data);
    fmb->setupBCFieldManagers(bcs,physics_blocks,*eqset_factory,cm_factory,bc_factory,closure_models,
                              *linear_object_factory,user_data);
    fmb->writeVolumeGraphvizDependencyFiles("Poisson", physics_blocks);


    // setup assembly engine
    /////////////////////////////////////////////////////////////

    panzer::AssemblyEngine_TemplateManager<panzer::Traits> assembly_engine_tm;
    panzer::AssemblyEngine_TemplateBuilder builder(fmb,linear_object_factory);
    assembly_engine_tm.buildObjects(builder);


    // Finalize construction of STK writer response library
    /////////////////////////////////////////////////////////////
    {
      user_data.set<int>("Workset Size",workset_size);
      error_response_library->buildResponseEvaluators(physics_blocks,
                                                      cm_factory,
                                                      closure_models,
                                                      user_data);
    }


    // assemble linear system
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::LinearObjContainer> ghostCont = linear_object_factory->buildGhostedLinearObjContainer();
    Teuchos::RCP<panzer::LinearObjContainer> container = linear_object_factory->buildLinearObjContainer();
    linear_object_factory->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                                      panzer::LinearObjContainer::F |
                                                      panzer::LinearObjContainer::Mat,*ghostCont);
    linear_object_factory->initializeContainer(panzer::LinearObjContainer::X |
                                               panzer::LinearObjContainer::F |
                                               panzer::LinearObjContainer::Mat,*container);
    ghostCont->initialize();
    container->initialize();

    panzer::AssemblyEngineInArgs input(ghostCont,container);
    input.alpha = 0;
    input.beta = 1;

    // evaluate physics: This does both the Jacobian and residual at once
    assembly_engine_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);



    // ---------------------------------------------------------------------------------
    // ---------------------------------- LINEAR SOLVER --------------------------------
    // ---------------------------------------------------------------------------------
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
      out_all << "Solution local length: " << tp_container->get_x()->getLocalLength() << std::endl;
      out_root << "Solution norm: " << tp_container->get_x()->norm2() << std::endl;
    }


    // ---------------------------------------------------------------------------------
    // ---------------------------------- REGION DRIVER --------------------------------
    // ---------------------------------------------------------------------------------
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::ArrayRCP;
      using Teuchos::TimeMonitor;
      using Teuchos::ParameterList;

      // =========================================================================
      // Convenient definitions
      // =========================================================================
      using STS = Teuchos::ScalarTraits<SC>;
      SC zero = STS::zero(), one = STS::one();
      using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
      using real_type = typename STS::coordinateType;
      using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;


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
        if(unstructuredRanks[idx] == my_rank) {useUnstructured = true;}
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
      RCP<panzer::TpetraLinearObjFactory<panzer::Traits,ST,LO,GO> > tp_object_factory = rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,ST,LO,GO>(comm, dof_manager));
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

      // std::cout << "p=" << my_rank << " | numSend=" << numSend << std::endl;
      // << ", numReceive=" << numReceive << std::endl;
      // std::cout << "p=" << my_rank << " | receiveGIDs: " << receiveGIDs << std::endl;
      // std::cout << "p=" << my_rank << " | receivePIDs: " << receivePIDs << std::endl;
      // std::cout << "p=" << my_rank << " | sendGIDs: " << sendGIDs << std::endl;
      // std::cout << "p=" << my_rank << " | sendPIDs: " << sendPIDs << std::endl;

      // Second we actually fill the send and receive arrays with appropriate data
      // which will allow us to compute the region and composite maps.
      // Now we can construct a list of GIDs that corresponds to rowMap
      Array<LO>  interfacesDimensions, interfacesLIDs;
      if(useUnstructured) {
        findInterface(num_dimensions, rNodesPerDim, boundaryConditions,
                      interfacesDimensions, interfacesLIDs);

        // std::cout << "p=" << my_rank << " | numLocalRegionNodes=" << numLocalRegionNodes
        //           << ", rNodesPerDim: " << rNodesPerDim << std::endl;
        // std::cout << "p=" << my_rank << " | boundaryConditions: " << boundaryConditions << std::endl
        //           << "p=" << my_rank << " | rNodesPerDim: " << rNodesPerDim << std::endl
        //           << "p=" << my_rank << " | interfacesDimensions: " << interfacesDimensions << std::endl
        //           << "p=" << my_rank << " | interfacesLIDs: " << interfacesLIDs << std::endl;
      }

      interfaceParams->set<Array<LO> >("interfaces: nodes per dimensions", interfacesDimensions); // nodesPerDimensions);
      interfaceParams->set<Array<LO> >("interfaces: interface nodes",      interfacesLIDs); // interfaceLIDs);

      // std::cout << "p=" << my_rank << " | compositeToRegionLIDs: " << compositeToRegionLIDs << std::endl;
      // std::cout << "p=" << my_rank << " | quasiRegionGIDs: " << quasiRegionGIDs << std::endl;
      // std::cout << "p=" << my_rank << " | interfaceGIDs: " << interfaceGIDs << std::endl;
      // std::cout << "p=" << my_rank << " | interfaceLIDsData: " << interfaceLIDsData << std::endl;
      // std::cout << "p=" << my_rank << " | interfaceLIDs: " << interfaceLIDs << std::endl;
      // std::cout << "p=" << my_rank << " | quasiRegionCoordGIDs: " << quasiRegionCoordGIDs() << std::endl;

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
      Teuchos::RCP<Import> rowImport;
      Teuchos::RCP<Import> colImport;
      rowImport = ImportFactory::Build(dofMap, rowMap);
      colImport = ImportFactory::Build(dofMap, colMap);
      Teuchos::RCP<Import> coordImporter = ImportFactory::Build(nodeMap, quasiRegCoordMap);

      comm->barrier();
      tmLocal = Teuchos::null;
      tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.3 - Import ghost GIDs")));

      Array<GO>  interfaceCompositeGIDs, interfaceRegionGIDs;
      ExtractListOfInterfaceRegionGIDs(revisedRowMap, interfaceLIDsData, interfaceRegionGIDs);

      Teuchos::RCP<Xpetra::MultiVector<LO, LO, GO, NO> > regionsPerGIDWithGhosts;
      Teuchos::RCP<Xpetra::MultiVector<GO, LO, GO, NO> > interfaceGIDsMV;
      MakeRegionPerGIDWithGhosts(nodeMap, revisedRowMap, rowImport,
                                 maxRegPerGID, numDofsPerNode,
                                 lNodesPerDim, sendGIDs, sendPIDs, interfaceLIDsData,
                                 regionsPerGIDWithGhosts, interfaceGIDsMV);

      Teuchos::ArrayRCP<LO> regionMatVecLIDs;
      Teuchos::RCP<Import> regionInterfaceImporter;
      SetupMatVec(interfaceGIDsMV, regionsPerGIDWithGhosts, revisedRowMap, rowImport,
                  regionMatVecLIDs, regionInterfaceImporter);

      comm->barrier();
      tmLocal = Teuchos::null;
      tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.4 - Build QuasiRegion Matrix")));

      std::cout << "About to create quasi region matrix" << std::endl;
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > quasiRegionMats;
      MakeQuasiregionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A),
                              regionsPerGIDWithGhosts, rowMap, colMap, rowImport,
                              quasiRegionMats, regionMatVecLIDs);
      std::cout << "Done creating quasi region matrix" << std::endl;

      comm->barrier();
      tmLocal = Teuchos::null;
      tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.5 - Build Region Matrix")));

      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats;
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
      Teuchos::RCP<MultiVector> regionNullspace;
      Teuchos::RCP<RealValuedMultiVector> regionCoordinates;

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
      Teuchos::RCP<ParameterList> coarseSolverData = rcp(new ParameterList());
      coarseSolverData->set<std::string>("coarse solver type", coarseSolverType);
      coarseSolverData->set<bool>("coarse solver rebalance", coarseSolverRebalance);
      coarseSolverData->set<int>("coarse rebalance num partitions", rebalanceNumPartitions);
      coarseSolverData->set<std::string>("amg xml file", coarseAmgXmlFile);
      coarseSolverData->set<std::string>("smoother xml file", coarseSmootherXMLFile);
      Teuchos::RCP<ParameterList> hierarchyData = rcp(new ParameterList());


      // Create MueLu Hierarchy Initially...
      // Read MueLu parameter list form xml file
      Teuchos::RCP<ParameterList> mueluParams = Teuchos::rcp(new ParameterList());
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
      Teuchos::RCP<Hierarchy> regHierarchy  = MueLu::CreateXpetraPreconditioner(regionMats, *mueluParams);

      {
        Teuchos::RCP<MueLu::Level> level = regHierarchy->GetLevel(0);
        level->Set<Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("rowImport",rowImport);
        level->Set<ArrayView<LocalOrdinal> > ("compositeToRegionLIDs", compositeToRegionLIDs() );
        level->Set<Teuchos::RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >("interfaceGIDs", interfaceGIDsMV);
        level->Set<Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >("regionsPerGIDWithGhosts", regionsPerGIDWithGhosts);
        level->Set<Teuchos::ArrayRCP<LocalOrdinal> >("regionMatVecLIDs", regionMatVecLIDs);
        level->Set<Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("regionInterfaceImporter", regionInterfaceImporter);
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
        Teuchos::RCP<MueLu::Level> level = regHierarchy->GetLevel(levelIdx);
        Teuchos::RCP<Xpetra::Import<LO, GO, NO> > regionInterfaceImport = level->Get<Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("regionInterfaceImporter");
        Teuchos::ArrayRCP<LO>            regionMatVecLIDs1     = level->Get<Teuchos::ArrayRCP<LO> >("regionMatVecLIDs");
        smootherParams[levelIdx]->set("Fast MatVec: interface LIDs",
                                      regionMatVecLIDs1);
        smootherParams[levelIdx]->set("Fast MatVec: interface importer",
                                      regionInterfaceImport);
      }

      // RCP<Teuchos::FancyOStream> fancy2 = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      // Teuchos::FancyOStream& out2 = *fancy2;
      // for(LO levelIdx = 0; levelIdx < numLevels; ++levelIdx) {
      //   out2 << "p=" << my_rank << " | regionMatVecLIDs on level " << levelIdx << std::endl;
      //   regionMatVecLIDsPerLevel[levelIdx]->describe(out2, Teuchos::VERB_EXTREME);
      // }

      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 5 - Solve with V-cycle")));

      {
        //    std::cout << my_rank << " | Running V-cycle ..." << std::endl;

        TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

        // We first use the non-level container variables to setup the fine grid problem.
        // This is ok since the initial setup just mimics the application and the outer
        // Krylov method.
        //
        // We switch to using the level container variables as soon as we enter the
        // recursive part of the algorithm.
        //

        // Composite residual vector
        Teuchos::RCP<Vector> compRes = VectorFactory::Build(dofMap, true);

        // transform composite vectors to regional layout
        Teuchos::RCP<Vector> quasiRegX;
        Teuchos::RCP<Vector> regX;
        compositeToRegional(X, quasiRegX, regX,
                            revisedRowMap, rowImport);

        Teuchos::RCP<Vector> quasiRegB;
        Teuchos::RCP<Vector> regB;
        compositeToRegional(B, quasiRegB, regB,
                            revisedRowMap, rowImport);
#ifdef DUMP_LOCALX_AND_A
        FILE *fp;
        char str[80];
        sprintf(str,"theMatrix.%d",my_rank);
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
        sprintf(str,"theX.%d",my_rank);
        fp = fopen(str,"w");
        ArrayRCP<SC> lX= regX->getDataNonConst(0);
        for (size_t kkk = 0; kkk < regionMats->getNodeNumRows(); kkk++) fprintf(fp, "%22.16e\n",lX[kkk]);
        fclose(fp);
#endif

        Teuchos::RCP<Vector> regRes;
        regRes = VectorFactory::Build(revisedRowMap, true);

        /////////////////////////////////////////////////////////////////////////
        // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
        /////////////////////////////////////////////////////////////////////////

        // Prepare output of residual norm to file
        Teuchos::RCP<std::ofstream> log;
        if (my_rank == 0)
        {
          log = rcp(new std::ofstream(convergenceLog.c_str()));
          (*log) << "# num procs = " << dofMap->getComm()->getSize() << "\n"
              << "# iteration | res-norm (scaled=" << scaleResidualHist << ")\n"
              << "#\n";
          *log << std::setprecision(16) << std::scientific;
        }

        // Print type of residual norm to the screen
        if (scaleResidualHist)
          out_root << "Using scaled residual norm." << std::endl;
        else
          out_root << "Using unscaled residual norm." << std::endl;


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
          Teuchos::RCP<MueLu::Level> level = regHierarchy->GetLevel(0);
          Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal> > regInterfaceScalings = level->Get<Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal> > >("regInterfaceScalings");
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
            out_root << cycle << "\t" << normRes << std::endl;
            if (my_rank == 0)
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
        out_root << "Number of iterations performed for this solve: " << cycle << std::endl;

        std::cout << std::setprecision(old_precision);
        std::cout.unsetf(std::ios::fixed | std::ios::scientific);
      }

      comm->barrier();
      tm = Teuchos::null;
      globalTimeMonitor = Teuchos::null;

      if (showTimerSummary)
      {
        Teuchos::RCP<ParameterList> reportParams = rcp(new ParameterList);
        const std::string filter = "";
        if (useStackedTimer) {
          Teuchos::StackedTimer::OutputOptions options;
          options.output_fraction = options.output_histogram = options.output_minmax = true;
          stacked_timer->report(out_root, comm, options);
        } else {
          std::ios_base::fmtflags ff(out_root.flags());
          TimeMonitor::report(comm.ptr(), out_root, filter, reportParams);
          out_root << std::setiosflags(ff);
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
      linear_object_factory->globalToGhostContainer(*container,*ghostCont, panzer::TpetraLinearObjContainer<ST,LO,GO>::X
                                                    | panzer::TpetraLinearObjContainer<ST,LO,GO>::DxDt);

      // get X Tpetra_Vector from ghosted container
      // TODO: there is some magic here with Tpetra objects that needs to be fixed
      //Teuchos::RCP<panzer::TpetraLinearObjContainer<ST,LO,GO> > tp_ghostCont = Teuchos::rcp_dynamic_cast<panzer::TpetraLinearObjContainer<ST,LO,GO> >(ghostCont);
      //panzer_stk::write_solution_data(*dof_manager,*mesh,*tp_ghostCont->get_x());

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

      Teuchos::RCP<panzer::ResponseBase> l2_resp = error_response_library->getResponse<panzer::Traits::Residual>("L2 Error");
      Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > l2_resp_func = Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(l2_resp);
      Teuchos::RCP<Thyra::VectorBase<double> > l2_respVec = Thyra::createMember(l2_resp_func->getVectorSpace());
      l2_resp_func->setVector(l2_respVec);


     Teuchos::RCP<panzer::ResponseBase> h1_resp = error_response_library->getResponse<panzer::Traits::Residual>("H1 Error");
     Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > h1_resp_func = Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(h1_resp);
     Teuchos::RCP<Thyra::VectorBase<double> > h1_respVec = Thyra::createMember(h1_resp_func->getVectorSpace());
     h1_resp_func->setVector(h1_respVec);


      error_response_library->addResponsesToInArgs<panzer::Traits::Residual>(respInput);
      error_response_library->evaluate<panzer::Traits::Residual>(respInput);

      out_root << "This is the Basis Order" << std::endl;
      out_root << "Basis Order = " << discretization_order << std::endl;
      out_root << "This is the L2 Error" << std::endl;
      out_root << "L2 Error = " << sqrt(l2_resp_func->value) << std::endl;
      out_root << "This is the H1 Error" << std::endl;
      out_root << "H1 Error = " << sqrt(h1_resp_func->value) << std::endl;
    }

    tm = Teuchos::null;
    globalTimeMonitor = Teuchos::null;

    if (showTimerSummary)
    {
      Teuchos::RCP<ParameterList> reportParams = rcp(new ParameterList);
      const std::string filter = "";
      if (useStackedTimer)
      {
        Teuchos::StackedTimer::OutputOptions options;
        options.output_fraction = options.output_histogram = options.output_minmax = true;
        stacked_timer->report(out_root, comm, options);
      }
      else
      {
        std::ios_base::fmtflags ff(out_root.flags());
        Teuchos::TimeMonitor::report(comm.ptr(), out_root, filter, reportParams);
        out_root << std::setiosflags(ff);
      }
    }

  } // Parallel initialization scope
  Tpetra::finalize();
  Kokkos::finalize();

  return 0;
} // main
