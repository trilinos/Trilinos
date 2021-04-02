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
#include "MueLu_CreateTpetraPreconditioner.hpp"
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

    // GH: comb back through everything and make sure I'm using my MPI comms properly when necessary
    Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
    Teuchos::RCP<const Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

    const int numRanks = comm->getSize();
    const int myRank = comm->getRank();

    if (myRank == 0) {
      *out << "Running MueLu region driver on " << numRanks << " ranks... \n";
    }

    // Parse command line arguments
    Teuchos::CommandLineProcessor clp(false);
    std::string exodusFileName        = "";                  clp.setOption("exodus-mesh",           &exodusFileName,          "Exodus hex mesh filename (overrides a pamgen-mesh if both specified)");
    std::string pamgenFileName        = "cylinder.rtp";      clp.setOption("pamgen-mesh",           &pamgenFileName,          "Pamgen hex mesh filename");
    std::string xmlFileName           = "";                  clp.setOption("xml",                   &xmlFileName,             "MueLu parameters file");
    int mesh_refinements              = 1;                   clp.setOption("mesh-refinements",      &mesh_refinements,        "Uniform mesh refinements");
    bool delete_parent_elements       = false;                clp.setOption("delete-parent-elements", "keep-parent-elements", &delete_parent_elements,"Save the parent elements in the perceptMesh");

    // debug options
    bool print_percept_mesh           = false;                clp.setOption("print-percept-mesh", "no-print-percept-mesh", &print_percept_mesh, "Calls perceptMesh's print_info routine");
    bool print_debug_info             = false;                clp.setOption("print-debug-info", "no-print-debug-info", &print_debug_info, "Print more debugging information");

    // timer options
    bool useStackedTimer              = false;                clp.setOption("stacked-timer","no-stacked-timer", &useStackedTimer, "use stacked timer");
    bool showTimerSummary             = false;                clp.setOption("show-timer-summary", "no-show-timer-summary", &showTimerSummary, "Switch on/off the timer summary at the end of the run.");

    clp.recogniseAllOptions(true);
    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }
    TEUCHOS_ASSERT(mesh_refinements >= 0); // get around the expected 7 arguments error for clp.setOption(...unsigned int...) for now

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
    /******************************* MESH AND WORKSETS ********************************/
    /**********************************************************************************/


    Teuchos::RCP<Teuchos::Time> meshTimer = Teuchos::TimeMonitor::getNewCounter("Step 1: Mesh generation");
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

      if(mesh_refinements>0)
      {
        mesh_pl->set("Levels of Uniform Refinement",mesh_refinements); // this multiplies the number of elements by 2^(dimension*level)
        mesh_pl->set("Keep Percept Data",true); // this is necessary to gather mesh hierarchy information
        mesh_pl->set("Keep Percept Parent Elements",!delete_parent_elements); // this is necessary to gather mesh hierarchy information
      }
    }
    else if(pamgenFileName.length())
    {
      // set the filename and type
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory);
      mesh_pl->set("File Name",pamgenFileName);
      mesh_pl->set("File Type","Pamgen");

      if(mesh_refinements>0)
      {
        mesh_pl->set("Levels of Uniform Refinement",mesh_refinements); // this multiplies the number of elements by 2^(dimension*level)
        mesh_pl->set("Keep Percept Data",true); // this is necessary to gather mesh hierarchy information
        mesh_pl->set("Keep Percept Parent Elements",!delete_parent_elements); // this is necessary to gather mesh hierarchy information
      }
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

    unsigned int dimension = mesh->getDimension();
    if(print_debug_info)
      std::cout << "Using dimension = " << dimension << std::endl;

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


    /**********************************************************************************/
    /********************************** CONSTRUCT REGIONS *****************************/
    /**********************************************************************************/

    unsigned int children_per_element = 1 << (dimension*mesh_refinements);
    if(print_debug_info)
      std::cout << "Number of mesh children = " << children_per_element << std::endl;

    // initialize data here that we will use for the region MG solver
    std::vector<GO> child_element_gids; // these don't always start at 0, and changes I'm making to Panzer keep changing this, so I'll store them for now
    std::vector<GO> child_element_region_gids;

    // do not run region MG if we delete parent elements or if we do not refine the mesh regularly
    if(mesh_refinements>0 && !delete_parent_elements)
    {
      // get the Percept mesh from Panzer
      Teuchos::RCP<percept::PerceptMesh> refinedMesh = mesh->getRefinedMesh();
      if(print_percept_mesh)
        refinedMesh->print_info(*out,"",1,true);

      //    Teuchos::RCP<percept::URP_Heterogeneous_3D> breakPattern = Teuchos::rcp(new percept::URP_Heterogeneous_3D(*refinedMesh));
      //    percept::UniformRefiner breaker(*refinedMesh,*breakPattern);
      //    breaker.deleteParentElements();

      // ids are linear within stk, but we need an offset because the original mesh info comes first
      size_t node_id_start = 0;
      {
        const stk::mesh::BucketVector & buckets = refinedMesh->get_bulk_data()->buckets(refinedMesh->node_rank());
        stk::mesh::Bucket & bucket = **buckets.begin() ;
        node_id_start = refinedMesh->id(bucket[0]);
        if(print_debug_info)
          std::cout << "Starting node id = " << node_id_start << std::endl;
      }

      size_t elem_id_start = 0;
      {
        const stk::mesh::BucketVector & buckets = refinedMesh->get_bulk_data()->buckets(refinedMesh->element_rank());
        stk::mesh::Bucket & bucket = **buckets.begin() ;
        elem_id_start = refinedMesh->id(bucket[0]);
        if(print_debug_info)
          std::cout << "Starting element id = " << elem_id_start << std::endl;
      }


      {
        const stk::mesh::BucketVector & buckets = refinedMesh->get_bulk_data()->buckets(refinedMesh->element_rank());
        int npar=0;
        int nchild=0;
        for (stk::mesh::BucketVector::const_iterator k = buckets.begin(); k != buckets.end(); ++k)
        {
          stk::mesh::Bucket & bucket = **k ;
          if(print_debug_info)
            std::cout << "New bucket" << std::endl;

          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
          {
            stk::mesh::Entity element = bucket[iElement];
            if (!refinedMesh->isParentElement(element, false))
            {
              ++nchild;

              // this is the important part here. take the id of the element and the id of the element's root
              child_element_gids.push_back(refinedMesh->id(element));
              child_element_region_gids.push_back(refinedMesh->id(refinedMesh->rootOfTree(element)));

              if(print_debug_info)
                std::cout << "Stk Element = " << element << std::endl;

              percept::MyPairIterRelation elem_nodes ( *refinedMesh, element,  stk::topology::NODE_RANK);

              for (unsigned i_node = 0; i_node < elem_nodes.size(); i_node++)
              {
                stk::mesh::Entity node = elem_nodes[i_node].entity();
                if(print_debug_info)
                  std::cout << "Stk Node = " << node << std::endl;

                //std::cout << "node " << node << std::endl;
                //        double *f_data = refinedMesh->field_data(field, node, null_u);
                //        for (int i_stride=0; i_stride < refinedMesh->get_spatial_dim(); i_stride++)
                //        {
                //        std::cout << f_data[i_stride] << " ";
                //      }
                //          std::cout << std::endl;
              }
            }
            else
            {
              std::cout << "parent= " << refinedMesh->id(element) << std::endl;
            }
          }
        }
      }
    }

    if(print_debug_info)
    {
      for(unsigned int i=0; i<child_element_gids.size(); ++i)
      {
        std::cout << "child= " << child_element_gids[i] << " parent= " << child_element_region_gids[i] << std::endl;
      }
    }

    // next we need to get the LIDs

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

      // GH: comment out H1 errors for now while I work on other aspects
      /*
      builder.comm = MPI_COMM_WORLD;
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

      // GH: comment out H1 errors for now while I work on other aspects
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


    // solve the linear system
    /////////////////////////////////////////////////////////////

    // convert generic linear object container to tpetra container
    Teuchos::RCP<panzer::TpetraLinearObjContainer<ST,LO,GO> > tp_container = Teuchos::rcp_dynamic_cast<panzer::TpetraLinearObjContainer<ST,LO,GO> >(container);

    Teuchos::RCP<MueLu::TpetraOperator<ST,LO,GO,NT> > mueLuPreconditioner;

    if(xmlFileName.size())
      mueLuPreconditioner = MueLu::CreateTpetraPreconditioner(Teuchos::rcp_dynamic_cast<Tpetra::Operator<ST,LO,GO,NT> >(tp_container->get_A()), xmlFileName);
    else
    {
      Teuchos::ParameterList mueLuParamList;
      mueLuParamList.set("verbosity", "low");
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
    std::cout << "Solution local length: " << tp_container->get_x()->getLocalLength() << std::endl;
    std::cout << "Solution norm: " << tp_container->get_x()->norm2() << std::endl;

    // output data (optional)
    /////////////////////////////////////////////////////////////

    // write the solution to matrix
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

    // compute the error of the finite element solution
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
