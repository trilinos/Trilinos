#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"

#include "Panzer_NodeType.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_InitialCondition_Builder.hpp"
#include "Panzer_CheckBCConsistency.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"
#include "Panzer_Response_Functional.hpp"

#include "Panzer_STK_MeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_SetupLOWSFactory.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_IOClosureModel_Factory_TemplateBuilder.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"

// includes specific to this example
#include "Thyra_BlockedLinearOpWithSolveBase.hpp"
#include "Panzer_STK_ParameterListCallbackBlocked.hpp"
#include "PanzerMiniEM_config.hpp"
#include "MiniEM_EquationSetFactory.hpp"
#include "MiniEM_BCStrategy_Factory.hpp"
#include "MiniEM_ClosureModel_Factory_TemplateBuilder.hpp"
#include "MiniEM_AddFieldsToMesh.hpp"
#include "MiniEM_OperatorRequestCallback.hpp"
#include "MiniEM_FullMaxwellPreconditionerFactory.hpp"
#include "MiniEM_DiscreteGradient.hpp"

#include <string>
#include <iostream>


void createExodusFile(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                      Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory,
                      Teuchos::RCP<panzer_stk::STK_Interface> mesh,
                      const bool & exodus_out);

Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >
buildSTKIOResponseLibrary(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
    const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linObjFactory,
    const Teuchos::RCP<panzer::WorksetContainer> & wkstContainer,
    const Teuchos::RCP<panzer::UniqueGlobalIndexerBase> & globalIndexer,
    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
    const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
    const Teuchos::ParameterList & closure_model_pl);

template <class Scalar>
void writeToExodus(double time_stamp,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & x,
    const panzer::ModelEvaluator<Scalar> & model,
    panzer::ResponseLibrary<panzer::Traits> & stkIOResponseLibrary,
    panzer_stk::STK_Interface & mesh);


void updateParams(const std::string & xml,
                  Teuchos::RCP<Teuchos::ParameterList> pl,
                  const Teuchos::RCP<const Teuchos::MpiComm<int> > comm,
                  const Teuchos::RCP<Teuchos::FancyOStream> out) {
  *out << "Loading solver config from " << xml << std::endl;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xml,pl.ptr(),*comm);
}


/********************************************************************************
 * Sets up an electromagetics problem driven by a simple Gaussian current pulse
 * on the domain [0,1]^3. First order Maxwell equations with edge-face
 * discretization for E,B. Backward Euler time-stepping with fixed CFL. Linear
 * systems solved with Belos GMRES using augmentation based block preconditioner
 * through Teko with multigrid subsolves from MueLu.
 *
 * This is meant to test the components of the Tpetra linear solver stack
 * required by EMPIRE-EM
 *
 * Command-line arguments:
 *   x-elements, y-elements,z-elements: # elements in each direction
 *   basis-order: order of finite element bases
 *   cfl: CFL with respect to speed of light. CFL = c*dt/min{dx,dy,dz}
 *   workset-size: size of worksets
 *
 * ******************************************************************************
 */

enum solverType {
  AUGMENTATION,
  MUELU_REFMAXWELL,
  ML_REFMAXWELL,
  CG
};

enum linearAlgebraType {
  linAlgTpetra,
  linAlgEpetra
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class blockedLinObjFactory>
int main_(Teuchos::CommandLineProcessor &clp, int argc,char * argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
  Teuchos::RCP<const Teuchos::MpiComm<int> > comm
    = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());
  if (comm->getSize() > 1) {
    out->setOutputToRootOnly(0);
  }

  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
  bool use_stacked_timer;

  {
    // defaults for command-line options
    int x_elements=-1,y_elements=-1,z_elements=-1,basis_order=1;
    double dt=0.0;
    std::string meshFile = "";
    bool exodus_output = false;
    bool matrix_output = false;
    std::string input_file = "maxwell.xml";
    std::string xml = "";
    solverType solverValues[4] = {AUGMENTATION, MUELU_REFMAXWELL, ML_REFMAXWELL, CG};
    const char * solverNames[4] = {"Augmentation", "MueLu-RefMaxwell", "ML-RefMaxwell", "CG"};
    solverType solver = MUELU_REFMAXWELL;
    int numTimeSteps = 1;
    bool doSolveTimings = false;
    int numReps = 0;
    linearAlgebraType linAlgebraValues[2] = {linAlgTpetra, linAlgEpetra};
    const char * linAlgebraNames[2] = {"Tpetra", "Epetra"};
    linearAlgebraType linAlgebra = linAlgTpetra;
    clp.setOption<linearAlgebraType>("linAlgebra",&linAlgebra,2,linAlgebraValues,linAlgebraNames);
    use_stacked_timer = true;
    clp.setOption("x-elements",&x_elements);
    clp.setOption("y-elements",&y_elements);
    clp.setOption("z-elements",&z_elements);
    clp.setOption("basis-order",&basis_order);
    clp.setOption("dt",&dt,"Override \"dt\" specified in input file.");
    clp.setOption("meshFile",&meshFile,"Override input mesh file specified in input file.");
    clp.setOption("exodus-output","no-exodus-output",&exodus_output);
    clp.setOption("matrix-output","no-matrix-output",&matrix_output);
    clp.setOption("inputFile",&input_file,"XML file with the problem definitions");
    clp.setOption("solverFile",&xml,"XML file with the solver params");
    clp.setOption<solverType>("solver",&solver,4,solverValues,solverNames,"Solver that is used");
    clp.setOption("numTimeSteps",&numTimeSteps);
    clp.setOption("doSolveTimings","no-doSolveTimings",&doSolveTimings,"repeat the first solve \"numTimeSteps\" times");
    clp.setOption("stacked-timer","no-stacked-timer",&use_stacked_timer,"Run with or without stacked timer output");

    // parse command-line argument
    clp.recogniseAllOptions(true);
    clp.throwExceptions(false);
    const Teuchos::CommandLineProcessor::EParseCommandLineReturn parseResult = clp.parse (argc, argv);
    switch (parseResult) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    if (use_stacked_timer)
      stacked_timer = rcp(new Teuchos::StackedTimer("Mini-EM"));
    Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: Total Time")));

    if (doSolveTimings) {
      numReps = numTimeSteps;
      numTimeSteps = 1;
    }

    Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
    Scalar one = Teuchos::ScalarTraits<Scalar>::one();

    // data container for auxiliary linear operators used in preconditioning (mass matrix and gradient)
    Teuchos::RCP<panzer::GlobalEvaluationDataContainer> auxGlobalData = Teuchos::rcp(new panzer::GlobalEvaluationDataContainer);

    // factory definitions
    Teuchos::RCP<mini_em::EquationSetFactory> eqset_factory =
      Teuchos::rcp(new mini_em::EquationSetFactory(auxGlobalData)); // where the maxwell equations are defined
    mini_em::BCFactory bc_factory;                                  // where boundary conditions are defined

    // GlobalData sets ostream and parameter interface to physics
    Teuchos::RCP<panzer::GlobalData> globalData = panzer::createGlobalData();


    // Parse the input file and broadcast to other processes
    Teuchos::RCP<Teuchos::ParameterList> input_params = Teuchos::rcp(new Teuchos::ParameterList("User_App Parameters"));
    Teuchos::updateParametersFromXmlFileAndBroadcast(input_file, input_params.ptr(), *comm);
    Teuchos::ParameterList & mesh_pl                 = input_params->sublist("Mesh");
    RCP<Teuchos::ParameterList> physicsBlock_pl      = rcp(new Teuchos::ParameterList(input_params->sublist("Physics Blocks")));
    Teuchos::ParameterList & assembly_pl             = input_params->sublist("Assembly");
    Teuchos::ParameterList & block_to_physics_pl     = input_params->sublist("Block ID to Physics ID Mapping");
    Teuchos::ParameterList & block_to_aux_physics_pl = input_params->sublist("Block ID to Auxiliary Physics ID Mapping");
    Teuchos::ParameterList & bcs_pl                  = input_params->sublist("Boundary Conditions");
    Teuchos::ParameterList & aux_bcs_pl              = input_params->sublist("Auxiliary Boundary Conditions");
    Teuchos::ParameterList & closure_models          = input_params->sublist("Closure Models");
    Teuchos::ParameterList responses                 = input_params->sublist("Responses");
    Teuchos::ParameterList & user_data               = input_params->sublist("User Data");

    RCP<panzer_stk::STK_Interface> mesh;
    int dim = 3;
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
    {
      Teuchos::TimeMonitor tMmesh(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: build mesh")));

      if (mesh_pl.get<std::string>("Source") == "Exodus File" || meshFile != "") { // Exodus file reader...
        RCP<Teuchos::ParameterList> pl;
        RCP<Teuchos::ParameterList> input_pl = rcp(new Teuchos::ParameterList(mesh_pl.sublist("Exodus File")));
        if (meshFile == "") {
          pl = rcp(new Teuchos::ParameterList(input_pl->sublist("Exodus Parameters")));
        } else {
          pl = rcp(new Teuchos::ParameterList());
          pl->set("File Name",meshFile);
        }
        mesh_factory = Teuchos::RCP<panzer_stk::STK_MeshFactory>(new panzer_stk::STK_ExodusReaderFactory());
        mesh_factory->setParameterList(pl);
        // build mesh
        mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
        //get dt
        if (dt <=0.)
          dt = input_pl->get<double>("dt");
      } else if (mesh_pl.get<std::string>("Source") == "Inline Mesh") { // Inline mesh generator
        // set mesh factory parameters
        Teuchos::ParameterList & inline_gen_pl = mesh_pl.sublist("Inline Mesh");
        RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList(inline_gen_pl.sublist("Mesh Factory Parameter List")));
        dim = inline_gen_pl.get<int>("Mesh Dimension");

        // overrides from command line
        if (x_elements > 0)
          pl->set<int>("X Elements",x_elements);
        if (y_elements > 0)
          pl->set<int>("Y Elements",y_elements);
        if (dim == 3 && z_elements > 0)
          pl->set<int>("Z Elements",z_elements);

        // build mesh
        if (dim == 3) {
          if (inline_gen_pl.get<std::string>("Mesh Type") == "tet")
            mesh_factory = rcp(new panzer_stk::CubeTetMeshFactory());
          else if (inline_gen_pl.get<std::string>("Mesh Type") == "quad")
            mesh_factory = rcp(new panzer_stk::CubeHexMeshFactory());
          else
            return EXIT_FAILURE;
        } else if (dim == 2) {
          if (inline_gen_pl.get<std::string>("Mesh Type") == "tet")
            mesh_factory = rcp(new panzer_stk::SquareTriMeshFactory());
          else if (inline_gen_pl.get<std::string>("Mesh Type") == "quad")
            mesh_factory = rcp(new panzer_stk::SquareQuadMeshFactory());
          else
            return EXIT_FAILURE;
        }
        mesh_factory->setParameterList(pl);
        mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);

        // set dt
        if (dt <=0.) {
          if (inline_gen_pl.isType<double>("dt"))
            dt = inline_gen_pl.get<double>("dt");
          else {
            double cfl = inline_gen_pl.get<double>("CFL");
            x_elements = pl->get<int>("X Elements");
            y_elements = pl->get<int>("Y Elements");
            double min_dx;
            if (dim == 3) {
              z_elements = pl->get<int>("Z Elements");
              min_dx = 1.0/std::max(x_elements,std::max(y_elements,z_elements));
            } else
              min_dx = std::max(x_elements,y_elements);
            // This is only correct when epsilon==epsilon0, mu==mu0
            double c = 299792458.0;
            dt = cfl * min_dx / c;
          }
        }
      } else
        return EXIT_FAILURE;
    }

    *out << "dt = " << dt << std::endl << std::endl;
    if (dt <= 0.0)
      return EXIT_FAILURE;

    RCP<Teuchos::ParameterList> lin_solver_pl = Teuchos::rcp(new Teuchos::ParameterList("Linear Solver"));
    if (xml == "") {
      // Load a solver configuration
      // This input deck choice depends on
      // * chosen solver
      // * linear algebra library
      // * spatial dimension
      // * node type
      if (solver == AUGMENTATION)
        if (linAlgebra == linAlgTpetra)
          updateParams("solverAugmentation.xml", lin_solver_pl, comm, out);
        else
          updateParams("solverAugmentationEpetra.xml", lin_solver_pl, comm, out);
      else if (solver == CG)
        if (linAlgebra == linAlgTpetra)
          updateParams("solverCG.xml", lin_solver_pl, comm, out);
        else
          return EXIT_FAILURE;
      else if (solver == ML_REFMAXWELL)
        updateParams("solverMLRefMaxwell.xml", lin_solver_pl, comm, out);
      else if (solver == MUELU_REFMAXWELL) {
        if (linAlgebra == linAlgTpetra) {
          updateParams("solverMueLuRefMaxwell.xml", lin_solver_pl, comm, out);

          if (dim == 2)
            updateParams("solverMueLuRefMaxwell2D.xml", lin_solver_pl, comm, out);

#ifdef KOKKOS_ENABLE_OPENMP
          if (typeid(panzer::TpetraNodeType).name() == typeid(Kokkos::Compat::KokkosOpenMPWrapperNode).name()) {
            if (linAlgebra == linAlgTpetra)
              updateParams("solverMueLuRefMaxwellOpenMP.xml", lin_solver_pl, comm, out);
            else {
              std::cout << std::endl
                        << "WARNING" << std::endl
                        << "MueLu RefMaxwell + Epetra + OpenMP does currently not work." << std::endl
                        << "The Xpetra-Epetra interface is missing \"setAllValues\" with kokkos views." << std::endl << std::endl;
              return EXIT_FAILURE;
            }
          }
#endif
#ifdef KOKKOS_ENABLE_CUDA
          if (typeid(panzer::TpetraNodeType).name() == typeid(Kokkos::Compat::KokkosCudaWrapperNode).name())
            updateParams("solverMueLuRefMaxwellCuda.xml", lin_solver_pl, comm, out);
#endif
        } else {
          updateParams("solverMueLuRefMaxwellEpetra.xml", lin_solver_pl, comm, out);

          if (dim == 2)
            updateParams("solverMueLuRefMaxwell2D.xml", lin_solver_pl, comm, out);
        }
      }
    } else
      updateParams(xml, lin_solver_pl, comm, out);
    lin_solver_pl->print(*out);

    // The curl-curl term needs to be scaled by dt, the RefMaxwell augmentation needs 1/dt
    closure_models.sublist("electromagnetics_aux").sublist("dt").set<double>("Value",dt);
    closure_models.sublist("electromagnetics_aux").sublist("1/dt").set<double>("Value",1.0/dt);
    closure_models.sublist("electromagnetics").sublist("1/dt").set<double>("Value",1.0/dt);
    if (closure_models.sublist("electromagnetics").sublist("CURRENT").get<std::string>("Type") == "GAUSSIAN PULSE")
      closure_models.sublist("electromagnetics").sublist("CURRENT").set<double>("dt",dt); // set pulse width such that dt resolves it

    std::map<std::string,std::string> block_ids_to_physics_ids, block_ids_to_aux_physics_ids;
    panzer::buildBlockIdToPhysicsIdMap(block_ids_to_physics_ids, block_to_physics_pl);
    panzer::buildBlockIdToPhysicsIdMap(block_ids_to_aux_physics_ids, block_to_aux_physics_pl);

    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
    for(auto itr=block_ids_to_physics_ids.begin();itr!=block_ids_to_physics_ids.end();itr++)
      block_ids_to_cell_topo[itr->first] = mesh->getCellTopology(itr->first);
    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > aux_block_ids_to_cell_topo;
    for(auto itr=block_ids_to_aux_physics_ids.begin();itr!=block_ids_to_aux_physics_ids.end();itr++)
      aux_block_ids_to_cell_topo[itr->first] = mesh->getCellTopology(itr->first);

    // define physics block parameter list and boundary conditions
    std::vector<panzer::BC> bcs;
    panzer::buildBCs(bcs,bcs_pl,globalData);

    std::vector<panzer::BC> aux_bcs;
    panzer::buildBCs(aux_bcs,aux_bcs_pl,globalData);

    int workset_size = assembly_pl.get("Workset Size",20);

    // build the physics blocks objects
    std::vector<RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      Teuchos::TimeMonitor tMphysics(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: build physics blocks")));

      bool build_transient_support = true;
      // Can be overridden by the equation set
      int default_integration_order = 2;
      std::vector<std::string> tangentParamNames;
      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                 block_ids_to_cell_topo,
                                 physicsBlock_pl,
                                 default_integration_order,
                                 workset_size,
                                 eqset_factory,
                                 globalData,
                                 build_transient_support,
                                 physicsBlocks,
                                 tangentParamNames);
    }

    // build the auxiliary physics blocks objects
    std::vector<RCP<panzer::PhysicsBlock> > auxPhysicsBlocks;
    {
      Teuchos::TimeMonitor tMaux_physics(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: build auxiliary physics blocks")));

      bool build_transient_support = false;
      // Can be overridden by the equation set
      int default_integration_order = 2;
      std::vector<std::string> tangentParamNames;
      panzer::buildPhysicsBlocks(block_ids_to_aux_physics_ids,
                                 aux_block_ids_to_cell_topo,
                                 physicsBlock_pl,
                                 default_integration_order,
                                 workset_size,
                                 eqset_factory,
                                 globalData,
                                 build_transient_support,
                                 auxPhysicsBlocks,
                                 tangentParamNames);
    }

    // Add fields to the mesh data base (this is a peculiarity of how STK classic requires the
    // fields to be setup)
    createExodusFile(physicsBlocks, mesh_factory, mesh, exodus_output);

    // build worksets
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
      = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
      = Teuchos::rcp(new panzer::WorksetContainer);
    Teuchos::RCP<panzer::WorksetContainer> auxWkstContainer  // attach it to a workset container (uses lazy evaluation)
      = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    auxWkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physicsBlocks.size();i++) {
        wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
        auxWkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),auxPhysicsBlocks[i]->getWorksetNeeds());
    }
    wkstContainer->setWorksetSize(workset_size);
    auxWkstContainer->setWorksetSize(workset_size);

    // build DOF Managers and linear object factories

    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager<LocalOrdinal,GlobalOrdinal> >
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager<GlobalOrdinal>(mesh));

    // blocked degree of freedom manager
    panzer::BlockedDOFManagerFactory<LocalOrdinal,GlobalOrdinal> globalIndexerFactory;
    std::string fieldOrder = assembly_pl.get<std::string>("Field Order");
    RCP<panzer::UniqueGlobalIndexerBase > dofManager = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager,fieldOrder);

    // auxiliary dof manager
    std::string auxFieldOrder = assembly_pl.get<std::string>("Auxiliary Field Order");
    RCP<panzer::UniqueGlobalIndexerBase > auxDofManager = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),auxPhysicsBlocks,conn_manager,auxFieldOrder);

    // construct some linear algebra objects, build object to pass to evaluators
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory = Teuchos::rcp(new blockedLinObjFactory(comm,rcp_dynamic_cast<panzer::BlockedDOFManager<LocalOrdinal,GlobalOrdinal> >(dofManager,true)));
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > auxLinObjFactory = Teuchos::rcp(new blockedLinObjFactory(comm,rcp_dynamic_cast<panzer::BlockedDOFManager<LocalOrdinal,GlobalOrdinal> >(auxDofManager,true)));

    // Assign the dof managers to worksets
    wkstContainer->setGlobalIndexer(dofManager);
    auxWkstContainer->setGlobalIndexer(auxDofManager);

    // setup closure model
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    mini_em::ClosureModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    // add full maxwell solver to teko
    RCP<Teko::Cloneable> clone = rcp(new Teko::AutoClone<mini_em::FullMaxwellPreconditionerFactory>());
    Teko::PreconditionerFactory::addPreconditionerFactory("Full Maxwell Preconditioner",clone);

    // add callbacks to request handler. these are for requesting auxiliary operators and for providing
    // coordinate information to MueLu
    Teuchos::RCP<Teko::RequestHandler> req_handler = Teuchos::rcp(new Teko::RequestHandler());
    req_handler->addRequestCallback(Teuchos::rcp(new mini_em::OperatorRequestCallback(auxGlobalData, matrix_output)));
    Teuchos::RCP<panzer_stk::ParameterListCallbackBlocked<LocalOrdinal,GlobalOrdinal> > callback =
      rcp(new panzer_stk::ParameterListCallbackBlocked<LocalOrdinal,GlobalOrdinal>(
                                                                                   rcp_dynamic_cast<panzer_stk::STKConnManager<GlobalOrdinal> >(conn_manager,true),
                                                                                   rcp_dynamic_cast<panzer::BlockedDOFManager<LocalOrdinal,GlobalOrdinal> >(dofManager,true),
                                                                                   rcp_dynamic_cast<panzer::BlockedDOFManager<LocalOrdinal,GlobalOrdinal> >(auxDofManager,true)));
    req_handler->addRequestCallback(callback);

    // add discrete gradient
    {
      Teuchos::TimeMonitor tMdiscGrad(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: add discrete gradient")));
      addDiscreteGradientToRequestHandler(auxLinObjFactory,req_handler);
    }

    // build linear solver
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > lowsFactory
      = panzer_stk::buildLOWSFactory(true, dofManager, conn_manager,
                                     Teuchos::as<int>(mesh->getDimension()),
                                     comm, lin_solver_pl,req_handler, false, false, auxDofManager);

    //setup model evaluators
    RCP<panzer::ModelEvaluator<Scalar> > physics = rcp(new panzer::ModelEvaluator<Scalar> (linObjFactory, lowsFactory, globalData, true, 0.0));
    RCP<panzer::ModelEvaluator<Scalar> > auxPhysics = rcp(new panzer::ModelEvaluator<Scalar> (auxLinObjFactory, lowsFactory, globalData, false, 0.0));

    // add a volume response functionals
    std::map<int,std::string> responseIndexToName;
    for(Teuchos::ParameterList::ConstIterator itr=responses.begin();itr!=responses.end();++itr) {
      const std::string name = responses.name(itr);
      TEUCHOS_ASSERT(responses.entry(itr).isList());
      Teuchos::ParameterList & lst = Teuchos::getValue<Teuchos::ParameterList>(responses.entry(itr));


      // parameterize the builder
      panzer::FunctionalResponse_Builder<int,int> builder;
      builder.comm = MPI_COMM_WORLD; // good enough
      builder.cubatureDegree = 2;
      builder.requiresCellIntegral = lst.isType<bool>("Requires Cell Integral") ? lst.get<bool>("Requires Cell Integral"): false;
      builder.quadPointField = lst.get<std::string>("Field Name");

      // add the respone
      std::vector<std::string> eblocks;
      panzer::StringTokenizer(eblocks,lst.get<std::string>("Element Blocks"),",",true);

      std::vector<panzer::WorksetDescriptor> wkst_descs;
      for(std::size_t i=0;i<eblocks.size();i++)
        wkst_descs.push_back(panzer::blockDescriptor(eblocks[i]));

      int respIndex = physics->addResponse(name,wkst_descs,builder);
      responseIndexToName[respIndex] = name;
    }


    {
      Teuchos::TimeMonitor tMphysicsEval(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: setup physics model evaluator")));
      physics->setupModel(wkstContainer,physicsBlocks,bcs,
                          *eqset_factory,
                          bc_factory,
                          cm_factory,
                          cm_factory,
                          closure_models,
                          user_data,false,"");

      // add auxiliary data to model evaluator
      for(panzer::GlobalEvaluationDataContainer::const_iterator itr=auxGlobalData->begin();itr!=auxGlobalData->end();++itr)
        physics->addNonParameterGlobalEvaluationData(itr->first,itr->second);
    }

    {
      Teuchos::TimeMonitor tMauxphysicsEval(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: setup auxiliary physics model evaluator")));

      auxPhysics->setupModel(auxWkstContainer,auxPhysicsBlocks,aux_bcs,
                             *eqset_factory,
                             bc_factory,
                             cm_factory,
                             cm_factory,
                             closure_models,
                             user_data,false,"");

      // evaluate the auxiliary model to obtain auxiliary operators
      for(panzer::GlobalEvaluationDataContainer::const_iterator itr=auxGlobalData->begin();itr!=auxGlobalData->end();++itr)
        auxPhysics->addNonParameterGlobalEvaluationData(itr->first,itr->second);
    }

    Thyra::ModelEvaluatorBase::InArgs<Scalar> auxInArgs = auxPhysics->getNominalValues();
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> auxOutArgs = auxPhysics->createOutArgs();
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > aux_W_op = auxPhysics->create_W_op();
    auxOutArgs.set_W_op(aux_W_op);
    auxPhysics->evalModel(auxInArgs, auxOutArgs);

    // setup a response library to write to the mesh
    RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
      = buildSTKIOResponseLibrary(physicsBlocks,linObjFactory,wkstContainer,dofManager,cm_factory,mesh,
                                  closure_models);



    // set up the solution vector, jacobian, and residual
    RCP<Thyra::VectorBase<Scalar> > solution_vec = Thyra::createMember(physics->get_x_space());
    Thyra::assign(solution_vec.ptr(),zero);
    RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian = physics->create_W();

    RCP<Thyra::VectorBase<Scalar> > residual = Thyra::createMember(physics->get_f_space());

    // set up the model evaluator
    Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = physics->createInArgs();
    inArgs.set_alpha(one/dt);
    inArgs.set_beta(one);

    // initial condition is zero, define x_dot accordingly
    RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
    RCP<Thyra::VectorBase<Scalar> > x_dot = Thyra::createMember(physics->get_x_space());
    Thyra::V_StV(x_dot.ptr(),one/dt,*x);
    inArgs.set_x_dot(x_dot);
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs = physics->createOutArgs();
    outArgs.set_f(residual);
    outArgs.set_W(jacobian);

    // compute the jacobian matrix only once
    Kokkos::fence();
    physics->evalModel(inArgs,outArgs);
    outArgs.set_W(RCP<Thyra::LinearOpWithSolveBase<Scalar> >(NULL));

    // take time-steps with Backward Euler
    if (exodus_output)
      writeToExodus<Scalar>(0,solution_vec,*physics,*stkIOResponseLibrary,*mesh);

    RCP<Thyra::VectorBase<Scalar> > correction_vec = Thyra::createMember(physics->get_x_space());
    Thyra::assign(correction_vec.ptr(),zero);

    {
      Teuchos::TimeMonitor tMts(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: timestepper")));
      for(int ts = 1; ts < numTimeSteps+1; ts++)
        {
          (*out) << std::endl;
          (*out) << "**************************************************" << std::endl;
          (*out) << "* starting time step " << ts << std::endl;

          RCP<Thyra::VectorBase<Scalar> > x_old = solution_vec->clone_v();

          inArgs.set_t(dt*ts);

          // start Newton loop (nonlinear case)
          // for() until convergence

          Thyra::V_StVpStV(x_dot.ptr(),one/dt,*solution_vec,-one/dt,*x_old);
          inArgs.set_x(solution_vec);
          inArgs.set_x_dot(x_dot);

          // construct the residual
          physics->evalModel(inArgs,outArgs);

          // solve
          if (doSolveTimings)
            for (int rep = 0; rep < numReps; rep++) {
              Thyra::assign(correction_vec.ptr(),zero);
              jacobian->solve(Thyra::NOTRANS,*residual,correction_vec.ptr());
            }
          else
            jacobian->solve(Thyra::NOTRANS,*residual,correction_vec.ptr());
          Thyra::V_StVpStV(solution_vec.ptr(),one,*solution_vec,-one,*correction_vec);

          // end for()
          // end Newton loop (nonlinear case)
          
          // print out responses
          {
            Thyra::ModelEvaluatorBase::InArgs<Scalar> respInArgs = physics->createInArgs();
            Thyra::ModelEvaluatorBase::OutArgs<Scalar> respOutArgs = physics->createOutArgs();

            respInArgs.set_x(solution_vec);

            for(int i=0;i<respOutArgs.Ng();i++) {
              Teuchos::RCP<Thyra::VectorBase<Scalar> > response = Thyra::createMember(*physics->get_g_space(i));
              respOutArgs.set_g(i,response);
            }

            physics->evalModel(respInArgs, respOutArgs);

            for (auto elem: responseIndexToName) {
              Teuchos::RCP<Thyra::VectorBase<Scalar> > g = respOutArgs.get_g(elem.first);
              *out << elem.second << " = " << Thyra::get_ele(*g,0) << std::endl;
            }
          }
          // for checking correctness:
          //    the current pulse is designed to peak at 10 time-steps
          //    consequently, EM energy increases quickly up to 10 time-steps and then growth slows
          //    asymptoting around 20 time-steps when the pulse is negligible.
          //    since the pulse depends on the time-step, so does the energy.
          //    energy normalized by 1/dt^2 is approximately the same across different values of dt.

          // write to an exodus file
          if (exodus_output)
            {
              Teuchos::TimeMonitor tMexodus(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: timestepper: writeToExodus")));
              writeToExodus<Scalar>(dt*ts,solution_vec,*physics,*stkIOResponseLibrary,*mesh);
            }

          (*out) << std::endl;
          (*out) << "* finished time step " << ts << std::endl;
          (*out) << "**************************************************" << std::endl << std::endl;
        }
    }
  }

  // Output timer data
  if (use_stacked_timer) {
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stacked_timer->report(*out, comm, options);
  } else
    Teuchos::TimeMonitor::summarize(*out,false,true,false,Teuchos::Union,"",true);

  return EXIT_SUCCESS;
}

int main(int argc,char * argv[]){
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  Kokkos::initialize(argc, argv);

  Teuchos::CommandLineProcessor clp(false);
  linearAlgebraType linAlgebraValues[2] = {linAlgTpetra, linAlgEpetra};
  const char * linAlgebraNames[2] = {"Tpetra", "Epetra"};
  linearAlgebraType linAlgebra = linAlgTpetra;
  clp.setOption<linearAlgebraType>("linAlgebra",&linAlgebra,2,linAlgebraValues,linAlgebraNames);
  solverType solverValues[4] = {AUGMENTATION, MUELU_REFMAXWELL, ML_REFMAXWELL, CG};
  const char * solverNames[4] = {"Augmentation", "MueLu-RefMaxwell", "ML-RefMaxwell", "CG"};
  solverType solver = MUELU_REFMAXWELL;
  clp.setOption<solverType>("solver",&solver,4,solverValues,solverNames,"Solver that is used");
  // bool useComplex = false;
  // clp.setOption("complex","real",&useComplex);
  clp.recogniseAllOptions(false);
  switch (clp.parse(argc, argv, NULL)) {
    case Teuchos::CommandLineProcessor::PARSE_ERROR:                return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:         break;
  }

  if (solver == ML_REFMAXWELL) {
    TEUCHOS_ASSERT(linAlgebra == linAlgEpetra);
    // TEUCHOS_ASSERT(!useComplex);
  }

  int retVal;
  if (linAlgebra == linAlgTpetra) {
    // if (useComplex) {
// #if defined(HAVE_TPETRA_COMPLEX_DOUBLE)
//       typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,std::complex<double>,int,panzer::Ordinal64> blockedLinObjFactory;
//       retVal = main_<std::complex<double>,int,panzer::Ordinal64,blockedLinObjFactory>(clp, argc, argv);
// #else
//       std::cout << std::endl
//                 << "WARNING" << std::endl
//                 << "Tpetra was compiled without Scalar=std::complex<double>." << std::endl << std::endl;
//       return EXIT_FAILURE;
// #endif
//     } else {
      typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64> blockedLinObjFactory;
      retVal = main_<double,int,panzer::Ordinal64,blockedLinObjFactory>(clp, argc, argv);
//    }
  } else if (linAlgebra == linAlgEpetra) {
    // TEUCHOS_ASSERT(!useComplex);
    typedef typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> blockedLinObjFactory;
    retVal = main_<double,int,int,blockedLinObjFactory>(clp, argc, argv);
  } else
    TEUCHOS_ASSERT(false);

  Kokkos::finalize();

  return retVal;
}

void createExodusFile(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                      Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory,
                      Teuchos::RCP<panzer_stk::STK_Interface> mesh,
                      const bool & exodus_out) {
  for(std::size_t i=0;i<physicsBlocks.size();i++) {
    Teuchos::RCP<panzer::PhysicsBlock> pb = physicsBlocks[i]; // we are assuming only one physics block

    const std::vector<panzer::StrPureBasisPair> & blockFields = pb->getProvidedDOFs();

    // insert all fields into a set
    std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp> fieldNames;
    fieldNames.insert(blockFields.begin(),blockFields.end());

    // build string for modifiying vectors
    std::vector<std::string> dimenStr(3);
    dimenStr[0] = "X"; dimenStr[1] = "Y"; dimenStr[2] = "Z";

    // add basis to DOF manager: block specific
    std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp>::const_iterator fieldItr;
    for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {
      Teuchos::RCP<const panzer::PureBasis> basis = fieldItr->second;
      if(basis->getElementSpace()==panzer::PureBasis::HGRAD)
        mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
      else if(basis->getElementSpace()==panzer::PureBasis::CONST )
        mesh->addCellField(fieldItr->first,pb->elementBlockID());
      else if(basis->getElementSpace()==panzer::PureBasis::HCURL ||
          basis->getElementSpace()==panzer::PureBasis::HDIV    ) {
        for(int i=0;i<basis->dimension();i++)
          mesh->addCellField(fieldItr->first+dimenStr[i],pb->elementBlockID());
      }
    }

    std::vector<std::string> block_names;
    mesh->getElementBlockNames(block_names);

    Teuchos::ParameterList output_pl("Output");
    output_pl.sublist("Cell Average Quantities");
    Teuchos::ParameterList& cell_avg_v = output_pl.sublist("Cell Average Vectors");
    cell_avg_v.set(block_names[0],"CURRENT");
    output_pl.sublist("Cell Quantities");
    output_pl.sublist("Nodal Quantities");
    output_pl.sublist("Allocate Nodal Quantities");
    mini_em::addFieldsToMesh(*mesh,output_pl);
  }
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

  if (exodus_out)
    mesh->setupExodusFile("mesh_output.exo");
}

Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >
buildSTKIOResponseLibrary(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
    const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linObjFactory,
    const Teuchos::RCP<panzer::WorksetContainer> & wkstContainer,
    const Teuchos::RCP<panzer::UniqueGlobalIndexerBase> & globalIndexer,
    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
    const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
    const Teuchos::ParameterList & closure_model_pl) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
  = rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,globalIndexer,linObjFactory));

  // get a vector of all the element blocks
  std::vector<std::string> eBlocks;
  mesh->getElementBlockNames(eBlocks);

  panzer_stk::RespFactorySolnWriter_Builder builder;
  builder.mesh = mesh;

  stkIOResponseLibrary->addResponse("Main Field Output",eBlocks,builder);

  std::vector<std::string> block_names;
  mesh->getElementBlockNames(block_names);

  // this automatically adds in the nodal fields
  Teuchos::ParameterList output_pl("Output");
  output_pl.sublist("Cell Average Quantities");
  Teuchos::ParameterList& cell_avg_v = output_pl.sublist("Cell Average Vectors");
  cell_avg_v.set(block_names[0],"CURRENT");
  output_pl.sublist("Cell Quantities");
  output_pl.sublist("Nodal Quantities");
  output_pl.sublist("Allocate Nodal Quantities");
  panzer_stk::IOClosureModelFactory_TemplateBuilder<panzer::Traits> io_cm_builder(cm_factory,mesh,
                                                                                  output_pl);
  panzer::ClosureModelFactory_TemplateManager<panzer::Traits> io_cm_factory;
  io_cm_factory.buildObjects(io_cm_builder);

  stkIOResponseLibrary->buildResponseEvaluators(physicsBlocks,
      io_cm_factory,
      closure_model_pl,
      Teuchos::ParameterList());

  return stkIOResponseLibrary;
}

template <class Scalar>
void writeToExodus(double time_stamp,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & x,
    const panzer::ModelEvaluator<Scalar> & model,
    panzer::ResponseLibrary<panzer::Traits> & stkIOResponseLibrary,
    panzer_stk::STK_Interface & mesh)
{
  // fill STK mesh objects
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = model.createInArgs();
  inArgs.set_x(x);
  inArgs.set_t(time_stamp);

  panzer::AssemblyEngineInArgs respInput;
  model.setupAssemblyInArgs(inArgs,respInput);

  stkIOResponseLibrary.addResponsesToInArgs<panzer::Traits::Residual>(respInput);
  stkIOResponseLibrary.evaluate<panzer::Traits::Residual>(respInput);

  mesh.writeToExodus(time_stamp);
}
