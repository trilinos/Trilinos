// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#include "Kokkos_View_Fad.hpp"
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include "Panzer_NodeType.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
# include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#endif
#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"
// #include "Panzer_InitialCondition_Builder.hpp"
// #include "Panzer_CheckBCConsistency.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"

#include "Panzer_STK_MeshFactory.hpp"
#include "Panzer_STK_SetupLOWSFactory.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_ParameterListCallbackBlocked.hpp"

// includes specific to this example
#include "PanzerMiniEM_config.hpp"
#include "MiniEM_EquationSetFactory.hpp"
#include "MiniEM_BCStrategy_Factory.hpp"
#include "MiniEM_ClosureModel_Factory_TemplateBuilder.hpp"
#include "MiniEM_OperatorRequestCallback.hpp"
#include "MiniEM_FullMaxwellPreconditionerFactory.hpp"
#include "MiniEM_HigherOrderMaxwellPreconditionerFactory.hpp"
#include "MiniEM_FullMaxwellPreconditionerFactory_Augmentation.hpp"
#include "MiniEM_FullDarcyPreconditionerFactory.hpp"
#include "MiniEM_Interpolation.hpp"
#include "MiniEM_helpers.hpp"

#include <string>
#include <iostream>


template <class Scalar>
void writeToExodus(double time_stamp,
                   const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & x,
                   const panzer::ModelEvaluator<Scalar> & model,
                   panzer::ResponseLibrary<panzer::Traits> & stkIOResponseLibrary,
                   panzer_stk::STK_Interface & mesh);

/********************************************************************************
 * This driver sets up either
 * - first order Maxwell equations with edge-face discretization for E, B,
 * - mixed form Darcy flow.
 * We use backward Euler time-stepping with fixed CFL.
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

using namespace mini_em;

using mini_em::physicsType, mini_em::MAXWELL, mini_em::DARCY;
using mini_em::solverType, mini_em::AUGMENTATION, mini_em::MUELU, mini_em::ML, mini_em::CG, mini_em::GMRES;
using mini_em::linearAlgebraType, mini_em::linAlgTpetra, mini_em::linAlgEpetra;


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class blockedLinObjFactory, bool useTpetra>
int main_(Teuchos::CommandLineProcessor &clp, int argc,char * argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
  Teuchos::RCP<const Teuchos::MpiComm<int> > comm
    = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());
  if (comm->getSize() > 1) {
    out->setOutputToRootOnly(0);
  }

  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
  bool use_stacked_timer;
  std::string test_name = "MiniEM 3D RefMaxwell";

  // Figure of merit data for acceptance testing
  bool print_fom;
  size_t fom_num_cells;

  {
    // defaults for command-line options
    int x_elements=-1,y_elements=-1,z_elements=-1,basis_order=1;
    int workset_size=2000;
    std::string pCoarsenScheduleStr = "1";
    double dt=0.0;
    std::string meshFile = "";
    bool exodus_output = false;
    bool matrix_output = false;
    std::string input_file = "maxwell.xml";
    std::string xml = "";
    solverType solverValues[5] = {AUGMENTATION, MUELU, ML, CG, GMRES};
    const char * solverNames[5] = {"Augmentation", "MueLu", "ML", "CG", "GMRES"};
    bool preferTPLs = false;
    bool truncateMueLuHierarchy = false;
    solverType solver = MUELU;
    int numTimeSteps = 1;
    double finalTime = -1.;
    bool resetSolver = false;
    bool doSolveTimings = false;
    bool matrixFree = false;
    int numReps = 0;
    linearAlgebraType linAlgebraValues[2] = {linAlgTpetra, linAlgEpetra};
    const char * linAlgebraNames[2] = {"Tpetra", "Epetra"};
    linearAlgebraType linAlgebra = linAlgTpetra;
    clp.setOption<linearAlgebraType>("linAlgebra",&linAlgebra,2,linAlgebraValues,linAlgebraNames);
    use_stacked_timer = true;
    print_fom = true;
    clp.setOption("x-elements",&x_elements);
    clp.setOption("y-elements",&y_elements);
    clp.setOption("z-elements",&z_elements);
    clp.setOption("basis-order",&basis_order);
    clp.setOption("workset-size",&workset_size);
    clp.setOption("pCoarsenSchedule",&pCoarsenScheduleStr);
    clp.setOption("dt",&dt,"Override \"dt\" specified in input file.");
    clp.setOption("meshFile",&meshFile,"Override input mesh file specified in input file.");
    clp.setOption("exodus-output","no-exodus-output",&exodus_output);
    clp.setOption("matrix-output","no-matrix-output",&matrix_output);
    clp.setOption("inputFile",&input_file,"XML file with the problem definitions");
    clp.setOption("solverFile",&xml,"XML file with the solver params");
    clp.setOption<solverType>("solver",&solver,5,solverValues,solverNames,"Solver that is used");
    clp.setOption("tpl", "no-tpl", &preferTPLs, "Prefer TPL usage over fused kernels");
    clp.setOption("truncateMueLuHierarchy", "no-truncateMueLuHierarchy", &truncateMueLuHierarchy, "Truncate the MueLu hierarchy");
    clp.setOption("numTimeSteps",&numTimeSteps);
    clp.setOption("finalTime",&finalTime);
    clp.setOption("matrixFree","no-matrixFree",&matrixFree,"matrix-free operators");
    clp.setOption("resetSolver","no-resetSolver",&resetSolver,"update the solver in every timestep");
    clp.setOption("doSolveTimings","no-doSolveTimings",&doSolveTimings,"repeat the first solve \"numTimeSteps\" times");
    clp.setOption("stacked-timer","no-stacked-timer",&use_stacked_timer,"Run with or without stacked timer output");
    clp.setOption("test-name", &test_name, "Name of test (for Watchr output)");
    clp.setOption("print-fom","no-print-fom",&print_fom,"print the figure of merit for acceptance testing");
#ifdef HAVE_TEUCHOS_STACKTRACE
    bool stacktrace = false;
    clp.setOption("stacktrace", "nostacktrace", &stacktrace, "display stacktrace");
#endif

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

#ifdef HAVE_TEUCHOS_STACKTRACE
    if (stacktrace)
      Teuchos::print_stack_on_segfault();
#endif


    if (use_stacked_timer) {
      stacked_timer = rcp(new Teuchos::StackedTimer("Mini-EM"));
      Teuchos::RCP<Teuchos::FancyOStream> verbose_out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
      verbose_out->setShowProcRank(true);
      stacked_timer->setVerboseOstream(verbose_out);
    }
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
    // RCP<Teuchos::ParameterList> physicsBlock_pl      = rcp(new Teuchos::ParameterList(input_params->sublist("Physics Blocks")));
    Teuchos::ParameterList &physicsBlock_pl          = input_params->sublist("Physics Blocks");
    Teuchos::ParameterList & assembly_pl             = input_params->sublist("Assembly");
    Teuchos::ParameterList & block_to_physics_pl     = input_params->sublist("Block ID to Physics ID Mapping");
    Teuchos::ParameterList & block_to_aux_physics_pl = input_params->sublist("Block ID to Auxiliary Physics ID Mapping");
    Teuchos::ParameterList & ops_pl                  = input_params->sublist("Operators");
    Teuchos::ParameterList & aux_ops_pl              = input_params->sublist("Auxiliary Operators");
    Teuchos::ParameterList & bcs_pl                  = input_params->sublist("Boundary Conditions");
    Teuchos::ParameterList & aux_bcs_pl              = input_params->sublist("Auxiliary Boundary Conditions");
    Teuchos::ParameterList & closure_models          = input_params->sublist("Closure Models");
    Teuchos::ParameterList responses                 = input_params->sublist("Responses");
    Teuchos::ParameterList & user_data               = input_params->sublist("User Data");

    // Set basis order
    Teuchos::ParameterList& physicsEqSet = physicsBlock_pl.sublist("Maxwell Physics").sublist("Maxwell Physics");
    const std::string physicsTypeStr = physicsEqSet.get<std::string>("Type");
    physicsType physics;
    if (physicsTypeStr == "Maxwell")
      physics = MAXWELL;
    else if (physicsTypeStr == "Darcy")
      physics = DARCY;
    else
      TEUCHOS_ASSERT(false);
    basis_order = physicsEqSet.get("Basis Order", basis_order);
    physicsEqSet.set("Integration Order", 2*basis_order);

    RCP<panzer_stk::STK_Interface> mesh;
    int dim;
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
    {
      Teuchos::TimeMonitor tMmesh(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: build mesh")));
      double mesh_size = 0.;
      mini_em::getMesh(mesh_pl, meshFile, x_elements, y_elements, z_elements, basis_order, comm, mesh, mesh_factory, mesh_size);
      dim = Teuchos::as<int>(mesh->getDimension());

      // set dt
      if (dt <= 0.) {
        if (mesh_pl.get<std::string>("Source") == "Exodus File" || meshFile != "") {
          RCP<Teuchos::ParameterList> input_pl = rcp(new Teuchos::ParameterList(mesh_pl.sublist("Exodus File")));
          if (input_pl->isType<double>("dt"))
            dt = input_pl->get<double>("dt");
          if (input_pl->isType<int>("num time steps"))
            numTimeSteps = input_pl->get<int>("num time steps");
          if (input_pl->isType<double>("final time"))
            finalTime = input_pl->get<double>("final time");
        } else if (mesh_pl.get<std::string>("Source") ==  "Pamgen Mesh") {
          Teuchos::ParameterList & pamgen_pl = mesh_pl.sublist("Pamgen Mesh");
          dt = pamgen_pl.get<double>("dt");
        } else if (mesh_pl.get<std::string>("Source") == "Inline Mesh") {
          Teuchos::ParameterList & inline_gen_pl = mesh_pl.sublist("Inline Mesh");
          if (inline_gen_pl.isType<double>("final time"))
            finalTime = inline_gen_pl.get<double>("final time");
          if (inline_gen_pl.isType<double>("dt"))
            dt = inline_gen_pl.get<double>("dt");
          else {
            if (physics == MAXWELL) {
              // This is only correct when epsilon==epsilon0, mu==mu0
              double cfl = inline_gen_pl.get<double>("CFL");
              double c = 299792458.0;
              dt = cfl * mesh_size / basis_order / c;
            } else if (physics == DARCY) {
              // dt = mesh_size * mesh_size / basis_order / basis_order;
              dt = mesh_size / basis_order;
            }
          }
        }
      }

      if (finalTime <= 0.)
        finalTime = numTimeSteps*dt;
      else {
        numTimeSteps = std::max(Teuchos::as<int>(std::ceil(finalTime/dt)), 1);
        dt = finalTime/numTimeSteps;
      }

      (*out) << "dt = " << dt << std::endl << std::endl;
      if (dt <= 0.0)
        throw;
    }

    RCP<Teuchos::ParameterList> lin_solver_pl = mini_em::getSolverParameters(linAlgebra, physics, solver, dim, comm, out, xml, basis_order, preferTPLs, truncateMueLuHierarchy);

    if (lin_solver_pl->sublist("Preconditioner Types").isSublist("Teko") &&
        lin_solver_pl->sublist("Preconditioner Types").sublist("Teko").isSublist("Inverse Factory Library")) {
      if (physics == MAXWELL)
        lin_solver_pl->sublist("Preconditioner Types").sublist("Teko").set("Inverse Type", "Maxwell");
      else if (physics == DARCY)
        lin_solver_pl->sublist("Preconditioner Types").sublist("Teko").set("Inverse Type", "Darcy");
      if (lin_solver_pl->sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").isSublist("Maxwell"))
        lin_solver_pl->sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("Maxwell").set("dt",dt);
      if (lin_solver_pl->sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").isSublist("Darcy"))
        lin_solver_pl->sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("Darcy").set("dt",dt);
    }
    lin_solver_pl->print(*out);

    std::string modelID = physicsEqSet.get<std::string>("Model ID");
    std::string auxModelID;
    setClosureParameters(physics, physicsEqSet, closure_models, dt, auxModelID);

    std::map<std::string,std::string> block_ids_to_physics_ids, block_ids_to_aux_physics_ids;
    panzer::buildBlockIdToPhysicsIdMap(block_ids_to_physics_ids, block_to_physics_pl);
    panzer::buildBlockIdToPhysicsIdMap(block_ids_to_aux_physics_ids, block_to_aux_physics_pl);

    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
    for(auto itr=block_ids_to_physics_ids.begin();itr!=block_ids_to_physics_ids.end();itr++)
      block_ids_to_cell_topo[itr->first] = mesh->getCellTopology(itr->first);
    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > aux_block_ids_to_cell_topo;
    for(auto itr=block_ids_to_aux_physics_ids.begin();itr!=block_ids_to_aux_physics_ids.end();itr++)
      aux_block_ids_to_cell_topo[itr->first] = mesh->getCellTopology(itr->first);

    std::string auxFieldOrder;
    setAuxiliaryOperatorParameters(physics, solver, basis_order, pCoarsenScheduleStr, matrixFree, *input_params, *lin_solver_pl, auxFieldOrder);

    // define physics block parameter list and boundary conditions
    std::vector<panzer::BC> bcs;
    panzer::buildBCs(bcs,bcs_pl,globalData);

    std::vector<panzer::BC> aux_bcs;
    panzer::buildBCs(aux_bcs,aux_bcs_pl,globalData);

    workset_size = assembly_pl.get("Workset Size", workset_size);

    // build the physics blocks objects
    std::vector<RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      Teuchos::TimeMonitor tMphysics(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: build physics blocks")));

      bool build_transient_support = true;
      // Can be overridden by the equation set
      int default_integration_order = 2*basis_order;
      std::vector<std::string> tangentParamNames;
      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                 block_ids_to_cell_topo,
                                 rcpFromRef(physicsBlock_pl),
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
      int default_integration_order = 2*basis_order;
      std::vector<std::string> tangentParamNames;
      panzer::buildPhysicsBlocks(block_ids_to_aux_physics_ids,
                                 aux_block_ids_to_cell_topo,
                                 rcpFromRef(physicsBlock_pl),
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
    createExodusFile(physicsBlocks, mesh_factory, mesh, exodus_output, comm, physics);

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
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    // blocked degree of freedom manager
    panzer::BlockedDOFManagerFactory globalIndexerFactory;
    std::string fieldOrder = assembly_pl.get<std::string>("Field Order");
    RCP<panzer::GlobalIndexer > dofManager = globalIndexerFactory.buildGlobalIndexer(comm->getRawMpiComm(),physicsBlocks,conn_manager,fieldOrder);

    // auxiliary dof manager
    RCP<panzer::GlobalIndexer > auxDofManager = globalIndexerFactory.buildGlobalIndexer(comm->getRawMpiComm(),auxPhysicsBlocks,conn_manager,auxFieldOrder);

    // construct some linear algebra objects, build object to pass to evaluators
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory = Teuchos::rcp(new blockedLinObjFactory(comm,rcp_dynamic_cast<panzer::BlockedDOFManager>(dofManager,true)));
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > auxLinObjFactory = Teuchos::rcp(new blockedLinObjFactory(comm,rcp_dynamic_cast<panzer::BlockedDOFManager>(auxDofManager,true)));

    // Assign the dof managers to worksets
    wkstContainer->setGlobalIndexer(dofManager);
    auxWkstContainer->setGlobalIndexer(auxDofManager);

    // setup closure model
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    mini_em::ClosureModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    // add full maxwell preconditioner to teko
    RCP<Teko::Cloneable> clone = rcp(new Teko::AutoClone<mini_em::FullMaxwellPreconditionerFactory>());
    Teko::PreconditionerFactory::addPreconditionerFactory("Full Maxwell Preconditioner",clone);

    // add full darcy preconditioner to teko
    RCP<Teko::Cloneable> cloneDarcy = rcp(new Teko::AutoClone<mini_em::FullDarcyPreconditionerFactory>());
    Teko::PreconditionerFactory::addPreconditionerFactory("Full Darcy Preconditioner",cloneDarcy);

    // add higher-order maxwell preconditioner to teko
    RCP<Teko::Cloneable> cloneHO = rcp(new Teko::AutoClone<mini_em::HigherOrderMaxwellPreconditionerFactory>());
    Teko::PreconditionerFactory::addPreconditionerFactory("Higher Order Maxwell Preconditioner",cloneHO);

    // add augmentation preconditioner to teko
    RCP<Teko::Cloneable> cloneAug = rcp(new Teko::AutoClone<mini_em::FullMaxwellPreconditionerFactory_Augmentation>());
    Teko::PreconditionerFactory::addPreconditionerFactory("Full Maxwell Preconditioner: Augmentation",cloneAug);

    // add callbacks to request handler. these are for requesting auxiliary operators and for providing
    // coordinate information to MueLu
    Teuchos::RCP<Teko::RequestHandler> req_handler = Teuchos::rcp(new Teko::RequestHandler());
    req_handler->addRequestCallback(Teuchos::rcp(new mini_em::OperatorRequestCallback(auxGlobalData, matrix_output)));
    Teuchos::RCP<panzer_stk::ParameterListCallbackBlocked> callback =
      rcp(new panzer_stk::ParameterListCallbackBlocked(rcp_dynamic_cast<panzer_stk::STKConnManager>(conn_manager,true),
                                                       rcp_dynamic_cast<panzer::BlockedDOFManager>(dofManager,true),
                                                       rcp_dynamic_cast<panzer::BlockedDOFManager>(auxDofManager,true)));
    req_handler->addRequestCallback(callback);

    {
      if (physics == MAXWELL) {
        // add discrete curl
        ops_pl.sublist("Discrete Curl").set("Source", "E_edge");
        ops_pl.sublist("Discrete Curl").set("Target", "B_face");
        ops_pl.sublist("Discrete Curl").set("Op", "curl");
        ops_pl.sublist("Discrete Curl").set("matrix-free", matrixFree);
      }

      if (physics == DARCY) {
        // add discrete div
        ops_pl.sublist("Discrete Div").set("Source", "sigma");
        ops_pl.sublist("Discrete Div").set("Target", "u");
        ops_pl.sublist("Discrete Div").set("Op", "div");
        ops_pl.sublist("Discrete Div").set("matrix-free", matrixFree);
      }

      // add request handlers for all interpolation type operators
      // (discrete grad, curl, div and interpolations between spaces of different orders)
      std::vector<std::pair<Teuchos::ParameterList,
                            Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > > > opLists = {{ops_pl, linObjFactory},
                                                                                                   {aux_ops_pl, auxLinObjFactory}};
      for (auto p = opLists.begin(); p != opLists.end(); ++p) {
        for (auto it = p->first.begin(); it != p->first.end(); ++it) {
          std::string name = it->first;
          Teuchos::ParameterList pl = it->second.getValue<Teuchos::ParameterList>(0);
          const std::string src = pl.get<std::string>("Source");
          const std::string tgt = pl.get<std::string>("Target");
          const std::string op = pl.get<std::string>("Op","value");
          const bool waitForRequest = pl.get<bool>("waitForRequest",true);
          const bool dump = pl.get<bool>("dump",false);
          const bool useMatrixFree = pl.get<bool>("matrix-free",matrixFree);
          Intrepid2::EOperator eOp;
          if (op == "value")
            eOp = Intrepid2::OPERATOR_VALUE;
          else if (op == "grad")
            eOp = Intrepid2::OPERATOR_GRAD;
          else if (op == "curl")
            eOp = Intrepid2::OPERATOR_CURL;
          else if (op == "div")
            eOp = Intrepid2::OPERATOR_DIV;
          else
            TEUCHOS_ASSERT(false);
          addInterpolationToRequestHandler(name, p->second, req_handler, src, tgt, eOp, waitForRequest, dump, workset_size, useMatrixFree);
        }
      }
    }

    // build linear solver
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > lowsFactory
      = panzer_stk::buildLOWSFactory(true, dofManager, conn_manager,
                                     Teuchos::as<int>(mesh->getDimension()),
                                     comm, lin_solver_pl,req_handler, false, false, auxDofManager);

    //setup model evaluators
    RCP<panzer::ModelEvaluator<Scalar> > physicsME = rcp(new panzer::ModelEvaluator<Scalar> (linObjFactory, lowsFactory, globalData, true, 0.0));
    RCP<panzer::ModelEvaluator<Scalar> > auxPhysicsME = rcp(new panzer::ModelEvaluator<Scalar> (auxLinObjFactory, lowsFactory, globalData, false, 0.0));
    physicsME->template disableEvaluationType<panzer::Traits::Tangent>();
    auxPhysicsME->template disableEvaluationType<panzer::Traits::Tangent>();

    // add a volume response functionals
    std::map<int,std::string> responseIndexToName;
    for(Teuchos::ParameterList::ConstIterator itr=responses.begin();itr!=responses.end();++itr) {
      const std::string name = responses.name(itr);
      TEUCHOS_ASSERT(responses.entry(itr).isList());
      Teuchos::ParameterList & lst = Teuchos::getValue<Teuchos::ParameterList>(responses.entry(itr));


      // parameterize the builder
      panzer::FunctionalResponse_Builder<int,int> builder;
      builder.comm = (*comm->getRawMpiComm())();
      builder.cubatureDegree = 2*basis_order;
      builder.requiresCellIntegral = lst.isType<bool>("Requires Cell Integral") ? lst.get<bool>("Requires Cell Integral"): false;
      builder.quadPointField = lst.get<std::string>("Field Name");

      // add the response
      std::vector<std::string> eblocks;
      panzer::StringTokenizer(eblocks,lst.get<std::string>("Element Blocks"),",",true);

      std::vector<panzer::WorksetDescriptor> wkst_descs;
      for(std::size_t i=0;i<eblocks.size();i++)
        wkst_descs.push_back(panzer::blockDescriptor(eblocks[i]));

      int respIndex = physicsME->addResponse(name,wkst_descs,builder);
      responseIndexToName[respIndex] = name;
    }


    {
      Teuchos::TimeMonitor tMphysicsEval(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: setup physics model evaluator")));
      physicsME->setupModel(wkstContainer,physicsBlocks,bcs,
                            *eqset_factory,
                            bc_factory,
                            cm_factory,
                            cm_factory,
                            closure_models,
                            user_data,false,"");

      // add auxiliary data to model evaluator
      for(panzer::GlobalEvaluationDataContainer::const_iterator itr=auxGlobalData->begin();itr!=auxGlobalData->end();++itr)
        physicsME->addNonParameterGlobalEvaluationData(itr->first,itr->second);
    }

    {
      Teuchos::TimeMonitor tMauxphysicsEval(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: setup auxiliary physics model evaluator")));

      auxPhysicsME->setupModel(auxWkstContainer,auxPhysicsBlocks,aux_bcs,
                               *eqset_factory,
                               bc_factory,
                               cm_factory,
                               cm_factory,
                               closure_models,
                               user_data,false,"");

      // evaluate the auxiliary model to obtain auxiliary operators
      for(panzer::GlobalEvaluationDataContainer::const_iterator itr=auxGlobalData->begin();itr!=auxGlobalData->end();++itr)
        auxPhysicsME->addNonParameterGlobalEvaluationData(itr->first,itr->second);
    }

    {
      Teuchos::TimeMonitor tMauxphysicsEval(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: eval auxiliary physics model evaluator")));
      Thyra::ModelEvaluatorBase::InArgs<Scalar> auxInArgs = auxPhysicsME->getNominalValues();
      Thyra::ModelEvaluatorBase::OutArgs<Scalar> auxOutArgs = auxPhysicsME->createOutArgs();
      Teuchos::RCP<Thyra::LinearOpBase<Scalar> > aux_W_op = auxPhysicsME->create_W_op();

      auxOutArgs.set_W_op(aux_W_op);
      auxPhysicsME->evalModel(auxInArgs, auxOutArgs);
    }

    // setup a response library to write to the mesh
    RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
      = buildSTKIOResponseLibrary(physicsBlocks,linObjFactory,wkstContainer,dofManager,cm_factory,mesh,
                                  closure_models, physics);

    // set up the model evaluator
    Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = physicsME->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs = physicsME->createOutArgs();

    // We are setting up a model evaluator for
    // f(x_dot, x, t) = M*x_dot + A*x - rhs
    //
    // Eval with
    // alpha = 1/dt, beta = 1
    // W = 1/dt*M + A
    //
    // Eval with
    // alpha = 1/dt, beta = 1
    // x = solution_vector, x_dot = 0
    // residual = A*solution_vector - rhs
    //
    // Solve
    // W*c = residual
    // Set
    // new_solution_vector = solution_vector - c
    //
    // This is a backward Euler scheme
    // 1/dt*M*(new_solution_vector - solution_vector) + A*new_solution_vector - rhs = -(1/dt*M + A)*c + A*solution_vector - rhs = -W*c + residual = 0

    inArgs.set_alpha(one/dt);
    inArgs.set_beta(one);

    // set up the jacobian
    RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian = physicsME->create_W();
    outArgs.set_W(jacobian);

    physicsME->evalModel(inArgs,outArgs);
    if (!resetSolver)
      // compute the jacobian matrix only once
      outArgs.set_W(Teuchos::null);

    // set up the solution vector, and residual
    RCP<Thyra::VectorBase<Scalar> > solution_vec = Thyra::createMember(physicsME->get_x_space());
    // initial condition is zero
    Thyra::assign(solution_vec.ptr(), zero);
    inArgs.set_x(solution_vec);
    RCP<Thyra::VectorBase<Scalar> > correction_vec = Thyra::createMember(physicsME->get_x_space());
    Thyra::assign(correction_vec.ptr(), zero);
    RCP<Thyra::VectorBase<Scalar> > residual = Thyra::createMember(physicsME->get_f_space());
    outArgs.set_f(residual);

    // set up responses
    Thyra::ModelEvaluatorBase::InArgs<Scalar> respInArgs = physicsME->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> respOutArgs = physicsME->createOutArgs();
    respInArgs.set_x(solution_vec);
    for(int i=0;i<respOutArgs.Ng();i++) {
      Teuchos::RCP<Thyra::VectorBase<Scalar> > response = Thyra::createMember(*physicsME->get_g_space(i));
      respOutArgs.set_g(i,response);
    }


    if (exodus_output)
      writeToExodus<Scalar>(0,solution_vec,*physicsME,*stkIOResponseLibrary,*mesh);

    {
      Teuchos::TimeMonitor tMts(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: timestepper")));
      auto time_step_timer = Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: Advance Time Step"));
      auto response_timer = Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: Compute responses"));
      for(int ts = 1; ts < numTimeSteps+1; ts++)
        {
          Teuchos::TimeMonitor adv_time_step_timer(*time_step_timer);
          {
            std::stringstream strStream;
            strStream << std::endl;
            strStream << "**************************************************" << std::endl;
            strStream << "* starting time step " << ts << std::endl;
            (*out) << strStream.str();
          }

          double time = dt*ts;
          inArgs.set_t(time);

          // construct the residual, and also the jacobian if (resetSolver == true)
          physicsME->evalModel(inArgs, outArgs);

          // solve
          if (doSolveTimings)
            for (int rep = 0; rep < numReps; rep++) {
              Thyra::assign(correction_vec.ptr(), zero);
              jacobian->solve(Thyra::NOTRANS, *residual, correction_vec.ptr());
            }
          else
            jacobian->solve(Thyra::NOTRANS, *residual, correction_vec.ptr());
          Thyra::V_StVpStV(solution_vec.ptr(), one, *solution_vec, -one, *correction_vec);

          // compute responses
          if (respOutArgs.Ng() > 0) {
            Teuchos::TimeMonitor response_tm(*response_timer);

            respInArgs.set_t(time);
            physicsME->evalModel(respInArgs, respOutArgs);

            std::stringstream strStream;
            for (auto elem: responseIndexToName) {
              Teuchos::RCP<Thyra::VectorBase<Scalar> > g = respOutArgs.get_g(elem.first);
              std::string responseName = elem.second;
              std::transform(responseName.begin(), responseName.end(), responseName.begin(), ::toupper);
              if ((responseName.find("ERROR") != std::string::npos) ||
                  (responseName.find("NORM") != std::string::npos)) {
                strStream << elem.second << " = " << std::sqrt(Thyra::get_ele(*g,0)) << std::endl;
                if (elem.second == "L2 Error E maxwell - analyticSolution")
                  TEUCHOS_ASSERT_INEQUALITY(std::sqrt(Thyra::get_ele(*g,0)), <, 0.065);
              } else
                strStream << elem.second << " = " << Thyra::get_ele(*g,0) << std::endl;
            }
            (*out) << strStream.str();
          }

          // write to an exodus file
          if (exodus_output) {
            Teuchos::TimeMonitor tMexodus(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: timestepper: writeToExodus")));
            writeToExodus<Scalar>(time,solution_vec,*physicsME,*stkIOResponseLibrary,*mesh);
          }

          {
            std::stringstream strStream;
            strStream << std::endl;
            strStream << "* finished time step " << ts << ", t = " << time << std::endl;
            strStream << "**************************************************" << std::endl << std::endl;
            (*out) << strStream.str();
          }
        }
    }

    // Collect FOM data before everything goes out of scope
    fom_num_cells = mesh->getEntityCounts(dim);
  }

  // Output timer data
  if (use_stacked_timer) {
    stacked_timer->stop("Mini-EM");
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stacked_timer->report(*out, comm, options);
    auto xmlOut = stacked_timer->reportWatchrXML(test_name + ' ' + std::to_string(comm->getSize()) + " ranks", comm);
    if(xmlOut.length())
      std::cout << "\nAlso created Watchr performance report " << xmlOut << '\n';

    if ( print_fom && (comm->getRank() == 0) ) {
      std::string fom_timer_name = "Mini-EM@Mini-EM: Total Time@Mini-EM: timestepper@Mini-EM: Advance Time Step@Stratimikos: BelosLOWS";
      double fom_time = stacked_timer->getMpiAverageTime(fom_timer_name);
      double fom_count = stacked_timer->getMpiAverageCount(fom_timer_name);

      *out << "\n=================================\n";
      *out << "FOM Calculation\n";
      *out << "=================================\n";
      *out << "  Number of cells = " << fom_num_cells << std::endl;
      *out << "  Time for Belos Linear Solve = " << fom_time << " seconds" <<std::endl;
      *out << "  Number of Time Steps (one linear solve per step) = " << fom_count << std::endl;
      if (fom_time > 0.0)
        *out << "  FOM ( num_cells * num_steps / solver_time / 1000) = "
             << double(fom_num_cells) * fom_count / fom_time / 1000.0
             << " k-cell-steps per second \n";
      *out << "=================================\n\n";
    }

  } else {
    Teuchos::TimeMonitor::summarize(*out,false,true,false,Teuchos::Union,"",true);
  }

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
  solverType solverValues[5] = {AUGMENTATION, MUELU, ML, CG, GMRES};
  const char * solverNames[5] = {"Augmentation", "MueLu", "ML", "CG", "GMRES"};
  solverType solver = MUELU;
  clp.setOption<solverType>("solver",&solver,5,solverValues,solverNames,"Solver that is used");
  // bool useComplex = false;
  // clp.setOption("complex","real",&useComplex);
  clp.recogniseAllOptions(false);
  switch (clp.parse(argc, argv, NULL)) {
    case Teuchos::CommandLineProcessor::PARSE_ERROR:                return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:         break;
  }

  if (solver == ML) {
    TEUCHOS_ASSERT(linAlgebra == linAlgEpetra);
    // TEUCHOS_ASSERT(!useComplex);
  }

  int retVal;
  if (linAlgebra == linAlgTpetra) {
    // if (useComplex) {
// #if defined(HAVE_TPETRA_COMPLEX_DOUBLE)
//       typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,std::complex<double>,int,panzer::GlobalOrdinal> blockedLinObjFactory;
//       retVal = main_<std::complex<double>,int,panzer::GlobalOrdinal,blockedLinObjFactory,true>(clp, argc, argv);
// #else
//       std::cout << std::endl
//                 << "WARNING" << std::endl
//                 << "Tpetra was compiled without Scalar=std::complex<double>." << std::endl << std::endl;
//       return EXIT_FAILURE;
// #endif
//     } else {
      typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal> blockedLinObjFactory;
      retVal = main_<double,int,panzer::GlobalOrdinal,blockedLinObjFactory,true>(clp, argc, argv);
//    }
#ifdef PANZER_HAVE_EPETRA_STACK
  } else if (linAlgebra == linAlgEpetra) {
    // TEUCHOS_ASSERT(!useComplex);
    typedef typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> blockedLinObjFactory;
    retVal = main_<double,int,int,blockedLinObjFactory,false>(clp, argc, argv);
#endif
  } else
    TEUCHOS_ASSERT(false);

  Kokkos::finalize();

  return retVal;
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
