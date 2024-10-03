// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <iostream>
#include <string>

// NOX
#include "NOX_Thyra.H"

// Panzer
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_NodeType.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_IOClosureModel_Factory_TemplateBuilder.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"
#include "Panzer_STK_SetupLOWSFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_String_Utilities.hpp"

// Seacass
#include <Ioss_SerializeIO.h>

// Teuchos
#include "Teuchos_as.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Files specifically for this example.
#include "myBCStrategyFactory.hpp"
#include "myClosureModelFactory_TemplateBuilder.hpp"
#include "myEquationSetFactory.hpp"

/**
 * 	\brief Build the STK I/O Response Library.
 *
 * 	Create a `ResponseLibrary` and add a response for the field to output, then
 * 	build the corresponding response `Evaluators`.
 *
 * 	\param[in] physicsBlocks  A list of all the physics blocks for our problem
 * 	                          (of which there is only one).
 * 	\param[in] linObjFactory  The linear object factory, which keeps track of
 * 	                          all the linear algebra pieces (vectors/matrices)
 * 	                          of our problem.
 * 	\param[in] wkstContainer  A container holding all of our worksets, which
 * 	                          are groups of cells (elements) that live on the
 * 	                          same process.
 * 	\param[in] globalIndexer  Think of this as our degree of freedom manager,
 * 	                          which handles all the indexing of our unknowns.
 * 	\param[in] cmFactory      The closure model factory, through which we'll
 * 	                          specify our source term.
 * 	\param[in] mesh           Our mesh data structure; that is, our concrete
 * 	                          implementation of our connection manager.
 * 	\param[in] closureModelPl The "Closure Models" `ParameterList` from the
 * 	                          input XML file.
 * 	\param[in] userData       The "User Data" `ParameterList`, which holds the
 * 	                          MPI communicator.
 *
 *  \returns The `ResponseLibrary`, which will allow us to write out the
 *           results at the end of the run.
 */
Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>>
buildSTKIOResponseLibrary(
  const std::vector<Teuchos::RCP<panzer::PhysicsBlock>>&        physicsBlocks,
  const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits>>& linObjFactory,
  const Teuchos::RCP<panzer::WorksetContainer>&                 wkstContainer,
  const Teuchos::RCP<panzer::GlobalIndexer>&          globalIndexer,
  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cmFactory,
  const Teuchos::RCP<panzer_stk::STK_Interface>&                mesh,
  const Teuchos::ParameterList&                                 closureModelPl,
  const Teuchos::ParameterList&                                 userData);

/**
 * 	\brief Write the results to an Exodus file.
 *
 * 	Get the responses, evaluate the model, and write out the results to an
 * 	Exodus file, which can be viewed with ParaView (www.paraview.org).
 *
 * 	\param[in]     x                    The solution vector.
 * 	\param[in]     model                The `ModelEvaluator` representing the
 * 	                                    problem we're solving.
 * 	\param[in/out] stkIOResponseLibrary Our response library, which is able to
 * 	                                    evaluate the fields we want to output.
 * 	\param[in/out] mesh                 Our mesh database, which does the
 * 	                                    actual writing to file.
 */
void
writeToExodus(
  const Teuchos::RCP<const Thyra::VectorBase<double>>& x,
  const panzer::ModelEvaluator<double>&                model,
  panzer::ResponseLibrary<panzer::Traits>&             stkIOResponseLibrary,
  panzer_stk::STK_Interface&                           mesh);

///////////////////////////////////////////////////////////////////////////////
//
//  main()
//
///////////////////////////////////////////////////////////////////////////////
int
main(
  int   argc,
  char* argv[])
{
  using   panzer::BC;
  using   panzer::buildBCs;
  using   panzer::buildBlockIdToPhysicsIdMap;
  using   panzer::buildPhysicsBlocks;
  using   panzer::ClosureModelFactory_TemplateManager;
  using   panzer::ConnManager;
  using   panzer::createGlobalData;
  using   panzer::DOFManagerFactory;
  using   panzer::BlockedEpetraLinearObjFactory;
  using   panzer::GlobalData;
  using   panzer::LinearObjFactory;
  using   panzer::PhysicsBlock;
  using   panzer::PureBasis;
  using   panzer::ResponseLibrary;
  using   panzer::StrPureBasisComp;
  using   panzer::StrPureBasisPair;
  using   panzer::Traits;
  using   panzer::GlobalIndexer;
  using   panzer::WorksetContainer;
  using   panzer_stk::buildLOWSFactory;
  using   panzer_stk::SquareQuadMeshFactory;
  using   panzer_stk::STKConnManager;
  using   panzer_stk::STK_Interface;
  using   panzer_stk::STK_MeshFactory;
  using   panzer_stk::WorksetFactory;
  using   shards::CellTopology;
  using   std::cout;
  using   std::endl;
  using   std::exception;
  using   std::map;
  using   std::runtime_error;
  using   std::set;
  using   std::size_t;
  using   std::string;
  using   std::vector;
  using   Teuchos::as;
  using   Teuchos::Comm;
  using   Teuchos::CommandLineProcessor;
  using   Teuchos::DefaultComm;
  using   Teuchos::FancyOStream;
  using   Teuchos::GlobalMPISession;
  using   Teuchos::oblackholestream;
  using   Teuchos::MpiComm;
  using   Teuchos::null;
  using   Teuchos::opaqueWrapper;
  using   Teuchos::ParameterList;
  using   Teuchos::RCP;
  using   Teuchos::rcp;
  using   Teuchos::rcp_dynamic_cast;
  using   Teuchos::Time;
  using   Teuchos::TimeMonitor;
  using   Teuchos::updateParametersFromXmlFileAndBroadcast;
  using   Thyra::assign;
  using   Thyra::createMember;
  using   Thyra::LinearOpWithSolveBase;
  using   Thyra::LinearOpWithSolveFactoryBase;
  using   Thyra::scale;
  using   Thyra::VectorBase;
  typedef panzer::ModelEvaluator<double>             PME;
  typedef Thyra::ModelEvaluatorBase::InArgs<double>  InArgs;
  typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
  int status(0);

  // Initialize Kokkos/MPI.
  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc, &argv, &blackhole);
  Kokkos::initialize(argc, argv);

  // Set up the fancy output stream.
  RCP<FancyOStream> out = rcp(new FancyOStream(rcp(&cout, false)));
  if (mpiSession.getNProc() > 1)
  {
    out->setShowProcRank(true);
    out->setOutputToRootOnly(0);
  }

  // Now on to the main part of our program.
  try
  {
    // Create a TimeMonitor to track all our timing information.
    RCP<Time> totalTime = TimeMonitor::getNewTimer("Total Time");
    TimeMonitor timer(*totalTime);

    // Get the MPI communicator.
    RCP<const MpiComm<int>> comm =
      rcp_dynamic_cast<const MpiComm<int>>(DefaultComm<int>::getComm());

    // Parse the command-line arguments.
    string inputFileName("input.xml");
    {
      CommandLineProcessor clp;
      clp.setOption("i", &inputFileName, "Input XML filename");
      CommandLineProcessor::EParseCommandLineReturn parseReturn =
        clp.parse(argc, argv, &std::cerr);
      TEUCHOS_TEST_FOR_EXCEPTION(
        parseReturn != CommandLineProcessor::PARSE_SUCCESSFUL, runtime_error,
        "Failed to parse the command line!");
    } // end of "Parse the command-line arguments"

    // Parse the input file and broadcast to the other processes.
    RCP<ParameterList> inputParams =
      rcp(new ParameterList("Input Parameters"));
    updateParametersFromXmlFileAndBroadcast(inputFileName, inputParams.ptr(),
      *comm);

    // Separate out the various blocks of the input parameters.
    RCP<ParameterList> meshPl          =
      rcp(new ParameterList(inputParams->sublist("Mesh")));
    RCP<ParameterList> physicsBlocksPl =
      rcp(new ParameterList(inputParams->sublist("Physics Blocks")));
    RCP<ParameterList> linSolverPl     =
      rcp(new ParameterList(inputParams->sublist("Linear Solver")));
    ParameterList& blockToPhysicsPl    =
      inputParams->sublist("Block ID to Physics ID Mapping");
    ParameterList& bcsPl               =
      inputParams->sublist("Boundary Conditions");
    ParameterList& closureModelsPl     =
      inputParams->sublist("Closure Models");
    ParameterList& userDataPl          =
      inputParams->sublist("User Data");

    // Store the MPI communicator in the "User Data" ParameterList.
    userDataPl.set<RCP<const Comm<int>>>("Comm", comm);

    // Create what's needed for the equation sets, boundary conditions, and
    // closure models.
    RCP<GlobalData> globalData = createGlobalData();
    RCP<MyEquationSetFactory> eqSetFactory = rcp(new MyEquationSetFactory);
    MyBCStrategyFactory bcFactory;
    MyClosureModelFactory_TemplateBuilder cmBuilder;
    ClosureModelFactory_TemplateManager<Traits> cmFactory;
    cmFactory.buildObjects(cmBuilder);

    // Read in the mesh database and build the un-committed mesh.
    RCP<STK_MeshFactory> meshFactory = rcp(new SquareQuadMeshFactory);
    meshFactory->setParameterList(meshPl);
    RCP<STK_Interface> mesh = meshFactory->buildUncommitedMesh(MPI_COMM_WORLD);

    // Read in the element block to physics block mapping, and create a
    // corresponding element block to cell topology mapping.
    map<string, string> blockIdsToPhysicsIds;
    buildBlockIdToPhysicsIdMap(blockIdsToPhysicsIds, blockToPhysicsPl);
    map<string, RCP<const CellTopology>> blockIdsToCellTopo;
    for (auto itr(blockIdsToPhysicsIds.begin());
      itr != blockIdsToPhysicsIds.end(); ++itr)
      blockIdsToCellTopo[itr->first] = mesh->getCellTopology(itr->first);

    // Build the physics blocks, using some reasonable defaults in case values
    // are omitted in the input parameters.
    vector<RCP<PhysicsBlock>> physicsBlocks;
    int worksetSize(20), defaultIntegrationOrder(2);
    bool buildTransientSupport(false);
    vector<string> tangentParamNames;
    buildPhysicsBlocks(blockIdsToPhysicsIds, blockIdsToCellTopo,
      physicsBlocksPl, defaultIntegrationOrder, worksetSize, eqSetFactory,
      globalData, buildTransientSupport, physicsBlocks, tangentParamNames);

    // Loop over the physics blocks.
    for (size_t i(0); i < physicsBlocks.size(); ++i)
    {
      // Get the degrees of freedom from the current physics block and insert
      // them all into a set to enforce uniqueness.
      RCP<PhysicsBlock> pb = physicsBlocks[i];
      const vector<StrPureBasisPair>& blockFields = pb->getProvidedDOFs();
      set<StrPureBasisPair, StrPureBasisComp>
        fieldNames(blockFields.begin(), blockFields.end());

      // Loop over the degrees of freedom, adding their bases to the mesh
      // database.
      set<StrPureBasisPair, StrPureBasisComp>::const_iterator field;
      for (field = fieldNames.begin(); field != fieldNames.end(); ++field)
        mesh->addSolutionField(field->first, pb->elementBlockID());
      meshFactory->completeMeshConstruction(*mesh, MPI_COMM_WORLD);
    } // end loop over physicsBlocks

    // Build the connection manager, which is an abstraction of the mesh.
    const RCP<ConnManager> connManager =
      rcp(new STKConnManager(mesh));

    // Build the degree of freedom manager and the linear object factory, which
    // creates the linear algebra objects (vectors/matrices) for the problem.
    RCP<GlobalIndexer> dofManager;
    RCP<LinearObjFactory<Traits>> linObjFactory;
    {
      DOFManagerFactory<int, int> globalIndexerFactory;
      dofManager = globalIndexerFactory.buildGlobalIndexer(
        opaqueWrapper(MPI_COMM_WORLD), physicsBlocks, connManager);
      linObjFactory =
        rcp(new BlockedEpetraLinearObjFactory<Traits, int>(comm, dofManager));
    }

    // Build the STK workset factory and attach it to a workset container.
    RCP<WorksetFactory> wkstFactory = rcp(new WorksetFactory(mesh));
    RCP<WorksetContainer> wkstContainer = rcp(new WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for (size_t i(0); i < physicsBlocks.size(); ++i) 
      wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),
        physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(worksetSize);
    wkstContainer->setGlobalIndexer(dofManager);

    // Build the linear solver we'll use to solve the system.
    RCP<LinearOpWithSolveFactoryBase<double>> lowsFactory =
      buildLOWSFactory(false, dofManager, connManager,
      as<int>(mesh->getDimension()), comm, linSolverPl, null);

    // Set up the boundary conditions.
    vector<BC> bcs;
    buildBCs(bcs, bcsPl, globalData);

    // Build and set up the ModelEvaluator.
    RCP<PME> physics = rcp(new PME(linObjFactory, lowsFactory, globalData,
      buildTransientSupport, 0.0));
    physics->setupModel(wkstContainer, physicsBlocks, bcs, *eqSetFactory,
      bcFactory, cmFactory, cmFactory, closureModelsPl, userDataPl, false, "");

    // Set up a response library to write to the mesh.
    RCP<ResponseLibrary<Traits>> stkIOResponseLibrary =
      buildSTKIOResponseLibrary(physicsBlocks, linObjFactory, wkstContainer,
      dofManager, cmFactory, mesh, closureModelsPl, userDataPl);

    // Allocate the vectors and matrix for the linear solve.
    RCP<VectorBase<double>> solutionVec = createMember(physics->get_x_space());
    assign(solutionVec.ptr(), 0.0);
    RCP<VectorBase<double>> residual = createMember(physics->get_f_space());
    RCP<LinearOpWithSolveBase<double>> jacobian = physics->create_W();

    // Do the problem assembly (this is where the evaluators are called and the
    // directed acyclic graph is execueted).
    {
      InArgs inArgs = physics->createInArgs();
      inArgs.set_x(solutionVec);
      OutArgs outArgs = physics->createOutArgs();
      outArgs.set_f(residual);
      outArgs.set_W(jacobian);

      // Construct the residual and Jacobian.
      physics->evalModel(inArgs, outArgs);
    }

    // Do a linear solve.
    auto status(jacobian->solve(Thyra::NOTRANS, *residual, solutionVec.ptr()));
    scale(-1.0, solutionVec.ptr());

    // Write the results to an Exodus file.
    writeToExodus(solutionVec, *physics, *stkIOResponseLibrary, *mesh);
  } // end of try

  // Catch any exceptions that were thrown along the way.
  catch (exception& e)
  {
    *out << "*********** Caught Exception: Begin Error Report ***********"
         << endl << e.what() << endl
         << "************ Caught Exception: End Error Report ************"
         << endl;
    status = -1;
  } // end of catch (exception& e)
  catch (string& msg)
  {
    *out << "*********** Caught Exception: Begin Error Report ***********"
         << endl << msg << endl
         << "************ Caught Exception: End Error Report ************"
         << endl;
    status = -1;
  } // end of catch (string& msg)
  catch (...)
  {
    *out << "*********** Caught Exception: Begin Error Report ***********"
         << endl << "Caught UNKOWN exception" << endl
         << "************ Caught Exception: End Error Report ************"
         << endl;
    status = -1;
  } // end of catch (...)

  // Print out the timing information.
  Teuchos::TimeMonitor::summarize(*out, false, true, false);

  // Print out if we were successful.
  if (status == 0)
    *out << "Run completed." << endl;

  // Shut things down.
  Kokkos::finalize_all();
  return status;
} // end of main()

///////////////////////////////////////////////////////////////////////////////
//
//  buildSTKIOResponseLibrary()
//
///////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>>
buildSTKIOResponseLibrary(
  const std::vector<Teuchos::RCP<panzer::PhysicsBlock>>&        physicsBlocks,
  const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits>>& linObjFactory,
  const Teuchos::RCP<panzer::WorksetContainer>&                 wkstContainer,
  const Teuchos::RCP<panzer::GlobalIndexer>&          globalIndexer,
  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cmFactory,
  const Teuchos::RCP<panzer_stk::STK_Interface>&                mesh,
  const Teuchos::ParameterList&                                 closureModelPl,
  const Teuchos::ParameterList&                                 userData)
{
  using panzer::ClosureModelFactory_TemplateManager;
  using panzer::ResponseLibrary;
  using panzer::Traits;
  using panzer_stk::IOClosureModelFactory_TemplateBuilder;
  using panzer_stk::RespFactorySolnWriter_Builder;
  using std::map;
  using std::string;
  using std::vector;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Create the ResponseLibrary and add a response for the field output.
  RCP<ResponseLibrary<Traits>> stkIOResponseLibrary = rcp(new
    ResponseLibrary<Traits>(wkstContainer, globalIndexer, linObjFactory));
  vector<string> eBlocks;
  mesh->getElementBlockNames(eBlocks);
  RespFactorySolnWriter_Builder builder;
  builder.mesh = mesh;
  stkIOResponseLibrary->addResponse("Main Field Output", eBlocks, builder);

  // Build the response Evaluators.
  map<string, vector<string>> nodalFields, cellFields;
  IOClosureModelFactory_TemplateBuilder<Traits> ioCmBuilder(cmFactory, mesh,
    nodalFields, cellFields);
  ClosureModelFactory_TemplateManager<Traits> ioCmFactory;
  ioCmFactory.buildObjects(ioCmBuilder);
  stkIOResponseLibrary->buildResponseEvaluators(physicsBlocks, ioCmFactory,
    closureModelPl, userData);
  return stkIOResponseLibrary;
} // end of buildSTKIOResponseLibrary()

///////////////////////////////////////////////////////////////////////////////
//
//  writeToExodus()
//
///////////////////////////////////////////////////////////////////////////////
void
writeToExodus(
  const Teuchos::RCP<const Thyra::VectorBase<double>>& x,
  const panzer::ModelEvaluator<double>&                model,
  panzer::ResponseLibrary<panzer::Traits>&             stkIOResponseLibrary,
  panzer_stk::STK_Interface&                           mesh)
{
  using   panzer::AssemblyEngineInArgs;
  typedef panzer::Traits::Residual Residual;
  typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;

  // Get the responses.
  InArgs inArgs = model.createInArgs();
  inArgs.set_x(x);
  AssemblyEngineInArgs respInput;
  model.setupAssemblyInArgs(inArgs, respInput);

  // Evaluate the model and write out the results.
  stkIOResponseLibrary.addResponsesToInArgs<Residual>(respInput);
  stkIOResponseLibrary.evaluate<Residual>(respInput);
  mesh.writeToExodus("output.exo");
} // end of writeToExodus()

// end of main.cpp
