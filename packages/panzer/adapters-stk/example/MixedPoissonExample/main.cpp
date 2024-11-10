// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StackedTimer.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"
#include "Panzer_Response_Functional.hpp"

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"
#include "Panzer_HierarchicParallelism.hpp"

#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Epetra_MpiComm.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "Panzer_STK_Utilities.hpp"
#include "AztecOO.h"
#endif

#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Ifpack2_Factory.hpp"

#include "Example_BCStrategy_Factory.hpp"
#include "Example_ClosureModel_Factory_TemplateBuilder.hpp"
#include "Example_EquationSetFactory.hpp"


#include <sstream>

using Teuchos::RCP;
using Teuchos::rcp;

void testInitialization(const int hgrad_basis_order,
                        const int hdiv_basis_order,
                        const Teuchos::RCP<Teuchos::ParameterList>& ipb,
                        const std::vector<std::string>& eBlockNames,
                        std::vector<panzer::BC>& bcs);

void solveEpetraSystem(panzer::LinearObjContainer & container);
void solveTpetraSystem(panzer::LinearObjContainer & container);

// calls MPI_Init and MPI_Finalize
int main(int argc,char * argv[])
{
   using std::endl;
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using panzer::StrPureBasisPair;
   using panzer::StrPureBasisComp;

   //panzer::HP::inst().overrideSizes(1,1,1);

   {
     const auto stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Panzer MixedPoisson Test"));
     Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
     stackedTimer->start("Mixed Poisson");


     Teuchos::GlobalMPISession mpiSession(&argc,&argv);
     Kokkos::initialize(argc,argv);
     Teuchos::RCP<const Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
     Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
     out.setOutputToRootOnly(0);
     out.setShowProcRank(true);

     // Build command line processor
     ////////////////////////////////////////////////////

     bool useTpetra = true;
     int x_blocks = 1;
     int x_elements=10,y_elements=10,z_elements=10,hgrad_basis_order=1,hdiv_basis_order=1;
     std::string celltype = "Hex"; // or "Tet"
     std::string output_filename="output_";
     int workset_size = 20;
     bool use_shared_mem_for_ad = true;
     bool check_order_for_shared_mem = true;
     bool stacked_timer_output = true;
     bool time_monitor_output = false;

     Teuchos::CommandLineProcessor clp;
     clp.throwExceptions(false);
     clp.setDocString("This example solves mixed Poisson problem (div,grad) with Hex and Tet inline mesh with high order.\n");

     clp.setOption("cell",&celltype); // this is place holder for tet (for now hex only)
     clp.setOption("use-tpetra","use-epetra",&useTpetra);
     clp.setOption("x-blocks",&x_blocks);
     clp.setOption("x-elements",&x_elements);
     clp.setOption("y-elements",&y_elements);
     clp.setOption("z-elements",&z_elements);
     clp.setOption("hgrad-basis-order",&hgrad_basis_order);
     clp.setOption("hdiv-basis-order",&hdiv_basis_order);
     clp.setOption("output-filename",&output_filename);
     clp.setOption("workset-size",&workset_size);
     clp.setOption("use-shared-mem-for-ad","no-use-shared-mem-for-ad",&use_shared_mem_for_ad);
     clp.setOption("check-order","no-check-order",&check_order_for_shared_mem);
     clp.setOption("stacked-timer-output","no-stacked-timer-output",&stacked_timer_output);
     clp.setOption("time-monitor-output","no-time-monitor-output",&time_monitor_output);

     // parse commandline argument
     Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );
     if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return  0;
     if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

     if(!useTpetra){
#ifndef PANZER_HAVE_EPETRA_STACK
        throw std::runtime_error("Trying to run Panzer MixedPoisson Test with Epetra, but Epetra is disabled!");
#endif
     }

     // cuda optimizations
     ////////////////////////////////////////////////////
     // Always use shared memory optimization for residual
     // evaluation. For AD with FADs, we run out of shared memory on
     // this problem for basis order 3 or higher. Disable shared mem
     // use in this case.
     if (check_order_for_shared_mem) {
       if (std::max(hdiv_basis_order,hgrad_basis_order) > 2)
	 use_shared_mem_for_ad = false;
     }
     panzer::HP::inst().setUseSharedMemory(true,use_shared_mem_for_ad);
     if (use_shared_mem_for_ad)
       out << "Use Shared Memory for AD: ON" << std::endl;
     else
       out << "Use Shared Memory for AD: OFF" << std::endl;

     // variable declarations
     ////////////////////////////////////////////////////

     // factory definitions
     Teuchos::RCP<Example::EquationSetFactory> eqset_factory =
       Teuchos::rcp(new Example::EquationSetFactory); // where poison equation is defined
     Example::BCStrategyFactory bc_factory;    // where boundary conditions are defined

     Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
     if      (celltype == "Hex") mesh_factory = Teuchos::rcp(new panzer_stk::CubeHexMeshFactory);
     else if (celltype == "Tet") mesh_factory = Teuchos::rcp(new panzer_stk::CubeTetMeshFactory);
     else
       throw std::runtime_error("not supported celltype argument: try Hex or Tet");

     // construction of uncommitted (no elements) mesh
     ////////////////////////////////////////////////////////

     // set mesh factory parameters
     RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
     pl->set("X Blocks",x_blocks);
     pl->set("Y Blocks",1);
     pl->set("Z Blocks",1);
     pl->set("X Elements",x_elements/x_blocks);
     pl->set("Y Elements",y_elements);
     pl->set("Z Elements",z_elements);
     mesh_factory->setParameterList(pl);

     RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);

     // construct input physics and physics block
     ////////////////////////////////////////////////////////

     Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
     std::vector<panzer::BC> bcs;
     std::vector<RCP<panzer::PhysicsBlock> > physicsBlocks;
     {
        bool build_transient_support = false;

        std::vector<std::string> eBlockNames;
        mesh->getElementBlockNames(eBlockNames);

        testInitialization(hgrad_basis_order,
                           hdiv_basis_order,
                           ipb, eBlockNames, bcs);

        const panzer::CellData volume_cell_data(workset_size, mesh->getCellTopology("eblock-0_0_0"));

        // GobalData sets ostream and parameter interface to physics
        Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

        // Can be overridden by the equation set
        int default_integration_order = 1;

        // the physics block nows how to build and register evaluator with the field manager
        for (const auto& block : eBlockNames) {
          RCP<panzer::PhysicsBlock> pb
            = rcp(new panzer::PhysicsBlock(ipb, block,
                                           default_integration_order,
                                           volume_cell_data,
                                           eqset_factory,
                                           gd,
                                           build_transient_support));

          // we can have more than one physics block, one per element block
          physicsBlocks.push_back(pb);
        }
     }

     // finish building mesh, set required field variables and mesh bulk data
     ////////////////////////////////////////////////////////////////////////

     for (const auto& pb : physicsBlocks) {
        const std::vector<StrPureBasisPair> & blockFields = pb->getProvidedDOFs();

        // insert all fields into a set
        std::set<StrPureBasisPair,StrPureBasisComp> fieldNames;
        fieldNames.insert(blockFields.begin(),blockFields.end());

        // build string for modifiying vectors
        std::vector<std::string> dimenStr(3);
        dimenStr[0] = "X"; dimenStr[1] = "Y"; dimenStr[2] = "Z";

        // add basis to DOF manager: block specific
        std::set<StrPureBasisPair,StrPureBasisComp>::const_iterator fieldItr;
        for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {
           Teuchos::RCP<const panzer::PureBasis> basis = fieldItr->second;
           if(basis->getElementSpace()==panzer::PureBasis::HGRAD)
              mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
           else if(basis->getElementSpace()==panzer::PureBasis::HCURL) {
              for(int i=0;i<basis->dimension();i++)
                 mesh->addCellField(fieldItr->first+dimenStr[i],pb->elementBlockID());
           }
           else if(basis->getElementSpace()==panzer::PureBasis::HDIV) {
              for(int i=0;i<basis->dimension();i++)
                 mesh->addCellField(fieldItr->first+dimenStr[i],pb->elementBlockID());
           }
        }
     }
     mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

     // build DOF Manager and linear object factory
     /////////////////////////////////////////////////////////////

     RCP<panzer::GlobalIndexer> dofManager;
     Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory;

     // build the connection manager
     if(!useTpetra) {
#ifdef PANZER_HAVE_EPETRA_STACK
       const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

       panzer::DOFManagerFactory globalIndexerFactory;
       RCP<panzer::GlobalIndexer> dofManager_int
             = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);
       dofManager = dofManager_int;

       // construct some linear algebra object, build object to pass to evaluators
       linObjFactory = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(comm.getConst(),dofManager_int));
#endif
     }
     else {
       const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

       panzer::DOFManagerFactory globalIndexerFactory;
       RCP<panzer::GlobalIndexer> dofManager_long
             = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);
       dofManager = dofManager_long;

       // construct some linear algebra object, build object to pass to evaluators
       linObjFactory = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal>(comm,dofManager_long));
     }

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

     // Setup STK response library for writing out the solution fields
     ////////////////////////////////////////////////////////////////////////
     Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
        = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,dofManager,linObjFactory));

     {
        // get a vector of all the element blocks
        std::vector<std::string> eBlocks;
        mesh->getElementBlockNames(eBlocks);

        panzer_stk::RespFactorySolnWriter_Builder builder;
        builder.mesh = mesh;
        stkIOResponseLibrary->addResponse("Main Field Output",eBlocks,builder);
     }

     // Setup response library for checking the error in this manufactered solution
     ////////////////////////////////////////////////////////////////////////

     Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > errorResponseLibrary
        = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,dofManager,linObjFactory));

     {
       const int integration_order = 2 * std::max(hgrad_basis_order,hdiv_basis_order) + 1;

       std::vector<std::string> eBlocks;
       mesh->getElementBlockNames(eBlocks);

       panzer::FunctionalResponse_Builder<int,int> builder;
       builder.comm = MPI_COMM_WORLD;
       builder.cubatureDegree = integration_order;
       builder.requiresCellIntegral = true;
       builder.quadPointField = "PHI_ERROR";

       errorResponseLibrary->addResponse("PHI L2 Error",eBlocks,builder);
     }

     // setup closure model
     /////////////////////////////////////////////////////////////

     // Add in the application specific closure model factory
     panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
     Example::ClosureModelFactory_TemplateBuilder cm_builder;
     cm_factory.buildObjects(cm_builder);

     Teuchos::ParameterList closure_models("Closure Models");
     {
        closure_models.sublist("solid").sublist("SOURCE").set<std::string>("Type","SINE SOURCE");
          // SOURCE field is required by the MixedPoissonEquationSet

       // required for error calculation
       closure_models.sublist("solid").sublist("PHI_ERROR").set<std::string>("Type","ERROR_CALC");
       closure_models.sublist("solid").sublist("PHI_ERROR").set<std::string>("Field A","PHI");
       closure_models.sublist("solid").sublist("PHI_ERROR").set<std::string>("Field B","PHI_EXACT");

       closure_models.sublist("solid").sublist("PHI_EXACT").set<std::string>("Type","PHI_EXACT");
     }

     Teuchos::ParameterList user_data("User Data"); // user data can be empty here

     // setup field manager builder
     /////////////////////////////////////////////////////////////

     Teuchos::RCP<panzer::FieldManagerBuilder> fmb =
           Teuchos::rcp(new panzer::FieldManagerBuilder);
     fmb->setWorksetContainer(wkstContainer);
     fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
     fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,
                               *linObjFactory,user_data);

     // setup assembly engine
     /////////////////////////////////////////////////////////////

     // build assembly engine: The key piece that brings together everything and
     //                        drives and controls the assembly process. Just add
     //                        matrices and vectors
     panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
     panzer::AssemblyEngine_TemplateBuilder builder(fmb,linObjFactory);
     ae_tm.buildObjects(builder);

     // Finalize construcition of STK writer response library
     /////////////////////////////////////////////////////////////
     {
        user_data.set<int>("Workset Size",workset_size);
        stkIOResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                          cm_factory,
                                          closure_models,
                                          user_data);

        user_data.set<int>("Workset Size",workset_size);
        errorResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                                      cm_factory,
                                                      closure_models,
                                                      user_data);
     }

     // assemble linear system
     /////////////////////////////////////////////////////////////

     // build linear algebra objects: Ghost is for parallel assembly, it contains
     //                               local element contributions summed, the global IDs
     //                               are not unique. The non-ghosted or "global"
     //                               container will contain the sum over all processors
     //                               of the ghosted objects. The global indices are unique.
     RCP<panzer::LinearObjContainer> ghostCont = linObjFactory->buildGhostedLinearObjContainer();
     RCP<panzer::LinearObjContainer> container = linObjFactory->buildLinearObjContainer();
     linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                               panzer::LinearObjContainer::F |
                                               panzer::LinearObjContainer::Mat,*ghostCont);
     linObjFactory->initializeContainer(panzer::LinearObjContainer::X |
                                        panzer::LinearObjContainer::F |
                                        panzer::LinearObjContainer::Mat,*container);
     ghostCont->initialize();
     container->initialize();

     // Actually evaluate
     /////////////////////////////////////////////////////////////

     panzer::AssemblyEngineInArgs input(ghostCont,container);
     input.alpha = 0;
     input.beta = 1;

     // evaluate physics: This does both the Jacobian and residual at once
     ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

     // solve linear system
     /////////////////////////////////////////////////////////////

     if(useTpetra)
        solveTpetraSystem(*container);
     else
        solveEpetraSystem(*container);

     // output data (optional)
     /////////////////////////////////////////////////////////////

     // write out solution
     if(true) {
        // fill STK mesh objects
        panzer::AssemblyEngineInArgs respInput(ghostCont,container);
        respInput.alpha = 0;
        respInput.beta = 1;

        stkIOResponseLibrary->addResponsesToInArgs<panzer::Traits::Residual>(respInput);
        stkIOResponseLibrary->evaluate<panzer::Traits::Residual>(respInput);

        // write to exodus
        // ---------------
        // Due to multiple instances of this test being run at the
        // same time (one for each order), we need to differentiate
        // output to prevent race conditions on output file. Multiple
        // runs for the same order are ok as they are staged one after
        // another in the ADD_ADVANCED_TEST cmake macro.
        std::ostringstream filename;
        filename << output_filename << hgrad_basis_order << "_" << hdiv_basis_order << ".exo";
        mesh->writeToExodus(filename.str());
     }

     // compute error norm
     /////////////////////////////////////////////////////////////

     if(true) {
        Teuchos::FancyOStream lout(Teuchos::rcpFromRef(std::cout));
        lout.setOutputToRootOnly(0);

        panzer::AssemblyEngineInArgs respInput(ghostCont,container);
        respInput.alpha = 0;
        respInput.beta = 1;

        Teuchos::RCP<panzer::ResponseBase> resp = errorResponseLibrary->getResponse<panzer::Traits::Residual>("PHI L2 Error");
        Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > resp_func =
               Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(resp);
        Teuchos::RCP<Thyra::VectorBase<double> > respVec = Thyra::createMember(resp_func->getVectorSpace());
        resp_func->setVector(respVec);

        errorResponseLibrary->addResponsesToInArgs<panzer::Traits::Residual>(respInput);
        errorResponseLibrary->evaluate<panzer::Traits::Residual>(respInput);

        lout << "HGrad Basis Order = " << hgrad_basis_order << std::endl;
        lout << "HDiv Basis Order = " << hdiv_basis_order << std::endl;
        lout << "Error = " << sqrt(resp_func->value) << std::endl;
     }

     stackedTimer->stop("Mixed Poisson");
     stackedTimer->stopBaseTimer();
     if (stacked_timer_output) {
       Teuchos::StackedTimer::OutputOptions options;
       options.output_fraction = true;
       options.output_minmax = true;
       options.output_histogram = true;
       options.num_histogram = 5;
       stackedTimer->report(std::cout, Teuchos::DefaultComm<int>::getComm(), options);
     }
     if (time_monitor_output) {
       Teuchos::TimeMonitor::summarize(out,false,true,false,Teuchos::Union);
     }

     // all done!
     /////////////////////////////////////////////////////////////

     if(useTpetra)
        out << "ALL PASSED: Tpetra" << std::endl;
     else
        out << "ALL PASSED: Epetra" << std::endl;
   }

   return 0;
}

void solveEpetraSystem(panzer::LinearObjContainer & container)
{
#ifdef PANZER_HAVE_EPETRA_STACK
   // convert generic linear object container to epetra container
   panzer::EpetraLinearObjContainer & ep_container
         = Teuchos::dyn_cast<panzer::EpetraLinearObjContainer>(container);

   // Setup the linear solve: notice A is used directly
   Epetra_LinearProblem problem(&*ep_container.get_A(),&*ep_container.get_x(),&*ep_container.get_f());

   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setShowProcRank(true);
   out.setOutputToRootOnly(-1);
   out << "SIZE = " << ep_container.get_x()->MyLength() << "/"
                    << ep_container.get_x()->GlobalLength() << std::endl;

   // build the solver
   AztecOO solver(problem);
   solver.SetAztecOption(AZ_solver,AZ_gmres); // we don't push out dirichlet conditions
   solver.SetAztecOption(AZ_kspace,150); // something else might be better but there is a
                                         // chance to much memory is allocated
   solver.SetAztecOption(AZ_output,1);
   solver.SetAztecOption(AZ_precond,AZ_Jacobi);

   // solve the linear system
   solver.Iterate(1000,1e-9);

   // we have now solved for the residual correction from
   // zero in the context of a Newton solve.
   //     J*e = -r = -(f - J*0) where f = J*u
   // Therefore we have  J*e=-J*u which implies e = -u
   // thus we will scale the solution vector
   ep_container.get_x()->Scale(-1.0);
#endif
}

void solveTpetraSystem(panzer::LinearObjContainer & container)
{
  typedef panzer::TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> LOC;

  LOC & tp_container = Teuchos::dyn_cast<LOC>(container);

  // do stuff
  // Wrap the linear problem to solve in a Belos::LinearProblem
  // object.  The "X" argument of the LinearProblem constructor is
  // only copied shallowly and will be overwritten by the solve, so we
  // make a deep copy here.  That way we can compare the result
  // against the original X_guess.
  typedef Tpetra::MultiVector<double,int,panzer::GlobalOrdinal> MV;
  typedef Tpetra::Operator<double,int,panzer::GlobalOrdinal> OP;
  typedef Belos::LinearProblem<double,MV, OP> ProblemType;
  Teuchos::RCP<ProblemType> problem(new ProblemType(tp_container.get_A(), tp_container.get_x(), tp_container.get_f()));
  auto prec = Ifpack2::Factory::create<Tpetra::RowMatrix<double,int,panzer::GlobalOrdinal>>("RELAXATION",tp_container.get_A(),4);
  Teuchos::ParameterList precParams;
  precParams.set("relaxation: type", "Gauss-Seidel");
  prec->setParameters(precParams);
  prec->initialize();
  prec->compute();
  problem->setRightPrec(prec);
  TEUCHOS_ASSERT(problem->setProblem());

  typedef Belos::PseudoBlockGmresSolMgr<double,MV,OP> SolverType;

  Teuchos::ParameterList belosList;
  belosList.set( "Num Blocks", 3000 );            // Maximum number of blocks in Krylov factorization
  belosList.set( "Block Size", 1 );              // Blocksize to be used by iterative solver
  belosList.set( "Maximum Iterations", 3000 );       // Maximum number of iterations allowed
  belosList.set( "Maximum Restarts", 1 );      // Maximum number of restarts allowed
  belosList.set( "Convergence Tolerance", 1e-9 );         // Relative convergence tolerance requested
  // belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
  belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails );
  belosList.set( "Output Frequency", 1 );
  belosList.set( "Output Style", 1 );

  SolverType solver(problem, Teuchos::rcpFromRef(belosList));

  Belos::ReturnType result = solver.solve();
  if (result == Belos::Converged)
    std::cout << "Result: Converged." << std::endl;
  else {
    TEUCHOS_ASSERT(false); // FAILURE!
  }

  // scale by -1
  tp_container.get_x()->scale(-1.0);

  tp_container.get_A()->resumeFill(); // where does this go?
}

void testInitialization(const int hgrad_basis_order,
                        const int hdiv_basis_order,
                        const Teuchos::RCP<Teuchos::ParameterList>& ipb,
                        const std::vector<std::string>& eBlockNames,
                        std::vector<panzer::BC>& bcs)
{
  {
    const int integration_order = 2 * std::max(hgrad_basis_order,hdiv_basis_order) + 1;
    Teuchos::ParameterList& p = ipb->sublist("MixedPoisson Physics");
    p.set("Type","MixedPoisson");
    p.set("Model ID","solid");
    p.set("HGrad Basis Order",hgrad_basis_order);
    p.set("HDiv Basis Order",hdiv_basis_order);
    p.set("Integration Order",integration_order);
  }

  std::size_t bc_id = 0;
  for (const auto& block : eBlockNames) {
    panzer::BCType bctype = panzer::BCT_Dirichlet;
    std::string sideset_id = "left";
    std::string element_block_id = block;
    std::string dof_name = "PHI";
    std::string strategy = "Constant";
    double value = 0.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name,
		  strategy, p);
    bcs.push_back(bc);
  }


  for (const auto& block : eBlockNames) {
    panzer::BCType bctype = panzer::BCT_Dirichlet;
    std::string sideset_id = "top";
    std::string element_block_id = block;
    std::string dof_name = "PHI";
    std::string strategy = "Constant";
    double value = 0.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name,
		  strategy, p);
    bcs.push_back(bc);
  }

  for (const auto& block : eBlockNames) {
    panzer::BCType bctype = panzer::BCT_Dirichlet;
    std::string sideset_id = "right";
    std::string element_block_id = block;
    std::string dof_name = "PHI";
    std::string strategy = "Constant";
    double value = 0.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name,
		  strategy, p);
    bcs.push_back(bc);
  }

  for (const auto& block : eBlockNames) {
    panzer::BCType bctype = panzer::BCT_Dirichlet;
    std::string sideset_id = "bottom";
    std::string element_block_id = block;
    std::string dof_name = "PHI";
    std::string strategy = "Constant";
    double value = 0.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name,
		  strategy, p);
    bcs.push_back(bc);
  }

  for (const auto& block : eBlockNames) {
    panzer::BCType bctype = panzer::BCT_Dirichlet;
    std::string sideset_id = "front";
    std::string element_block_id = block;
    std::string dof_name = "PHI";
    std::string strategy = "Constant";
    double value = 0.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name,
		  strategy, p);
    bcs.push_back(bc);
  }

  for (const auto& block : eBlockNames) {
    panzer::BCType bctype = panzer::BCT_Dirichlet;
    std::string sideset_id = "back";
    std::string element_block_id = block;
    std::string dof_name = "PHI";
    std::string strategy = "Constant";
    double value = 0.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name,
		  strategy, p);
    bcs.push_back(bc);
  }
}
