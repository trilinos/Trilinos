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

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
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
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"

#include "Example_BCStrategy_Factory.hpp"
#include "Example_ClosureModel_Factory_TemplateBuilder.hpp"
#include "Example_EquationSetFactory.hpp"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosMultiVecTraits_Tpetra.hpp"
#include "BelosOperatorTraits_Tpetra.hpp"

#include "MatrixMarket_Tpetra.hpp"

#include <sstream>
#include <fstream>
#include <cmath>

using Teuchos::RCP;
using Teuchos::rcp;

void testInitialization(const int basis_order,
                        const Teuchos::RCP<Teuchos::ParameterList>& ipb,
                        std::vector<panzer::BC>& bcs,
                        const bool curvilinear);

// calls MPI_Init and MPI_Finalize
int main(int argc,char * argv[])
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using panzer::StrPureBasisPair;
   using panzer::StrPureBasisComp;

   const auto stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Panzer MixedPoisson Test"));
   Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
   stackedTimer->start("Mixed Poisson");

   Teuchos::GlobalMPISession mpiSession(&argc,&argv);
   Kokkos::initialize(argc,argv);
   RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setOutputToRootOnly(0);
   out.setShowProcRank(true);

   // Build command line processor
   ////////////////////////////////////////////////////

   int x_elements=10,y_elements=10,basis_order=1;
   std::string celltype = "Quad"; // or "Tri"
   std::string mesh_name = "";
   std::string problem_name = "rectangle";
   Teuchos::CommandLineProcessor clp;
   clp.throwExceptions(false);
   clp.setDocString("This example solves a Poisson problem with Quad and Tri inline mesh on a square domain"
       " or with a user-supplied mesh on an annular domain.\n");

   clp.setOption("cell",&celltype);         // ignored if mesh file is supplied
   clp.setOption("x-elements",&x_elements); // ignored if mesh file is supplied
   clp.setOption("y-elements",&y_elements); // ignored if mesh file is supplied
   clp.setOption("basis-order",&basis_order);
   clp.setOption("mesh-filename",&mesh_name);
   clp.setOption("problem",&problem_name,"Problem type: rectangle or annulus. Defaults to rectangle.");

   // parse commandline argument
   Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );
   if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return  0;
   if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

   // two problems are supported -- one on a rectangular domain, one on an annular domain
   bool curvilinear = false;
   if      (problem_name == "rectangle") {}
   else if (problem_name == "annulus") {curvilinear = true;}
   else
     throw std::runtime_error("Problem not supported: try rectangle or annulus");

   // variable declarations
   ////////////////////////////////////////////////////

   // factory definitions
   Teuchos::RCP<Example::EquationSetFactory> eqset_factory =
     Teuchos::rcp(new Example::EquationSetFactory); // where poison equation is defined
   Example::BCStrategyFactory bc_factory;    // where boundary conditions are defined

   Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
   if (mesh_name == "") {
    if      (celltype == "Quad") mesh_factory = Teuchos::rcp(new panzer_stk::SquareQuadMeshFactory);
    else if (celltype == "Tri")  mesh_factory = Teuchos::rcp(new panzer_stk::SquareTriMeshFactory);
    else
      throw std::runtime_error("not supported celltype argument: try Quad or Tri");
    // construction of uncommitted (no elements) mesh
    ////////////////////////////////////////////////////////

    // set mesh factory parameters
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",x_elements);
    pl->set("Y Elements",y_elements);
    mesh_factory->setParameterList(pl);
   } else {
    mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory(mesh_name));
   }

   // other declarations
   const std::size_t workset_size = 2000;

   RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);

   // construct input physics and physics block
   ////////////////////////////////////////////////////////

   Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
   std::vector<panzer::BC> bcs;
   std::vector<RCP<panzer::PhysicsBlock> > physicsBlocks;
   {
      bool build_transient_support = false;

      testInitialization(basis_order, ipb, bcs, curvilinear);

      const panzer::CellData volume_cell_data(workset_size, mesh->getCellTopology("eblock-0_0"));

      // GobalData sets ostream and parameter interface to physics
      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      // Can be overridden by the equation set
      int default_integration_order = 1;

      // the physics block nows how to build and register evaluator with the field manager
      RCP<panzer::PhysicsBlock> pb
	= rcp(new panzer::PhysicsBlock(ipb,
				       "eblock-0_0",
				       default_integration_order,
				       volume_cell_data,
				       eqset_factory,
				       gd,
				       build_transient_support));

      // we can have more than one physics block, one per element block
      physicsBlocks.push_back(pb);
   }

   // finish building mesh, set required field variables and mesh bulk data
   ////////////////////////////////////////////////////////////////////////

   {
      RCP<panzer::PhysicsBlock> pb = physicsBlocks[0]; // we are assuming only one physics block

      const std::vector<StrPureBasisPair> & blockFields = pb->getProvidedDOFs();

      // insert all fields into a set
      std::set<StrPureBasisPair,StrPureBasisComp> fieldNames;
      fieldNames.insert(blockFields.begin(),blockFields.end());

      // add basis to DOF manager: block specific
      std::set<StrPureBasisPair,StrPureBasisComp>::const_iterator fieldItr;
      for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr)
         mesh->addSolutionField(fieldItr->first,pb->elementBlockID());

      mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);
   }

   // build DOF Manager and linear object factory
   /////////////////////////////////////////////////////////////

   // build the connection manager
   const Teuchos::RCP<panzer::ConnManager> conn_manager =
     Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

   panzer::DOFManagerFactory globalIndexerFactory;
   RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

   // construct some linear algebra object, build object to pass to evaluators
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
     = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,panzer::LocalOrdinal,panzer::GlobalOrdinal>(tComm,dofManager));

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

   // Setup response library for checking the error in this manufactured solution and computing area
   ////////////////////////////////////////////////////////////////////////

   Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > exampleResponseLibrary
      = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,dofManager,linObjFactory));

   {
     const int integration_order = 10;

     std::vector<std::string> eBlocks;
     mesh->getElementBlockNames(eBlocks);

     panzer::FunctionalResponse_Builder<int,int> builder;
     builder.comm = MPI_COMM_WORLD;
     builder.cubatureDegree = integration_order;
     builder.requiresCellIntegral = true;

     builder.quadPointField = "TEMPERATURE_L2_ERROR";
     exampleResponseLibrary->addResponse("L2 Error",eBlocks,builder);

     builder.quadPointField = "TEMPERATURE_H1_ERROR";
     exampleResponseLibrary->addResponse("H1 Error",eBlocks,builder);

     builder.quadPointField = "AREA";
     exampleResponseLibrary->addResponse("Area",eBlocks,builder);
   }

   // setup closure model
   /////////////////////////////////////////////////////////////

   // Add in the application specific closure model factory
   panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
   Example::ClosureModelFactory_TemplateBuilder cm_builder;
   cm_factory.buildObjects(cm_builder);

   Teuchos::ParameterList closure_models("Closure Models");
   {
     // the exact solution is different depending on rectangular or annular domain
     closure_models.set<bool>("Curvilinear",curvilinear);
     // compute area (since we are in 2D)
     closure_models.sublist("solid").sublist("AREA").set<double>("Value",1.0);
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

   Teuchos::RCP<panzer::FieldManagerBuilder> fmb =
         Teuchos::rcp(new panzer::FieldManagerBuilder);
   fmb->setWorksetContainer(wkstContainer);
   fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
   fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,
                             *linObjFactory,user_data);

   fmb->writeVolumeGraphvizDependencyFiles("Poisson", physicsBlocks);

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

      exampleResponseLibrary->buildResponseEvaluators(physicsBlocks,
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

   panzer::AssemblyEngineInArgs input(ghostCont,container);
   input.alpha = 0;
   input.beta = 1;

   // evaluate physics: This does both the Jacobian and residual at once
   ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);

   // solve linear system
   /////////////////////////////////////////////////////////////

   Teuchos::ParameterList belosList;
   belosList.set("Num Blocks", 1000);              // Maximum number of blocks in Krylov factorization
   belosList.set("Block Size", 1);                 // Blocksize to be used by iterative solver
   belosList.set("Maximum Iterations", 1000);      // Maximum number of iterations allowed
   belosList.set("Maximum Restarts", 1);           // Maximum number of restarts allowed
   belosList.set("Convergence Tolerance", 1.0e-9); // Relative convergence tolerance requested
   belosList.set("Timer Label", "Belos Init");     // Label used by the timers in this solver
   const bool verbose = false;
   if (verbose) {
     belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
                    Belos::TimingDetails + Belos::StatusTestDetails );
     const int frequency = 1;
     if (frequency > 0)
       belosList.set( "Output Frequency", frequency );
   }
   else
     belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

   using ST = double;
   using LO = panzer::LocalOrdinal;
   using GO = panzer::GlobalOrdinal;
   using NT = panzer::TpetraNodeType;
   using OP  = typename Tpetra::Operator<ST,LO,GO,NT>;
   using MV  = typename Tpetra::MultiVector<ST,LO,GO,NT>;
   // using MT = typename Teuchos::ScalarTraits<ST>::magnitudeType;

   using TPLOC = panzer::TpetraLinearObjContainer<ST,LO,GO>;
   auto tp_container = rcp_dynamic_cast<TPLOC>(container);
   Belos::LinearProblem<ST,MV,OP> initProblem(tp_container->get_A(), tp_container->get_x(), tp_container->get_f());
   TEUCHOS_ASSERT(initProblem.setProblem());

   RCP< Belos::SolverManager<ST,MV,OP> > solver
     = rcp(new Belos::PseudoBlockGmresSolMgr<ST,MV,OP>(rcp(&initProblem,false), rcp(&belosList,false)));

   Belos::ReturnType belos_solve_status = solver->solve();
   TEUCHOS_ASSERT(belos_solve_status == Belos::Converged);

   out << "Linear Solver Converged, achieved tol=" << solver->achievedTol() << ", num iters=" << solver->getNumIters() << std::endl;

   // we have now solved for the residual correction from
   // zero in the context of a Newton solve.
   //     J*e = -r = -(f - J*0) where f = J*u
   // Therefore we have  J*e=-J*u which implies e = -u
   // thus we will scale the solution vector
   tp_container->get_x()->scale(-1.0);

   // output data (optional)
   /////////////////////////////////////////////////////////////

   // write out linear system
   if (false) {
     stackedTimer->start("Write Linear System");
     const auto crsOp = *tp_container->get_A();
     Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile("a_op_rowmap.mm",*(crsOp.getRowMap()));
     Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile("a_op_colmap.mm",*(crsOp.getColMap()));
     Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile("a_op_domainmap.mm",*(crsOp.getDomainMap()));
     Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile("a_op_rangemap.mm",*(crsOp.getRangeMap()));
     Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeSparseFile("a_op.mm",crsOp);
     Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeDenseFile("x_vec.mm",*(tp_container->get_x()));
     Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeDenseFile("b_vec.mm",*(tp_container->get_f()));
     stackedTimer->stop("Write Linear System");
   }

   // write out solution to matrix
   if (true) {
     stackedTimer->start("Write Solution to Exodus");
     panzer::AssemblyEngineInArgs respInput(ghostCont,container);
     respInput.alpha = 0;
     respInput.beta = 1;

     stkIOResponseLibrary->addResponsesToInArgs<panzer::Traits::Residual>(respInput);
     stkIOResponseLibrary->evaluate<panzer::Traits::Residual>(respInput);

     // Due to multiple instances of this test being run at the same
     // time (one for each celltype and each order), we need to
     // differentiate output to prevent race conditions on output
     // file. Multiple runs for different mesh refinement levels for
     // the same celltype/order are ok as they are staged one after
     // another in the ADD_ADVANCED_TEST cmake macro.
     std::ostringstream filename;
     if (curvilinear) filename << "annulus_";
     filename << "output_" << celltype << "_p" << basis_order << ".exo";
     mesh->writeToExodus(filename.str());
     stackedTimer->stop("Write Solution to Exodus");
   }

   // compute error norm
   /////////////////////////////////////////////////////////////

   if (true) {
      stackedTimer->start("Compute Responses");
      Teuchos::FancyOStream lout(Teuchos::rcpFromRef(std::cout));
      lout.setOutputToRootOnly(0);

      panzer::AssemblyEngineInArgs respInput(ghostCont,container);
      respInput.alpha = 0;
      respInput.beta = 1;

      Teuchos::RCP<panzer::ResponseBase> area_resp = exampleResponseLibrary->getResponse<panzer::Traits::Residual>("Area");
      Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > area_resp_func =
             Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(area_resp);
      Teuchos::RCP<Thyra::VectorBase<double> > area_respVec = Thyra::createMember(area_resp_func->getVectorSpace());
      area_resp_func->setVector(area_respVec);

      Teuchos::RCP<panzer::ResponseBase> l2_resp = exampleResponseLibrary->getResponse<panzer::Traits::Residual>("L2 Error");
      Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > l2_resp_func =
             Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(l2_resp);
      Teuchos::RCP<Thyra::VectorBase<double> > l2_respVec = Thyra::createMember(l2_resp_func->getVectorSpace());
      l2_resp_func->setVector(l2_respVec);

      Teuchos::RCP<panzer::ResponseBase> h1_resp = exampleResponseLibrary->getResponse<panzer::Traits::Residual>("H1 Error");
      Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > h1_resp_func =
             Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(h1_resp);
      Teuchos::RCP<Thyra::VectorBase<double> > h1_respVec = Thyra::createMember(h1_resp_func->getVectorSpace());
      h1_resp_func->setVector(h1_respVec);

      exampleResponseLibrary->addResponsesToInArgs<panzer::Traits::Residual>(respInput);
      exampleResponseLibrary->evaluate<panzer::Traits::Residual>(respInput);

      double area_exact = 1.;
      if (curvilinear)
        area_exact = M_PI * 1.0 * 1.0 - M_PI * .5 * .5;

      lout << "This is the Basis Order" << std::endl;
      lout << "Basis Order = " << basis_order << std::endl;
      lout << "This is the L2 Error" << std::endl;
      lout << "L2 Error = " << sqrt(l2_resp_func->value) << std::endl;
      lout << "This is the H1 Error" << std::endl;
      lout << "H1 Error = " << sqrt(h1_resp_func->value) << std::endl;
      lout << "This is the error in area" << std::endl;
      lout << "Area Error = " << std::abs(area_resp_func->value - area_exact) << std::endl;
      stackedTimer->stop("Compute Responses");
   }

   stackedTimer->stop("Mixed Poisson");
   stackedTimer->stopBaseTimer();
   if (true) {
     std::ostringstream filename;
     filename << "timing_" << problem_name << "_" << celltype << "_p" << basis_order
              << "_" << x_elements << "x"<< y_elements << ".log";

     std::fstream timing_stream(filename.str().c_str(),std::fstream::out|std::fstream::trunc);
     Teuchos::StackedTimer::OutputOptions options;
     options.output_fraction = true;
     options.output_minmax = false;
     options.output_histogram = false;
     options.num_histogram = 5;
     stackedTimer->report(timing_stream, Teuchos::DefaultComm<int>::getComm(), options);
   }

   // all done!
   /////////////////////////////////////////////////////////////

   out << "ALL PASSED" << std::endl;

   return 0;
}

void testInitialization(const int basis_order,
                        const Teuchos::RCP<Teuchos::ParameterList>& ipb,
                        std::vector<panzer::BC>& bcs,
                        const bool curvilinear)
{
  {
    const int integration_order = 10;
    Teuchos::ParameterList& p = ipb->sublist("Poisson Physics");
    p.set("Type","Poisson");
    p.set("Model ID","solid");
    p.set("Basis Type","HGrad");
    p.set("Basis Order",basis_order);
    p.set("Integration Order",integration_order);
  }

  if (curvilinear) { // Boundary conditions for annulus are T = 5 on inner circle, T = 1 on outer
   {
      std::size_t bc_id = 0;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "inner";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name,
  		    strategy, p);
      bcs.push_back(bc);
   }

   {
      std::size_t bc_id = 1;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "outer";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 1.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name,
  		    strategy, p);
      bcs.push_back(bc);
   }
  } else { // Boundary conditions for rectangular domain are 0 Dirichlet everywhere
   {
      std::size_t bc_id = 0;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 0.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name,
  		    strategy, p);
      bcs.push_back(bc);
   }

   {
      std::size_t bc_id = 1;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 0.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name,
  		    strategy, p);
      bcs.push_back(bc);
   }

   {
      std::size_t bc_id = 2;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 0.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name,
  		    strategy, p);
      bcs.push_back(bc);
   }

   {
      std::size_t bc_id = 3;
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = "bottom";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 0.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name,
  		    strategy, p);
      bcs.push_back(bc);
   }
  }
}
