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

#include "Panzer_NodeType.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_InitialCondition_Builder.hpp"
#include "Panzer_CheckBCConsistency.hpp"

#include "Panzer_STK_CubeHexMeshFactory.hpp"
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
//#include "MiniEM_RefMaxwellPreconditionerFactory.hpp"
#include "MiniEM_DiscreteGradient.hpp"

#include <string>
#include <iostream>

Teuchos::RCP<Teuchos::ParameterList> maxwellParameterList(const int basis_order);
std::vector<panzer::BC> homogeneousBoundaries();
std::vector<panzer::BC> auxiliaryBoundaries();
Teuchos::RCP<Teuchos::ParameterList> auxOpsParameterList(const int basis_order);
Teuchos::RCP<Teuchos::ParameterList> maxwellSolverParameterList(const bool use_ilu, const bool use_refmaxwell, const bool print_diagnostics);
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
void writeToExodus(double time_stamp,
    const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
    const panzer::ModelEvaluator<double> & model,
    panzer::ResponseLibrary<panzer::Traits> & stkIOResponseLibrary,
    panzer_stk::STK_Interface & mesh);

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
 *   use-ilu: use ILU as smoother for E-field Schur complement
 *
 * ******************************************************************************
 */

int main(int argc,char * argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
  Kokkos::initialize(argc, argv);

  { // Put main inside braces to delete all views prior to finalize of device

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
    if (mpiSession.getNProc() > 1) {
      out->setShowProcRank(true);
      out->setOutputToRootOnly(0);
    }

    Teuchos::RCP<const Teuchos::MpiComm<int> > comm
      = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

    {
      Teuchos::RCP<Teuchos::TimeMonitor> tM = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: Total Time"))));
  
      // defaults for command-line options
      int x_elements=10,y_elements=10,z_elements=10,x_procs=-1,y_procs=1,z_procs=1,basis_order=1;
      double cfl=4.0;
      std::size_t workset_size = 20;
      bool exodus_output = false;
      bool matrix_output = false;
      bool use_ilu = false;
      bool use_refmaxwell = false;
      bool print_diagnostics = false;
      Teuchos::CommandLineProcessor clp;
      clp.setOption("x-elements",&x_elements);
      clp.setOption("y-elements",&y_elements);
      clp.setOption("z-elements",&z_elements);
      clp.setOption("x-procs",&x_procs);
      clp.setOption("y-procs",&y_procs);
      clp.setOption("z-procs",&z_procs);
      clp.setOption("basis-order",&basis_order);
      clp.setOption("cfl",&cfl);
      clp.setOption("workset-size",&workset_size);
      clp.setOption("exodus-output","no-exodus-output",&exodus_output);
      clp.setOption("matrix-output","no-matrix-output",&matrix_output);
      clp.setOption("use-ilu","no-ilu",&use_ilu);
      clp.setOption("use-refmaxwell","use-augmentation",&use_refmaxwell);
      clp.setOption("subsolve-diagnostics","no-subsolve-diagnostics",&print_diagnostics);
  
      // parse command-line argument
      TEUCHOS_ASSERT(clp.parse(argc,argv)==Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL);
  
      // compute dt from cfl
      double c  = std::sqrt(1.0/8.854187817e-12/1.2566370614e-6);
      double min_dx = 1.0/std::max(x_elements,std::max(y_elements,z_elements));
      double dt = cfl*min_dx/c;
  
      // set mesh factory parameters
      panzer_stk::CubeHexMeshFactory mesh_factory;
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Blocks",1);
      pl->set("Y Blocks",1);
      pl->set("Z Blocks",1);
      pl->set("X Elements",x_elements);
      pl->set("Y Elements",y_elements);
      pl->set("Z Elements",z_elements);
      pl->set("X Procs",x_procs);
      pl->set("Y Procs",y_procs);
      pl->set("Z Procs",z_procs);
  
      // periodic boundaries
      //      Teuchos::ParameterList& per_pl = pl->sublist("Periodic BCs");
      //      per_pl.set("Count", 3);
      //      per_pl.set("Periodic Condition 1", "xy-all 1e-8: front;back");
      //      per_pl.set("Periodic Condition 2", "xz-all 1e-8: top;bottom");
      //      per_pl.set("Periodic Condition 3", "yz-all 1e-8: left;right");
  
      mesh_factory.setParameterList(pl);
  
      // build mesh
      RCP<panzer_stk::STK_Interface> mesh = mesh_factory.buildUncommitedMesh(MPI_COMM_WORLD);
  
      // data container for auxiliary linear operators used in preconditioning (mass matrix and gradient)
      Teuchos::RCP<panzer::GlobalEvaluationDataContainer> auxGlobalData = Teuchos::rcp(new panzer::GlobalEvaluationDataContainer);
  
      // factory definitions
      Teuchos::RCP<mini_em::EquationSetFactory> eqset_factory = 
        Teuchos::rcp(new mini_em::EquationSetFactory(auxGlobalData)); // where the maxwell equations are defined
      mini_em::BCFactory bc_factory;                                  // where boundary conditions are defined 
  
      // GobalData sets ostream and parameter interface to physics
      Teuchos::RCP<panzer::GlobalData> globalData = panzer::createGlobalData();
  
      // define physics block parameter list and boundary conditions
      Teuchos::RCP<Teuchos::ParameterList> physicsBlock_pl = maxwellParameterList(basis_order);
      std::vector<panzer::BC> bcs = homogeneousBoundaries();
      std::vector<panzer::BC> aux_bcs;// = auxiliaryBoundaries();
  
      // build the physics blocks objects
      std::vector<RCP<panzer::PhysicsBlock> > physicsBlocks;
      {
        bool build_transient_support = true;
  
        const panzer::CellData volume_cell_data(workset_size, mesh->getCellTopology("eblock-0_0_0"));
  
        // Can be overridden by the equation set
        int default_integration_order = 2;
        
        // the physics block knows how to build and register evaluator with the field manager
        RCP<panzer::PhysicsBlock> pb 
  	= rcp(new panzer::PhysicsBlock(physicsBlock_pl,
  				       "eblock-0_0_0", 
  				       default_integration_order,
  				       volume_cell_data,
  				       eqset_factory,
  				       globalData,
  				       build_transient_support));
  
        // we can have more than one physics block, one per element block
        physicsBlocks.push_back(pb);
      }
  
      // build the auxiliary physics blocks objects
      Teuchos::RCP<Teuchos::ParameterList> auxPhysicsBlock_pl = auxOpsParameterList(basis_order);
      std::vector<RCP<panzer::PhysicsBlock> > auxPhysicsBlocks;
      {
        bool build_transient_support = false;
  
        const panzer::CellData volume_cell_data(workset_size, mesh->getCellTopology("eblock-0_0_0"));
  
        // Can be overridden by the equation set
        int default_integration_order = 2;
        
        // the physics block knows how to build and register evaluator with the field manager
        RCP<panzer::PhysicsBlock> pb 
  	= rcp(new panzer::PhysicsBlock(auxPhysicsBlock_pl,
  				       "eblock-0_0_0", 
  				       default_integration_order,
  				       volume_cell_data,
  				       eqset_factory,
  				       globalData,
  				       build_transient_support));
  
        // we can have more than one physics block, one per element block
        auxPhysicsBlocks.push_back(pb);
      }
  
      // Add fields to the mesh data base (this is a peculiarity of how STK classic requires the
          // fields to be setup)
      createExodusFile(physicsBlocks, Teuchos::rcpFromRef(mesh_factory), mesh, exodus_output);
  
      // build worksets
      Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
         = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
      Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
         = Teuchos::rcp(new panzer::WorksetContainer);
      Teuchos::RCP<panzer::WorksetContainer> auxWkstContainer  // attach it to a workset container (uses lazy evaluation)
         = Teuchos::rcp(new panzer::WorksetContainer);
      wkstContainer->setFactory(wkstFactory);
      auxWkstContainer->setFactory(wkstFactory);
      for(size_t i=0;i<physicsBlocks.size();i++) 
      {
        wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
        auxWkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),auxPhysicsBlocks[i]->getWorksetNeeds());
      }
      wkstContainer->setWorksetSize(workset_size);
      auxWkstContainer->setWorksetSize(workset_size);
  
      // build DOF Managers and linear object factories
   
      // build the connection manager 
      const Teuchos::RCP<panzer::ConnManager<int,panzer::Ordinal64> > 
        conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager<panzer::Ordinal64>(mesh));
  
      // blocked degree of freedom manager
      panzer::BlockedDOFManagerFactory<int,panzer::Ordinal64> globalIndexerFactory;
      std::string fieldOrder = "blocked: B_face E_edge";
      RCP<panzer::UniqueGlobalIndexerBase > dofManager = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager,fieldOrder);
  
      // auxiliary dof manager
      std::string auxFieldOrder = "blocked: AUXILIARY_NODE AUXILIARY_EDGE";
      RCP<panzer::UniqueGlobalIndexerBase > auxDofManager = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),auxPhysicsBlocks,conn_manager,auxFieldOrder);
  
      // construct some linear algebra objects, build object to pass to evaluators
      Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
           = Teuchos::rcp(new panzer::BlockedTpetraLinearObjFactory<panzer::Traits, double, int, panzer::Ordinal64>(comm,rcp_dynamic_cast<panzer::BlockedDOFManager<int,panzer::Ordinal64> >(dofManager,true)));
      Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > auxLinObjFactory
           = Teuchos::rcp(new panzer::BlockedTpetraLinearObjFactory<panzer::Traits, double, int, panzer::Ordinal64>(comm,rcp_dynamic_cast<panzer::BlockedDOFManager<int,panzer::Ordinal64> >(auxDofManager,true)));
  
      // Assign the dof managers to worksets
      wkstContainer->setGlobalIndexer(dofManager);
      auxWkstContainer->setGlobalIndexer(auxDofManager);
  
      // setup closure model
      panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory; 
      mini_em::ClosureModelFactory_TemplateBuilder cm_builder;
      cm_factory.buildObjects(cm_builder);
      Teuchos::ParameterList closure_models("Closure Models");
      {
        closure_models.sublist("electromagnetics").sublist("CURRENT").set<std::string>("Type","GAUSSIAN PULSE"); // a gaussian current source
        closure_models.sublist("electromagnetics").sublist("CURRENT").set<double>("dt",dt); // set pulse width such that dt resolves it
      }
  
      Teuchos::ParameterList user_data("User Data"); // user data can be empty here
  
      // add full maxwell solver to teko
      RCP<Teko::Cloneable> clone = rcp(new Teko::AutoClone<mini_em::FullMaxwellPreconditionerFactory>());
      Teko::PreconditionerFactory::addPreconditionerFactory("Full Maxwell Preconditioner",clone);
      // add refMaxwell solver to teko
      //clone = rcp(new Teko::AutoClone<mini_em::RefMaxwellPreconditionerFactory>());
      //Teko::PreconditionerFactory::addPreconditionerFactory("RefMaxwell Preconditioner",clone);
  
      // add callbacks to request handler. these are for requesting auxiliary operators and for providing
      // coordinate information to MueLu
      Teuchos::RCP<Teko::RequestHandler> req_handler = Teuchos::rcp(new Teko::RequestHandler());
      req_handler->addRequestCallback(Teuchos::rcp(new mini_em::OperatorRequestCallback(auxGlobalData, matrix_output)));
      Teuchos::RCP<panzer_stk::ParameterListCallbackBlocked<int,panzer::Ordinal64> > callback =
                  rcp(new panzer_stk::ParameterListCallbackBlocked<int,panzer::Ordinal64>(
                             rcp_dynamic_cast<panzer_stk::STKConnManager<panzer::Ordinal64> >(conn_manager,true),
                             rcp_dynamic_cast<panzer::BlockedDOFManager<int,panzer::Ordinal64> >(dofManager,true),
                             rcp_dynamic_cast<panzer::BlockedDOFManager<int,panzer::Ordinal64> >(auxDofManager,true)));
      req_handler->addRequestCallback(callback);

      // add discrete gradient
      addDiscreteGradientToRequestHandler(auxLinObjFactory,req_handler);
  
      // build linear solver
      RCP<Teuchos::ParameterList> lin_solver_pl = maxwellSolverParameterList(use_ilu,use_refmaxwell,print_diagnostics);
      RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory
      = panzer_stk::buildLOWSFactory(true, dofManager, conn_manager,
          Teuchos::as<int>(mesh->getDimension()),
          comm, lin_solver_pl,req_handler, false, false, auxDofManager);
  
      //setup model evaluators
      RCP<panzer::ModelEvaluator<double> > physics = rcp(new panzer::ModelEvaluator<double> (linObjFactory, lowsFactory, globalData, true, 0.0));
      RCP<panzer::ModelEvaluator<double> > auxPhysics = rcp(new panzer::ModelEvaluator<double> (auxLinObjFactory, lowsFactory, globalData, false, 0.0));
  
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
  
      Thyra::ModelEvaluatorBase::InArgs<double> auxInArgs = auxPhysics->getNominalValues();
      Thyra::ModelEvaluatorBase::OutArgs<double> auxOutArgs = auxPhysics->createOutArgs();
      Teuchos::RCP<Thyra::LinearOpBase<double> > aux_W_op = auxPhysics->create_W_op();
      auxOutArgs.set_W_op(aux_W_op);
      auxPhysics->evalModel(auxInArgs, auxOutArgs);
  
      // setup a response library to write to the mesh
      RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
      = buildSTKIOResponseLibrary(physicsBlocks,linObjFactory,wkstContainer,dofManager,cm_factory,mesh,
          closure_models);
  
      // set up the solution vector, jacobian, and residual
      RCP<Thyra::VectorBase<double> > solution_vec = Thyra::createMember(physics->get_x_space());
      Thyra::assign(solution_vec.ptr(),0.0);
      RCP<Thyra::LinearOpWithSolveBase<double> > jacobian = physics->create_W();
  
      RCP<Thyra::VectorBase<double> > residual = Thyra::createMember(physics->get_f_space());
  
      // set up the model evaluator
      Thyra::ModelEvaluatorBase::InArgs<double> inArgs = physics->createInArgs();
      inArgs.set_alpha(1.0/dt);
      inArgs.set_beta(1.0);
  
      // initial condition is zero, define x_dot accordingly
      RCP<const Thyra::VectorBase<double> > x = inArgs.get_x();
      RCP<Thyra::VectorBase<double> > x_dot = Thyra::createMember(physics->get_x_space());
      Thyra::V_StV(x_dot.ptr(),1.0/dt,*x);
      inArgs.set_x_dot(x_dot);
      Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = physics->createOutArgs();
      outArgs.set_f(residual);
      outArgs.set_W(jacobian);
  
      // compute the jacobian matrix only once
      Kokkos::fence();
      physics->evalModel(inArgs,outArgs);
      outArgs.set_W(RCP<Thyra::LinearOpWithSolveBase<double> >(NULL));
  
      // take 20 time-steps with Backward Euler
      if (exodus_output)
        writeToExodus(0,solution_vec,*physics,*stkIOResponseLibrary,*mesh);

      {
        Teuchos::RCP<Teuchos::TimeMonitor> tM = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: timestepper"))));
        for(int ts = 1; ts < 21; ts++)
        {
          RCP<Thyra::VectorBase<double> > x_old = solution_vec->clone_v();
    
          inArgs.set_t(dt*ts);
          Thyra::V_StVpStV(x_dot.ptr(),1.0/dt,*x,-1.0/dt,*x_old);
          inArgs.set_x_dot(x_dot);
    
          // construct the residual
          physics->evalModel(inArgs,outArgs);
    
          // solve
          jacobian->solve(Thyra::NOTRANS,*residual,solution_vec.ptr());
          Thyra::V_StVpStV(solution_vec.ptr(),1.0,*x_old,-1.0,*solution_vec);
    
          // write to an exodus file
          if (exodus_output)
          {
            Teuchos::RCP<Teuchos::TimeMonitor> tM = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: timestepper: writeToExodus"))));
            writeToExodus(dt*ts,solution_vec,*physics,*stkIOResponseLibrary,*mesh);
          }
    
          (*out) << "finished time step " << ts << std::endl;
        }
      }
    }
    Teuchos::TimeMonitor::summarize(*out,false,true,false,Teuchos::Union);

  }
  Kokkos::finalize();
  return 0;
}

//! Create a parameter list defining the Maxwell equations physics block
Teuchos::RCP<Teuchos::ParameterList> maxwellParameterList(const int basis_order)
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList("Physics Block"));
  const int integration_order = 2*basis_order;

  Teuchos::ParameterList& p = pl->sublist("Maxwell Physics");
  p.set("Type","Maxwell");
  p.set("Model ID","electromagnetics");
  p.set("Basis Order",basis_order);
  p.set("Integration Order",integration_order);
  //p.set("Epsilon",1.0);
  //p.set("Mu",1.0);

  return pl;
}

//! Create BCs for E x n = 0 and B . n = 0 on all boundaries
std::vector<panzer::BC> homogeneousBoundaries()
{
  std::vector<panzer::BC> bcs;

  std::string sidesets[6] = {"left","right","top","bottom","back","front"};
  std::string dofs[2]     = {"E_edge","B_face"};

  std::size_t bc_id = 0;
  for (int s = 0; s < 6; s++)
    for (int d = 0; d < 0; d++)
    {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = sidesets[s];
      std::string element_block_id = "eblock-0_0_0";
      std::string dof_name = dofs[d];
      std::string strategy = "Constant";
      Teuchos::ParameterList p;
      p.set("Value X",0.0);
      p.set("Value Y",0.0);
      p.set("Value Z",0.0);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
                    strategy, p);
      bcs.push_back(bc);
      bc_id++;
    }
 
  return bcs;
}

//! Create BCs for auxiliary operators
std::vector<panzer::BC> auxiliaryBoundaries()
{
  std::vector<panzer::BC> bcs;

  std::string sidesets[6] = {"left","right","top","bottom","back","front"};
  std::string eq_sets[2]  = {"Mass Matrix AUXILIARY_NODE","Weak Gradient"};
  std::string dofs[2]     = {"AUXILIARY_NODE","AUXILIARY_EDGE"};

  std::size_t bc_id = 0;
  for (int s = 0; s < 6; s++)
    for (int d = 0; d < 2; d++)
    {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = sidesets[s];
      std::string element_block_id = "eblock-0_0_0";
      std::string dof_name = eq_sets[d];
      std::string strategy = "AuxConstant";
      Teuchos::ParameterList p;
      p.set("Field Name",dofs[d]);
      p.set("Value",0.0);
      panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
                    strategy, p);
      bcs.push_back(bc);
      bc_id++;
    }
 
  return bcs;
}

//! Create parameter list defining nodal mass matrix and node-edge weak gradient
Teuchos::RCP<Teuchos::ParameterList> auxOpsParameterList(const int basis_order)
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList("Aux Physics Block"));
  const int integration_order = 2*basis_order;

  {
    Teuchos::ParameterList& p = pl->sublist("Auxiliary Mass Physics");
    p.set("Type","Auxiliary Mass Matrix");
    p.set("DOF Name","AUXILIARY_NODE");
    p.set("Basis Type","HGrad");
    p.set("Basis Order",basis_order);
    p.set("Integration Order",integration_order);
  }

  {
    Teuchos::ParameterList& p = pl->sublist("Auxiliary Gradient Physics");
    p.set("Type","Auxiliary Weak Gradient");
    p.set("Vector Name","AUXILIARY_EDGE");
    p.set("Scalar Name","AUXILIARY_NODE");
    p.set("Basis Order",basis_order);
    p.set("Integration Order",integration_order);
  }

  return pl;
}

//! Create parameter list defining linear solver
Teuchos::RCP<Teuchos::ParameterList> maxwellSolverParameterList(const bool use_ilu, const bool use_refmaxwell, const bool print_diagnostics)
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList("Linear Solver"));
  pl->set("Linear Solver Type","Belos");
  {
    Teuchos::ParameterList& ls_types = pl->sublist("Linear Solver Types");
    {
      Teuchos::ParameterList& belos = ls_types.sublist("Belos");
      belos.set("Solver Type","Pseudo Block GMRES");
      {
        Teuchos::ParameterList& solver_types = belos.sublist("Solver Types");
        {
          Teuchos::ParameterList& gmres = solver_types.sublist("Pseudo Block GMRES");
          gmres.set("Convergence Tolerance", 1.0e-5);
          gmres.set("Output Frequency",      1);
          gmres.set("Output Style",          1);
          gmres.set("Verbosity",             33);
          gmres.set("Maximum Iterations",    100);
          gmres.set("Block Size",            1);
          gmres.set("Num Blocks",            100);
        }
      }
      {
        Teuchos::ParameterList& verb = belos.sublist("VerboseObject");
        verb.set("Output File", "none");
        verb.set("Verbosity Level", "medium");
      }
    }
  }
  pl->set("Preconditioner Type", "Teko");
  {
    Teuchos::ParameterList& prec_types = pl->sublist("Preconditioner Types");
    {
      Teuchos::ParameterList& teko = prec_types.sublist("Teko");
      teko.set("Inverse Type", "Maxwell");
      {
        Teuchos::ParameterList& inv_fac = teko.sublist("Inverse Factory Library");
        { // Maxwell block preconditioner settings
          Teuchos::ParameterList& maxwell = inv_fac.sublist("Maxwell");
          maxwell.set("Type", "Full Maxwell Preconditioner");
          maxwell.set("Use ILU", use_ilu);
          maxwell.set("Use refMaxwell", use_refmaxwell);
          maxwell.set("Print Diagnostics", print_diagnostics);
        }
      }
    }
  }

  return pl;
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

    Teuchos::ParameterList output_pl("Output");
    Teuchos::ParameterList& cell_avg_q = output_pl.sublist("Cell Average Quantities");
    Teuchos::ParameterList& cell_avg_v = output_pl.sublist("Cell Average Vectors");
    cell_avg_v.set("eblock-0_0_0","CURRENT");
    Teuchos::ParameterList& cell_q = output_pl.sublist("Cell Quantities");
    Teuchos::ParameterList& nodal_q = output_pl.sublist("Nodal Quantities");
    Teuchos::ParameterList& a_nodal_q = output_pl.sublist("Allocate Nodal Quantities");
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

  // this automatically adds in the nodal fields
  Teuchos::ParameterList output_pl("Output");
  Teuchos::ParameterList& cell_avg_q = output_pl.sublist("Cell Average Quantities");
  Teuchos::ParameterList& cell_avg_v = output_pl.sublist("Cell Average Vectors");
  cell_avg_v.set("eblock-0_0_0","CURRENT");
  Teuchos::ParameterList& cell_q = output_pl.sublist("Cell Quantities");
  Teuchos::ParameterList& nodal_q = output_pl.sublist("Nodal Quantities");
  Teuchos::ParameterList& a_nodal_q = output_pl.sublist("Allocate Nodal Quantities");
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

void writeToExodus(double time_stamp,
    const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
    const panzer::ModelEvaluator<double> & model,
    panzer::ResponseLibrary<panzer::Traits> & stkIOResponseLibrary,
    panzer_stk::STK_Interface & mesh)
{
  // fill STK mesh objects
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = model.createInArgs();
  inArgs.set_x(x);

  panzer::AssemblyEngineInArgs respInput;
  model.setupAssemblyInArgs(inArgs,respInput);

  stkIOResponseLibrary.addResponsesToInArgs<panzer::Traits::Residual>(respInput);
  stkIOResponseLibrary.evaluate<panzer::Traits::Residual>(respInput);

  mesh.writeToExodus(time_stamp);
}

