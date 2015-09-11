// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Phalanx_KokkosUtilities.hpp"

#include "Panzer_config.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_GlobalData.hpp"

#include "Panzer_STK_config.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_Utilities.hpp"

#include "Epetra_MpiComm.h"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"

#include "Example_BCStrategy_Factory.hpp"
#include "Example_ClosureModel_Factory_TemplateBuilder.hpp"
#include "Example_EquationSetFactory.hpp"

#include "AztecOO.h"

using Teuchos::RCP;
using Teuchos::rcp;

static void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
                               std::vector<panzer::BC>& bcs);

static void solve_Ax_eq_b(Epetra_CrsMatrix& A, Epetra_Vector& b, Epetra_Vector& x)
{
  // Setup the linear solve: notice A is used directly 
  Epetra_LinearProblem problem(&A, &x, &b);

  // build the solver
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_solver,AZ_gmres); // we don't push out dirichlet conditions
  solver.SetAztecOption(AZ_precond,AZ_none);
  solver.SetAztecOption(AZ_kspace,300);
  solver.SetAztecOption(AZ_output,0);
  solver.SetAztecOption(AZ_precond,AZ_Jacobi);

  // solve the linear system
  solver.Iterate(1000,1e-5);

  // we have now solved for the residual correction from
  // zero in the context of a Newton solve.
  //     J*e = -r = -(f - J*0) where f = J*u
  // Therefore we have  J*e=-J*u which implies e = -u
  // thus we will scale the solution vector 
  x.Scale(-1.0);
}

static void
assembleAndSolve(panzer::AssemblyEngine_TemplateManager<panzer::Traits>& ae_tm,
                 const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& linObjFactory,
                 Teuchos::RCP<panzer::EpetraLinearObjContainer>& ep_container,
                 Teuchos::RCP<panzer::LinearObjContainer>& ghostCont,
                 const double rtol = 1e-10)
{
  double bnorm;
  Teuchos::RCP<Epetra_Vector> x, dx;
  Teuchos::RCP<Epetra_CrsMatrix> A;
  for (int it = 0; it < 500; ++it) {
    // assemble linear system
    // build linear algebra objects: Ghost is for parallel assembly, it contains
    //                               local element contributions summed, the global IDs
    //                               are not unique. The non-ghosted or "global"
    //                               container will contain the sum over all processors
    //                               of the ghosted objects. The global indices are unique.
    ghostCont = linObjFactory->buildGhostedLinearObjContainer();
    const Teuchos::RCP<panzer::LinearObjContainer> container = linObjFactory->buildLinearObjContainer();
    linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                              panzer::LinearObjContainer::F |
                                              panzer::LinearObjContainer::Mat, *ghostCont);
    linObjFactory->initializeContainer(panzer::LinearObjContainer::X |
                                       panzer::LinearObjContainer::F |
                                       panzer::LinearObjContainer::Mat, *container);
    ghostCont->initialize();
    container->initialize();

    // Convert generic linear object container to epetra container.
    ep_container = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(container);
    if (x.is_null()) {
      x = ep_container->get_x();
      dx = Teuchos::rcp(new Epetra_Vector(ep_container->get_x()->Map()));
    } else
      ep_container->set_x(x);

    panzer::AssemblyEngineInArgs input(ghostCont, ep_container);
    input.alpha = 0;
    input.beta = 1;

    // Evaluate physics. Evaluate the Residual by itself because the Jacobian
    // gather/scatter evaluators don't work yet.
    ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
    Teuchos::RCP<Epetra_Vector> r;
    if (it == 0) {
      r = Teuchos::rcp(new Epetra_Vector(*ep_container->get_f()));
      // Get the inexact Jacobian. This is a linear problem, so need to evaluate
      // it just once.
      ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);
      {
        A = Teuchos::rcp(new Epetra_CrsMatrix(*ep_container->get_A()));
        // The Jacobian is rank deficient at the moment. It's wrong, anyway. For
        // now, just add I to the diag.
        Epetra_Vector diag(r->Map());
        A->ExtractDiagonalCopy(diag);
        for (int i = 0; i < diag.MyLength(); ++i)
          diag[i] += 1;
        A->ReplaceDiagonalValues(diag);
        //EpetraExt::RowMatrixToMatrixMarketFile("A.mm", *A);
      }
      // Set f to its correct value. Some side effect makes evaluating Jacobian
      // then Residual wrong.
      ep_container->set_f(r);
    } else
      r = ep_container->get_f();

    double rnorm;
    r->Norm2(&rnorm);
    if (it == 0) bnorm = rnorm;
    if (it % 10 == 0) printf("it %2d rnorm %1.3e\n", it, rnorm);
    if (rnorm <= rtol*bnorm) break;

    solve_Ax_eq_b(*A, *r, *dx);

    x->Update(1, *dx, 1);
  }
}

// calls MPI_Init and MPI_Finalize
int main(int argc,char * argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using panzer::StrPureBasisPair;
  using panzer::StrPureBasisComp;

  PHX::InitializeKokkosDevice();

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  out.setShowProcRank(true);

  // variable declarations
  ////////////////////////////////////////////////////

  // factory definitions
  Teuchos::RCP<Example::EquationSetFactory> eqset_factory = 
    Teuchos::rcp(new Example::EquationSetFactory); // where poison equation is defined
  Example::BCStrategyFactory bc_factory;    // where boundary conditions are defined 

  panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

  // other declarations
  const std::size_t workset_size = 20;

  // construction of uncommitted (no elements) mesh 
  ////////////////////////////////////////////////////////

  // set mesh factory parameters
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  pl->set("X Blocks",2);
  pl->set("Y Blocks",1);
  pl->set("X Elements",10); // per block
  pl->set("Y Elements",10); // per block
  mesh_factory.setParameterList(pl);

  RCP<panzer_stk_classic::STK_Interface> mesh = mesh_factory.buildUncommitedMesh(MPI_COMM_WORLD);

  // construct input physics and physics block
  ////////////////////////////////////////////////////////

  const Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
  std::vector<panzer::BC> bcs;
  std::vector<RCP<panzer::PhysicsBlock> > physicsBlocks;
  {
    testInitialization(ipb, bcs);

    std::map<std::string,std::string> block_ids_to_physics_ids;
    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;

    block_ids_to_physics_ids["eblock-0_0"] = "Poisson Physics Left";
    block_ids_to_physics_ids["eblock-1_0"] = "Poisson Physics Right";

    block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
    block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");
      
    // GobalData sets ostream and parameter interface to physics
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

    // Can be overridden by the equation set
    int default_integration_order = 1;
      
    // the physics block knows how to build and register evaluator with the field manager
    panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                               block_ids_to_cell_topo,
                               ipb,
                               default_integration_order,
                               workset_size,
                               eqset_factory,
                               gd,
                               false,
                               physicsBlocks);
  }


  // finish building mesh, set required field variables and mesh bulk data
  ////////////////////////////////////////////////////////////////////////
  {
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator physIter;
    for(physIter=physicsBlocks.begin();physIter!=physicsBlocks.end();++physIter) {
      Teuchos::RCP<const panzer::PhysicsBlock> pb = *physIter;
      const std::vector<StrPureBasisPair> & blockFields = pb->getProvidedDOFs();

      // insert all fields into a set
      std::set<StrPureBasisPair,StrPureBasisComp> fieldNames;
      fieldNames.insert(blockFields.begin(),blockFields.end());

      // add basis to DOF manager: block specific
      std::set<StrPureBasisPair,StrPureBasisComp>::const_iterator fieldItr;
      for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {

        mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
      }
    }

    mesh_factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  }

  // build worksets
  ////////////////////////////////////////////////////////

  Teuchos::RCP<panzer_stk_classic::WorksetFactory> wkstFactory
    = Teuchos::rcp(new panzer_stk_classic::WorksetFactory(mesh)); // build STK workset factory
  Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
    = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,physicsBlocks,workset_size));

  std::vector<std::string> elementBlockNames;
  mesh->getElementBlockNames(elementBlockNames);
  std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > volume_worksets;
  panzer::getVolumeWorksetsFromContainer(*wkstContainer,elementBlockNames,volume_worksets);

  // build DOF Manager and linear object factory
  /////////////////////////////////////////////////////////////
 
  // build the connection manager 
  const Teuchos::RCP<panzer::ConnManager<int,int> > 
    conn_manager = Teuchos::rcp(new panzer_stk_classic::STKConnManager<int>(mesh));

  panzer::DOFManagerFactory<int,int> globalIndexerFactory;
  RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
    = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

  // construct some linear algebra object, build object to pass to evaluators
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
    = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(Comm.getConst(),dofManager));

  // setup closure model
  /////////////////////////////////////////////////////////////
 
  // Add in the application specific closure model factory
  panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory; 
  Example::ClosureModelFactory_TemplateBuilder cm_builder;
  cm_factory.buildObjects(cm_builder);

  Teuchos::ParameterList closure_models("Closure Models");
  closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE1").set<double>("Value",0.0); // a constant source
  closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE2").set<double>("Value",0.0); // a constant source
  // SOURCE_TEMPERATURE field is required by the PoissonEquationSet

  Teuchos::ParameterList user_data("User Data"); // user data can be empty here

  // setup field manager builder
  /////////////////////////////////////////////////////////////

  Teuchos::RCP<panzer::FieldManagerBuilder> fmb = 
    Teuchos::rcp(new panzer::FieldManagerBuilder);
  fmb->setWorksetContainer(wkstContainer);
  fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
  fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,
                            *linObjFactory,user_data);
  fmb->writeVolumeGraphvizDependencyFiles("volume", physicsBlocks);
  fmb->writeBCGraphvizDependencyFiles("bc");

  // setup assembly engine
  /////////////////////////////////////////////////////////////

  // build assembly engine: The key piece that brings together everything and 
  //                        drives and controls the assembly process. Just add
  //                        matrices and vectors
  panzer::AssemblyEngine_TemplateManager<panzer::Traits> ae_tm;
  panzer::AssemblyEngine_TemplateBuilder builder(fmb,linObjFactory);
  ae_tm.buildObjects(builder);

  // assemble and solve
  /////////////////////////////////////////////////////////////
  Teuchos::RCP<panzer::EpetraLinearObjContainer> ep_container;
  Teuchos::RCP<panzer::LinearObjContainer> ghostCont;
  assembleAndSolve(ae_tm, linObjFactory, ep_container, ghostCont);
  
  // output data (optional)
  /////////////////////////////////////////////////////////////

  // write out linear system
  if(false) {
    EpetraExt::RowMatrixToMatrixMarketFile("a_op.mm",*ep_container->get_A());
    EpetraExt::VectorToMatrixMarketFile("x_vec.mm",*ep_container->get_x());
    EpetraExt::VectorToMatrixMarketFile("b_vec.mm",*ep_container->get_f());
  }

  // write out solution to matrix
  if(true) {
    // redistribute solution vector to ghosted vector
    linObjFactory->globalToGhostContainer(*ep_container,*ghostCont, panzer::EpetraLinearObjContainer::X 
                                          | panzer::EpetraLinearObjContainer::DxDt); 

    // get X Epetra_Vector from ghosted container
    RCP<panzer::EpetraLinearObjContainer> ep_ghostCont = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ghostCont);
    panzer_stk_classic::write_solution_data(*dofManager,*mesh,*ep_ghostCont->get_x());
    mesh->writeToExodus("output.exo");
  }

  // all done!
  /////////////////////////////////////////////////////////////

  PHX::FinalizeKokkosDevice();

  return 0;
}

void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
                        std::vector<panzer::BC>& bcs)
{
  {
    Teuchos::ParameterList& p = ipb->sublist("Poisson Physics Left").sublist("Equation Set 1");
    p.set("Type","Poisson");
    p.set("Model ID","solid");
    p.set("Basis Type","HGrad");
    p.set("Basis Order",1);
    p.set("Integration Order",2);
    p.set("Suffix","1");
  }
  {
    Teuchos::ParameterList& p = ipb->sublist("Poisson Physics Right").sublist("Equation Set 2");
    p.set("Type","Poisson");
    p.set("Model ID","solid");
    p.set("Basis Type","HGrad");
    p.set("Basis Order",1);
    p.set("Integration Order",2);
    p.set("Suffix","2");
  }
  
  {
    std::size_t bc_id = 0;
    panzer::BCType bctype = panzer::BCT_Dirichlet;
    std::string sideset_id = "left";
    std::string element_block_id = "eblock-0_0";
    std::string dof_name = "TEMPERATURE1";
    std::string strategy = "Constant";
    double value = 0.5;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
                  strategy, p);
    bcs.push_back(bc);
  }
  {
    std::size_t bc_id = 3;
    panzer::BCType bctype = panzer::BCT_Dirichlet;
    std::string sideset_id = "right";
    std::string element_block_id = "eblock-1_0";
    std::string dof_name = "TEMPERATURE2";
    std::string strategy = "Constant";
    double value = -0.5;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
                  strategy, p);
    bcs.push_back(bc);
  }

  const bool doit = true;
  if ( ! doit) {
    std::size_t bc_id = 1;
    panzer::BCType bctype = panzer::BCT_Neumann;
    std::string sideset_id = "vertical_0";
    std::string element_block_id = "eblock-0_0";
    std::string dof_name = "TEMPERATURE1";
    std::string strategy = "Neumann Constant";
    double value = 1.4;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
                  strategy, p);
    bcs.push_back(bc);
  }
  if (true) {
    std::size_t bc_id = 2;
    panzer::BCType bctype = panzer::BCT_Neumann;
    std::string sideset_id = "vertical_0";
    std::string element_block_id = "eblock-1_0";
    std::string dof_name = "TEMPERATURE2";
    std::string strategy = "Neumann Constant";
    double value = 1.4;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, bctype, sideset_id, element_block_id, dof_name, 
                  strategy, p);
    bcs.push_back(bc);
  }
  if ( ! doit) return;
  for (int di = 0; di < 1 /* just NeumannMatch atm */; ++di) {
    std::size_t bc_id = 4 + di;
    Teuchos::ParameterList p;
    p.set("Type", "Interface");
    p.set("Sideset ID", "vertical_0");
    p.set("Element Block ID", di == 0 ? "eblock-0_0" : "eblock-1_0");
    p.set("Equation Set Name", di == 0 ? "TEMPERATURE1" : "TEMPERATURE2");
    p.set("Strategy", di == 0 ? "Neumann Match Interface" : "Weak Dirichlet Match Interface");
    p.set("Element Block ID2", di == 0 ? "eblock-1_0" : "eblock-0_0");
    p.set("Equation Set Name2", di == 0 ? "TEMPERATURE2" : "TEMPERATURE1");
    bcs.push_back(panzer::BC(bc_id, p));
  }
}
