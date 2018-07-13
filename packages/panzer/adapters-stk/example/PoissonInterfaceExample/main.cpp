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
#include <Teuchos_CommandLineProcessor.hpp>

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
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_DOFManager.hpp"
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
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_Utilities.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"

#include "Example_BCStrategy_Factory.hpp"
#include "Example_ClosureModel_Factory_TemplateBuilder.hpp"
#include "Example_EquationSetFactory.hpp"

#include "AztecOO.h"

using Teuchos::RCP;
using Teuchos::rcp;

struct ProblemOptions {
  std::string mesh_filename;
  int nxelem, nyelem;
  bool is3d;

  bool generate_mesh_only;

  std::vector<std::string>
    dof_names, // nominally A1, A2, B, C1, C2
    eb_names,  // two element block names
    ss_names;  // left, interface, right sideset names

  int integration_order;
  // Use a nonlinear version of the Robin interface condition.
  bool nonlinear_Robin;

  // Wrap an outer iteration around the linear solve.
  bool outer_iteration;
  // Test against finite-difference Jacobian.
  bool test_Jacobian;
  // Convergence tolerance.
  double rtol;
  bool check_error;
  Teuchos::RCP<Epetra_CrsMatrix> J;
};

std::ostream& operator<< (std::ostream& os, const ProblemOptions& po) {
  os << "ProblemOptions:\n  DOF Names:";
  for (size_t i = 0; i < po.dof_names.size(); ++i) os << " " << po.dof_names[i];
  os << "\n"
     << "  nxelem " << po.nxelem << "\n"
     << "  nyelem " << po.nyelem << "\n"
     << "  nonlinear_Robin: " << po.nonlinear_Robin << "\n"
     << "  outer_iteration: " << po.outer_iteration << "\n"
     << "  rtol: " << po.rtol << "\n"
     << "  test_Jacobian: " << po.test_Jacobian << "\n";
  return os;
}

std::string strint (const std::string& str, const int i) {
  std::stringstream ss;
  ss << str << i;
  return ss.str();
}

bool solve_Ax_eq_b (Epetra_CrsMatrix& A, Epetra_Vector& b, Epetra_Vector& x, const double rtol)
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
  const int nits = 1000;
  solver.Iterate(nits, rtol);
  return solver.NumIters() < nits;
}

void assembleAndSolve (panzer::AssemblyEngine_TemplateManager<panzer::Traits>& ae_tm,
                       const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& linObjFactory,
                       Teuchos::RCP<panzer::EpetraLinearObjContainer>& ep_container,
                       Teuchos::RCP<panzer::LinearObjContainer>& ghost_container,
                       const ProblemOptions& po)
{
  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  out.setShowProcRank(true);

  ghost_container = linObjFactory->buildGhostedLinearObjContainer();
  const Teuchos::RCP<panzer::LinearObjContainer> container = linObjFactory->buildLinearObjContainer();
  // Convert generic linear object container to epetra container.
  ep_container = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(container);
  linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                            panzer::LinearObjContainer::F |
                                            panzer::LinearObjContainer::Mat, *ghost_container);

  double bnorm(0);
  Teuchos::RCP<Epetra_Vector> x, dx;
  Teuchos::RCP<Epetra_CrsMatrix> A;
  // Newton iteration.
  for (int it = 0; it < 10; ++it) {
    linObjFactory->initializeContainer(panzer::LinearObjContainer::X |
                                       panzer::LinearObjContainer::F |
                                       panzer::LinearObjContainer::Mat, *container);
    container->initialize();
    ghost_container->initialize();

    if (x.is_null()) {
      x = ep_container->get_x();
      dx = Teuchos::rcp(new Epetra_Vector(ep_container->get_x()->Map()));
    } else
      ep_container->set_x(x);

    panzer::AssemblyEngineInArgs input(ghost_container, ep_container);
    input.alpha = 0;
    input.beta = 1;

    // Evaluate physics.
    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);
    if (it == 0 && Teuchos::nonnull(po.J))
      A = po.J;
    else
      A = ep_container->get_A();

    double rnorm, dxnorm;
    ep_container->get_f()->Norm2(&rnorm);
    if (it == 0) bnorm = rnorm;
    dx->Norm2(&dxnorm);
    const bool done = rnorm <= po.rtol*bnorm;
    out << "it " << it << " norm(r) " << rnorm << " norm(dx) " << dxnorm << "\n";
    if (done) break;

    if ( ! solve_Ax_eq_b(*A, *ep_container->get_f(), *dx, 1e-2*po.rtol))
      out << "  Linear solver did not converge; continuing outer iteration.\n";

    x->Update(-1, *dx, 1);
  }
}

struct SparseTriple {
  int m, n;
  std::vector<int> i, j;
  std::vector<double> v;
};

// Serial only. Very straightforward: one eval per column. Could speed up a lot
// using graph coloring.
double
testJacobian (panzer::AssemblyEngine_TemplateManager<panzer::Traits>& ae_tm,
              const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& linObjFactory,
              const Teuchos::RCP<const Epetra_Vector> x0 = Teuchos::null)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  const RCP<panzer::LinearObjContainer> ghost_con = linObjFactory->buildGhostedLinearObjContainer();
  const RCP<panzer::LinearObjContainer> con = linObjFactory->buildLinearObjContainer();
  linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                            panzer::LinearObjContainer::F |
                                            panzer::LinearObjContainer::Mat, *ghost_con);
  const RCP<panzer::EpetraLinearObjContainer>
    ep_con = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(con);

  SparseTriple t;
  { // Finite-difference Jacobian. Should use graph coloring to greatly reduce
    // number of evaluations. But brute force for now.
    RCP<Epetra_Vector> x, f0;
    for (int i = -1; ; ++i) {
      linObjFactory->initializeContainer(panzer::LinearObjContainer::X | panzer::LinearObjContainer::F |
                                         panzer::LinearObjContainer::Mat, *con);
      ghost_con->initialize();
      con->initialize();
      const Epetra_Map& row_map = ep_con->get_A()->RowMap(), col_map = ep_con->get_A()->ColMap();
      if (i == -1) {
        if (Teuchos::nonnull(x0))
          x = rcp(new Epetra_Vector(*x0));
        else
          x = rcp(new Epetra_Vector(ep_con->get_x()->Map()));
        t.m = t.n = x->GlobalLength();
      }
      // For a linear problem, could make delta 1 to remove cancellation
      // error. But I want to look at a nonlinear Robin condition, so do a true
      // finite difference.
      const double delta = 1e-6;
      const bool i_mine = row_map.MyGID(i);
      double x_prev = 0;
      int i_lid = 0;
      if (i_mine && i >= 0) {
        i_lid = col_map.LID(i);
        x_prev = (*x)[i_lid];
        (*x)[i_lid] += delta;
      }
      ep_con->set_x(x);
      panzer::AssemblyEngineInArgs input(ghost_con, ep_con);
      input.alpha = 0;
      input.beta = 1;
      ae_tm.getAsObject<panzer::Traits::Residual>()->evaluate(input);
      if (i == -1)
        f0 = ep_con->get_f();
      else {
        const Epetra_Vector& f = *ep_con->get_f();
        if (i_mine) (*x)[i_lid] = x_prev;
        for (int k = 0; k < f.MyLength(); ++k) {
          const double d = f[k] - (*f0)[k];
          if (d == 0) continue;
          t.i.push_back(row_map.GID(k));
          t.j.push_back(i);
          t.v.push_back(d/delta);
        }
      }
      if (i + 1 == t.m) break;
    }
  }

  { // AD Jacobian.
    linObjFactory->initializeContainer(panzer::LinearObjContainer::X | panzer::LinearObjContainer::F |
                                       panzer::LinearObjContainer::Mat, *con);
    ghost_con->initialize();
    con->initialize();
    if (Teuchos::nonnull(x0))
      ep_con->set_x(rcp(new Epetra_Vector(*x0)));
    else
      ep_con->get_x()->PutScalar(0);
    panzer::AssemblyEngineInArgs input(ghost_con, ep_con);
    input.alpha = 0;
    input.beta = 1;
    ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);
    EpetraExt::RowMatrixToMatrixMarketFile("A_ad.mm", *ep_con->get_A());
  }

  RCP<Epetra_CrsMatrix> A_fd;
  {
    const Epetra_Map& row_map = ep_con->get_A()->RowMap(), col_map = ep_con->get_A()->ColMap();
    A_fd = rcp(new Epetra_CrsMatrix(Copy, row_map, 40));
    for (int i = 0; i < static_cast<int>(t.v.size()); ++i)
      A_fd->InsertGlobalValues(t.i[i], 1, &t.v[i], &t.j[i]);
    A_fd->FillComplete();
  }
  EpetraExt::RowMatrixToMatrixMarketFile("A_fd.mm", *A_fd);

  { // Check error.
    const double A_fd_inf_norm = A_fd->NormInf();
    RCP<Epetra_CrsMatrix> D = rcp(new Epetra_CrsMatrix(Copy, ep_con->get_A()->Graph()));
    D->FillComplete();
    Epetra_CrsMatrix* D_ptr = D.get();
    EpetraExt::MatrixMatrix::Add(*ep_con->get_A(), false, -1, *A_fd, false, 1, D_ptr);
    const double D_inf_norm = D->NormInf();
    return D_inf_norm / A_fd_inf_norm;
  }
}

bool hasInterfaceCondition (const std::vector<panzer::BC>& bcs)
{
  for (std::vector<panzer::BC>::const_iterator bcit = bcs.begin(); bcit != bcs.end(); ++bcit)
    if (bcit->bcType() == panzer::BCT_Interface)
      return true;
  return false;
}

Teuchos::RCP<panzer_stk::STKConnManager<int> >
getSTKConnManager (const Teuchos::RCP<panzer::ConnManagerBase<int> >& conn_mgr)
{
  const Teuchos::RCP<panzer_stk::STKConnManager<int> > stk_conn_mgr =
    Teuchos::rcp_dynamic_cast<panzer_stk::STKConnManager<int> >(conn_mgr);
  TEUCHOS_TEST_FOR_EXCEPTION(stk_conn_mgr.is_null(), std::logic_error,
                             "There are interface conditions, but the connection manager"
                             " does not support the necessary connections.");
  return stk_conn_mgr;
}

void buildInterfaceConnections (const std::vector<panzer::BC>& bcs,
                                const Teuchos::RCP<panzer::ConnManagerBase<int> >& conn_mgr)
{
  const Teuchos::RCP<panzer_stk::STKConnManager<int> >
    stk_conn_mgr = getSTKConnManager(conn_mgr);
  for (std::vector<panzer::BC>::const_iterator bcit = bcs.begin(); bcit != bcs.end(); ++bcit)
    if (bcit->bcType() == panzer::BCT_Interface)
      stk_conn_mgr->associateElementsInSideset(bcit->sidesetID());
}

void checkInterfaceConnections (const Teuchos::RCP<panzer::ConnManagerBase<int> >& conn_mgr,
                                const Teuchos::RCP<Teuchos::Comm<int> >& comm)
{
  const Teuchos::RCP<panzer_stk::STKConnManager<int> >
    stk_conn_mgr = getSTKConnManager(conn_mgr);
  std::vector<std::string> sidesets = stk_conn_mgr->checkAssociateElementsInSidesets(*comm);
  if ( ! sidesets.empty()) {
    std::stringstream ss;
    ss << "Sideset IDs";
    for (std::size_t i = 0; i < sidesets.size(); ++i)
      ss << " " << sidesets[i];
    ss << " did not yield associations, but these sidesets correspond to BCT_Interface BCs.";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, ss.str());
  }
}

void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb, std::vector<panzer::BC>& bcs,
                        const ProblemOptions& po);

// calls MPI_Init and MPI_Finalize
int main (int argc, char* argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using panzer::StrPureBasisPair;
  using panzer::StrPureBasisComp;

  Kokkos::initialize(argc,argv);

  {

    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
    out.setOutputToRootOnly(0);
    out.setShowProcRank(true);
  
    const std::size_t workset_size = 20;
  
    ProblemOptions po;
    {
      // Set up this problem with two discontinuous (A, C) and one continuous (B)
      // fields.
      //   On fields A are imposed Neumann and weak Dirichlet matching interface conditions.
      //   On fields C are imposed Robin interface conditions, with one-way
      // coupling to field B.
      //   If the Robin condition is linear, then the default setup is such that
      // A, B, C all converge to the same solution for which the Solution
      // evaluator provides the exact expression. A response function reports the
      // error so a convergence test can be wrapped around multiple runs of this
      // program.
      //   If the Robin condition is nonlinear, then the source is 0 and the
      // solution is two planes with a jump of 0.4 at the interface.
  
      Teuchos::CommandLineProcessor clp;
      po.nxelem = 10;
      clp.setOption("nx", &po.nxelem, "Number of elements in x direction");
      po.nonlinear_Robin = false;
      clp.setOption("nonlinear", "linear", &po.nonlinear_Robin,
                    "Use a nonlinear Robin interface condition");
      po.rtol = 1e-10;
      clp.setOption("rtol", &po.rtol, "Tolerance on residual norm");
      po.is3d = false;
      clp.setOption("3d", "2d", &po.is3d, "3D test instead of 2D");
      po.mesh_filename = "";
      clp.setOption("mesh-filename", &po.mesh_filename, "Optionally read from an Exodus mesh");
      po.test_Jacobian = false;
      clp.setOption("test-jacobian", "dont-test-jacobian", &po.test_Jacobian,
                    "Test Jacobian using finite differences.");
      po.generate_mesh_only = false;
      clp.setOption("generate-mesh-only", "dont-generate-mesh-only", &po.generate_mesh_only,
                    "Generate mesh, save, and quit.");
      try {
        clp.parse(argc, argv);
      } catch (...) {
        Kokkos::finalize_all();
        return -1;
      }
  
      po.nyelem = po.nxelem;
      po.dof_names.push_back("A");
      po.dof_names.push_back("B");
      po.dof_names.push_back("C");
      po.ss_names.push_back("left");
      po.ss_names.push_back("vertical_0");
      po.ss_names.push_back("right");
      po.outer_iteration = true;
      po.check_error = true;
  
      out << po << "\n";
    }
    bool pass = true;
  
    // Can be overridden by the equation set.
    po.integration_order = 2;
  
    // Construct mesh.
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
    if ( ! po.mesh_filename.empty()) {
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory(po.mesh_filename));
    } else {
      if (po.is3d)
        mesh_factory = Teuchos::rcp(new panzer_stk::CubeHexMeshFactory);
      else
        mesh_factory = Teuchos::rcp(new panzer_stk::SquareQuadMeshFactory);
    }
  
    if (po.mesh_filename.empty()) {
      // set mesh factory parameters
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Blocks",2);
      pl->set("Y Blocks",1);
      if (po.is3d) pl->set("Z Blocks",1);
      pl->set("X Elements", po.nxelem); // per block
      pl->set("Y Elements", po.nyelem);
      if (po.is3d) {
        pl->set("Z Elements", po.nyelem);
        pl->set("Build Interface Sidesets", true);
      }
      { // If np is even, put ranks in both x and y directions; if not, go with
        // default, which is x direction only. The x direction is the harder case.
        const int np = mpiSession.getNProc();
        if (np % 2 == 0 && np >= 4) {
          const int nxp = np/2, nyp = 2;
          pl->set("X Procs", nxp);
          pl->set("Y Procs", nyp);
        }
      }
      mesh_factory->setParameterList(pl);
    }
  
    RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
    if (po.generate_mesh_only) {
      mesh_factory->completeMeshConstruction(*mesh, MPI_COMM_WORLD);
      mesh->writeToExodus("output.exo");
      out << "Stopping after writing mesh because --generate-mesh-only was requested.\n";
      Kokkos::finalize_all();
      return 0;
    }
  
    //todo mesh->getDimension() may not be right if mesh_factory is the Exodus
    // reader.
    po.is3d = mesh->getMetaData()->spatial_dimension() == 3;
  
    if (po.is3d) {
      po.eb_names.push_back("eblock-0_0_0");
      po.eb_names.push_back("eblock-1_0_0");
    } else {
      po.eb_names.push_back("eblock-0_0");
      po.eb_names.push_back("eblock-1_0");
    }
  
    // construct input physics and physics block
    ////////////////////////////////////////////////////////
  
    // factory definitions
    Teuchos::RCP<Example::EquationSetFactory> eqset_factory = 
      Teuchos::rcp(new Example::EquationSetFactory); // where poisson equation is defined
    Example::BCStrategyFactory bc_factory;    // where boundary conditions are defined 
  
    const Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    std::vector<RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      testInitialization(ipb, bcs, po);
  
      std::map<std::string,std::string> block_ids_to_physics_ids;
      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
  
      block_ids_to_physics_ids[po.eb_names[0]] = "Poisson Physics Left";
      block_ids_to_physics_ids[po.eb_names[1]] = "Poisson Physics Right";
  
      block_ids_to_cell_topo[po.eb_names[0]] = mesh->getCellTopology(po.eb_names[0]);
      block_ids_to_cell_topo[po.eb_names[1]] = mesh->getCellTopology(po.eb_names[1]);
        
      // GobalData sets ostream and parameter interface to physics
      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
        
      // the physics block knows how to build and register evaluator with the field manager
      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                 block_ids_to_cell_topo,
                                 ipb,
                                 po.integration_order,
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
        for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr)
          mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
      }
  
      mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);
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
  
    std::vector<std::string> elementBlockNames;
    mesh->getElementBlockNames(elementBlockNames);
    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > volume_worksets;
    panzer::getVolumeWorksetsFromContainer(*wkstContainer,elementBlockNames,volume_worksets);
  
    // build DOF Manager and linear object factory
    /////////////////////////////////////////////////////////////
  
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager;
    {
      const Teuchos::RCP<panzer::ConnManager<int,int> >
        conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager<int>(mesh)); 
      const bool has_interface_condition = hasInterfaceCondition(bcs);
      if (has_interface_condition)
        buildInterfaceConnections(bcs, conn_manager);
      panzer::DOFManagerFactory<int,int> globalIndexerFactory;
      globalIndexerFactory.setUseNeighbors(has_interface_condition);
      dofManager = globalIndexerFactory.buildUniqueGlobalIndexer(
        Teuchos::opaqueWrapper(MPI_COMM_WORLD), physicsBlocks, conn_manager, "");
      if (has_interface_condition)
        checkInterfaceConnections(conn_manager, dofManager->getComm());
    }
  
    // construct some linear algebra object, build object to pass to evaluators
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
      = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));
  
    std::vector<std::string> names;
    std::vector<std::vector<std::string> > eblocks;
    const int c_name_start = 3;
    {
      for (int i = 1; i <= 2; ++i) {
        names.push_back(strint(po.dof_names[0], i));
        eblocks.push_back(std::vector<std::string>());
        eblocks.back().push_back(po.eb_names[i-1]);
      }
      names.push_back(po.dof_names[1]);
      eblocks.push_back(std::vector<std::string>());
      eblocks.back().push_back(po.eb_names[0]);
      eblocks.back().push_back(po.eb_names[1]);
      if (po.dof_names.size() >= 3)
        for (int i = 1; i <= 2; ++i) {
          names.push_back(strint(po.dof_names[2], i));
          eblocks.push_back(std::vector<std::string>());
          eblocks.back().push_back(po.eb_names[i-1]);
        }
    }
    
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > errorResponseLibrary
      = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer, dofManager, linObjFactory));
    {
      for (std::size_t i = po.nonlinear_Robin ? c_name_start : 0; i < names.size(); ++i) {
        panzer::FunctionalResponse_Builder<int,int> builder;
        builder.comm = MPI_COMM_WORLD;
        builder.cubatureDegree = po.integration_order;
        builder.requiresCellIntegral = true;
        builder.quadPointField = names[i] + "_ERROR";
        errorResponseLibrary->addResponse(names[i] + " L2 Error", eblocks[i], builder);
      }
    }
  
    // setup closure model
    /////////////////////////////////////////////////////////////
   
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory; 
    Example::ClosureModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);
  
    Teuchos::ParameterList closure_models("Closure Models");
    {    
      Teuchos::ParameterList& s = closure_models.sublist("solid");
      for (std::vector<std::string>::const_iterator it = names.begin(); it != names.end(); ++it) {
        if (po.nonlinear_Robin)
          s.sublist(std::string("SOURCE_") + *it).set<double>("Value", 0.0);
        else
          s.sublist(std::string("SOURCE_") + *it).set<std::string>("Type", "SIMPLE SOURCE");
      }
      if (po.check_error)
        for (std::size_t i = po.nonlinear_Robin ? c_name_start : 0; i < names.size(); ++i) {
          const std::string err = names[i] + "_ERROR";
          s.sublist(err).set<std::string>("Type", "ERROR_CALC");
          s.sublist(err).set<std::string>("Field A", names[i]);
          s.sublist(err).set<std::string>("Field B", "EXACT");
        }
      if (po.check_error)
        s.sublist("EXACT").set<std::string>("Type", po.nonlinear_Robin ? "EXACT nonlinear Robin" : "EXACT");
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
  
    user_data.set<int>("Workset Size", workset_size);
    if (po.check_error)
      errorResponseLibrary->buildResponseEvaluators(physicsBlocks, cm_factory, closure_models, user_data);
  
    // assemble and solve
    /////////////////////////////////////////////////////////////
    Teuchos::RCP<panzer::EpetraLinearObjContainer> ep_container;
    Teuchos::RCP<panzer::LinearObjContainer> ghost_container;
    if ( ! po.outer_iteration) {
      // Straightfoward solve
      // build linear algebra objects: Ghost is for parallel assembly, it contains
      //                               local element contributions summed, the global IDs
      //                               are not unique. The non-ghosted or "global"
      //                               container will contain the sum over all processors
      //                               of the ghosted objects. The global indices are unique.
      ghost_container = linObjFactory->buildGhostedLinearObjContainer();
      RCP<panzer::LinearObjContainer> container = linObjFactory->buildLinearObjContainer();
      linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                                panzer::LinearObjContainer::F |
                                                panzer::LinearObjContainer::Mat,*ghost_container);
      linObjFactory->initializeContainer(panzer::LinearObjContainer::X |
                                         panzer::LinearObjContainer::F |
                                         panzer::LinearObjContainer::Mat,*container);
      ghost_container->initialize();
      container->initialize();
  
      panzer::AssemblyEngineInArgs input(ghost_container,container);
      input.alpha = 0;
      input.beta = 1;
  
      // evaluate physics: This does both the Jacobian and residual at once
      ae_tm.getAsObject<panzer::Traits::Jacobian>()->evaluate(input);
  
      // solve linear system
      /////////////////////////////////////////////////////////////
      // convert generic linear object container to epetra container
      ep_container = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(container);
  
      // Setup the linear solve: notice A is used directly 
      Epetra_LinearProblem problem(&*ep_container->get_A(),&*ep_container->get_x(),&*ep_container->get_f()); 
  
      // build the solver
      AztecOO solver(problem);
      solver.SetAztecOption(AZ_solver,AZ_gmres); // we don't push out dirichlet conditions
      solver.SetAztecOption(AZ_precond,AZ_none);
      solver.SetAztecOption(AZ_kspace,300);
      solver.SetAztecOption(AZ_output,10);
      solver.SetAztecOption(AZ_precond,AZ_Jacobi);
  
      // solve the linear system
      solver.Iterate(1000,1e-5);
  
      // we have now solved for the residual correction from
      // zero in the context of a Newton solve.
      //     J*e = -r = -(f - J*0) where f = J*u
      // Therefore we have  J*e=-J*u which implies e = -u
      // thus we will scale the solution vector 
      ep_container->get_x()->Scale(-1.0);
    } else { // Some analysis and an outer iteration if necessary.
      Teuchos::RCP<Epetra_CrsMatrix> J_fd;
      assembleAndSolve(ae_tm, linObjFactory, ep_container, ghost_container, po);
      if (po.test_Jacobian) {
        const double nwre = testJacobian(ae_tm, linObjFactory, ep_container->get_x());
        out << "TEST JACOBIAN " << nwre << "\n";
        if (nwre < 0 || nwre > 1e-5) pass = false;
      }
    }
    
    // output data (optional)
    /////////////////////////////////////////////////////////////
  
    // write linear system
    if (false) {
      EpetraExt::RowMatrixToMatrixMarketFile("a_op.mm",*ep_container->get_A());
      EpetraExt::VectorToMatrixMarketFile("x_vec.mm",*ep_container->get_x());
      EpetraExt::VectorToMatrixMarketFile("b_vec.mm",*ep_container->get_f());
    }
  
    if (po.check_error) {
      std::vector<Teuchos::RCP<panzer::Response_Functional<panzer::Traits::Residual> > > rfs(names.size());
      for (std::size_t i = po.nonlinear_Robin ? c_name_start : 0; i < names.size(); ++i) {
        Teuchos::RCP<panzer::ResponseBase>
          resp = errorResponseLibrary->getResponse<panzer::Traits::Residual>(names[i] + " L2 Error");
        rfs[i] = Teuchos::rcp_dynamic_cast<panzer::Response_Functional<panzer::Traits::Residual> >(resp);
        Teuchos::RCP<Thyra::VectorBase<double> > respVec = Thyra::createMember(rfs[i]->getVectorSpace());
        rfs[i]->setVector(respVec);
      }
  
      panzer::AssemblyEngineInArgs respInput(ghost_container, ep_container);
      respInput.alpha = 0;
      respInput.beta = 1;
      errorResponseLibrary->addResponsesToInArgs<panzer::Traits::Residual>(respInput);
      errorResponseLibrary->evaluate<panzer::Traits::Residual>(respInput);
  
      // Record a max error so we can use convergence_rate.py.
      double max_err = -1;
      for (std::size_t i = po.nonlinear_Robin ? c_name_start : 0; i < names.size(); ++i) {
        const double err = sqrt(rfs[i]->value);
        max_err = std::max(max_err, err);
        out << names[i] << " ERROR = " << err << "\n";
        if (err < 0 || err > (po.nonlinear_Robin ? 1e-10 : 0.03/(po.nxelem*po.nxelem/25.0)))
          pass = false;
      }
      out << "Error = " << max_err << "\n";
    }
  
    // Write solution except in the special case of a generated 3D mesh and #rank
    // > 1. In that case, something in the mesh-gen and rebalance code is causing
    // a failure in IossBridge::write_side_data_to_ioss.
    if ( ! (po.is3d && mpiSession.getNProc() > 1 && po.mesh_filename.empty())) {
      // redistribute solution vector to ghosted vector
      linObjFactory->globalToGhostContainer(*ep_container,*ghost_container,
                                            panzer::EpetraLinearObjContainer::X 
                                            | panzer::EpetraLinearObjContainer::DxDt); 
  
      // get X Epetra_Vector from ghosted container
      RCP<panzer::EpetraLinearObjContainer> ep_ghost_container =
        rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ghost_container);
      panzer_stk::write_solution_data(*dofManager,*mesh,*ep_ghost_container->get_x());
      mesh->writeToExodus("output.exo");
    }
  
    // all done!
    /////////////////////////////////////////////////////////////
    out << (pass ? "PASS" : "FAIL") << " BASICS\n";
  }

  return 0;
}

void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb, std::vector<panzer::BC>& bcs,
                        const ProblemOptions& po)
{
  const char* const pb_strings[2] = {"Poisson Physics Left", "Poisson Physics Right"};
  const double dleft = 0.5, dright = -0.3;
  size_t bc_id = 0;

  { // A discontinuous field with matching value (in a weak sense) and
    // derivative at the interface.
    for (int i = 0; i < 2; ++i) {
      Teuchos::ParameterList& p = ipb->sublist(pb_strings[i]).sublist("Equation Set 1");
      p.set("Type", "Poisson");
      p.set("Model ID", "solid");
      p.set("Basis Type", "HGrad");
      p.set("Basis Order", 1);
      p.set("Integration Order", po.integration_order);
      p.set("DOF Name", strint(po.dof_names[0], i+1));
    }

    const std::string t1_name = strint(po.dof_names[0], 1), t2_name = strint(po.dof_names[0], 2);
    {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = po.ss_names[0];
      std::string element_block_id = po.eb_names[0];
      std::string dof_name = t1_name;
      std::string strategy = "Constant";
      double value = dleft;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name, strategy, p);
      bcs.push_back(bc);
    }
    {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = po.ss_names[2];
      std::string element_block_id = po.eb_names[1];
      std::string dof_name = t2_name;
      std::string strategy = "Constant";
      double value = dright;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name, strategy, p);
      bcs.push_back(bc);
    }

    // There is no reason to turn off either of these, but in early development
    // it was useful.
    const bool
      neumann_match = 1,
      weak_dirichlet_match = 1;

    if (neumann_match) {
      Teuchos::ParameterList p;
      p.set("Type", "Interface");
      p.set("Sideset ID", po.ss_names[1]);
      p.set("Element Block ID" , po.eb_names[0]);
      p.set("Element Block ID2", po.eb_names[1]);
      p.set("Equation Set Name" , t1_name);
      p.set("Equation Set Name2", t2_name);
      p.set("Strategy", "Neumann Match Interface");
      bcs.push_back(panzer::BC(bc_id++, p));
    } else {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = po.ss_names[1];
      std::string element_block_id = po.eb_names[0];
      std::string dof_name = t1_name;
      std::string strategy = "Constant";
      double value = -0.4;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name, strategy, p);
      bcs.push_back(bc);
    }

    if (weak_dirichlet_match) {
      Teuchos::ParameterList p;
      p.set("Type", "Interface");
      p.set("Sideset ID", po.ss_names[1]);
      p.set("Element Block ID" , po.eb_names[1]);
      p.set("Element Block ID2", po.eb_names[0]);
      p.set("Equation Set Name" , t2_name);
      p.set("Equation Set Name2", t1_name);
      p.set("Strategy", "Weak Dirichlet Match Interface");
      bcs.push_back(panzer::BC(bc_id++, p));
    } else {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = po.ss_names[1];
      std::string element_block_id = po.eb_names[1];
      std::string dof_name = t2_name;
      std::string strategy = "Constant";
      double value = 0.4;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name, strategy, p);
      bcs.push_back(bc);
    }
  }

  { // A field that is continuous across the interface. It has 1 DOF per
    // interface node.
    for (int ipbs = 0; ipbs < 2; ++ipbs) {
      Teuchos::ParameterList& p = ipb->sublist(pb_strings[ipbs]).sublist("Equation Set 2");
      p.set("Type","Poisson");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",po.integration_order);
      p.set("DOF Name", po.dof_names[1]);
    }
    {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = po.ss_names[0];
      std::string element_block_id = po.eb_names[0];
      std::string dof_name = po.dof_names[1];
      std::string strategy = "Constant";
      double value = dleft;
      Teuchos::ParameterList p;
      p.set("Value", value);
      panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name, strategy, p);
      bcs.push_back(bc);
    }
    {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = po.ss_names[2];
      std::string element_block_id = po.eb_names[1];
      std::string dof_name = po.dof_names[1];
      std::string strategy = "Constant";
      double value = dright;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name, strategy, p);
      bcs.push_back(bc);
    }
  }

  if (po.dof_names.size() >= 3) {
    // Another discontinuous field that will couple to the continuous field at
    // the interface.
    for (int i = 0; i < 2; ++i) {
      Teuchos::ParameterList& p = ipb->sublist(pb_strings[i]).sublist("Equation Set 3");
      p.set("Type", "Poisson");
      p.set("Model ID", "solid");
      p.set("Basis Type", "HGrad");
      p.set("Basis Order", 1);
      p.set("Integration Order", po.integration_order);
      p.set("DOF Name", strint(po.dof_names[2], i+1));
    }

    const std::string t1_name = strint(po.dof_names[2], 1), t2_name = strint(po.dof_names[2], 2);
    {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = po.ss_names[0];
      std::string element_block_id = po.eb_names[0];
      std::string dof_name = t1_name;
      std::string strategy = "Constant";
      double value = dleft;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name, strategy, p);
      bcs.push_back(bc);
    }
    {
      panzer::BCType bctype = panzer::BCT_Dirichlet;
      std::string sideset_id = po.ss_names[2];
      std::string element_block_id = po.eb_names[1];
      std::string dof_name = t2_name;
      std::string strategy = "Constant";
      double value = dright;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id++, bctype, sideset_id, element_block_id, dof_name, strategy, p);
      bcs.push_back(bc);
    }

    for (int ibc = 0; ibc < 2; ++ibc) {
      Teuchos::ParameterList p;
      p.set("Type", "Interface");
      p.set("Sideset ID", po.ss_names[1]);
      p.set("Element Block ID" , po.eb_names[ibc]);
      p.set("Element Block ID2", po.eb_names[(ibc + 1) % 2]);
      p.set("Equation Set Name" , ibc == 0 ? t1_name : t2_name);
      p.set("Equation Set Name2", ibc == 0 ? t2_name : t1_name);
      p.set("Strategy", "Robin Interface");
      Teuchos::ParameterList d;
      d.set("Coupling DOF Name", po.dof_names[1]);
      if (po.nonlinear_Robin) {
        // There's no real meaning to these cofficients; I just solved for them
        // to get the solution I want.
        const double c[] = {1.1, 0.1, 1.66};
        d.set("a", c[0]);
        d.set("b", c[1]);
        d.set("c", c[2]);
      } else {
        const double c[] = {1, -1, 1};
        d.set("a", c[0]);
        // Switch the two so that the same coefficients are applied to C_left and
        // C_right. Switch sign since the normal has opposite direction.
        d.set("b", (ibc == 0 ? c[1] : -c[2]));
        d.set("c", (ibc == 0 ? c[2] : -c[1]));
      };
      d.set("Nonlinear", po.nonlinear_Robin);
      p.set("Data", d);
      bcs.push_back(panzer::BC(bc_id++, p));
    }
  }
}
