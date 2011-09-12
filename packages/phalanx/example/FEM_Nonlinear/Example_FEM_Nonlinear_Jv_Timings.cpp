// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_FancyOStream.hpp"

// User defined objects
#include "Element_Linear2D.hpp"
#include "Workset.hpp"
#include "Traits.hpp"
#include "FactoryTraits.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "MeshBuilder.hpp"
#include "LinearObjectFactory.hpp"

// Linear solver
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// Preconditioner
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "Ifpack.h"
#include "ml_epetra_preconditioner.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void printVector(std::string filename_prefix, const Epetra_Vector& vector, 
		 int newton_step)
{
  std::stringstream ss;
  ss << filename_prefix << "_" << newton_step << ".dat";
  ofstream file( ss.str().c_str(), ios::out | ios::app );
  vector.Print(file);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void printMatrix(std::string filename_prefix, const Epetra_CrsMatrix& matrix,
		 int newton_step)
{
  std::stringstream ss;
  ss << filename_prefix << "_" << newton_step << ".dat";
  ofstream file( ss.str().c_str(), ios::out | ios::app );
  matrix.Print(file);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void checkConvergence(int my_pid,
		      const int num_newton_steps, const Epetra_Vector& f, 
		      const Epetra_Vector& delta_x, bool& converged)
{
  double norm_f = 0.0;
  f.Norm2(&norm_f);
  double norm_delta_x = 0.0;
  delta_x.Norm2(&norm_delta_x);

  if ( (norm_f < 1.0e-8) && (norm_delta_x < 1.0e-5) )
    converged = true;
  else
    converged = false;

  if (my_pid == 0)
    cout << "Step: " << num_newton_steps << ", ||f|| = "  << norm_f
	 << ", ||delta x|| = " << norm_delta_x << endl;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void applyBoundaryConditions(const double& value,
			     const Epetra_Vector& x, Epetra_CrsMatrix& Jac, 
			     Epetra_Vector& f, const MeshBuilder mb)
{
  // This function assumes a serial run.  It will not check process rank.
  // Assumes two equations.

  // Left edge includes global nodes 0 and 1, dofs 0,1,2,3
  int num_nonzero_entries_in_row = 0;
  int* row_indices;
  double* row_values;

  int num_eq = 2;

  const std::vector<int>& left_nodes = mb.leftNodeSetGlobalIds();
  for (std::size_t node = 0; node < left_nodes.size();  ++node) {
    
    for (int eq = 0; eq < num_eq; ++eq) {

      int lid = Jac.RowMap().LID(left_nodes[node] * num_eq + eq);

      Jac.ExtractMyRowView(lid, num_nonzero_entries_in_row, row_values, 
			   row_indices);
      
      for (int col=0; col < num_nonzero_entries_in_row; ++col) {
	if (row_indices[col] == lid)
	  row_values[col] = 1.0;
	else
	  row_values[col] = 0.0;
      }
      f[lid] = value - x[lid];
    
    }

  }
  
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI  
  MPI_Init(&argc, &argv);
#endif

  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    RCP<Time> residual_eval_time = 
      TimeMonitor::getNewTimer("Residual Evaluation Time");
    RCP<Time> jacobian_eval_time = 
      TimeMonitor::getNewTimer("Jacobian Evaluation Time");
    RCP<Time> linear_solve_time = 
      TimeMonitor::getNewTimer("Linear Solve Time");
    RCP<Time> nonlinear_solve_time = 
      TimeMonitor::getNewTimer("Nonlinear Solve Time");
    RCP<Time> preconditioner_solve_time = 
      TimeMonitor::getNewTimer("Preconditioner Time");
    RCP<Time> setup_time = 
      TimeMonitor::getNewTimer("Setup Time (not scalable)");
    RCP<Time> jv_eval_time = 
      TimeMonitor::getNewTimer("Jv (AD)");
    RCP<Time> matvec = 
      TimeMonitor::getNewTimer("Matvec");

    setup_time->start();

    bool print_debug_info = false;

#ifdef HAVE_MPI
    RCP<Epetra_Comm> comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    RCP<Epetra_Comm> comm = rcp(new Epetra_SerialComm);
#endif
    
    Teuchos::basic_FancyOStream<char> os(rcp(&std::cout,false));
    os.setShowProcRank(true);
    os.setProcRankAndSize(comm->MyPID(), comm->NumProc());

    if (comm->MyPID() == 0)
      cout << "\nStarting FEM_Nonlinear Example!\n" << endl;

    // *********************************************************
    // * Build the Finite Element data structures
    // *********************************************************

    // Problem dimension - a 2D problem
    const static int dim = 2;

    // Create the mesh
    MeshBuilder mb(comm, 10, 3, 1.0, 1.0, 8);

    if (print_debug_info) 
      os << mb;

    std::vector<Element_Linear2D>& cells = *(mb.myElements());

    // Divide mesh into workset blocks
    const std::size_t workset_size = 15;
    std::vector<MyWorkset> worksets;
    {
      std::vector<Element_Linear2D>::iterator cell_it = cells.begin();
      std::size_t count = 0;
      MyWorkset w;
      w.local_offset = cell_it->localElementIndex();
      w.begin = cell_it;
      for (; cell_it != cells.end(); ++cell_it) {
	++count;
	std::vector<Element_Linear2D>::iterator next = cell_it;
	++next;
	
	if ( count == workset_size || next == cells.end()) {
	  w.end = next;
	  w.num_cells = count;
	  worksets.push_back(w);
	  count = 0;
	  
	  if (next != cells.end()) {
	    w.local_offset = next->localElementIndex();
	    w.begin = next;
	  }
	}
      }
    }
    
    if (print_debug_info) {
      cout << "Printing Element Information" << endl;
      for (std::size_t i = 0; i < worksets.size(); ++i) {
	std::vector<Element_Linear2D>::iterator it = worksets[i].begin;
	for (; it != worksets[i].end; ++it)
	  cout << *it << endl;
      }
    }
    
    if (print_debug_info) {
      for (std::size_t i = 0; i < worksets.size(); ++i) {
	cout << "Printing Workset Information" << endl;
	cout << "worksets[" << i << "]" << endl;
	cout << "  num_cells =" << worksets[i].num_cells << endl;
	cout << "  local_offset =" << worksets[i].local_offset << endl;
	std::vector<Element_Linear2D>::iterator it = worksets[i].begin;
	for (; it != worksets[i].end; ++it)
	  cout << "  cell_local_index =" << it->localElementIndex() << endl;
      }
      cout << endl;
    }

    // *********************************************************
    // * Build the Newton solver data structures
    // *********************************************************

    // Setup Nonlinear Problem (build Epetra_Vector and Epetra_CrsMatrix)
    // Newton's method: J delta_x = -f
    const std::size_t num_eq = 2;
    
    LinearObjectFactory lof(mb, comm, num_eq);

    if (print_debug_info) {
      ofstream file("OwnedGraph.dat", ios::out | ios::app);
      Teuchos::basic_FancyOStream<char> p(rcp(&file,false)); 
      p.setShowProcRank(true); 
      p.setProcRankAndSize(comm->MyPID(), comm->NumProc()); 	
      lof.ownedGraph()->Print(p);
    }

    Epetra_Map owned_map = *(lof.ownedMap());
    Epetra_Map overlapped_map = *(lof.overlappedMap());
    Epetra_CrsGraph owned_graph = *(lof.ownedGraph());
    Epetra_CrsGraph overlapped_graph = *(lof.overlappedGraph());

    // Solution vector x
    RCP<Epetra_Vector> owned_x = rcp(new Epetra_Vector(owned_map,true));
    RCP<Epetra_Vector> overlapped_x =  
      rcp(new Epetra_Vector(overlapped_map,true));

    // Update vector x
    RCP<Epetra_Vector> owned_delta_x = rcp(new Epetra_Vector(owned_map,true));

    // Residual vector f
    RCP<Epetra_Vector> owned_f = rcp(new Epetra_Vector(owned_map,true));
    RCP<Epetra_Vector> overlapped_f =  
      rcp(new Epetra_Vector(overlapped_map,true));

    // Jacobian Matrix
    Epetra_DataAccess copy = ::Copy;
    RCP<Epetra_CrsMatrix> owned_jac = 
      rcp(new Epetra_CrsMatrix(copy, owned_graph));
    RCP<Epetra_CrsMatrix> overlapped_jac = 
      rcp(new Epetra_CrsMatrix(copy, overlapped_graph));

    // Import/export
    RCP<Epetra_Import> importer = 
      rcp(new Epetra_Import(overlapped_map, owned_map));
    RCP<Epetra_Export> exporter = 
      rcp(new Epetra_Export(overlapped_map, owned_map));

    // Sets bc for initial guess
    applyBoundaryConditions(1.0, *owned_x, *owned_jac, *owned_f, mb);
    
    // *********************************************************
    // * Build the FieldManager
    // *********************************************************

    RCP< vector<string> > dof_names = rcp(new vector<string>(num_eq));
    (*dof_names)[0] = "Temperature";
    (*dof_names)[1] = "Velocity X";

    RCP<DataLayout> qp_scalar = 
      rcp(new MDALayout<Cell,QuadPoint>(workset_size,4));
    RCP<DataLayout> node_scalar = 
      rcp(new MDALayout<Cell,Node>(workset_size,4));
    
    RCP<DataLayout> qp_vec = 
      rcp(new MDALayout<Cell,QuadPoint,Dim>(workset_size,4,dim));
    RCP<DataLayout> node_vec = 
      rcp(new MDALayout<Cell,Node,Dim>(workset_size,4,dim));

    RCP<DataLayout> dummy = rcp(new MDALayout<Cell>(0));
    
    map<string, RCP<ParameterList> > evaluators_to_build;
    
    { // Gather Solution
      RCP<ParameterList> p = rcp(new ParameterList);
      int type = MyFactoryTraits<MyTraits>::id_gather_solution;
      p->set<int>("Type", type);
      p->set< RCP< vector<string> > >("Solution Names", dof_names);
      p->set< RCP<Epetra_Vector> >("Solution Vector", overlapped_x);
      p->set< RCP<DataLayout> >("Data Layout", node_scalar);
      evaluators_to_build["Gather Solution"] = p;
    }
    
    { // FE Interpolation - Temperature
      RCP<ParameterList> p = rcp(new ParameterList);
      
      int type = MyFactoryTraits<MyTraits>::id_feinterpolation;
      p->set<int>("Type", type);
      
      p->set<string>("Node Variable Name", "Temperature");
      p->set<string>("QP Variable Name", "Temperature");
      p->set<string>("Gradient QP Variable Name", "Temperature Gradient");
      
      p->set< RCP<DataLayout> >("Node Data Layout", node_scalar);
      p->set< RCP<DataLayout> >("QP Scalar Data Layout", qp_scalar);
      p->set< RCP<DataLayout> >("QP Vector Data Layout", qp_vec);
      
      evaluators_to_build["FE Interpolation Temperature"] = p;
    }
     
    { // FE Interpolation - Velocity X
      RCP<ParameterList> p = rcp(new ParameterList);
      
      int type = MyFactoryTraits<MyTraits>::id_feinterpolation;
      p->set<int>("Type", type);
      
      p->set<string>("Node Variable Name", "Velocity X");
      p->set<string>("QP Variable Name", "Velocity X");
      p->set<string>("Gradient QP Variable Name", "Velocity X Gradient");
      
      p->set< RCP<DataLayout> >("Node Data Layout", node_scalar);
      p->set< RCP<DataLayout> >("QP Scalar Data Layout", qp_scalar);
      p->set< RCP<DataLayout> >("QP Vector Data Layout", qp_vec);
      
      evaluators_to_build["FE Interpolation Velocity X"] = p;
    }

    { // Evaluate residual
      RCP<ParameterList> p = rcp(new ParameterList);
      int type = MyFactoryTraits<MyTraits>::id_equations;
      p->set<int>("Type", type);
      p->set< RCP< vector<string> > >("Solution Names", dof_names);
      p->set< RCP<DataLayout> >("Node Data Layout", node_scalar);
      p->set< RCP<DataLayout> >("QP Data Layout", qp_scalar);
      p->set< RCP<DataLayout> >("Gradient QP Data Layout", qp_vec);
      evaluators_to_build["Equations"] = p;
    }
 
    { // Scatter Solution
      RCP<ParameterList> p = rcp(new ParameterList);
      int type = MyFactoryTraits<MyTraits>::id_scatter_residual;
      p->set<int>("Type", type);

      RCP< vector<string> > res_names = rcp(new vector<string>(num_eq));
      (*res_names)[0] = "Residual Temperature";
      (*res_names)[1] = "Residual Velocity X";

      p->set< RCP< vector<string> > >("Residual Names", res_names);
      p->set< RCP<Epetra_Vector> >("Residual Vector", overlapped_f);
      p->set< RCP<Epetra_CrsMatrix> >("Jacobian Matrix", overlapped_jac);
      p->set< RCP<DataLayout> >("Dummy Data Layout", dummy);
      p->set< RCP<DataLayout> >("Data Layout", node_scalar);
      evaluators_to_build["Scatter Residual"] = p;
    }
    
    // Build Field Evaluators for each evaluation type
    EvaluatorFactory<MyTraits,MyFactoryTraits<MyTraits> > factory;
    RCP< vector< RCP<Evaluator_TemplateManager<MyTraits> > > > 
      evaluators;
    evaluators = factory.buildEvaluators(evaluators_to_build);
    
    // Create a FieldManager
    FieldManager<MyTraits> fm;
      
    // Register all Evaluators 
    registerEvaluators(evaluators, fm);

    // Request quantities to assemble RESIDUAL PDE operators
    {
      typedef MyTraits::Residual::ScalarT ResScalarT;
      Tag<ResScalarT> res_tag("Scatter", dummy);
      fm.requireField<MyTraits::Residual>(res_tag);
      
      // Request quantities to assemble JACOBIAN PDE operators
      typedef MyTraits::Jacobian::ScalarT JacScalarT;
      Tag<JacScalarT> jac_tag("Scatter", dummy);
      fm.requireField<MyTraits::Jacobian>(jac_tag);

      // Request quantities to assemble Jv operators
      typedef MyTraits::Jv::ScalarT JvScalarT;
      Tag<JvScalarT> jv_tag("Scatter", dummy);
      fm.requireField<MyTraits::Jv>(jv_tag);
    }

    {
      RCP<Time> registration_time = 
	TimeMonitor::getNewTimer("Post Registration Setup Time");
      {
	TimeMonitor t(*registration_time);
	fm.postRegistrationSetupForType<MyTraits::Residual>(NULL);
	fm.postRegistrationSetupForType<MyTraits::Jacobian>(NULL);
	fm.postRegistrationSetupForType<MyTraits::Jv>(NULL);
      }
    }

    if (print_debug_info)
      cout << fm << endl;

    // *********************************************************
    // * Evaluate Jacobian and Residual (required for ML to 
    // * to be constructed properly
    // *********************************************************    
    {
      //TimeMonitor t(*jacobian_eval_time);

      overlapped_x->Import(*owned_x, *importer, Insert);
      
      owned_f->PutScalar(0.0);
      overlapped_f->PutScalar(0.0);
      owned_jac->PutScalar(0.0);
      overlapped_jac->PutScalar(0.0);
      
      for (std::size_t i = 0; i < worksets.size(); ++i)
	fm.evaluateFields<MyTraits::Jacobian>(worksets[i]);
      
      owned_f->Export(*overlapped_f, *exporter, Add);
      owned_jac->Export(*overlapped_jac, *exporter, Add);
      
      applyBoundaryConditions(1.0, *owned_x, *owned_jac, *owned_f, mb);
    }

    // *********************************************************
    // * Build Preconditioner (Ifpack or ML)
    // *********************************************************    
    bool use_ml = false;

    RCP<Belos::EpetraPrecOp> belosPrec;
    RCP<Ifpack_Preconditioner> ifpack_prec;
    RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec;

    if (!use_ml) {
      Ifpack Factory;
      std::string PrecType = "ILU";
      int OverlapLevel = 0;
      ifpack_prec = 
	Teuchos::rcp( Factory.Create(PrecType,owned_jac.get(),OverlapLevel) );
      ParameterList ifpackList;
      ifpackList.set("fact: drop tolerance", 1e-9);
      ifpackList.set("fact: level-of-fill", 1);
      ifpackList.set("schwarz: combine mode", "Add");
      IFPACK_CHK_ERR(ifpack_prec->SetParameters(ifpackList));
      IFPACK_CHK_ERR(ifpack_prec->Initialize());
      belosPrec = rcp( new Belos::EpetraPrecOp( ifpack_prec ) );
    }
    else {
      ParameterList ml_params;
      ML_Epetra::SetDefaults("SA",ml_params);
      //ml_params.set("Base Method Defaults", "SA");
      ml_params.set("ML label", "Phalanx_Test");
      ml_params.set("ML output", 10);
      ml_params.set("print unused", 1);
      ml_params.set("max levels", 4);
      ml_params.set("PDE equations",2);
      ml_params.set("prec type","MGV");
      ml_params.set("increasing or decreasing","increasing");
      // ml_params.set("aggregation: nodes per aggregate",50);
      ml_params.set("aggregation: type","Uncoupled");
      ml_params.set("aggregation: damping factor", 0.0);
      ml_params.set("coarse: type","Amesos-KLU");
      //ml_params.set("coarse: type","IFPACK");
      ml_params.set("coarse: max size", 1000);
      //ml_params.set("smoother: type","IFPACK");
      ml_params.set("smoother: type","block Gauss-Seidel");
      ml_params.set("smoother: ifpack type","ILU");
      ml_params.set("smoother: ifpack overlap",1);
      ml_params.sublist("smoother: ifpack list").set("fact: level-of-fill",1);
      ml_params.sublist("smoother: ifpack list").set("schwarz: reordering type","rcm");
      ml_prec = rcp( new ML_Epetra::MultiLevelPreconditioner(*owned_jac, 
							     ml_params, 
							     true) );
      belosPrec = rcp( new Belos::EpetraPrecOp( ml_prec ) );
    }

    // *********************************************************
    // * Build linear solver (Belos)
    // *********************************************************
    
    // Linear solver parameters
    typedef double                            ST;
    typedef Teuchos::ScalarTraits<ST>        SCT;
    typedef SCT::magnitudeType                MT;
    typedef Epetra_MultiVector                MV;
    typedef Epetra_Operator                   OP;
    typedef Belos::MultiVecTraits<ST,MV>     MVT;
    typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
    
    RCP<ParameterList> belosList = rcp(new ParameterList);
    belosList->set<int>("Num Blocks", 400);
    belosList->set<int>("Block Size", 1);
    belosList->set<int>("Maximum Iterations", 400);
    belosList->set<int>("Maximum Restarts", 0);
    belosList->set<MT>( "Convergence Tolerance", 1.0e-4);
    int verbosity = Belos::Errors + Belos::Warnings;
    if (false) {
      verbosity += Belos::TimingDetails + Belos::StatusTestDetails;
      belosList->set<int>( "Output Frequency", -1);
    }
    if (print_debug_info) {
      verbosity += Belos::Debug;
      belosList->set<int>( "Output Frequency", -1);
    }
    belosList->set( "Verbosity", verbosity );
    
    RCP<Epetra_MultiVector> F = 
      Teuchos::rcp_implicit_cast<Epetra_MultiVector>(owned_f);
    
    RCP<Epetra_MultiVector> DX = 
      Teuchos::rcp_implicit_cast<Epetra_MultiVector>(owned_delta_x);
    
    RCP<Belos::LinearProblem<double,MV,OP> > problem =
      rcp(new Belos::LinearProblem<double,MV,OP>(owned_jac, DX, F) );

    problem->setRightPrec( belosPrec );

    RCP< Belos::SolverManager<double,MV,OP> > gmres_solver = 
      rcp( new Belos::BlockGmresSolMgr<double,MV,OP>(problem, belosList) );
    
    setup_time->stop();










    // Timing for Mat-Vec using sacado
    RCP<Epetra_Vector> owned_v = rcp(new Epetra_Vector(owned_map,true));
    RCP<Epetra_Vector> overlapped_v = 
      rcp(new Epetra_Vector(overlapped_map,true));
    RCP<Epetra_Vector> owned_Jv = rcp(new Epetra_Vector(owned_map,true));
    RCP<Epetra_Vector> overlapped_Jv = 
      rcp(new Epetra_Vector(overlapped_map,true));

    owned_x->PutScalar(1.0);
    owned_v->PutScalar(1.0);

    int iter = 0;
    while (iter != 30) {

      overlapped_x->Import(*owned_x, *importer, Insert);

      overlapped_v->Import(*owned_v, *importer, Insert);

      owned_f->PutScalar(0.0);
      overlapped_f->PutScalar(0.0);
      owned_jac->PutScalar(0.0);
      overlapped_jac->PutScalar(0.0);
      owned_Jv->PutScalar(0.0);
      overlapped_Jv->PutScalar(0.0);

      for (std::size_t i = 0; i < worksets.size(); ++i) {
	worksets[i].v = overlapped_v;
	worksets[i].Jv = overlapped_Jv;
      }
      {
	TimeMonitor t(*residual_eval_time);
	for (std::size_t i = 0; i < worksets.size(); ++i)
	  fm.evaluateFields<MyTraits::Residual>(worksets[i]);
      }
      {
	TimeMonitor t(*jacobian_eval_time);
	for (std::size_t i = 0; i < worksets.size(); ++i)
	  fm.evaluateFields<MyTraits::Jacobian>(worksets[i]);
      }
      {
	TimeMonitor t(*jv_eval_time);
	for (std::size_t i = 0; i < worksets.size(); ++i)
	  fm.evaluateFields<MyTraits::Jv>(worksets[i]);
      }

      owned_f->Export(*overlapped_f, *exporter, Add);
      owned_jac->Export(*overlapped_jac, *exporter, Add);
      owned_Jv->Export(*overlapped_Jv, *exporter, Add);

      applyBoundaryConditions(1.0, *owned_x, *owned_jac, *owned_f, mb);

      {
	TimeMonitor t(*matvec);
	owned_jac->Apply(*owned_x, *owned_f);
      }

      ++iter;
    }

    //owned_jac->Print(std::cout);
    //overlapped_Jv->Print(std::cout);
    //owned_Jv->Print(std::cout);
    




    


    // NOTE: in the future we can set up the nonlinear solver below to
    // do Jacobian-Free Newton-Krylov solves to test the matvec


    /*


    // *********************************************************
    // * Solve the system
    // *********************************************************

    // Set initial guess
    owned_x->PutScalar(1.0);

    // Evaluate Residual
    {
      TimeMonitor t(*residual_eval_time);

      overlapped_x->Import(*owned_x, *importer, Insert);

      owned_f->PutScalar(0.0);
      overlapped_f->PutScalar(0.0);

      for (std::size_t i = 0; i < worksets.size(); ++i)
	fm.evaluateFields<MyTraits::Residual>(worksets[i]);
      
      owned_f->Export(*overlapped_f, *exporter, Add);

      applyBoundaryConditions(1.0, *owned_x, *owned_jac, *owned_f, mb);
    }
    
    if (print_debug_info) {
      printVector("x_owned", *owned_x, -1);
      printVector("f_owned", *owned_f, -1);
    }

    // Newton Loop
    bool converged = false;
    std::size_t num_newton_steps = 0;
    std::size_t num_gmres_iterations = 0;
    checkConvergence(comm->MyPID(), num_newton_steps, *owned_f, 
		     *owned_delta_x, converged);

    while (!converged && num_newton_steps < 20) {
      
      TimeMonitor t(*nonlinear_solve_time);

      // Evaluate Residual and Jacobian
      {
	TimeMonitor t(*jacobian_eval_time);

	overlapped_x->Import(*owned_x, *importer, Insert);

	owned_f->PutScalar(0.0);
	overlapped_f->PutScalar(0.0);
	owned_jac->PutScalar(0.0);
	overlapped_jac->PutScalar(0.0);
	
	for (std::size_t i = 0; i < worksets.size(); ++i)
	  fm.evaluateFields<MyTraits::Jacobian>(worksets[i]);

	owned_f->Export(*overlapped_f, *exporter, Add);
	owned_jac->Export(*overlapped_jac, *exporter, Add);

	applyBoundaryConditions(1.0, *owned_x, *owned_jac, *owned_f, mb);
      }

      if (print_debug_info) {
	printVector("x_owned", *owned_x, num_newton_steps);
	printVector("x_overlapped", *overlapped_x, num_newton_steps);
	printVector("f_owned", *owned_f, num_newton_steps);
	printMatrix("jacobian_owned", *owned_jac, num_newton_steps);
      }
      
      owned_f->Scale(-1.0);

      // Solve linear problem
      {
	TimeMonitor t(*linear_solve_time);
	
	owned_delta_x->PutScalar(0.0);

	{
	  TimeMonitor tp(*preconditioner_solve_time);
	  
	  if (use_ml)
	    ml_prec->ReComputePreconditioner();
	  else
	    IFPACK_CHK_ERR(ifpack_prec->Compute());
	}

	problem->setProblem();

	Belos::ReturnType ret = gmres_solver->solve();

	int num_iters = gmres_solver->getNumIters();
	num_gmres_iterations += num_iters; 
	//if (print_debug_info) 
	if (comm->MyPID() == 0)
	  std::cout << "Number of gmres iterations performed for this solve: " 
		    << num_iters << std::endl;
	
	if (ret!=Belos::Converged && comm->MyPID() == 0) {
	  std::cout << std::endl << "WARNING:  Belos did not converge!" 
		    << std::endl;
	}
	
      }
      
      owned_x->Update(1.0, *owned_delta_x, 1.0);

      { // Evaluate Residual Only
	TimeMonitor t(*residual_eval_time);
	
	overlapped_x->Import(*owned_x, *importer, Insert);

	owned_f->PutScalar(0.0);
	overlapped_f->PutScalar(0.0);

	for (std::size_t i = 0; i < worksets.size(); ++i)
	  fm.evaluateFields<MyTraits::Residual>(worksets[i]);
	
	owned_f->Export(*overlapped_f, *exporter, Add);

	applyBoundaryConditions(1.0, *owned_x, *owned_jac, *owned_f, mb);
      }

      num_newton_steps += 1;

      checkConvergence(comm->MyPID(), num_newton_steps, *owned_f, 
		       *owned_delta_x, converged);
    }
     
    if (print_debug_info)
      printVector("f_owned", *owned_f, num_newton_steps);

    if (comm->MyPID() == 0) {
      if (converged)
	cout << "\n\nNewton Iteration Converged!\n" << endl;
      else
	cout << "\n\nNewton Iteration Failed to Converge!\n" << endl;
    }

    RCP<Time> file_io = 
      TimeMonitor::getNewTimer("File IO");
    {
      TimeMonitor t(*file_io);
      
      // Create a list of node coordinates
      std::map<int, std::vector<double> > coordinates;
      Teuchos::RCP< std::vector<Element_Linear2D> > cells = mb.myElements();
      for (std::vector<Element_Linear2D>::iterator cell = cells->begin();
	   cell != cells->end(); ++cell) {
	
	const shards::Array<double,shards::NaturalOrder,Node,Dim>& coords = 
	  cell->nodeCoordinates();

	for (int node=0; node < cell->numNodes(); ++node) {
	  coordinates[cell->globalNodeId(node)].resize(dim);
	  coordinates[cell->globalNodeId(node)][0] = coords(node,0);
	  coordinates[cell->globalNodeId(node)][1] = coords(node,1);
	}	
      }


      {
	std::vector< RCP<ofstream> > files; 
	for (std::size_t eq = 0; eq < num_eq; ++eq) {
	  std::stringstream ost;
	  ost << "upper_DOF" << eq << "_PID" << comm->MyPID() << ".dat";
	  files.push_back( rcp(new std::ofstream(ost.str().c_str()), 
			       ios::out | ios::trunc) );
	  files[eq]->precision(10);
	}
	
	const std::vector<int>& node_list = mb.topNodeSetGlobalIds();
	for (std::size_t node = 0; node < node_list.size(); ++node) {
	  int lid = owned_x->Map().LID(node_list[node] * num_eq);
	  for (std::size_t eq = 0; eq < num_eq; ++eq) {
	    int dof_index = lid + eq;
	    *(files[eq]) << coordinates[node_list[node]][0] << "   " 
			 << (*owned_x)[dof_index] << endl;
	  }
	}
      }

      {
	std::vector< RCP<ofstream> > files; 
	for (std::size_t eq = 0; eq < num_eq; ++eq) {
	  std::stringstream ost;
	  ost << "lower_DOF" << eq << "_PID" << comm->MyPID() << ".dat";
	  files.push_back( rcp(new std::ofstream(ost.str().c_str()), 
			       ios::out | ios::trunc) );
	  files[eq]->precision(10);
	}
	
	const std::vector<int>& node_list = mb.bottomNodeSetGlobalIds();
	for (std::size_t node = 0; node < node_list.size(); ++node) {
	  int lid = owned_x->Map().LID(node_list[node] * num_eq);
	  for (std::size_t eq = 0; eq < num_eq; ++eq) {
	    int dof_index = lid + eq;
	    *(files[eq]) << coordinates[node_list[node]][0] << "   " 
			 << (*owned_x)[dof_index] << endl;
	  }
	}
      }

    }

    TEST_FOR_EXCEPTION(!converged, std::runtime_error,
		       "Problem failed to converge!");

    TEST_FOR_EXCEPTION(num_newton_steps != 10, std::runtime_error,
		       "Incorrect number of Newton steps!");

    // Only check num gmres steps in serial
#ifndef HAVE_MPI
    TEST_FOR_EXCEPTION(num_gmres_iterations != 10, std::runtime_error,
		       "Incorrect number of GMRES iterations!");
#endif


    */

    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    if (comm->MyPID() == 0)
      std::cout << "\nRun has completed successfully!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
