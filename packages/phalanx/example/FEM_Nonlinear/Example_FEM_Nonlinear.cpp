// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

// User defined objects
#include "Element_Linear2D.hpp"
#include "Workset.hpp"
#include "Traits.hpp"
#include "FactoryTraits.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

// Linear solver
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// Preconditioner
#include "Ifpack.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void checkConvergence(const int num_newton_steps, const Epetra_Vector& f, 
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

  cout << "Step: " << num_newton_steps << ", ||f|| = "  << norm_f
       << ", ||delta x|| = " << norm_delta_x << endl;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void applyBoundaryConditions(const double& T_left, const double& Vx_left,
			     const Epetra_Vector& x, Epetra_CrsMatrix& Jac, 
			     Epetra_Vector& f)
{
  // This function assumes a serial run.  It will not check process rank.
  // Assumes two equations.

  // Left edge includes global nodes 0 and 1, dofs 0,1,2,3
  int num_nonzero_entries_in_row = 0;
  int* row_indices;
  double* row_values;

  // Temp @ node 0
  Jac.ExtractMyRowView(0, num_nonzero_entries_in_row, row_values, row_indices);
  for (int i=0; i < num_nonzero_entries_in_row; ++i) {
    if (row_indices[i] == 0)
      row_values[i] = 1.0;
    else
      row_values[i] = 0.0;
  }
  f[0] = T_left - x[0];

  // Vx @ node 0
  Jac.ExtractMyRowView(1, num_nonzero_entries_in_row, row_values, row_indices);
  for (int i=0; i < num_nonzero_entries_in_row; ++i) {
    if (row_indices[i] == 1)
      row_values[i] = 1.0;
    else
      row_values[i] = 0.0;
  }
  f[1] = Vx_left - x[1];

  // Temp @ node 1
  Jac.ExtractMyRowView(2, num_nonzero_entries_in_row, row_values, row_indices);
  for (int i=0; i < num_nonzero_entries_in_row; ++i) {
    if (row_indices[i] == 2)
      row_values[i] = 1.0;
    else
      row_values[i] = 0.0;
  }
  f[2] = T_left - x[2];

  // Vx @ node 1
  Jac.ExtractMyRowView(3, num_nonzero_entries_in_row, row_values, row_indices);
  for (int i=0; i < num_nonzero_entries_in_row; ++i) {
    if (row_indices[i] == 3)
      row_values[i] = 1.0;
    else
      row_values[i] = 0.0;
  }
  f[3] = Vx_left - x[3];

}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    bool print_debug_info = false;

    
    cout << "\nStarting FEM_Nonlinear Example!\n" << endl;


    // *********************************************************
    // * Build the Finite Element data structures
    // *********************************************************

    // Problem dimension - a 2D problem
    const static int dim = 2;

    // Create the mesh, one strip of 2D elements.
    const std::size_t num_local_cells = 1005;
    
    double domain_length = 1.0;
    double dx = domain_length / static_cast<double>(num_local_cells);
    std::vector<unsigned> global_id(4);
    global_id[0] = 0;
    global_id[1] = 2;
    global_id[2] = 3;
    global_id[3] = 1;
    std::vector<double> x_coords(4);
    std::vector<double> y_coords(4);
    std::vector<Element_Linear2D> cells;
    for (std::size_t i = 0; i < num_local_cells; ++i) {
      
      x_coords[0] = static_cast<double>(i) * dx;
      x_coords[1] = x_coords[0] + dx;
      x_coords[2] = x_coords[0] + dx;
      x_coords[3] = static_cast<double>(i) * dx;
      y_coords[0] = 0.0;
      y_coords[1] = 0.0;
      y_coords[2] = 1.0;
      y_coords[3] = 1.0;
      
      Element_Linear2D e(global_id, i, x_coords, y_coords);
      cells.push_back(e);
      
      // update global indices for next element
      for (std::size_t i=0; i < global_id.size(); ++i)
	global_id[i] += 2;
      
    }
    
    // Divide mesh into workset blocks
    const std::size_t workset_size = 50;
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
    const std::size_t num_nodes = 2 * (num_local_cells +1);
    const std::size_t num_dof = num_nodes * num_eq;
    RCP<Epetra_Vector> x;
    RCP<Epetra_Vector> delta_x;
    RCP<Epetra_Vector> f;
    RCP<Epetra_CrsMatrix> Jac;
    {
      Epetra_SerialComm comm;
      Epetra_Map map(num_dof, num_dof, 0, comm);      
      x = rcp(new Epetra_Vector(map, true));         
      delta_x = rcp(new Epetra_Vector(map, true)); 
      f = rcp(new Epetra_Vector(map, true));

      Epetra_DataAccess d = ::Copy;
      const std::size_t approximate_indices_per_row = 6 * num_eq;
      Epetra_CrsGraph graph(d, map, approximate_indices_per_row);
      std::vector<Element_Linear2D>::iterator cell_it = cells.begin();
      for (; cell_it != cells.end(); ++cell_it) {

	  std::vector<int> col_indices(0);
	  for (std::size_t node = 0; node < cell_it->numNodes(); ++node)
	    for (std::size_t eq = 0; eq < num_eq; ++eq)
	      col_indices.push_back(cell_it->globalNodeId(node)*num_eq + eq);
	  
	  for (std::size_t node = 0; node < cell_it->numNodes(); ++node)
	    for (std::size_t eq = 0; eq < num_eq; ++eq)
	      graph.InsertGlobalIndices(cell_it->globalNodeId(node)*num_eq+eq, 
					col_indices.size(), &col_indices[0]);

      }
      graph.FillComplete();
      Jac = rcp(new Epetra_CrsMatrix(d,graph));
      
    }

    //if (print_debug_info) {
    if (false) {
      x->Print(cout);
      Jac->Print(cout);
    }

    // Sets bc for initial guess
    applyBoundaryConditions(1.0, 1.0, *x, *Jac, *f);
    
    // *********************************************************
    // * Build the FieldManager
    // *********************************************************

    RCP< vector<string> > dof_names = rcp(new vector<string>(num_eq));
    (*dof_names)[0] = "Temperature";
    (*dof_names)[1] = "Velocity X";

    RCP<DataLayout> qp_scalar = rcp(new MDALayout<QuadPoint>(4));
    RCP<DataLayout> node_scalar = rcp(new MDALayout<Node>(4));
    
    RCP<DataLayout> qp_vec = rcp(new MDALayout<QuadPoint,Dim>(4,dim));
    RCP<DataLayout> node_vec = rcp(new MDALayout<Node,Dim>(4,dim));

    RCP<DataLayout> dummy = rcp(new FlatLayout("Dummy",0));
    
    map<string, RCP<ParameterList> > evaluators_to_build;
    
    { // Gather Solution
      RCP<ParameterList> p = rcp(new ParameterList);
      int type = MyFactoryTraits<MyTraits>::id_gather_solution;
      p->set<int>("Type", type);
      p->set< RCP< vector<string> > >("Solution Names", dof_names);
      p->set< RCP<Epetra_Vector> >("Solution Vector", x);
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
      p->set< RCP<Epetra_Vector> >("Solution Vector", x);
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
      p->set< RCP<Epetra_Vector> >("Residual Vector", f);
      p->set< RCP<Epetra_CrsMatrix> >("Jacobian Matrix", Jac);
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
    }

    {
      RCP<Time> registration_time = 
	TimeMonitor::getNewTimer("Post Registration Setup Time");
      {
	TimeMonitor t(*registration_time);
	fm.postRegistrationSetup(workset_size);
      }
    }

    if (print_debug_info)
      cout << fm << endl;

    // *********************************************************
    // * Build Preconditioner (Ifpack)
    // *********************************************************
    
    Ifpack Factory;
    std::string PrecType = "ILU";
    int OverlapLevel = 1;
    RCP<Ifpack_Preconditioner> Prec = 
      Teuchos::rcp( Factory.Create(PrecType, &*Jac, OverlapLevel) );
    ParameterList ifpackList;
    ifpackList.set("fact: drop tolerance", 1e-9);
    ifpackList.set("fact: level-of-fill", 1);
    ifpackList.set("schwarz: combine mode", "Add");
    IFPACK_CHK_ERR(Prec->SetParameters(ifpackList));
    IFPACK_CHK_ERR(Prec->Initialize());
    RCP<Belos::EpetraPrecOp> belosPrec = 
      rcp( new Belos::EpetraPrecOp( Prec ) );

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
    belosList->set<int>("Num Blocks", num_dof);
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
      Teuchos::rcp_implicit_cast<Epetra_MultiVector>(f);
    
    RCP<Epetra_MultiVector> DX = 
      Teuchos::rcp_implicit_cast<Epetra_MultiVector>(delta_x);
    
    RCP<Belos::LinearProblem<double,MV,OP> > problem =
      rcp(new Belos::LinearProblem<double,MV,OP>(Jac, DX, F) );
    
    problem->setRightPrec( belosPrec );

    RCP< Belos::SolverManager<double,MV,OP> > gmres_solver = 
      rcp( new Belos::BlockGmresSolMgr<double,MV,OP>(problem, belosList) );
    
    // *********************************************************
    // * Solve the system
    // *********************************************************

    RCP<Time> residual_eval_time = 
      TimeMonitor::getNewTimer("Residual Evaluation Time");
    RCP<Time> jacobian_eval_time = 
      TimeMonitor::getNewTimer("Jacobian Evaluation Time");
    RCP<Time> linear_solve_time = 
      TimeMonitor::getNewTimer("Linear Solve Time");
    RCP<Time> nonlinear_solve_time = 
      TimeMonitor::getNewTimer("Nonlinear Solve Time");

    // Set initial guess
    x->PutScalar(1.0);

    // Evaluate Residual
    {
      TimeMonitor t(*residual_eval_time);

      f->PutScalar(0.0);

      for (std::size_t i = 0; i < worksets.size(); ++i)
	fm.evaluateFields<MyTraits::Residual>(worksets[i]);
     
      applyBoundaryConditions(1.0, 1.0, *x, *Jac, *f);
    }
    
    if (print_debug_info) {
      cout << "Printing Initial Residual" << endl;
      f->Print(cout);
      cout << endl;
    }

    // Newton Loop
    bool converged = false;
    std::size_t num_newton_steps = 0;
    std::size_t num_gmres_iterations = 0;
    checkConvergence(num_newton_steps, *f, *delta_x, converged);

    while (!converged && num_newton_steps < 21) {
      
      TimeMonitor t(*nonlinear_solve_time);

      // Evaluate Residual and Jacobian
      {
	TimeMonitor t(*jacobian_eval_time);

	f->PutScalar(0.0);
	Jac->PutScalar(0.0);

	for (std::size_t i = 0; i < worksets.size(); ++i)
	  fm.evaluateFields<MyTraits::Jacobian>(worksets[i]);

	applyBoundaryConditions(1.0, 1.0, *x, *Jac, *f);
      }

      if (print_debug_info) {
	cout << "Residual AND Jacobian for Newton step: " << num_newton_steps
	     << endl;
	f->Print(cout);
	Jac->Print(cout);
	cout << endl;
      }
      
      f->Scale(-1.0);

      // Solve linear problem
      {
	TimeMonitor t(*linear_solve_time);
	
	delta_x->PutScalar(0.0);

	IFPACK_CHK_ERR(Prec->Compute());

	problem->setProblem();

	Belos::ReturnType ret = gmres_solver->solve();

	int num_iters = gmres_solver->getNumIters();
	num_gmres_iterations += num_iters; 
	//if (print_debug_info) 
	if (true)
	  std::cout << "Number of gmres iterations performed for this solve: " 
		    << num_iters << std::endl;
	
	if (ret!=Belos::Converged) {
	  std::cout << std::endl << "WARNING:  Belos did not converge!" 
		    << std::endl;
	}
	
      }
      
      x->Update(1.0, *delta_x, 1.0);

      { // Evaluate Residual Only
	TimeMonitor t(*residual_eval_time);
	
	f->PutScalar(0.0);

	for (std::size_t i = 0; i < worksets.size(); ++i)
	  fm.evaluateFields<MyTraits::Residual>(worksets[i]);
	
	applyBoundaryConditions(1.0, 1.0, *x, *Jac, *f);
      }

      num_newton_steps += 1;

      checkConvergence(num_newton_steps, *f, *delta_x, converged);
    }
     
    RCP<Time> file_io = 
      TimeMonitor::getNewTimer("File IO");
    {
      TimeMonitor t(*file_io);
      
      ofstream file_upper_temp("upper_temp.dat", ios::out | ios::trunc);
      ofstream file_upper_vx("upper_vx.dat", ios::out | ios::trunc);
      ofstream file_lower_temp("lower_temp.dat", ios::out | ios::trunc);
      ofstream file_lower_vx("lower_vx.dat", ios::out | ios::trunc);
      
      for (std::size_t i = 0; i < worksets.size(); ++i) {
	std::vector<Element_Linear2D>::iterator element = worksets[i].begin;
	for (; element != worksets[i].end; ++element) {
	  const PHX::Array<double,NaturalOrder,Node,Dim>& coords = 
	    element->nodeCoordinates();
	  
	  int first_DOF = element->globalNodeId(0) * num_eq;
	  file_lower_temp << coords(0,0) << "   " << (*x)[first_DOF] << endl;
	  file_lower_vx   << coords(0,0) << "   " << (*x)[first_DOF+1] << endl;
	  first_DOF = element->globalNodeId(1) * num_eq;
	  file_lower_temp << coords(1,0) << "   " << (*x)[first_DOF] << endl;
	  file_lower_vx   << coords(1,0) << "   " << (*x)[first_DOF+1] << endl;
	  first_DOF = element->globalNodeId(3) * num_eq;
	  file_upper_temp << coords(3,0) << "   " << (*x)[first_DOF] << endl;
	  file_upper_vx   << coords(3,0) << "   " << (*x)[first_DOF+1] << endl;
	  first_DOF = element->globalNodeId(2) * num_eq;
	  file_upper_temp << coords(2,0) << "   " << (*x)[first_DOF] << endl;
	  file_upper_vx   << coords(2,0) << "   " << (*x)[first_DOF+1] << endl;
	}
      }
    }

    TEST_FOR_EXCEPTION(!converged, std::runtime_error,
		       "Problem failed to converge!");

    TEST_FOR_EXCEPTION(num_newton_steps != 4, std::runtime_error,
		       "Incorrect number of Newton steps!");

    TEST_FOR_EXCEPTION(num_gmres_iterations != 4, std::runtime_error,
		       "Incorrect number of GMRES iterations!");


    // *********************************************************************
    // Finished all testing
    // *********************************************************************
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
    
  return 0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
