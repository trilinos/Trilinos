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
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"

// User defined objects
#include "Element_Linear2D.hpp"
#include "Workset.hpp"
#include "Traits.hpp"
#include "FactoryTraits.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_VbrRowMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Util.h"

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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

#ifdef HAVE_MPI  
  MPI_Init(&argc, &argv);
#endif

  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    bool print_debug_info = true;

    cout << "\nStarting Epetra_VBR_Test Example!\n" << endl;


    // *********************************************************
    // * Build the Finite Element data structures
    // *********************************************************

    // Create the mesh, one strip of 2D elements.
    const std::size_t num_local_cells = 5;
    
    double domain_length = 1.0;
    double dx = domain_length / static_cast<double>(num_local_cells);
    std::vector<int> global_id(4);
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
      
      Element_Linear2D e(global_id, i, i, x_coords, y_coords);
      cells.push_back(e);
      
      // update global indices for next element
      for (std::size_t i=0; i < global_id.size(); ++i)
	global_id[i] += 2;
      
    }
    
    // Divide mesh into workset blocks
    const std::size_t workset_size = 5;
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
    RCP<Epetra_VbrRowMatrix> Jac;
    {
      Epetra_SerialComm comm;
      Epetra_BlockMap map(static_cast<int>(num_nodes), static_cast<int>(num_eq), 0, comm);
      Epetra_DataAccess copy = ::Copy;
      Epetra_CrsGraph graph(copy, map, 3);

      std::vector<Element_Linear2D>::iterator e = cells.begin();
      for (; e != cells.end(); ++e) {
	for (int row = 0; row < e->numNodes(); ++row) {
	  for (int col = 0; col < e->numNodes(); ++col) {
	    int global_row = e->globalNodeId(row);
	    int global_col = e->globalNodeId(col);
	    graph.InsertGlobalIndices(global_row, 1, &global_col);
	  }
	}
      }
      graph.FillComplete();
      graph.Print(cout);
      
      Epetra_SerialDenseMatrix block_matrix(2,2);
      Epetra_SerialDenseMatrix diag_block_matrix(2,2);

      RCP<Epetra_VbrMatrix> Jac_vbr = rcp(new Epetra_VbrMatrix(copy,graph));

      Epetra_Util util;

      e = cells.begin();
      for (; e != cells.end(); ++e) {
	for (int row = 0; row < e->numNodes(); ++row) {
	  
	  int global_row = e->globalNodeId(row);
	  
	  block_matrix(0,0) = util.RandomDouble();
	  block_matrix(0,1) = util.RandomDouble();
	  block_matrix(1,0) = util.RandomDouble();
	  block_matrix(1,1) = util.RandomDouble();

	  diag_block_matrix(0,0) = 100.0*util.RandomDouble();
	  diag_block_matrix(0,1) = util.RandomDouble();
	  diag_block_matrix(1,0) = util.RandomDouble();
	  diag_block_matrix(1,1) = 100.0*util.RandomDouble();
	  
	  for (int col = 0; col < e->numNodes(); ++col) {
	    int global_col = e->globalNodeId(col);
	    Jac_vbr->BeginReplaceMyValues(global_row, 1, &global_col);
	    if (global_row==global_col) 
	      Jac_vbr->SubmitBlockEntry(diag_block_matrix);
	    else
	      Jac_vbr->SubmitBlockEntry(block_matrix);
	    Jac_vbr->EndSubmitEntries();
	  }
	}
      }
      Jac_vbr->FillComplete();
      x = rcp(new Epetra_Vector(map));
      delta_x = rcp(new Epetra_Vector(map));
      f = rcp(new Epetra_Vector(map));
      
      x->PutScalar(1.0);
      Jac_vbr->Apply(*x,*f);

      Jac = 
	rcpWithEmbeddedObjPostDestroy(new Epetra_VbrRowMatrix(Jac_vbr.get()),
				      Jac_vbr);

    }

    if (print_debug_info) {
      x->Print(cout);
      Jac->Print(cout);
      f->Print(cout); 
    }

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
    IFPACK_CHK_ERR(Prec->Compute());
    IFPACK_CHK_ERR(Prec->Condest());
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

    RCP<Time> linear_solve_time = 
      TimeMonitor::getNewTimer("Linear Solve Time");

    
    std::size_t num_gmres_iterations = 0;

    {
      TimeMonitor t(*linear_solve_time);
      
      delta_x->PutScalar(0.0);
      
      IFPACK_CHK_ERR(Prec->Compute());
      
      problem->setProblem();
      
      Belos::ReturnType ret = gmres_solver->solve();
      
      int num_iters = gmres_solver->getNumIters();
      num_gmres_iterations += num_iters; 
      if (print_debug_info)
	std::cout << "Number of gmres iterations performed for this solve: " 
		  << num_iters << std::endl;
      
      if (ret!=Belos::Converged) {
	std::cout << std::endl << "WARNING:  Belos did not converge!" 
		  << std::endl;
      }
      
    }
      
    delta_x->Print(cout);

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

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
