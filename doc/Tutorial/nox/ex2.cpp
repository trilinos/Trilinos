
// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// Trilinos Tutorial
// -----------------
// Simple nonlinear PDE problem.
// This file shows how to solve the nonlinear problem
//
// -\Delta u + \lambda e^u = 0  in \Omega = (0,1) \times (0,1)
//                       u = 0  on \partial \Omega
//
// using NOX 

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "NOX.H"
#include "NOX_Epetra_Interface.H"
#include "NOX_Epetra_Group.H"

// this is required to know the number of lower, upper, left and right
// node for each node of the Cartesian grid (composed by nx \timex ny 
// elements)

static void  get_neighbours( const int i, const int nx, const int ny,
			     int & left, int & right, 
			     int & lower, int & upper) 
{

  int ix, iy;
  ix = i%nx;
  iy = (i - ix)/nx;

  if( ix == 0 ) 
    left = -1;
  else 
    left = i-1;
  if( ix == nx-1 ) 
    right = -1;
  else
    right = i+1;
  if( iy == 0 ) 
    lower = -1;
  else
    lower = i-nx;
  if( iy == ny-1 ) 
    upper = -1;
  else
    upper = i+nx;

  return;

}

// This function creates a CrsMatrix, whose elements corresponds
// to the discretization of a Laplacian over a Cartesian grid,
// with nx grid point along the x-axis and and ny grid points 
// along the y-axis. For the sake of simplicity, I suppose that
// all the nodes in the matrix are internal nodes (Dirichlet
// boundary nodes are supposed to have been already condensated)

Epetra_CrsMatrix * CreateLaplacian( const int nx, const int ny,
				    const Epetra_Comm * Comm)
{

  int NumGlobalElements = nx * ny;
    
  // create a map
  Epetra_Map * Map = new Epetra_Map(NumGlobalElements,0,*Comm);
  // local number of rows
  int NumMyElements = Map->NumMyElements();
  // get update list
  int * MyGlobalElements = Map->MyGlobalElements();

  double hx = 1.0/(nx-1);
  double hy = 1.0/(ny-1);
  double off_left  = -1.0/(hx*hx);
  double off_right = -1.0/(hx*hx);
  double off_lower = -1.0/(hy*hy);
  double off_upper = -1.0/(hy*hy);
  double diag      =  2.0/(hx*hx) + 2.0/(hy*hy);
  
  int left, right, lower, upper;
    
  // a bit overestimated the nonzero per row
  
  Epetra_CrsMatrix * A = new Epetra_CrsMatrix(Copy,*Map,5);
    
  // Add  rows one-at-a-time
    
  double *Values = new double[4];
  int *Indices = new int[4];
  int NumEntries;
    
  for( int i=0 ; i<NumMyElements; ++i ) {
    int NumEntries=0;
    get_neighbours(  MyGlobalElements[i], nx, ny, 
		     left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = off_left;
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = off_right;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = off_lower;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = off_upper;
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
				 Values, Indices)==0);
    // Put in the diagonal entry
    assert(A->InsertGlobalValues(MyGlobalElements[i], 1, 
				 &diag, MyGlobalElements+i)==0);
  }

  // put matrix in local ordering
  A->TransformToLocal();

  return A;
  
} /* createJacobian */

// ==========================================================================
// This class contians the main definition of the nonlinear problem at
// hand. A method is provided to compute F(x) for a given x, and another
// method to update the entries of the Jacobian matrix, for a given x.
// As the Jacobian matrix J can be written as
//    J = L + diag(lambda*exp(x[i])),
// where L corresponds to the discretization of a Laplacian, and diag
// is a diagonal matrix with lambda*exp(x[i]). Basically, to update
// the jacobian we simply update the diagonal entries. Similarly, to compute
// F(x), we reset J to be equal to L, then we multiply it by the
// (distributed) vector x, then we add the diagonal contribution
// ==========================================================================

class PDEProblem {

public:

  // constructor. Requires the number of nodes along the x-axis
  // and y-axis, the value of lambda, and the Epetra_Communicator
  // (to define a Map, which is a linear map in this case)
  PDEProblem(const int nx, const int ny, const double lambda,
	     const Epetra_Comm * Comm) :
    nx_(nx), ny_(ny), lambda_(lambda) 
  {
    hx_ = 1.0/(nx_-1);
    hy_ = 1.0/(ny_-1);
    Matrix_ = CreateLaplacian(nx_,ny_,Comm);
  }

  // destructor
  ~PDEProblem() 
  {
    delete Matrix_;
  }

  // compute F(x)
  void ComputeF(const Epetra_Vector & x, Epetra_Vector & f) 
  {
    // reset diagonal entries
    double diag      =  2.0/(hx_*hx_) + 2.0/(hy_*hy_);
  
    int NumMyElements = Matrix_->Map().NumMyElements();
    // get update list
    int * MyGlobalElements = Matrix_->Map().MyGlobalElements( );
    
    int NumEntries;
    
    for( int i=0 ; i<NumMyElements; ++i ) {
      // Put in the diagonal entry
      assert(Matrix_->ReplaceGlobalValues(MyGlobalElements[i], 1, 
					  &diag, MyGlobalElements+i)==0);
    }
    // matrix-vector product (intra-processes communication occurs
    // in this call)
    Matrix_->Multiply(false,x,f);

    // add diagonal contributions
    for( int i=0 ; i<NumMyElements; ++i ) {
      // Put in the diagonal entry
      f[i] += lambda_*exp(x[i]);
    }
  }

  // update the Jacobian matrix for a given x
  void UpdateJacobian(const Epetra_Vector & x) 
  {
    
    double diag      =  2.0/(hx_*hx_) + 2.0/(hy_*hy_);
  
    int NumMyElements = Matrix_->Map().NumMyElements();
    // get update list
    int * MyGlobalElements = Matrix_->Map().MyGlobalElements( );
  
    int NumEntries;
    
    for( int i=0 ; i<NumMyElements; ++i ) {
      // Put in the diagonal entry
      double newdiag = diag + lambda_*x[i];
      Matrix_->ReplaceGlobalValues(MyGlobalElements[i], 1, 
				   &newdiag, MyGlobalElements+i);
    }

  }

  // returns a pointer to the internally stored matrix
  Epetra_CrsMatrix * GetMatrix()
  {
    return Matrix_;
  }
  
private:
  
  int nx_, ny_;
  double hx_, hy_;
  Epetra_CrsMatrix * Matrix_;
  double lambda_;

}; /* class PDEProblem */
 
// ==========================================================================
// This is the main NOX class for this example. Here we define
// the interface between the nonlinear problem at hand, and NOX.
// The constructor accepts a PDEProblem object. Using a pointer
// to this object, we can update the Jacobian and compute F(x),
// using the definition of our problem. This interface is bit
// crude: For instance, no PrecMatrix nor Preconditioner is specified.
// ==========================================================================

class SimpleProblemInterface : public NOX::Epetra::Interface {

public:
 
  //! Constructor
  SimpleProblemInterface( PDEProblem * Problem ) :
    Problem_(Problem) {};

  //! Destructor
  ~SimpleProblemInterface() 
  {
  };

  bool computeF(const Epetra_Vector & x, Epetra_Vector & f,
                NOX::Epetra::Interface::FillType F )
  {
    Problem_->ComputeF(x,f);
    return true;
  };
  
  bool computeJacobian(const Epetra_Vector & x, Epetra_Operator & Jac)
  {

    Problem_->UpdateJacobian(x);
    return true;
  }

  bool computePrecMatrix(const Epetra_Vector & x, Epetra_RowMatrix & M) 
  {
    cout << "*ERR* SimpleProblem::preconditionVector()\n";
    cout << "*ERR* don't use explicit preconditioning" << endl;
    throw 1;
  }  
  
  bool computePreconditioner(const Epetra_Vector & x, Epetra_Operator & O)
  {
    cout << "*ERR* SimpleProblem::preconditionVector()\n";
    cout << "*ERR* don't use explicit preconditioning" << endl;
    throw 1;
  }  

private:
  
  PDEProblem * Problem_;
  
}; /* class SimpleProblemInterface */

// =========== //
// main driver //
// =========== //

int main( int argc, char **argv )
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // define the parameters of the nonlinear PDE problem
  int nx = 5;
  int ny = 6;  
  double lambda = 1.0;

  PDEProblem Problem(nx,ny,lambda,&Comm);
 
  // starting solution, here a zero vector
  Epetra_Vector InitialGuess(Problem.GetMatrix()->Map());
  InitialGuess.PutScalar(0.0);

  // Set up the problem interface
  SimpleProblemInterface Interface(&Problem);
  
  // Create the top level parameter list
  NOX::Parameter::List nlParams;

  // Set the nonlinear solver method
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", Comm.MyPID()); 
  printParams.setParameter("Output Precision", 3);
  printParams.setParameter("Output Processor", 0);
  printParams.setParameter("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);

  // start definition of nonlinear solver parameters
  // Sublist for line search 
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  searchParams.setParameter("Method", "Full Step");

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  dirParams.setParameter("Method", "Newton");

  NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
  newtonParams.setParameter("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "GMRES");  
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Output Frequency", 50);    
  lsParams.setParameter("Aztec Preconditioner", "ilu"); 

  NOX::Epetra::Group grp(printParams, lsParams, Interface, InitialGuess, 
                         dynamic_cast<Epetra_RowMatrix&>(*(Problem.GetMatrix()))); 

  // Set up the status tests
  NOX::StatusTest::NormF statusTestA(grp, 1.0e-4);
  NOX::StatusTest::MaxIters statusTestB(20);
  // this will be the convergence test to be used
  NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

  // Create the list of solver parameters
  NOX::Parameter::List solverParameters;

  // Set the level of output (this is the default)
  solverParameters.sublist("Printing").setParameter("Output Information", 
						    NOX::Utils::OuterIteration);

  // Set the solver (this is the default)
  solverParameters.setParameter("Nonlinear Solver", "Line Search Based");

  // Create the line search parameters sublist
  NOX::Parameter::List& lineSearchParameters = solverParameters.sublist("Line Search");

  // Set the line search method
  lineSearchParameters.setParameter("Method","More'-Thuente");

  // Create the solver
  NOX::Solver::Manager solver(grp, statusTestsCombo, solverParameters);

  // Solve the nonlinesar system
  NOX::StatusTest::StatusType status = solver.solve();

  // Print the answer
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver.getParameterList().print(cout);

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group & finalGroup = 
    dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector & finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  if( Comm.MyPID() == 0 ) cout << "Computed solution : " << endl;
  cout << finalSolution;

} /* main */
