//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                                                
#include "NOX_Common.H"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "DennisSchnabel.H"

// Constructor - creates the Epetra objects (maps and vectors) 
DennisSchnabel::DennisSchnabel(int numGlobalElements, Epetra_Comm& comm) :
  GenericProblem(numGlobalElements, comm),
  Comm(&comm),
  NumGlobalElements(numGlobalElements)
{

  // Commonly used variables
  int i;
  MyPID = Comm->MyPID();      // Process ID
  NumProc = Comm->NumProc();  // Total number of processes

  // Construct a Source Map that puts approximately the same 
  // Number of equations on each processor in uniform global ordering
  StandardMap = new Epetra_Map(NumGlobalElements, 0, *Comm);

  // Get the number of elements owned by this processor
  NumMyElements = StandardMap->NumMyElements();

  // Construct an overlaped map for the fill calls **********************
  /* The overlap map is needed for multiprocessor jobs.  The unknowns 
   * are stored in a distributed vector where each node owns one unknown.  
   * During a function or Jacobian evaluation, each processor needs to have 
   * both of the unknown values.  The overlapped vector can get this data 
   * by importing the owned values from the other node to this overlapped map. 
   * Actual solves must be done using the Standard map where everything is 
   * distributed.
   */
  // For single processor jobs, the overlap and standard map are the same
  if (NumProc == 1) {
    OverlapMap = new Epetra_Map(*StandardMap);
  } 
  else {

    int OverlapNumMyElements = 2;
    int OverlapMyGlobalElements[2];
    
    for (i = 0; i < OverlapNumMyElements; i ++) 
      OverlapMyGlobalElements[i] = i;
    
    OverlapMap = new Epetra_Map(-1, OverlapNumMyElements, 
				OverlapMyGlobalElements, 0, *Comm);
  } // End Overlap map construction *************************************

  // Construct Linear Objects  
  Importer = new Epetra_Import(*OverlapMap, *StandardMap);
  initialSolution = new Epetra_Vector(*StandardMap);
  AA = new Epetra_CrsGraph(Copy, *StandardMap, 5);

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  generateGraph(*AA);

  // Use the graph AA to create a Matrix.
  A = new Epetra_CrsMatrix (Copy, *AA);

  // Transform the global matrix coordinates to local so the matrix can 
  // be operated upon.
  A->TransformToLocal();

  // Create the solver parameter list
  createSolverParameters();
}

// Destructor
DennisSchnabel::~DennisSchnabel()
{
  delete AA;
  delete A;
  delete initialSolution;
  delete Importer;
  delete OverlapMap;
  delete StandardMap;
}

void DennisSchnabel::initializeSolution()
{
  if (MyPID==0) {
    (*initialSolution)[0]=2.0;
    if (NumProc==1) 
      (*initialSolution)[1]=0.5;
  } 
  else 
    (*initialSolution)[0]=0.5;
}

// Matrix and Residual Fills
bool DennisSchnabel::evaluate(FillType f, 
			      const Epetra_Vector* soln, 
			      Epetra_Vector* tmp_rhs, 
			      Epetra_RowMatrix* tmp_matrix)
{
  flag = f;

  // Set the incoming linear objects
  if (flag == F_ONLY) {
    rhs = tmp_rhs;
  } 
  else if (flag == MATRIX_ONLY) {
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } 
  else if (flag == ALL) { 
    rhs = tmp_rhs;
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } 
  else {
    cout << "ERROR: DennisSchnabel::fillMatrix() - No such flag as " 
	 << flag << endl;
    throw;
  }

  // Create the overlapped solution
  Epetra_Vector u(*OverlapMap);

  // Export Solution to Overlap vector so we have all unknowns required
  // for function and Jacobian evaluations.
  u.Import(*soln, *Importer, Insert);

  // Declare required variables
  int i,ierr;
  int OverlapMinMyGID;
  if (MyPID==0) OverlapMinMyGID = StandardMap->MinMyGID();
  else OverlapMinMyGID = StandardMap->MinMyGID()-1;

  // Begin F fill
  if((flag == F_ONLY) || (flag == ALL)) {

    // Zero out the F vector
    i=rhs->PutScalar(0.0);

    // Processor 0 always fills the first equation.
    if (MyPID==0) { 
      (*rhs)[0]=(u[0]*u[0] + u[1]*u[1] - 2.);

      // If it's a single processor job, fill the second equation on proc 0.
      if (NumProc==1) 
	(*rhs)[1]=(exp(u[0]-1.) + u[1]*u[1]*u[1] - 2.);
    } 
    // Multiprocessor job puts the second equation on processor 1.
    else { 
      (*rhs)[0]=(exp(u[0]-1.) + u[1]*u[1]*u[1] - 2.);
    }
  }

  
  int* column = new int[2];
  double* jac = new double[2];

  // The matrix is 2 x 2 and will always be 0 and 1 regardless of 
  // the coordinates being local or global.
  column[0] = 0; 
  column[1] = 1;

  // Begin Jacobian fill
  if((flag == MATRIX_ONLY) || (flag == ALL)) {

    // Zero out Jacobian
    i=A->PutScalar(0.0);

    if (MyPID==0) {
      // Processor 0 always fills the first equation.
      jac[0] = 2.*u[0];
      jac[1] = 2.*u[1];
      ierr=A->ReplaceGlobalValues(0, 2, jac, column);

      // If it's a single processor job, fill the second equation on proc 0.
      if (NumProc==1) {	 
	jac[0] = exp(u[0]-1.);
	jac[1] = 3.*u[1]*u[1];
	ierr=A->ReplaceGlobalValues(1, 2, jac, column);
      }
    } 
    // Multiprocessor job puts the second equation on processor 1.
    else {
      jac[0] = exp(u[0]-1.);
      jac[1] = 3.*u[1]*u[1];
      ierr=A->ReplaceGlobalValues(1, 2, jac, column);
    }
  } 

  delete [] column;
  delete [] jac;

  // Sync up processors to be safe
  Comm->Barrier();
 
  // Transform matrix so it can be operated upon.
  A->TransformToLocal();

  return true;
}

bool DennisSchnabel::computePrecMatrix(const Epetra_Vector& x, 
                                             Epetra_RowMatrix& M)
{
  // For now we will compute the entire Jacobian for the preconditioner
  Epetra_RowMatrix* precMatrix = dynamic_cast<Epetra_RowMatrix*>(&M);
  if (precMatrix == NULL) {
    cout << "ERROR: DennisSchnabel::computePrecMatrix() - The supplied"
         << "Epetra_Operator is NOT an Epetra_RowMatrix!" << endl;
    throw;
  }
  return evaluate(MATRIX_ONLY, &x, NULL, precMatrix);
}
 
bool DennisSchnabel::computePreconditioner(const Epetra_Vector& x, 
                                                 Epetra_Operator& M)
{
  // Here, we throw an exception, but this method could be expanded
  // in the future to provide an example of generic user-provided
  // preconditioning.
  cout << "ERROR: DennisSchnabel::computePreconditioner"
       << " - Use Explicit Jacobian only for this test problem!" << endl;
  throw;
}

Epetra_Vector& DennisSchnabel::getSolution()
{
  return *initialSolution;
}
  
Epetra_CrsMatrix& DennisSchnabel::getJacobian()
{
  return *A;
}

NOX::Parameter::List& DennisSchnabel::getParameters()
{
  return nlParams;
}

NOX::Parameter::List& DennisSchnabel::getlsParameters()
{
  return lsParams;
}

Epetra_CrsGraph& DennisSchnabel::generateGraph(Epetra_CrsGraph& AA)
{
  
  int* index = new int[2];

  if (MyPID==0) {
    index[0]=0;
    index[1]=1;
    AA.InsertGlobalIndices(0, 2, index);
  
    if (NumProc==1) {
      index[0]=0;
      index[1]=1;
      AA.InsertGlobalIndices(1, 2, index);
    }
  } else {
    index[0]=0;
    index[1]=1;
    AA.InsertGlobalIndices(1, 2, index);
  }
  
  delete [] index;
  
  AA.TransformToLocal();
//   AA.SortIndices();
//   AA.RemoveRedundantIndices();
  return AA;
}

void DennisSchnabel::createSolverParameters()
{

  // Create the top level parameter list

  // Set the nonlinear solver method
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");
  //nlParams.setParameter("Nonlinear Solver", "Trust Region Based");

  // Set the printing parameters in the "Printing" sublist
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", MyPID); 
  printParams.setParameter("Output Precision", 3);
  printParams.setParameter("Output Processor", 0);
  printParams.setParameter("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);

  // Sublist for line search 
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  //searchParams.setParameter("Method", "Full Step");
  //searchParams.setParameter("Method", "Interval Halving");
  searchParams.setParameter("Method", "Polynomial");
  //searchParams.setParameter("Method", "NonlinearCG");
  //searchParams.setParameter("Method", "Quadratic");
  //searchParams.setParameter("Method", "More'-Thuente");

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  dirParams.setParameter("Method", "Newton");
  NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    newtonParams.setParameter("Forcing Term Method", "Constant");
    //newtonParams.setParameter("Forcing Term Method", "Type 1");
    //newtonParams.setParameter("Forcing Term Method", "Type 2");
    //newtonParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-4);
    //newtonParams.setParameter("Forcing Term Maximum Tolerance", 0.1);
    //NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  //dirParams.setParameter("Method", "Steepest Descent");
  //NOX::Parameter::List& sdParams = dirParams.sublist("Steepest Descent");
    //NOX::Parameter::List& lsParams = sdParams.sublist("Linear Solver");
    //sdParams.setParameter("Scaling Type", "None");
    //sdParams.setParameter("Scaling Type", "2-Norm");
    //sdParams.setParameter("Scaling Type", "Quadratic Model Min");
  //dirParams.setParameter("Method", "NonlinearCG");
  //NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
    //nlcgParams.setParameter("Restart Frequency", 2000);
    //nlcgParams.setParameter("Precondition", "On");
    //nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
    //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

  // Sublist for linear solver
  //lsParams.setParameter("Aztec Solver", "GMRES");  
  //lsParams.setParameter("Max Iterations", 800);  
  //lsParams.setParameter("Tolerance", 1e-4);
  //lsParams.setParameter("Output Frequency", 50);    
  //lsParams.setParameter("Preconditioning", "None");   
  //lsParams.setParameter("Preconditioning", "AztecOO: Jacobian Matrix");   
  //lsParams.setParameter("Preconditioning", "AztecOO: User RowMatrix"); 
  //lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");
  //lsParams.setParameter("Aztec Preconditioner", "ilu"); 
  //lsParams.setParameter("Overlap", 2);  
  //lsParams.setParameter("Graph Fill", 2); 
  //lsParams.setParameter("Aztec Preconditioner", "ilut"); 
  //lsParams.setParameter("Overlap", 2);   
  //lsParams.setParameter("Fill Factor", 2);   
  //lsParams.setParameter("Drop Tolerance", 1.0e-12);   
  //lsParams.setParameter("Aztec Preconditioner", "Polynomial"); 
  //lsParams.setParameter("Polynomial Order", 6); 
}
