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
#include "NOX_Epetra.H"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include <fstream>

#include "FiniteElementProblem.H"

// Basis Constructor
FiniteElementProblem::Basis::Basis() {
  phi = new double[2];
  dphide = new double[2];
}

// Destructor
FiniteElementProblem::Basis::~Basis() {
  delete [] phi;
  delete [] dphide;
}

// Calculates a linear 1D basis
void FiniteElementProblem::Basis::getBasis(int gp, double *x, double *u) {
  int N = 2;
  if (gp==0) {eta=-1.0/sqrt(3.0); wt=1.0;}
  if (gp==1) {eta=1.0/sqrt(3.0); wt=1.0;}

  // Calculate basis function and derivatives at nodel pts
  phi[0]=(1.0-eta)/2.0;
  phi[1]=(1.0+eta)/2.0;
  dphide[0]=-0.5;
  dphide[1]=0.5;
  
  // Caculate basis function and derivative at GP.
  dx=0.5*(x[1]-x[0]);
  xx=0.0;
  uu=0.0;
  duu=0.0;
  for (int i=0; i < N; i++) {
    xx += x[i] * phi[i];
    uu += u[i] * phi[i];
    duu += u[i] * dphide[i];
  }

  return;
}


// Constructor - creates the Epetra objects (maps and vectors) 
FiniteElementProblem::FiniteElementProblem(int numGlobalElements, Epetra_Comm& comm) :
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

  // Construct an overlaped map for the finite element fill *************
  // For single processor jobs, the overlap and standard map are the same
  if (NumProc == 1) {
    OverlapMap = new Epetra_Map(*StandardMap);
  } else {

    int OverlapNumMyElements;
    int OverlapMinMyGID;
    OverlapNumMyElements = NumMyElements + 2;
    if ((MyPID == 0) || (MyPID == NumProc - 1)) 
      OverlapNumMyElements --;
    
    if (MyPID==0) 
      OverlapMinMyGID = StandardMap->MinMyGID();
    else 
      OverlapMinMyGID = StandardMap->MinMyGID() - 1;
    
    int* OverlapMyGlobalElements = new int[OverlapNumMyElements];
    
    for (i = 0; i < OverlapNumMyElements; i ++) 
      OverlapMyGlobalElements[i] = OverlapMinMyGID + i;
    
    OverlapMap = new Epetra_Map(-1, OverlapNumMyElements, 
			    OverlapMyGlobalElements, 0, *Comm);

    delete [] OverlapMyGlobalElements;

  } // End Overlap map construction *************************************

  // Construct Linear Objects  
  Importer = new Epetra_Import(*OverlapMap, *StandardMap);
  initialSolution = new Epetra_Vector(*StandardMap);
  AA = new Epetra_CrsGraph(Copy, *StandardMap, 5);

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  generateGraph(*AA);

  // Create a second matrix using graph of first matrix - this creates a 
  // static graph so we can refill the new matirx after TransformToLocal()
  // is called.
  A = new Epetra_CrsMatrix (Copy, *AA);
  A->TransformToLocal();

  // Create the solver parameter list
  createSolverParameters();
}

// Destructor
FiniteElementProblem::~FiniteElementProblem()
{
  delete AA;
  delete A;
  delete initialSolution;
  delete Importer;
  delete OverlapMap;
  delete StandardMap;
}

void FiniteElementProblem::initializeSolution()
{
  // Initialize Solution
  initialSolution->PutScalar(1.0);
}

// Matrix and Residual Fills
bool FiniteElementProblem::evaluate(FillType f, 
				    const Epetra_Vector* soln, 
				    Epetra_Vector* tmp_rhs, 
				    Epetra_RowMatrix* tmp_matrix)
{
  flag = f;

  // Set the incoming linear objects
  if (flag == F_ONLY) {
    rhs = tmp_rhs;
  } else if (flag == MATRIX_ONLY) {
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } else if (flag == ALL) { 
    rhs = tmp_rhs;
    A = dynamic_cast<Epetra_CrsMatrix*> (tmp_matrix);
  } else {
    cout << "ERROR: FiniteElementProblem::fillMatrix() - FillType flag is broken" << endl;
    throw;
  }

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector x(*OverlapMap);

  // Export Solution to Overlap vector
  u.Import(*soln, *Importer, Insert);

  // Declare required variables
  int i,j,ierr;
  int OverlapNumMyElements = OverlapMap->NumMyElements();

  int OverlapMinMyGID;
  if (MyPID==0) OverlapMinMyGID = StandardMap->MinMyGID();
  else OverlapMinMyGID = StandardMap->MinMyGID()-1;

  int row, column;
  double factor=1000.0;
  double jac;
  double xx[2];
  double uu[2];
  Basis basis;

  // Create the nodal coordinates
  double Length=1.0;
  double dx=Length/((double) NumGlobalElements-1);
  for (i=0; i < OverlapNumMyElements; i++) {
    x[i]=dx*((double) OverlapMinMyGID+i);
  }
  
  // Zero out the objects that will be filled
  if ((flag == MATRIX_ONLY) || (flag == ALL)) i=A->PutScalar(0.0);
  if ((flag == F_ONLY)    || (flag == ALL)) i=rhs->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyElements-1; ne++) {
    
    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++) {
      // Get the solution and coordinates at the nodes 
      xx[0]=x[ne];
      xx[1]=x[ne+1];
      uu[0]=u[ne];
      uu[1]=u[ne+1];
      // Calculate the basis function at the gauss point
      basis.getBasis(gp, xx, uu);
	            
      // Loop over Nodes in Element
      for (i=0; i< 2; i++) {
	row=OverlapMap->GID(ne+i);
	//printf("Proc=%d GlobalRow=%d LocalRow=%d Owned=%d\n",
	//     MyPID, row, ne+i,StandardMap.MyGID(row));
	if (StandardMap->MyGID(row)) {
	  if ((flag == F_ONLY)    || (flag == ALL)) {
	    (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))]+=
	      +basis.wt*basis.dx
	      *((1.0/(basis.dx*basis.dx))*basis.duu*
		basis.dphide[i]+factor*basis.uu*basis.uu*basis.phi[i]);
	  }
	}
	// Loop over Trial Functions
	if ((flag == MATRIX_ONLY) || (flag == ALL)) {
	  for(j=0;j < 2; j++) {
	    if (StandardMap->MyGID(row)) {
	      column=OverlapMap->GID(ne+j);
	      jac=basis.wt*basis.dx*((1.0/(basis.dx*basis.dx))*
				     basis.dphide[j]*basis.dphide[i]
				     +2.0*factor*basis.uu*basis.phi[j]*
				     basis.phi[i]);  
	      ierr=A->SumIntoGlobalValues(row, 1, &jac, &column);
	    }
	  }
	}
      }
    }
  } 

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(0)=1
  if (MyPID==0) {
    if ((flag == F_ONLY)    || (flag == ALL)) 
      (*rhs)[0]= (*soln)[0] - 1.0;
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      column=0;
      jac=1.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column=1;
      jac=0.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
    }
  }

  // Sync up processors to be safe
  Comm->Barrier();
 
  A->TransformToLocal();

  return true;
}

bool FiniteElementProblem::computePrecMatrix(const Epetra_Vector& x, 
                                             Epetra_RowMatrix& M)
{
  // For now we will compute the entire Jacobian for the preconditioner
  Epetra_RowMatrix* precMatrix = dynamic_cast<Epetra_RowMatrix*>(&M);
  if (precMatrix == NULL) {
    cout << "ERROR: FiniteElementProblem::computePrecMatrix() - The supplied"
         << "Epetra_Operator is NOT an Epetra_RowMatrix!" << endl;
    throw;
  }
  return evaluate(MATRIX_ONLY, &x, NULL, precMatrix);
}
 
bool FiniteElementProblem::computePreconditioner(const Epetra_Vector& x, 
                                                 Epetra_Operator& M)
{
  // Here, we throw an exception, but this method could be expanded
  // in the future to provide an example of generic user-provided
  // preconditioning.
  cout << "ERROR: FiniteElementProblem::computePreconditioner"
       << " - Use Explicit Jacobian only for this test problem!" << endl;
  throw;
}

Epetra_Vector& FiniteElementProblem::getSolution()
{
  return *initialSolution;
}
  
Epetra_CrsMatrix& FiniteElementProblem::getJacobian()
{
  return *A;
}

NOX::Parameter::List& FiniteElementProblem::getParameters()
{
  return nlParams;
}

NOX::Parameter::List& FiniteElementProblem::getlsParameters()
{
  return lsParams;
}

Epetra_CrsGraph& FiniteElementProblem::generateGraph(Epetra_CrsGraph& AA)
{
  
  // Declare required variables
  int i,j;
  int row, column;
  int OverlapNumMyElements = OverlapMap->NumMyElements();
  int OverlapMinMyGID;
  if (MyPID==0) OverlapMinMyGID = StandardMap->MinMyGID();
  else OverlapMinMyGID = StandardMap->MinMyGID()-1;
  
  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyElements-1; ne++) {
          
    // Loop over Nodes in Element
    for (i=0; i< 2; i++) {
      row=OverlapMap->GID(ne+i);
      
      // Loop over Trial Functions
      for(j=0;j < 2; j++) {
	
	// If this row is owned by current processor, add the index
	if (StandardMap->MyGID(row)) {
	  column=OverlapMap->GID(ne+j);
	  AA.InsertGlobalIndices(row, 1, &column);
	}
      } 	
    }
  }
  AA.TransformToLocal();
//   AA.SortIndices();
//   AA.RemoveRedundantIndices();
  return AA;
}

void FiniteElementProblem::createSolverParameters()
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
  //dirParams.setParameter("Method", "Newton");
  //NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    //newtonParams.setParameter("Forcing Term Method", "Constant");
    //newtonParams.setParameter("Forcing Term Method", "Type 1");
    //newtonParams.setParameter("Forcing Term Method", "Type 2");
    //newtonParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-4);
    //newtonParams.setParameter("Forcing Term Maximum Tolerance", 0.1);
    //NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  dirParams.setParameter("Method", "Steepest Descent");
  NOX::Parameter::List& sdParams = dirParams.sublist("Steepest Descent");
    NOX::Parameter::List& lsParams = sdParams.sublist("Linear Solver");
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


void FiniteElementProblem::outputResults(NOX::Solver::Manager& solver, 
                   NOX::Parameter::List& printParams)
{
  // Output the parameter list
  NOX::Utils utils(printParams);
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << endl << "Final Parameters" << endl
	 << "****************" << endl;
    solver.getParameterList().print(cout);
    cout << endl;
  }

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>
        (finalGroup.getX())).getEpetraVector();

  // Print solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = finalSolution.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d",finalSolution.Map().Comm().MyPID());
  ifp = fopen(file_name, "w");
  for (int i=0; i<NumMyElements; i++)
    fprintf(ifp,"%d  %E\n",finalSolution.Map().MinMyGID()+i,finalSolution[i]);
  fclose(ifp);

}

