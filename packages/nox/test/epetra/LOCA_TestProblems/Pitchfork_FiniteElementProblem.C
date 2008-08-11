// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Common.H"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Basis.H"

#include "Pitchfork_FiniteElementProblem.H"

// Constructor - creates the Epetra objects (maps and vectors) 
Pitchfork_FiniteElementProblem::Pitchfork_FiniteElementProblem(
						       int numGlobalElements, 
						       Epetra_Comm& comm) :
  flag(F_ONLY),
  StandardMap(NULL),
  OverlapMap(NULL),
  Importer(NULL),
  initialSolution(NULL),
  rhs(NULL),
  AA(NULL),
  A(NULL),
  Comm(&comm),
  MyPID(0),
  NumProc(0),
  NumMyElements(0),
  NumGlobalElements(numGlobalElements),
  lambda(0.0),
  alpha(0.0),
  beta(0.0)
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
  // static graph so we can refill the new matirx after FillComplete()
  // is called.
  A = new Epetra_CrsMatrix (Copy, *AA);
  A->FillComplete();

  // Set default bifurcation values
  lambda = -2.25;
  alpha = 1.0;
  beta = 0.0;
}

// Destructor
Pitchfork_FiniteElementProblem::~Pitchfork_FiniteElementProblem()
{
  delete AA;
  delete A;
  delete initialSolution;
  delete Importer;
  delete OverlapMap;
  delete StandardMap;
}

// Matrix and Residual Fills
bool 
Pitchfork_FiniteElementProblem::evaluate(FillType f, 
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
    cout << "ERROR: Pitchfork_FiniteElementProblem::fillMatrix() - FillType flag is broken" << endl;
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
  double jac;
  double xx[2];
  double uu[2];
  Basis basis;

  // Create the nodal coordinates
  double Length=2.0;
  double dx=Length/((double) NumGlobalElements-1);
  for (i=0; i < OverlapNumMyElements; i++) {
    x[i]=-1.0 + dx*((double) OverlapMinMyGID+i);
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
	      *((-1.0/(basis.dx*basis.dx))*basis.duu*
		basis.dphide[i]-source_term(basis.uu)*basis.phi[i]);
	  }
	}
	// Loop over Trial Functions
	if ((flag == MATRIX_ONLY) || (flag == ALL)) {
	  for(j=0;j < 2; j++) {
	    if (StandardMap->MyGID(row)) {
	      column=OverlapMap->GID(ne+j);
	      jac=basis.wt*basis.dx*((-1.0/(basis.dx*basis.dx))*
				     basis.dphide[j]*basis.dphide[i]
				     -source_deriv(basis.uu)*basis.phi[j]*
				     basis.phi[i]);  
	      ierr=A->SumIntoGlobalValues(row, 1, &jac, &column);
	    }
	  }
	}
      }
    }
  } 

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(-1)=beta
  if (MyPID==0) {
    if ((flag == F_ONLY)    || (flag == ALL)) 
      (*rhs)[0]= (*soln)[0] - beta;
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      column=0;
      jac=1.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column=1;
      jac=0.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
    }
  }

  // U(1)=beta
  if ( StandardMap->LID(StandardMap->MaxAllGID()) >= 0 ) {
    int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
    if ((flag == F_ONLY)    || (flag == ALL)) 
       (*rhs)[lastDof] = (*soln)[lastDof] - beta;
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
       int row = StandardMap->MaxAllGID();
       column = row;
       jac=1.0;
       A->ReplaceGlobalValues(row, 1, &jac, &column);
       column=row-1;
       jac=0.0;
       A->ReplaceGlobalValues(row, 1, &jac, &column);
    }
  }
  // Sync up processors to be safe
  Comm->Barrier();
 
  A->FillComplete();

  return true;
}

Epetra_Vector& 
Pitchfork_FiniteElementProblem::getSolution()
{
  return *initialSolution;
}
  
Epetra_CrsMatrix& 
Pitchfork_FiniteElementProblem::getJacobian()
{
  return *A;
}

bool 
Pitchfork_FiniteElementProblem::setParameter(string label, double value)
{
  if (label == "lambda")
    lambda = value;
  else if (label == "alpha")
    alpha = value;
  else if (label == "beta")
    beta = value;
  else if (label == "Homotopy Continuation Parameter") {
    // do nothing for now
  }
  else {
    // do nothing (may be a constraint parameter that we don't know about)
  }
  return true;
}

Epetra_CrsGraph& 
Pitchfork_FiniteElementProblem::generateGraph(Epetra_CrsGraph& AAA)
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
	  AAA.InsertGlobalIndices(row, 1, &column);
	}
      } 	
    }
  }
  AAA.FillComplete();
  return AAA;
}

double
Pitchfork_FiniteElementProblem::source_term(double x) {
  return lambda*x - alpha*x*x + beta*x*x*x;
}

double
Pitchfork_FiniteElementProblem::source_deriv(double x) {
  return lambda - 2.0*alpha*x + 3.0*beta*x*x;
}
