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
#include "Basis.H"

#include "Equation_A.H"

// Constructor - creates the Epetra objects (maps and vectors) 
Equation_A::Equation_A(Epetra_Comm& comm, int numGlobalNodes) :
  GenericEpetraProblem(comm, numGlobalNodes),
  xmin(0.0),
  xmax(1.0),
  dt(1.0e-1)
{

  // Create mesh and solution vectors

  // We first initialize the mesh and then the solution since the latter
  // can depend on the mesh.
  xptr = new Epetra_Vector(*StandardMap);
  double Length= xmax - xmin;
  dx=Length/((double) NumGlobalNodes-1);
  for (int i=0; i < NumMyNodes; i++) {
    (*xptr)[i]=xmin + dx*((double) StandardMap->MinMyGID()+i);
  }

  // Create extra vector needed for this transient problem
  oldSolution = new Epetra_Vector(*StandardMap);

  // Create extra vector needed as data from other problem (B)
  auxSolution = new Epetra_Vector(*StandardMap);

  // Next we create and initialize the solution vector
  initialSolution = new Epetra_Vector(*StandardMap);
  initializeSolution();

  // Allocate the memory for a matrix dynamically (i.e. the graph is dynamic).
  AA = new Epetra_CrsGraph(Copy, *StandardMap, 0);
  generateGraph();

#ifdef DEBUG
  AA->Print(cout);
#endif

  // Create a matrix using the graph just created - this creates a
  // static graph so we can refill the new matirx after TransformToLocal()
  // is called.
  A = new Epetra_CrsMatrix (Copy, *AA);
  A->TransformToLocal();

  // Create the Importer needed for FD coloring
  ColumnToOverlapImporter = new Epetra_Import(A->ColMap(),*OverlapMap);
}

// Destructor
Equation_A::~Equation_A()
{
  delete A; A = 0;
  delete AA; AA = 0;
  delete xptr; xptr = 0;
  delete oldSolution; oldSolution = 0;
  delete auxSolution; auxSolution = 0;
  delete initialSolution; initialSolution = 0;
  delete ColumnToOverlapImporter; ColumnToOverlapImporter = 0;
}

// Reset function
void Equation_A::reset(const Epetra_Vector& x)
{
  *oldSolution = x;
}

// Empty Reset function
void Equation_A::reset()
{
  cout << "WARNING: reset called without passing any update vector !!" 
       << endl;
}

// Set initialSolution to desired initial condition
void Equation_A::initializeSolution()
{
  // Aliases for convenience
  Epetra_Vector& soln = *initialSolution;
  Epetra_Vector& x = *xptr;

  // Here we do a sinusoidal perturbation of the unstable
  // steady state.

  double pi = 4.*atan(1.0);

  for (int i=0; i<x.MyLength(); i++)
    soln[i] = 0.6 + 1.e-1*sin(1.0*pi*x[i]);
  
  *oldSolution = soln;

/*
  // Fill Equation_B solution from file (obtained from Brusselator)
  char file_name[] = "EqB_soln.dat";
  int idum;
  float ddum1, ddum2, value;
  FILE *ifp;
  ifp = fopen(file_name, "r");
  for (int i=0; i<x.MyLength(); i++) {
    fscanf(ifp, "%d %f %f %f", &idum, &ddum1, &ddum2, &value);
    (*auxSolution)[i] = value;
  }
  fclose(ifp);

//  cout << *auxSolution << endl;
*/
} 

// Matrix and Residual Fills
bool Equation_A::evaluate(FillType f, 
				    const Epetra_Vector* soln, 
				    Epetra_Vector* tmp_rhs, 
				    Epetra_RowMatrix* tmp_matrix,
                                    NOX::Epetra::Interface::FillType fillFlag)
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
    cout << "ERROR: Equation_A::fillMatrix() - FillType flag is broken" << endl;
    throw;
  }

  // Create the overlapped solution and position vectors
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector uold(*OverlapMap);
  Epetra_Vector aux(*OverlapMap);
  Epetra_Vector xvec(*OverlapMap);

  // Export Solution to Overlap vector
  // If the vector to be used in the fill is already in the Overlap form,
  // we simply need to map on-processor from column-space indices to
  // OverlapMap indices. Note that the old solution is simply fixed data that
  // needs to be sent to an OverlapMap (ghosted) vector.  The conditional
  // treatment for the current soution vector arises from use of
  // FD coloring in parallel.
  uold.Import(*oldSolution, *Importer, Insert);
  aux.Import(*auxSolution, *Importer, Insert);
  xvec.Import(*xptr, *Importer, Insert);
  if( fillFlag == NOX::Epetra::Interface::FiniteDifferenceF)
    // Overlap vector for solution received from FD coloring, so simply reorder
    // on processor
    u.Export(*soln, *ColumnToOverlapImporter, Insert);
  else // Communication to Overlap vector is needed
    u.Import(*soln, *Importer, Insert);

  // Declare required variables
  int i,j,ierr;
  int OverlapNumMyNodes = OverlapMap->NumMyElements();

  int OverlapMinMyNodeGID;
  if (MyPID==0) OverlapMinMyNodeGID = StandardMap->MinMyGID();
  else OverlapMinMyNodeGID = StandardMap->MinMyGID()-1;

  int row, column;
  double Dcoeff = 0.025;
  double alpha = 0.6;
  double beta = 2.0;
  double jac;
  double xx[2];
  double uu[2]; 
  double uuold[2];
  double aaux[2];
  Basis basis;

  
  // Zero out the objects that will be filled
  if ((flag == MATRIX_ONLY) || (flag == ALL)) i=A->PutScalar(0.0);
  if ((flag == F_ONLY)    || (flag == ALL)) i=rhs->PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes-1; ne++) {
    
    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++) {
      // Get the solution and coordinates at the nodes 
      xx[0]=xvec[ne];
      xx[1]=xvec[ne+1];
      uu[0] = u[ne];
      uu[1] = u[ne+1];
      uuold[0] = uold[ne];
      uuold[1] = uold[ne+1];
      aaux[0] = aux[ne];
      aaux[1] = aux[ne+1];
      // Calculate the basis function at the gauss point
      basis.getBasis(gp, xx, uu, uuold, aaux);

      // Loop over Nodes in Element
      for (i=0; i< 2; i++) {
	row=OverlapMap->GID(ne+i);
	if (StandardMap->MyGID(row)) {
	  if ((flag == F_ONLY)    || (flag == ALL)) {
	    (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))]+=
	      +basis.wt*basis.dx
	      *((basis.uu - basis.uuold)/dt * basis.phi[i] 
              +(1.0/(basis.dx*basis.dx))*Dcoeff*basis.duu*basis.dphide[i]
              + basis.phi[i] * ( -alpha + (beta+1.0)*basis.uu
                - basis.uu*basis.uu*basis.aaux) );
	  }
	}
	// Loop over Trial Functions
	if ((flag == MATRIX_ONLY) || (flag == ALL)) {
	  for(j=0;j < 2; j++) {
	    if (StandardMap->MyGID(row)) {
	      column=OverlapMap->GID(ne+j);
	      jac=basis.wt*basis.dx*(
                      basis.phi[j]/dt*basis.phi[i] 
                      +(1.0/(basis.dx*basis.dx))*Dcoeff*basis.dphide[j]*
                                                        basis.dphide[i]
                      + basis.phi[i] * ( (beta+1.0)*basis.phi[j]
                      - 2.0*basis.uu*basis.phi[j]*basis.aaux) );  
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
      (*rhs)[0]= (*soln)[0] - 0.6;
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      int column=0;
      double jac=1.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
      column=1;
      jac=0.0;
      A->ReplaceGlobalValues(0, 1, &jac, &column);
    }
  }
  // U(1)=1
  if ( StandardMap->LID(StandardMap->MaxAllGID()) >= 0 ) {
    int lastDof = StandardMap->LID(StandardMap->MaxAllGID());
    if ((flag == F_ONLY)    || (flag == ALL))
      (*rhs)[lastDof] = (*soln)[lastDof] - 0.6;
    if ((flag == MATRIX_ONLY) || (flag == ALL)) {
      int row=StandardMap->MaxAllGID();
      int column = row;
      double jac = 1.0;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
      jac=0.0;
      column--;
      A->ReplaceGlobalValues(row, 1, &jac, &column);
    }
  }

  // Sync up processors to be safe
  Comm->Barrier();
 
  A->TransformToLocal();

#ifdef DEBUG
  A->Print(cout);

  if( (flag == F_ONLY) || (flag == ALL) )
    cout << "For residual fill :" << endl << *rhs << endl;

  if( (flag == MATRIX_ONLY) || (flag == ALL) ) {
    cout << "For jacobian fill :" << endl;
    A->Print(cout);
  }

#endif


  return true;
}

Epetra_Vector& Equation_A::getOldSoln()
{
  return *oldSolution;
} 
  
double Equation_A::getdt()
{
  return dt;
}

void Equation_A::generateGraph()
{
  
  // Declare required variables
  int i,j;
  int row, column;
  int OverlapNumMyNodes = OverlapMap->NumMyElements();
  int OverlapMinMyNodeGID;
  if (MyPID==0) OverlapMinMyNodeGID = StandardMap->MinMyGID();
  else OverlapMinMyNodeGID = StandardMap->MinMyGID()-1;
  
  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyNodes-1; ne++) {
          
    // Loop over Nodes in Element
    for (i=0; i<2; i++) {

      // If this node is owned by current processor, add indices
      if (StandardMap->MyGID(OverlapMap->GID(ne+i))) {

        // Loop over unknowns in Node
        row=OverlapMap->GID(ne+i);

        // Loop over supporting nodes
        for(j=0; j<2; j++) {

          // Loop over unknowns at supporting nodes
          column=OverlapMap->GID(ne+j);
          //printf("\t\tWould like to insert -> (%d, %d)\n",row,column);
          AA->InsertGlobalIndices(row, 1, &column);
        }
      }
    }
  }
  AA->TransformToLocal();
//   AA->SortIndices();
//   AA->RemoveRedundantIndices();
  
  return;
}
