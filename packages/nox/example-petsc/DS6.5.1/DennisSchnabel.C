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
#include "petscsles.h"

#include "DennisSchnabel.H"

// Constructor - creates the Petsc objects (maps and vectors) 
DennisSchnabel::DennisSchnabel(int numGlobalElements) :
  NumGlobalElements(numGlobalElements)
{

  // Commonly used variables
  int i, ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);//CHKERRQ(ierr);// Process ID
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&NumProc);//CHKERRQ(ierr);// # Procs 

  // Construct a Source Map that puts approximately the same 
  // Number of equations on each processor in uniform global ordering
  //StandardMap = new Epetra_Map(NumGlobalElements, 0, *Comm);

  // Get the number of elements owned by this processor
  //NumMyElements = StandardMap->NumMyElements();

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
    //OverlapMap = new Epetra_Map(*StandardMap);
  } 
  else {

    //int OverlapNumMyElements = 2;
    //int OverlapMyGlobalElements[2];
    
    //for (i = 0; i < OverlapNumMyElements; i ++) 
    //  OverlapMyGlobalElements[i] = i;
    
    //OverlapMap = new Epetra_Map(-1, OverlapNumMyElements, 
    //				OverlapMyGlobalElements, 0, *Comm);
  } // End Overlap map construction *************************************

  // Construct Linear Objects  
  initialSolution = new Vec;
  ierr = VecCreate(PETSC_COMM_WORLD, initialSolution);//CHKERRQ(ierr);
  ierr = VecSetSizes(*initialSolution, PETSC_DECIDE, 2);//CHKERRQ(ierr);
  ierr = VecSetFromOptions(*initialSolution);//CHKERRQ(ierr);

  A = new Mat;
  MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2, 2, A);//CHKERRQ(ierr);
  ierr = MatSetFromOptions(*A);//CHKERRQ(ierr);

  // Create Mapping for overlap solution vector using Petsc IS
  overlapSolution = new Vec;
  VecCreateSeq(PETSC_COMM_SELF,2,overlapSolution);
  int indexMap[] = {0, 1};
  IS from, to;
  ISCreateGeneral(PETSC_COMM_WORLD,2,indexMap,&from);
  ISCreateGeneral(PETSC_COMM_WORLD,2,indexMap,&to);
  petscMap = new VecScatter;
  VecScatterCreate(*initialSolution, from, *overlapSolution, to, petscMap);
  
}

// Destructor
DennisSchnabel::~DennisSchnabel()
{
  delete A;
  delete initialSolution;
}

// Matrix and Residual Fills
bool DennisSchnabel::evaluate(FillType f, 
			      const Vec* soln, 
			      Vec* tmp_rhs, 
			      Mat* tmp_matrix)
{
  flag = f;

  // Set the incoming linear objects
  if (flag == RHS_ONLY) {
    rhs = tmp_rhs;
  } 
  else if (flag == MATRIX_ONLY) {
    A = tmp_matrix;
  } 
  else if (flag == ALL) { 
    rhs = tmp_rhs;
    A = tmp_matrix;
  } 
  else {
    cout << "ERROR: DennisSchnabel::fillMatrix() - No such flag as " 
	 << flag << endl;
    throw;
  }

  // Create the overlapped solution
  VecScatterBegin(*soln, *overlapSolution, INSERT_VALUES, SCATTER_FORWARD, *petscMap);
  VecScatterEnd(*soln, *overlapSolution, INSERT_VALUES, SCATTER_FORWARD, *petscMap);
  // Export Solution to Overlap vector so we have all unknowns required
  // for function and Jacobian evaluations.
  double* u;
  VecGetArray(*overlapSolution, &u);
//  for(int i=0; i<2; i++)
//    cout << i << " :\t" << u[i] << endl;
//  cin.get();


  // Declare required variables
  int i,ierr;
  double zero = 0.0;
  double resid[2];
  int* column = new int[2];
  column[0] = 0; 
  column[1] = 1;

  // Begin RHS fill
  if((flag == RHS_ONLY) || (flag == ALL)) {

    // Zero out the RHS vector
    //ierr = VecSet(&zero, *rhs);CHKERRQ(ierr);
    //ierr = VecAssemblyBegin(*rhs);CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(*rhs);CHKERRQ(ierr);

    // Processor 0 always fills the first equation.
    if(NumProc==1) {
      resid[0] = u[0]*u[0] + u[1]*u[1] - 2.;
      resid[1] = exp(u[0]-1.) + u[1]*u[1]*u[1] - 2.;
      ierr = VecSetValues(*rhs,2,column,resid,INSERT_VALUES);CHKERRQ(ierr);
    }
    else {  // NumProc==2
      double value;
      if (MyPID==0) { 
        value = u[0]*u[0] + u[1]*u[1] - 2.;
      }
      else {
        value = exp(u[0]-1.) + u[1]*u[1]*u[1] - 2.;
      }
      ierr = VecSetValues(*rhs,1,&column[MyPID],&value,
                                              INSERT_VALUES);CHKERRQ(ierr);
    }
    /*
       Assemble vector, using the 2-step process:
       VecAssemblyBegin(), VecAssemblyEnd()
       Computations can be done while messages are in transition,
       by placing code between these two statements.
    */
    ierr = VecAssemblyBegin(*rhs);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*rhs);CHKERRQ(ierr);

// Check resid values from Vec rhs
//  Vec mask;
//  double dmask[] = {1.0, 0.0};
//  ierr = VecCreate(PETSC_COMM_WORLD, &mask);CHKERRQ(ierr);
//  ierr = VecSetSizes(mask, PETSC_DECIDE, 2);CHKERRQ(ierr);
//  ierr = VecSetFromOptions(mask);CHKERRQ(ierr);
//  ierr = VecSetValues(mask,2,column,dmask, INSERT_VALUES);CHKERRQ(ierr);
//  ierr = VecAssemblyBegin(*rhs);CHKERRQ(ierr);
//  ierr = VecAssemblyEnd(*rhs);CHKERRQ(ierr);
//  
//  cout << "Checking residuals .... " << endl;
//  double val1, val2;
//  VecDot(*rhs, mask, &val1);
//  dmask[0] = 0.0;
//  dmask[1] = 1.0;
//  ierr = VecSetValues(mask,2,column,dmask, INSERT_VALUES);CHKERRQ(ierr);
//  ierr = VecAssemblyBegin(*rhs);CHKERRQ(ierr);
//  ierr = VecAssemblyEnd(*rhs);CHKERRQ(ierr);
//  VecDot(*rhs, mask, &val2);
//  cout << "rhs vals:\t" << val1 << "\t" << val2 << endl;
//  cin.get();

  }

  // Here we fill a Jacobian double array on all active procs and
  // insert values  
  double* jac = new double[4];

  // The matrix is 2 x 2 and will always be 0 and 1 regardless of 
  // the coordinates being local or global.

  // Begin Jacobian fill
  if((flag == MATRIX_ONLY) || (flag == ALL)) {

    // Zero out Jacobian
    //i=A->PutScalar(0.0);

    jac[0] = 2.*u[0];
    jac[1] = 2.*u[1];
    jac[2] = exp(u[0]-1.);
    jac[3] = 3.*u[1]*u[1];
    ierr = MatSetValues(*A,2,column,2,column,jac,INSERT_VALUES);CHKERRQ(ierr);
    /* 
       Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd()
       Computations can be done while messages are in transition
       by placing code between these two statements.
    */
    ierr = MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    delete [] column;
    delete [] jac;
  } 

  VecRestoreArray(*overlapSolution, &u);

  return true;
}

Vec& DennisSchnabel::getSolution()
{
  return *initialSolution;
}
  
Mat& DennisSchnabel::getJacobian()
{
  return *A;
}

