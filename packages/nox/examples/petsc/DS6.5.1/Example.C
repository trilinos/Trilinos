// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//  Simple 2 equation test for quadratic and cubic line searches
//  from Dennis & Schnabel's book, chp 6.  The test problem is from
//  Example 6.5.1
/*
 *    U0**2 + U1**2 - 2 = 0
 *    exp(U0-1) + U1**3 -2 = 0
 *
 *  NOTE: To reproduce the results from the reference, the linesearch option
 *        must be set to polynomial.  This can be done either from the
 *        command line using "-nox_linesearch_type polynomial" or by
 *        placing the same option specification in an input file,
 *        eg ${HOME}/.petscrc
 *        Also, the solver must be -nox_linesearch_based with
 *        "-nox_direction_type newton", which can be used by not
 *        explicitly setting either of these since they are the defaults.
 */

static char help[] =
       "Solves Dennis & Schnabel example problem in parallel.\n\n";


// Petsc Objects
#include "petscksp.h"
/*
  Include "petscsles.h" so that we can use SLES solvers.  Note that this file
  automatically includes:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - system routines       petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/

// NOX Library
#include "NOX.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Petsc_Vector.H"
#include "NOX_Petsc_SharedJacobian.H"
#include "NOX_Petsc_Group.H"
#include "NOX_Petsc_Options.H"

// User's application specific files
#include "Problem_Interface.H" // Interface file to NOX
#include "DennisSchnabel.H"

using namespace std;

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char *argv[])
{
  int ierr = 0;

  // Initialize MPI
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);

  // Get the process ID and the total number of processors
  int MyPID, NumProc;
#ifdef HAVE_MPI
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&NumProc);CHKERRQ(ierr);
#else
  MyPID = 0;
  NumProc = 1;
#endif

  int NumGlobalElements = 2;  // Hardcoded for D&S Example problem

  // A maximum of 2 procesors is allowed since there are only 2 equations
  if (NumProc >= 3) {
    std::cout << "ERROR: Maximum number of processors is 2!" << std::endl;
    exit(1);
  }

  // Create the Problem class.  This creates all required
  // Petsc objects for the problem and allows calls to the
  // function (RHS) and Jacobian evaluation routines.
  DennisSchnabel Problem(NumGlobalElements);

  // Get the vector from the Problem
  Vec& soln = Problem.getSolution();

  // Initialize Solution. For simplicity, this is done on all (both) procs.
  int globalIndex[2];
  double doubleArray[2];
  globalIndex[0] = 0;
  globalIndex[1] = 1;
  doubleArray[0] = 2.0;
  doubleArray[1] = 0.5;
  ierr = VecSetValues(soln,2,globalIndex,doubleArray,
                                         INSERT_VALUES);CHKERRQ(ierr);

  ierr = VecAssemblyBegin(soln);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(soln);CHKERRQ(ierr);

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  NOX::Petsc::Options optionHandler(nlParams, MyPID);

  // Warn if not using Polynimial linesearch
  std::string lsMethod = nlParams.sublist("Line Search").get("Method", "Full Step");
  if( "Polynomial" != lsMethod )
    if (MyPID==0)
      std::cout << "Line Search Method is set to \"" << lsMethod << "\".  This test is designed to "
           << "work with \"Polynomial\".  You can set this using \"-nox_linesearch_type polynomial\" "
           << "on the command line." << std::endl;

  // Create the interface between the test problem and the nonlinear solver
  // This is created using inheritance of the abstract base class:
  // NOX::Petsc::Interface
  Problem_Interface interface(Problem);

  // Get a reference to the Petsc_RowMatrix created in Problem_Interface.
  Mat& A = Problem.getJacobian();

  // Create the Group
  Teuchos::RCP<NOX::Petsc::Group> grp = Teuchos::rcp( new NOX::Petsc::Group(interface, soln, A) );
  grp->computeF(); // Needed to establish the initial convergence state

  // Create the method and solve
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(grp, optionHandler.getStatusTest(), nlParamsPtr);
  NOX::StatusTest::StatusType status = solver->solve();

  if (status != NOX::StatusTest::Converged)
  {
    if (MyPID==0)
      std::cout << "Nonlinear solver failed to converge!" << std::endl;
  if( "Polynomial" != lsMethod )
      if (MyPID==0)
        std::cout << "\nLine Search Method is set to \"" << lsMethod << "\".  This test is designed to "
             << "work with \"Polynomial\".  You can set this using \"-nox_linesearch_type polynomial\" "
             << "on the command line." << std::endl;
  }


  // Get the Petsc_Vector with the final solution from the solver
  const NOX::Petsc::Group& finalGroup =
      dynamic_cast<const NOX::Petsc::Group&>(solver->getSolutionGroup());
  const Vec& finalSolution = (dynamic_cast<const NOX::Petsc::Vector&>
        (finalGroup.getX())).getPetscVector();

  // End Nonlinear Solver **************************************

  // Print solution
  if(MyPID==0)
    printf("Final solution :\n");
  char file_name[25];
  FILE *ifp;
  (void) sprintf(file_name, "output.%d",MyPID);
  double* finalVals;
  ifp = fopen(file_name, "w");
  VecGetArray( finalSolution, &finalVals );
  for(int i=0; i<3-NumProc; i++) {
    printf("(proc %d)\t%d\t%e\n",MyPID,i,finalVals[i]);
    fprintf(ifp, "%d  %E\n",i, finalVals[i]);
  }
  VecRestoreArray( finalSolution, &finalVals );
  fclose(ifp);

  ierr = PetscFinalize();CHKERRQ(ierr);

return 0;

} // end main
