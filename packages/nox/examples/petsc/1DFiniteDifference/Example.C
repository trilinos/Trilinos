// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//  Problem description goes here.

static char help[] =
       "1D Finite Difference example problem in parallel.\n\n";


// Petsc Objects
/*
   Include "petscda.h" so that we can use distributed arrays (DAs).
   Include "petscdraw.h" so that we can use PETSc drawing routines.
   Include "petscsnes.h" so that we can use SNES solvers.  Note that this
   file automatically includes:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - system routines       petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
     petscsles.h   - linear solvers
*/

#include "petscdm.h"
#include "petscdmda.h"
#include "petscsnes.h"

// NOX Library
#include "NOX.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Petsc_Vector.H"
#include "NOX_Petsc_SharedJacobian.H"
#include "NOX_Petsc_Group.H"
#include "NOX_Petsc_Options.H"

// User's application specific files
#include "Problem_Interface.H" // Interface to NOX
#include "FiniteDifference.H"  // The PDE class used for fills

/*
   User-defined routines.  Note that immediately before each routine below,
   we define the macro __FUNCT__ to be a std::string containing the routine name.
   If defined, this macro is used in the PETSc error handlers to provide a
   complete traceback of routine names.  All PETSc library routines use this
   macro, and users can optionally employ it as well in their application
   codes.  Note that users can get a traceback of PETSc errors regardless of
   whether they define __FUNCT__ in application codes; this macro merely
   provides the added traceback detail of the application routine names.
*/
int FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
int FormFunction(SNES,Vec,Vec,void*);
int FormInitialGuess(Vec);
int Monitor(SNES,int,PetscReal,void *);
int StepCheck(SNES,void *,Vec,PetscBool *);

/*
   User-defined context for monitoring
*/
typedef struct {
   PetscViewer viewer;
} MonitorCtx;

/*
   User-defined context for checking candidate iterates that are
   determined by line search methods
*/
typedef struct {
   Vec       last_step;  /* previous iterate */
   PetscReal tolerance;  /* tolerance for changes between successive iterates */} StepCheckCtx;

using namespace std;

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char *argv[])
{
  int ierr = 0; // Error inidicator used with CHKERRQ macro

  ApplicationCtx ctx;                  /* user-defined context */

  // Initialize MPI
  int N = 5;
  PetscInitialize(&argc,&argv,(char *)0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&ctx.rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&ctx.size);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&N,PETSC_NULL);CHKERRQ(ierr);
  ctx.h = 1.0/(N-1);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create a dummy nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  SNES snes; // Needed as a dummy arg to FromResidual, etc.
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create vector data structures; set function evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create distributed array (DA) to manage parallel grid and vectors
  */
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,N,1,1,PETSC_NULL,&ctx.da);CHKERRQ(ierr);

  /*
     Extract global and local vectors from DA; then duplicate for remaining
     vectors that are the same types
  */
  Vec x, r;
  ierr = DMCreateGlobalVector(ctx.da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

  // Create the Problem class.  For now, this is simply a holder for
  // objects needed to interface to native FormXXX routines.
  FiniteDifference Problem(&snes,&ctx);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix data structures; set Jacobian evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  Mat J;
  ierr = MatCreate( PETSC_COMM_SELF, &J );
  ierr = MatSetSizes( J, PETSC_DECIDE, PETSC_DECIDE, N, N);
  ierr = MatSetFromOptions(J);CHKERRQ(ierr);
  ierr = MatSetUp(J);

  // This should be replaced by NOX interface setup
  //ierr = SNESSetJacobian(snes,J,J,FormJacobian,&ctx);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize application:
     Store right-hand-side of PDE and exact solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Get local grid boundaries (for 1-dimensional DA):
       xs, xm - starting grid index, width of local grid (no ghost points)
  */
  int xs, xm;
  ierr = DMDAGetCorners(ctx.da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess; then solve nonlinear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscReal one = 1.0;
  ierr = VecSet( x, one );CHKERRQ(ierr);

// Additional NOX setup

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());


  // Allow options to be set from command line or from file
  NOX::Petsc::Options optionHandler(nlParams, ctx.rank);

  // Create the interface between the test problem and the nonlinear solver
  // This is created using inheritance of the abstract base class:
  // NOX::Petsc::Interface
  Problem_Interface interface(Problem);

  // Create the Group
  //NOX::Petsc::Group* grp = new NOX::Petsc::Group(interface, x, J);
  Teuchos::RCP<NOX::Petsc::Group> grp = Teuchos::rcp( new NOX::Petsc::Group(interface, x, J) );
  grp->computeF(); // Needed to establish the initial convergence state

  // Create the method and solve
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(grp, optionHandler.getStatusTest(), nlParamsPtr);
  NOX::StatusTest::StatusType status = solver->solve();

  if (status != NOX::StatusTest::Converged)
    std::cout << "Nonlinear solver failed to converge!" << std::endl;

  // Get the Petsc_Vector with the final solution from the solver
  const NOX::Petsc::Group& finalGroup =
      dynamic_cast<const NOX::Petsc::Group&>(solver->getSolutionGroup());
  const Vec& finalSolution = (dynamic_cast<const NOX::Petsc::Vector&>
        (finalGroup.getX())).getPetscVector();

  // End Nonlinear Solver **************************************

  printf("\n\tHere is the solution:\n\n");
  VecView(finalSolution, PETSC_VIEWER_STDOUT_WORLD);

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx.da);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}

/* -----------------------------------------------------
 // end main
--------------------------------------------------*/
