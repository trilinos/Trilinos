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

#include "petscda.h"
#include "petscsnes.h"

// NOX Library
#include "NOX.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Petsc_Vector.H"
#include "NOX_Petsc_SharedJacobian.H"
#include "NOX_Petsc_Group.H"

// User's application specific files 
#include "Problem_Interface.H" // Interface file to NOX
#include "FiniteDifference.H"              

/*
   User-defined routines.  Note that immediately before each routine below,
   we define the macro __FUNCT__ to be a string containing the routine name.
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
int StepCheck(SNES,void *,Vec,PetscTruth *);

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
  ierr = DACreate1d(PETSC_COMM_WORLD,DA_NONPERIODIC,N,1,1,PETSC_NULL,&ctx.da);CHKERRQ(ierr);

  /*
     Extract global and local vectors from DA; then duplicate for remaining
     vectors that are the same types
  */
  Vec x, r;
  ierr = DACreateGlobalVector(ctx.da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

  // Create the Problem class.  For now, this is simply a holder for
  // objects needed to interface to native FormXXX routines.
  FiniteDifference Problem(&snes,&ctx);

  // This should be replaced by NOX interface setup
  //ierr = SNESSetFunction(snes,r,FormFunction,&ctx);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix data structures; set Jacobian evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  Mat J;
  ierr = MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,N,N,&J);CHKERRQ(ierr);
  ierr = MatSetFromOptions(J);CHKERRQ(ierr);

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
  ierr = DAGetCorners(ctx.da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess; then solve nonlinear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscReal zero = 0.0, pfive = 0.5, one = 1.0;
  ierr = VecSet(&one,x);CHKERRQ(ierr);

// Additional NOX setup

  // Create the top level parameter list
  NOX::Parameter::List nlParams;

  // Specify nonlinear solver method
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");
  //nlParams.setParameter("Nonlinear Solver", "Trust Region Based");

  // Set the printing parameters in the "Printing" sublist
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", ctx.rank);
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
  searchParams.setParameter("Method", "Full Step");
  //searchParams.setParameter("Method", "Interval Halving");
  //searchParams.setParameter("Method", "Polynomial");
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
  // OR
  //dirParams.setParameter("Method", "Steepest Descent");
  //NOX::Parameter::List& sdParams = dirParams.sublist("Steepest Descent");
    //sdParams.setParameter("Scaling Type", "None");
    //sdParams.setParameter("Scaling Type", "2-Norm");
    //sdParams.setParameter("Scaling Type", "Quadratic Model Min");
  // OR
  //dirParams.setParameter("Method", "NonlinearCG");
  //NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
    //nlcgParams.setParameter("Restart Frequency", 2000);
    //nlcgParams.setParameter("Precondition", "On");
    //nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
    //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

  // Sublist for linear solver
  // Note that preconditioning options as well as the following can be
  // specified on the command line or via the file .petscrc
  // See Petsc documentation for more info.
  NOX::Parameter::List& lsParams = dirParams.sublist("Linear Solver");
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Iteration Output Frequency", 50);    
  lsParams.setParameter("Preconditioning Matrix Type", "None"); 

  // Create the interface between the test problem and the nonlinear solver
  // This is created using inheritance of the abstract base class:
  // NOX::Petsc::Interface
  Problem_Interface interface(Problem);

  // Create the Group
  NOX::Petsc::Group* grp = new NOX::Petsc::Group(lsParams, interface, x, J);
//  NOX::Petsc::Group* grp2 = new NOX::Petsc::Group(*grp);
//  cout << "\n\n\t\tGroup created .....";
//  grp->computeF();
//  grp2->computeF();
//  delete grp;
//  delete grp2;
//  cout << " and destroyed !!" << endl;
/*
  A test to demonstrate memory leaks

  for(int j=0; j<1000000; j++)
  {
    grp = new NOX::Petsc::Group(lsParams, interface, x, J);
    delete grp;
    if( j % 10000 == 0 )
      cout << "  Finished : " << j << endl;
  }
  exit(0);
*/

  grp->computeF(); // Needed to establish the initial convergence state

  // Create the convergence tests
  NOX::StatusTest::NormF testNormF(1.0e-6);
  NOX::StatusTest::MaxIters testMaxIters(10);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, testNormF, testMaxIters);

  //ierr = SNESSolve(snes,x,&it);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF,"Newton iterations = %d\n\n",it);CHKERRQ(ierr);

  // Create the method and solve
  NOX::Solver::Manager solver(*grp, combo, nlParams);
// 
  //Cut the solve step for now to debug....

  NOX::StatusTest::StatusType status = solver.solve();

  if (status != NOX::StatusTest::Converged)
    cout << "Nonlinear solver failed to converge!" << endl;

  // Get the Petsc_Vector with the final solution from the solver
  const NOX::Petsc::Group& finalGroup = 
      dynamic_cast<const NOX::Petsc::Group&>(solver.getSolutionGroup());
  const Vec& finalSolution = (dynamic_cast<const NOX::Petsc::Vector&>
        (finalGroup.getX())).getPetscVector();

  // End Nonlinear Solver **************************************

  Vec& nonconst_x = const_cast<Vec&>(x);

  printf("\n\tHere is the solution:\n\n");
//  VecView(nonconst_x, PETSC_VIEWER_STDOUT_WORLD);
  VecView(finalSolution, PETSC_VIEWER_STDOUT_WORLD);
//

//  NOX::Petsc::Group* grp2 = new NOX::Petsc::Group(*grp);
//  delete grp; grp = 0;
//  delete grp2;
//  delete grp; grp = 0;

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);
  ierr = MatDestroy(J);CHKERRQ(ierr);
  ierr = SNESDestroy(snes);CHKERRQ(ierr);
  ierr = DADestroy(ctx.da);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}

/* -----------------------------------------------------
 // end main
--------------------------------------------------*/
