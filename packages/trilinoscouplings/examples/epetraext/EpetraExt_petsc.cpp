// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ml_config.h"
#include "EpetraExt_config.h"
#if defined(HAVE_PETSC) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO) && defined(HAVE_MPI)
#include "petscksp.h"
#include "EpetraExt_PETScAIJMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"

/* 
   This example demonstrates how to apply a Trilinos preconditioner to a PETSc
   linear system.

   For information on configuring and building Trilinos with the PETSc aij
   interface enabled, please see EpetraExt's doxygen documentation at
   http://trilinos.sandia.gov/packages/epetraext, development version
   or release 9.0 or later.

   The PETSc matrix is a 2D 5-point Laplace operator stored in AIJ format.
   This matrix is wrapped as an Epetra_PETScAIJMatrix, and an ML AMG
   preconditioner is created for it.  The associated linear system is solved twice,
   the first time using AztecOO's preconditioned CG, the second time using PETSc's.

   To invoke this example, use something like:

       mpirun -np 5 ./TrilinosCouplings_petsc.exe -m 150 -n 150 -petsc_smoother -ksp_monitor_true_residual
*/

static char help[] = "Demonstrates how to solve a PETSc linear system with KSP\
and a Trilinos AMG preconditioner.  In particular, it shows how to wrap a PETSc\
AIJ matrix as an Epetra matrix, create the AMG preconditioner, and wrap it as a\
shell preconditioner for a PETSc Krylov method.\
Input parameters include:\n\
  -random_exact_sol : use a random exact solution vector\n\
  -view_exact_sol   : write exact solution vector to stdout\n\
  -m <mesh_x>       : number of mesh points in x-direction\n\
  -n <mesh_n>       : number of mesh points in y-direction\n\n";

#if (PETSC_VERSION_MAJOR < 3) || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 1)
extern PetscErrorCode ShellApplyML(void*,Vec,Vec);
#else
extern PetscErrorCode ShellApplyML(PC,Vec,Vec);
#endif

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            x,b,u;  /* approx solution, RHS, exact solution */
  Mat            A;        /* linear system matrix */
  KSP            ksp;     /* linear solver context */
  KSP            kspSmoother=0;  /* solver context for PETSc fine grid smoother */
  PC pc;
  PetscRandom    rctx;     /* random number generator context */
  PetscReal      norm;     /* norm of solution error */
  PetscInt       i,j,Ii,J,Istart,Iend,its;
  PetscInt       m = 50,n = 50; /* #mesh points in x & y directions, resp. */
  PetscErrorCode ierr;
  PetscBool     flg;
  PetscScalar    v,one = 1.0,neg_one = -1.0;
  PetscInt rank=0;
  MPI_Comm comm;

  PetscInitialize(&argc,&args,(char *)0,help);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n);CHKERRQ(ierr);
  ierr = MatSetType(A, MATAIJ);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,5,PETSC_NULL,5,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  PetscObjectGetComm( (PetscObject)A, &comm);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (!rank) printf("Matrix has %d (%dx%d) rows\n",m*n,m,n);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

  for (Ii=Istart; Ii<Iend; Ii++) { 
    v = -1.0; i = Ii/n; j = Ii - i*n;  
    if (i>0)   {J = Ii - n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
    if (i<m-1) {J = Ii + n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
    if (j>0)   {J = Ii - 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
    if (j<n-1) {J = Ii + 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
    v = 4.0; ierr = MatSetValues(A,1,&Ii,1,&Ii,&v,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Wrap the PETSc matrix as an Epetra_PETScAIJMatrix. This is lightweight,
     i.e., no deep data copies. */
  Epetra_PETScAIJMatrix epA(A);

  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = VecSetSizes(u,PETSC_DECIDE,m*n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&b);CHKERRQ(ierr); 

  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL,"-random_exact_sol",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
    ierr = VecSetRandom(u,rctx);CHKERRQ(ierr);
    ierr = PetscRandomDestroy(&rctx);CHKERRQ(ierr);
  } else {
    ierr = VecSet(u,one);CHKERRQ(ierr);
  }
  ierr = MatMult(A,u,b);CHKERRQ(ierr);
  ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
  if (rank==0) printf("||b|| = %f\n",norm);

  /* Copy the PETSc vectors to Epetra vectors. */
  Epetra_Vector epu(epA.Map()), epb(epA.Map());
  PetscScalar *vals;
  ierr = VecGetArray(u,&vals);CHKERRQ(ierr);
  PetscInt length;
  ierr = VecGetLocalSize(u,&length);CHKERRQ(ierr);
  for (int j=0; j<length; j++) epu[j] = vals[j];
  ierr = VecRestoreArray(u,&vals);CHKERRQ(ierr);
  epA.Multiply(false, epu, epb);

  /* Check norms of the Epetra and PETSc vectors. */
  epu.Norm2(&norm);
  if (rank == 0) printf("||epetra u||_2 = %f\n",norm);
  ierr = VecNorm(u,NORM_2,&norm);CHKERRQ(ierr);
  if (rank == 0) printf("||petsc u||_2  = %f\n",norm);
  epb.Norm2(&norm);
  if (rank == 0) printf("||epetra b||_2 = %f\n",norm);
  ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
  if (rank == 0) printf("||petsc b||_2  = %f\n",norm);

  /* Create the ML AMG preconditioner. */

  /* Parameter list that holds options for AMG preconditioner. */
  Teuchos::ParameterList mlList;
  /* Set recommended defaults for Poisson-like problems. */
  ML_Epetra::SetDefaults("SA",mlList);
  /* Specify how much information ML prints to screen.
     0 is the minimum (no output), 10 is the maximum. */
  mlList.set("ML output",10);
  /* Set the fine grid smoother.  PETSc will be much faster for any
     smoother requiring row access, e.g., SOR.  For any smoother whose
     kernel is a matvec, Trilinos/PETSc performance should be comparable,
     as Trilinos simply calls the PETSc matvec.

     To use a PETSc smoother, create a KSP object, set the KSP type to
     KSPRICHARDSON, and set the desired smoother as the KSP preconditioner.
     It is important that you call KSPSetInitialGuessNonzero.  Otherwise, the
     post-smoother phase will incorrectly ignore the current approximate
     solution.  The KSP pointer must be cast to void* and passed to ML via
     the parameter list.

     You are responsible for freeing the KSP object.
  */
  ierr = PetscOptionsHasName(PETSC_NULL,"-petsc_smoother",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = KSPCreate(comm,&kspSmoother);CHKERRQ(ierr);
    ierr = KSPSetOperators(kspSmoother,A,A);CHKERRQ(ierr);
    ierr = KSPSetType(kspSmoother,KSPRICHARDSON);CHKERRQ(ierr);
    ierr = KSPSetTolerances(kspSmoother, 1e-12, 1e-50, 1e7,1);
    ierr = KSPSetInitialGuessNonzero(kspSmoother,PETSC_TRUE);CHKERRQ(ierr);
    ierr = KSPGetPC(kspSmoother,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc, PCSOR);CHKERRQ(ierr);
    ierr = PCSORSetSymmetric(pc,SOR_LOCAL_SYMMETRIC_SWEEP);CHKERRQ(ierr);
    ierr = PCSetFromOptions(pc);CHKERRQ(ierr);
    ierr = KSPSetUp(kspSmoother);CHKERRQ(ierr);
    mlList.set("smoother: type (level 0)","petsc");
    mlList.set("smoother: petsc ksp (level 0)",(void*)kspSmoother);
  } else {
    mlList.set("smoother: type (level 0)","symmetric Gauss-Seidel");
  }

  /* how many fine grid pre- or post-smoothing sweeps to do */
  mlList.set("smoother: sweeps (level 0)",2);

  ML_Epetra::MultiLevelPreconditioner *Prec = new ML_Epetra::MultiLevelPreconditioner(epA,mlList);

  /* Trilinos CG */
  epu.PutScalar(0.0);
  Epetra_LinearProblem Problem(&epA, &epu, &epb);
  AztecOO solver(Problem);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);
  solver.SetAztecOption(AZ_conv, AZ_noscaled);
  solver.SetPrecOperator(Prec);
  solver.Iterate(30, 1e-12);

  /* PETSc CG */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1e-12,1.e-50,PETSC_DEFAULT,
                          PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);

  /* Wrap the ML preconditioner as a PETSc shell preconditioner. */
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,ShellApplyML);CHKERRQ(ierr);
  ierr = PCShellSetContext(pc,(void*)Prec);CHKERRQ(ierr);
  ierr = PCShellSetName(pc,"ML AMG");CHKERRQ(ierr);

  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %e iterations %D\n",
                     norm,its);CHKERRQ(ierr);

  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);  ierr = MatDestroy(&A);CHKERRQ(ierr);

  if (kspSmoother) {ierr = KSPDestroy(&kspSmoother);CHKERRQ(ierr);}
  if (Prec) delete Prec;

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
} /*main*/

/* ***************************************************************** */

#if (PETSC_VERSION_MAJOR < 3) || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 1)
PetscErrorCode ShellApplyML(void *ctx,Vec x,Vec y)
#else
PetscErrorCode ShellApplyML(PC pc,Vec x,Vec y)
{
  PetscErrorCode  ierr;
  ML_Epetra::MultiLevelPreconditioner *mlp = 0;
  void* ctx;

  ierr = PCShellGetContext(pc,&ctx); CHKERRQ(ierr);  
  mlp = (ML_Epetra::MultiLevelPreconditioner*)ctx;
#endif

  /* Wrap x and y as Epetra_Vectors. */
  // NOTE: We should not be const casting xvals, but we do promise not to modify the entries.
  const PetscScalar *xvals;
  PetscScalar *nonconst_xvals, *yvals;
  ierr = VecGetArrayRead(x,&xvals);CHKERRQ(ierr);
  nonconst_xvals = const_cast<PetscScalar*>(xvals);
  Epetra_Vector epx(View,mlp->OperatorDomainMap(),nonconst_xvals);
  ierr = VecGetArray(y,&yvals);CHKERRQ(ierr);
  Epetra_Vector epy(View,mlp->OperatorRangeMap(),yvals);

  /* Apply ML. */
  mlp->ApplyInverse(epx,epy);
  
  /* Clean up and return. */
  ierr = VecRestoreArrayRead(x,&xvals);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&yvals);CHKERRQ(ierr);
  return 0;
} /*ShellApplyML*/
#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  printf("Please configure EpetraExt with:");
#if !defined(HAVE_PETSC)
  printf("--enable-petsc");
#endif
#if !defined(HAVE_ML_EPETRA)
  printf("--enable-epetra");
#endif
#if !defined(HAVE_ML_TEUCHOS)
  printf("--enable-teuchos");
#endif
#if !defined(HAVE_ML_AZTECOO)
  printf("--enable-aztecoo");
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return(EXIT_SUCCESS);
}
#endif
