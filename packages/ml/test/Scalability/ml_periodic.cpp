/////////////////////////////////////////////////////////////
// The goal of this code is to test the multiple right hand sides
// capabilities of ML
/////////////////////////////////////////////////////////////
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA)
#include "Epetra_ConfigDefs.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Time.h"

#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_operator.h"
#include "ml_epetra_utils.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "AztecOO.h"
#include <fstream>
using namespace Teuchos;

int main(int argc, char *argv[]) {

  int myPid = 0;
  int numProc = 1;

#ifdef ML_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPid);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  int numElePerDirection = 101;

  if (argc > 1) sscanf(argv[1],"%d",&numElePerDirection);
   
  int numNodes = (numElePerDirection - 1)*(numElePerDirection - 1);

  // Determine the number of processors to use in each coordinate 
  // direction. The general idea is to try and use about the same
  // number of processors in each direction.

  int numProcXDirection = (int)ceil(sqrt((double)(numProc)));
  int numProcYDirection = numProcXDirection;
  while ( numProcYDirection*numProcXDirection != numProc ) {
     numProcXDirection--; 
     numProcYDirection = numProc/numProcXDirection; 
  }
  int nx = numElePerDirection-1; 
  int ny = nx;

  // Figure out how many local points in this processor and determine the
  // corresponding global indices.

  int PerProcSmallXDirection= (int) (((double) nx)/((double) numProcXDirection));
  int PerProcSmallYDirection= (int) (((double) ny)/((double) numProcYDirection));
  int NBigXDirection = nx - PerProcSmallXDirection*numProcXDirection;
  int NBigYDirection = ny - PerProcSmallYDirection*numProcYDirection;

  int startx, starty, endx, endy;
  int xpid = myPid % numProcXDirection;
  int ypid = myPid / numProcXDirection;

  if    (xpid < NBigXDirection) startx = xpid*(PerProcSmallXDirection+1);
  else  startx = NBigXDirection*(PerProcSmallXDirection+1) + (xpid-NBigXDirection)*PerProcSmallXDirection;
  endx = startx + PerProcSmallXDirection;
  if ( xpid < NBigXDirection) endx++;

  if    (ypid < NBigYDirection) starty = ypid*(PerProcSmallYDirection+1);
  else  starty = NBigYDirection*(PerProcSmallYDirection+1) + (ypid-NBigYDirection)*PerProcSmallYDirection;
  endy = starty + PerProcSmallYDirection;
  if ( ypid < NBigYDirection) endy++;

  int NumMyElements = (endx - startx) * (endy - starty);
  vector<int> MyGlobalElements(NumMyElements);
  int count = 0;

  double *x_coords, *y_coords;
  x_coords= (double *) malloc(sizeof(double)*(endx-startx+2)*(endy-starty+2));
  y_coords= (double *) malloc(sizeof(double)*(endx-startx+2)*(endy-starty+2));

  for (int i = startx ; i < endx ; ++i) {
    for (int j = starty ; j < endy ; ++j) {
      x_coords[count] = ((double) i)*.001;
      y_coords[count] = (double) j;
      MyGlobalElements[count++] = i + j * nx;
    }
  }

  Epetra_Map Map(numNodes,  NumMyElements, &MyGlobalElements[0], 0, Comm);


  // Define the matrix A from the discretization of Laplace equation 
  // with homogeneous Dirichlet condition (using square Q1 elements)

  Epetra_CrsMatrix A(Copy, Map, 0);
  for (int iY = 0; iY < numElePerDirection - 1; ++iY) {
    for (int iX = 0; iX < numElePerDirection - 1; ++iX) {

      int node = iX + iY*(numElePerDirection-1);
      if (Map.LID(node) == -1)
        continue;

      double val = 12.01;
      int pos = node;
      A.InsertGlobalValues(node, 1, &val, &pos);
      if (iY > 0) {
        val = -2.0;
        pos = iX + (iY-1)*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
      }
      else {
        val = -2.0;
        pos = iX + (numElePerDirection-2)*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
      }
      if (iY < numElePerDirection-2) {
        val = -2.0;
        pos = iX + (iY+1)*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
      }
      else {
        val = -2.0;
        pos = iX + (0)*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
      }

      // Insert values on the left side of stencil
      if (iX > 0) {
        val = -2.0;
        pos = iX-1 + iY*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
        if (iY > 0) {
          val = -1.0;
          pos = iX-1 + (iY-1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        else {
          val = -1.0;
          pos = iX-1 + (numElePerDirection-2)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        if (iY < numElePerDirection-2) {
          val = -1.0;
          pos = iX-1 + (iY+1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        else {
          val = -1.0;
          pos = iX-1 + (0)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
      }
      else {
        val = -2.0;
        pos = numElePerDirection-2 + iY*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
        if (iY > 0) {
          val = -1.0;
          pos = numElePerDirection-2 + (iY-1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        else {
          val = -1.0;
          pos = numElePerDirection-2 + (numElePerDirection-2)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        if (iY < numElePerDirection-2) {
          val = -1.0;
          pos = numElePerDirection-2 + (iY+1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        else {
          val = -1.0;
          pos = numElePerDirection-2 + (0)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
      }

      // Insert values on the right side of stencil
      if (iX < numElePerDirection - 2) {
        val = -2.0;
        pos = iX+1 + iY*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
        if (iY > 0) {
          val = -1.0;
          pos = iX+1 + (iY-1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        else {
          val = -1.0;
          pos = iX+1 + (numElePerDirection-2)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        if (iY < numElePerDirection-2) {
          val = -1.0;
          pos = iX+1 + (iY+1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        else {
          val = -1.0;
          pos = iX+1 + (0)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
      }
      else {
        val = -2.0;
        pos = 0 + iY*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
        if (iY > 0) {
          val = -1.0;
          pos = 0 + (iY-1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        else {
          val = -1.0;
          pos = 0 + (numElePerDirection-2)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        if (iY < numElePerDirection-2) {
          val = -1.0;
          pos = 0 + (iY+1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        else {
          val = -1.0;
          pos = 0 + (0)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
      }

    } // for (int iX = 0; iX < numElePerDirection - 1; ++iX)
  } // for (int iY = 0; iY < numElePerDirection - 1; ++iY)
  A.FillComplete();
  A.OptimizeStorage();

  Epetra_Vector LHS(Map);       // solution vector
  Epetra_Vector LHSexact(Map);  // exact solution, check later
  Epetra_Vector RHS(Map);       // right-hand side
  LHS.PutScalar(0.0);                       // zero starting solution
  LHSexact.Random();                        // random exact solution

  // create a parameter list for ML options
  ParameterList MLList;

  ML_Epetra::SetDefaults("SA", MLList);
  MLList.set("output", 10);
  MLList.set("max levels", 10);
  MLList.set("smoother: type", "symmetric Gauss-Seidel");
  MLList.set("smoother: sweeps", 1);
  MLList.set("smoother: pre or post", "both");
  MLList.set("coarse: max size",200);

  MLList.set("x-coordinates", x_coords);
  MLList.set("y-coordinates", y_coords);
//  MLList.set("aggregation: aux: enable", true);
//  MLList.set("aggregation: aux: threshold", 0.2);
  MLList.set("repartition: enable", 1);
  MLList.set("repartition: max min ratio", 1.327);
  MLList.set("repartition: min per proc", 10);
  MLList.set("repartition: Zoltan dimensions", 2);


  ML_Epetra::MultiLevelPreconditioner* MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(A, MLList);

  // =========================== end of ML part =============================

  A.Multiply(false,LHSexact,RHS);

  Epetra_LinearProblem Problem(&A,&LHS,&RHS);
  AztecOO solver(Problem);

  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);
  solver.SetAztecOption(AZ_kspace, 3);
  solver.SetPrecOperator(MLPrec);

  // solve with 500 iterations and 1e-12 as tolerance on the
  // relative residual  

  solver.Iterate(500, 1e-8);
 // delete the preconditioner. Do it BEFORE calling MPI_Finalize
  delete MLPrec;


#ifdef ML_MPI
  MPI_Finalize();
#endif

  return 0;

}
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

  puts("Please configure ML with --enable-epetra ");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

#endif /* #if defined(HAVE_ML_EPETRA) */
