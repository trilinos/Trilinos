/////////////////////////////////////////////////////////////
// The goal of this code is to test the multiple right hand sides
// capabilities of ML
/////////////////////////////////////////////////////////////

#include "ml_config.h"
#if defined(HAVE_ML_EPETRA)
#include "Epetra_ConfigDefs.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Time.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "ml_include.h"
#include "ml_epetra_operator.h"
#include "ml_epetra_utils.h"
#include "ml_utils.h"

int main(int argc, char *argv[]) {

  int myPid = 0;
  int numProc = 1;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPid);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  int numElePerDirection = 5;
  int numNodes = (numElePerDirection - 1)*(numElePerDirection - 1);

  Epetra_Map Map(numNodes, 0, Comm);

  // Define the matrix A from the discretization of Laplace equation 
  // with homogeneous Dirichlet condition (using square Q1 elements)

  Epetra_CrsMatrix A(Copy, Map, 0);
  for (int iY = 0; iY < numElePerDirection - 1; ++iY) {
    for (int iX = 0; iX < numElePerDirection - 1; ++iX) {

      int node = iX + iY*(numElePerDirection-1);
      if (Map.LID(node) == -1)
        continue;

      double val = 4.0;
      int pos = node;
      A.InsertGlobalValues(node, 1, &val, &pos);
      if (iY > 0) {
        val = -2.0;
        pos = iX + (iY-1)*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
      }
      if (iY < numElePerDirection-2) {
        val = -2.0;
        pos = iX + (iY+1)*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
      }

      // Insert values on the left side of stencil
      if (iX > 0) {
        val = -2.0;
        pos = iX-1 + iY*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
        if (iY > 0) {
          val = 1.0;
          pos = iX-1 + (iY-1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        if (iY < numElePerDirection-2) {
          val = 1.0;
          pos = iX-1 + (iY+1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
      }

      // Insert values on the right side of stencil
      if (iX < numElePerDirection - 2) {
        val = -2.0;
        pos = iX+1 + iY*(numElePerDirection-1);
        A.InsertGlobalValues(node, 1, &val, &pos);
        if (iY > 0) {
          val = 1.0;
          pos = iX+1 + (iY-1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
        if (iY < numElePerDirection-2) {
          val = 1.0;
          pos = iX+1 + (iY+1)*(numElePerDirection-1);
          A.InsertGlobalValues(node, 1, &val, &pos);
        }
      }

    } // for (int iX = 0; iX < numElePerDirection - 1; ++iX)
  } // for (int iY = 0; iY < numElePerDirection - 1; ++iY)
  A.FillComplete();

  // Define the ML preconditioner
  int AMG_NLevels = 10;
  ML *ml_handle;
  ML_Create(&ml_handle, AMG_NLevels);

  EpetraMatrix2MLMatrix(ml_handle, 0, &A);

  ML_Aggregate *agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_Set_MaxCoarseSize(agg_object, 1);
  ML_Aggregate_Set_Threshold(agg_object, 0.0);
  ML_Aggregate_Set_Threshold(agg_object, 0.0);
  ML_Aggregate_Set_DampingFactor(agg_object, 1.0);
  AMG_NLevels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0, ML_INCREASING, agg_object);

  // Set a smoother for the AMG method
  int degree = 2;
  for (int j = 0; j < AMG_NLevels; j++)
    ML_Gen_Smoother_MLS(ml_handle, j, ML_BOTH, 30., degree);

  ML_Gen_Solver(ml_handle, ML_MGV, 0, AMG_NLevels-1);

  ML_Epetra::MultiLevelOperator Prec(ml_handle, Comm, Map, Map);
  
  Epetra_CrsMatrix *Mat;
  
  double time1 = 0.0;
  int numNZ = 0;
  
  ML_Operator2EpetraCrsMatrix(ml_handle->Rmat, Mat, numNZ, true, time1);
  
  cout << (*Mat);

  ML_Operator_Print(ml_handle->Rmat, "Rmat");

#ifdef EPETRA_MPI
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

  puts("Please configure ML with --enable-epetra.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
#endif //if defined(HAVE_ML_EPETRA)

