#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_ZOLTAN)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"

// includes required by ML
#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "ml_epetra_utils.h"
#include "ml_MultiLevelOperator.h"
#include "Teuchos_ParameterList.hpp"

using namespace ML_Epetra;
using namespace Teuchos;
using namespace Trilinos_Util;

extern "C" {
int ML_DecomposeGraph_with_Zoltan(ML_Operator *Amatrix,
				  int N_parts,
				  int graph_decomposition[],
				  double bdry_nodes[], double old_x[], 
				  double old_y[], double old_z[],
				  int current_level);
}


int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  int nx = 4;
  int ny = nx;
  Epetra_Time Time(Comm);
  
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size",nx * ny);

  // get pointer to the linear system matrix
  Epetra_CrsMatrix * A = Gallery.GetMatrix();

  // get a pointer to the map
  const Epetra_Map * Map = Gallery.GetMap();

  double* x = 0;
  double* y = 0;
  double* z = 0;
  Gallery.GetCartesianCoordinates(x,y,z);

  ML_Comm* MLComm;
  ML_Comm_Create(&MLComm);
  ML_Operator* M = ML_Operator_Create(MLComm);
  Epetra2MLMatrix(A,M);

  vector<int> graph_decomposition(A->NumMyRows());

  ML_DecomposeGraph_with_Zoltan(M,Comm.NumProc(),
                                &graph_decomposition[0],
                                0, x, y, z, 0);
  for (int i = 0 ; i < A->NumMyRows() ; ++i) {
    printf("(%f,%f) ", x[i], y[i]);
    cout << "graph[" << A->Map().GID(i) << "] = " << graph_decomposition[i] << endl;
  }

  ML_Operator_Destroy(&M);
  ML_Comm_Destroy(&MLComm);

  delete [] x;
  delete [] y;
  delete [] z;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return(0);
  
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */

