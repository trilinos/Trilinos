/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_config.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#include "Galeri_Utils.h"

using namespace Teuchos;
using namespace Galeri;
using namespace ML_Epetra;

//
// \author Marzio Sala, SNL 9214
//

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int nx;
  if (argc > 1)
    nx = (int) strtol(argv[1],NULL,10);
  else
    nx = 256;
  int ny = nx * Comm.NumProc(); // each subdomain is a square

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  int NumNodes = nx*ny;
  int NumPDEEqns = 2;

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* CrsA = CreateCrsMatrix("Laplace2D", Map, GaleriList);
  Epetra_VbrMatrix* A = CreateVbrMatrix(CrsA, NumPDEEqns);

  Epetra_Vector LHS(A->DomainMap()); LHS.PutScalar(0);
  Epetra_Vector RHS(A->DomainMap()); RHS.Random();
  Epetra_LinearProblem Problem(A, &LHS, &RHS);
  AztecOO solver(Problem);
  double *x_coord = 0, *y_coord = 0, *z_coord = 0;
  
  Epetra_MultiVector *coords = CreateCartesianCoordinates("2D", &(CrsA->Map()),
                                                          GaleriList);

  double **ttt;
  if (!coords->ExtractView(&ttt)) {
    x_coord = ttt[0];
    y_coord = ttt[1];
  } else {
    printf("Error extracting coordinate vectors\n");
#   ifdef HAVE_MPI
    MPI_Finalize() ;
#   endif
    exit(EXIT_FAILURE);
  }

  ParameterList MLList;
  SetDefaults("SA",MLList);
  MLList.set("ML output",10);
  MLList.set("max levels",10);
  MLList.set("increasing or decreasing","increasing");
  MLList.set("smoother: type", "Chebyshev");
  MLList.set("smoother: sweeps", 3);

  // *) if a low number, it will use all the available processes
  // *) if a big number, it will use only processor 0 on the next level
  MLList.set("aggregation: next-level aggregates per process", 1);

  MLList.set("aggregation: type (level 0)", "Zoltan");
  MLList.set("aggregation: type (level 1)", "Uncoupled");
  MLList.set("aggregation: type (level 2)", "Zoltan");
  MLList.set("aggregation: smoothing sweeps", 2);

  MLList.set("x-coordinates", x_coord);
  MLList.set("y-coordinates", y_coord);
  MLList.set("z-coordinates", z_coord);

  // specify the reduction with respect to the previous level
  // (very small values can break the code)
  int ratio = 16;
  MLList.set("aggregation: global aggregates (level 0)", 
             NumNodes / ratio);
  MLList.set("aggregation: global aggregates (level 1)", 
             NumNodes / (ratio * ratio));
  MLList.set("aggregation: global aggregates (level 2)", 
             NumNodes / (ratio * ratio * ratio));

  MultiLevelPreconditioner* MLPrec = 
    new MultiLevelPreconditioner(*A, MLList, true);

  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 1);
  solver.Iterate(100, 1e-12);
  
  // compute the real residual
  Epetra_Vector Residual(A->DomainMap());
  //1.0 * RHS + 0.0 * RHS - 1.0 * (A * LHS)
  A->Apply(LHS,Residual);
  Residual.Update(1.0, RHS, 0.0, RHS, -1.0);
  double rn;
  Residual.Norm2(&rn);
  
  if (Comm.MyPID() == 0 )
    std::cout << "||b-Ax||_2 = " << rn << endl;

  if (Comm.MyPID() == 0 && rn > 1e-5) {
    std::cout << "TEST FAILED!!!!" << endl;
#   ifdef HAVE_MPI
    MPI_Finalize() ;
#   endif
    exit(EXIT_FAILURE);
  }

  delete MLPrec;
  delete coords;
  delete Map;
  delete CrsA;
  delete A;

  if (Comm.MyPID() == 0)
    std::cout << "TEST PASSED" << endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  exit(EXIT_SUCCESS);
  
}
