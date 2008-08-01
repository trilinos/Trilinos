/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_include.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Galeri;

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int nx;
  if (argc > 2) nx = (int) strtol(argv[1],NULL,10);
  else          nx = 256;

  nx = nx*(int)sqrt((double)Comm.NumProc());
  int ny = nx;

  char xmlFile[80];
  bool readXML=false;
  if (argc > 2) {strcpy(xmlFile,argv[2]); readXML = true;}
  else sprintf(xmlFile,"params.xml");

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);

  Epetra_MultiVector *coords = CreateCartesianCoordinates("2D", Map,GaleriList);
  double *x_coord=0,*y_coord=0,*z_coord=0;
  double **ttt;
  if (!coords->ExtractView(&ttt)) {
    x_coord = ttt[0];
    y_coord = ttt[1];
  } else {
    if (!Comm.MyPID()) printf("Error extracting coordinate vectors\n");
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }
    
  Epetra_Vector LHS(*Map); LHS.Random();
  Epetra_Vector RHS(*Map); RHS.PutScalar(0.0);
  Epetra_LinearProblem Problem(A, &LHS, &RHS);
  AztecOO solver(Problem);

  // =========================== begin of ML part ===========================
  
  ParameterList MLList;

  if (argc > 1) {
    MLList.set("read XML",true);
    MLList.set("XML input file",argv[1]);
  }
  else {
    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("smoother: type","Chebyshev");
    MLList.set("smoother: sweeps",3);
    MLList.set("coarse: max size",1);
  }

  MLList.set("x-coordinates", x_coord);
  MLList.set("y-coordinates", y_coord);
  MLList.set("z-coordinates", z_coord);
  MLList.set("ML output",10);

  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // =========================== end of ML part =============================
  
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);
  const ML* mlptr = MLPrec->GetML();
  mlptr->Amat[0].apply_without_comm_time = 0;
  mlptr->Amat[0].apply_time = 0;
  solver.Iterate(50, 1e-8);

  // destroy the preconditioner
  delete MLPrec;
  
  // compute the real residual

  double residual;
  LHS.Norm2(&residual);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
  }

  delete A;
  delete Map;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
