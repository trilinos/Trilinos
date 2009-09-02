
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

// Goal of this example is to present the basic usage of
// the ML_Epetra::MultiLevelPreconditioner class.
// The example builds a simple matrix and solves the corresponding
// linear system using AztecOO and ML as a preconditioner. It finally
// checks the accuracy of the computed solution.
//
// \author Marzio Sala, ETHZ/COLAB
//
// \data Last modified on 28-Oct-05

#include "ml_include.h"

// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), requires Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example also
// requires --enable-galeri (for the definition of the linear systems)
// and --enable-aztecoo (to solve the linear system)

// epetra objects
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
// required to build the example matrix
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
// required by the linear system solver
#include "AztecOO.h"
// required by ML
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Galeri;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

#define ML_SCALING
#ifdef ML_SCALING
   const int ntimers=4;
   enum {total, probBuild, precBuild, solve};
   ml_DblLoc timeVec[ntimers], maxTime[ntimers], minTime[ntimers];

  for (int i=0; i<ntimers; i++) timeVec[i].rank = Comm.MyPID();
  timeVec[total].value = MPI_Wtime();
#endif

  int nx;
  if (argc > 1) nx = (int) strtol(argv[1],NULL,10);
  else          nx = 256;

  if (nx < 1) nx = 256; // input a nonpositive integer if you want to specify
                        // the XML input file name.
  nx = nx*(int)sqrt((double)Comm.NumProc());
  int ny = nx;

  printf("nx = %d\nny = %d\n",nx,ny);
  fflush(stdout);

  char xmlFile[80];
  bool readXML=false;
  if (argc > 2) {strcpy(xmlFile,argv[2]); readXML = true;}
  else sprintf(xmlFile,"%s","params.xml");

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);

#ifdef ML_SCALING
  timeVec[probBuild].value = MPI_Wtime();
#endif
  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);

  if (!Comm.MyPID()) printf("nx = %d, ny = %d, mx = %d, my = %d\n",nx,ny,GaleriList.get("mx",-1),GaleriList.get("my",-1));
  fflush(stdout);
  //avoid potential overflow
  double numMyRows = A->NumMyRows();
  double numGlobalRows;
  Comm.SumAll(&numMyRows,&numGlobalRows,1);
  if (!Comm.MyPID()) printf("# global rows = %1.0f\n",numGlobalRows);
  //printf("pid %d: #rows = %d\n",Comm.MyPID(),A->NumMyRows());
  fflush(stdout);

  Epetra_MultiVector *coords = CreateCartesianCoordinates("2D", Map,GaleriList);
  double *x_coord=0,*y_coord=0,*z_coord=0;
  double **ttt;
  if (!coords->ExtractView(&ttt)) {
    x_coord = ttt[0];
    y_coord = ttt[1];
  } else {
    if (!Comm.MyPID()) printf("Error extracting coordinate vectors\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
    
  Epetra_Vector LHS(*Map); LHS.Random();
  Epetra_Vector RHS(*Map); RHS.PutScalar(0.0);
  Epetra_LinearProblem Problem(A, &LHS, &RHS);
  AztecOO solver(Problem);

#ifdef ML_SCALING
  timeVec[probBuild].value = MPI_Wtime() - timeVec[probBuild].value;
#endif

  // =========================== begin of ML part ===========================
  
#ifdef ML_SCALING
  timeVec[precBuild].value = MPI_Wtime();
#endif
  ParameterList MLList;

  if (readXML) {
    MLList.set("read XML",true);
    MLList.set("XML input file",xmlFile);
  }
  else {
    cout << "here" << endl;
    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("smoother: type","Chebyshev");
    MLList.set("smoother: sweeps",3);
    MLList.set("coarse: max size",1);
  }

  MLList.set("x-coordinates", x_coord);
  MLList.set("y-coordinates", y_coord);
  MLList.set("z-coordinates", z_coord);

/*
RCP<std::vector<int> >
   m_smootherAztecOptions = rcp(new std::vector<int>(AZ_OPTIONS_SIZE));
RCP<std::vector<double> >
   m_smootherAztecParams = rcp(new std::vector<double>(AZ_PARAMS_SIZE));
//int             m_smootherAztecOptions[AZ_OPTIONS_SIZE];
//double          m_smootherAztecParams[AZ_PARAMS_SIZE];
  
std::string smootherType("Aztec");
AZ_defaults(&(*m_smootherAztecOptions)[0],&(*m_smootherAztecParams)[0]);
(*m_smootherAztecOptions)[AZ_precond]         = AZ_dom_decomp;
(*m_smootherAztecOptions)[AZ_subdomain_solve] = AZ_icc;
bool smootherAztecAsSolver = true;
*/
  MLList.set("ML output",10);

  MLList.set("repartition: enable",1);
  MLList.set("repartition: max min ratio",1.3);
  MLList.set("repartition: min per proc",200);
  MLList.set("repartition: partitioner","Zoltan");
  MLList.set("repartition: Zoltan dimensions",2);
  MLList.set("repartition: put on single proc",1);
  MLList.set("repartition: Zoltan type","hypergraph");
  MLList.set("repartition: estimated iterations",13);

/*
MLList.set("smoother: Aztec options",m_smootherAztecOptions);
MLList.set("smoother: Aztec params",m_smootherAztecParams);
MLList.set("smoother: type",smootherType.c_str());
MLList.set("smoother: Aztec as solver", smootherAztecAsSolver); 

MLList.set("ML print initial list", 0); 
MLList.set("ML print final list", 0); 
*/

  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // verify unused parameters on process 0 (put -1 to print on all
  // processes)
  MLPrec->PrintUnused(0);
#ifdef ML_SCALING
  timeVec[precBuild].value = MPI_Wtime() - timeVec[precBuild].value;
#endif

  // =========================== end of ML part =============================
  
#ifdef ML_SCALING
  timeVec[solve].value = MPI_Wtime();
#endif
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);
  solver.Iterate(500, 1e-12);
#ifdef ML_SCALING
  timeVec[solve].value = MPI_Wtime() - timeVec[solve].value;
#endif

  // destroy the preconditioner
  delete MLPrec;
  
  // compute the real residual

  double residual;
  LHS.Norm2(&residual);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
  }

  // for testing purposes
  if (residual > 1e-5)
    exit(EXIT_FAILURE);

  delete A;
  delete Map;
  delete coords;

#ifdef ML_SCALING
  timeVec[total].value = MPI_Wtime() - timeVec[total].value;

  //avg
  double dupTime[ntimers],avgTime[ntimers];
  for (int i=0; i<ntimers; i++) dupTime[i] = timeVec[i].value;
  MPI_Reduce(dupTime,avgTime,ntimers,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  for (int i=0; i<ntimers; i++) avgTime[i] = avgTime[i]/Comm.NumProc();
  //min
  MPI_Reduce(timeVec,minTime,ntimers,MPI_DOUBLE_INT,MPI_MINLOC,0,MPI_COMM_WORLD);
  //max
  MPI_Reduce(timeVec,maxTime,ntimers,MPI_DOUBLE_INT,MPI_MAXLOC,0,MPI_COMM_WORLD);

  if (Comm.MyPID() == 0) {
    printf("timing :  max (pid)  min (pid)  avg\n");
    printf("Problem build         :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[probBuild].value,maxTime[probBuild].rank,
             minTime[probBuild].value,minTime[probBuild].rank,
             avgTime[probBuild]);
    printf("Preconditioner build  :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[precBuild].value,maxTime[precBuild].rank,
             minTime[precBuild].value,minTime[precBuild].rank,
             avgTime[precBuild]);
    printf("Solve                 :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[solve].value,maxTime[solve].rank,
             minTime[solve].value,minTime[solve].rank,
             avgTime[solve]);
    printf("Total                 :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[total].value,maxTime[total].rank,
             minTime[total].value,minTime[total].rank,
             avgTime[total]);
  }
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
