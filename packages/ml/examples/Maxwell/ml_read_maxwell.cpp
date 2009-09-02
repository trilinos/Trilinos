/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
//@HEADER
/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Branch:           $Branch$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

//#define CurlCurlAndMassAreSeparate

/*
   Sample driver for Maxwell equation AMG solver in the ML package. This
   example reads in data from a file.  All data must be in the MatrixMarket
   format.  (The EpetraExt documentation describes the format for maps.)  This
   example can be compiled and used in two different ways:

   -----------------------------------------------------------------------
   USAGE CASE 1 (default): curl,curl and mass are provided as one matrix
   -----------------------------------------------------------------------

   By default, it's assumed that the edge FE input matrix is curl,curl + mass.
   In this case, invoke the example as follows:

      ml_read_maxwell.exe Ke T Kn [edge map] [node map]
  
   where

      Ke is the edge FE matrix (curlcurl + mass)
      T is the topological gradient matrix
      Kn is the nodal FE matrix
      edge map (optional)
      node map (optional)

   -----------------------------------------------------------------------
   USAGE CASE 2:  curl,curl and mass matrices are separate
   -----------------------------------------------------------------------

   If the macro CurlCurlAndMassAreSeparate is defined, then the edge FE
   (curl,curl) and mass matrices are read in separately.  In this case, invoke
   the example as follows:

      ml_read_maxwell.exe S M T Kn [edge map] [node map]
  
   where

      S is the curl curl matrix
      M is the mass matrix
      T is the discrete gradient matrix
      Kn is the nodal FE matrix
      edge map (optional)
      node map (optional)
  
   Matrices are read from file via the local function 

     int MatrixMarketFileToCrsMatrix(const char *filename,
                                  const Epetra_Map & rowMap,
                                  const Epetra_Map& rangeMap,
                                  const Epetra_Map& domainMap,
                                  Epetra_CrsMatrix * & A)

   and the row maps with

       EpetraExt::MatrixMarketFileToMap(datafile, Comm, nodeMap).

   In this function the domainMap is calculated on the fly in order to ensure
   the same number of non zeros independent of the number of processors used.

*/

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Export.h"
#include "AztecOO.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

//read in Crs Matrix and calculate domainMap on the fly
int MatrixMarketFileToCrsMatrix(const char *filename,
                                const Epetra_Map & rowMap,
                                const Epetra_Map& rangeMap,
                                const Epetra_Map& domainMap,
                                Epetra_CrsMatrix * & A)
{
  A = new Epetra_CrsMatrix(Copy, rowMap, 0);
  return(EpetraExt::MatrixMarketFileToCrsMatrixHandle(filename, A->Comm(), A,
                                                      &rowMap,NULL,
                                                      &rangeMap, &domainMap));
} 

ML_Comm *mlcomm;

int main(int argc, char *argv[])
{

#ifdef ML_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  ML_Comm_Create(&mlcomm);

  char *datafile;

#ifdef CurlCurlAndMassAreSeparate
  if (argc != 5 && argc != 7) {
    if (Comm.MyPID() == 0) {
      cout << "usage: ml_maxwell.exe <S> <M> <T> <Kn> [edge map] [node map]"
           << endl;
      cout << "        S = edge stiffness matrix file" << endl;
      cout << "        M = edge mass matrix file" << endl;
      cout << "        T = discrete gradient file" << endl;
      cout << "       Kn = auxiliary nodal FE matrix file" << endl;
      cout << " edge map = edge distribution over processors" << endl;
      cout << " node map = node distribution over processors" << endl;
      cout << argc << endl;
    }
#else //ifdef CurlCurlAndMassAreSeparate
  if (argc != 4 && argc != 6) {
    if (Comm.MyPID() == 0) {
      cout << "usage: ml_maxwell.exe <A> <T> <Kn> [edge map] [node map]" <<endl;
      cout << "        A = edge element matrix file" << endl;
      cout << "        T = discrete gradient file" << endl;
      cout << "       Kn = auxiliary nodal FE matrix file" << endl;
      cout << " edge map = edge distribution over processors" << endl;
      cout << " node map = node distribution over processors" << endl;
      cout << argc << endl;
    }
#endif //ifdef CurlCurlAndMassAreSeparate
#ifdef ML_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  Epetra_Map *edgeMap, *nodeMap;
  Epetra_CrsMatrix *CCplusM=NULL, *CurlCurl=NULL, *Mass=NULL, *T=NULL, *Kn=NULL;

  // ================================================= //
  // READ IN MAPS FROM FILE                            //
  // ================================================= //
  // every processor reads this in
#ifdef CurlCurlAndMassAreSeparate
  if (argc > 5)
#else
  if (argc > 4)
#endif
  {
    datafile = argv[5];
    if (Comm.MyPID() == 0) {
      printf("Reading in edge map from %s ...\n",datafile);
      fflush(stdout);
    }
    EpetraExt::MatrixMarketFileToMap(datafile, Comm, edgeMap);
    datafile = argv[6];
    if (Comm.MyPID() == 0) {
      printf("Reading in node map from %s ...\n",datafile);
      fflush(stdout);
    }
    EpetraExt::MatrixMarketFileToMap(datafile, Comm, nodeMap);
  }
  else { // linear maps
         // Read the T matrix to determine the map sizes
         // and then construct linear maps

    if (Comm.MyPID() == 0)
      printf("Using linear edge and node maps ...\n");

    const int lineLength = 1025;
    char line[lineLength];
    FILE *handle;
    int M,N,NZ;
#ifdef CurlCurlAndMassAreSeparate
     handle = fopen(argv[3],"r"); 
#else
     handle = fopen(argv[2],"r"); 
#endif
    if (handle == 0) EPETRA_CHK_ERR(-1); // file not found
    // Strip off header lines (which start with "%")
    do {
       if(fgets(line, lineLength, handle)==0) {if (handle!=0) fclose(handle);}
    } while (line[0] == '%');
    // Get problem dimensions: M, N, NZ
    if(sscanf(line,"%d %d %d", &M, &N, &NZ)==0) {if (handle!=0) fclose(handle);}
    fclose(handle);
    edgeMap = new Epetra_Map(M,0,Comm);
    nodeMap = new Epetra_Map(N,0,Comm); 
  }

  // ===================================================== //
  // READ IN MATRICES FROM FILE                            //
  // ===================================================== //
#ifdef CurlCurlAndMassAreSeparate
  for (int i = 1; i <5; i++) {
    datafile = argv[i];
    if (Comm.MyPID() == 0) {
      printf("reading %s ....\n",datafile); fflush(stdout);
    }
    switch (i) {
    case 1: //Curl
      MatrixMarketFileToCrsMatrix(datafile,*edgeMap,*edgeMap,*edgeMap,CurlCurl);
      break;
    case 2: //Mass
      MatrixMarketFileToCrsMatrix(datafile, *edgeMap, *edgeMap, *edgeMap, Mass);
      break;
    case 3: //Gradient
      MatrixMarketFileToCrsMatrix(datafile, *edgeMap, *edgeMap,*nodeMap, T);
      break;
    case 4: //Auxiliary nodal matrix
      MatrixMarketFileToCrsMatrix(datafile, *nodeMap,*nodeMap, *nodeMap, Kn);
      break;
    } //switch
  } //for (int i = 1; i <5; i++)

#else
  for (int i = 1; i <4; i++) {
    datafile = argv[i];
    if (Comm.MyPID() == 0) {
      printf("reading %s ....\n",datafile); fflush(stdout);
    }
    switch (i) {
    case 1: //Edge element matrix
      MatrixMarketFileToCrsMatrix(datafile,*edgeMap,*edgeMap,*edgeMap,CCplusM);
      break;
    case 2: //Gradient
      MatrixMarketFileToCrsMatrix(datafile, *edgeMap, *edgeMap, *nodeMap, T);
      break;
    case 3: //Auxiliary nodal matrix
      MatrixMarketFileToCrsMatrix(datafile, *nodeMap, *nodeMap, *nodeMap, Kn);
      break;
    } //switch
  } //for (int i = 1; i <4; i++)
#endif //ifdef CurlCurlAndMassAreSeparate

  // ==================================================== //
  // S E T U P   O F    M L   P R E C O N D I T I O N E R //
  // ==================================================== //

  Teuchos::ParameterList MLList;
  int *options    = new int[AZ_OPTIONS_SIZE];
  double *params  = new double[AZ_PARAMS_SIZE];
  ML_Epetra::SetDefaults("maxwell", MLList, options, params);

  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("coarse: max size", 30);
  MLList.set("aggregation: threshold", 0.0);
  //MLList.set("negative conductivity",true);
  //MLList.set("smoother: type", "Jacobi");
  MLList.set("subsmoother: type", "symmetric Gauss-Seidel");
  //MLList.set("max levels", 2);

  // coarse level solve
  MLList.set("coarse: type", "Amesos-KLU");
  //MLList.set("coarse: type", "Hiptmair");
  //MLList.set("coarse: type", "Jacobi");

  //MLList.set("dump matrix: enable", true);

#ifdef CurlCurlAndMassAreSeparate
  //Create the matrix of interest.
  CCplusM = Epetra_MatrixAdd(CurlCurl,Mass,1.0);
#endif

#ifdef CurlCurlAndMassAreSeparate
  ML_Epetra::MultiLevelPreconditioner * MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(*CurlCurl, *Mass, *T, *Kn, MLList);
  // Comment out the line above and uncomment the next one if you have
  // mass and curl separately but want to precondition as if they are added
  // together.
/*
  ML_Epetra::MultiLevelPreconditioner * MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(*CCplusM, *T, *Kn, MLList);
*/
#else
  ML_Epetra::MultiLevelPreconditioner * MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(*CCplusM, *T, *Kn, MLList);
#endif //ifdef CurlCurlAndMassAreSeparate

  MLPrec->PrintUnused(0);

  // ========================================================= //
  // D E F I N I T I O N   O F   A Z T E C O O   P R O B L E M //
  // ========================================================= //

  // create left-hand side and right-hand side, and populate them with
  // data from file. Both vectors are defined on the domain map of the
  // edge matrix.
  // Epetra_Vectors can be created in View mode, to accept pointers to
  // double vectors.

  if (Comm.MyPID() == 0)
    cout << "Putting in a zero initial guess and random rhs (in the range of S+M)" << endl;
  Epetra_Vector x(CCplusM->DomainMap());
  x.Random();
  Epetra_Vector rhs(CCplusM->DomainMap());
  CCplusM->Multiply(false,x,rhs);
  x.PutScalar(0.0);

  double vecnorm;
  rhs.Norm2(&vecnorm);
  if (Comm.MyPID() == 0) cout << "||rhs|| = " << vecnorm << endl;
  x.Norm2(&vecnorm);
  if (Comm.MyPID() == 0) cout << "||x|| = " << vecnorm << endl;

  // for AztecOO, we need an Epetra_LinearProblem
  Epetra_LinearProblem Problem(CCplusM,&x,&rhs);
  // AztecOO Linear problem
  AztecOO solver(Problem);
  // set MLPrec as precondititoning operator for AztecOO linear problem
  //cout << "no ml preconditioner!!!" << endl;
  solver.SetPrecOperator(MLPrec);
  //solver.SetAztecOption(AZ_precond, AZ_Jacobi);

  // a few options for AztecOO
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);

  solver.Iterate(75, 1e-12);

  // =============== //
  // C L E A N   U P //
  // =============== //

  delete MLPrec;    // destroy phase prints out some information
  delete CurlCurl;
  delete CCplusM;
  delete Mass;
  delete T;
  delete Kn;

  ML_Comm_Destroy(&mlcomm);
#ifdef ML_MPI
  MPI_Finalize();
#endif
        
  return 0;
        
} //main

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

  puts("Please configure ML with:");
#if !defined(HAVE_ML_EPETRA)
  puts("--enable-epetra");
#endif
#if !defined(HAVE_ML_TEUCHOS)
  puts("--enable-teuchos");
#endif
#if !defined(HAVE_ML_EPETRAEXT)
  puts("--enable-epetraext");
#endif
#if !defined(HAVE_ML_AZTECOO)
  puts("--enable-aztecoo");
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

#endif
