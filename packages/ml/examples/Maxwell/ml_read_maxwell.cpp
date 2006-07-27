//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Branch:           $Branch$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

/* Sample driver for Maxwell equation AMG solver in the ML package. 
   This example reads in the matrices from a file via the local function 

   int MatrixMarketFileToCrsMatrix(const char *filename,
                                const Epetra_Map & rowMap,
                                const Epetra_Map& rangeMap,
                                const Epetra_Map& domainMap,
				Epetra_CrsMatrix * & A)

   and the row maps with  EpetraExt::MatrixMarketFileToMap(datafile, Comm, nodeMap);
   In this function the domainMap is calculated on the fly in order to assure the same
   number of non zeros independent of the number of processors used.

   The files are all written in the applicaion (fermaXX) with: 

   EpetraExt::RowMatrixToMatrixMarketFile("K.mtx", *K); (K = A - sigma*M)

   EpetraExt::RowMatrixToMatrixMarketFile("Y.mtx", *Y); (Y = gradient matrix)

   EpetraExt::RowMatrixToMatrixMarketFile("H.mtx", *H); (H poisson matrix)

   EpetraExt::BlockMapToMatrixMarketFile("A_rowmap.txt", A->RowMap());
   EpetraExt::BlockMapToMatrixMarketFile("H_rowmap.txt", H->RowMap());

  
   Usage: mpirun -np 8 ../ml_read_maxwell K.mtx Y.mtx H.mtx A_rowmap.txt H_rowmap.txt | tee p1.out


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

/*
    read in Crs Matrix and calculate domainMap on the fly
    this functionality is missing in Epetra_Ext!
*/
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

  if (argc != 4 && argc != 6) {
    if (Comm.MyPID() == 0) {
      cout << "usage: ml_maxwell.exe <A> <T> <Kn> [edge map] [node map]" << endl;
      cout << "        A = edge element matrix file" << endl;
      cout << "        T = discrete gradient file" << endl;
      cout << "       Kn = auxiliary nodal FE matrix file" << endl;
      cout << " edge map = edge distribution over processors" << endl;
      cout << " node map = node distribution over processors" << endl;
      cout << argc << endl;
    }
#ifdef ML_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  Epetra_Map *edgeMap=NULL, *nodeMap=NULL;
  Epetra_CrsMatrix *CCplusM=NULL, *Mass=NULL, *T=NULL, *Kn=NULL;

  // ================================================= //
  // READ IN MAPS FROM FILE                            //
  // ================================================= //
  // every processor reads this in
  if (argc > 4) {
    datafile = argv[4];
    if (Comm.MyPID() == 0) {
      printf("Reading in edge map from %s ...\n",datafile);
      fflush(stdout);
    }
    EpetraExt::MatrixMarketFileToMap(datafile, Comm, edgeMap);
    datafile = argv[5];
    if (Comm.MyPID() == 0) {
      printf("Reading in node map from %s\n",datafile);
      fflush(stdout);
    }
    EpetraExt::MatrixMarketFileToMap(datafile, Comm, nodeMap);
  }
  else { //use linear maps

    edgeMap = new Epetra_Map(2450,0,Comm);
    nodeMap = new Epetra_Map(873,0,Comm);
  }

  // ===================================================== //
  // READ IN MATRICES FROM FILE                            //
  // ===================================================== //

  for (int i = 1; i <4; i++) {
    datafile = argv[i];
    if (Comm.MyPID() == 0) {
      printf("reading %s ....\n",datafile); fflush(stdout);
    }
    switch (i) {
    case 1: //Edge element matrix
      MatrixMarketFileToCrsMatrix(datafile, *edgeMap, *edgeMap, *edgeMap, CCplusM);
      break;
    case 2: //Gradient
      MatrixMarketFileToCrsMatrix(datafile, *edgeMap, *edgeMap, *nodeMap, T);
      break;
    case 3: //Auxiliary nodal matrix
      MatrixMarketFileToCrsMatrix(datafile, *nodeMap, *nodeMap, *nodeMap, Kn);
      break;
    } //switch
  } 

  // ==================================================== //
  // S E T U P   O F    M L   P R E C O N D I T I O N E R //
  // ==================================================== //

  Teuchos::ParameterList MLList;
  int *options    = new int[AZ_OPTIONS_SIZE];
  double *params  = new double[AZ_PARAMS_SIZE];
  ML_Epetra::SetDefaults("maxwell", MLList, options, params);

  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("coarse: max size", 30);
  MLList.set("aggregation: damping factor",1.3333);
  MLList.set("subsmoother: type", "MLS");
  //MLList.set("subsmoother: type", "symmetric Gauss-Seidel");
  MLList.set("output",10);

  ML_Set_PrintLevel(10);

  ML_Epetra::MultiLevelPreconditioner * MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(*CCplusM, *T, *Kn, MLList);

  MLPrec->PrintList(0);
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

  //EpetraExt::RowMatrixToMatrixMarketFile("checkingA.mm",*CCplusM,"curlcurl plus mass");

  // for AztecOO, we need an Epetra_LinearProblem
  //Epetra_CrsMatrix *Combined = Epetra_MatrixAdd(CCplusM,Mass,1.0);
  //EpetraExt::RowMatrixToMatrixMarketFile("CplusM.mm",*Combined,"curlcurl plus mass");
  Epetra_LinearProblem Problem(CCplusM,&x,&rhs);
  // AztecOO Linear problem
  AztecOO solver(Problem);
  // set MLPrec as precondititoning operator for AztecOO linear problem
  solver.SetPrecOperator(MLPrec);

  // a few options for AztecOO
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);

  solver.Iterate(150, 1e-10);

  // =============== //
  // C L E A N   U P //
  // =============== //

  delete CCplusM;
  delete Mass;
  delete T;
  delete Kn;
  delete edgeMap;
  delete nodeMap;
  delete [] options;
  delete [] params;
  ML_Comm_Destroy(&mlcomm);

  delete MLPrec;    // destroy phase prints out some information
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
