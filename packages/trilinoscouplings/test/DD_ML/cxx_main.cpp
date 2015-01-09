
// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
//                 Copyright (2001) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Ifpack_ConfigDefs.h"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Ifpack_CrsIct.h"
#include "Ifpack.h"
#include "Teuchos_RefCountPtr.hpp"
#ifdef JJH_PRINTMATRIX
#include "/data/home/jhu/Trilinos/development-branch/Trilinos/packages/epetraext/src/inout/EpetraExt_RowMatrixOut.h"
#endif

#ifdef HAVE_IFPACK_ML
#include "ml_MultiLevelPreconditioner.h"
#endif

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  int nx, ny, mx, my, overlap;
  if (argc > 1) overlap = (int) strtol(argv[1],NULL,10);
  else          overlap = 0;
  if (argc > 2) nx = (int) strtol(argv[2],NULL,10);
  else          nx = 30;
  if (argc > 3) ny = (int) strtol(argv[3],NULL,10);
  else          ny = nx * Comm.NumProc();
  if (argc > 4) mx = (int) strtol(argv[4],NULL,10);
  else          mx = 1;
  if (argc > 5) my = (int) strtol(argv[5],NULL,10);
  else          my = Comm.NumProc();

  printf("\n--------\noverlap=%d\nnx=%d,ny=%d\nmx=%d,my=%d\n--------\n",
          overlap,nx,ny,mx,my);

  int MyPID = Comm.MyPID();
  bool verbose = false; 
  if (MyPID==0) verbose = true;

  // The problem is defined on a 2D grid, global size is nx * nx.
  Teuchos::ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", mx);
  GaleriList.set("my", my);
/*
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());
*/
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Cartesian2D", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );
  Teuchos::RefCountPtr<Epetra_MultiVector> LHS = Teuchos::rcp( new Epetra_MultiVector(*Map, 1) );
  Teuchos::RefCountPtr<Epetra_MultiVector> RHS = Teuchos::rcp( new Epetra_MultiVector(*Map, 1) );
  //LHS->PutScalar(0.0); RHS->Random();
  LHS->PutScalar(0.0); RHS->PutScalar(1.0);

  if (!MyPID) cout << "#global rows in A = " << A->NumGlobalRows() << endl;

#ifdef JJH_PRINT_LOCAL_MATRIX
  A->Print(cout);
#endif
  
#ifdef JJH_PRINTMATRIX
  char filename[80];
  sprintf(filename,"Acheck_%dproc",Comm.NumProc());
  EpetraExt::RowMatrixToMatlabFile(filename,*A);
#endif

  // ============================ //
  // Construct ML preconditioner  //
  // ---------------------------- //

  int Niters = 500;
  AztecOO solver;
  Ifpack Factory;
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec =
                Teuchos::rcp( Factory.Create("ML", &*A, overlap) );

  Teuchos::ParameterList List;
  ML_Epetra::SetDefaults("SA",List);
  List.set("ML output",0);
  List.set("max levels",1);
  List.set("XML input file","ml_parameters.xml");
  List.set("read XML",true);
  IFPACK_CHK_ERR(Prec->SetParameters(List));
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());

  // Here we create an AztecOO object
  LHS->PutScalar(0.0);

  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_cg);
  solver.SetAztecOption(AZ_conv,AZ_noscaled);
  solver.SetPrecOperator(&*Prec);
  solver.SetAztecOption(AZ_output, 1); 
  solver.Iterate(Niters, 1.0e-8);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
