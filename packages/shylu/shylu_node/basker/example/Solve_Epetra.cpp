// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_MpiComm.h>
#include <EpetraExt_CrsMatrixIn.h>

#include "Amesos2.hpp"
#include "Amesos2_Meta.hpp"

#include <Kokkos_Core.hpp>

using namespace std;

int main(int argc, char *argv[])
{
  typedef int            LO;
  typedef int            GO;
  typedef double         VAL;

  typedef Kokkos::OpenMP         Node;
  typedef Epetra_CrsMatrix       MAT;
  typedef Epetra_MultiVector     VEC;


  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const Epetra_MpiComm comm (MPI_COMM_WORLD);


  Kokkos::initialize(argc,argv);

  printf(" Basker Solver_Epetra: After init \n");
  
  std::string solver_name = string(argv[1]);
  char *matrix_name = argv[2];

  //Does it exist?
  if(!Amesos2::query(solver_name))
  {
    printf("Does not exist \n");
    return 0;
  }

  printf(" Basker Solver_Epetra: After query\n");

  MAT* A;
  int ret = EpetraExt::MatrixMarketFileToCrsMatrix(matrix_name, comm, A, false, false);

  const Epetra_Map dmnmap = A->DomainMap();
  const Epetra_Map rngmap = A->RangeMap();
  
  //Create random x
  Teuchos::RCP<VEC> X = Teuchos::rcp(new VEC(dmnmap,1));
  X->Random();
  Teuchos::RCP<VEC> B = Teuchos::rcp(new VEC(rngmap,1));
  B->Random();

  printf(" Basker Solver_Epetra: After load\n");
  
  Teuchos::RCP< Amesos2::Solver<MAT,VEC> > solver;
  try
  {
    solver = Amesos2::create<MAT,VEC>("Basker", Teuchos::rcp(A),
  }
  catch(std::invalid_argument e)
  {
    std::cout << e.what() << std::endl;
    return 0;
  }

  printf(" Basker Solver_Epetra: After solver create\n");

  solver->solve();

  printf("Done \n");
  
  Kokkos::finalize();

}
