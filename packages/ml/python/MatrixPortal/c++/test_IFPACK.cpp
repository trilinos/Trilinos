#include "Ifpack_ConfigDefs.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"
#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::ParameterList GaleriList;

  // The problem is defined on a 2D grid, global size is nx * nx.
  int nx = 30; 
  GaleriList.set("n", nx * nx);
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx);
  Epetra_Map* Map = Galeri::CreateMap("Linear", Comm, GaleriList);
  Epetra_RowMatrix* A = Galeri::CreateCrsMatrix("Laplace2D", Map, GaleriList);

  // =============================================================== //
  // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
  // =============================================================== //

  Teuchos::ParameterList List;

  // allocates an IFPACK factory. No data is associated 
  // to this object (only method Create()).
  Ifpack Factory;

  // create the preconditioner. For valid PrecType values,
  // please check the documentation
  string PrecType = "ILU"; // incomplete LU
  int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1,
                        // it is ignored.

  Ifpack_Preconditioner* Prec = Factory.Create(PrecType, A, OverlapLevel);
  assert(Prec != 0);

  // specify parameters for ILU
  List.set("fact: drop tolerance", 1e-9);
  List.set("fact: level-of-fill", 1);
  // the combine mode is on the following:
  // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
  // Their meaning is as defined in file Epetra_CombineMode.h   
  List.set("schwarz: combine mode", "Add");
  // sets the parameters
  IFPACK_CHK_ERR(Prec->SetParameters(List));

  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  IFPACK_CHK_ERR(Prec->Initialize());

  // Builds the preconditioners, by looking for the values of 
  // the matrix.
  IFPACK_CHK_ERR(Prec->Compute());

  // =================================================== //
  // E N D   O F   I F P A C K   C O N S T R U C T I O N //
  // =================================================== //

  // At this point, we need some additional objects
  // to define and solve the linear system.

  // defines LHS and RHS
  Epetra_Vector LHS(A->OperatorDomainMap());
  Epetra_Vector RHS(A->OperatorDomainMap());

  // solution is constant
  LHS.PutScalar(1.0);
  // now build corresponding RHS
  A->Apply(LHS,RHS);

  // now randomize the solution
  RHS.Random();

  // need an Epetra_LinearProblem to define AztecOO solver
  Epetra_LinearProblem Problem(A,&LHS,&RHS);

  // now we can allocate the AztecOO solver
  AztecOO Solver(Problem);

  // specify solver
  Solver.SetAztecOption(AZ_solver,AZ_gmres);
  Solver.SetAztecOption(AZ_output,32);

  // HERE WE SET THE IFPACK PRECONDITIONER
  Solver.SetPrecOperator(Prec);

  // .. and here we solve
  Solver.Iterate(1550,1e-8);

  cout << *Prec;

  // delete memory
  delete Prec;
  delete Map;
  delete A;
  
#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  return(EXIT_SUCCESS);
}
