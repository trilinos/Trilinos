/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_include.h"
#if defined(HAVE_ML_IFPACK) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRA)
#include "ml_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"

#include "Epetra_Map.h" 
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "ml_ifpack.h"
#include "ml_ifpack_wrap.h"
#include "ml_RowMatrix.h"
#include "Ifpack.h"

#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

using namespace ML_Epetra;

// ================================================ ====== ==== ==== == =

int ML_Ifpack_Gen(ML *ml, int curr_level, int * options,
		  double * params, void ** Ifpack_Handle)
{

#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  ML_Operator *Ke = &(ml->Amat[curr_level]);

  RowMatrix* Ifpack_Matrix = new RowMatrix(Ke, Comm);
  assert (Ifpack_Matrix != 0);

  Ifpack Factory;
  Ifpack_Preconditioner* Prec;

  // set these values in the List
  Teuchos::ParameterList List;

  string Type;
  string AmesosType;

  // local preconditioner
  switch (options[ML_IFPACK_TYPE]) {
  case ML_IFPACK_AMESOS:
    Type = "Amesos";
    break;
  case ML_IFPACK_JACOBI:
    Type = "point relaxation";
    List.set("point: type","Jacobi");
    break;
  case ML_IFPACK_GS:
    Type = "point relaxation";
    List.set("point: type","Gauss-Seidel");
    break;
  case ML_IFPACK_SGS:
    Type = "point relaxation";
    List.set("point: type","symmetric Gauss-Seidel");
    break;
  case ML_IFPACK_BLOCK_JACOBI:
  case ML_IFPACK_BLOCK_JACOBI_AMESOS:
    Type = "block relaxation";
    List.set("block: type","Jacobi");
    break;
  case ML_IFPACK_BLOCK_GS:
  case ML_IFPACK_BLOCK_GS_AMESOS:
    Type = "block relaxation";
    List.set("block: type","Gauss-Seidel");
    break;
  case ML_IFPACK_BLOCK_SGS:
  case ML_IFPACK_BLOCK_SGS_AMESOS:
    Type = "block relaxation";
    List.set("block: type","symmetric Gauss-Seidel");
    break;
  case ML_IFPACK_ICT:
    Type = "ICT";
    break;
  case ML_IFPACK_RILUK:
    Type = "RILUK";
    break;
  default:
    cerr << "Value for options[ML_IFPACK_TYPE] not recognized." << endl;
    exit(EXIT_FAILURE);
  }

  // overlap
  int Overlap = options[ML_IFPACK_OVERLAP];

  // overlap among blocks (only for Ifpack_Jacobi)
  int BlockOverlap = options[ML_IFPACK_BLOCK_OVERLAP];

  // level-of-fill
  int LevelOfFill = options[ML_IFPACK_LEVEL_OF_FILL];

  // number of local blocks
  int LocalParts = options[ML_IFPACK_LOCAL_PARTS];

  // sweeps
  int NumSweeps = options[ML_IFPACK_SWEEPS];

  // damping factor
  double DampingFactor = params[ML_IFPACK_DAMPING_FACTOR];

  List.set("partitioner: local parts", LocalParts);
  List.set("partitioner: overlap", BlockOverlap);
  List.set("point: damping factor", LocalParts);
  List.set("point: sweeps", LocalParts);
  List.set("block: damping factor", LocalParts);
  List.set("block: sweeps", LocalParts);
  List.set("fact: level-of-fill", LocalParts);

  // create the preconditioner 
  Prec = Factory.Create(Type, Ifpack_Matrix, Overlap);
  assert (Prec != 0);

  Prec->SetParameters(List);
  ML_CHK_ERR(Prec->Initialize());
  ML_CHK_ERR(Prec->Compute());

  *Ifpack_Handle = (void *)Prec;

  return 0;
  
} /* ML_Ifpack_Gen */

// ================================================ ====== ==== ==== == =

int ML_Ifpack_Solve(void * Ifpack_Handle, double * x, double * rhs )
{

  Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner *)Ifpack_Handle;

  Epetra_Vector Erhs(View, Prec->OperatorRangeMap(), rhs);
  Epetra_Vector Ex(View, Prec->OperatorDomainMap(), x);
  Prec->ApplyInverse(Erhs,Ex); 

  return 0;

} /* ML_Ifpack_Solve */

// ================================================ ====== ==== ==== == =

void ML_Ifpack_Destroy(void * Ifpack_Handle)
{

  Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner *)Ifpack_Handle;
  delete &(Prec->Matrix());
  delete Prec;

} /* ML_Ifpack_Destroy */

#else

#endif /* #ifdef HAVE_ML_IFPACK */
