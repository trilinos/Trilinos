#include "Ifpack_ConfigDefs.h"
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Jacobi.h"
#include "Ifpack_GaussSeidel.h"
#include "Ifpack_SymGaussSeidel.h"
#include "Ifpack_SOR.h"
#include "Ifpack_SSOR.h"
#include "Ifpack_BlockJacobi.h"
#include "Ifpack_BlockGaussSeidel.h"
#include "Ifpack_BlockSymGaussSeidel.h"
#include "Ifpack_Ict.h"
#include "Ifpack_Riluk.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_SparseContainer.h"
#ifdef HAVE_IFPACK_AMESOS
#include "Ifpack_Amesos.h"
#endif

Ifpack_Preconditioner* Ifpack::Create(const string PrecType,
				      Epetra_RowMatrix* Matrix,
				      const int Overlap)
{

  if (PrecType == "Jacobi") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_Jacobi>(Matrix, Overlap));
  }
  else if (PrecType == "Gauss-Seidel") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_GaussSeidel>(Matrix,Overlap));
  }
  else if (PrecType == "symmetric Gauss-Seidel") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_SymGaussSeidel>(Matrix,Overlap));
  }
  if (PrecType == "SOR") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_SOR>(Matrix, Overlap));
  }
  if (PrecType == "SSOR") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_SSOR>(Matrix, Overlap));
  }
  else if (PrecType == "block Jacobi") {
    return(new Ifpack_AdditiveSchwarz<
	    Ifpack_BlockJacobi<Ifpack_DenseContainer> >(Matrix,Overlap));
  }
  else if (PrecType == "block Gauss-Seidel") {
    return(new Ifpack_AdditiveSchwarz<
	    Ifpack_BlockGaussSeidel<Ifpack_DenseContainer> >(Matrix,Overlap));
  }
  else if (PrecType == "block symmetric Gauss-Seidel") {
    return(new Ifpack_AdditiveSchwarz<
	    Ifpack_BlockSymGaussSeidel<Ifpack_DenseContainer> >(Matrix,Overlap));
  }
  else if (PrecType == "block Jacobi (Amesos)") {
    return(new Ifpack_AdditiveSchwarz<
	    Ifpack_BlockJacobi<Ifpack_SparseContainer<Ifpack_Amesos> > >(Matrix,Overlap));
  }
  else if (PrecType == "block Gauss-Seidel (Amesos)") {
    return(new Ifpack_AdditiveSchwarz<
	    Ifpack_BlockGaussSeidel<Ifpack_SparseContainer<Ifpack_Amesos> > >(Matrix,Overlap));
  }
  else if (PrecType == "block symmetric Gauss-Seidel (Amesos)") {
    return(new Ifpack_AdditiveSchwarz<
	    Ifpack_BlockSymGaussSeidel<Ifpack_SparseContainer<Ifpack_Amesos> > >(Matrix,Overlap));
  }
  else if (PrecType == "Amesos") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(Matrix,Overlap));
  }
  else if (PrecType == "ICT") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_Ict>(Matrix,Overlap));

  } 
  else if (PrecType == "RILUK") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_Riluk>(Matrix,Overlap));
  }
  else
    return(0);

}
