#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_TEUCHOS)
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_gIct.h"
#include "Ifpack_vIct.h"
#include "Ifpack_gRiluk.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_SparseContainer.h"
#ifdef HAVE_IFPACK_AMESOS
#include "Ifpack_Amesos.h"
#endif

//==============================================================================
Ifpack_Preconditioner* Ifpack::Create(const string PrecType,
                                      Epetra_RowMatrix* Matrix,
                                      const int Overlap)
{

  if (PrecType == "point relaxation") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation>(Matrix, Overlap));
  }
  else if (PrecType == "block relaxation") {
    return(new Ifpack_AdditiveSchwarz<
            Ifpack_BlockRelaxation<Ifpack_DenseContainer> >(Matrix,Overlap));
  }
#ifdef HAVE_IFPACK_AMESOS
  else if (PrecType == "block relaxation (Amesos)") {
    return(new Ifpack_AdditiveSchwarz<
            Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > >(Matrix,Overlap));
  }
  else if (PrecType == "Amesos") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(Matrix,Overlap));
  }
#endif
  else if (PrecType == "gICT") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_gIct>(Matrix,Overlap));

  } 
  else if (PrecType == "vICT") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_vIct>(Matrix,Overlap));

  } 
  else if (PrecType == "gRILUK") {
    return(new Ifpack_AdditiveSchwarz<Ifpack_gRiluk>(Matrix,Overlap));
  }
  else
    return(0);

}

#endif // HAVE_IFPACK_TEUCHOS
