#ifndef ML_REFMAXWELL_UTILS_H
#define ML_REFMAXWELL_UTILS_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_Array.hpp"

namespace ML_Epetra
{

  class CoordPack {
  public:
    Teuchos::Array<double> x, y, z, material;
  };

  // Output an Epetra_IntVector
  void IVOUT(const Epetra_IntVector & A, const char *of);

  // Getrow that replaces all non-zero entries with ones
  int CSR_getrow_ones(ML_Operator *data, int N_requested_rows, int requested_rows[],
		      int allocated_space, int columns[], double values[], int row_lengths[]);

  // Migrate Coordinates from DomainMap to ColumnMap
  int RefMaxwell_SetupCoordinates(ML_Operator* A, 
                                  const double *in_coordx, const double *in_coordy, const double *in_coordz, const double* in_material,
                                  double *&coordx, double *&coordy, double *&coordz, double*& material);

  // Aggregation
  int RefMaxwell_Aggregate_Nodes(const Epetra_CrsMatrix & A, Teuchos::ParameterList & List, ML_Comm * ml_comm, std::string PrintMsg,
				 ML_Aggregate_Struct *& MLAggr,ML_Operator *&P, int &NumAggregates, CoordPack & pack);

  // Edge Nullspace
  Epetra_MultiVector* Build_Edge_Nullspace(const Epetra_CrsMatrix & D0Clean_Matrix, const Teuchos::ArrayRCP<int> BCedges, Teuchos::ParameterList & List, bool verbose);

  // Pi Matrix
  Epetra_CrsMatrix *Build_Pi_Matrix(const Epetra_CrsMatrix & D0Clean_Matrix, const Teuchos::ArrayRCP<int> BCedges, Teuchos::ParameterList & List, bool verbose);
}





#endif

#endif
