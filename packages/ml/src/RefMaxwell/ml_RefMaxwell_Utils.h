#ifndef ML_REFMAXWELL_UTILS_H
#define ML_REFMAXWELL_UTILS_H

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"

namespace ML_Epetra
{

  // Output an Epetra_IntVector
  void IVOUT(const Epetra_IntVector & A, const char *of);
    
  // Getrow that replaces all non-zero entries with ones
  int CSR_getrow_ones(ML_Operator *data, int N_requested_rows, int requested_rows[],
		      int allocated_space, int columns[], double values[], int row_lengths[]);

  // Aggregation
  int RefMaxwell_Aggregate_Nodes(const Epetra_CrsMatrix & A, Teuchos::ParameterList & List, ML_Comm * ml_comm, std::string PrintMsg,
				 ML_Aggregate_Struct *& MLAggr,ML_Operator *&P, int &NumAggregates);
}

    


#endif

#endif
