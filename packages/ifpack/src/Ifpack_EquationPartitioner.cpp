#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_EquationPartitioner.h"
#include "Ifpack_Graph.h"

#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
int Ifpack_EquationPartitioner::ComputePartitions()
{
  
  int mod = NumMyRows() / NumLocalParts_;
  if (mod)
    IFPACK_CHK_ERR(-1); // rows must be multiples of number of equations

  for (int i = 0 ; i < NumMyRows() ; ++i) {
    // Dirichlet nodes get -1
    if (Mask_[i] == -1)
      Partition_[i] = -1;
    else {
      Partition_[i] = i % mod;
    }
  }

  return(0);
}
#endif // HAVE_IFPACK_TEUCHOS
