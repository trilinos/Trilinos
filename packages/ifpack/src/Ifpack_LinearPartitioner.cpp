#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_LinearPartitioner.h"
#include "Ifpack_Graph.h"

#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
int Ifpack_LinearPartitioner::ComputePartitions()
{
  
  int mod = NumMyRows() / NumLocalParts_;
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    // Dirichlet nodes get -1
    if (Mask_[i] == -1)
      Partition_[i] = -1;
    else {
      Partition_[i] = i / mod;
      if (Partition_[i] >= NumLocalParts_)
	Partition_[i] = NumLocalParts_ - 1;
    }
  }

  return(0);
}
#endif // HAVE_IFPACK_TEUCHOS
