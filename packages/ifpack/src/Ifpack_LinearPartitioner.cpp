#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_LinearPartitioner.h"

//==============================================================================
int Ifpack_LinearPartitioner::ComputePartitions()
{
  
  int mod = NumMyRows() / NumLocalParts_;
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    Partition_[i] = i / mod;
    if (Partition_[i] >= NumLocalParts_)
      Partition_[i] = NumLocalParts_ - 1;
  }

  return(0);
}
#endif // HAVE_IFPACK_TEUCHOS
