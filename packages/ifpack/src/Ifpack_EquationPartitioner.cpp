#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_EquationPartitioner.h"
#include "Epetra_CrsGraph.h"

//==============================================================================
int Ifpack_EquationPartitioner::ComputePartitions()
{
  
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    Partition_[i] = i % NumLocalParts_;
  }

  return(0);
}
#endif // HAVE_IFPACK_TEUCHOS
