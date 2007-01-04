#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_UserPartitioner.h"
#include "Epetra_CrsGraph.h"
#include <vector>

//==============================================================================
int Ifpack_UserPartitioner::ComputePartitions()
{
  
  if (Map_ == 0)
    IFPACK_CHK_ERR(-1);

  // simply copy user's vector
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    Partition_[i] = Map_[i];
  }

  // put together all partitions composed by 1 one vertex
  // (if any)
  vector<int> singletons(NumLocalParts());
  for (unsigned int i = 0 ; i < singletons.size() ; ++i) {
    singletons[i] = 0;
  }

#if 0
  // may want to uncomment the following to ensure that no
  // partitions are in fact singletons
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    ++singletons[Partition_[i]];
  }
  
  int count = 0;
  for (unsigned int i = 0 ; i < singletons.size() ; ++i) {
    if (singletons[i] == 1)
      ++count;
  }

  int index = -1;
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    int j = Partition_[i];
    if (singletons[j] == 1) {
      if (index == -1)
        index = j;
      else
        Partition_[i] = index;
    }
  }
#endif

  return(0);
}
