#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_GreedyPartitioner.h"
#include "Ifpack_Graph.h"

#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
int Ifpack_GreedyPartitioner::ComputePartitions()
{
  vector<int> ElementsPerPart;
  ElementsPerPart.resize(NumLocalParts());

  vector<int> count;
  count.resize(NumLocalParts());

  // define how many nodes have to be put on each part

  int div = NumMyRows() / NumLocalParts();
  int mod = NumMyRows() % NumLocalParts();

  for (int i = 0 ; i < NumLocalParts() ; ++i) {
    count[i] = 0;
    ElementsPerPart[i] = div;
    if (i < mod) ElementsPerPart[i]++;
  }

  for( int i=0 ; i<NumMyRows() ; ++i ) {
    Partition_[i] = -1;
  }

  int CrsNumEntries;

  int CurrentPart = 0;
  vector<int> Indices;
  Indices.resize(MaxNumEntries());
  
  bool ok = true;

  // pick up a root node. Check that this is not a Dirichlet node.

  int RootNode = RootNode_;

  // start from row 0, assigned to domain 0
  Partition_[RootNode] = 0;      

  // cycle over all non-Dirichlet nodes to perform the aggregation.

  while (ok == true) {

    int ierr = Graph_->ExtractMyRowCopy(RootNode, MaxNumEntries(),
				CrsNumEntries, &Indices[0]);

    IFPACK_CHK_ERR(ierr);

    ok = false;

    for (int j = 0 ; j < CrsNumEntries ; ++j) {

      // go to the next indices if Dirichlet node
      // or off-processor one.

      if (Indices[j] >= NumMyRows()) 
	continue;

      // filter for Dirichlet
      int col = Indices[j];

      if (count[CurrentPart] == ElementsPerPart[CurrentPart]) {
	CurrentPart++;
      }

      if (Partition_[col] == -1) {
	Partition_[col] = CurrentPart;
	if (ok == false) {
	  ok = true;
	  RootNode = col;
	}
	count[CurrentPart]++;
      }
    }

    // if ok == false, it means that no neighboring node
    // was included in the previous part. This may signify that all
    // neighboring nodes are already assigned, or that the graph is
    // disconnected. In both cases, we look for a new root node, starting from
    // 0. 

    if (ok == false) {
      for (int j = 0 ; j < NumMyRows() ; ++j) {

	if (Partition_[j] == -1 ) {
	  RootNode = j;
	  ok = true;
	  break;
	}
      }
    }
  }

  return(0);
}

#endif // HAVE_IFPACK_TEUCHOS
