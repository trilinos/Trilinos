#ifndef IFPACK_GREEDYPARTITIONER_H
#define IFPACK_GREEDYPARTITIONER_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Teuchos_ParameterList.hpp"
#include <vector>
class Epetra_Comm;
class Ifpack_Graph;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;

//! Ifpack_GreedyPartitioner: A class to decompose Ifpack_Graph's using a simple greedy algorithm.

class Ifpack_GreedyPartitioner : public Ifpack_OverlappingPartitioner {

public:

  //! Constructor.
  Ifpack_GreedyPartitioner(const Ifpack_Graph* Graph) :
    Ifpack_OverlappingPartitioner(Graph),
    RootNode_(0)
  {}

  //! Destructor.
  ~Ifpack_GreedyPartitioner() {};

  //! Sets all the parameters for the partitioner (root node).
  int SetPartitionParameters(Teuchos::ParameterList& List)
  {
    RootNode_ = List.get("partitioner: root node", RootNode_);

    return(0);
  }

  //! Computes the partitions. Returns 0 if successful.
  int ComputePartitions();

private:

  int RootNode_;

}; // class Ifpack_GreedyPartitioner

#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_GREEDYPARTITIONER_H
