#ifndef IFPACK_METISPARTITIONER_H
#define IFPACK_METISPARTITIONER_H

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

//! Ifpack_METISPartitioner: A class to decompose overlapping and non-overlapping Ifpack_Graph's.

class Ifpack_METISPartitioner : public Ifpack_OverlappingPartitioner {

public:

  Ifpack_METISPartitioner(const Ifpack_Graph* Graph) :
    Ifpack_OverlappingPartitioner(Graph)
  {}

  ~Ifpack_METISPartitioner() {};

  //! Sets all the parameters for the partitioner.
  int SetPartitionParameters(Teuchos::ParameterList& List)
  {
    return(0);
  }

  //! Computes the partitions. Returns 0 if successful.
  int ComputePartitions();

}; // class Ifpack_METISPartitioner

#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_METISPARTITIONER_H
