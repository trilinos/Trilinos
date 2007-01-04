#ifndef IFPACK_USERPARTITIONER_H
#define IFPACK_USERPARTITIONER_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Teuchos_ParameterList.hpp"
class Epetra_Comm;
class Ifpack_Graph;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;

//! Ifpack_UserPartitioner: A class to define linear partitions.

class Ifpack_UserPartitioner : public Ifpack_OverlappingPartitioner {

public:

  //! Constructor.
  Ifpack_UserPartitioner(const Ifpack_Graph* Graph) :
    Ifpack_OverlappingPartitioner(Graph),
    Map_(0)
  {}

  //! Destructor.
  virtual ~Ifpack_UserPartitioner() {};

  //! Sets all the parameters for the partitioner (none for linear partioning).
  int SetPartitionParameters(Teuchos::ParameterList& List)
  {
    Map_ = List.get("partitioner: map",Map_);
    if (Map_ == 0)
      IFPACK_CHK_ERR(-1);

    return(0);
  }

  //! Computes the partitions. Returns 0 if successful.
  int ComputePartitions();

private:
  int* Map_;
};

#endif // IFPACK_USERPARTITIONER_H
