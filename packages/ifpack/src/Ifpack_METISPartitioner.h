#ifndef IFPACK_METISPARTITIONER_H
#define IFPACK_METISPARTITIONER_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Teuchos_ParameterList.hpp"
class Ifpack_Graph;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;

//! Ifpack_METISPartitioner: A class to decompose Ifpack_Graph's using METIS.
/*!
Class Ifpack_METISPartitioner enables the decomposition of the local
Ifpack_Graph's using METIS. In order to work properly, this class requires
IFPACK to be configured with option \c --enable-ifpack-metis.
Otherwise, this class will always create one partition.

\date Last modified: Oct-04.
*/

class Ifpack_METISPartitioner : public Ifpack_OverlappingPartitioner {

public:

  //! Constructor.
  Ifpack_METISPartitioner(const Ifpack_Graph* Graph) :
    Ifpack_OverlappingPartitioner(Graph),
    UseSymmetricGraph_(true)
  {}

  //! Destructor.
  virtual ~Ifpack_METISPartitioner() {};

  //! Sets all the parameters for the partitioner (none at moment).
  int SetPartitionParameters(Teuchos::ParameterList& List)
  {
    // `true' is the safest option, as singletons may still be
    // in the system, even after Ifpack_SingletonFilter (think
    // for example to an upper triangular matrix, which has one
    // singleton, and still has one after the elimination of 
    // the first one...)
    UseSymmetricGraph_ = List.get("partitioner: use symmetric graph", 
				  UseSymmetricGraph_);

    return(0);
  }

  //! Computes the partitions. Returns 0 if successful.
  int ComputePartitions();

private:
  bool UseSymmetricGraph_;

}; // class Ifpack_METISPartitioner

#endif // IFPACK_METISPARTITIONER_H
