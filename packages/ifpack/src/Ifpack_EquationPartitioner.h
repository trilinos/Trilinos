#ifndef IFPACK_EQUATIONPARTITIONER_H
#define IFPACK_EQUATIONPARTITIONER_H

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

//! Ifpack_EquationPartitioner: A class to decompose an Ifpack_Graph so that each block will contain all the rows for a different equation.

/*!
Ifpack_EquationPartitioner enables a decomposition into blocks of equations.
Suppose that the input Ifpack_Graph is based on an Epetra_RowMatrix, whose
rows represent (U_i,V_i,P_i) for each grid node i. This partitioner
will decompose the graph into three subgraphs, each of them containing
the rows of U, then V, than P.

The number of equations is set as the number of local partitions.

\note It is required that NumRows % NumLocalParts() = 0.

\date Sep-04.
*/
class Ifpack_EquationPartitioner : public Ifpack_OverlappingPartitioner {

public:

  //! Constructor.
  Ifpack_EquationPartitioner(const Ifpack_Graph* Graph) :
    Ifpack_OverlappingPartitioner(Graph)
  {}

  //! Destructor.
  ~Ifpack_EquationPartitioner() {};

  //! Sets all the parameters for the partitioner.
  int SetPartitionParameters(Teuchos::ParameterList& List)
  {
    return(0);
  }

  //! Computes the partitions. Returns 0 if successful.
  int ComputePartitions();

private:

};

#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_EQUATIONPARTITIONER_H
