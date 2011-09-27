
#ifndef _ZOLTAN2_PARTITIONINGPROBLEM_HPP_
#define _ZOLTAN2_PARTITIONINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>

/*! \file Zoltan2_PartitioningProblem.hpp

  This file contains the PartitioningProblem class, which derives from 
  the Problem class.
*/

using Teuchos::RCP;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
CONSISTENT_CLASS_TEMPLATE_LINE
class PartitioningProblem : public Problem<CONSISTENT_TEMPLATE_PARAMS>
{
protected:
  void createPartitioningProblem();

  //TODO RCP<PartitioningSolution> _solution;

public:

  // Destructor
  virtual ~PartitioningProblem() {}

  //! Constructor with Tpetra Matrix interface.
  PartitioningProblem(
    Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &A,
    Teuchos::ParameterList &p
  ) : Problem<CONSISTENT_TEMPLATE_PARAMS>(A, p) 
  {
    HELLO;
    createPartitioningProblem();
  }

  virtual void solve();
};

CONSISTENT_FN_TEMPLATE_LINE
void PartitioningProblem<CONSISTENT_TEMPLATE_PARAMS>::solve()
{
  HELLO;
}

//! createPartitioningProblem 
CONSISTENT_FN_TEMPLATE_LINE
void PartitioningProblem<CONSISTENT_TEMPLATE_PARAMS>::createPartitioningProblem()
{
  HELLO;
  //KDDZoltan2::InputAdapterType adapterType = _inputAdapter->inputAdapterType();

  // Select Model based on parameters and InputAdapter type

}

}
#endif
