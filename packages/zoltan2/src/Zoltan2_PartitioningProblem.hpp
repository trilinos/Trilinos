
#ifndef _ZOLTAN2_PARTITIONINGPROBLEM_HPP_
#define _ZOLTAN2_PARTITIONINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningAlgorithms.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

/*! \file Zoltan2_PartitioningProblem.hpp

  This file contains the PartitioningProblem class, which derives from 
  the Problem class.
*/

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
template<Z2CLASS_TEMPLATE>
class PartitioningProblem : public Problem<Z2PARAM_TEMPLATE>
{
protected:
  void createPartitioningProblem();

  RCP<PartitioningSolution<Z2PARAM_TEMPLATE> > solution_;

public:

  // Destructor
  virtual ~PartitioningProblem() {}

  //! Constructor with Tpetra Matrix interface.
  PartitioningProblem(
    Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &A,
    ParameterList &p
  ) : Problem<Z2PARAM_TEMPLATE>(A, p) 
  {
    HELLO;
    createPartitioningProblem();
  }

  // Other methods
  virtual void solve();
  virtual void redistribute();

};

////////////////////////////////////////////////////////////////////////
template <Z2FN_TEMPLATE>
void PartitioningProblem<Z2PARAM_TEMPLATE>::solve()
{
  HELLO;
  // Determine which algorithm to use based on defaults and parameters.
  // For now, assuming Scotch graph partitioning.
  // Need some exception handling here, too.

  AlgScotch<Z2PARAM_TEMPLATE> alg(this->model_, this->solution_, this->params_);
}

////////////////////////////////////////////////////////////////////////
template <Z2FN_TEMPLATE>
void PartitioningProblem<Z2PARAM_TEMPLATE>::redistribute()
{
  HELLO;
}

////////////////////////////////////////////////////////////////////////
//! createPartitioningProblem 
//  Method with common functionality for creating a PartitioningProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <Z2FN_TEMPLATE>
void PartitioningProblem<Z2PARAM_TEMPLATE>::createPartitioningProblem()
{
  HELLO;
  cout << __func__ << " input adapter type " 
       << this->inputAdapter_->inputAdapterType() << " " 
       << this->inputAdapter_->inputAdapterName() << endl;

  InputAdapterType adapterType = this->inputAdapter_->inputAdapterType();

  // Determine which parameters are relevant here.
  // For now, assume parameters similar to Zoltan:
  //   MODEL = graph, hypergraph, geometric, ids
  //   APPROACH = partition, repartition
  //   ALGORITHM = metis, parmetis, scotch, ptscotch, patoh, 
  //               phg, rcb, rib, hsfc, block, cyclic, random
  // TODO: I will need help from Lee Ann understanding how to use the parameter
  // functionality in Zoltan2.  For now, I will set a few parameters and
  // continue computing.
  ModelType model = GraphModelType;

  // Select Model based on parameters and InputAdapter type
  switch (model) {
  case GraphModelType:
    switch (adapterType) {
    case MatrixAdapterType:
      cout << __func__ << " Matrix adapter switch" << endl;
      //TODO model_ = 
      //TODO    new GraphModel<MatrixInput, Z2PARAM_TEMPLATE>(this->inputAdapter_);
      break;
    case GraphAdapterType:
      cout << __func__ << " Graph adapter switch" << endl;
      //TODO model_ = 
      //TODO    new GraphModel<GraphInput, Z2PARAM_TEMPLATE>(this->inputAdapter_);
      break;
    case MeshAdapterType:
    case CoordAdapterType:
    case IdAdapterType:
      cout << __func__ 
           << " PartitioningProblem not yet implemented for this input adapter "
           << this->inputAdapter_->inputAdapterName() << endl;
      break;
    default:
      cout << "Invalid adapter type; this condition should never happen." 
           << endl;
      break;
    }
    break;
  case HypergraphModelType:
  case GeometryModelType:
  case IdModelType:
    cout << __func__ << " Model type " << model << " not yet supported." 
         << endl;
    break;
  default:
    cout << __func__ << " Invalid model" << model << endl;
    break;
  }
}

}
#endif
