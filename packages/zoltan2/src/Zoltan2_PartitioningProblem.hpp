
#ifndef _ZOLTAN2_PARTITIONINGPROBLEM_HPP_
#define _ZOLTAN2_PARTITIONINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningAlgorithms.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <Zoltan2_GraphModel.hpp>

/*! \file Zoltan2_PartitioningProblem.hpp

  This file contains the PartitioningProblem class, which derives from 
  the Problem class.
*/

using Teuchos::rcp_dynamic_cast;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
template<typename Adapter>
class PartitioningProblem : public Problem<Adapter>
{
protected:
  void createPartitioningProblem();

  RCP<PartitioningSolution<Adapter> > solution_;

public:

  // Destructor
  virtual ~PartitioningProblem() {};

#if 0  // KDDKDD Don't know how to use shortcut with Adapter template
  //! Constructor with Tpetra Matrix interface.
  PartitioningProblem(Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &A,
    ParameterList &p
  ) : Problem<Adapter>(A, p) 
  {
    HELLO;
    createPartitioningProblem();
  }
#endif

  //! Constructor with InputAdapter Interface
  PartitioningProblem(Adapter *A, Teuchos::ParameterList *p) 
                      : Problem<Adapter>(A, p) 
  {
    HELLO;
    createPartitioningProblem();
  };

  // Other methods
  virtual void solve();
  virtual void redistribute();

  PartitioningSolution<Adapter> *getSolution() {
    return solution_.getRawPtr();
  };
};

////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void PartitioningProblem<Adapter>::solve()
{
  HELLO;

  this->solution_ = rcp(new PartitioningSolution<Adapter>);

  // Determine which algorithm to use based on defaults and parameters.
  // For now, assuming Scotch graph partitioning.
  // Need some exception handling here, too.

  AlgPTScotch<Adapter>(this->graphModel_, this->solution_, this->params_,
                       this->comm_);
}

////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void PartitioningProblem<Adapter>::redistribute()
{
  HELLO;
}

////////////////////////////////////////////////////////////////////////
//! createPartitioningProblem 
//  Method with common functionality for creating a PartitioningProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <typename Adapter>
void PartitioningProblem<Adapter>::createPartitioningProblem()
{
  HELLO;
  cout << __func__ << " input adapter type " 
       << this->inputAdapter_->inputAdapterType() << " " 
       << this->inputAdapter_->inputAdapterName() << endl;

  // Determine which parameters are relevant here.
  // For now, assume parameters similar to Zoltan:
  //   MODEL = graph, hypergraph, geometric, ids
  //   APPROACH = partition, repartition
  //   ALGORITHM = metis, parmetis, scotch, ptscotch, patoh, 
  //               phg, rcb, rib, hsfc, block, cyclic, random
  // TODO: I will need help from Lee Ann understanding how to use the parameter
  // functionality in Zoltan2.  For now, I will set a few parameters and
  // continue computing.
  ModelType modelType = GraphModelType;

  // Select Model based on parameters and InputAdapter type
  switch (modelType) {

  case GraphModelType:

    cout << __func__ << "Xpetra matrix adapter switch" << endl;

    this->graphModel_ = RCP<GraphModel<Adapter> > 
                        (new GraphModel<Adapter>(this->inputAdapter_,
                                                 this->comm_, this->env_));
    break;

  case HypergraphModelType:
  case GeometryModelType:
  case IdModelType:
    cout << __func__ << " Model type " << modelType << " not yet supported." 
         << endl;
    break;

  default:
    cout << __func__ << " Invalid model" << modelType << endl;
    break;
  }
}

}
#endif
