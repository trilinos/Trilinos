
#ifndef _ZOLTAN2_PARTITIONINGPROBLEM_HPP_
#define _ZOLTAN2_PARTITIONINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
// #include <Zoltan2_PartitioningAlgorithms.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <Zoltan2_GraphModel.hpp>

/*! \file Zoltan2_PartitioningProblem.hpp

  This file contains the PartitioningProblem class, which derives from 
  the Problem class.
*/

using Teuchos::rcp_dynamic_cast;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
template<typename User>
class PartitioningProblem : public Problem<User>
{
protected:
  void createPartitioningProblem();

  RCP<PartitioningSolution<User> > solution_;

public:

  // Destructor
  virtual ~PartitioningProblem() {};

#if 0  // KDDKDD Don't know how to use shortcut with User template
  //! Constructor with Tpetra Matrix interface.
  PartitioningProblem(Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> &A,
    ParameterList &p
  ) : Problem<User>(A, p) 
  {
    HELLO;
    createPartitioningProblem();
  }
#endif

  //! Constructor with InputAdapter Interface
  PartitioningProblem(InputAdapter<User> *A, Teuchos::ParameterList *p) 
                      : Problem<User>(A, p) 
  {
    HELLO;
    createPartitioningProblem();
  };

  // Other methods
  virtual void solve();
  virtual void redistribute();

};

////////////////////////////////////////////////////////////////////////
template <typename User>
void PartitioningProblem<User>::solve()
{
  HELLO;
  // Determine which algorithm to use based on defaults and parameters.
  // For now, assuming Scotch graph partitioning.
  // Need some exception handling here, too.

  // AlgScotch<User> alg(this->model_, this->solution_, this->params_);
}

////////////////////////////////////////////////////////////////////////
template <typename User>
void PartitioningProblem<User>::redistribute()
{
  HELLO;
}

////////////////////////////////////////////////////////////////////////
//! createPartitioningProblem 
//  Method with common functionality for creating a PartitioningProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <typename User>
void PartitioningProblem<User>::createPartitioningProblem()
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
  ModelType modelType = GraphModelType;

  // Select Model based on parameters and InputAdapter type
  switch (modelType) {

  case GraphModelType:

    switch (adapterType) {

    case XpetraCrsMatrixAdapterType: {
      cout << __func__ << "Xpetra matrix adapter switch" << endl;

      GraphModel<XpetraCrsMatrixInput<User> > *model =
         new GraphModel<XpetraCrsMatrixInput<User> >(
                        rcp_dynamic_cast<const XpetraCrsMatrixInput<User> >
                                        (this->inputAdapter_, true),
                        this->comm_,
                        this->env_);
      delete model;  //KDDKDD REMOVE LATER WHEN HAVE RCPs WORKING.

      //KDDKDD NOT WORKING YET RCP<GraphModel<XpetraCrsMatrixInput<User> > > rcpmodel = rcp(model);

      //KDDKDD NOT WORKING YET this->model_ = rcp_dynamic_cast<GraphModel<InputAdapter<User> > >(rcpmodel, true);
                       
      // KDDKDD Question:  Not clear on use of RCPs here; compiler complains
      // KDDKDD that it cannot do conversions when I use GraphModel and 
      // KDDKDD XpetraCrsMatrixInput.
      break;
    }

    case MatrixAdapterType: {
      cout << __func__ << " Matrix adapter switch" << endl;
      GraphModel<MatrixInput<User> > *model = 
           new GraphModel<MatrixInput<User> >(//KDD No GraphModel 
                                              //KDD constructor yet.
                                              //KDD this->inputAdapter_,
                                              //KDD this->comm_,
                                              //KDD this->env_
                                              );
      //KDDKDD NOT WORKING YET this->model_ = RCP<GraphModel<MatrixInput<User> > >(model);
      break;
    }

    case GraphAdapterType:
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
