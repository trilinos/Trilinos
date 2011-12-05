
#ifndef _ZOLTAN2_PARTITIONINGPROBLEM_HPP_
#define _ZOLTAN2_PARTITIONINGPROBLEM_HPP_

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Ptr.hpp>

#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningAlgorithms.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <Zoltan2_GraphModel.hpp>
#ifdef HAVE_OVIS
#include <ovis.h>
#endif

/*! \file Zoltan2_PartitioningProblem.hpp

  This file contains the PartitioningProblem class, which derives from 
  the Problem class.
*/

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
template<typename Adapter>
class PartitioningProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lid_t lid_t;
  typedef typename Adapter::lno_t lno_t;

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
  PartitioningProblem(Adapter *A, Teuchos::ParameterList *p): 
    Problem<Adapter>(A,p), generalParams_(), partitioningParams_(),solution_()
  {
    HELLO;
    createPartitioningProblem();
  };

  // Other methods
  //   LRIESEN - Do we restate virtual in the concrete class?  I
  //    don't think I've seen this style before.
  virtual void solve();

  PartitioningSolution<gid_t, lid_t, lno_t> &getSolution() {
    return *(solution_.getRawPtr());
  };

private:
  void createPartitioningProblem();

  Teuchos::Ptr<const Teuchos::ParameterList> generalParams_;
  Teuchos::Ptr<const Teuchos::ParameterList> partitioningParams_;
  RCP<PartitioningSolution<gid_t, lid_t, lno_t> > solution_;

};

////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void PartitioningProblem<Adapter>::solve()
{
  HELLO;

  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lid_t lid_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

  size_t nVtx = this->graphModel_->getLocalNumVertices();

  std::string algorithm = partitioningParams_->get<std::string>(
   string("algorithm"));

  size_t numGlobalParts = partitioningParams_->get<int>(
    string("num_global_parts"));

  size_t *tmp = new size_t [nVtx];
  Z2_GLOBAL_MEMORY_ASSERTION(*(this->env_), nVtx, !nVtx||tmp);

  ArrayRCP<size_t> parts(tmp, 0, nVtx, true);

  try {
    // Determine which algorithm to use based on defaults and parameters.
    // For now, assuming Scotch graph partitioning.

    AlgPTScotch<base_adapter_t>(this->envConst_, this->comm_, 
      this->graphModel_, numGlobalParts, parts.persistingView(0, numGlobalParts));
  }
  Z2_FORWARD_EXCEPTIONS;

  // Create the solution.   TODO add exception handling

  this->solution_ = 
    rcp(new PartitioningSolution<gid_t,lid_t,lno_t>( numGlobalParts, nVtx,
      this->inputAdapter_->haveLocalIds() ? nVtx : 0));

  this->solution_->getPartsRCP() = parts;

  typedef IdentifierMap<lid_t,gid_t,lno_t,gno_t> idmap_t;
  const RCP<const idmap_t> idMap = this->graphModel_->getIdentifierMap();

  ArrayView<const gno_t> vtxID;
  ArrayView<const scalar_t> xyz, vtxWt;
  this->graphModel_->getVertexList(vtxID, xyz, vtxWt);
  ArrayRCP<gno_t> gnos = arcpFromArrayView(av_const_cast<gno_t>(vtxID));

  if (idMap->gnosAreGids()){
    this->solution_->getGidsRCP() = gnos;   // TODO may need a cast
  }
  else{
    ArrayRCP<gid_t> gids(nVtx);
    idMap->gidTranslate(gids, gnos, TRANSLATE_LIB_TO_APP);
    this->solution_->getGidsRCP() = gids;
  }

  if (this->inputAdapter_->haveLocalIds()){
    ArrayRCP<lid_t> lids(nVtx);
    idMap->lidTranslate(lids, gnos, TRANSLATE_LIB_TO_APP);
    this->solution_->getLidsRCP() = lids;
  }
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

#ifdef HAVE_OVIS
  ovis_enabled(this->comm_->getRank());
#endif

  // The problem communicator is this->comm_.  
  // Set the application communicator to MPI_COMM_WORLD.

  this->env_->setCommunicator(DefaultComm<int>::getComm());

  // Finalize parameters.  If the Problem wants to set or
  // change any parameters, do it before this call.

  this->env_->commitParameters();

  // Get the parameters

  generalParams_ = Teuchos::Ptr<const Teuchos::ParameterList>(
    &this->env_->getParams());

  partitioningParams_ = Teuchos::Ptr<const Teuchos::ParameterList>(
    &this->env_->getPartitioningParams());

  // Determine which parameters are relevant here.
  // For now, assume parameters similar to Zoltan:
  //   MODEL = graph, hypergraph, geometric, ids
  //   APPROACH = partition, repartition
  //   ALGORITHM = metis, parmetis, scotch, ptscotch, patoh, 
  //               phg, rcb, rib, hsfc, block, cyclic, random
  // TODO: I will need help from Lee Ann understanding how to use the parameter
  // functionality in Zoltan2.  For now, I will set a few parameters and
  // continue computing.

  ModelType modelType = GraphModelType;  // make this a class variable
  typedef typename Adapter::base_adapter_t base_adapter_t;

  // TODO: This doesn't work.  baseInputAdapter_.getRawPtr() is NULL.
  //
  // RCP<const base_adapter_t> baseInputAdapter_ = 
  //   rcp_implicit_cast<const base_adapter_t>(this->inputAdapter_);
  //
  // So to pass the InputAdapter to the Model we use a raw pointer.
  // Since the Problem creates the Model and will destroy it when
  // done, the Model doesn't really need an RCP to the InputAdapter.
  // But it would be nice if that worked.

  const base_adapter_t *baseAdapter = this->inputAdapter_.getRawPtr();

  // Select Model based on parameters and InputAdapter type
  switch (modelType) {

  case GraphModelType:
    this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
      baseAdapter, this->env_, false, true));
    break;

  case HypergraphModelType:
  case GeometryModelType:
  case IdentifierModelType:
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
