// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_PartitioningProblem.hpp

  This file contains the PartitioningProblem class, which derives from 
  the Problem class.
*/

#ifndef _ZOLTAN2_PARTITIONINGPROBLEM_HPP_
#define _ZOLTAN2_PARTITIONINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningAlgorithms.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_IdentifierModel.hpp>

#include <Teuchos_Ptr.hpp>

#ifdef HAVE_OVIS
#include <ovis.h>
#endif

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
template<typename Adapter>
class PartitioningProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::gid_t gid_t;
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
    Problem<Adapter>(A,p), generalParams_(), partitioningParams_(),solution_(),
    inputType_(InvalidAdapterType), modelType_(InvalidModel), algorithm_()
  {
    HELLO;
    createPartitioningProblem();
  };

  // Other methods
  //   LRIESEN - Do we restate virtual in the concrete class?  I
  //    don't think I've seen this style before.
  virtual void solve();

  PartitioningSolution<gid_t, lno_t> &getSolution() {
    return *(solution_.getRawPtr());
  };

private:
  void createPartitioningProblem();

  Teuchos::Ptr<Teuchos::ParameterList> generalParams_;
  Teuchos::Ptr<Teuchos::ParameterList> partitioningParams_;
  RCP<PartitioningSolution<gid_t, lno_t> > solution_;

  InputAdapterType inputType_;
  ModelType modelType_;
  std::string algorithm_;
};

////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void PartitioningProblem<Adapter>::solve()
{
  HELLO;

  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

  size_t nObj = this->generalModel_->getLocalNumObjects();

  size_t numGlobalParts = 
    partitioningParams_->get<int>(string("num_global_parts"));

  // Create the solution.   TODO add exception handling

  solution_ = rcp(new PartitioningSolution<gid_t,lno_t>(numGlobalParts, nObj));

  ArrayRCP<gid_t> &solnGids = solution_->getGidsRCP();
  ArrayRCP<size_t> &solnParts = solution_->getPartsRCP();

  // Call the algorithm

  try {
    if (algorithm_ == string("scotch")){
      AlgPTScotch<base_adapter_t>(this->envConst_, this->comm_, 
        this->graphModel_, numGlobalParts, solnParts(0, nObj));
    }
  }
  Z2_FORWARD_EXCEPTIONS;

  // Write User's GIDs (in the same order in which they appears in the
  //   input adapter's "get" method) to the solution.

  typedef IdentifierMap<gid_t,lno_t,gno_t> idmap_t;
  const RCP<const idmap_t> idMap = this->generalModel_->getIdentifierMap();

  ArrayView<const gno_t> vtxGNO;
  this->generalModel_->getGlobalObjectIds(vtxGNO);

  ArrayRCP<gno_t> gnos = arcpFromArrayView(av_const_cast<gno_t>(vtxGNO));

  if (idMap->gnosAreGids()){
    solnGids = arcp_reinterpret_cast<gid_t>(gnos);
  }
  else{
    idMap->gidTranslate(solnGids(0, nObj), gnos(0, nObj), 
      TRANSLATE_LIB_TO_APP);
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
  using std::string;

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

  generalParams_ = Teuchos::Ptr<Teuchos::ParameterList>(
    &this->env_->getParamsNonConst());

  partitioningParams_ = Teuchos::Ptr<Teuchos::ParameterList>(
    &this->env_->getPartitioningParamsNonConst());

  string paramNotSet("unset");

  // What type of input did the User provide.  If they didn't
  //   specify a model and/or an algorithm, we can use the input
  //   type to choose some defaults.

  inputType_ = this->inputAdapter_->inputAdapterType();

  // Did user specify a computational model?

  string &model = partitioningParams_->get<string>(string("model"), 
    paramNotSet);

  // What type of partitioning algorithm, if any, is the user asking for?

  string &algorithm = partitioningParams_->get<string>(string("algorithm"), 
    paramNotSet);

  // Does the algorithm require consecutive global IDs?

  bool needConsecutiveGlobalIds = false;

  // Do matrix and graph algorithms require that diagonal entries
  //    or self-edges be removed?

  bool removeSelfEdges= false;

  // Determine algorithm, model, and algorithm requirements.  This
  // is a first pass.  Feel free to change this and add to it.
  
  if (algorithm != paramNotSet){
    // Figure out the model required by the algorithm
    if (algorithm == string("block") ||
        algorithm == string("random") ||
        algorithm == string("cyclic") ){

      modelType_ = IdentifierModelType;
      algorithm_ = algorithm;
      needConsecutiveGlobalIds = true;
    }
    else if (algorithm == string("rcb") ||
             algorithm == string("rib") ||
             algorithm == string("hsfc")){

      modelType_ = GeometryModelType;
      algorithm_ = algorithm;
      needConsecutiveGlobalIds = true;
    }
    else if (algorithm == string("metis") ||
             algorithm == string("parmetis") ||
             algorithm == string("scotch") ||
             algorithm == string("ptscotch")){

      modelType_ = GraphModelType;
      algorithm_ = algorithm;
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
    }
    else if (algorithm == string("patoh") ||
             algorithm == string("phg")){

      if ((modelType_ != GraphModelType) &&
          (modelType_ != HypergraphModelType) ){
        modelType_ = HypergraphModelType;
      }
      algorithm_ = algorithm;
      needConsecutiveGlobalIds = true;
    }
    else{
      // Parameter list should ensure this does not happen.
      throw std::logic_error("parameter list algorithm is invalid");
    }
  }
  else if (model != paramNotSet){
    // Figure out the algorithm suggested by the model.
    if (model == string("hypergraph")){
      modelType_ = HypergraphModelType;
      if (this->comm_->getSize() > 1)
        algorithm_ = string("phg"); 
      else
        algorithm_ = string("patoh"); 
      needConsecutiveGlobalIds = true;
    }
    else if (model == string("graph")){
      modelType_ = GraphModelType;
#ifdef HAVE_SCOTCH
      if (this->comm_->getSize() > 1)
        algorithm_ = string("ptscotch"); 
      else
        algorithm_ = string("scotch"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#else
#ifdef HAVE_PARMETIS
      if (this->comm_->getSize() > 1)
        algorithm_ = string("parmetis"); 
      else
        algorithm_ = string("metis"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#else
      if (this->comm_->getSize() > 1)
        algorithm_ = string("phg"); 
      else
        algorithm_ = string("patoh"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#endif
#endif
    }
    else if (model == string("geometry")){
      modelType_ = GeometryModelType;
      algorithm_ = string("rib");
      needConsecutiveGlobalIds = true;
    }
    else if (model == string("ids")){
      modelType_ = IdentifierModelType;
      algorithm_ = string("block");
      needConsecutiveGlobalIds = true;
    }
    else{
      // Parameter list should ensure this does not happen.
      throw std::logic_error("parameter list model type is invalid");
    }
  }
  else{   
    // Determine an algorithm and model suggested by the input type.
    //   TODO: this is a good time to use the time vs. quality parameter
    //     in choosing an algorithm, and setting some parameters

    if (inputType_ == MatrixAdapterType){
      modelType_ = HypergraphModelType;
      if (this->comm_->getSize() > 1)
        algorithm_ = string("phg"); 
      else
        algorithm_ = string("patoh"); 
    }
    else if (inputType_ == GraphAdapterType ||
        inputType_ == MeshAdapterType){
      modelType_ = GraphModelType;
      if (this->comm_->getSize() > 1)
        algorithm_ = string("phg"); 
      else
        algorithm_ = string("patoh"); 
    }
    else if (inputType_ == CoordinateAdapterType){
      modelType_ = GeometryModelType;
      algorithm_ = string("rib");
    }
    else if (inputType_ == VectorAdapterType ||
             inputType_ == MultiVectorAdapterType ||
             inputType_ == IdentifierAdapterType){
      modelType_ = IdentifierModelType;
      algorithm_ = string("block");
    }
    else{
      // This should never happen
      throw std::logic_error("input type is invalid");
    }
  }

  // TODO: This doesn't work.  baseInputAdapter_.getRawPtr() is NULL.
  //
  // RCP<const base_adapter_t> baseInputAdapter_ = 
  //   rcp_implicit_cast<const base_adapter_t>(this->inputAdapter_);
  //
  // So to pass the InputAdapter to the Model we use a raw pointer.
  // Since the Problem creates the Model and will destroy it when
  // done, the Model doesn't really need an RCP to the InputAdapter.
  // But it would be nice if that worked.

  // TODO - check for exceptions

  typedef typename Adapter::base_adapter_t base_adapter_t;
  const base_adapter_t *baseAdapter = this->inputAdapter_.getRawPtr();

  // Create the computational model.

  switch (modelType_) {

  case GraphModelType:
    this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
      baseAdapter, this->envConst_, needConsecutiveGlobalIds, removeSelfEdges));

    this->generalModel_ = rcp_implicit_cast<Model<base_adapter_t> >(
      this->graphModel_);

    break;

  case HypergraphModelType:
    break;

  case GeometryModelType:
    break;

  case IdentifierModelType:
    this->identifierModel_ = rcp(new IdentifierModel<base_adapter_t>(
      baseAdapter, this->envConst_, needConsecutiveGlobalIds));

    this->generalModel_ = rcp_implicit_cast<Model<base_adapter_t> >(
      this->identifierModel_);
    break;

  default:
    cout << __func__ << " Invalid model" << modelType_ << endl;
    break;
  }
}

}
#endif
