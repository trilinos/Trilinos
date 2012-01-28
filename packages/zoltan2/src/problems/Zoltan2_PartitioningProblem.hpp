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

#ifdef HAVE_ZOLTAN2_OVIS
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
  typedef typename Adapter::user_t user_t;

  // Destructor
  virtual ~PartitioningProblem() {};

#ifdef HAVE_ZOLTAN2_MPI
  //! Constructor for MPI builds
  PartitioningProblem(Adapter *A, Teuchos::ParameterList *p, 
    MPI_Comm comm=MPI_COMM_WORLD); 
#else
  //! Constructor for serial builds
  PartitioningProblem(Adapter *A, Teuchos::ParameterList *p) ;
#endif

  // Other methods
  void solve();

  PartitioningSolution<user_t> &getSolution() {
    return *(solution_.getRawPtr());
  };

  /*! Set relative sizes for the partitions that Zoltan2 will create.
   *
   *  \param len  The size of the partIds and partSizes lists
   *  \param partIds   A list of len partition identifiers.  Partition
   *           identifiers range from zero to one less than the global
   *           number of identifiers.  
   *  \param partSizes  A list of len relative sizes corresponding to
   *           the partIds.
   *  \param makeCopy  If true, Zoltan2 will make a copy of the ids and sizes
   *      that are provided in this call.  If false, Zoltan2 will just save
   *      the pointers to to the caller's lists.  If the pointers will remain
   *      remain valid throughout the lifetime of the PartitioningProblem,
   *      and memory use is an issue, then set makeCopy to false.  By default,
   *      Zoltan2 will copy the caller's list of ids and sizes.
   *
   * A given partid should only be provided once across the application.
   * Duplicate partIds will generate a std::runtime_error exception when
   * the PartitioningSolution is created.  Part
   * ids that are omitted will be assigned the average of the sizes that
   * have been specified.
   *
   * Subsequent calls to setPartSizes will replace the list of part ids
   * and part sizes provided previously.
   * 
   * If the application has set multiple weights per object, then the
   * part sizes supplied in this method are applied to the first weight.
   *
   * Zoltan2 assumes that uniform partition sizes are desired by the caller,
   * unless specified otherwise in a call to setPartSizes or 
   * setPartSizesForCritiera.
   */

  void setPartSizes(int len, size_t *partIds, float *partSizes, 
    bool makeCopy=true) 
  { 
    setPartSizesForCritiera(0, len, partIds, partSizes, makeCopy);
  }

  /*! Set the relative sizes (per weight) for the partitions that Zoltan2 will 
   *    create.
   *
   *  \param criteria the criteria (weight dimension) for which these 
   *      part sizes apply.  Criteria range from zero to one less than
   *     the number of weights per object specified in the 
   *     caller's InputAdapter.
   *  \param len  The size of the partIds and partSizes lists
   *  \param partIds   A list of len partition identifiers.  Partition
   *           identifiers range from zero to one less than the global
   *           number of identifiers.  
   *  \param partSizes  A list of len relative sizes corresponding to
   *           the partIds.
   *  \param makeCopy  If true, Zoltan2 will make a copy of the ids and sizes
   *      that are provided in this call.  If false, Zoltan2 will just save
   *      the pointers to to the caller's lists.  If the pointers will remain
   *      remain valid throughout the lifetime of the PartitioningProblem,
   *      and memory use is an issue, then set makeCopy to false.  By default,
   *      Zoltan2 will copy the caller's list of ids and sizes.
   *
   * A given partid should only be provided once across the application.
   * Duplicate partIds will generate a std::runtime_error exception when
   * the PartitioningSolution is created.  Part
   * ids that are omitted will be assigned the average of the sizes that
   * have been specified.
   *
   * Subsequent calls to setPartSizes for the same criteria will replace 
   * the list of part ids and part sizes provided for that criteria previously.
   *
   * Zoltan2 assumes that uniform partition sizes are desired by the caller,
   * unless specified otherwise in a call to setPartSizes or 
   * setPartSizesForCritiera.
   */

  void setPartSizesForCritiera(int criteria, int len, size_t *partIds, 
    float *partSizes, bool makeCopy=true) ;

private:
  void createPartitioningProblem();

  Teuchos::Ptr<Teuchos::ParameterList> generalParams_;
  Teuchos::Ptr<Teuchos::ParameterList> partitioningParams_;
  RCP<PartitioningSolution<user_t> > solution_;

  InputAdapterType inputType_;
  ModelType modelType_;
  std::string algorithm_;

  int numberOfWeights_;

  // Suppose Array<size_t> partIds = partIds_[w].  If partIds.size() > 0
  // then the user supplied part sizes for weight index "w", and the sizes
  // corresponding to the Ids in partIds are partSizes[w].
  //
  // If numberOfWeights_ >= 0, then there is an Id and Sizes array for
  // for each weight.  Otherwise the user did not supply object weights,
  // but they can still specify part sizes. 
  // So numberOfCriteria_ is numberOfWeights_ or one, whichever is greater.

  ArrayRCP<ArrayRCP<size_t> > partIds_;
  ArrayRCP<ArrayRCP<float> > partSizes_;
  int numberOfCriteria_;

  // Number of parts to be computed at each level in hierarchical partitioning.
  
  ArrayRCP<int> numberOfParts_;
  bool hierarchical_;
};
////////////////////////////////////////////////////////////////////////

#ifdef HAVE_ZOLTAN2_MPI
template <typename Adapter>
  PartitioningProblem<Adapter>::PartitioningProblem(Adapter *A, 
    ParameterList *p, MPI_Comm comm):
      Problem<Adapter>(A,p,comm), 
      generalParams_(), partitioningParams_(),solution_(),
      inputType_(InvalidAdapterType), modelType_(InvalidModel), algorithm_(),
      numberOfWeights_(), partIds_(), partSizes_(), 
      numberOfCriteria_(), numberOfParts_(), hierarchical_(false)
#else
template <typename Adapter>
  PartitioningProblem<Adapter>::PartitioningProblem(Adapter *A, 
    ParameterList *p):
      Problem<Adapter>(A,p), 
      generalParams_(), partitioningParams_(),solution_(),
      inputType_(InvalidAdapterType), modelType_(InvalidModel), algorithm_(),
      numberOfWeights_(), 
      partIds_(), partSizes_(), numberOfCriteria_(), 
      numberOfParts_(), hierarchical_(false)
#endif
{
  HELLO;
  createPartitioningProblem();
}

template <typename Adapter>
  void PartitioningProblem<Adapter>::setPartSizesForCritiera(
    int criteria, int len, size_t *partIds, float *partSizes, bool makeCopy) 
{
  Z2_LOCAL_INPUT_ASSERTION(*this->env_, "invalid length", len>= 0, BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*this->env_, "invalid criteria", 
    criteria >= 0 && criteria < numberOfWeights_, BASIC_ASSERTION);

  if (len == 0){
    partIds_[criteria] = ArrayRCP<size_t>();
    partSizes_[criteria] = ArrayRCP<float>();
    return;
  }

  Z2_LOCAL_INPUT_ASSERTION(*this->env_, "invalid arrays", partIds && partSizes, 
    BASIC_ASSERTION);

  // The global validity of the partIds and partSizes arrays is performed
  // by the PartitioningSolution, which computes global part distribution and
  // part sizes.

  size_t *z2_partIds = partIds;
  float *z2_partSizes = partSizes;
  bool own_memory = false;

  if (makeCopy){
    z2_partIds = NULL;
    z2_partIds = new size_t [len];
    Z2_LOCAL_MEMORY_ASSERTION(*this->env_, len, z2_partIds);
    z2_partSizes = NULL;
    z2_partSizes = new float [len];
    Z2_LOCAL_MEMORY_ASSERTION(*this->env_, len, z2_partSizes);
    bool own_memory = true;
  }

  partIds_[criteria] = arcp(z2_partIds, 0, len, own_memory);
  partSizes_[criteria] = arcp(z2_partSizes, 0, len, own_memory);
}

template <typename Adapter>
void PartitioningProblem<Adapter>::solve()
{
  HELLO;

  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

  // TODO: If hierarchical_

  // Create the solution. The algorithm will query the Solution
  //   for part and weight information. The algorithm will
  //   update the solution with part assignments and quality
  //   metrics.  The Solution object itself will convert our internal
  //   global numbers back to application global Ids.

  int weightDim = this->generalModel_->getNumWeights();

  RCP<const IdentifierMap<user_t> > idMap = 
    this->generalModel_->getIdentifierMap();

  solution_ = rcp(new PartitioningSolution<user_t>( this->envConst_,
    this->comm_, idMap, weightDim, partIds_.view(0, numberOfCriteria_), 
    partSizes_.view(0, numberOfCriteria_)));

  // Call the algorithm

  try {
    if (algorithm_ == string("scotch")){
      AlgPTScotch<base_adapter_t>(this->envConst_, this->comm_, 
        this->graphModel_, solution_);
    }
    else if (algorithm_ == string("block")){
      AlgPTBlock<base_adapter_t>(this->envConst_, this->comm_, 
        this->identifierModel_, solution_);
    }
    else{
      throw std::logic_error("partitioning algorithm not supported yet");
    }
  }
  Z2_FORWARD_EXCEPTIONS;
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

#ifdef HAVE_ZOLTAN2_OVIS
  ovis_enabled(this->comm_->getRank());
#endif

  // Finalize parameters.  If the Problem wants to set or
  // change any parameters, do it before this call.

  this->env_->commitParameters();

  // Get the parameters

  generalParams_ = Teuchos::Ptr<Teuchos::ParameterList>(
    &this->env_->getParamsNonConst());

  partitioningParams_ = Teuchos::Ptr<Teuchos::ParameterList>(
    &this->env_->getPartitioningParamsNonConst());

  string paramNotSet("unset");

  // Hierarchical partitioning?

  ParameterEntry *topo = partitioningParams_->getEntryPtr("topology");
  if (topo){
    Array<int> *values = NULL;
    Array<int> &valueList = topo->getValue(values);
    if (!Zoltan2::noValuesAreInRangeList(valueList)){
      int *n = new int [valueList.size() + 1];
      numberOfParts_ = arcp(n, 0, valueList.size() + 1, true);
      int procsPerNode = 1;
      for (int i=0; i < valueList.size(); i++){
        numberOfParts_[i+1] = valueList[i];
        procsPerNode *= valueList[i];
      }
      numberOfParts_[0] = this->env_->numProcs_ / procsPerNode;

      if (this->env_->numProcs_ % procsPerNode > 0)
        numberOfParts_[0]++;
    }
  }

  hierarchical_ = numberOfParts_.size() > 0;

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
      // TODO: add a parameter by which user tells us there are
      // no self edges to be removed.
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
#ifdef HAVE_ZOLTAN2_SCOTCH
      if (this->comm_->getSize() > 1)
        algorithm_ = string("ptscotch"); 
      else
        algorithm_ = string("scotch"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#else
#ifdef HAVE_ZOLTAN2_PARMETIS
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
             inputType_ == IdentifierAdapterType){
      modelType_ = IdentifierModelType;
      algorithm_ = string("block");
    }
    else{
      // This should never happen
      throw std::logic_error("input type is invalid");
    }
  }

  // Create the computational model.
  // Models are instantiated for base input adapter types (mesh,
  // matrix, graph, and so on).  We pass a pointer to the input
  // adapter, cast as the base input type.

  // TODO - check for exceptions

  typedef typename Adapter::base_adapter_t base_adapter_t;

  switch (modelType_) {

  case GraphModelType:
    this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
      this->baseInputAdapter_, this->envConst_, this->comm_, 
      needConsecutiveGlobalIds, removeSelfEdges));

    this->generalModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
      this->graphModel_);

    break;

  case HypergraphModelType:
    break;

  case GeometryModelType:
    break;

  case IdentifierModelType:
    this->identifierModel_ = rcp(new IdentifierModel<base_adapter_t>(
      this->baseInputAdapter_, this->envConst_, this->comm_,
      needConsecutiveGlobalIds));

    this->generalModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
      this->identifierModel_);
    break;

  default:
    cout << __func__ << " Invalid model" << modelType_ << endl;
    break;
  }

  // The Caller can specify part sizes in setPartSizes().  If he/she
  // does not, the part size arrays are empty.

  numberOfWeights_ = this->generalModel_->getNumWeights();

  numberOfCriteria_ = (numberOfWeights_ > 1) ? numberOfWeights_ : 1;

  ArrayRCP<size_t> *noIds = new ArrayRCP<size_t> [numberOfCriteria_];
  ArrayRCP<float> *noSizes = new ArrayRCP<float> [numberOfCriteria_];

  partIds_ = arcp(noIds, 0, numberOfCriteria_, true);
  partSizes_ = arcp(noSizes, 0, numberOfCriteria_, true);
}

}
#endif
