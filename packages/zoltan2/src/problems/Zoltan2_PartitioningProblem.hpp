// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_PartitioningProblem.hpp
    \brief Defines the PartitioningProblem class.
*/

#ifndef _ZOLTAN2_PARTITIONINGPROBLEM_HPP_
#define _ZOLTAN2_PARTITIONINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningAlgorithms.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningSolutionQuality.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_IntegerRangeList.hpp>
#include <Zoltan2_MachineRepresentation.hpp>
#include <Zoltan2_TaskMapping.hpp>

#ifndef _WIN32
#include <unistd.h>
#else
#include <process.h>
#define NOMINMAX
#include <windows.h>
#endif

#ifdef HAVE_ZOLTAN2_OVIS
#include <ovis.h>
#endif

namespace Zoltan2{

/*! \brief PartitioningProblem sets up partitioning problems for the user.
 *
 *  The PartitioningProblem is the core of the Zoltan2 partitioning API.
 *  Based on the the user's input and parameters, the PartitioningProblem
 *  sets up a computational Model, and a Solution object.  When the user
 *  calls the solve() method, the PartitioningProblem runs the algorithm,
 *  after which the Solution object may be obtained by the user.
 *  \todo include pointers to examples
 *
 *  The template parameter is the InputAdapter containing the data that
 *  is to be partitioned.
 *
 *  \todo hierarchical partitioning
 *  \todo repartition given an initial solution
 *  \todo follow partitioning with global or local ordering
 *  \todo allow unsetting of part sizes by passing in null pointers
 *  \todo add a parameter by which user tells us there are no self 
 *        edges to be removed.
 *  \todo - Should Problems and Solution have interfaces for returning
 *          views and for returning RCPs?  Or just one?  At a minimum, 
 *          we should have the word "View" in function names that return views.
 */

template<typename Adapter>
class PartitioningProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::zgid_t zgid_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

#ifdef HAVE_ZOLTAN2_MPI
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
  /*! \brief Constructor where MPI communicator can be specified
   */
  PartitioningProblem(Adapter *A, ParameterList *p, MPI_Comm comm):
      Problem<Adapter>(A,p,comm), solution_(),
      problemComm_(), problemCommConst_(),
      inputType_(InvalidAdapterType), 
      graphFlags_(), idFlags_(), coordFlags_(), algName_(),
      numberOfWeights_(), partIds_(), partSizes_(), 
      numberOfCriteria_(), levelNumberParts_(), hierarchical_(false), 
      timer_(), metricsRequested_(false), metrics_()
  {
    for(int i=0;i<MAX_NUM_MODEL_TYPES;i++) modelAvail_[i]=false;
    initializeProblem();
  }
#endif

  //! \brief Constructor where communicator is the Teuchos default.
  PartitioningProblem(Adapter *A, ParameterList *p):
      Problem<Adapter>(A,p), solution_(),
      problemComm_(), problemCommConst_(),
      inputType_(InvalidAdapterType), 
      graphFlags_(), idFlags_(), coordFlags_(), algName_(),
      numberOfWeights_(), 
      partIds_(), partSizes_(), numberOfCriteria_(), 
      levelNumberParts_(), hierarchical_(false), timer_(),
      metricsRequested_(false), metrics_()
  {
    for(int i=0;i<MAX_NUM_MODEL_TYPES;i++) modelAvail_[i]=false;
    initializeProblem();
  }

  //! \brief Constructor where Teuchos communicator is specified
  PartitioningProblem(Adapter *A, ParameterList *p,
                      RCP<const Teuchos::Comm<int> > &comm):
      Problem<Adapter>(A,p,comm), solution_(),
      problemComm_(), problemCommConst_(),
      inputType_(InvalidAdapterType), 
      graphFlags_(), idFlags_(), coordFlags_(), algName_(),
      numberOfWeights_(), 
      partIds_(), partSizes_(), numberOfCriteria_(), 
      levelNumberParts_(), hierarchical_(false), timer_(),
      metricsRequested_(false), metrics_()
  {
    for(int i=0;i<MAX_NUM_MODEL_TYPES;i++) modelAvail_[i]=false;
    initializeProblem();
  }

  /*! \brief Destructor
   */
  ~PartitioningProblem() {};

  //!  \brief Direct the problem to create a solution.
  //
  //    \param updateInputData   If true this indicates that either
  //          this is the first attempt at solution, or that we
  //          are computing a new solution and the input data has
  //          changed since the previous solution was computed.  
  //          By input data we mean coordinates, topology, or weights.
  //          If false, this indicates that we are computing a
  //          new solution using the same input data was used for
  //          the previous solution, even though the parameters
  //          may have been changed.
  //
  //  For the sake of performance, we ask the caller to set \c updateInputData
  //  to false if he/she is computing a new solution using the same input data,
  //  but different problem parameters, than that which was used to compute 
  //  the most recent solution.

  void solve(bool updateInputData=true );
 
  //!  \brief Return the part overlapping a given point in space; 
  //          when a point lies on a part boundary, the lowest part
  //          number on that boundary is returned.
  //          Note that not all partitioning algorithms will support
  //          this method.
  //
  //   \param dim : the number of dimensions specified for the point in space
  //   \param point : the coordinates of the point in space; array of size dim
  //   \return the part number of a part overlapping the given point
  //
  //   TODO:  This method might be more appropriate in a partitioning solution
  //   TODO:  But then the solution would need to know about the algorithm.
  //   TODO:  Consider moving the algorithm to the partitioning solution.
  part_t pointAssign(int dim, scalar_t *point) 
  {
    part_t p;
    try {
      if (this->algorithm_ == Teuchos::null)
        throw std::logic_error("no partitioning algorithm has been run yet");

      p = this->algorithm_->pointAssign(dim, point); 
    }
    Z2_FORWARD_EXCEPTIONS
    return p;
  }

  //!  \brief Return an array of all parts overlapping a given box in space.
  //   This method allocates memory for the return argument, but does not
  //   control that memory.  The user is responsible for freeing the 
  //   memory.
  //
  //   \param dim : (in) the number of dimensions specified for the box
  //   \param lower : (in) the coordinates of the lower corner of the box; 
  //                   array of size dim
  //   \param upper : (in) the coordinates of the upper corner of the box; 
  //                   array of size dim
  //   \param nPartsFound : (out) the number of parts overlapping the box
  //   \param partsFound :  (out) array of parts overlapping the box
  //   TODO:  This method might be more appropriate in a partitioning solution
  //   TODO:  But then the solution would need to know about the algorithm.
  //   TODO:  Consider moving the algorithm to the partitioning solution.
  void boxAssign(int dim, scalar_t *lower, scalar_t *upper,
                 size_t &nPartsFound, part_t **partsFound) 
  {
    try {
      if (this->algorithm_ == Teuchos::null)
        throw std::logic_error("no partitioning algorithm has been run yet");

      this->algorithm_->boxAssign(dim, lower, upper, nPartsFound, partsFound); 
    }
    Z2_FORWARD_EXCEPTIONS
  }

  //!  \brief Get the solution to the problem.
  //
  //   \return  a reference to the solution to the most recent solve().

  const PartitioningSolution<Adapter> &getSolution() {
    return *(solution_.getRawPtr());
  };

  /*! \brief Returns the imbalance of the solution.
   *   \param idx If there are multiple weights per object,
   *      specify the index for which the imbalance
   *      is desired, ranging from zero to one less then
   *     number of weights per object. 
   *   Imbalance was only computed if user requested
   *   metrics with a parameter.
   */
  const scalar_t getWeightImbalance(int idx=0) const {
    scalar_t imb = 0;
    if (!metrics_.is_null())
      metrics_->getWeightImbalance(imb, idx);

    return imb;
  }

  /*! \brief Get the array of metrics
   *   Metrics were only computed if user requested
   *   metrics with a parameter.
   */
  ArrayRCP<const MetricValues<scalar_t> > getMetrics() const {
   if (metrics_.is_null()){
      ArrayRCP<const MetricValues<scalar_t> > emptyMetrics;
      return emptyMetrics;
    }
    else
      return metrics_->getMetrics();
  }

  /*! \brief Print the array of metrics
   *   \param os the output stream for the report.
   *   Metrics were only computed if user requested
   *   metrics with a parameter.
   */
  void printMetrics(std::ostream &os) const {
    if (metrics_.is_null())
      os << "No metrics available." << endl;
    else
      metrics_->printMetrics(os);
  };

  /*! \brief Set or reset relative sizes for the parts that Zoltan2 will create.
   *
   *  \param len  The size of the \c partIds and \c partSizes lists
   *  \param partIds   A list of \c len part identifiers.  Part
   *           identifiers range from zero to one less than the global
   *           number of identifiers.  
   *  \param partSizes  A list of \c len relative sizes corresponding to
   *           the \c partIds.
   *  \param makeCopy  If true, Zoltan2 will make a copy of the ids and sizes
   *      that are provided in this call.  If false, Zoltan2 will just save
   *      the pointers to to the caller's lists.  If the pointers will remain
   *      remain valid throughout the lifetime of the PartitioningProblem,
   *      and memory use is an issue, then set makeCopy to false.  By default,
   *      Zoltan2 will copy the caller's list of ids and sizes.
   *
   * A given partid should be provided only once across all ranks.
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
   * Zoltan2 assumes that uniform part sizes are desired by the caller,
   * unless specified otherwise in a call to setPartSizes or 
   * setPartSizesForCriteria.
   *
   * \todo A user should be able to give us one set of part sizes
   *            that applies to all weight indices.  Right now
   *            for each weight index that does not have
   *            uniform part sizes, the user has to give us the
   *            part sizes once for each.
   */

  void setPartSizes(int len, part_t *partIds, scalar_t *partSizes, 
    bool makeCopy=true) 
  { 
    setPartSizesForCriteria(0, len, partIds, partSizes, makeCopy);
  }

  /*! \brief Set or reset the relative sizes (per weight) for the parts
   *    that Zoltan2 will create.
   *
   *  \param criteria the criteria for which these 
   *      part sizes apply.  Criteria range from zero to one less than
   *     the number of weights per object specified in the 
   *     caller's InputAdapter.
   *  \param len  The size of the \c partIds and \c partSizes lists
   *  \param partIds   A list of \c len part identifiers.  Part
   *           identifiers range from zero to one less than the global
   *           number of identifiers.  
   *  \param partSizes  A list of \c len relative sizes corresponding to
   *           the \c partIds.
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
   * Zoltan2 assumes that uniform part sizes are desired by the caller,
   * unless specified otherwise in a call to setPartSizes or 
   * setPartSizesForCriteria.
   */

  void setPartSizesForCriteria(int criteria, int len, part_t *partIds,
    scalar_t *partSizes, bool makeCopy=true) ;
/*
  void setMachine(MachineRepresentation<typename Adapter::base_adapter_t::scalar_t> *machine);
*/
  /*! \brief Reset the list of parameters
   */
  void resetParameters(ParameterList *params)
  {
    Problem<Adapter>::resetParameters(params);  // creates new environment
    if (timer_.getRawPtr() != NULL)
      this->env_->setTimer(timer_);
  }

  /*! \brief Get the current Environment.
   *   Useful for testing.
   */
  
  const RCP<const Environment> & getEnvironment() const 
  {
    return this->envConst_;
  }

private:
  void initializeProblem();

  void createPartitioningProblem(bool newData);

  RCP<PartitioningSolution<Adapter> > solution_;

  RCP<MachineRepresentation <typename Adapter::base_adapter_t::scalar_t>  > machine_;

  RCP<Comm<int> > problemComm_;
  RCP<const Comm<int> > problemCommConst_;

  BaseAdapterType inputType_;

  //ModelType modelType_;
  bool modelAvail_[MAX_NUM_MODEL_TYPES];

  modelFlag_t graphFlags_;
  modelFlag_t idFlags_;
  modelFlag_t coordFlags_;
  std::string algName_;

  int numberOfWeights_;

  // Suppose Array<part_t> partIds = partIds_[w].  If partIds.size() > 0
  // then the user supplied part sizes for weight index "w", and the sizes
  // corresponding to the Ids in partIds are partSizes[w].
  //
  // If numberOfWeights_ >= 0, then there is an Id and Sizes array for
  // for each weight.  Otherwise the user did not supply object weights,
  // but they can still specify part sizes. 
  // So numberOfCriteria_ is numberOfWeights_ or one, whichever is greater.

  ArrayRCP<ArrayRCP<part_t> > partIds_;
  ArrayRCP<ArrayRCP<scalar_t> > partSizes_;
  int numberOfCriteria_;

  // Number of parts to be computed at each level in hierarchical partitioning.
  
  ArrayRCP<int> levelNumberParts_;
  bool hierarchical_;

  // Create a Timer if the user asked for timing stats.

  RCP<TimerManager> timer_;

  // Did the user request metrics?

  bool metricsRequested_;
  RCP<const PartitioningSolutionQuality<Adapter> > metrics_;
};
////////////////////////////////////////////////////////////////////////

/*
template <typename Adapter>
void PartitioningProblem<Adapter>::setMachine(MachineRepresentation<typename Adapter::base_adapter_t::scalar_t> *machine){
  this->machine_ = RCP<MachineRepresentation<typename Adapter::base_adapter_t::scalar_t> > (machine, false);
}
*/


template <typename Adapter>
  void PartitioningProblem<Adapter>::initializeProblem()
{
  HELLO;

  this->env_->debug(DETAILED_STATUS, "PartitioningProblem::initializeProblem");

  if (getenv("DEBUGME")){
#ifndef _WIN32
    std::cout << getpid() << std::endl;
    sleep(15);
#else
    std::cout << _getpid() << std::endl;
    Sleep(15000);
#endif
  }

#ifdef HAVE_ZOLTAN2_OVIS
  ovis_enabled(this->comm_->getRank());
#endif

  // Create a copy of the user's communicator.  

  problemComm_ = this->comm_->duplicate();
  problemCommConst_ = rcp_const_cast<const Comm<int> > (problemComm_);

  machine_ = RCP <Zoltan2::MachineRepresentation<typename Adapter::scalar_t> >(new Zoltan2::MachineRepresentation<typename Adapter::scalar_t>(problemComm_));

  // Number of criteria is number of user supplied weights if non-zero.
  // Otherwise it is 1 and uniform weight is implied.

  numberOfWeights_ = this->inputAdapter_->getNumWeightsPerID();

  numberOfCriteria_ = (numberOfWeights_ > 1) ? numberOfWeights_ : 1;

  inputType_ = this->inputAdapter_->adapterType();

  // The Caller can specify part sizes in setPartSizes().  If he/she
  // does not, the part size arrays are empty.

  ArrayRCP<part_t> *noIds = new ArrayRCP<part_t> [numberOfCriteria_];
  ArrayRCP<scalar_t> *noSizes = new ArrayRCP<scalar_t> [numberOfCriteria_];

  partIds_ = arcp(noIds, 0, numberOfCriteria_, true);
  partSizes_ = arcp(noSizes, 0, numberOfCriteria_, true);

  if (this->env_->getDebugLevel() >= DETAILED_STATUS){
    std::ostringstream msg;
    msg << problemComm_->getSize() << " procs,"
      << numberOfWeights_ << " user-defined weights\n";
    this->env_->debug(DETAILED_STATUS, msg.str());
  }

  this->env_->memory("After initializeProblem");
}

template <typename Adapter>
  void PartitioningProblem<Adapter>::setPartSizesForCriteria(
    int criteria, int len, part_t *partIds, scalar_t *partSizes, bool makeCopy) 
{
  this->env_->localInputAssertion(__FILE__, __LINE__, "invalid length", 
    len>= 0, BASIC_ASSERTION);

  this->env_->localInputAssertion(__FILE__, __LINE__, "invalid criteria", 
    criteria >= 0 && criteria < numberOfCriteria_, BASIC_ASSERTION);

  if (len == 0){
    partIds_[criteria] = ArrayRCP<part_t>();
    partSizes_[criteria] = ArrayRCP<scalar_t>();
    return;
  }

  this->env_->localInputAssertion(__FILE__, __LINE__, "invalid arrays", 
    partIds && partSizes, BASIC_ASSERTION);

  // The global validity of the partIds and partSizes arrays is performed
  // by the PartitioningSolution, which computes global part distribution and
  // part sizes.

  part_t *z2_partIds = NULL;
  scalar_t *z2_partSizes = NULL;
  bool own_memory = false;

  if (makeCopy){
    z2_partIds = new part_t [len];
    z2_partSizes = new scalar_t [len];
    this->env_->localMemoryAssertion(__FILE__, __LINE__, len, z2_partSizes);
    memcpy(z2_partIds, partIds, len * sizeof(part_t));
    memcpy(z2_partSizes, partSizes, len * sizeof(scalar_t));
    own_memory=true;
  }
  else{
    z2_partIds = partIds;
    z2_partSizes = partSizes;
  }

  partIds_[criteria] = arcp(z2_partIds, 0, len, own_memory);
  partSizes_[criteria] = arcp(z2_partSizes, 0, len, own_memory);
}

template <typename Adapter>
void PartitioningProblem<Adapter>::solve(bool updateInputData)
{
  HELLO;
  this->env_->debug(DETAILED_STATUS, "Entering solve");

  // Create the computational model.

  this->env_->timerStart(MACRO_TIMERS, "create problem");

  createPartitioningProblem(updateInputData);

  this->env_->timerStop(MACRO_TIMERS, "create problem");

  // TODO: If hierarchical_

  // Create the solution. The algorithm will query the Solution
  //   for part and weight information. The algorithm will
  //   update the solution with part assignments and quality
  //   metrics.  The Solution object itself will convert our internal
  //   global numbers back to application global Ids if needed.

  RCP<const IdentifierMap<user_t> > idMap = 
    this->baseModel_->getIdentifierMap();

  PartitioningSolution<Adapter> *soln = NULL;

  this->env_->timerStart(MACRO_TIMERS, "create solution");

  try{
    soln = new PartitioningSolution<Adapter>( 
      this->envConst_, problemCommConst_, idMap, numberOfWeights_, 
      partIds_.view(0, numberOfCriteria_), 
      partSizes_.view(0, numberOfCriteria_));
  }
  Z2_FORWARD_EXCEPTIONS;

  solution_ = rcp(soln);

  this->env_->timerStop(MACRO_TIMERS, "create solution");

  this->env_->memory("After creating Solution");

  // Call the algorithm

  try {
    if (algName_ == std::string("scotch")) {

      this->algorithm_ = rcp(new AlgPTScotch<Adapter>(this->envConst_,
                                            problemComm_,
                                            this->graphModel_));
      this->algorithm_->partition(solution_);
    }

    else if (algName_ == std::string("parmetis")) {

      this->algorithm_ = rcp(new AlgParMETIS<Adapter>(this->envConst_,
                                            problemComm_,
                                            this->graphModel_));
      this->algorithm_->partition(solution_);
    }

    else if (algName_ == std::string("block")) {

      this->algorithm_ = rcp(new AlgBlock<Adapter>(this->envConst_,
                                         problemComm_, this->identifierModel_));
      this->algorithm_->partition(solution_);
    }

    else if (algName_ == std::string("rcb")) {

      this->algorithm_ = rcp(new AlgRCB<Adapter>(this->envConst_, problemComm_,
                                                 this->coordinateModel_));
      this->algorithm_->partition(solution_);
    }

    else if (algName_ == std::string("multijagged")) {

      this->algorithm_ = rcp(new Zoltan2_AlgMJ<Adapter>(this->envConst_,
                                              problemComm_,
                                              this->coordinateModel_));
      this->algorithm_->partition(solution_);
    }

    else if (algName_ == std::string("wolf")) {

      this->algorithm_ = rcp(new AlgWolf<Adapter>(this->envConst_,
                                        problemComm_,this->graphModel_,
                                        this->coordinateModel_));

      // need to add coordModel, make sure this is built
      this->algorithm_->partition(solution_);
    }
    else {
      throw std::logic_error("partitioning algorithm not supported");
    }
  }
  Z2_FORWARD_EXCEPTIONS;

  //if mapping is requested
  const Teuchos::ParameterEntry *pe = this->envConst_->getParameters().getEntryPtr("mapping_type");
  int mapping_type = -1;
  if (pe){
    mapping_type = pe->getValue(&mapping_type);
  }
  //if mapping is 0 -- coordinate mapping

  if (mapping_type == 0){

    //part_t *task_communication_xadj = NULL, *task_communication_adj = NULL;

    Zoltan2::CoordinateTaskMapper <Adapter, part_t> *ctm=
                  new Zoltan2::CoordinateTaskMapper<Adapter,part_t>(
                          problemComm_.getRawPtr(),
                          machine_.getRawPtr(),
                          this->coordinateModel_.getRawPtr(),
                          solution_.getRawPtr(),
                          this->envConst_.getRawPtr()
                          //,task_communication_xadj,
                          //task_communication_adj
                          );
    //for now just delete the object.
    delete ctm;
  }
  else if (mapping_type == 1){
    //if mapping is 1 -- graph mapping
  }

  if (metricsRequested_){
    typedef PartitioningSolution<Adapter> ps_t;
    typedef PartitioningSolutionQuality<Adapter> psq_t;

    psq_t *quality = NULL;
    RCP<const ps_t> solutionConst = rcp_const_cast<const ps_t>(solution_);
    RCP<const Adapter> adapter = rcp(this->inputAdapter_, false);

    try{
      quality = new psq_t(this->envConst_, problemCommConst_, adapter, 
        solutionConst);
    }
    Z2_FORWARD_EXCEPTIONS

    metrics_ = rcp(quality);
  }

  this->env_->debug(DETAILED_STATUS, "Exiting solve");
}

template <typename Adapter>
void PartitioningProblem<Adapter>::createPartitioningProblem(bool newData)
{
  HELLO;
  this->env_->debug(DETAILED_STATUS, 
    "PartitioningProblem::createPartitioningProblem");

  using std::string;
  using Teuchos::ParameterList;

  // A Problem object may be reused.  The input data may have changed and
  // new parameters or part sizes may have been set.
  //
  // Save these values in order to determine if we need to create a new model.

  //ModelType previousModel = modelType_;
  bool prevModelAvail[MAX_NUM_MODEL_TYPES];
  for(int i=0;i<MAX_NUM_MODEL_TYPES;i++)
  {
    prevModelAvail[i] = modelAvail_[i];
  }


  modelFlag_t previousGraphModelFlags = graphFlags_;
  modelFlag_t previousIdentifierModelFlags = idFlags_;
  modelFlag_t previousCoordinateModelFlags = coordFlags_;

  //modelType_ = InvalidModel;
  for(int i=0;i<MAX_NUM_MODEL_TYPES;i++)
  {
    modelAvail_[i] = false;
  }

  graphFlags_.reset();
  idFlags_.reset();
  coordFlags_.reset();

  ////////////////////////////////////////////////////////////////////////////
  // It's possible at this point that the Problem may want to
  // add problem parameters to the parameter list in the Environment. 
  //
  // Since the parameters in the Environment have already been
  // validated in its constructor, a new Environment must be created:
  ////////////////////////////////////////////////////////////////////////////
  // Teuchos::RCP<const Teuchos::Comm<int> > oldComm = this->env_->comm_;
  // const ParameterList &oldParams = this->env_->getUnvalidatedParameters();
  // 
  // ParameterList newParams = oldParams;
  // newParams.set("new_parameter", "new_value");
  // 
  // ParameterList &newPartParams = newParams.sublist("partitioning");
  // newPartParams.set("new_partitioning_parameter", "its_value");
  // 
  // this->env_ = rcp(new Environment(newParams, oldComm));
  ////////////////////////////////////////////////////////////////////////////

  this->env_->debug(DETAILED_STATUS, "    parameters");
  Environment &env = *(this->env_);
  ParameterList &pl = env.getParametersNonConst();

  std::string defString("default");

  // Did the user ask for computation of quality metrics?

  int yesNo=0;
  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("compute_metrics");
  if (pe){
    yesNo = pe->getValue<int>(&yesNo);
    metricsRequested_ = true;
  }

  // Did the user specify a computational model?

  std::string model(defString);
  pe = pl.getEntryPtr("model");
  if (pe)
    model = pe->getValue<std::string>(&model);

  // Did the user specify an algorithm?

  std::string algorithm(defString);
  pe = pl.getEntryPtr("algorithm");
  if (pe)
    algorithm = pe->getValue<std::string>(&algorithm);

  // Possible algorithm requirements that must be conveyed to the model:

  bool needConsecutiveGlobalIds = false;
  bool removeSelfEdges= false;

  ///////////////////////////////////////////////////////////////////
  // Determine algorithm, model, and algorithm requirements.  This
  // is a first pass.  Feel free to change this and add to it.

  if (algorithm != defString)
  {

    // Figure out the model required by the algorithm
    if (algorithm == std::string("block") ||
        algorithm == std::string("random") ||
        algorithm == std::string("cyclic") ){

      //modelType_ = IdentifierModelType;
      modelAvail_[IdentifierModelType] = true;

      algName_ = algorithm;
      needConsecutiveGlobalIds = true;
    }
    else if (algorithm == std::string("rcb") ||
             algorithm == std::string("rib") ||
             algorithm == std::string("multijagged") ||
             algorithm == std::string("hsfc"))
    {
      //modelType_ = CoordinateModelType;
      modelAvail_[CoordinateModelType]=true;
    
      algName_ = algorithm;
    }
    else if (algorithm == std::string("metis") ||
             algorithm == std::string("parmetis") ||
             algorithm == std::string("scotch") ||
             algorithm == std::string("ptscotch"))
    {

      //modelType_ = GraphModelType;
      modelAvail_[GraphModelType]=true;

      algName_ = algorithm;
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
    }
    else if (algorithm == std::string("patoh") ||
             algorithm == std::string("phg"))
    {
      // if ((modelType_ != GraphModelType) &&
      //     (modelType_ != HypergraphModelType) )
      if ((modelAvail_[GraphModelType]==false) &&
          (modelAvail_[HypergraphModelType]==false) )
      {
        //modelType_ = HypergraphModelType;
        modelAvail_[HypergraphModelType]=true;
      }
      algName_ = algorithm;
      needConsecutiveGlobalIds = true;
    }
    else if (algorithm == std::string("wolf"))
    {
      modelAvail_[GraphModelType]=true;
      modelAvail_[CoordinateModelType]=true;
      algName_ = algorithm;
    }
    else
    {
      // Parameter list should ensure this does not happen.
      throw std::logic_error("parameter list algorithm is invalid");
    }
  }
  else if (model != defString)
  {
    // Figure out the algorithm suggested by the model.
    if (model == std::string("hypergraph"))
    {      
      //modelType_ = HypergraphModelType;
      modelAvail_[HypergraphModelType]=true;

      if (problemComm_->getSize() > 1)
        algName_ = std::string("phg"); 
      else
        algName_ = std::string("patoh"); 
      needConsecutiveGlobalIds = true;
    }
    else if (model == std::string("graph"))
    {
      //modelType_ = GraphModelType;
      modelAvail_[GraphModelType]=true;

#ifdef HAVE_ZOLTAN2_SCOTCH
      if (problemComm_->getSize() > 1)
        algName_ = std::string("ptscotch"); 
      else
        algName_ = std::string("scotch"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#else
#ifdef HAVE_ZOLTAN2_PARMETIS
      if (problemComm_->getSize() > 1)
        algName_ = std::string("parmetis"); 
      else
        algName_ = std::string("metis"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#else
      if (problemComm_->getSize() > 1)
        algName_ = std::string("phg"); 
      else
        algName_ = std::string("patoh"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#endif
#endif
    }
    else if (model == std::string("geometry"))
    {
      //modelType_ = CoordinateModelType;
      modelAvail_[CoordinateModelType]=true;

      algName_ = std::string("rcb");
    }
    else if (model == std::string("ids"))
    {
      //modelType_ = IdentifierModelType;
      modelAvail_[IdentifierModelType]=true;

      algName_ = std::string("block");
      needConsecutiveGlobalIds = true;
    }
    else
    {
      // Parameter list should ensure this does not happen.
      env.localBugAssertion(__FILE__, __LINE__, 
        "parameter list model type is invalid", 1, BASIC_ASSERTION);
    }
  }
  else
  {   
    // Determine an algorithm and model suggested by the input type.
    //   TODO: this is a good time to use the time vs. quality parameter
    //     in choosing an algorithm, and setting some parameters

    if (inputType_ == MatrixAdapterType)
    {
      //modelType_ = HypergraphModelType;
      modelAvail_[HypergraphModelType]=true;
      
      if (problemComm_->getSize() > 1)
        algName_ = std::string("phg"); 
      else
        algName_ = std::string("patoh"); 
    }
    else if (inputType_ == GraphAdapterType ||
        inputType_ == MeshAdapterType)
    {
      //modelType_ = GraphModelType;
      modelAvail_[GraphModelType]=true;

      if (problemComm_->getSize() > 1)
        algName_ = std::string("phg"); 
      else
        algName_ = std::string("patoh"); 
    }
    else if (inputType_ == CoordinateAdapterType)
    {
      //modelType_ = CoordinateModelType;
      modelAvail_[CoordinateModelType]=true;

      if(algName_ != std::string("multijagged"))
      algName_ = std::string("rcb");
    }
    else if (inputType_ == VectorAdapterType ||
             inputType_ == IdentifierAdapterType)
    {
      //modelType_ = IdentifierModelType;
      modelAvail_[IdentifierModelType]=true;

      algName_ = std::string("block");
    }
    else{
      // This should never happen
      throw std::logic_error("input type is invalid");
    }
  }

  // Hierarchical partitioning?

  Array<int> valueList;
  pe = pl.getEntryPtr("topology");

  if (pe){
    valueList = pe->getValue<Array<int> >(&valueList);

    if (!Zoltan2::noValuesAreInRangeList<int>(valueList)){
      int *n = new int [valueList.size() + 1];
      levelNumberParts_ = arcp(n, 0, valueList.size() + 1, true);
      int procsPerNode = 1;
      for (int i=0; i < valueList.size(); i++){
        levelNumberParts_[i+1] = valueList[i];
        procsPerNode *= valueList[i];
      }
      // Number of parts in the first level
      levelNumberParts_[0] = env.numProcs_ / procsPerNode;

      if (env.numProcs_ % procsPerNode > 0)
        levelNumberParts_[0]++;
    }
  }
  else{
    levelNumberParts_.clear();
  }

  hierarchical_ = levelNumberParts_.size() > 0;

  // Object to be partitioned? (rows, columns, etc)

  std::string objectOfInterest(defString);
  pe = pl.getEntryPtr("objects_to_partition");
  if (pe)
    objectOfInterest = pe->getValue<std::string>(&objectOfInterest);

  ///////////////////////////////////////////////////////////////////
  // Set model creation flags, if any.

  this->env_->debug(DETAILED_STATUS, "    models");
  //  if (modelType_ == GraphModelType)
  if (modelAvail_[GraphModelType]==true)
  {

    // Any parameters in the graph sublist?

    std::string symParameter(defString);
    pe = pl.getEntryPtr("symmetrize_graph");
    if (pe){
      symParameter = pe->getValue<std::string>(&symParameter);
      if (symParameter == std::string("transpose"))
        graphFlags_.set(SYMMETRIZE_INPUT_TRANSPOSE);
      else if (symParameter == std::string("bipartite"))
        graphFlags_.set(SYMMETRIZE_INPUT_BIPARTITE);
    } 

    int sgParameter = 0;
    pe = pl.getEntryPtr("subset_graph");
    if (pe)
      sgParameter = pe->getValue<int>(&sgParameter);

    if (sgParameter == 1)
        graphFlags_.set(GRAPH_IS_A_SUBSET_GRAPH);

    // Any special behaviors required by the algorithm?
    
    if (removeSelfEdges)
      graphFlags_.set(SELF_EDGES_MUST_BE_REMOVED);

    if (needConsecutiveGlobalIds)
      graphFlags_.set(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

    // How does user input map to vertices and edges?

    if (inputType_ == MatrixAdapterType){
      if (objectOfInterest == defString ||
          objectOfInterest == std::string("matrix_rows") )
        graphFlags_.set(VERTICES_ARE_MATRIX_ROWS);
      else if (objectOfInterest == std::string("matrix_columns"))
        graphFlags_.set(VERTICES_ARE_MATRIX_COLUMNS);
      else if (objectOfInterest == std::string("matrix_nonzeros"))
        graphFlags_.set(VERTICES_ARE_MATRIX_NONZEROS);
    }

    else if (inputType_ == MeshAdapterType){
      if (objectOfInterest == defString ||
          objectOfInterest == std::string("mesh_nodes") )
        graphFlags_.set(VERTICES_ARE_MESH_NODES);
      else if (objectOfInterest == std::string("mesh_elements"))
        graphFlags_.set(VERTICES_ARE_MESH_ELEMENTS);
    } 
  }
  //MMW is it ok to remove else?
  //  else if (modelType_ == IdentifierModelType)
  if (modelAvail_[IdentifierModelType]==true)
  {

    // Any special behaviors required by the algorithm?
    
    if (needConsecutiveGlobalIds)
      idFlags_.set(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  }
  //  else if (modelType_ == CoordinateModelType)
  if (modelAvail_[CoordinateModelType]==true)
  {

    // Any special behaviors required by the algorithm?
    
    if (needConsecutiveGlobalIds)
      coordFlags_.set(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  }


  if ( newData ||
       (modelAvail_[GraphModelType]!=prevModelAvail[GraphModelType]) ||
       (modelAvail_[HypergraphModelType]!=prevModelAvail[HypergraphModelType]) ||
       (modelAvail_[CoordinateModelType]!=prevModelAvail[CoordinateModelType]) ||
       (modelAvail_[IdentifierModelType]!=prevModelAvail[IdentifierModelType]) ||
	//       (modelType_ != previousModel) ||
       (graphFlags_ != previousGraphModelFlags) ||
       (coordFlags_ != previousCoordinateModelFlags) ||
       (idFlags_ != previousIdentifierModelFlags) ) 
  {

    // Create the computational model.
    // Models are instantiated for base input adapter types (mesh,
    // matrix, graph, and so on).  We pass a pointer to the input
    // adapter, cast as the base input type.

    //KDD Not sure why this shadow declaration is needed
    //KDD Comment out for now; revisit later if problems.
    //KDD const Teuchos::ParameterList pl = this->envConst_->getParameters();
    //bool exceptionThrow = true;

    if(modelAvail_[GraphModelType]==false && modelAvail_[HypergraphModelType]==false &&
       modelAvail_[CoordinateModelType]==false && modelAvail_[IdentifierModelType]==false)
    {
      cout << __func__zoltan2__ << " Invalid model"  << endl;
    }
    else
    {
      if(modelAvail_[GraphModelType]==true)
      {
        this->env_->debug(DETAILED_STATUS, "    building graph model");
        this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
          this->baseInputAdapter_, this->envConst_, problemComm_, graphFlags_));

        this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(this->graphModel_);
      }
      if(modelAvail_[HypergraphModelType]==true)
      {
	std::cout << "Hypergraph model not implemented yet..." << std::endl;
      }

      if(modelAvail_[CoordinateModelType]==true)
      {
      	this->env_->debug(DETAILED_STATUS, "    building coordinate model");
      	this->coordinateModel_ = rcp(new CoordinateModel<base_adapter_t>(
      				     this->baseInputAdapter_, this->envConst_, problemComm_, coordFlags_));

        this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(this->coordinateModel_);
      }

      if(modelAvail_[IdentifierModelType]==true)
      {
        this->env_->debug(DETAILED_STATUS, "    building identifier model");
        this->identifierModel_ = rcp(new IdentifierModel<base_adapter_t>(
                                     this->baseInputAdapter_, this->envConst_, problemComm_, idFlags_));

        this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(this->identifierModel_);
      }
  

    }



    this->env_->memory("After creating Model");
    this->env_->debug(DETAILED_STATUS, "createPartitioningProblem done");
  }

}

}  // namespace Zoltan2
#endif
