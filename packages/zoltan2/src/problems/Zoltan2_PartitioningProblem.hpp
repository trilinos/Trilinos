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
#include <Zoltan2_EvaluatePartition.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_IntegerRangeList.hpp>
#include <Zoltan2_MachineRepresentation.hpp>
#include <Zoltan2_AlgSerialGreedy.hpp>
#ifdef ZOLTAN2_TASKMAPPING_MOVE
#include <Zoltan2_TaskMapping.hpp>
#endif

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
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

  //! \brief Constructor where Teuchos communicator is specified
  PartitioningProblem(Adapter *A, ParameterList *p,
                      const RCP<const Teuchos::Comm<int> > &comm):
      Problem<Adapter>(A,p,comm), 
      solution_(),
      inputType_(InvalidAdapterType),
      graphFlags_(), idFlags_(), coordFlags_(),
      algName_(), numberOfWeights_(), partIds_(), partSizes_(),
      numberOfCriteria_(), levelNumberParts_(), hierarchical_(false)
  {
    for(int i=0;i<MAX_NUM_MODEL_TYPES;i++) modelAvail_[i]=false;
    initializeProblem();
  }

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor where MPI communicator can be specified
   */
  PartitioningProblem(Adapter *A, ParameterList *p, MPI_Comm mpicomm):
  PartitioningProblem(A, p, 
                      rcp<const Comm<int> >(new Teuchos::MpiComm<int>(
                                            Teuchos::opaqueWrapper(mpicomm))))
  {}
#endif

  //! \brief Constructor where communicator is the Teuchos default.
  PartitioningProblem(Adapter *A, ParameterList *p):
  PartitioningProblem(A, p, Tpetra::getDefaultComm())
  {}

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

  void solve(bool updateInputData=true);

  //!  \brief Get the solution to the problem.
  //
  //   \return  a reference to the solution to the most recent solve().

  const PartitioningSolution<Adapter> &getSolution() {
    return *(solution_.getRawPtr());
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

  /*! \brief Set up validators specific to this Problem
  */
  static void getValidParameters(ParameterList & pl)
  {
    Zoltan2_AlgMJ<Adapter>::getValidParameters(pl);
    AlgPuLP<Adapter>::getValidParameters(pl);
    AlgPTScotch<Adapter>::getValidParameters(pl);
    AlgSerialGreedy<Adapter>::getValidParameters(pl);
    AlgForTestingOnly<Adapter>::getValidParameters(pl);

    // This set up does not use tuple because we didn't have constructors
    // that took that many elements - Tuple will need to be modified and I
    // didn't want to have low level changes with this particular refactor
    // TO DO: Add more Tuple constructors and then redo this code to be
    //  Teuchos::tuple<std::string> algorithm_names( "rcb", "multijagged" ... );
    Array<std::string> algorithm_names(17);
    algorithm_names[0] = "rcb";
    algorithm_names[1] = "multijagged";
    algorithm_names[2] = "rib";
    algorithm_names[3] = "hsfc";
    algorithm_names[4] = "patoh";
    algorithm_names[5] = "phg";
    algorithm_names[6] = "metis";
    algorithm_names[7] = "parmetis";
    algorithm_names[8] = "pulp";
    algorithm_names[9] = "parma";
    algorithm_names[10] = "scotch";
    algorithm_names[11] = "ptscotch";
    algorithm_names[12] = "block";
    algorithm_names[13] = "cyclic";
    algorithm_names[14] = "random";
    algorithm_names[15] = "zoltan";
    algorithm_names[16] = "forTestingOnly";
    RCP<Teuchos::StringValidator> algorithm_Validator = Teuchos::rcp(
      new Teuchos::StringValidator( algorithm_names ));
    pl.set("algorithm", "random", "partitioning algorithm",
      algorithm_Validator);

    // bool parameter
    pl.set("keep_partition_tree", false, "If true, will keep partition tree",
      Environment::getBoolValidator());

    // bool parameter
    pl.set("rectilinear", false, "If true, then when a cut is made, all of the "
      "dots located on the cut are moved to the same side of the cut. The "
      "resulting regions are then rectilinear. The resulting load balance may "
      "not be as good as when the group of dots is split by the cut. ",
      Environment::getBoolValidator());

    RCP<Teuchos::StringValidator> partitioning_objective_Validator =
      Teuchos::rcp( new Teuchos::StringValidator(
       Teuchos::tuple<std::string>( "balance_object_count",
         "balance_object_weight", "multicriteria_minimize_total_weight",
         "multicriteria_minimize_maximum_weight",
         "multicriteria_balance_total_maximum", "minimize_cut_edge_count",
         "minimize_cut_edge_weight", "minimize_neighboring_parts",
         "minimize_boundary_vertices" )));
    pl.set("partitioning_objective", "balance_object_weight",
      "objective of partitioning", partitioning_objective_Validator);

    pl.set("imbalance_tolerance", 1.1, "imbalance tolerance, ratio of "
      "maximum load over average load", Environment::getAnyDoubleValidator());

    // num_global_parts >= 1
    RCP<Teuchos::EnhancedNumberValidator<int>> num_global_parts_Validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(
        1, Teuchos::EnhancedNumberTraits<int>::max()) ); // no maximum
    pl.set("num_global_parts", 1, "global number of parts to compute "
      "(0 means use the number of processes)", num_global_parts_Validator);

    // num_local_parts >= 0
    RCP<Teuchos::EnhancedNumberValidator<int>> num_local_parts_Validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(
        0, Teuchos::EnhancedNumberTraits<int>::max()) ); // no maximum
    pl.set("num_local_parts", 0, "number of parts to compute for this "
      "process (num_global_parts == sum of all num_local_parts)", 
      num_local_parts_Validator);

    RCP<Teuchos::StringValidator> partitioning_approach_Validator =
      Teuchos::rcp( new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "partition", "repartition",
          "maximize_overlap" )));
    pl.set("partitioning_approach", "partition", "Partition from scratch, "
      "partition incrementally from current partition, of partition from "
      "scratch but maximize overlap  with the current partition",
      partitioning_approach_Validator);

    RCP<Teuchos::StringValidator> objects_to_partition_Validator =
      Teuchos::rcp( new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "matrix_rows", "matrix_columns",
          "matrix_nonzeros", "mesh_elements", "mesh_nodes", "graph_edges",
        "graph_vertices", "coordinates", "identifiers" )));
    pl.set("objects_to_partition", "graph_vertices", "Objects to be partitioned",
      objects_to_partition_Validator);

    RCP<Teuchos::StringValidator> model_Validator = Teuchos::rcp(
      new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "hypergraph", "graph",
          "geometry", "ids" )));
    pl.set("model", "graph", "This is a low level parameter. Normally the "
      "library will choose a computational model based on the algorithm or "
      "objective specified by the user.", model_Validator);

    // bool parameter
    pl.set("remap_parts", false, "remap part numbers to minimize migration "
      "between old and new partitions", Environment::getBoolValidator() );

    pl.set("mapping_type", -1, "Mapping of solution to the processors. -1 No"
      " Mapping, 0 coordinate mapping.", Environment::getAnyIntValidator());

    RCP<Teuchos::EnhancedNumberValidator<int>> ghost_layers_Validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(1, 10, 1, 0) );
    pl.set("ghost_layers", 2, "number of layers for ghosting used in "
      "hypergraph ghost method", ghost_layers_Validator);
  }

private:
  void initializeProblem();

  void createPartitioningProblem(bool newData);

  RCP<PartitioningSolution<Adapter> > solution_;
#ifdef ZOLTAN2_TASKMAPPING_MOVE
  RCP<MachineRepresentation<scalar_t,part_t> > machine_;
#endif

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

#ifdef ZOLTAN2_TASKMAPPING_MOVE
  machine_ = RCP<MachineRepresentation<scalar_t,part_t> >(
                 new MachineRepresentation<scalar_t,part_t>(*(this->comm_)));
#endif

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
    msg << this->comm_->getSize() << " procs,"
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
  this->env_->debug(DETAILED_STATUS, "Entering solve");

  // Create the computational model.

  this->env_->timerStart(MACRO_TIMERS, "create problem");

  createPartitioningProblem(updateInputData);

  this->env_->timerStop(MACRO_TIMERS, "create problem");

  // TODO: If hierarchical_

  // Create the solution. The algorithm will query the Solution
  // for part and weight information. The algorithm will
  // update the solution with part assignments and quality metrics.

  // Create the algorithm
  try {
    if (algName_ == std::string("multijagged")) {
      this->algorithm_ = rcp(new Zoltan2_AlgMJ<Adapter>(this->envConst_,
                                              this->comm_,
                                              this->coordinateModel_));
    }
    else if (algName_ == std::string("zoltan")) {
      this->algorithm_ = rcp(new AlgZoltan<Adapter>(this->envConst_,
                                           this->comm_,
                                           this->baseInputAdapter_));
    }
    else if (algName_ == std::string("parma")) {
      this->algorithm_ = rcp(new AlgParMA<Adapter>(this->envConst_,
                                           this->comm_,
                                           this->baseInputAdapter_));
    }
    else if (algName_ == std::string("scotch")) {
      this->algorithm_ = rcp(new AlgPTScotch<Adapter>(this->envConst_,
                                            this->comm_,
                                            this->baseInputAdapter_));
    }
    else if (algName_ == std::string("parmetis")) {
      this->algorithm_ = rcp(new AlgParMETIS<Adapter>(this->envConst_,
                                            this->comm_,
                                            this->graphModel_));
    }
    else if (algName_ == std::string("pulp")) {
      this->algorithm_ = rcp(new AlgPuLP<Adapter>(this->envConst_,
                                            this->comm_,
                                            this->baseInputAdapter_));
    }
    else if (algName_ == std::string("block")) {
      this->algorithm_ = rcp(new AlgBlock<Adapter>(this->envConst_,
                                         this->comm_, this->identifierModel_));
    }
    else if (algName_ == std::string("phg") ||
             algName_ == std::string("patoh")) {
      // phg and patoh provided through Zoltan
      Teuchos::ParameterList &pl = this->env_->getParametersNonConst();
      Teuchos::ParameterList &zparams = pl.sublist("zoltan_parameters",false);
      if (numberOfWeights_ > 0) {
        char strval[20];
        sprintf(strval, "%d", numberOfWeights_);
        zparams.set("OBJ_WEIGHT_DIM", strval);
      }
      zparams.set("LB_METHOD", algName_.c_str());
      zparams.set("LB_APPROACH", "PARTITION"); 
      algName_ = std::string("zoltan");

      this->algorithm_ = rcp(new AlgZoltan<Adapter>(this->envConst_,
                                           this->comm_,
                                           this->baseInputAdapter_));
    }
    else if (algName_ == std::string("forTestingOnly")) {
      this->algorithm_ = rcp(new AlgForTestingOnly<Adapter>(this->envConst_,
                                           this->comm_,
                                           this->baseInputAdapter_));
    }
    // else if (algName_ == std::string("rcb")) {
    //  this->algorithm_ = rcp(new AlgRCB<Adapter>(this->envConst_,this->comm_,
    //                                             this->coordinateModel_));
    // }
    else {
      throw std::logic_error("partitioning algorithm not supported");
    }
  }
  Z2_FORWARD_EXCEPTIONS;

  // Create the solution
  this->env_->timerStart(MACRO_TIMERS, "create solution");
  PartitioningSolution<Adapter> *soln = NULL;

  try{
    soln = new PartitioningSolution<Adapter>(
      this->envConst_, this->comm_, numberOfWeights_,
      partIds_.view(0, numberOfCriteria_),
      partSizes_.view(0, numberOfCriteria_), this->algorithm_);
  }
  Z2_FORWARD_EXCEPTIONS;

  solution_ = rcp(soln);

  this->env_->timerStop(MACRO_TIMERS, "create solution");
  this->env_->memory("After creating Solution");

  // Call the algorithm

  try {
    this->algorithm_->partition(solution_);
  }
  Z2_FORWARD_EXCEPTIONS;

  //if mapping is requested
  const Teuchos::ParameterEntry *pe = this->envConst_->getParameters().getEntryPtr("mapping_type");
  int mapping_type = -1;
  if (pe){
    mapping_type = pe->getValue(&mapping_type);
  }
  //if mapping is 0 -- coordinate mapping

#if ZOLTAN2_TASKMAPPING_MOVE
  if (mapping_type == 0){

    //part_t *task_communication_xadj = NULL, *task_communication_adj = NULL;

    Zoltan2::CoordinateTaskMapper <Adapter, part_t> *ctm=
                  new Zoltan2::CoordinateTaskMapper<Adapter,part_t>(
                          this->comm_.getRawPtr(),
                          machine_.getRawPtr(),
                          this->coordinateModel_.getRawPtr(),
                          solution_.getRawPtr(),
                          this->envConst_.getRawPtr()
                          //,task_communication_xadj,
                          //task_communication_adj
                          );

    // KDD  For now, we would need to re-map the part numbers in the solution.
    // KDD  I suspect we'll later need to distinguish between part numbers and 
    // KDD  process numbers to provide separation between partitioning and
    // KDD  mapping.  For example, does this approach here assume #parts == #procs?
    // KDD  If we map k tasks to p processes with k > p, do we effectively reduce
    // KDD  the number of tasks (parts) in the solution?

#ifdef KDD_READY
    const part_t *oldParts = solution_->getPartListView();
    size_t nLocal = ia->getNumLocalIds();
    for (size_t i = 0; i < nLocal; i++) {
      // kind of cheating since oldParts is a view; probably want an interface in solution 
      // for resetting the PartList rather than hacking in like this.
      oldParts[i] = ctm->getAssignedProcForTask(oldParts[i]);  
    }
#endif 

    //for now just delete the object.
    delete ctm;
  }
#endif

  else if (mapping_type == 1){
    //if mapping is 1 -- graph mapping
  }

  this->env_->debug(DETAILED_STATUS, "Exiting solve");
}

template <typename Adapter>
void PartitioningProblem<Adapter>::createPartitioningProblem(bool newData)
{
  this->env_->debug(DETAILED_STATUS,
    "PartitioningProblem::createPartitioningProblem");

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

  // Did the user specify a computational model?

  std::string model(defString);
  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("model");
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
    }
    else if (algorithm == std::string("zoltan") ||
             algorithm == std::string("parma") ||
             algorithm == std::string("forTestingOnly"))
    {
      algName_ = algorithm;
    }
    else if (algorithm == std::string("rcb") ||
             algorithm == std::string("rib") ||
             algorithm == std::string("hsfc"))
    {
      // rcb, rib, hsfc provided through Zoltan
      Teuchos::ParameterList &zparams = pl.sublist("zoltan_parameters",false);
      zparams.set("LB_METHOD", algorithm);
      if (numberOfWeights_ > 0) {
        char strval[20];
        sprintf(strval, "%d", numberOfWeights_);
        zparams.set("OBJ_WEIGHT_DIM", strval);
      }
      algName_ = std::string("zoltan");
    }
    else if (algorithm == std::string("multijagged"))
    {
      //modelType_ = CoordinateModelType;
      modelAvail_[CoordinateModelType]=true;

      algName_ = algorithm;
    }
    else if (algorithm == std::string("metis") ||
             algorithm == std::string("parmetis"))
    {

      //modelType_ = GraphModelType;
      modelAvail_[GraphModelType]=true;
      algName_ = algorithm;
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
    }
    else if (algorithm == std::string("scotch") ||
             algorithm == std::string("ptscotch")) // BDD: Don't construct graph for scotch here
    {
      algName_ = algorithm;
    }
    else if (algorithm == std::string("pulp"))
    {
      algName_ = algorithm;
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

      algName_ = std::string("phg");
    }
    else if (model == std::string("graph"))
    {
      //modelType_ = GraphModelType;
      modelAvail_[GraphModelType]=true;

#ifdef HAVE_ZOLTAN2_SCOTCH
      modelAvail_[GraphModelType]=false; // graph constructed by AlgPTScotch
      if (this->comm_->getSize() > 1)
        algName_ = std::string("ptscotch");
      else
        algName_ = std::string("scotch");
#else
#ifdef HAVE_ZOLTAN2_PARMETIS
      if (this->comm_->getSize() > 1)
        algName_ = std::string("parmetis");
      else
        algName_ = std::string("metis");
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#else
#ifdef HAVE_ZOLTAN2_PULP
      // TODO: XtraPuLP
      //if (this->comm_->getSize() > 1)
      //  algName_ = std::string("xtrapulp");
      //else
      algName_ = std::string("pulp");
#else
      algName_ = std::string("phg");
#endif
#endif
#endif
    }
    else if (model == std::string("geometry"))
    {
      //modelType_ = CoordinateModelType;
      modelAvail_[CoordinateModelType]=true;

      algName_ = std::string("multijagged");
    }
    else if (model == std::string("ids"))
    {
      //modelType_ = IdentifierModelType;
      modelAvail_[IdentifierModelType]=true;

      algName_ = std::string("block");
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

      algName_ = std::string("phg");
    }
    else if (inputType_ == GraphAdapterType ||
        inputType_ == MeshAdapterType)
    {
      //modelType_ = GraphModelType;
      modelAvail_[GraphModelType]=true;

      algName_ = std::string("phg");
    }
    else if (inputType_ == VectorAdapterType)
    {
      //modelType_ = CoordinateModelType;
      modelAvail_[CoordinateModelType]=true;

      algName_ = std::string("multijagged");
    }
    else if (inputType_ == IdentifierAdapterType)
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

    bool sgParameter = false;
    pe = pl.getEntryPtr("subset_graph");
    if (pe)
      sgParameter = pe->getValue(&sgParameter);

    if (sgParameter == 1)
        graphFlags_.set(BUILD_SUBSET_GRAPH);

    // Any special behaviors required by the algorithm?

    if (removeSelfEdges)
      graphFlags_.set(REMOVE_SELF_EDGES);

    if (needConsecutiveGlobalIds)
      graphFlags_.set(GENERATE_CONSECUTIVE_IDS);

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

  }
  //  else if (modelType_ == CoordinateModelType)
  if (modelAvail_[CoordinateModelType]==true)
  {

    // Any special behaviors required by the algorithm?

  }


  if (newData ||
      (modelAvail_[GraphModelType]!=prevModelAvail[GraphModelType]) ||
      (modelAvail_[HypergraphModelType]!=prevModelAvail[HypergraphModelType])||
      (modelAvail_[CoordinateModelType]!=prevModelAvail[CoordinateModelType])||
      (modelAvail_[IdentifierModelType]!=prevModelAvail[IdentifierModelType])||
      // (modelType_ != previousModel) ||
      (graphFlags_ != previousGraphModelFlags) ||
      (coordFlags_ != previousCoordinateModelFlags) ||
      (idFlags_    != previousIdentifierModelFlags))
  {
    // Create the computational model.
    // Models are instantiated for base input adapter types (mesh,
    // matrix, graph, and so on).  We pass a pointer to the input
    // adapter, cast as the base input type.

    //KDD Not sure why this shadow declaration is needed
    //KDD Comment out for now; revisit later if problems.
    //KDD const Teuchos::ParameterList pl = this->envConst_->getParameters();
    //bool exceptionThrow = true;

    if(modelAvail_[GraphModelType]==true)
    {
      this->env_->debug(DETAILED_STATUS, "    building graph model");
      this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
            this->baseInputAdapter_, this->envConst_, this->comm_,
            graphFlags_));

      this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
            this->graphModel_);
    } 
    if(modelAvail_[HypergraphModelType]==true)
    {
      //KDD USING ZOLTAN FOR HYPERGRAPH FOR NOW
      //KDD std::cout << "Hypergraph model not implemented yet..." << std::endl;
    }

    if(modelAvail_[CoordinateModelType]==true)
    {
      this->env_->debug(DETAILED_STATUS, "    building coordinate model");
      this->coordinateModel_ = rcp(new CoordinateModel<base_adapter_t>(
            this->baseInputAdapter_, this->envConst_, this->comm_,
            coordFlags_));

      this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
            this->coordinateModel_);
    }

    if(modelAvail_[IdentifierModelType]==true)
    {
      this->env_->debug(DETAILED_STATUS, "    building identifier model");
      this->identifierModel_ = rcp(new IdentifierModel<base_adapter_t>(
            this->baseInputAdapter_, this->envConst_, this->comm_,
            idFlags_));

      this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
            this->identifierModel_);
    }

    this->env_->memory("After creating Model");
    this->env_->debug(DETAILED_STATUS, "createPartitioningProblem done");
  }
}

}  // namespace Zoltan2
#endif
