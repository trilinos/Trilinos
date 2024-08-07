// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_PartitioningSolution.hpp
    \brief Defines the PartitioningSolution class.
*/

#ifndef _ZOLTAN2_PARTITIONINGSOLUTION_HPP_
#define _ZOLTAN2_PARTITIONINGSOLUTION_HPP_

namespace Zoltan2 {
template <typename Adapter>
class PartitioningSolution;
}

#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Solution.hpp>
#include <Zoltan2_GreedyMWM.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_CoordinatePartitioningGraph.hpp>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <sstream>

#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif


namespace Zoltan2 {

/*! \brief A PartitioningSolution is a solution to a partitioning problem.

    It is initialized by a PartitioningProblem,
    written to by an algorithm, and may be read by the user or by
    a data migration routine in an input adapter.

    \todo Problem computes metrics using the Solution.  Should
  Solution have a pointer to the metrics, since it may persist after
  the Problem is gone?
    \todo save an RCB tree, so it can be used in repartitioning, and
                supplied to the caller.
    \todo doxyfy the comments in this file.
*/

template <typename Adapter>
  class PartitioningSolution : public Solution
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
#endif

/*! \brief Constructor when part sizes are not supplied.
 *
 *   The Solution constructor may require global communication.
 *   The rest of the Solution methods do not.
 *
 *    \param env the environment for the application
 *    \param comm the communicator for the problem associated with
 *             this solution
 *    \param nUserWeights  the number of weights supplied by the
 *         application for each object.
 *    \param algorithm  Algorithm, if any, used to compute the solution.
 *
 *   It is possible that part sizes were supplied on other processes,
 *   so this constructor does do a check to see if part sizes need
 *   to be globally calculated.
 */

  PartitioningSolution( const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    int nUserWeights, 
    const RCP<Algorithm<Adapter> > &algorithm = Teuchos::null);

/*! \brief Constructor when part sizes are supplied.
 *
 *   The Solution constructor may require global communication.
 *   The rest of the Solution methods do not.
 *
 *    \param env the environment for the application
 *    \param comm the communicator for the problem associated with
 *                        this solution
 *    \param nUserWeights  the number of weights supplied
 *                         by the application
 *    \param reqPartIds  reqPartIds[i] is a list of
 *          of part numbers for weight index i.
 *    \param reqPartSizes  reqPartSizes[i] is the list
 *          of part sizes for weight i corresponding to parts in
 *          reqPartIds[i]
 *    \param algorithm  Algorithm, if any, used to compute the solution.
 *
 *   If <tt>reqPartIds[i].size()</tt> and <tt>reqPartSizes[i].size()</tt>
 *           are zero for
 *   all processes, it is assumed that part sizes for weight
 *   dimension "i" are uniform.
 *
 *   If across the application there are some part numbers that are not
 *   included in the reqPartIds lists, then those part sizes are assumed
 *   to be 1.0.
 *
 *   \todo handle errors that may arise - like duplicate part numbers
 */

  PartitioningSolution(const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    int nUserWeights, ArrayView<ArrayRCP<part_t> > reqPartIds,
    ArrayView<ArrayRCP<scalar_t> > reqPartSizes,
    const RCP<Algorithm<Adapter> > &algorithm = Teuchos::null);

  ////////////////////////////////////////////////////////////////////
  // Information that the algorithm may wish to query.

/*! \brief Returns the global number of parts desired in the solution.
 */
  size_t getTargetGlobalNumberOfParts() const { return nGlobalParts_; }

/*! \brief Returns the actual global number of parts provided in setParts().
 */
  size_t getActualGlobalNumberOfParts() const { return nGlobalPartsSolution_; }

/*! \brief Returns the number of parts to be assigned to this process.
 */
  size_t getLocalNumberOfParts() const { return nLocalParts_; }

/*! \brief If parts are divided across processes, return the fraction of
             a part on this process.
    \return zero if parts are not split across processes, approximate
                fraction of the part otherwise.
    \todo More useful to get number of processes owning part?  Or
              not useful at all - remove this?
 */
  scalar_t getLocalFractionOfPart() const { return localFraction_; }

/*! \brief Is the part-to-process distribution is one-to-one.
     \return true if Process p owns part p for all p, and false if the part
                  to process distribution is more complex.

   If this is true, then getPartDistribution() and getProcDistribution()
   return NULL pointers.  If either of the latter two methods is non-NULL,
   then this method returns false.
 */

  bool oneToOnePartDistribution() const { return onePartPerProc_; }

/*! \brief Return a distribution by part.

    \return If any parts are
       divided across processes, then a mapping \c A is returned.
     \c A such that \c A[i] is the lowest numbered process
       owning part \c i.  The length of the array is one greater than the
       global number of parts.
       The value of the last element is the global number of processes.

     Parts are divided across processes only if there are fewer parts
     than processes and the caller did not define "num_local_parts" for
     each process.  In this case, parts are divided somewhat evenly
     across the processes.  This situation is more likely to arise in
     Zoltan2 algorithms than in user applications.

   If either oneToOnePartDistribution() is true or getProcDistribution() is
   non-NULL, then this method returns a NULL pointer.  The three are mutually
   exclusive and collective exhaustive.
 */
  const int *getPartDistribution() const {
    if (partDist_.size() > 0) return &partDist_[0];
    else return NULL;
  }

/*! \brief Return a distribution by process.

    \return If the mapping of parts to processes is not one-to-one, and
      if parts are not divided across processes, then the mapping
      \c A is returned. \c A such that \c A[i] is the first part
      assigned to process \c i.
      (Parts are assigned sequentially to processes.)
      However if <tt> A[i+1]</tt> is equal to \c A[i], there
      are no parts assigned to process \c i.  The length of the array is
      one greater than the number of processes.  The last element of
      the array is the global number of parts.

    If the mapping is one-to-one, or if parts are divided across processes,
    then this method returns NULL pointer, and either
    oneToOnePartDistribution() or getPartDistribution() describes the mapping.
 */
  const part_t *getProcDistribution() const {
    if (procDist_.size() > 0) return &procDist_[0];
    else return NULL;
  }

/*! \brief Get the number of criteria (object weights)
    \return the number of criteria for which the solution has part sizes.
 */
  int getNumberOfCriteria() const { return nWeightsPerObj_; }


/*! \brief Determine if balancing criteria has uniform
                part sizes.  (User can specify differing part sizes.)
    \param idx   A value from 0 to one less than the number of weights per
                   object.
    \return true if part sizes are uniform for this criteria.
 */
  bool criteriaHasUniformPartSizes(int idx) const { return pSizeUniform_[idx];}

/*! \brief Get the size for a given weight index and a given part.

    \param idx   A value from 0 to one less than the number of weights per
                       object.
    \param part  A value from 0 to one less than the global number of parts
                   to be computed
    \return   The size for that part.  Part sizes for a given weight
                    dimension sum to 1.0.

      \todo It would be useful to algorithms to get the sum of
           part sizes from a to b, or the sum or a list of parts.
 */
  scalar_t getCriteriaPartSize(int idx, part_t part) const {
    if (pSizeUniform_[idx])
      return 1.0 / nGlobalParts_;
    else if (pCompactIndex_[idx].size())
      return pSize_[idx][pCompactIndex_[idx][part]];
    else
      return pSize_[idx][part];
  }

/*! \brief Return true if the two weight indices have the same
 *          part size information.

    \param c1   A value from 0 through one less than the number of weights.
    \param c2   A value from 0 through one less than the number of weights.
    \return   If weight index \c c1 and weight index \c c2 have
        the same part size information, the \c true is returned, otherwise
        \c false is returned.

    It may be a problem for some algorithms if there are multiple weight
    dimensions with differing part size constraints to be satisfied.
 */

  bool criteriaHaveSamePartSizes(int c1, int c2) const;

  ////////////////////////////////////////////////////////////////////
  // Method used by the algorithm to set results.

  /*! \brief The algorithm uses setParts to set the solution.
   *
   *   \param gnoList  A list of global numbers.
   *
   *   \param partList  The part assigned to gnoList[i] by the algorithm
   *      should be in partList[i].  The partList is allocated and written
   *      by the algorithm.
   *
   *   \param dataDidNotMove The algorithm did not change the order of the
   *      data provided by the model; that is, it did not move the data
   *      to other processes or reorganize within the process.  Thus,
   *      the gnoList and partList are ordered in the same way as the
   *      array provided by the model.  Setting this flag to true avoids
   *      processing to remap the data to the original process that provided
   *      the data.
   *
   * The global numbers supplied by the algorithm do not need to be
   * those representing the global Ids of that process.  But
   * all global numbers should be assigned a part by exactly one
   * process.
   *
   * setParts() must be called by all processes in the problem, as
   * the part for each global identifier supplied by each process
   * in its InputAdapter is found and saved in this PartitioningSolution.
   */

  void setParts(ArrayRCP<part_t> &partList);

  ////////////////////////////////////////////////////////////////////

  /*! \brief Remap a new partition for maximum overlap with an input partition.
   *
   * Assumptions for this version:
   * input part assignment == processor rank for every local object.
   * assuming nGlobalParts <= num ranks
   * TODO:  Write a version that takes the input part number as input;
   *        this change requires input parts in adapters to be provided in
   *        the Adapter.
   * TODO:  For repartitioning, compare to old remapping results; see Zoltan1.
   */

  void RemapParts();

  ////////////////////////////////////////////////////////////////////
  /* Return the weight of objects staying with a given remap.
   * If remap is NULL, compute weight of objects staying with given partition
   */
  long measure_stays(part_t *remap, int *idx, part_t *adj, long *wgt,
                     part_t nrhs, part_t nlhs)
  {
    long staying = 0;
    for (part_t i = 0; i < nrhs; i++) { 
      part_t k = (remap ? remap[i] : i);
      for (part_t j = idx[k]; j < idx[k+1]; j++) { 
        if (i == (adj[j]-nlhs)) {
          staying += wgt[j];
          break;
        }
      }
    }
    return staying;
  }

  ////////////////////////////////////////////////////////////////////
  // Results that may be queried by the user, by migration methods,
  // or by metric calculation methods.
  // We return raw pointers so users don't have to learn about our
  // pointer wrappers.

  /*! \brief Return the communicator associated with the solution.
   */
  inline const RCP<const Comm<int> > &getCommunicator() const { return comm_;}

  /*! \brief Return the environment associated with the solution.
   */
  inline const RCP<const Environment> &getEnvironment() const { return env_;}

  /*! \brief Returns the part list corresponding to the global ID list.
   */
  const part_t *getPartListView() const {
    if (parts_.size() > 0) return parts_.getRawPtr();
    else                   return NULL;
  }

  /*! \brief Returns the process list corresponding to the global ID list.
      \return The return value is a NULL pointer if part IDs are
                synonomous with process IDs.
   */
  const int *getProcListView() const {
    if (procs_.size() > 0) return procs_.getRawPtr();
    else                   return NULL;
  }

  /*! \brief calculate if partition tree is binary.
   */
  virtual bool isPartitioningTreeBinary() const
  {
    if (this->algorithm_ == Teuchos::null)
      throw std::logic_error("no partitioning algorithm has been run yet");
    return this->algorithm_->isPartitioningTreeBinary();
  }

  /*! \brief get the partition tree - fill the relevant arrays
   */
  void getPartitionTree(part_t & numTreeVerts,
                        std::vector<part_t> & permPartNums,
                        std::vector<part_t> & splitRangeBeg,
                        std::vector<part_t> & splitRangeEnd,
                        std::vector<part_t> & treeVertParents) const {

    part_t numParts = static_cast<part_t>(getTargetGlobalNumberOfParts());

    if (this->algorithm_ == Teuchos::null)
      throw std::logic_error("no partitioning algorithm has been run yet");
    this->algorithm_->getPartitionTree(
      numParts, // may want to change how this is passed through
      numTreeVerts,
      permPartNums,
      splitRangeBeg,
      splitRangeEnd,
      treeVertParents);
  }

  /*! \brief returns the part box boundary list.
   */
  std::vector<Zoltan2::coordinateModelPartBox> &
  getPartBoxesView() const
  {
    return this->algorithm_->getPartBoxesView();
  }

  //!  \brief Return the part overlapping a given point in space; 
  //          when a point lies on a part boundary, the lowest part
  //          number on that boundary is returned.
  //          Note that not all partitioning algorithms will support
  //          this method.
  //
  //   \param dim : the number of dimensions specified for the point in space
  //   \param point : the coordinates of the point in space; array of size dim
  //   \return the part number of a part overlapping the given point
  part_t pointAssign(int dim, scalar_t *point) const
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
  void boxAssign(int dim, scalar_t *lower, scalar_t *upper,
                 size_t &nPartsFound, part_t **partsFound) const
  {
    try {
      if (this->algorithm_ == Teuchos::null)
        throw std::logic_error("no partitioning algorithm has been run yet");

      this->algorithm_->boxAssign(dim, lower, upper, nPartsFound, partsFound); 
    }
    Z2_FORWARD_EXCEPTIONS
  }


  /*! \brief returns communication graph resulting from geometric partitioning.
   */
  void getCommunicationGraph(ArrayRCP <part_t> &comXAdj,
                             ArrayRCP <part_t> &comAdj) const
  {
    try {
      if (this->algorithm_ == Teuchos::null)
        throw std::logic_error("no partitioning algorithm has been run yet");

      this->algorithm_->getCommunicationGraph(this, comXAdj, comAdj);
    }
    Z2_FORWARD_EXCEPTIONS
  }

  /*! \brief Get the parts belonging to a process.
   *  \param procId a process rank
   *  \param numParts on return will be set the number of parts belonging
   *                    to the process.
   *  \param partMin on return will be set to minimum part number
   *  \param partMax on return will be set to maximum part number
   *
   * Normally \c numParts is at least one. But if there are more processes
   * than parts, one of two things can happen.  Either there are processes
   * with no parts, and so \c numParts will be zero, or a part may be
   * split across more than one process, in which \c numParts will
   * be non-zero but less than 1.
   *
   * In the latter case, \c numParts is 1.0 divided by the number of
   * processes that share the part.
   */

  void getPartsForProc(int procId, double &numParts, part_t &partMin,
    part_t &partMax) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid process id",
      procId >= 0 && procId < comm_->getSize(), BASIC_ASSERTION);

    procToPartsMap(procId, numParts, partMin, partMax);
  }

  /*! \brief Get the processes containing a part.
   *  \param partId a part number from 0 to one less than the global number
   *                of parts.
   *  \param procMin on return will be set to minimum proc number
   *  \param procMax on return will be set to maximum proc number
   *
   * Normally \c procMin and \c procMax are the same value and a part
   * is assigned to one process.  But if there are more processes than
   * parts, it's possible that a part will be divided across more than
   * one process.
   */
  void getProcsForPart(part_t partId, part_t &procMin, part_t &procMax) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid part id",
      partId >= 0 && partId < nGlobalParts_, BASIC_ASSERTION);

    partToProcsMap(partId, procMin, procMax);
  }

private:
  void partToProc(bool doCheck, bool haveNumLocalParts, bool haveNumGlobalParts,
    int numLocalParts, int numGlobalParts);

  void procToPartsMap(int procId, double &numParts, part_t &partMin,
    part_t &partMax) const;

  void partToProcsMap(part_t partId, int &procMin, int &procMax) const;

  void setPartDistribution();

  void setPartSizes(ArrayView<ArrayRCP<part_t> > reqPartIds,
    ArrayView<ArrayRCP<scalar_t> > reqPartSizes);

  void computePartSizes(int widx, ArrayView<part_t> ids,
    ArrayView<scalar_t> sizes);

  void broadcastPartSizes(int widx);


  RCP<const Environment> env_;             // has application communicator
  const RCP<const Comm<int> > comm_;       // the problem communicator

  //part box boundaries as a result of geometric partitioning algorithm.
  RCP<std::vector<Zoltan2::coordinateModelPartBox> > partBoxes;

  part_t nGlobalParts_;// target global number of parts
  part_t nLocalParts_; // number of parts to be on this process

  scalar_t localFraction_; // approx fraction of a part on this process
  int nWeightsPerObj_;      // if user has no weights, this is 1  TODO:  WHY???

  // If process p is to be assigned part p for all p, then onePartPerProc_
  // is true. Otherwise it is false, and either procDist_ or partDist_
  // describes the allocation of parts to processes.
  //
  // If parts are never split across processes, then procDist_ is defined
  // as follows:
  //
  //   partId              = procDist_[procId]
  //   partIdNext          = procDist_[procId+1]
  //   globalNumberOfParts = procDist_[numProcs]
  //
  // meaning that the parts assigned to process procId range from
  // [partId, partIdNext).  If partIdNext is the same as partId, then
  // process procId has no parts.
  //
  // If the number parts is less than the number of processes, and the
  // user did not specify "num_local_parts" for each of the processes, then
  // parts are split across processes, and partDist_ is defined rather than
  // procDist_.
  //
  //   procId              = partDist_[partId]
  //   procIdNext          = partDist_[partId+1]
  //   globalNumberOfProcs = partDist_[numParts]
  //
  // which implies that the part partId is shared by processes in the
  // the range [procId, procIdNext).
  //
  // We use std::vector so we can use upper_bound algorithm

  bool             onePartPerProc_;   // either this is true...
  std::vector<int>      partDist_;      // or this is defined ...
  std::vector<part_t> procDist_;      // or this is defined.
  bool procDistEquallySpread_;        // if procDist_ is used and
                                      // #parts > #procs and
                                      // num_local_parts is not specified,
                                      // parts are evenly distributed to procs

  // In order to minimize the storage required for part sizes, we
  // have three different representations.
  //
  // If the part sizes for weight index w are all the same, then:
  //    pSizeUniform_[w] = true
  //    pCompactIndex_[w].size() = 0
  //    pSize_[w].size() = 0
  //
  // and the size for part p is 1.0 / nparts.
  //
  // If part sizes differ for each part in weight index w, but there
  // are no more than 64 distinct sizes:
  //    pSizeUniform_[w] = false
  //    pCompactIndex_[w].size() = number of parts
  //    pSize_[w].size() = number of different sizes
  //
  // and the size for part p is pSize_[pCompactIndex_[p]].
  //
  // If part sizes differ for each part in weight index w, and there
  // are more than 64 distinct sizes:
  //    pSizeUniform_[w] = false
  //    pCompactIndex_[w].size() = 0
  //    pSize_[w].size() = nparts
  //
  // and the size for part p is pSize_[p].
  //
  // NOTE: If we expect to have similar cases, i.e. a long list of scalars
  //   where it is highly possible that just a few unique values appear,
  //   then we may want to make this a class.  The advantage is that we
  //   save a long list of 1-byte indices instead of a long list of scalars.

  ArrayRCP<bool> pSizeUniform_;
  ArrayRCP<ArrayRCP<unsigned char> > pCompactIndex_;
  ArrayRCP<ArrayRCP<scalar_t> > pSize_;

  ////////////////////////////////////////////////////////////////
  // The algorithm sets these values upon completion.

  ArrayRCP<part_t> parts_;      // part number assigned to localid[i]

  bool haveSolution_;

  part_t nGlobalPartsSolution_; // global number of parts in solution

  ////////////////////////////////////////////////////////////////
  // The solution calculates this from the part assignments,
  // unless onePartPerProc_.

  ArrayRCP<int> procs_;       // process rank assigned to localid[i]

  ////////////////////////////////////////////////////////////////
  // Algorithm used to compute the solution; 
  // needed for post-processing with pointAssign or getCommunicationGraph
  const RCP<Algorithm<Adapter> > algorithm_;  // 
};

////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////

template <typename Adapter>
  PartitioningSolution<Adapter>::PartitioningSolution(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    int nUserWeights,
    const RCP<Algorithm<Adapter> > &algorithm)
    : env_(env), comm_(comm),
      partBoxes(),
      nGlobalParts_(0), nLocalParts_(0),
      localFraction_(0),  nWeightsPerObj_(),
      onePartPerProc_(false), partDist_(), procDist_(),
      procDistEquallySpread_(false),
      pSizeUniform_(), pCompactIndex_(), pSize_(),
      parts_(), haveSolution_(false), nGlobalPartsSolution_(0),
      procs_(), algorithm_(algorithm)
{
  nWeightsPerObj_ = (nUserWeights ? nUserWeights : 1);  // TODO:  WHY??  WHY NOT ZERO?

  setPartDistribution();

  // We must call setPartSizes() because part sizes may have
  // been provided by the user on other processes.

  ArrayRCP<part_t> *noIds = new ArrayRCP<part_t> [nWeightsPerObj_];
  ArrayRCP<scalar_t> *noSizes = new ArrayRCP<scalar_t> [nWeightsPerObj_];
  ArrayRCP<ArrayRCP<part_t> > ids(noIds, 0, nWeightsPerObj_, true);
  ArrayRCP<ArrayRCP<scalar_t> > sizes(noSizes, 0, nWeightsPerObj_, true);

  setPartSizes(ids.view(0, nWeightsPerObj_), sizes.view(0, nWeightsPerObj_));

  env_->memory("After construction of solution");
}

template <typename Adapter>
  PartitioningSolution<Adapter>::PartitioningSolution(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    int nUserWeights,
    ArrayView<ArrayRCP<part_t> > reqPartIds,
    ArrayView<ArrayRCP<scalar_t> > reqPartSizes,
    const RCP<Algorithm<Adapter> > &algorithm)
    : env_(env), comm_(comm),
      partBoxes(),
      nGlobalParts_(0), nLocalParts_(0),
      localFraction_(0),  nWeightsPerObj_(),
      onePartPerProc_(false), partDist_(), procDist_(),
      procDistEquallySpread_(false),
      pSizeUniform_(), pCompactIndex_(), pSize_(),
      parts_(), haveSolution_(false), nGlobalPartsSolution_(0),
      procs_(), algorithm_(algorithm)
{
  nWeightsPerObj_ = (nUserWeights ? nUserWeights : 1);  // TODO:  WHY?? WHY NOT ZERO?

  setPartDistribution();

  setPartSizes(reqPartIds, reqPartSizes);

  env_->memory("After construction of solution");
}

template <typename Adapter>
  void PartitioningSolution<Adapter>::setPartDistribution()
{
  // Did the caller define num_global_parts and/or num_local_parts?

  const ParameterList &pl = env_->getParameters();
  size_t haveGlobalNumParts=0, haveLocalNumParts=0;
  int numLocal=0, numGlobal=0;

  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("num_global_parts");

  if (pe){
    haveGlobalNumParts = 1;
    nGlobalParts_ = part_t(pe->getValue(&nGlobalParts_));
    numGlobal = nGlobalParts_;
  }

  pe = pl.getEntryPtr("num_local_parts");

  if (pe){
    haveLocalNumParts = 1;
    nLocalParts_ = part_t(pe->getValue(&nLocalParts_));
    numLocal = nLocalParts_;
  }

  try{
    // Sets onePartPerProc_, partDist_, and procDist_

    partToProc(true, haveLocalNumParts, haveGlobalNumParts,
      numLocal, numGlobal);
  }
  Z2_FORWARD_EXCEPTIONS

  int nprocs = comm_->getSize();
  int rank = comm_->getRank();

  if (onePartPerProc_){
    nGlobalParts_ = nprocs;
    nLocalParts_ = 1;
  }
  else if (partDist_.size() > 0){   // more procs than parts
    nGlobalParts_ = partDist_.size() - 1;
    int pstart = partDist_[0];
    for (part_t i=1; i <= nGlobalParts_; i++){
      int pend = partDist_[i];
      if (rank >= pstart && rank < pend){
        int numOwners = pend - pstart;
        nLocalParts_ = 1;
        localFraction_ = 1.0 / numOwners;
        break;
      }
      pstart = pend;
    }
  }
  else if (procDist_.size() > 0){  // more parts than procs
    nGlobalParts_ = procDist_[nprocs];
    nLocalParts_ = procDist_[rank+1] - procDist_[rank];
  }
  else {
    throw std::logic_error("partToProc error");
  }

}

template <typename Adapter>
  void PartitioningSolution<Adapter>::setPartSizes(
    ArrayView<ArrayRCP<part_t> > ids, ArrayView<ArrayRCP<scalar_t> > sizes)
{
  int widx = nWeightsPerObj_;
  bool fail=false;

  size_t *countBuf = new size_t [widx*2];
  ArrayRCP<size_t> counts(countBuf, 0, widx*2, true);

  fail = ((ids.size() != widx) || (sizes.size() != widx));

  for (int w=0; !fail && w < widx; w++){
    counts[w] = ids[w].size();
    if (ids[w].size() != sizes[w].size()) fail=true;
  }

  env_->globalBugAssertion(__FILE__, __LINE__, "bad argument arrays", fail==0,
    COMPLEX_ASSERTION, comm_);

  // Are all part sizes the same?  This is the common case.

  ArrayRCP<scalar_t> *emptySizes= new ArrayRCP<scalar_t> [widx];
  pSize_ = arcp(emptySizes, 0, widx);

  ArrayRCP<unsigned char> *emptyIndices= new ArrayRCP<unsigned char> [widx];
  pCompactIndex_ = arcp(emptyIndices, 0, widx);

  bool *info = new bool [widx];
  pSizeUniform_ = arcp(info, 0, widx);
  for (int w=0; w < widx; w++)
    pSizeUniform_[w] = true;

  if (nGlobalParts_ == 1){
    return;   // there's only one part in the whole problem
  }

  size_t *ptr1 = counts.getRawPtr();
  size_t *ptr2 = counts.getRawPtr() + widx;

  try{
    reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_MAX, widx, ptr1, ptr2);
  }
  Z2_THROW_OUTSIDE_ERROR(*env_);

  bool zero = true;

  for (int w=0; w < widx; w++)
    if (counts[widx+w] > 0){
      zero = false;
      pSizeUniform_[w] = false;
    }

  if (zero) // Part sizes for all criteria are uniform.
    return;

  // Compute the part sizes for criteria for which part sizes were
  // supplied.  Normalize for each criteria so part sizes sum to one.

  int nprocs = comm_->getSize();
  int rank = comm_->getRank();

  for (int w=0; w < widx; w++){
    if (pSizeUniform_[w]) continue;

    // Send all ids and sizes to one process.
    // (There is no simple gather method in Teuchos.)

    part_t length = ids[w].size();
    part_t *allLength = new part_t [nprocs];
    Teuchos::gatherAll<int, part_t>(*comm_, 1, &length,
      nprocs, allLength);

    if (rank == 0){
      int total = 0;
      for (int i=0; i < nprocs; i++)
        total += allLength[i];

      part_t *partNums = new part_t [total];
      scalar_t *partSizes = new scalar_t [total];

      ArrayView<part_t> idArray(partNums, total);
      ArrayView<scalar_t> sizeArray(partSizes, total);

      if (length > 0){
        for (int i=0; i < length; i++){
          *partNums++ = ids[w][i];
          *partSizes++ = sizes[w][i];
        }
      }

      for (int p=1; p < nprocs; p++){
        if (allLength[p] > 0){
          Teuchos::receive<int, part_t>(*comm_, p,
            allLength[p], partNums);
          Teuchos::receive<int, scalar_t>(*comm_, p,
            allLength[p], partSizes);
          partNums += allLength[p];
          partSizes += allLength[p];
        }
      }

      delete [] allLength;

      try{
        computePartSizes(w, idArray, sizeArray);
      }
      Z2_FORWARD_EXCEPTIONS

      delete [] idArray.getRawPtr();
      delete [] sizeArray.getRawPtr();
    }
    else{
      delete [] allLength;
      if (length > 0){
        Teuchos::send<int, part_t>(*comm_, length, ids[w].getRawPtr(), 0);
        Teuchos::send<int, scalar_t>(*comm_, length, sizes[w].getRawPtr(), 0);
      }
    }

    broadcastPartSizes(w);
  }
}

template <typename Adapter>
  void PartitioningSolution<Adapter>::broadcastPartSizes(int widx)
{
  env_->localBugAssertion(__FILE__, __LINE__, "preallocations",
    pSize_.size()>widx &&
    pSizeUniform_.size()>widx && pCompactIndex_.size()>widx,
    COMPLEX_ASSERTION);

  int rank = comm_->getRank();
  int nprocs = comm_->getSize();
  part_t nparts = nGlobalParts_;

  if (nprocs < 2)
    return;

  char flag=0;

  if (rank == 0){
    if (pSizeUniform_[widx] == true)
      flag = 1;
    else if (pCompactIndex_[widx].size() > 0)
      flag = 2;
    else
      flag = 3;
  }

  try{
    Teuchos::broadcast<int, char>(*comm_, 0, 1, &flag);
  }
  Z2_THROW_OUTSIDE_ERROR(*env_);

  if (flag == 1){
    if (rank > 0)
      pSizeUniform_[widx] = true;

    return;
  }

  if (flag == 2){

    // broadcast the indices into the size list

    unsigned char *idxbuf = NULL;

    if (rank > 0){
      idxbuf = new unsigned char [nparts];
      env_->localMemoryAssertion(__FILE__, __LINE__, nparts, idxbuf);
    }
    else{
      env_->localBugAssertion(__FILE__, __LINE__, "index list size",
        pCompactIndex_[widx].size() == nparts, COMPLEX_ASSERTION);
      idxbuf = pCompactIndex_[widx].getRawPtr();
    }

    try{
      // broadcast of unsigned char is not supported
      Teuchos::broadcast<int, char>(*comm_, 0, nparts,
        reinterpret_cast<char *>(idxbuf));
    }
    Z2_THROW_OUTSIDE_ERROR(*env_);

    if (rank > 0)
      pCompactIndex_[widx] = arcp(idxbuf, 0, nparts, true);

    // broadcast the list of different part sizes

    unsigned char maxIdx=0;
    for (part_t p=0; p < nparts; p++)
      if (idxbuf[p] > maxIdx) maxIdx = idxbuf[p];

    int numSizes = maxIdx + 1;

    scalar_t *sizeList = NULL;

    if (rank > 0){
      sizeList = new scalar_t [numSizes];
      env_->localMemoryAssertion(__FILE__, __LINE__, numSizes, sizeList);
    }
    else{
      env_->localBugAssertion(__FILE__, __LINE__, "wrong number of sizes",
        numSizes == pSize_[widx].size(), COMPLEX_ASSERTION);

      sizeList = pSize_[widx].getRawPtr();
    }

    try{
      Teuchos::broadcast<int, scalar_t>(*comm_, 0, numSizes, sizeList);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_);

    if (rank > 0)
      pSize_[widx] = arcp(sizeList, 0, numSizes, true);

    return;
  }

  if (flag == 3){

    // broadcast the size of each part

    scalar_t *sizeList = NULL;

    if (rank > 0){
      sizeList = new scalar_t [nparts];
      env_->localMemoryAssertion(__FILE__, __LINE__, nparts, sizeList);
    }
    else{
      env_->localBugAssertion(__FILE__, __LINE__, "wrong number of sizes",
        nparts == pSize_[widx].size(), COMPLEX_ASSERTION);

      sizeList = pSize_[widx].getRawPtr();
    }

    try{
      Teuchos::broadcast<int, scalar_t >(*comm_, 0, nparts, sizeList);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_);

    if (rank > 0)
      pSize_[widx] = arcp(sizeList, 0, nparts);

    return;
  }
}

template <typename Adapter>
  void PartitioningSolution<Adapter>::computePartSizes(int widx,
    ArrayView<part_t> ids, ArrayView<scalar_t> sizes)
{
  int len = ids.size();

  if (len == 0){
    pSizeUniform_[widx] = true;
    return;
  }

  env_->localBugAssertion(__FILE__, __LINE__, "bad array sizes",
    len>0 && sizes.size()==len, COMPLEX_ASSERTION);

  env_->localBugAssertion(__FILE__, __LINE__, "bad index",
    widx>=0 && widx<nWeightsPerObj_, COMPLEX_ASSERTION);

  env_->localBugAssertion(__FILE__, __LINE__, "preallocations",
    pSize_.size()>widx &&
    pSizeUniform_.size()>widx && pCompactIndex_.size()>widx,
    COMPLEX_ASSERTION);

  // Check ids and sizes and find min, max and average sizes.
  // If sizes are very close to uniform, call them uniform parts.

  part_t nparts = nGlobalParts_;
  unsigned char *buf = new unsigned char [nparts];
  env_->localMemoryAssertion(__FILE__, __LINE__, nparts, buf);
  memset(buf, 0, nparts);
  ArrayRCP<unsigned char> partIdx(buf, 0, nparts, true);

  scalar_t epsilon = 10e-5 / nparts;
  scalar_t min=sizes[0], max=sizes[0], sum=0;

  for (int i=0; i < len; i++){
    part_t id = ids[i];
    scalar_t size = sizes[i];

    env_->localInputAssertion(__FILE__, __LINE__, "invalid part id",
      id>=0 && id<nparts, BASIC_ASSERTION);

    env_->localInputAssertion(__FILE__, __LINE__, "invalid part size", size>=0,
      BASIC_ASSERTION);

    // TODO: we could allow users to specify multiple sizes for the same
    // part if we add a parameter that says what we are to do with them:
    // add them or take the max.

    env_->localInputAssertion(__FILE__, __LINE__,
      "multiple sizes provided for one part", partIdx[id]==0, BASIC_ASSERTION);

    partIdx[id] = 1;    // mark that we have a size for this part

    if (size < min) min = size;
    if (size > max) max = size;
    sum += size;
  }

  if (sum == 0){

    // User has given us a list of parts of size 0, we'll set
    // the rest to them to equally sized parts.

    scalar_t *allSizes = new scalar_t [2];
    env_->localMemoryAssertion(__FILE__, __LINE__, 2, allSizes);

    ArrayRCP<scalar_t> sizeArray(allSizes, 0, 2, true);

    int numNonZero = nparts - len;

    allSizes[0] = 0.0;
    allSizes[1] = 1.0 / numNonZero;

    for (part_t p=0; p < nparts; p++)
      buf[p] = 1;                 // index to default part size

    for (int i=0; i < len; i++)
      buf[ids[i]] = 0;            // index to part size zero

    pSize_[widx] = sizeArray;
    pCompactIndex_[widx] = partIdx;

    return;
  }

  if (max - min <= epsilon){
    pSizeUniform_[widx] = true;
    return;
  }

  // A size for parts that were not specified:
  scalar_t avg = sum / nparts;

  // We are going to merge part sizes that are very close.  This takes
  // computation time now, but can save considerably in the storage of
  // all part sizes on each process.  For example, a common case may
  // be some parts are size 1 and all the rest are size 2.

  scalar_t *tmp = new scalar_t [len];
  env_->localMemoryAssertion(__FILE__, __LINE__, len, tmp);
  memcpy(tmp, sizes.getRawPtr(), sizeof(scalar_t) * len);
  ArrayRCP<scalar_t> partSizes(tmp, 0, len, true);

  std::sort(partSizes.begin(), partSizes.end());

  // create a list of sizes that are unique within epsilon

  Array<scalar_t> nextUniqueSize;
  nextUniqueSize.push_back(partSizes[len-1]);   // largest
  scalar_t curr = partSizes[len-1];
  int avgIndex = len;
  bool haveAvg = false;
  if (curr - avg <= epsilon)
     avgIndex = 0;

  for (int i=len-2; i >= 0; i--){
    scalar_t val = partSizes[i];
    if (curr - val > epsilon){
      nextUniqueSize.push_back(val);  // the highest in the group
      curr = val;
      if (avgIndex==len && val > avg && val - avg <= epsilon){
        // the average would be in this group
        avgIndex = nextUniqueSize.size() - 1;
        haveAvg = true;
      }
    }
  }

  partSizes.clear();

  size_t numSizes = nextUniqueSize.size();
  int sizeArrayLen = numSizes;

  if (numSizes < 64){
    // We can store the size for every part in a compact way.

    // Create a list of all sizes in increasing order

    if (!haveAvg) sizeArrayLen++;   // need to include average

    scalar_t *allSizes = new scalar_t [sizeArrayLen];
    env_->localMemoryAssertion(__FILE__, __LINE__, sizeArrayLen, allSizes);
    ArrayRCP<scalar_t> sizeArray(allSizes, 0, sizeArrayLen, true);

    int newAvgIndex = sizeArrayLen;

    for (int i=numSizes-1, idx=0; i >= 0; i--){

      if (newAvgIndex == sizeArrayLen){

        if (haveAvg && i==avgIndex)
          newAvgIndex = idx;

        else if (!haveAvg && avg < nextUniqueSize[i]){
          newAvgIndex = idx;
          allSizes[idx++] = avg;
        }
      }

      allSizes[idx++] = nextUniqueSize[i];
    }

    env_->localBugAssertion(__FILE__, __LINE__, "finding average in list",
      newAvgIndex < sizeArrayLen, COMPLEX_ASSERTION);

    for (int i=0; i < nparts; i++){
      buf[i] = newAvgIndex;   // index to default part size
    }

    sum = (nparts - len) * allSizes[newAvgIndex];

    for (int i=0; i < len; i++){
      int id = ids[i];
      scalar_t size = sizes[i];
      int index;

      // Find the first size greater than or equal to this size.

      if (size < avg && avg - size <= epsilon)
        index = newAvgIndex;
      else{
        typename ArrayRCP<scalar_t>::iterator found =
          std::lower_bound(sizeArray.begin(), sizeArray.end(), size);

        env_->localBugAssertion(__FILE__, __LINE__, "size array",
          found != sizeArray.end(), COMPLEX_ASSERTION);

        index = found - sizeArray.begin();
      }

      buf[id] = index;
      sum += allSizes[index];
    }

    for (int i=0; i < sizeArrayLen; i++){
      sizeArray[i] /= sum;
    }

    pCompactIndex_[widx] = partIdx;
    pSize_[widx] = sizeArray;
  }
  else{
    // To have access to part sizes, we must store nparts scalar_ts on
    // every process.  We expect this is a rare case.

    tmp = new scalar_t [nparts];
    env_->localMemoryAssertion(__FILE__, __LINE__, nparts, tmp);

    sum += ((nparts - len) * avg);

    for (int i=0; i < nparts; i++){
      tmp[i] = avg/sum;
    }

    for (int i=0; i < len; i++){
      tmp[ids[i]] = sizes[i]/sum;
    }

    pSize_[widx] = arcp(tmp, 0, nparts);
  }
}

template <typename Adapter>
  void PartitioningSolution<Adapter>::setParts(ArrayRCP<part_t> &partList)
{
  env_->debug(DETAILED_STATUS, "Entering setParts");

  size_t len = partList.size();

  // Find the actual number of parts in the solution, which can
  // be more or less than the nGlobalParts_ target.
  // (We may want to compute the imbalance of a given solution with
  // respect to a desired solution.  This solution may have more or
  // fewer parts that the desired solution.)

  part_t lMax = -1;
  part_t lMin = (len > 0 ? std::numeric_limits<part_t>::max() : 0);
  part_t gMax, gMin;

  for (size_t i = 0; i < len; i++) {
    if (partList[i] < lMin) lMin = partList[i];
    if (partList[i] > lMax) lMax = partList[i];
  }
  Teuchos::reduceAll<int, part_t>(*comm_, Teuchos::REDUCE_MIN, 1, &lMin, &gMin);
  Teuchos::reduceAll<int, part_t>(*comm_, Teuchos::REDUCE_MAX, 1, &lMax, &gMax);

  nGlobalPartsSolution_ = gMax - gMin + 1;
  parts_ = partList;

  // Now determine which process gets each object, if not one-to-one.

  if (!onePartPerProc_){

    int *procs = new int [len];
    env_->localMemoryAssertion(__FILE__, __LINE__, len, procs);
    procs_ = arcp<int>(procs, 0, len);

    if (len > 0) {
      part_t *parts = partList.getRawPtr();
  
      if (procDist_.size() > 0){    // parts are not split across procs
  
        int procId;
        for (size_t i=0; i < len; i++){
          partToProcsMap(parts[i], procs[i], procId);
        }
      }
      else{  // harder - we need to split the parts across multiple procs
  
        lno_t *partCounter = new lno_t [nGlobalPartsSolution_];
        env_->localMemoryAssertion(__FILE__, __LINE__, nGlobalPartsSolution_,
          partCounter);
  
        int numProcs = comm_->getSize();
  
        //MD NOTE: there was no initialization for partCounter.
        //I added the line below, correct me if I am wrong.
        memset(partCounter, 0, sizeof(lno_t) * nGlobalPartsSolution_);
  
        for (typename ArrayRCP<part_t>::size_type i=0; i < partList.size(); i++)
          partCounter[parts[i]]++;
  
        lno_t *procCounter = new lno_t [numProcs];
        env_->localMemoryAssertion(__FILE__, __LINE__, numProcs, procCounter);
  
        int proc1;
        int proc2 = partDist_[0];
  
        for (part_t part=1; part < nGlobalParts_; part++){
          proc1 = proc2;
          proc2 = partDist_[part+1];
          int numprocs = proc2 - proc1;
  
          double dNum = partCounter[part];
          double dProcs = numprocs;
  
          //std::cout << "dNum:" << dNum << " dProcs:" << dProcs << std::endl;
          double each = floor(dNum/dProcs);
          double extra = fmod(dNum,dProcs);
  
          for (int proc=proc1, i=0; proc<proc2; proc++, i++){
            if (i < extra)
              procCounter[proc] = lno_t(each) + 1;
            else
              procCounter[proc] = lno_t(each);
          }
        }
  
        delete [] partCounter;
  
        for (typename ArrayRCP<part_t>::size_type i=0; i < partList.size(); i++){
          if (partList[i] >= nGlobalParts_){
            // Solution has more parts that targeted.  These
            // objects just remain on this process.
            procs[i] = comm_->getRank();
            continue;
          }
          part_t partNum = parts[i];
          proc1 = partDist_[partNum];
          proc2 = partDist_[partNum + 1];
  
          int proc;
          for (proc=proc1; proc < proc2; proc++){
            if (procCounter[proc] > 0){
              procs[i] = proc;
              procCounter[proc]--;
              break;
            }
          }
          env_->localBugAssertion(__FILE__, __LINE__, "part to proc",
            proc < proc2, COMPLEX_ASSERTION);
        }
  
        delete [] procCounter;
      }
    }
  }

  // Now that parts_ info is back on home process, remap the parts.
  // TODO:  The parts will be inconsistent with the proc assignments after
  // TODO:  remapping.  This problem will go away after we separate process
  // TODO:  mapping from setParts.  But for MueLu's use case, the part
  // TODO:  remapping is all that matters; they do not use the process mapping.
  bool doRemap = false;
  const Teuchos::ParameterEntry *pe =
                 env_->getParameters().getEntryPtr("remap_parts");
  if (pe) doRemap = pe->getValue(&doRemap);
  if (doRemap) RemapParts();

  haveSolution_ = true;

  env_->memory("After Solution has processed algorithm's answer");
  env_->debug(DETAILED_STATUS, "Exiting setParts");
}


template <typename Adapter>
  void PartitioningSolution<Adapter>::procToPartsMap(int procId,
    double &numParts, part_t &partMin, part_t &partMax) const
{
  if (onePartPerProc_){
    numParts = 1.0;
    partMin = partMax = procId;
  }
  else if (procDist_.size() > 0){
    partMin = procDist_[procId];
    partMax = procDist_[procId+1] - 1;
    numParts = procDist_[procId+1] - partMin;
  }
  else{
    // find the first p such that partDist_[p] > procId

    std::vector<int>::const_iterator entry;
    entry = std::upper_bound(partDist_.begin(), partDist_.end(), procId);

    size_t partIdx = entry - partDist_.begin();
    int numProcs = partDist_[partIdx] - partDist_[partIdx-1];
    partMin = partMax = int(partIdx) - 1;
    numParts = 1.0 / numProcs;
  }
}

template <typename Adapter>
  void PartitioningSolution<Adapter>::partToProcsMap(part_t partId,
    int &procMin, int &procMax) const
{
  if (partId >= nGlobalParts_){
    // setParts() may be given an initial solution which uses a
    // different number of parts than the desired solution.  It is
    // still a solution.  We keep it on this process.
    procMin = procMax = comm_->getRank();
  }
  else if (onePartPerProc_){
    procMin = procMax = int(partId);
  }
  else if (procDist_.size() > 0){
    if (procDistEquallySpread_) {
      // Avoid binary search.
      double fProcs = comm_->getSize();
      double fParts = nGlobalParts_;
      double each = fParts / fProcs;
      procMin = int(partId / each);
      while (procDist_[procMin] > partId) procMin--;
      while (procDist_[procMin+1] <= partId) procMin++;
      procMax = procMin;
    }
    else {
      // find the first p such that procDist_[p] > partId.
      // For now, do a binary search.

      typename std::vector<part_t>::const_iterator entry;
      entry = std::upper_bound(procDist_.begin(), procDist_.end(), partId);

      size_t procIdx = entry - procDist_.begin();
      procMin = procMax = int(procIdx) - 1;
    }
  }
  else{
    procMin = partDist_[partId];
    procMax = partDist_[partId+1] - 1;
  }
}

template <typename Adapter>
  bool PartitioningSolution<Adapter>::criteriaHaveSamePartSizes(
    int c1, int c2) const
{
  if (c1 < 0 || c1 >= nWeightsPerObj_ || c2 < 0 || c2 >= nWeightsPerObj_ )
    throw std::logic_error("criteriaHaveSamePartSizes error");

  bool theSame = false;

  if (c1 == c2)
    theSame = true;

  else if (pSizeUniform_[c1] == true && pSizeUniform_[c2] == true)
    theSame = true;

  else if (pCompactIndex_[c1].size() == pCompactIndex_[c2].size()){
    theSame = true;
    bool useIndex = pCompactIndex_[c1].size() > 0;
    if (useIndex){
      for (part_t p=0; theSame && p < nGlobalParts_; p++)
        if (pSize_[c1][pCompactIndex_[c1][p]] !=
            pSize_[c2][pCompactIndex_[c2][p]])
          theSame = false;
    }
    else{
      for (part_t p=0; theSame && p < nGlobalParts_; p++)
        if (pSize_[c1][p] != pSize_[c2][p])
          theSame = false;
    }
  }

  return theSame;
}

/*! \brief  Compute the assignment of parts to processes.
 *    \param doCheck  if true, do a global check to reconcile the
 *       numLocalParts and numGlobalParts values on all processes.
 *       If false, we assume numLocalParts and numGlobalParts are
 *       the same across processes.
 *   \param haveNumLocalParts true if this process set the num_local_parts
 *         parameter, false otherwise
 *   \param numLocalParts the number of local parts specified for this
 *        process as a parameter, if any
 *   \param haveNumGlobalParts true if this process set the num_global_parts
 *         parameter, false otherwise
 *   \param numGlobalParts the number of global parts specified by this
 *        process as a parameter, if any
 */

template <typename Adapter>
  void PartitioningSolution<Adapter>::partToProc(
    bool doCheck, bool haveNumLocalParts, bool haveNumGlobalParts,
    int numLocalParts, int numGlobalParts)
{
#ifdef _MSC_VER
  typedef SSIZE_T ssize_t;
#endif
  int nprocs = comm_->getSize();
  ssize_t reducevals[4];
  ssize_t sumHaveGlobal=0, sumHaveLocal=0;
  ssize_t sumGlobal=0, sumLocal=0;
  ssize_t maxGlobal=0, maxLocal=0;
  ssize_t vals[4] = {haveNumGlobalParts, haveNumLocalParts,
                     numGlobalParts, numLocalParts};

  partDist_.clear();
  procDist_.clear();

  if (doCheck){

    try{
      reduceAll<int, ssize_t>(*comm_, Teuchos::REDUCE_SUM, 4, vals, reducevals);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_);

    sumHaveGlobal = reducevals[0];
    sumHaveLocal = reducevals[1];
    sumGlobal = reducevals[2];
    sumLocal = reducevals[3];

    env_->localInputAssertion(__FILE__, __LINE__,
      "Either all procs specify num_global/local_parts or none do",
      (sumHaveGlobal == 0 || sumHaveGlobal == nprocs) &&
      (sumHaveLocal == 0 || sumHaveLocal == nprocs),
      BASIC_ASSERTION);
  }
  else{
    if (haveNumLocalParts)
      sumLocal = numLocalParts * nprocs;
    if (haveNumGlobalParts)
      sumGlobal = numGlobalParts * nprocs;

    sumHaveGlobal = haveNumGlobalParts ? nprocs : 0;
    sumHaveLocal = haveNumLocalParts ? nprocs : 0;

    maxLocal = numLocalParts;
    maxGlobal = numGlobalParts;
  }

  if (!haveNumLocalParts && !haveNumGlobalParts){
    onePartPerProc_ = true;   // default if user did not specify
    return;
  }

  if (haveNumGlobalParts){
    if (doCheck){
      vals[0] = numGlobalParts;
      vals[1] = numLocalParts;
      try{
        reduceAll<int, ssize_t>(
          *comm_, Teuchos::REDUCE_MAX, 2, vals, reducevals);
      }
      Z2_THROW_OUTSIDE_ERROR(*env_);

      maxGlobal = reducevals[0];
      maxLocal = reducevals[1];

      env_->localInputAssertion(__FILE__, __LINE__,
        "Value for num_global_parts is different on different processes.",
        maxGlobal * nprocs == sumGlobal, BASIC_ASSERTION);
    }

    if (sumLocal){
      env_->localInputAssertion(__FILE__, __LINE__,
        "Sum of num_local_parts does not equal requested num_global_parts",
        sumLocal == numGlobalParts, BASIC_ASSERTION);

      if (sumLocal == nprocs && maxLocal == 1){
        onePartPerProc_ = true;   // user specified one part per proc
        return;
      }
    }
    else{
      if (maxGlobal == nprocs){
        onePartPerProc_ = true;   // user specified num parts is num procs
        return;
      }
    }
  }

  // If we are here, we do not have #parts == #procs.

  if (sumHaveLocal == nprocs){
    //
    // We will go by the number of local parts specified.
    //

    try{
      procDist_.resize(nprocs+1);
    }
    catch (std::exception &){
      throw(std::bad_alloc());
    }

    part_t *procArray = &procDist_[0];

    try{
      part_t tmp = part_t(numLocalParts);
      gatherAll<int, part_t>(*comm_, 1, &tmp, nprocs, procArray + 1);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_);

    procArray[0] = 0;

    for (int proc=0; proc < nprocs; proc++)
      procArray[proc+1] += procArray[proc];
  }
  else{
    //
    // We will allocate global number of parts to the processes.
    //
    double fParts = numGlobalParts;
    double fProcs = nprocs;

    if (fParts < fProcs){

      try{
        partDist_.resize(size_t(fParts+1));
      }
      catch (std::exception &){
        throw(std::bad_alloc());
      }

      int *partArray = &partDist_[0];

      double each = floor(fProcs / fParts);
      double extra = fmod(fProcs, fParts);
      partDist_[0] = 0;

      for (part_t part=0; part < numGlobalParts; part++){
        int numOwners = int(each + ((part<extra) ? 1 : 0));
        partArray[part+1] = partArray[part] + numOwners;
      }

      env_->globalBugAssertion(__FILE__, __LINE__, "#parts != #procs",
        partDist_[numGlobalParts] == nprocs, COMPLEX_ASSERTION, comm_);
    }
    else if (fParts > fProcs){

      // User did not specify local number of parts per proc;
      // Distribute the parts evenly among the procs.

      procDistEquallySpread_ = true;

      try{
        procDist_.resize(size_t(fProcs+1));
      }
      catch (std::exception &){
        throw(std::bad_alloc());
      }

      part_t *procArray = &procDist_[0];

      double each = floor(fParts / fProcs);
      double extra = fmod(fParts, fProcs);
      procArray[0] = 0;

      for (int proc=0; proc < nprocs; proc++){
        part_t numParts = part_t(each + ((proc<extra) ? 1 : 0));
        procArray[proc+1] = procArray[proc] + numParts;
      }

      env_->globalBugAssertion(__FILE__, __LINE__, "#parts != #procs",
        procDist_[nprocs] == numGlobalParts, COMPLEX_ASSERTION, comm_);
    }
    else{
      env_->globalBugAssertion(__FILE__, __LINE__,
        "should never get here", 1, COMPLEX_ASSERTION, comm_);
    }
  }
}

////////////////////////////////////////////////////////////////////
// Remap a new part assignment vector for maximum overlap with an input
// part assignment.
//
// Assumptions for this version:
//   input part assignment == processor rank for every local object.
//   assuming nGlobalParts_ <= num ranks
// TODO:  Write a version that takes the input part number as input;
//        this change requires input parts in adapters to be provided in
//        the Adapter.
// TODO:  For repartitioning, compare to old remapping results; see Zoltan1.

template <typename Adapter>
void PartitioningSolution<Adapter>::RemapParts()
{
  size_t len = parts_.size();

  part_t me = comm_->getRank();
  int np = comm_->getSize();

  if (np < nGlobalParts_) {
    if (me == 0) {
      std::ostringstream msg; 
      msg << "Remapping not yet supported for "
           << "num_global_parts " << nGlobalParts_
           << " > num procs " << np << std::endl;
      env_->debug(DETAILED_STATUS, msg.str());
    }
    return;
  }
  // Build edges of a bipartite graph with np + nGlobalParts_ vertices,
  // and edges between a process vtx and any parts to which that process'
  // objects are assigned.
  // Weight edge[parts_[i]] by the number of objects that are going from
  // this rank to parts_[i].
  // We use a std::map, assuming the number of unique parts in the parts_ array
  // is small to keep the binary search efficient.
  // TODO We use the count of objects to move; should change to SIZE of objects
  // to move; need SIZE function in Adapter.

  std::map<part_t, long> edges;
  long lstaying = 0;  // Total num of local objects staying if we keep the
                      // current mapping. TODO:  change to SIZE of local objs
  long gstaying = 0;  // Total num of objects staying in the current partition

  for (size_t i = 0; i < len; i++) {
    edges[parts_[i]]++;                // TODO Use obj size instead of count
    if (parts_[i] == me) lstaying++;    // TODO Use obj size instead of count
  }

  Teuchos::reduceAll<int, long>(*comm_, Teuchos::REDUCE_SUM, 1,
                                &lstaying, &gstaying);
//TODO  if (gstaying == Adapter::getGlobalNumObjs()) return;  // Nothing to do

  part_t *remap = NULL;

  int nedges = edges.size();

  // Gather the graph to rank 0.
  part_t tnVtx = np + nGlobalParts_;  // total # vertices
  int *idx = NULL;    // Pointer index into graph adjacencies
  int *sizes = NULL;  // nedges per rank
  if (me == 0) {
    idx = new int[tnVtx+1];
    sizes = new int[np];
    sizes[0] = nedges;
  }
  if (np > 1)
    Teuchos::gather<int, int>(&nedges, 1, sizes, 1, 0, *comm_);

  // prefix sum to build the idx array
  if (me == 0) {
    idx[0] = 0;
    for (int i = 0; i < np; i++)
      idx[i+1] = idx[i] + sizes[i];
  }

  // prepare to send edges
  int cnt = 0;
  part_t *bufv = NULL;
  long *bufw = NULL;
  if (nedges) {
    bufv = new part_t[nedges];
    bufw = new long[nedges];
    // Create buffer with edges (me, part[i]) and weight edges[parts_[i]].
    for (typename std::map<part_t, long>::iterator it = edges.begin();
         it != edges.end(); it++) {
      bufv[cnt] = it->first;  // target part
      bufw[cnt] = it->second; // weight
      cnt++;
    }
  }

  // Prepare to receive edges on rank 0
  part_t *adj = NULL;
  long *wgt = NULL;
  if (me == 0) {
//SYM    adj = new part_t[2*idx[np]];  // need 2x space to symmetrize later
//SYM    wgt = new long[2*idx[np]];  // need 2x space to symmetrize later
    adj = new part_t[idx[np]];
    wgt = new long[idx[np]];
  }

  Teuchos::gatherv<int, part_t>(bufv, cnt, adj, sizes, idx, 0, *comm_);
  Teuchos::gatherv<int, long>(bufw, cnt, wgt, sizes, idx, 0, *comm_);
  delete [] bufv;
  delete [] bufw;

  // Now have constructed graph on rank 0.
  // Call the matching algorithm

  int doRemap;
  if (me == 0) {
    // We have the "LHS" vertices of the bipartite graph; need to create
    // "RHS" vertices.
    for (int i = 0; i < idx[np]; i++) {
      adj[i] += np;  // New RHS vertex number; offset by num LHS vertices
    }

    // Build idx for RHS vertices
    for (part_t i = np; i < tnVtx; i++) {
      idx[i+1] = idx[i];  // No edges for RHS vertices
    }

#ifdef KDDKDD_DEBUG
    std::cout << "IDX ";
    for (part_t i = 0; i <= tnVtx; i++) std::cout << idx[i] << " ";
    std::cout << std::endl;

    std::cout << "ADJ ";
    for (part_t i = 0; i < idx[tnVtx]; i++) std::cout << adj[i] << " ";
    std::cout << std::endl;

    std::cout << "WGT ";
    for (part_t i = 0; i < idx[tnVtx]; i++) std::cout << wgt[i] << " ";
    std::cout << std::endl;
#endif

    // Perform matching on the graph
    part_t *match = new part_t[tnVtx];
    for (part_t i = 0; i < tnVtx; i++) match[i] = i;
    part_t nmatches =
             Zoltan2::GreedyMWM<part_t, long>(idx, adj, wgt, tnVtx, match);

#ifdef KDDKDD_DEBUG
    std::cout << "After matching:  " << nmatches << " found" << std::endl;
    for (part_t i = 0; i < tnVtx; i++)
      std::cout << "match[" << i << "] = " << match[i]
           << ((match[i] != i &&
               (i < np && match[i] != i+np))
                  ? " *" : " ")
           << std::endl;
#endif

    // See whether there were nontrivial changes in the matching.
    bool nontrivial = false;
    if (nmatches) {
      for (part_t i = 0; i < np; i++) {
        if ((match[i] != i) && (match[i] != (i+np))) {
          nontrivial = true;
          break;
        }
      }
    }

    // Process the matches
    if (nontrivial) {
      remap = new part_t[nGlobalParts_];
      for (part_t i = 0; i < nGlobalParts_; i++) remap[i] = -1;

      bool *used = new bool[np];
      for (part_t i = 0; i < np; i++) used[i] = false;

      // First, process all matched parts
      for (part_t i = 0; i < nGlobalParts_; i++) {
        part_t tmp = i + np;
        if (match[tmp] != tmp) {
          remap[i] = match[tmp];
          used[match[tmp]] = true;
        }
      }

      // Second, process unmatched parts; keep same part number if possible
      for (part_t i = 0; i < nGlobalParts_; i++) {
        if (remap[i] > -1) continue;
        if (!used[i]) {
          remap[i] = i;
          used[i] = true;
        }
      }

      // Third, process unmatched parts; give them the next unused part
      for (part_t i = 0, uidx = 0; i < nGlobalParts_; i++) {
        if (remap[i] > -1) continue;
        while (used[uidx]) uidx++;
        remap[i] = uidx;
        used[uidx] = true;
      }
      delete [] used;
    }
    delete [] match;

#ifdef KDDKDD_DEBUG
    std::cout << "Remap vector: ";
    for (part_t i = 0; i < nGlobalParts_; i++) std::cout << remap[i] << " ";
    std::cout << std::endl;
#endif

    long newgstaying = measure_stays(remap, idx, adj, wgt,
                                                      nGlobalParts_, np);
    doRemap = (newgstaying > gstaying);
    std::ostringstream msg;
    msg << "gstaying " << gstaying << " measure(input) "
         << measure_stays(NULL, idx, adj, wgt, nGlobalParts_, np)
         << " newgstaying " << newgstaying
         << " nontrivial " << nontrivial
         << " doRemap " << doRemap << std::endl;
    env_->debug(DETAILED_STATUS, msg.str());
  }
  delete [] idx;
  delete [] sizes;
  delete [] adj;
  delete [] wgt;

  Teuchos::broadcast<int, int>(*comm_, 0, 1, &doRemap);

  if (doRemap) {
    if (me != 0) remap = new part_t[nGlobalParts_];
    Teuchos::broadcast<int, part_t>(*comm_, 0, nGlobalParts_, remap);
    for (size_t i = 0; i < len; i++) {
      parts_[i] = remap[parts_[i]];
    }
  }

  delete [] remap;  // TODO May want to keep for repartitioning as in Zoltan
}


}  // namespace Zoltan2

#endif
