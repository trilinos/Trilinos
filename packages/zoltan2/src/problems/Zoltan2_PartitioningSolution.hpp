// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER
//

/*! \file Zoltan2_PartitioningSolution.hpp
    \brief Defines the PartitioningSolution class.
*/

#ifndef _ZOLTAN2_PARTITIONINGSOLUTION_HPP_
#define _ZOLTAN2_PARTITIONINGSOLUTION_HPP_

#include <Zoltan2_IdentifierMap.hpp>
#include <Zoltan2_Solution.hpp>

#include <cmath>

namespace Zoltan2 {

/*! \brief A PartitioningSolution is a solution to a partitioning problem.

    It is initialized by a PartitioningProblem,
    written to by an algorithm, and may be read by the user or by
    a data migration routine in an input adapter.
    
    \todo handle the case where number of parts does not equal 
                     number of processes
    \todo handle more metrics
*/

template <typename User_t>
  class PartitioningSolution : public Solution
{
public:

  typedef typename InputTraits<User_t>::gno_t gno_t;
  typedef typename InputTraits<User_t>::lno_t lno_t;
  typedef typename InputTraits<User_t>::gid_t gid_t;

/*! \brief Constructor when part sizes are not supplied.
 *
 *   The Solution constructor may require global communication.
 *   The rest of the Solution methods do not.
 *
 *    \param env the environment for the application
 *    \param comm the communicator for the problem associated with this solution
 *    \param idMap  the IdentifierMap corresponding to the solution
 *    \param objWeightDim  the number of weights supplied by the application
 *                         for each object.
 *
 *   It is assumed that the parameters in the Environment have
 *   already been committed by the Problem.  
 *
 *   It is possible that part sizes were supplied on other processes,
 *   so this constructor does do a check to see if part sizes need
 *   to be globally calculated.
 */
 
  PartitioningSolution( RCP<const Environment> &env,
    RCP<const Comm<int> > &comm,
    RCP<const IdentifierMap<User_t> > &idMap,
    int userWeightDim);

/*! \brief Constructor when part sizes are supplied.
 *
 *   The Solution constructor may require global communication.
 *   The rest of the Solution methods do not.
 *
 *    \param env the environment for the application
 *    \param comm the communicator for the problem associated with this solution
 *    \param idMap  the IdentifierMap corresponding to the solution
 *    \param objWeightDim  the number of weights supplied by the application
 *                         for each object.
 *    \param reqPartIds  reqPartIds[i] is a list of
 *          of part numbers for weight dimension i.
 *    \param reqPartSizes  reqPartSizes[i] is the list
 *          of part sizes for weight i corresponding to parts in 
 *          reqPartIds[i]
 *
 *   It is assumed that the parameters in the Environment have
 *   already been committed by the Problem.  
 *
 *   If reqPartIds[i].size() and reqPartSizes[i].size() are zero for 
 *   all processes, it is assumed that part sizes for weight 
 *   dimension "i" are uniform.
 *
 *   If across the application there are some part numbers that are not
 *   included in the reqPartIds lists, then those part sizes are assumed
 *   to be the average of the supplied part sizes.
 */

  PartitioningSolution( RCP<const Environment> &env,
    RCP<const Comm<int> > &comm,
    RCP<const IdentifierMap<User_t> > &idMap,
    int userWeightDim, ArrayView<ArrayRCP<size_t> > reqPartIds,
    ArrayView<ArrayRCP<float> > reqPartSizes );
  
  ////////////////////////////////////////////////////////////////////
  // Information that is setup by the Problem for the algorithm.
  // It describes the parts to be created.  The algorithm may
  // query for this information.
  
/*! \brief Returns the global number of parts requested by the user
    \return the global number of parts
 */
  size_t getGlobalNumberOfParts() const { return nGlobalParts_; }
  
/*! \brief Returns the number of parts to be assigned to this process
    \return the number of parts
 */
  size_t getLocalNumberOfParts() const { return nLocalParts_; }
  
/*! \brief Return whether or not the part to process distribution is one-to-one.
     \return true if Process p owns part p for all p, and false if the part
                  to process distribution is more complex.
 */

  bool oneToOnePartDistribution() { return onePartPerProc_; }
  
/*! \brief If we do not have simply one part per process, then get the 
        distribution of parts across processes.

    \return an array A such that A[i] is process owning part i.  The
       length of the array is one greater than the global number of parts.
       The value of the last element is the number of processes that have parts.

     If Process p has part p for all p, then getPartDistribution returns NULL.
 */
  const int *getPartDistribution() const { return partDist_.getRawPtr(); }
  
/*! \brief If we do not have simply one part per process, then get the 
                 distribution of processes to parts.

    \return an array A such that A[i] is the first part 
      assigned to process i.  (Parts are assigned sequentially to processes.)
      However if A[i+1] is equal to A[i], there
      are no parts assigned to process i.  The length of the array is
      one greater than the number of processes.  The last element of
      the array is the global number of parts.

     If Process p has part p for all p, then getProcDistribution returns NULL.
 */
  const size_t *getProcDistribution() const { return procDist_.getRawPtr(); }
  
/*! \brief Determine if balancing criteria (weight dimension) has uniform
                parts.
    \param idx   A value from 0 to one less than the number of weights per 
                   object.
    \return true if part sizes are uniform for this criteria.
 */
  bool criteriaHasUniformPartSizes(int idx) const { return pSizeUniform_[idx];}
  
/*! \brief Get the size for a given weight dimension and a given part. 

    \param idx   A value from 0 to one less than the number of weights per 
                       object.
    \param part  A value from 0 to one less than the global number of parts 
                   to be computed
    \return   The size for that part.  Part sizes for a given weight 
                    dimension sum to 1.0.

      \todo It would be useful to algorithms to get the sum of
           part sizes from a to b, or the sum or a list of parts.
 */
  float getCriteriaPartSize(int idx, size_t part) const { 
    if (pSizeUniform_[idx]) 
      return 1.0 / nGlobalParts_;
    else if (pCompactIndex_[idx].size())
      return pSize_[idx][pCompactIndex_[idx][part]];
    else
      return pSize_[idx][part];
  }

  
  ////////////////////////////////////////////////////////////////////
  // Method used by the algorithm to set results.

  /*! \brief The algorithm uses setParts to set the solution.
   *
   *   \param gnoList  A view of the list of global numbers that was
   *     supplied to the algorithm by the model. If the model used
   *     the application's global IDs as our internal global numbers,
   *     this may be a view of the global IDs in the application's memory.
   *     Otherwise it is a view of the internal global numbers created
   *     by the IdentifierMap, or a re-ordered set of global numbers
   *     created during the course of the algorithm.
   *
   *   \param partList  The part assigned to gnoList[i] by the algorithm
   *      should be in partList[i].  The partList is allocated and written
   *      by the algorithm.
   *
   *   \param  imbalance  The algorithm computes one imbalance for each
   *      weight dimension.  It allocates and writes the imbalance array.
   *
   */
  
  void setParts(ArrayView<const gno_t> gnoList, ArrayRCP<size_t> &partList,
    ArrayRCP<float> &imbalance);
  
  ////////////////////////////////////////////////////////////////////
  // Results that may be queried by the user or by migration methods.
  // We return raw pointers so users don't have to learn about our
  // pointer wrappers.

  /*! \brief Returns the local number of Ids.
   */
  size_t getNumberOfIds() const { return gids_.size(); }

  /*! \brief Returns the user's global ID list.
   */
  const gid_t *getGlobalIdList() const { return gids_.getRawPtr(); }

  /*! \brief Returns the part list corresponding to the global ID list.
   */
  const size_t *getPartList() const { return parts_.getRawPtr();}

  /*! \brief Returns the imbalance of the solution.
   */
  const float *getImbalance() const { return imbalance_.getRawPtr(); }

private:

  void setPartDistribution();

  void setPartSizes(ArrayView<ArrayRCP<size_t> > reqPartIds,
    ArrayView<ArrayRCP<float> > reqPartSizes);

  void computePartSizes(int wdim, ArrayView<size_t> ids, 
    ArrayView<float> sizes);

  void broadcastPartSizes(int wdim);

  RCP<const Environment> env_;             // has application communicator
  RCP<const Comm<int> > comm_;             // the problem communicator
  RCP<const IdentifierMap<User_t> > idMap_;

  size_t nGlobalParts_;
  size_t nLocalParts_;
  int weightDim_;        // if user has no weights, this is 1

  // If process p is to be assigned part p for all p, then onePartPerProc_ 
  // is true.
  // Otherwise it is false, and partDist_ and procDist_ describe the
  // allocation of parts to processes.
  //
  // Process p's first part is procDist_[p].  (But if procDist_[p+1] equals
  // procDist_[p], then process p has no parts.) procDist[nprocs] is nparts.
  //
  // Part i belongs to process partDist_[i], and partDist_[nparts] is the
  // number of processes that have parts.

  bool             onePartPerProc_;
  ArrayRCP<int>    partDist_;
  ArrayRCP<size_t> procDist_;

  // In order to minimize the storage required for part sizes, we
  // have three different representations.
  //
  // If the part sizes for weight dimension w are all the same, then:
  //    pSizeUniform_[w] = true
  //    pCompactIndex_[w].size() = 0
  //    pSize_[w].size() = 0
  //
  // and the size for part p is 1.0 / nparts.
  //
  // If part sizes differ for each part in weight dimension w, but there
  // are no more than 64 distinct sizes:
  //    pSizeUniform_[w] = false
  //    pCompactIndex_[w].size() = number of parts
  //    pSize_[w].size() = number of different sizes
  //
  // and the size for part p is pSize_[pCompactIndex_[p]].
  //
  // If part sizes differ for each part in weight dimension w, and there
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
  ArrayRCP<ArrayRCP<float> > pSize_;

  ////////////////////////////////////////////////////////////////
  // The algorithm sets these values upon completion.

  ArrayRCP<const gid_t>  gids_; // User's global IDs from adapter "get" method
  ArrayRCP<size_t> parts_;      // part number assigned to gids_[i]
  ArrayRCP<float> imbalance_;  // weightDim_ imbalance measures
};

////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////

template <typename User_t>
  PartitioningSolution<User_t>::PartitioningSolution(
    RCP<const Environment> &env, 
    RCP<const Comm<int> > &comm,
    RCP<const IdentifierMap<User_t> > &idMap, int userWeightDim)
    : env_(env), comm_(comm), idMap_(idMap),
      nGlobalParts_(0), nLocalParts_(0), weightDim_(),
      onePartPerProc_(false), partDist_(), procDist_(), 
      pSizeUniform_(), pCompactIndex_(), pSize_(),
      gids_(), parts_(), imbalance_()
{
  weightDim_ = (userWeightDim ? userWeightDim : 1); 

  setPartDistribution();

  // We must call setPartSizes() because part sizes may have
  // been provided by the user on other processes.

  ArrayRCP<size_t> *noIds = new ArrayRCP<size_t> [weightDim_];
  ArrayRCP<float> *noSizes = new ArrayRCP<float> [weightDim_];
  ArrayRCP<ArrayRCP<size_t> > ids(noIds, 0, weightDim_, true);
  ArrayRCP<ArrayRCP<float> > sizes(noSizes, 0, weightDim_, true);

  setPartSizes(ids.view(0, weightDim_), sizes.view(0, weightDim_));
}

template <typename User_t>
  PartitioningSolution<User_t>::PartitioningSolution(
    RCP<const Environment> &env, 
    RCP<const Comm<int> > &comm,
    RCP<const IdentifierMap<User_t> > &idMap, int userWeightDim,
    ArrayView<ArrayRCP<size_t> > reqPartIds, 
    ArrayView<ArrayRCP<float> > reqPartSizes)
    : env_(env), comm_(comm), idMap_(idMap),
      nGlobalParts_(0), nLocalParts_(0), weightDim_(),
      onePartPerProc_(false), partDist_(), procDist_(), 
      pSizeUniform_(), pCompactIndex_(), pSize_(),
      gids_(), parts_(), imbalance_()
{
  weightDim_ = (userWeightDim ? userWeightDim : 1); 

  setPartDistribution();

  setPartSizes(reqPartIds, reqPartSizes);
}

template <typename User_t>
  void PartitioningSolution<User_t>::setPartDistribution()
{
  int rank = env_->myRank_;
  int nprocs = env_->numProcs_;
  
  // Did the caller define num_global_parts and/or num_local_parts?

  size_t haveGlobalNumParts=0, haveLocalNumParts=0;

  const ParameterList *params = NULL;
  if (env_->hasPartitioningParameters()){
    const ParameterList &allParams = env_->getParameters();
    const ParameterList &partitioningParams = allParams.sublist("partitioning");
    params = &partitioningParams;

    const double *entry1 = params->getPtr<double>(string("num_global_parts"));
    const double *entry2 = params->getPtr<double>(string("num_local_parts"));
  
    if (entry1){
      haveGlobalNumParts = 1;
      double val = *entry1;
      nGlobalParts_ = static_cast<size_t>(val);
    }
  
    if (entry2){
      haveLocalNumParts = 1;
      double val = *entry2;
      nLocalParts_ = static_cast<size_t>(val);
    }
  }

  size_t vals[4] = {haveGlobalNumParts, haveLocalNumParts, 
    nGlobalParts_, nLocalParts_};
  size_t reducevals[4];

  try{
    reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 4, vals, reducevals);
  }
  Z2_THROW_OUTSIDE_ERROR(*env_, e);

  size_t sumHaveGlobal = reducevals[0];
  size_t sumHaveLocal = reducevals[1];
  size_t sumGlobal = reducevals[2];
  size_t sumLocal = reducevals[3];

  env_->localInputAssertion(__FILE__, __LINE__,
    "Either all procs specify num_global/local_parts or none do",
    (sumHaveGlobal == 0 || sumHaveGlobal == nprocs) &&
    (sumHaveLocal == 0 || sumHaveLocal == nprocs), 
    BASIC_ASSERTION);

  if (sumHaveGlobal == 0 && sumHaveLocal == 0){
    onePartPerProc_ = true;   // default if user did not specify
    nGlobalParts_ = nprocs;
    nLocalParts_ = 1;
    return;
  }

  if (sumHaveGlobal == nprocs){
    vals[0] = nGlobalParts_;
    vals[1] = nLocalParts_;
    try{
      reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_MAX, 2, vals, reducevals);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_, e);

    size_t maxGlobal = reducevals[0];
    size_t maxLocal = reducevals[1];

    env_->localInputAssertion(__FILE__, __LINE__,
      "Value for num_global_parts is different on different processes.",
      maxGlobal * nprocs == sumGlobal, BASIC_ASSERTION);

    if (sumLocal){
      env_->localInputAssertion(__FILE__, __LINE__,
        "Sum of num_local_parts does not equal requested num_global_parts",
        sumLocal == nGlobalParts_, BASIC_ASSERTION);

      if (sumLocal == nprocs && maxLocal == 1){
        onePartPerProc_ = true;   // user specified one part per proc
        return;
      }
    }
    else{
      if (maxGlobal == nprocs){
        onePartPerProc_ = true;   // user specified num parts is num procs
        nLocalParts_ = 1;
        return;
      }
    }
  }

  size_t *nparts = new size_t [nprocs];
  env_->globalMemoryAssertion(__FILE__, __LINE__, nprocs, nparts, comm_);

  if (sumHaveLocal == nprocs){
    //
    // We will go by the number of local parts specified.
    //

    try{
      gatherAll<int, size_t>(*comm_, 1, &nLocalParts_, nprocs, nparts); 
    }
    Z2_THROW_OUTSIDE_ERROR(*env_, e);

    nGlobalParts_ = sumLocal;

  }
  else{
    //
    // We will divide the global number of parts across the processes.
    //

    float f1 = nGlobalParts_;
    float f2 = nprocs;

    float each = floor(f1 / f2);
    float extra = fmod(f1, f2);

    for (int p=0; p < nprocs; p++){
      if (p < extra)
        nparts[p] = size_t(each + 1);
      else
        nparts[p] = size_t(each);
    }

    nLocalParts_ = nparts[rank];
  }

  // partDist_[part] = process
  // procDist_[process] = part

  int *partArray = new int [nGlobalParts_+1];
  size_t *procArray = new size_t [nprocs+1];

  env_->globalMemoryAssertion(__FILE__, __LINE__, nprocs+nGlobalParts_+2,
    partArray && procArray, comm_);

  partDist_ = arcp(partArray, 0, nGlobalParts_+1, true);
  procDist_ = arcp(procArray, 0, nprocs+1, true);

  procDist_[0] = 0;

  for (int p=0; p < nprocs; p++){
    procDist_[p+1] = procDist_[p] + nparts[p];

    for (size_t part=procDist_[p]; part < procDist_[p+1]; part++)
      partDist_[part] = p;
  }

  delete [] nparts;

  partDist_[nGlobalParts_] = partDist_[nGlobalParts_-1] + 1;
}

template <typename User_t>
  void PartitioningSolution<User_t>::setPartSizes(
    ArrayView<ArrayRCP<size_t> > ids, ArrayView<ArrayRCP<float> > sizes)
{
  int wdim = weightDim_;
  bool fail=false;
  size_t *counts = new size_t [wdim*2];

  fail = ((ids.size() != wdim) || (sizes.size() != wdim));

  for (int w=0; !fail && w < wdim; w++){
    counts[w] = ids[w].size();
    if (ids[w].size() != sizes[w].size()) fail=true;
  }

  env_->globalBugAssertion(__FILE__, __LINE__, "bad argument arrays", fail==0, 
    COMPLEX_ASSERTION, comm_);

  // Are all part sizes the same?  This is the common case.

  ArrayRCP<float> *emptySizes= new ArrayRCP<float> [wdim];
  pSize_ = arcp(emptySizes, 0, wdim);

  ArrayRCP<unsigned char> *emptyIndices= new ArrayRCP<unsigned char> [wdim];
  pCompactIndex_ = arcp(emptyIndices, 0, wdim);

  bool *info = new bool [wdim];
  pSizeUniform_ = arcp(info, 0, wdim);
  for (int w=0; w < wdim; w++)
    pSizeUniform_[w] = true;

  if (nGlobalParts_ == 1){   
    return;   // there's only one part in the whole problem
  }

  try{
    reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_MAX, wdim, counts, 
      counts + wdim);
  }
  Z2_THROW_OUTSIDE_ERROR(*env_, e);

  bool zero = true;

  for (int w=0; w < wdim; w++)
    if (counts[wdim+w] > 0){
      zero = false;
      pSizeUniform_[w] = false;
    }

  if (zero) // Part sizes for all criteria are uniform.
    return; 

  // Compute the part sizes for criteria for which part sizes were
  // supplied.  Normalize for each criteria so part sizes sum to one.

  int nprocs = comm_->getSize();
  int rank = comm_->getRank();
  Array<long> sendCount(nprocs, 0);
  ArrayRCP<long> recvCount;
  ArrayRCP<float> recvSizes;
  ArrayRCP<size_t> recvIds;

  for (int w=0; w < wdim; w++){
    if (pSizeUniform_[w]) continue;
    
    // Send all ids and sizes to rank 0

    long len = ids[w].size();
    sendCount[0] = len;

    try{
      AlltoAllv<size_t, long>(*comm_, *env_, ids[w].view(0, len),
        sendCount.view(0, nprocs), recvIds, recvCount);
    }
    Z2_FORWARD_EXCEPTIONS

    try{
      AlltoAllv<float, long>(*comm_, *env_, sizes[w].view(0, len),
        sendCount.view(0, nprocs), recvSizes, recvCount);
    }
    Z2_FORWARD_EXCEPTIONS

    if (rank == 0){
      try{
        size_t numVals = recvIds.size();
        computePartSizes(w, recvIds.view(0,numVals), recvSizes.view(0,numVals));
      }
      Z2_FORWARD_EXCEPTIONS
    }

    broadcastPartSizes(w);
  } 
}

template <typename User_t>
  void PartitioningSolution<User_t>::broadcastPartSizes(int wdim)
{
  env_->localBugAssertion(__FILE__, __LINE__, "preallocations", 
    pSize_.size()>wdim && 
    pSizeUniform_.size()>wdim && pCompactIndex_.size()>wdim, 
    COMPLEX_ASSERTION);

  int rank = env_->myRank_;
  int nprocs = env_->numProcs_;
  size_t nparts = nGlobalParts_;

  if (nprocs < 2)
    return;

  char flag=0;

  if (rank == 0){
    if (pSizeUniform_[wdim] == true)
      flag = 1;
    else if (pCompactIndex_[wdim].size() > 0)
      flag = 2;
    else
      flag = 3;
  }

  try{
    Teuchos::broadcast<int, char>(*comm_, 0, 1, &flag);
  }
  Z2_THROW_OUTSIDE_ERROR(*env_, e);

  if (flag == 1){
    if (rank > 0)
      pSizeUniform_[wdim] = true;

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
        pCompactIndex_[wdim].size() == nparts, COMPLEX_ASSERTION);
      idxbuf = pCompactIndex_[wdim].getRawPtr();
    }

    try{
      // broadcast of unsigned char is not supported
      Teuchos::broadcast<int, char>(*comm_, 0, nparts, 
        reinterpret_cast<char *>(idxbuf));
    }
    Z2_THROW_OUTSIDE_ERROR(*env_, e);

    if (rank > 0)
      pCompactIndex_[wdim] = arcp(idxbuf, 0, nparts);

    // broadcast the list of different part sizes

    unsigned char maxIdx=0;
    for (size_t p=0; p < nparts; p++)
      if (idxbuf[p] > maxIdx) maxIdx = idxbuf[p];

    int numSizes = maxIdx + 1;
  
    float *sizeList = NULL;

    if (rank > 0){
      sizeList = new float [numSizes];
      env_->localMemoryAssertion(__FILE__, __LINE__, numSizes, sizeList);
    }
    else{
      env_->localBugAssertion(__FILE__, __LINE__, "wrong number of sizes", 
        numSizes == pSize_[wdim].size(), COMPLEX_ASSERTION);

      sizeList = pSize_[wdim].getRawPtr();
    }

    try{
      Teuchos::broadcast<int, float>(*comm_, 0, numSizes, sizeList);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_, e);

    if (rank > 0)
      pSize_[wdim] = arcp(sizeList, 0, numSizes);

    return;
  }

  if (flag == 3){

    // broadcast the size of each part

    float *sizeList = NULL;

    if (rank > 0){
      sizeList = new float [nparts];
      env_->localMemoryAssertion(__FILE__, __LINE__, nparts, sizeList);
    }
    else{
      env_->localBugAssertion(__FILE__, __LINE__, "wrong number of sizes", 
        nparts == pSize_[wdim].size(), COMPLEX_ASSERTION);

      sizeList = pSize_[wdim].getRawPtr();
    }

    try{
      Teuchos::broadcast<int, float>(*comm_, 0, nparts, sizeList);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_, e);

    if (rank > 0)
      pSize_[wdim] = arcp(sizeList, 0, nparts);

    return;
  }
}

template <typename User_t>
  void PartitioningSolution<User_t>::computePartSizes(int wdim,
    ArrayView<size_t> ids, ArrayView<float> sizes)
{
  int len = ids.size();

  if (len == 1){
    pSizeUniform_[wdim] = true;
    return;
  }

  env_->localBugAssertion(__FILE__, __LINE__, "bad array sizes", 
    len>0 && sizes.size()==len, COMPLEX_ASSERTION);

  env_->localBugAssertion(__FILE__, __LINE__, "bad index", 
    wdim>=0 && wdim<weightDim_, COMPLEX_ASSERTION);

  env_->localBugAssertion(__FILE__, __LINE__, "preallocations", 
    pSize_.size()>wdim && 
    pSizeUniform_.size()>wdim && pCompactIndex_.size()>wdim, 
    COMPLEX_ASSERTION);

  // Check ids and sizes and find min, max and average sizes.

  size_t nparts = nGlobalParts_;
  unsigned char *buf = new unsigned char [nparts];
  env_->localMemoryAssertion(__FILE__, __LINE__, nparts, buf);
  memset(buf, 0, nparts);
  ArrayRCP<unsigned char> partIdx(buf, 0, nparts, true);

  float min, max, sum;

  for (int i=0; i < len; i++){
    size_t id = ids[i];
    float size = sizes[i];

    env_->localInputAssertion(__FILE__, __LINE__, "invalid part id", 
      id>=0 && id<nparts, BASIC_ASSERTION);

    env_->localInputAssertion(__FILE__, __LINE__, "invalid part size", size>=0,
      BASIC_ASSERTION);

    // TODO: we could allow users to specify multiple sizes for the same
    // part if we add a parameter that says what we are to do with them:
    // add them or take the max.

    env_->localInputAssertion(__FILE__, __LINE__, 
      "multiple sizes provided for one part", partIdx[id]==0, BASIC_ASSERTION);

    partIdx[id] = 1;

    if (i==0 || size < min) min = size;
    if (i==0 || size > max) max = size;
    sum += size;
  }

  if (sum == 0){   // user has given us a list of parts of size 0
    
    float *allSizes = new float [2];
    env_->localMemoryAssertion(__FILE__, __LINE__, 2, allSizes);

    ArrayRCP<float> sizeArray(allSizes, 0, 2, true);

    allSizes[0] = 0.0;
    allSizes[1] = 1.0 / (nparts - len);

    for (size_t p=0; p < nparts; p++)
      buf[p] = 1;                 // index to default part size

    for (int i=0; i < len; i++)
      buf[ids[i]] = 0;            // index to part size zero
    
    pSize_[wdim] = sizeArray;
    pCompactIndex_[wdim] = partIdx;

    return;
  }

  double epsilon = (max - min) * 10e-5;  // to distinguish different sizes

  if (max - min <= epsilon){
    pSizeUniform_[wdim] = true;
    return;
  }

  float avg = sum / len;    // size for parts that were not specified

  // We are going to merge part sizes that are very close.  This takes
  // computation time now, but can save considerably in the storage of
  // all part sizes on each process.  For example, a common case may
  // be some parts are size 1 and all the rest are size 2.

  float *tmp = new float [len];
  env_->localMemoryAssertion(__FILE__, __LINE__, len, tmp);
  memcpy(tmp, sizes.getRawPtr(), sizeof(float) * len);
  ArrayRCP<float> partSizes(tmp, 0, len, true);

  std::sort(partSizes.begin(), partSizes.end());

  // create a list of sizes that are unique within epsilon

  Array<float> nextUniqueSize;
  nextUniqueSize.push_back(partSizes[len-1]);   // largest
  float curr = partSizes[len-1];
  int avgIndex = len;
  bool haveAvg = false;
  if (curr - avg <= epsilon)
     avgIndex = 0;

  for (int i=len-2; i >= 0; i--){
    float val = partSizes[i];
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
    
    float *allSizes = new float [sizeArrayLen];
    env_->localMemoryAssertion(__FILE__, __LINE__, sizeArrayLen, allSizes);
    ArrayRCP<float> sizeArray(allSizes, 0, sizeArrayLen, true);

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
      float size = sizes[i];
      int index;

      // Find the first size greater than or equal to this size.

      if (size < avg && avg - size <= epsilon)
        index = newAvgIndex;
      else{
        ArrayRCP<float>::iterator found = std::lower_bound(sizeArray.begin(),
          sizeArray.end(), size);

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

    pCompactIndex_[wdim] = partIdx;
    pSize_[wdim] = sizeArray;
  }
  else{
    // To have access to part sizes, we must store nparts floats on 
    // every process.  We expect this is a rare case.

    partIdx.clear();

    tmp = new float [nparts];
    env_->localMemoryAssertion(__FILE__, __LINE__, nparts, tmp);

    sum += ((nparts - len) * avg);

    for (int i=0; i < nparts; i++){
      tmp[i] = avg/sum;
    }

    for (int i=0; i < len; i++){
      tmp[ids[i]] = sizes[i]/sum;
    }

    pSize_[wdim] = arcp(tmp, 0, nparts);
  }
}


template <typename User_t>
  void PartitioningSolution<User_t>::setParts(
    ArrayView<const gno_t> gnoList, ArrayRCP<size_t> &partList,
    ArrayRCP<float> &imbalance)
{
  size_t ngnos = gnoList.size();
  
  if (ngnos){
  
    // Create list of application's global IDs: gids_
  
    if (idMap_->gnosAreGids()){
      gids_ = Teuchos::arcpFromArrayView<const gid_t>(gnoList);
    }
    else{
      ArrayView<gno_t> gnoView = av_const_cast<gno_t>(gnoList);
  
      gid_t *gidList = new gid_t [gnoList.size()];
      env_->localMemoryAssertion(__FILE__, __LINE__, ngnos, gidList);
  
      gids_ = arcp<const gid_t>(gidList, 0, ngnos);
   
      ArrayView<gid_t> gidView(gidList, ngnos);
      
      // TODO test for exception
      idMap_->gidTranslate(gidView , gnoView, TRANSLATE_LIB_TO_APP);
    }
  
    // Create list of corresponding parts: parts_
  
    size_t *parts = new size_t [ngnos];
    env_->localMemoryAssertion(__FILE__, __LINE__, ngnos, parts);
  
    memcpy(parts, partList.getRawPtr(), sizeof(size_t) * ngnos);
  
    parts_ = arcp<size_t>(parts, 0, ngnos);
  }

  // Create imbalance list: one for each weight dimension: weights_
  
  float *imbList = new float [weightDim_];
  env_->localMemoryAssertion(__FILE__, __LINE__, weightDim_, imbList);
  memcpy(imbList, imbalance.getRawPtr(), sizeof(float) * weightDim_);
  
  imbalance_ = arcp<float>(imbList, 0, weightDim_);
}

}  // namespace Zoltan2

#endif
