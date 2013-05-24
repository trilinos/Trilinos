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

/*! \file Zoltan2_PartMapping.hpp
    \brief Defines the PartMapping class.
*/

#ifndef _ZOLTAN2_PARTITIONMAPPING_HPP_
#define _ZOLTAN2_PARTITIONMAPPING_HPP_

namespace Zoltan2 {

long measure_stays(partId_t *, int *, partId_t *, long *, partId_t, partId_t);


/*! \brief PartitionMapping maps a solution or an input distribution to ranks.
*/

template <typename Adapter>
  class PartitionMapping
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::user_t user_t;
#endif

/*! \brief Constructor 
 */

  PartitionMapping( );

// TODO  I don't know whether the following methods are needed, but they
// should not be in the PartitioningSolution, so I am leaving them here for now.

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
  const partId_t *getProcDistribution() const {
    if (procDist_.size() > 0) return &procDist_[0];
    else return NULL;
  }

  /*! \brief Returns the process list corresponding to the global ID list.
      \return The return value is a NULL pointer if part IDs are
                synonomous with process IDs.
   */
  const int *getProcList() const { return procs_.getRawPtr();}

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

  void getPartsForProc(int procId, double &numParts, partId_t &partMin,
    partId_t &partMax) const
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
  void getProcsForPart(partId_t partId, int &procMin, int &procMax) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid part id",
      partId >= 0 && partId < nGlobalParts_, BASIC_ASSERTION);

    partToProcsMap(partId, procMin, procMax);
  }

private:
  void partToProc(bool doCheck, bool haveNumLocalParts, bool haveNumGlobalParts,
    int numLocalParts, int numGlobalParts);

  void procToPartsMap(int procId, double &numParts, partId_t &partMin,
    partId_t &partMax) const;

  void partToProcsMap(partId_t partId, int &procMin, int &procMax) const;

  void setPartDistribution();

  scalar_t localFraction_; // approx fraction of a part on this process

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
  std::vector<partId_t> procDist_;      // or this is defined.
  bool procDistEquallySpread_;        // if procDist_ is used and
                                      // #parts > #procs and
                                      // num_local_parts is not specified,
                                      // parts are evenly distributed to procs

  ////////////////////////////////////////////////////////////////
  // The solution calculates this from the part assignments,
  // unless onePartPerProc_.

  ArrayRCP<int> procs_;       // process rank assigned to gids_[i]
};


template <typename Adapter>
  void PartitioningSolution<Adapter>::setPartDistribution()
{
  // Did the caller define num_global_parts and/or num_local_parts?

  const ParameterList &pl = env_->getParameters();
  size_t haveGlobalNumParts=0, haveLocalNumParts=0;
  int numLocal=0, numGlobal=0;
  double val;

  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("num_global_parts");

  if (pe){
    val = pe->getValue<double>(&val);  // TODO: KDD Skip this double get
    haveGlobalNumParts = 1;            // TODO: KDD Should be unnecessary once
    numGlobal = static_cast<int>(val); // TODO: KDD paramlist handles long long.
    nGlobalParts_ = partId_t(numGlobal); // TODO: KDD  also do below.
  }

  pe = pl.getEntryPtr("num_local_parts");

  if (pe){
    val = pe->getValue<double>(&val);
    haveLocalNumParts = 1;
    numLocal = static_cast<int>(val);
    nLocalParts_ = partId_t(numLocal);
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
    for (partId_t i=1; i <= nGlobalParts_; i++){
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
    throw logic_error("partToProc error");
  }

}



template <typename Adapter>
  void PartitioningSolution<Adapter>::procToPartsMap(int procId,
    double &numParts, partId_t &partMin, partId_t &partMax) const
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
  void PartitioningSolution<Adapter>::partToProcsMap(partId_t partId,
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

      std::vector<partId_t>::const_iterator entry;
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
    catch (std::exception &e){
      throw(std::bad_alloc());
    }

    int *procArray = &procDist_[0];

    try{
      partId_t tmp = partId_t(numLocalParts);
      gatherAll<int, partId_t>(*comm_, 1, &tmp, nprocs, procArray + 1);
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
      catch (std::exception &e){
        throw(std::bad_alloc());
      }

      int *partArray = &partDist_[0];

      double each = floor(fProcs / fParts);
      double extra = fmod(fProcs, fParts);
      partDist_[0] = 0;

      for (partId_t part=0; part < numGlobalParts; part++){
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
      catch (std::exception &e){
        throw(std::bad_alloc());
      }

      int *procArray = &procDist_[0];

      double each = floor(fParts / fProcs);
      double extra = fmod(fParts, fProcs);
      procArray[0] = 0;

      for (int proc=0; proc < nprocs; proc++){
        partId_t numParts = partId_t(each + ((proc<extra) ? 1 : 0));
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



}  // namespace Zoltan2

#endif
