// @HEADER
// ***********************************************************************
//                Copyright message goes here.   
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_PartToProc.hpp
 *  \brief Map part numbers to process numbers for partitioning.
 */

#ifndef ZOLTAN2_PARTTOPROC_HPP
#define ZOLTAN2_PARTTOPROC_HPP

#include <Zoltan2_Standards.hpp>

namespace Zoltan2{

/*! \brief  Compute the assignment of parts to processes.
 *    \param env an Enviroment for error handling
 *    \param comm the communicator
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
 *   \param oneToOne on return is true if there is a one-to-one assignment
 *       of parts to processes.  Process \c p gets part \c p for all \c p.
 *   \param partDist On return, if !oneToOne and if there are more 
 *        processes than 
 *       parts , then part ID \c p is divided across processes beginning
 *       with <tt> partDist[p] </tt> up to but not including process
 *      <tt> partDist[p+1] </tt>.  
 *   \param procDist On return, if !oneToOne and there are more
 *      parts than processes, the process \c p owns parts
 *        <tt> procDist[p] </tt> up to but not including 
 *       <tt> procDist[p+1] </tt>.
 */

void partToProc(
  const RCP<const Environment> &env,
  const RCP<const Teuchos::Comm<int> > &comm,
  bool doCheck, bool haveNumLocalParts, bool haveNumGlobalParts,
  partId_t numLocalParts, partId_t numGlobalParts, 
  bool &oneToOne, vector<int> &partDist, vector<partId_t> procDist)
{
  int nprocs = comm->getSize();
  int rank = comm->getRank();
  size_t vals[4], reducevals[4];
  size_t sumHaveGlobal, sumHaveLocal;
  size_t sumGlobal, sumLocal;
  size_t maxGlobal, maxLocal;

  partDist.clear();
  procDist.clear();

  if (doCheck){
    vals[4] = {haveNumGlobalParts, haveNumLocalParts,
      numGlobalParts, numLocalParts};

    try{
      reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 4, vals, reducevals);
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
    oneToOne = true;   // default if user did not specify
    return;
  }

  if (haveNumGlobalParts){
    if (doCheck){
      vals[0] = numGlobalParts;
      vals[1] = numLocalParts;
      try{
        reduceAll<int, size_t>(
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
        oneToOne = true;   // user specified one part per proc
        return;
      }
    }
    else{
      if (maxGlobal == nprocs){
        oneToOne = true;   // user specified num parts is num procs
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
      procDist.resize(nprocs+1);
    }
    catch (std::exception &e){
      throw(std::bad_alloc());
    }

    int *procArray = &procDist[0];

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
        partDist.resize(size_t(fParts+1));
      }
      catch (std::exception &e){
        throw(std::bad_alloc());
      }

      int *partArray = &partDist[0];

      double each = floor(fProcs / fParts);
      double extra = fmod(fProcs, fParts);
      partDist[0] = 0;     

      for (partId_t part=0; part < numGlobalParts; part++){
        int numOwners = int(each + ((part<extra) ? 1 : 0));
        partArray[part+1] = partArray[part] + numOwners;
      }

      env_->globalBugAssertion(__FILE__, __LINE__, "#parts != #procs", 
        partDist[numGlobalParts] == nprocs, COMPLEX_ASSERTION, comm_);
    }
    else if (fParts > fProcs){ 

      try{
        procDist.resize(size_t(fProcs+1));
      }
      catch (std::exception &e){
        throw(std::bad_alloc());
      }

      int *procArray = &procDist[0];

      double each = floor(fParts / fProcs);
      double extra = fmod(fParts, fProcs);
      procArray[0] = 0;     

      for (int proc=0; proc < nprocs; proc++){
        partId_t numParts = partId_t(each + ((proc<extra) ? 1 : 0));
        procArray[proc+1] = procArray[proc] + numParts;
      }

      env_->globalBugAssertion(__FILE__, __LINE__, "#parts != #procs", 
        procDist[nprocs] == numGlobalParts, COMPLEX_ASSERTION, comm_);
    }
    else{
      env_->globalBugAssertion(__FILE__, __LINE__, 
        "should never get here", 1, COMPLEX_ASSERTION, comm_);
    }
  }
}

}  // namespace Zoltan2

#endif


