// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER
//

/*! \file Zoltan2_PartitioningSolution.hpp

    \brief A solution to a partitioning problem
*/


#ifndef _ZOLTAN2_PARTITIONINGSOLUTION_HPP_
#define _ZOLTAN2_PARTITIONINGSOLUTION_HPP_

#include <Zoltan2_IdentifierMap.hpp>
#include <Zoltan2_Solution.hpp>

namespace Zoltan2 {

/*! \brief A PartitioningSolution is a solution to a partitioning problem.

    It is initialized by a PartitioningProblem,
    written to by an algorithm, and may be read by the user or by
    a data migration routine in an input adapter.
*/

// TODO can we template on User?

template <typename User_t>
  class PartitioningSolution : public Solution
{
public:

  typedef typename InputTraits<User_t>::gno_t gno_t;
  typedef typename InputTraits<User_t>::lno_t lno_t;
  typedef typename InputTraits<User_t>::gid_t gid_t;

/*! Constructor
 *
 *   The Solution constructor may require global communication.
 *   The rest of the Solution methods do not.
 *
 *    \param env the environment for the application
 *    \param idMap  the IdentifierMap corresponding to the solution
 *    \param objWeightDim  the number of weights supplied by the application
 *                         for each object.
 *
 *   It is assumed that the parameters in the Environment have
 *   already been committed by the Problem.  
 *
 *   It is assumed that part sizes are uniform.
 */
 
  PartitioningSolution( RCP<const Environment> &env,
    RCP<const IdentifierMap<User_t> > &idMap,
    int userWeightDim);

/*! Constructor
 *
 *   The Solution constructor may require global communication.
 *   The rest of the Solution methods do not.
 *
 *    \param env the environment for the application
 *    \param idMap  the IdentifierMap corresponding to the solution
 *    \param objWeightDim  the number of weights supplied by the application
 *                         for each object.
 *    \param reqPartIds  if the user specified part sizes using
 *          a Problem set method, the reqPartIds[i] is a list of
 *          of part numbers for weight dimension i.
 *    \param reqPartSizes  if the user specified part sizes using
 *          a Problem set method, then reqPartSizes[i] is the list
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
    RCP<const IdentifierMap<User_t> > &idMap,
    int userWeightDim, ArrayView<Array<size_t> > reqPartIds,
    ArrayView<Array<float> > reqPartSizes );
  
  ////////////////////////////////////////////////////////////////////
  // Information that is setup by the Problem for the algorithm.
  // It describes the parts to be created.  The algorithm may
  // query for this information.
  
  size_t getGlobalNumberOfParts() const { return nGlobalParts_; }
  
  size_t getLocalNumberOfParts() const { return nLocalParts_; }
  
  const int *getPartDistribution() const { return partDist_.getRawPtr(); }
  
  const size_t *getProcsParts() const { return procDist_.getRawPtr(); }
  
  bool criteriaHasUniformPartSizes(int idx) const { return !partSizes_[idx].size();}
  
  float *getCriteriaPartSizes(int idx) const { return partSizes_[idx].getRawPtr();}
  
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
   *  TODO: add edge cuts and other metrics of interest as well.
   *      
   */
  
  void setParts(ArrayView<const gno_t> gnoList, ArrayRCP<size_t> &partList,
    ArrayRCP<float> &imbalance);
  
  ////////////////////////////////////////////////////////////////////
  // Results that may be queried by the user or by migration methods.
  // We return raw pointers so users don't have to learn about our
  // pointer wrappers.

  size_t getNumberOfIds() const { return gids_.size(); }

  const gid_t *getGlobalIdList() const { return gids_.getRawPtr(); }

  const size_t *getPartList() const { return parts_.getRawPtr();}

  const float *getImbalance() const { return imbalance_.getRawPtr(); }

private:

  RCP<const Environment> env_;
  RCP<const IdentifierMap<User_t> > idMap_;

  size_t nGlobalParts_;
  size_t nLocalParts_;
  int weightDim_;        // if user has no weights, this is 1

  // proc p's first part is partDist_[p], if p has no parts, partDist[p] is -1,
  // partDist[nprocs] is nparts

  ArrayRCP<int>    partDist_;

  // part i belongs to process procDist_[i], procDist_[nparts] is nprocs

  ArrayRCP<size_t> procDist_;

  // If part sizes were provided for weight i, then partSizes_[i] is a
  //   list of nGlobalParts_ relative sizes.  Sizes sum to 1.0.
  //   If part sizes for weight i are uniform, partSizes_[i] has length zero.
  //
  // NOTE: Although in Zoltan we store part size information for all parts on
  // each process, we should avoid that in Zoltan2.  Perhaps part size
  // information needs to be in a Tpetra_Vector, with remote information
  // obtained if necessary at limited synchronization points in the algorithm.

  ArrayRCP<Array<float> > partSizes_;  // nGlobalParts_ * weightDim_ sizes

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
    RCP<const IdentifierMap<User_t> > &idMap, int userWeightDim)
    : env_(env), idMap_(idMap),
      nGlobalParts_(), nLocalParts_(), weightDim_(),
      partDist_(), procDist_(), partSizes_(), gids_(), parts_(), imbalance_()
{
  using std::string;
  nGlobalParts_ = env_->numProcs_;
  nLocalParts_ = 1;
  weightDim_ = (userWeightDim ? userWeightDim : 1); 

  // TODO:
  // Figure out part/proc assignments from parameters.
  // Set part size information to uniform.
}

template <typename User_t>
  PartitioningSolution<User_t>::PartitioningSolution(
    RCP<const Environment> &env, 
    RCP<const IdentifierMap<User_t> > &idMap, int userWeightDim,
    ArrayView<Array<size_t> > reqPartIds, 
    ArrayView<Array<float> > reqPartSizes)
    : env_(env), idMap_(idMap),
      nGlobalParts_(), nLocalParts_(), weightDim_(),
      partDist_(), procDist_(), partSizes_(), gids_(), parts_(), imbalance_()
{
  using std::string;
  nGlobalParts_ = env_->numProcs_;
  nLocalParts_ = 1;
  weightDim_ = (userWeightDim ? userWeightDim : 1); 

  // TODO:
  // Compute number of parts, part/proc assignments and part size info.
  // Take this from Zoltan function Zoltan_LB.
  // This requires global information.
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
      Z2_LOCAL_MEMORY_ASSERTION(*env_, ngnos, gidList);
  
      gids_ = arcp<const gid_t>(gidList, 0, ngnos);
   
      ArrayView<gid_t> gidView(gidList, ngnos);
      
      // TODO test for exception
      idMap_->gidTranslate(gidView , gnoView, TRANSLATE_LIB_TO_APP);
    }
  
    // Create list of corresponding parts: parts_
  
    size_t *parts = new size_t [ngnos];
    Z2_LOCAL_MEMORY_ASSERTION(*env_, ngnos, parts);
  
    memcpy(parts, partList.getRawPtr(), sizeof(size_t) * ngnos);
  
    parts_ = arcp<size_t>(parts, 0, ngnos);
  }

  // Create imbalance list: one for each weight dimension: weights_
  
  float *imbList = new float [weightDim_];
  Z2_LOCAL_MEMORY_ASSERTION(*env_, weightDim_, imbList);
  memcpy(imbList, imbalance.getRawPtr(), sizeof(float) * weightDim_);
  
  imbalance_ = arcp<float>(imbList, 0, weightDim_);
}

}  // namespace Zoltan2

#endif
