// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Model.hpp

    \brief The abstract interface for a computational model.
*/


#ifndef _ZOLTAN2_PARTITIONINGSOLUTION_HPP_
#define _ZOLTAN2_PARTITIONINGSOLUTION_HPP_

#include <Zoltan2_IdentifierMap.hpp>
#include <Zoltan2_Model.hpp>
#include <Teuchos_Ptr.hpp>

namespace Zoltan2 {

/*! Zoltan2::PartitioningSolution
    Just a placeholder for now.
*/

// TODO can we template on User?
template <typename gid_t, typename lno_t, typename gno_t>
  class PartitioningSolution : public Solution<gid_t, lno_t, gno_t>
{
public:

/*! Constructor
 *
 *   The Solution constructor can require global communication.
 *   The rest of the Solution methods must not.
 *
 *    \param env the environment for the application
 *    \param model  the computational model for the problem
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
 *   dimension "i" are all equal.
 *
 *   If across the application there are some part numbers that are not
 *   included in the reqPartIds lists, then those part sizes are assumed
 *   to be the average of the supplied part sizes.
 */
 
PartitioningSolution(
  RCP<const Environment> &env,
  RPC<const Model> &model,
  ArrayView<Array<size_t> > &reqPartIds,
  ArrayView<Array<double> > &reqPartSizes 
)
  : env_(env), idMap_(model->getIdentifierMap()), 
    nGlobalParts_(), nLocalParts_(), weightDim_(),
    partDist_(), procDist_(), partSizes_(), gids_(), parts_(), imbalance_()
{
  using std::string;
  numGlobalParts_ = env_->numProcs_;
  numLocalParts_ = 1;
  int userWeightDim = model->getNumWeights();
  weightDim_ = (userWeightDim ? userWeightDim : 1); 

  // TODO:
  // Compute number of parts, part/proc assignments and part size info.
  // Take this from Zoltan function Zoltan_LB.
  // This requires global information.
}

  ////////////////////////////////////////////////////////////////////
  // Information that is setup by the Problem for the algorithm.
  // It describes the parts to be created.  The algorithm may
  // query for this information.
  
  size_t getGlobalNumberOfParts() { return nGlobalParts_; }
  
  size_t getLocalNumberOfParts() { return nLocalParts_; }
  
  const int *getPartDistribution() { return partDist_.getRawPtr(); }
  
  const size_t *getProcsParts() { return procDist_.getRawPtr(); }
  
  bool criteriaHasUniformPartSizes(int idx) { return !partSizes_[idx].size();}
  
  double *getCriteriaPartSizes(int idx) { return partSizes_[idx].getRawPtr();}

  ////////////////////////////////////////////////////////////////////
  // Results that are set by the algorithm.

void setParts(ArrayView<const gno_t> &gnoList, ArrayView<size_t> &partList,
  ArrayView<double> &imbalance)
{
  size_t ngnos = gnoList.size();
  
  if (ngnos){
  
    // Create list of application's global IDs: gids_
  
    if (idMap_->gnosAreGids()){
      gids_ = arcp_reinterpret_cast<const gid_t>(gnoList);
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
  
  double *imbList = new double [weightDim_];
  Z2_LOCAL_MEMORY_ASSERTION(*env_, weightDim_, imbList);
  memcpy(imbList, imbalance.getRawPtr(), sizeof(double) * weightDim_);
  
  imbalance_ = arcp<double>(imbList, 0, weightDim_);
}
  ////////////////////////////////////////////////////////////////////
  // Results that may be queried by the user or by migration methods.
  // We return raw pointers so users don't have to learn about our
  // pointer wrappers.

  size_t getNumberOfIds() { return gids_.size(); }

  const gid_t *getGlobalIdList() { return gids_.getRawPtr(); }

  const size_t *getPartList() { return parts_.getRawPtr();}

  double *getImbalance() { return imbalance_.getRawPtr(); }

private:

  ////////////////////////////////////////////////////////////////
  // The Problem gives us the Environment and the Model in the
  // constructor.  We get the IdentifierMap from the Model.
  
  RCP<const Environment> env_;
  RCP<const IdentifierMap> idMap_;  
  
  ////////////////////////////////////////////////////////////////
  // We set the following information in the constructor. Some
  // of it we get from the parameters (which are in the Environment).  
  // We get the number of weights from the Model.  We get the part size 
  // information, if any from arguments to the constructor. Part size 
  // information can be set by the user with Problem set methods.
  
  size_t nGlobalParts_;
  size_t nLocalParts_;
  int weightDim_;
  
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
  
  ArrayRCP<Array<double> > partSizes_;  // nGlobalParts_ * weightDim_ sizes
  
  ////////////////////////////////////////////////////////////////
  // The algorithm sets these values upon completion.
  
  ArrayRCP<const gid_t>  gids_; // User's global IDs from adapter "get" method
  ArrayRCP<size_t> parts_;      // part number assigned to gids_[i]
  ArrayRCP<double> imbalance_;  // weightDim_ imbalance measures
};

}

#endif
