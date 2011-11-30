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

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Solution.hpp>

namespace Zoltan2 {

/*! Zoltan2::PartitioningSolution
    Just a placeholder for now.
*/

template <typename gid_t, typename lid_t, typename lno_t>
  class PartitioningSolution : public Solution<gid_t, lid_t, lno_t>
{
public:

  //////////////////////////////////////////////
  // Constructor allocates memory for the solution.
  PartitioningSolution(
    size_t nparts,
    size_t ngids,
    size_t nlids
  )
  {
    HELLO;
    nParts_ = nparts;
    gids_   = ArrayRCP<gid_t>(ngids);
    if (nlids) lids_   = ArrayRCP<lid_t>(nlids);
    else       lids_   = ArrayRCP<lid_t>(Teuchos::null);
    parts_  = ArrayRCP<size_t>(ngids);
  }

  //////////////////////////////////////////////
  // Accessor functions, allowing algorithms to get ptrs to solution memory.
  // Algorithms can then load the memory.  
  // Non-RCP versions are provided for applications to use.
  inline size_t getNumParts() {return nParts_;}

  inline ArrayRCP<gid_t>  getGidsRCP()  {return gids_;}
  inline ArrayRCP<lid_t>  getLidsRCP()  {return lids_;}
  inline ArrayRCP<size_t> getPartsRCP() {return parts_;}

  inline gid_t  *getGids(size_t *length) {
    *length = gids_.size(); return gids_.getRawPtr();
  }
  inline lid_t  *getLids(size_t *length) {
    *length = lids_.size();
    return (lids_.is_null() ? (lid_t*) NULL : lids_.getRawPtr());
  }
  inline size_t *getParts(size_t *length) {
    *length = parts_.size(); return parts_.getRawPtr();
  }

protected:
  // Partitioning solution consists of GIDs, LIDs, and part assignments.
  size_t nParts_;
  ArrayRCP<gid_t>  gids_;
  ArrayRCP<lid_t>  lids_;
  ArrayRCP<size_t> parts_;
};

}

#endif
