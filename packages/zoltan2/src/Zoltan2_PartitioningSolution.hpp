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

template <typename gid_t, typename lno_t>
  class PartitioningSolution : public Solution<gid_t, lno_t>
{
public:

  //////////////////////////////////////////////
  // Constructor allocates memory for the solution.
  PartitioningSolution(
    size_t nparts,
    size_t ngids
  )
  {
    HELLO;
    nParts_ = nparts;
    gids_   = ArrayRCP<gid_t>(ngids);
    parts_  = ArrayRCP<size_t>(ngids);
  }

  //////////////////////////////////////////////
  // Accessor functions, allowing algorithms to get ptrs to solution memory.
  // Algorithms can then load the memory.
  // Non-RCP versions are provided for applications to use.
  inline size_t getNumParts() {return nParts_;}

  inline ArrayRCP<gid_t>  &getGidsRCP()  {return gids_;}
  inline ArrayRCP<size_t> &getPartsRCP() {return parts_;}

  inline ArrayRCP<gid_t>  &getGidsRCPConst()  const
  {
    return const_cast<ArrayRCP<gid_t>& > (gids_);
  }
  inline ArrayRCP<size_t> &getPartsRCPConst() const
  {
    return const_cast<ArrayRCP<size_t>& > (parts_);
  }

  inline gid_t  *getGids(size_t *length)
  {
    *length = gids_.size();
    return gids_.getRawPtr();
  }
  inline size_t *getParts(size_t *length)
  {
    *length = parts_.size();
    return parts_.getRawPtr();
  }

protected:
  // Partitioning solution consists of User's global IDs, listed in
  // the order in which they appeared in the input adapters "get" method,
  // and the part assigned to each global ID.

  size_t nParts_;
  ArrayRCP<gid_t>  gids_;    // User's global IDs from adapter "get" method
  ArrayRCP<size_t> parts_;   // part number assigned to gids_[i]
};

}

#endif
