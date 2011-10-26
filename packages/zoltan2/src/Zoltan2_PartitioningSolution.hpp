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

template <typename Adapter>
  class PartitioningSolution : public Solution<Adapter>
{
public:
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lid_t lid_t;

  void setPartitioningSolution(
    size_t length,   // Length of arrays
    gid_t *gids,     // GIDs
    lid_t *lids,     // LIDs
    size_t *parts    // Part assignments
  )
  {
    HELLO;
    gids_ = ArrayView<gid_t>(gids, length);

    if (lids != NULL)
      lids_ = ArrayView<lid_t>(lids, length);
    else     // lids may be NULL
      lids_ = ArrayView<lid_t>(Teuchos::null);

    parts_ = ArrayView<size_t>(parts, length);
  }

  void getPartitioningSolution(
    size_t *length,   // returned: Length of arrays
    gid_t **gids,     // returned: GIDs
    lid_t **lids,     // returned: LIDs
    size_t **parts    // returned: Part assignments
  )
  {
    *length = gids_.size();
    *gids   = gids_.getRawPtr();

    if (lids_ != Teuchos::null) *lids = lids_.getRawPtr();
    else                        *lids = (lid_t *) NULL;

    *parts  = parts_.getRawPtr();
  }

protected:
  // Partitioning solution consists of GIDs, LIDs, and part assignments.
  ArrayView<gid_t>  gids_;
  ArrayView<lid_t>  lids_;
  ArrayView<size_t> parts_;
};

}

#endif
