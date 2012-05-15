// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Util.hpp
 *  \brief A gathering of useful namespace methods.
 *  \todo Should each class of utility functions be in a separate source file
 *         instead of having a source file with the unhelpful name of Util?
 */

#ifndef ZOLTAN2_UTIL_HPP
#define ZOLTAN2_UTIL_HPP

#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_DefaultComm.hpp>

namespace Zoltan2{

long getProcessKilobytes();

template <typename scalar_t>
  inline bool outsideRegion(scalar_t val, scalar_t mark, double epsilon){
    return ((val < mark-epsilon) || (val > mark+epsilon));
}

#ifdef HAVE_ZOLTAN2_MPI

/*! \brief Convert an MPI communicator to a MpiComm object.
 */

template <typename Ordinal>
  RCP<MpiComm<Ordinal> >
    getTeuchosMpiComm(const MPI_Comm &comm)
{
  RCP<Teuchos::OpaqueWrapper<MPI_Comm> >handle = 
    Teuchos::opaqueWrapper<MPI_Comm>(comm);
  RCP<MpiComm<Ordinal> > tcommPtr(new MpiComm<Ordinal>(handle));

  return tcommPtr;
}
#endif

} // namespace Zoltan2

#endif
