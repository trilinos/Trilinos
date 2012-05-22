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

#include <Zoltan2_Standards.hpp>
#include <Teuchos_DefaultComm.hpp>

namespace Zoltan2{

long getProcessKilobytes();

template <typename scalar_t>
  inline bool outsideRegion(scalar_t val, scalar_t mark, double epsilon){
    return ((val < mark-epsilon) || (val > mark+epsilon));
}

#ifdef HAVE_ZOLTAN2_MPI

RCP<Teuchos::MpiComm<int> > MPI2Teuchos(const MPI_Comm &comm);

RCP<const Teuchos::MpiComm<int> > MPI2TeuchosConst(const MPI_Comm &comm);

MPI_Comm Teuchos2MPI(const RCP<Comm<int> > &comm);

MPI_Comm TeuchosConst2MPI(const RCP<const Comm<int> > &comm);

#endif

} // namespace Zoltan2

#endif
