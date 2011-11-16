#ifndef _ZOLTAN2_ALGRCM_HPP_
#define _ZOLTAN2_ALGRCM_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgRCM.hpp
//! \brief RCM ordering of a graph (serial)


// Placeholder for real error handling.
#define KDD_HANDLE_ERROR {\
    cout << __func__ << ":" << __LINE__ << " KDDERROR" << endl;\
    }

namespace Zoltan2{

template <typename Adapter>
int AlgRCM(
  const RCP<GraphModel<Adapter> > &model, 
  const RCP<OrderingSolution<Adapter> > &solution,
  const RCP<Teuchos::ParameterList> &pl,
  const RCP<const Teuchos::Comm<int> > &comm
) 
{
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lid_t lid_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr= 0;

  HELLO; // Test

  return ierr;
}

}
#endif
