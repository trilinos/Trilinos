// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_RCP_STD_SHAREDPTR_CONVERSIONS_HPP
#define TEUCHOS_RCP_STD_SHAREDPTR_CONVERSIONS_HPP

#include "Teuchos_RCPStdSharedPtrConversionsDecl.hpp"
#include "Teuchos_RCP.hpp"


template<class T>
Teuchos::RCP<T>
Teuchos::rcp( const std::shared_ptr<T> &sptr )
{
  if (sptr.get()) {
    // First, see if the RCP is in the shared_ptr deleter object
    const StdSharedPtrRCPDeleter<T>
      *rcpd = std::get_deleter<StdSharedPtrRCPDeleter<T> >(sptr);
    if (rcpd) {
      return rcpd->ptr();
    }
#ifdef TEUCHOS_DEBUG
    // Second, see if the an RCP node pointing to this type already exists
    // from being wrapped already from a prior call to this function where the
    // add_new_RCPNode(...) function could have been called already!.
    RCPNode* existingRCPNode = RCPNodeTracer::getExistingRCPNode(sptr.get());
    if (existingRCPNode) {
      return RCP<T>(sptr.get(), RCPNodeHandle(existingRCPNode, RCP_STRONG, false));
    }
#endif
    // Lastly, we just create a new RCP and RCPNode ...
    return rcpWithDealloc(sptr.get(), DeallocStdSharedPtr<T>(sptr), true);
  }
  return null;
}


template<class T>
std::shared_ptr<T>
Teuchos::get_shared_ptr( const RCP<T> &rcp )
{
  if (nonnull(rcp)) {
    Ptr<const DeallocStdSharedPtr<T> >
      dbsp = get_optional_dealloc<DeallocStdSharedPtr<T> >(rcp);
    if (nonnull(dbsp)) {
      return dbsp->ptr();
    }
    return std::shared_ptr<T>(rcp.get(), StdSharedPtrRCPDeleter<T>(rcp));
  }
  return std::shared_ptr<T>();
}


#endif	// TEUCHOS_RCP_STD_SHAREDPTR_CONVERSIONS_HPP
