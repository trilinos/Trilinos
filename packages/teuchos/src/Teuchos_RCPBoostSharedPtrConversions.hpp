// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_RCP_SHAREDPTR_CONVERSIONS_HPP
#define TEUCHOS_RCP_SHAREDPTR_CONVERSIONS_HPP

#include "Teuchos_RCPBoostSharedPtrConversionsDecl.hpp"
#include "Teuchos_RCP.hpp"


template<class T>
Teuchos::RCP<T>
Teuchos::rcp( const boost::shared_ptr<T> &sptr )
{
  if (sptr.get()) {
    // First, see if the RCP is in the shared_ptr deleter object
    const RCPDeleter<T>
      *rcpd = boost::get_deleter<RCPDeleter<T> >(sptr);
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
    return rcp(sptr.get(), DeallocBoostSharedPtr<T>(sptr), true);
  }
  return null;
}


template<class T>
boost::shared_ptr<T>
Teuchos::shared_pointer( const RCP<T> &rcp )
{
  if (nonnull(rcp)) {
    Ptr<const DeallocBoostSharedPtr<T> >
      dbsp = get_optional_dealloc<DeallocBoostSharedPtr<T> >(rcp);
    if (nonnull(dbsp))
      return dbsp->ptr();
    return boost::shared_ptr<T>(rcp.get(), RCPDeleter<T>(rcp));
  }
  return boost::shared_ptr<T>();
}


#endif	// TEUCHOS_RCP_SHAREDPTR_CONVERSIONS_HPP
