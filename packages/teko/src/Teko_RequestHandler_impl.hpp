/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// Export of this program may require a license from the United States
// Government.
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#ifndef __Teko_RequestHandler_impl_hpp__
#define __Teko_RequestHandler_impl_hpp__

template <typename DataT>
DataT RequestHandler::request(const RequestMesg & rd) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   // type of call back to use
   typedef RequestCallback<DataT> CallbackT;

   // distribute message over all callbacks
   std::vector<RCP<RequestCallbackBase> >::iterator itr;
   for(itr=callbacks_.begin();itr!=callbacks_.end();++itr) {
      RCP<CallbackT> cb = rcp_dynamic_cast<CallbackT>(*itr);

      // call back not right type
      if(cb==Teuchos::null)
         continue;

      // does call back handle this request?
      if(cb->handlesRequest(rd)) 
         return cb->request(rd);
   }

   TEST_FOR_EXCEPTION(true,std::runtime_error,
         "RequestHandler::request could not find appropriate callback: " << rd);
}

template <typename DataT>
void RequestHandler::preRequest(const RequestMesg & rd) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   // type of call back to use
   typedef RequestCallback<DataT> CallbackT;

   // distribute message over all callbacks
   std::vector<RCP<RequestCallbackBase> >::iterator itr;
   for(itr=callbacks_.begin();itr!=callbacks_.end();++itr) {
      RCP<CallbackT> cb = rcp_dynamic_cast<CallbackT>(*itr);

      // call back not right type
      if(cb==Teuchos::null)
         continue;

      // does call back handle this request?
      if(cb->handlesRequest(rd)) 
         return cb->preRequest(rd);
   }

   TEST_FOR_EXCEPTION(true,std::runtime_error,
         "RequestHandler::preRequest could not find appropriate callback: " << rd);
}

#endif
