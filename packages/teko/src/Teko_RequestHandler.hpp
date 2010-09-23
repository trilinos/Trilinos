/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#ifndef __Teko_RequestHandler_hpp__
#define __Teko_RequestHandler_hpp__

#include <vector>

#include "Teko_RequestCallback.hpp"

#include "Teuchos_RCP.hpp"

namespace Teko {

/** Classes that handles and distrubutes requests.
  * This has two types of users. Those that register
  * callbacks to handle the requests, and those that
  * make requests.
  *
  * This class is passive-aggressive.  It calls required
  * data a "request", however it actually requires 
  * a handler to satisfy the request.
  *
  * \note The implemntation is based on the listener pattern.
  */
class RequestHandler {
public:
   explicit RequestHandler();

   /** Add a call back object to handle requests
     *
     * \param[in] callback Ref-count-pointer to a call back object.
     */
   void addRequestCallback(const Teuchos::RCP<RequestCallbackBase> & callback);

   /** Get the data for a particular request.
     *
     * \param[in] rm The message describing the request.
     */
   template <typename DataT>
   DataT request(const RequestMesg & rm) const;

   /** Get the data for a particular request. Short hand version.
     *
     * \param[in] rm The message describing the request.
     */
   template <typename DataT>
   inline DataT request(const std::string & rm) const
   { return request<DataT>(RequestMesg(rm)); }

   /** Send a pre-request message to the callback
     * allowing them to do some work up front and
     * ahead of time. This is meant to be called
     * at construction/initialization.
     *
     * \param[in] rm The message describing the request.
     */
   template <typename DataT>
   void preRequest(const RequestMesg & rm) const;

   /** Send a pre-request message to the callback
     * allowing them to do some work up front and
     * ahead of time. This is meant to be called
     * at construction/initialization. Short hand version.
     *
     * \param[in] rm The message describing the request.
     */
   template <typename DataT>
   inline void preRequest(const std::string & rm) const
   { preRequest<DataT>(RequestMesg(rm)); }

private:
   // stores the callbacks to be used by this handler.
   mutable std::vector<Teuchos::RCP<RequestCallbackBase> > callbacks_;

   // hidden from the user
   RequestHandler(const RequestHandler & rh);
   RequestHandler & operator=(const RequestHandler &);
   const RequestHandler & operator=(const RequestHandler &) const;
};

#include "Teko_RequestHandler_impl.hpp"

}

#endif
