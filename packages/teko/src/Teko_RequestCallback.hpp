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

#ifndef __Teko_RequestCallback_hpp__
#define __Teko_RequestCallback_hpp__

#include "Teko_RequestMesg.hpp"

namespace Teko {

/** Base class that basically describes a sub
  * set of functonality for a listener.  The application
  * should derive off <code>RequestCallback</code>.
  */
class RequestCallbackBase {
public:
   virtual ~RequestCallbackBase() {}

   /** Does this object satisfy the request?
     * 
     * \returns True if the object can handle the request.
     */ 
   virtual bool handlesRequest(const RequestMesg &) = 0;

   /** Let the object know, that there is a request pending. This
     * permits the application code to ``prepare'' to satisfy the
     * request. This only must be called on construction of the 
     * preconditioner factory. If the result of the request
     * must be constructed multiple times, that is incombant on
     * the developer of the callback.
     *
     * \note No assumption of <code>handlesRequest</code> being true
     *       for this message is made. It is up to the application to
     *       verify this.
     */
   virtual void preRequest(const RequestMesg &) = 0;
};

/** Primary reference point for the application to 
  * derive off for a request callback. This one 
  * has a request function that returns the specified
  * type of data. Note that the type may be a reference
  * counted pointer (<code>Teuchos::RCP</code> for instance)
  * if a persistant value is needed.
  */
template <typename DataT>
class RequestCallback : public RequestCallbackBase {
public:
   virtual ~RequestCallback() {}

   /** Satisfy a request. Similar to <code>preRequest</code>
     * no assumption on the response to <code>handlesRequest</code>
     * is made. Again its up to the user to satisfy this.
     */
   virtual DataT request(const RequestMesg &) = 0;

   /** Does this object satisfy the request?
     * 
     * \returns True if the object can handle the request.
     */ 
   virtual bool handlesRequest(const RequestMesg &) = 0;
};

}

#endif
