// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

}  // namespace Teko

#endif
