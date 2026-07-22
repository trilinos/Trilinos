// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_RequestHandlerContainer_hpp__
#define __Teko_RequestHandlerContainer_hpp__

#include "Teko_RequestHandler.hpp"

namespace Teko {

/** Pure virtual interface class that
 * provides a mechanism for setting the
 * request handler.
 */
class RequestHandlerContainer {
 public:
  virtual ~RequestHandlerContainer() {}

  //! Set the request handler with pointers to the appropriate callbacks
  virtual void setRequestHandler(const Teuchos::RCP<RequestHandler>& rh) = 0;

  //! Get the request handler with pointers to the appropriate callbacks
  virtual Teuchos::RCP<RequestHandler> getRequestHandler() const = 0;
};

}  // end namespace Teko

#endif
