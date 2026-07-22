// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_RequestHandler.hpp"

namespace Teko {

RequestHandler::RequestHandler() {}

void RequestHandler::addRequestCallback(const Teuchos::RCP<RequestCallbackBase>& callback) {
  callbacks_.push_back(callback);
}

}  // end namespace Teko
