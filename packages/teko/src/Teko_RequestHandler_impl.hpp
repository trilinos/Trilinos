// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_RequestHandler_impl_hpp__
#define __Teko_RequestHandler_impl_hpp__

template <typename DataT>
DataT RequestHandler::request(const RequestMesg& rd) const {
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // type of call back to use
  typedef RequestCallback<DataT> CallbackT;

  // distribute message over all callbacks
  std::vector<RCP<RequestCallbackBase> >::iterator itr;
  for (itr = callbacks_.begin(); itr != callbacks_.end(); ++itr) {
    RCP<CallbackT> cb = rcp_dynamic_cast<CallbackT>(*itr);

    // call back not right type
    if (cb == Teuchos::null) continue;

    // does call back handle this request?
    if (cb->handlesRequest(rd)) return cb->request(rd);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                             "RequestHandler::request could not find appropriate callback: " << rd);
}

template <typename DataT>
void RequestHandler::preRequest(const RequestMesg& rd) const {
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // type of call back to use
  typedef RequestCallback<DataT> CallbackT;

  // distribute message over all callbacks
  std::vector<RCP<RequestCallbackBase> >::iterator itr;
  for (itr = callbacks_.begin(); itr != callbacks_.end(); ++itr) {
    RCP<CallbackT> cb = rcp_dynamic_cast<CallbackT>(*itr);

    // call back not right type
    if (cb == Teuchos::null) continue;

    // does call back handle this request?
    if (cb->handlesRequest(rd)) return cb->preRequest(rd);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error,
      "RequestHandler::preRequest could not find appropriate callback: " << rd);
}

#endif
