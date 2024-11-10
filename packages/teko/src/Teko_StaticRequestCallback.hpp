// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_StaticRequestCallback_hpp__
#define __Teko_StaticRequestCallback_hpp__

#include "Teko_RequestCallback.hpp"

namespace Teko {

/** A simple request interface that takes
 * a static bit of data to return to Teko.
 * This is meant primarily as a testing and
 * early stage development tool.
 *
 * The constructor takes an object of the same
 * type the class was templated on. It also
 * takes a string to be used to match the
 * request.
 */
template <typename DataT>
class StaticRequestCallback : public RequestCallback<DataT> {
 public:
  StaticRequestCallback(const std::string& name, DataT data) : name_(name), data_(data) {}

  DataT request(const RequestMesg& rm);
  void preRequest(const RequestMesg& rm);
  bool handlesRequest(const RequestMesg& rm);

 private:
  std::string name_;
  DataT data_;
};

template <typename DataT>
DataT StaticRequestCallback<DataT>::request(const RequestMesg& rm) {
  TEUCHOS_ASSERT(handlesRequest(rm));

  return data_;
}

template <typename DataT>
void StaticRequestCallback<DataT>::preRequest(const RequestMesg& rm) {
  TEUCHOS_ASSERT(handlesRequest(rm));
}

template <typename DataT>
bool StaticRequestCallback<DataT>::handlesRequest(const RequestMesg& rm) {
  return name_ == rm.getName();
}

}  // namespace Teko

#endif
