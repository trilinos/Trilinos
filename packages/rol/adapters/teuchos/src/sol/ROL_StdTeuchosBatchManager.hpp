// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROLSTDTEUCHOSBATCHMANAGER_HPP
#define ROLSTDTEUCHOSBATCHMANAGER_HPP

#include "ROL_TeuchosBatchManager.hpp"
#include "ROL_StdVector.hpp"

namespace ROL {

template<class Real, class Ordinal>
class StdTeuchosBatchManager : public TeuchosBatchManager<Real,Ordinal> {
public:
  StdTeuchosBatchManager(const ROL::Ptr<const Teuchos::Comm<int>> &comm)
    : TeuchosBatchManager<Real,Ordinal>(comm) {}

  using TeuchosBatchManager<Real,Ordinal>::sumAll;
  void sumAll(Vector<Real> &input, Vector<Real> &output) {
    std::vector<Real> &idata = *dynamic_cast<StdVector<Real>&>(input).getVector();
    std::vector<Real> &odata = *dynamic_cast<StdVector<Real>&>(output).getVector();
    int size = idata.size();
    ROL_TEST_FOR_EXCEPTION(size != static_cast<int>(odata.size()), std::invalid_argument,
      ">>> (ROL::StdTeuchosBatchManager::SumAll): Dimension mismatch!");
    TeuchosBatchManager<Real,Ordinal>::sumAll(&idata[0],&odata[0],size);
  }
};

}

#endif
