// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROLSINGLETONTEUCHOSBATCHMANAGER_HPP
#define ROLSINGLETONTEUCHOSBATCHMANAGER_HPP

#include "ROL_TeuchosBatchManager.hpp"
#include "ROL_SingletonVector.hpp"

namespace ROL {

template<class Real, class Ordinal>
class SingletonTeuchosBatchManager : public TeuchosBatchManager<Real,Ordinal> {
public:
  SingletonTeuchosBatchManager(const ROL::Ptr<const Teuchos::Comm<int> > &comm)
    : TeuchosBatchManager<Real,Ordinal>(comm) {}

  using TeuchosBatchManager<Real,Ordinal>::sumAll;

  void sumAll(Vector<Real> &input, Vector<Real> &output) {
    Real inval = dynamic_cast<SingletonVector<Real>&>(input).getValue();
    Real outval(0);
    TeuchosBatchManager<Real,Ordinal>::sumAll(&inval,&outval,1);
    dynamic_cast<SingletonVector<Real>&>(output).setValue(outval);
  }
};

}

#endif
