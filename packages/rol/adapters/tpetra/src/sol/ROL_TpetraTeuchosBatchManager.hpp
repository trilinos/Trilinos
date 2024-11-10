// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TPETRATEUCHOSBATCHMANAGER_HPP
#define ROL_TPETRATEUCHOSBATCHMANAGER_HPP

#include "ROL_TeuchosBatchManager.hpp"
#include "ROL_TpetraMultiVector.hpp"

namespace ROL {

template<class Real,
         class LO=Tpetra::Map<>::local_ordinal_type,
         class GO=Tpetra::Map<>::global_ordinal_type,
         class Node=Tpetra::Map<>::node_type > 
class TpetraTeuchosBatchManager : public TeuchosBatchManager<Real,GO> {
  typedef Tpetra::MultiVector<Real,LO,GO,Node> Tpetra_Vector;
  typedef TpetraMultiVector<Real,LO,GO,Node> OptVector;

public:
  TpetraTeuchosBatchManager(const ROL::Ptr<const Teuchos::Comm<int> > &comm)
    : TeuchosBatchManager<Real,GO>(comm) {}

  void sumAll(Vector<Real> &input, Vector<Real> &output) {
    ROL::Ptr<Tpetra_Vector> ivec = dynamic_cast<OptVector&>(input).getVector();
    ROL::Ptr<Tpetra_Vector> ovec = dynamic_cast<OptVector&>(output).getVector();

    size_t ilength = ivec->getLocalLength(), olength = ovec->getLocalLength();
    ROL_TEST_FOR_EXCEPTION(ilength != olength, std::invalid_argument,
      ">>> (TpetraTeuchosBatchManager::sumAll): Inconsistent local lengths!");

    size_t invec = ivec->getNumVectors(), onvec = ovec->getNumVectors();
    ROL_TEST_FOR_EXCEPTION(invec != onvec, std::invalid_argument,
      ">>> (TpetraTeuchosBatchManager::sumAll): Inconsistent number of vectors!");

    for (size_t i = 0; i < invec; ++i) {
      TeuchosBatchManager<Real,GO>::sumAll((ivec->getDataNonConst(i)).getRawPtr(),
                                           (ovec->getDataNonConst(i)).getRawPtr(),
                                           ilength);
    }
  }
};

}

#endif
