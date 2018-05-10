// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
  TpetraTeuchosBatchManager(const ROL::Ptr<const Teuchos::Comm<GO> > &comm)
    : TeuchosBatchManager<Real,GO>(comm) {}

  void sumAll(Vector<Real> &input, Vector<Real> &output) {
    ROL::Ptr<Tpetra_Vector> ivec = dynamic_cast<OptVector&>(input).getVector();
    ROL::Ptr<Tpetra_Vector> ovec = dynamic_cast<OptVector&>(output).getVector();

    size_t ilength = ivec->getLocalLength(), olength = ovec->getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION(ilength != olength, std::invalid_argument,
      ">>> (TpetraTeuchosBatchManager::sumAll): Inconsistent local lengths!");

    size_t invec = ivec->getNumVectors(), onvec = ovec->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(invec != onvec, std::invalid_argument,
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
