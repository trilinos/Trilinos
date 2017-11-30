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

#ifndef ROLTEUCHOSBATCHMANAGER_HPP
#define ROLTEUCHOSBATCHMANAGER_HPP

#include "ROL_BatchManager.hpp"
#include "ROL_Ptr.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace ROL {

template<class Real, class Ordinal>
class TeuchosBatchManager : public BatchManager<Real> {
private:
  const ROL::Ptr<const Teuchos::Comm<Ordinal> > comm_;

public:
  TeuchosBatchManager(const ROL::Ptr<const Teuchos::Comm<Ordinal> > &comm)
    : comm_(comm) {}

  int batchID(void) {
    return Teuchos::rank<Ordinal>(*comm_);
  }

  int numBatches(void) {
    return Teuchos::size<Ordinal>(*comm_);
  }

  void reduceAll(Real* input, Real* output, int dim,
                 const Elementwise::ReductionOp<Real> &r) {
    int nB = this->numBatches();
    std::vector<Real> receiveBuffer(nB);
    Teuchos::gather<Ordinal,Real>(input,1,&receiveBuffer[0],1,0,*comm_);
    output[0] = r.initialValue();
    for (int i = 0; i < nB; i++) {
      r.reduce(receiveBuffer[i],output[0]);
    }
    Teuchos::broadcast<Ordinal,Real>(*comm_,0,1,output);
  }

  void minAll(Real* input, Real* output, int dim) {
    Teuchos::reduceAll<Ordinal,Real>(*comm_,Teuchos::REDUCE_MIN,
      dim, input, output);
  }

  void maxAll(Real* input, Real* output, int dim) {
    Teuchos::reduceAll<Ordinal,Real>(*comm_,Teuchos::REDUCE_MAX,
      dim, input, output);
  }

  void sumAll(Real* input, Real* output, int dim) {
    Teuchos::reduceAll<Ordinal,Real>(*comm_,Teuchos::REDUCE_SUM,
      dim, input, output);
  }

  void broadcast(Real* input, int cnt, int root) {
    Teuchos::broadcast<Ordinal,Real>(*comm_,root,cnt,input);
  }

  virtual void sumAll(Vector<Real> &input, Vector<Real> &output) {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::TeuchosBatchManager): sumAll(Vector<Real> &input, Vector<Real> &output) is not implemented");
  }

  void barrier(void) {
    Teuchos::barrier<Ordinal>(*comm_); 
  }
};

}

#endif
