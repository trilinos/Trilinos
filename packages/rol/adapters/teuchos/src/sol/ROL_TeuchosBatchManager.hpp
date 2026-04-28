// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  const ROL::Ptr<const Teuchos::Comm<int> > comm_;

public:
  TeuchosBatchManager(const ROL::Ptr<const Teuchos::Comm<int> > &comm)
    : comm_(comm) {}

  int batchID(void) {
    return Teuchos::rank<int>(*comm_);
  }

  int numBatches(void) {
    return Teuchos::size<int>(*comm_);
  }

  void reduceAll(Real* input, Real* output, int dim,
                 const Elementwise::ReductionOp<Real> &r) {
    int nB = this->numBatches();
    std::vector<Real> receiveBuffer(nB);
    Teuchos::gather<int,Real>(input,1,&receiveBuffer[0],1,0,*comm_);
    output[0] = r.initialValue();
    for (int i = 0; i < nB; i++) {
      r.reduce(receiveBuffer[i],output[0]);
    }
    Teuchos::broadcast<int,Real>(*comm_,0,1,output);
  }

  void minAll(Real* input, Real* output, int dim) {
    Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_MIN,
      dim, input, output);
  }

  void maxAll(Real* input, Real* output, int dim) {
    Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_MAX,
      dim, input, output);
  }

  void sumAll(Real* input, Real* output, int dim) {
    Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,
      dim, input, output);
  }

  void gatherAll(const Real* send, const int ssize, Real *receive, const int rsize) const {
    Teuchos::gatherAll<int,Real>(*comm_,ssize,send,rsize,receive);
  }

  void gather(const Real* send, int ssize, Real *receive, int rsize, int root) const {
    std::vector<int> rcnt, disp;
    if (comm_->getRank()==root) rcnt.resize(comm_->getSize());
    Teuchos::gather<int,int>(&ssize,1,&rcnt[0],1,root,*comm_);
    if (comm_->getRank()==root) {
      disp.resize(rcnt.size(),0);
      for (unsigned i=1u; i<rcnt.size(); ++i) disp[i] = disp[i-1]+rcnt[i-1];
    }
    Teuchos::gatherv<int,Real>(send,ssize,receive,rcnt.data(),disp.data(),root,*comm_);
  }

  void broadcast(Real* input, int cnt, int root) {
    Teuchos::broadcast<int,Real>(*comm_,root,cnt,input);
  }

  virtual void sumAll(Vector<Real> &input, Vector<Real> &output) {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::TeuchosBatchManager): sumAll(Vector<Real> &input, Vector<Real> &output) is not implemented");
  }

  void barrier(void) {
    Teuchos::barrier<int>(*comm_); 
  }
};

}

#endif
