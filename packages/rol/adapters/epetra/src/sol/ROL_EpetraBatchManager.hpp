// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_EPETRABATCHMANAGER_HPP
#define ROL_EPETRABATCHMANAGER_HPP

#include "ROL_BatchManager.hpp"
#include "Epetra_Comm.h"

namespace ROL {

template<class Real> 
class EpetraBatchManager : public BatchManager<Real> {
protected: 
  ROL::Ptr<Epetra_Comm> comm_;

public:
  virtual ~EpetraBatchManager() {}

  EpetraBatchManager(ROL::Ptr<Epetra_Comm> &comm) : comm_(comm) {}

  int batchID(void) { return comm_->MyPID(); }

  int numBatches(void) { return comm_->NumProc(); }

  void reduceAll(Real* input, Real* output, int dim,
                 const Elementwise::ReductionOp<Real> &r) {
    int nB = this->numBatches();
    std::vector<Real> receiveBuffer(nB);
    comm_->GatherAll(input,&receiveBuffer[0],1);
    output[0] = r.initialValue();
    for (int i = 0; i < nB; i++) {
      r.reduce(receiveBuffer[i],output[0]);
    }
    comm_->Broadcast(output,1,0);
  }

  void sumAll(Real* input, Real* output, int dim) { comm_->SumAll(input,output,dim); }

  void broadcast(Real* input, int cnt, int root) {
    comm_->Broadcast(input, cnt, root);
  }

  void barrier(void) {
    comm_->Barrier();
  }
};

}

#endif
