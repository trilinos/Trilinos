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

#ifndef ROL_EPETRABATCHMANAGER_HPP
#define ROL_EPETRABATCHMANAGER_HPP

#include "ROL_BatchManager.hpp"
#include "Epetra_Comm.h"

namespace ROL {

template<class Real> 
class EpetraBatchManager : public BatchManager<Real> {
protected: 
  Teuchos::RCP<Epetra_Comm> comm_;

public:
  virtual ~EpetraBatchManager() {}

  EpetraBatchManager(Teuchos::RCP<Epetra_Comm> &comm) : comm_(comm) {}

  int batchID(void) { return comm_->MyPID(); }

  int numBatches(void) { return comm_->NumProc(); }

  void reduceAll(Real* input, Real* output,
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
