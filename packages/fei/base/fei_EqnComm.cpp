/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include "fei_EqnComm.hpp"
#include "fei_sstream.hpp"

namespace fei {

EqnComm::EqnComm(MPI_Comm comm, int numLocalEqns)
 : comm_(comm),
   globalOffsets_(2)
{
#ifndef FEI_SER

  int numProcs, localProc;
  MPI_Comm_size(comm, &numProcs);
  MPI_Comm_rank(comm, &localProc);

  std::vector<int> local(numProcs*2, 0);
  int* global = &local[0] + numProcs;

  if (numLocalEqns < 0) {
    throw std::runtime_error("fei::EqnComm ERROR, negative numLocalEqns not allowed.");
  }

  local[localProc] = numLocalEqns;

  MPI_Allreduce(&local[0], global, numProcs, MPI_INT, MPI_MAX, comm_);

  globalOffsets_.resize(numProcs+1);

  int offset = 0;
  for(int i=0; i<numProcs; ++i) {
    globalOffsets_[i] = offset;
    offset += global[i];
  }
  globalOffsets_[numProcs] = offset;

#else
  globalOffsets_[0] = 0;
  globalOffsets_[1] = numLocalEqns;
#endif
}
  
EqnComm::EqnComm(MPI_Comm comm, int numLocalEqns, const std::vector<int>& globalOffsets)
 : comm_(comm),
   globalOffsets_(globalOffsets)
{
}
  
EqnComm::~EqnComm()
{
}

const std::vector<int>&
EqnComm::getGlobalOffsets() const
{
  return(globalOffsets_);
}

int
EqnComm::getOwnerProc(int eqn) const
{
//  std::vector<int>::const_iterator
//   iter = std::lower_bound(globalOffsets_.begin(), globalOffsets_.end(),
//                           eqn);
//  int proc = iter - globalOffsets_.begin() - 1;
//  if (*iter==eqn) ++proc;
  int proc = 0;
  for(unsigned p=1; p<globalOffsets_.size(); ++p) {
    if (eqn < globalOffsets_[p]) {
      proc = p-1;
      break;
    }
  }

#ifndef NDEBUG
  int numProcs = globalOffsets_.size()-1;
  if (proc >= numProcs) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::EqnComm::getOwnerProc: input eqn="<<eqn<<", proc="<<proc
      << ", ERROR, proc should be in [0.."<<numProcs-1<<"].";
    throw std::runtime_error(std::string(osstr.str().c_str()));
  }
#endif

  return((int)proc);
}

}//namespace fei

