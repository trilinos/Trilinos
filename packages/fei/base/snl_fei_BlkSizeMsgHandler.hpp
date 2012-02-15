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


#ifndef _snl_fei_BlkSizeMsgHandler_hpp_
#define _snl_fei_BlkSizeMsgHandler_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Graph.hpp>
#include <fei_CommUtils.hpp>

namespace snl_fei {

/** MessageHandler implementation for block-size data. */
class BlkSizeMsgHandler : public fei::MessageHandler<int> {
 public:
  /** constructor */
  BlkSizeMsgHandler(fei::VectorSpace* vspace,
		    fei::Graph* graph,
		    MPI_Comm comm);
  /** destructor */
  virtual ~BlkSizeMsgHandler();

  /** clumsy method to launch the data exchange. */
  int do_the_exchange();

  /** Get list of procs to send to. */
  std::vector<int>& getSendProcs();
  /** Get list of procs to recv from. */
  std::vector<int>& getRecvProcs();

  /** Get length of message for specified destination proc. */
  int getSendMessageLength(int destProc, int& messageLength);
  /** Get message to send to specified destination proc. */
  int getSendMessage(int destProc, std::vector<int>& message);
  /** process message received from specified source proc. */
  int processRecvMessage(int srcProc, std::vector<int>& message);

 private:
  fei::comm_map* remote_colIndices_;
  fei::comm_map* local_colIndices_;
  fei::VectorSpace* vecSpace_;
  snl_fei::PointBlockMap* ptBlkMap_;
  fei::Graph* graph_;
  MPI_Comm comm_;
  std::vector<int> sendProcs_;
  std::vector<int> recvProcs_;

  bool firstExchange_;
};

} // namespace snl_fei

#endif

