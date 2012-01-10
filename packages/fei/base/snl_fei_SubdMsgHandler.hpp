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


#ifndef _snl_fei_SubdMsgHandler_hpp_
#define _snl_fei_SubdMsgHandler_hpp_

#include <fei_macros.hpp>
#include <fei_CommUtils.hpp>
#include <fei_fwd.hpp>

#include <vector>

namespace fei {
  template<typename T> class SharedIDs;
}

namespace snl_fei {
  /** implementation of MessageHandler for subdomain data */
  class SubdMsgHandler : public fei::MessageHandler<int> {
  public:
    /** constructor */
    SubdMsgHandler(RecordCollection* recordCollection,
		   fei::SharedIDs<int>* sharedIDTable,
		   fei::SharedIDs<int>* subdomainIDTable);
    /** destructor */
    virtual ~SubdMsgHandler();

    /** get list of processors to be sent to */
    std::vector<int>& getSendProcs();

    /** get list of processors to be recvd from */
    std::vector<int>& getRecvProcs();

    /** get length of message to be sent to specified proc */
    int getSendMessageLength(int destProc, int& messageLength);

    /** get message to be sent to specified proc */
    int getSendMessage(int destProc, std::vector<int>& message);

    /** process message received from specified recv proc */
    int processRecvMessage(int srcProc, std::vector<int>& message);

    /** set pattern describing procs to be sent to */
    void setSendPattern(fei::comm_map* pattern)
      { sendPattern_ = pattern; }

    /** set pattern describing procs to be recvd from */
    void setRecvPattern(fei::comm_map* pattern)
      { recvPattern_ = pattern; }

  private:
    fei::comm_map* sendPattern_;
    fei::comm_map* recvPattern_;
    RecordCollection* recordCollection_;
    fei::SharedIDs<int>* sharedIDTable_;
    fei::SharedIDs<int>* subdomainIDTable_;

    std::vector<int> sendProcs_;
    std::vector<int> recvProcs_;
  };//class SubdMsgHandler
}//namespace snl_fei

#endif // _snl_fei_SubdMsgHandler_hpp_

