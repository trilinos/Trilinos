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


#include <fei_macros.hpp>

#include <snl_fei_SubdMsgHandler.hpp>
#include <snl_fei_RecordCollection.hpp>
#include <fei_SharedIDs.hpp>
#include <fei_TemplateUtils.hpp>

#undef fei_file
#define fei_file "snl_fei::SubdMsgHandler"
#include <fei_ErrMacros.hpp>

snl_fei::SubdMsgHandler::SubdMsgHandler(RecordCollection* recordCollection,
					fei::SharedIDs<int>* sharedIDTable,
					fei::SharedIDs<int>* subdomainIDTable)
  :sendPattern_(NULL),
   recvPattern_(NULL),
   recordCollection_(recordCollection),
   sharedIDTable_(sharedIDTable),
   subdomainIDTable_(subdomainIDTable),
   sendProcs_(0, 32),
   recvProcs_(0, 32)
{
}

snl_fei::SubdMsgHandler::~SubdMsgHandler()
{
}

std::vector<int>& snl_fei::SubdMsgHandler::getSendProcs()
{
  fei::copyKeysToVector(sendPattern_->getMap(), sendProcs_);
  return(sendProcs_);
}

std::vector<int>& snl_fei::SubdMsgHandler::getRecvProcs()
{
  fei::copyKeysToVector(recvPattern_->getMap(), recvProcs_);
  return(recvProcs_);
}

int snl_fei::SubdMsgHandler::getSendMessageLength(int destProc,
						  int& messageLength)
{
  fei::comm_map::row_type* list = sendPattern_->getRow(destProc);
  if (list == NULL) {
    FEI_COUT << "SubdMsdHandler: destProc: " << destProc << ", list is NULL."<<FEI_ENDL;
    return(-1);
  }

  messageLength = list->size();

  return( 0 );
}

int snl_fei::SubdMsgHandler::getSendMessage(int destProc,
					  std::vector<int>& message)
{
  int len = sendPattern_->getRow(destProc)->size();
  message.resize(len);
  int* msgPtr = &message[0];

  fei::comm_map::row_type* ids = sendPattern_->getRow(destProc);
  fei::comm_map::row_type::const_iterator
    id_iter = ids->begin(),
    id_end = ids->end();

  int offset = 0;
  for(; id_iter != id_end; ++id_iter) {
    int ID = *id_iter;
    fei::Record<int>* rec = recordCollection_->getRecordWithID(ID);
    if (rec == NULL) {
      ERReturn(-1);
    }

    if (rec->isInLocalSubdomain_) {
      msgPtr[offset++] = 1;
    }
    else {
      msgPtr[offset++] = 0;
    }
  }

  return(0);
}

int snl_fei::SubdMsgHandler::processRecvMessage(int srcProc,
						std::vector<int>& message)
{
  fei::comm_map::row_type* ids = recvPattern_->getRow(srcProc);
  fei::comm_map::row_type::const_iterator
    id_iter = ids->begin(),
    id_end = ids->end();

  int* msgPtr = &message[0];

  if (message.size() != ids->size()) {
    ERReturn(-1);
  }

  int offset = 0;
  for(; id_iter != id_end; ++id_iter) {
    int ID = *id_iter;

    bool isInRemoteSubdomain = msgPtr[offset++] > 1 ? true : false;

    if (isInRemoteSubdomain) {
      subdomainIDTable_->addSharedID(ID, 1, &srcProc);
    }
  }

  return(0);
}

