// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#include "Tpetra_Details_DistributorPlan.hpp"

namespace Tpetra {
namespace Details {

std::string
DistributorSendTypeEnumToString (EDistributorSendType sendType)
{
  if (sendType == DISTRIBUTOR_ISEND) {
    return "Isend";
  }
  else if (sendType == DISTRIBUTOR_RSEND) {
    return "Rsend";
  }
  else if (sendType == DISTRIBUTOR_SEND) {
    return "Send";
  }
  else if (sendType == DISTRIBUTOR_SSEND) {
    return "Ssend";
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Invalid "
      "EDistributorSendType enum value " << sendType << ".");
  }
}

std::string
DistributorHowInitializedEnumToString (EDistributorHowInitialized how)
{
  switch (how) {
  case Details::DISTRIBUTOR_NOT_INITIALIZED:
    return "Not initialized yet";
  case Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS:
    return "By createFromSends";
  case Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_RECVS:
    return "By createFromRecvs";
  case Details::DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS_N_RECVS:
    return "By createFromSendsAndRecvs";
  case Details::DISTRIBUTOR_INITIALIZED_BY_REVERSE:
    return "By createReverseDistributor";
  case Details::DISTRIBUTOR_INITIALIZED_BY_COPY:
    return "By copy constructor";
  default:
    return "INVALID";
  }
}

DistributorPlan::DistributorPlan(Teuchos::RCP<const Teuchos::Comm<int>> comm)
  : comm_(comm),
    howInitialized_(DISTRIBUTOR_NOT_INITIALIZED),
    sendType_(DISTRIBUTOR_SEND),
    barrierBetweenRecvSend_(barrierBetween_default),
    useDistinctTags_(useDistinctTags_default),
    sendMessageToSelf_(false),
    numSendsToOtherProcs_(0),
    maxSendLength_(0),
    numReceives_(0),
    totalReceiveLength_(0)
{ }

DistributorPlan::DistributorPlan(const DistributorPlan& otherPlan)
  : comm_(otherPlan.comm_),
    howInitialized_(DISTRIBUTOR_INITIALIZED_BY_COPY),
    sendType_(otherPlan.sendType_),
    barrierBetweenRecvSend_(otherPlan.barrierBetweenRecvSend_),
    useDistinctTags_(otherPlan.useDistinctTags_),
    sendMessageToSelf_(otherPlan.sendMessageToSelf_),
    numSendsToOtherProcs_(otherPlan.numSendsToOtherProcs_),
    procIdsToSendTo_(otherPlan.procIdsToSendTo_),
    startsTo_(otherPlan.startsTo_),
    lengthsTo_(otherPlan.lengthsTo_),
    maxSendLength_(otherPlan.maxSendLength_),
    indicesTo_(otherPlan.indicesTo_),
    numReceives_(otherPlan.numReceives_),
    totalReceiveLength_(otherPlan.totalReceiveLength_),
    lengthsFrom_(otherPlan.lengthsFrom_),
    procsFrom_(otherPlan.procsFrom_),
    startsFrom_(otherPlan.startsFrom_),
    indicesFrom_(otherPlan.indicesFrom_)
{ }

}
}
