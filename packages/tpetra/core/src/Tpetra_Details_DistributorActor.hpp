// @HEADER
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

#ifndef TPETRA_DETAILS_DISTRIBUTOR_ACTOR_HPP
#define TPETRA_DETAILS_DISTRIBUTOR_ACTOR_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Time.hpp"

namespace Tpetra {
namespace Details {

class DistributorActor {

public:
  DistributorActor();
  DistributorActor(const DistributorActor& otherActor);

  /// \brief Communication requests associated with nonblocking receives and sends.
  ///
  /// \note To implementers: Distributor uses requests_.size() as
  ///   the number of outstanding nonblocking receives and sends.
  ///   This means you should always resize to zero after completing
  ///   receive and send requests.
  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest<int>>> requests_;

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::RCP<Teuchos::Time> timer_doPosts3TA_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4TA_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3KV_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4KV_;
  Teuchos::RCP<Teuchos::Time> timer_doWaits_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3TA_recvs_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4TA_recvs_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3TA_barrier_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4TA_barrier_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3TA_sends_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4TA_sends_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3TA_sends_slow_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4TA_sends_slow_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3TA_sends_fast_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4TA_sends_fast_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3KV_recvs_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4KV_recvs_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3KV_barrier_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4KV_barrier_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3KV_sends_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4KV_sends_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3KV_sends_slow_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4KV_sends_slow_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts3KV_sends_fast_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4KV_sends_fast_;

  //! Make the instance's timers.  (Call only in constructor.)
  void makeTimers();
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS
};

}
}

#endif
