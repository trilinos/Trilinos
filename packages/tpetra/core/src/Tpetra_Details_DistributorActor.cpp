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

#include "Tpetra_Details_DistributorActor.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Tpetra {
namespace Details {

  DistributorActor::DistributorActor() {
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    makeTimers();
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  }

  DistributorActor::DistributorActor(const DistributorActor& otherActor)
    : requests_(otherActor.requests_)
  {
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    makeTimers();
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  }

  void DistributorActor::doWaits(const DistributorPlan& plan) {
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMon (*timer_doWaits_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    if (requests_.size() > 0) {
      Teuchos::waitAll(*plan.getComm(), requests_());

      // Restore the invariant that requests_.size() is the number of
      // outstanding nonblocking communication requests.
      requests_.resize(0);
    }
  }

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  void DistributorActor::makeTimers () {
    timer_doWaits_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doWaits");

    timer_doPosts3KV_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(3) KV");
    timer_doPosts4KV_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(4) KV");

    timer_doPosts3KV_recvs_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(3): recvs KV");
    timer_doPosts4KV_recvs_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(4): recvs KV");

    timer_doPosts3KV_barrier_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(3): barrier KV");
    timer_doPosts4KV_barrier_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(4): barrier KV");

    timer_doPosts3KV_sends_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(3): sends KV");
    timer_doPosts4KV_sends_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(4): sends KV");
    timer_doPosts3KV_sends_slow_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(3): sends KV SLOW");
    timer_doPosts4KV_sends_slow_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(4): sends KV SLOW");
    timer_doPosts3KV_sends_fast_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(3): sends KV FAST");
    timer_doPosts4KV_sends_fast_ = Teuchos::TimeMonitor::getNewTimer (
                           "Tpetra::Distributor: doPosts(4): sends KV FAST");
  }
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS
}
}
