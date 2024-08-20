// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_DistributorActor.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Tpetra {
namespace Details {

  DistributorActor::DistributorActor()
    : mpiTag_(DEFAULT_MPI_TAG)
  {
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    makeTimers();
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  }

  DistributorActor::DistributorActor(const DistributorActor& otherActor)
    : mpiTag_(otherActor.mpiTag_),
      requests_(otherActor.requests_)
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

  bool DistributorActor::isReady() const {
    bool result = true;
    for (auto& request : requests_) {
      result &= request->isReady();
    }
    return result;
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
