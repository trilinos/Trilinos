// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_DistributorActor.hpp"

namespace Tpetra::Details {

DistributorActor::DistributorActor()
  : mpiTag_(DEFAULT_MPI_TAG) {}

void DistributorActor::doWaits(const DistributorPlan& plan) {
  doWaitsRecv(plan);
  doWaitsSend(plan);
}

void DistributorActor::doWaitsRecv(const DistributorPlan& plan) {
  if (requestsRecv_.size() > 0) {
    ProfilingRegion wr("Tpetra::Distributor::doWaitsRecv");

    Teuchos::waitAll(*plan.getComm(), requestsRecv_());

    // Restore the invariant that requests_.size() is the number of
    // outstanding nonblocking communication requests.
    requestsRecv_.resize(0);
  }

  doWaitsIalltofewv(plan);
}

void DistributorActor::doWaitsSend(const DistributorPlan& plan) {
  if (requestsSend_.size() > 0) {
    ProfilingRegion ws("Tpetra::Distributor::doWaitsSend");

    Teuchos::waitAll(*plan.getComm(), requestsSend_());

    // Restore the invariant that requests_.size() is the number of
    // outstanding nonblocking communication requests.
    requestsSend_.resize(0);
  }
}

void DistributorActor::doWaitsIalltofewv(const DistributorPlan& plan) {
#ifdef HAVE_TPETRA_MPI
  if (ialltofewv_.req) {
    ProfilingRegion ws("Tpetra::Distributor::doWaitsIalltofewv");
    ialltofewv_.impl.wait(*ialltofewv_.req);

    ialltofewv_.sendcounts.reset();
    ialltofewv_.sdispls.reset();
    ialltofewv_.recvcounts.reset();
    ialltofewv_.rdispls.reset();
    ialltofewv_.req = std::nullopt;
    ialltofewv_.roots.clear();
  }
#endif
}

bool DistributorActor::isReady() const {
  bool result = true;
  for (auto& request : requestsRecv_) {
    result &= request->isReady();
  }
  for (auto& request : requestsSend_) {
    result &= request->isReady();
  }

  // isReady just calls MPI_Test and returns flag != 0
  // don't use test because these are for a collective, and not
  // all ranks may call test, so progress may not be possible
#ifdef HAVE_TPETRA_MPI
  if (ialltofewv_.req) {
    int flag;
    ialltofewv_.impl.get_status(*ialltofewv_.req, &flag, MPI_STATUS_IGNORE);
    result &= flag;
  }
#endif

  return result;
}
}  // namespace Tpetra::Details
