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

  DistributorActor::DistributorActor(const DistributorActor& otherActor)
    : mpiTag_(otherActor.mpiTag_),
    requestsRecv_(otherActor.requestsRecv_),
    requestsSend_(otherActor.requestsSend_) {}

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

  bool DistributorActor::isReady() const {
    bool result = true;
    for (auto& request : requestsRecv_) {
      result &= request->isReady();
    }
    for (auto& request : requestsSend_) {
      result &= request->isReady();
    }
    return result;
  }
}
