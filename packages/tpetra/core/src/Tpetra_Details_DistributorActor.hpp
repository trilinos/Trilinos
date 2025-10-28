// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_DISTRIBUTOR_ACTOR_HPP
#define TPETRA_DETAILS_DISTRIBUTOR_ACTOR_HPP

#include <optional>
#include <vector>

#include "Teuchos_Assert.hpp"
#include "Tpetra_Details_DistributorPlan.hpp"
#include "Tpetra_Util.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include "Tpetra_Details_MpiTypeTraits.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Teuchos_RCP.hpp"

#include "Kokkos_TeuchosCommAdapters.hpp"

#ifdef HAVE_TPETRA_MPI
#include "mpi.h"
#include "Tpetra_Details_Ialltofewv.hpp"
#endif

namespace Tpetra::Details {

template <class View>
constexpr bool isKokkosView = Kokkos::is_view<View>::value;

template <class View1, class View2>
constexpr bool areKokkosViews = Kokkos::is_view<View1>::value&& Kokkos::is_view<View2>::value;

class DistributorActor {
  using IndexView     = DistributorPlan::IndexView;
  using SubViewLimits = DistributorPlan::SubViewLimits;

 public:
  static constexpr int DEFAULT_MPI_TAG = 1;

  DistributorActor();
  DistributorActor(const DistributorActor& otherActor) = default;

  template <class ExpView, class ImpView>
  void doPostsAndWaits(const DistributorPlan& plan,
                       const ExpView& exports,
                       size_t numPackets,
                       const ImpView& imports);

  template <class ExpView, class ImpView>
  void doPostsAndWaits(const DistributorPlan& plan,
                       const ExpView& exports,
                       const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                       const ImpView& imports,
                       const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

  template <class ImpView>
  void doPostRecvs(const DistributorPlan& plan,
                   size_t numPackets,
                   const ImpView& imports);

  template <class ImpView>
  void doPostRecvs(const DistributorPlan& plan,
                   const ImpView& imports,
                   const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

  template <class ExpView, class ImpView>
  void doPostSends(const DistributorPlan& plan,
                   const ExpView& exports,
                   size_t numPackets,
                   const ImpView& imports);

  template <class ExpView, class ImpView>
  void doPostSends(const DistributorPlan& plan,
                   const ExpView& exports,
                   const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                   const ImpView& imports,
                   const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

  template <class ExpView, class ImpView>
  void doPosts(const DistributorPlan& plan,
               const ExpView& exports,
               size_t numPackets,
               const ImpView& imports);

  template <class ExpView, class ImpView>
  void doPosts(const DistributorPlan& plan,
               const ExpView& exports,
               const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
               const ImpView& imports,
               const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

  void doWaits(const DistributorPlan& plan);

  void doWaitsRecv(const DistributorPlan& plan);

  void doWaitsSend(const DistributorPlan& plan);

  void doWaitsIalltofewv(const DistributorPlan& plan);

  bool isReady() const;

  int getMpiTag() const { return mpiTag_; };

 private:
  template <class ImpView>
  void doPostRecvsImpl(const DistributorPlan& plan,
                       const ImpView& imports,
                       const SubViewLimits& totalPacketsFrom);

  template <class ExpView, class ImpView>
  void doPostSendsImpl(const DistributorPlan& plan,
                       const ExpView& exports,
                       const SubViewLimits& exportSubViewLimits,
                       const ImpView& imports,
                       const SubViewLimits& importSubViewLimits);

#ifdef HAVE_TPETRA_MPI
  template <class ExpView, class ImpView>
  void doPostsAllToAllImpl(const DistributorPlan& plan,
                           const ExpView& exports,
                           const SubViewLimits& exportSubViewLimits,
                           const ImpView& imports,
                           const SubViewLimits& importSubViewLimits);

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  template <class ExpView, class ImpView>
  void doPostsNbrAllToAllVImpl(const DistributorPlan& plan,
                               const ExpView& exports,
                               const SubViewLimits& exportSubViewLimits,
                               const ImpView& imports,
                               const SubViewLimits& importSubViewLimits);
#endif  // HAVE_TPETRACORE_MPI_ADVANCE

  template <typename ExpView, typename ImpView>
  void doPostsIalltofewvImpl(const DistributorPlan& plan,
                             const ExpView& exports,
                             const SubViewLimits& exportSubViewLimits,
                             const ImpView& imports,
                             const SubViewLimits& importSubViewLimits);

  // ialltofewv members
  struct {
    Details::Ialltofewv impl;
    std::optional<Details::Ialltofewv::Req> req;
    Teuchos::RCP<std::vector<int>> sendcounts;
    Teuchos::RCP<std::vector<int>> sdispls;
    Teuchos::RCP<std::vector<int>> recvcounts;
    Teuchos::RCP<std::vector<int>> rdispls;
    std::vector<int> roots;
  } ialltofewv_;

#endif  // HAVE_TPETRA_MPI

  int mpiTag_;

  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest<int>>> requestsRecv_;
  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest<int>>> requestsSend_;
};

template <class ExpView, class ImpView>
void DistributorActor::doPosts(const DistributorPlan& plan,
                               const ExpView& exports,
                               size_t numPackets,
                               const ImpView& imports) {
  doPostRecvs(plan, numPackets, imports);
  doPostSends(plan, exports, numPackets, imports);
}

template <class ExpView, class ImpView>
void DistributorActor::doPostsAndWaits(const DistributorPlan& plan,
                                       const ExpView& exports,
                                       size_t numPackets,
                                       const ImpView& imports) {
  static_assert(areKokkosViews<ExpView, ImpView>,
                "Data arrays for DistributorActor::doPostsAndWaits must be Kokkos::Views");
  doPosts(plan, exports, numPackets, imports);
  doWaits(plan);
}

template <class ExpView, class ImpView>
void DistributorActor::doPosts(const DistributorPlan& plan,
                               const ExpView& exports,
                               const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                               const ImpView& imports,
                               const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID) {
  doPostRecvs(plan, imports, numImportPacketsPerLID);
  doPostSends(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
}

template <class ExpView, class ImpView>
void DistributorActor::doPostsAndWaits(const DistributorPlan& plan,
                                       const ExpView& exports,
                                       const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                                       const ImpView& imports,
                                       const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID) {
  static_assert(areKokkosViews<ExpView, ImpView>,
                "Data arrays for DistributorActor::doPostsAndWaits must be Kokkos::Views");
  doPosts(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
  doWaits(plan);
}

template <typename ViewType>
using HostAccessibility = Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, typename ViewType::memory_space>;

template <typename DstViewType, typename SrcViewType>
using enableIfHostAccessible = std::enable_if_t<HostAccessibility<DstViewType>::accessible &&
                                                HostAccessibility<SrcViewType>::accessible>;

template <typename DstViewType, typename SrcViewType>
using enableIfNotHostAccessible = std::enable_if_t<!HostAccessibility<DstViewType>::accessible ||
                                                   !HostAccessibility<SrcViewType>::accessible>;

template <typename DstViewType, typename SrcViewType>
enableIfHostAccessible<DstViewType, SrcViewType>
packOffset(const DstViewType& dst,
           const SrcViewType& src,
           const size_t dst_offset,
           const size_t src_offset,
           const size_t size) {
  memcpy((void*)(dst.data() + dst_offset), src.data() + src_offset, size * sizeof(typename DstViewType::value_type));
}

template <typename DstViewType, typename SrcViewType>
enableIfNotHostAccessible<DstViewType, SrcViewType>
packOffset(const DstViewType& dst,
           const SrcViewType& src,
           const size_t dst_offset,
           const size_t src_offset,
           const size_t size) {
  Kokkos::Compat::deep_copy_offset(dst, src, dst_offset, src_offset, size);
}

#ifdef HAVE_TPETRA_MPI

template <class ExpView, class ImpView>
void DistributorActor::doPostsIalltofewvImpl(const DistributorPlan& plan,
                                             const ExpView& exports,
                                             const SubViewLimits& exportSubViewLimits,
                                             const ImpView& imports,
                                             const SubViewLimits& importSubViewLimits) {
  using size_type   = Teuchos::Array<size_t>::size_type;
  using ExportValue = typename ExpView::non_const_value_type;

  ProfilingRegion pr("Tpetra::Distributor::doPostsIalltofewvImpl");

  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Ialltofewv\" only works for fast-path communication.");

  TEUCHOS_TEST_FOR_EXCEPTION(
      bool(ialltofewv_.req), std::runtime_error,
      "This actor has an active Ialltofewv already");

  TEUCHOS_TEST_FOR_EXCEPTION(
      bool(ialltofewv_.sendcounts), std::runtime_error,
      "This actor has an active Ialltofewv already");

  TEUCHOS_TEST_FOR_EXCEPTION(
      bool(ialltofewv_.sdispls), std::runtime_error,
      "This actor has an active Ialltofewv already");

  TEUCHOS_TEST_FOR_EXCEPTION(
      bool(ialltofewv_.recvcounts), std::runtime_error,
      "This actor has an active Ialltofewv already");

  TEUCHOS_TEST_FOR_EXCEPTION(
      bool(ialltofewv_.rdispls), std::runtime_error,
      "This actor has an active Ialltofewv already");

  auto comm = plan.getComm();

  const auto& [importStarts, importLengths] = importSubViewLimits;
  const auto& [exportStarts, exportLengths] = exportSubViewLimits;

  ialltofewv_.roots      = plan.getRoots();
  const int nroots       = ialltofewv_.roots.size();
  const int* roots       = ialltofewv_.roots.data();
  ialltofewv_.req        = std::make_optional<Details::Ialltofewv::Req>();
  ialltofewv_.sendcounts = Teuchos::rcp(new std::vector<int>(nroots));
  ialltofewv_.sdispls    = Teuchos::rcp(new std::vector<int>(nroots));
  ialltofewv_.recvcounts = Teuchos::rcp(new std::vector<int>);
  ialltofewv_.rdispls    = Teuchos::rcp(new std::vector<int>);

  for (int rootIdx = 0; rootIdx < nroots; ++rootIdx) {
    const int root = roots[rootIdx];

    // if we can't find the root proc index in our plan, it just means we send 0
    //  also make sure root only appears once in getProcsTo()
    size_type rootProcIndex = plan.getProcsTo().size();  // sentinel value -> not found
    for (size_type pi = 0; pi < plan.getProcsTo().size(); ++pi) {
      if (plan.getProcsTo()[pi] == root) {
        rootProcIndex = pi;
        break;
      }
    }

    // am I sending to root?
    int sendcount = 0;
    if (rootProcIndex != plan.getProcsTo().size()) {
      sendcount = exportLengths[rootProcIndex];
    }
    (*ialltofewv_.sendcounts)[rootIdx] = sendcount;

    int sdispl = 0;
    if (0 != sendcount) {
      sdispl = exportStarts[rootProcIndex];
    }
    (*ialltofewv_.sdispls)[rootIdx] = sdispl;

    // If I happen to be this root, set up my receive metadata
    if (comm->getRank() == root) {
      // don't recv anything from anywhere by default
      ialltofewv_.recvcounts->resize(comm->getSize());
      std::fill(ialltofewv_.recvcounts->begin(), ialltofewv_.recvcounts->end(), 0);
      ialltofewv_.rdispls->resize(comm->getSize());
      std::fill(ialltofewv_.rdispls->begin(), ialltofewv_.rdispls->end(), 0);

      const size_type actualNumReceives =
          Teuchos::as<size_type>(plan.getNumReceives()) +
          Teuchos::as<size_type>(plan.hasSelfMessage() ? 1 : 0);

      for (size_type i = 0; i < actualNumReceives; ++i) {
        const int src                  = plan.getProcsFrom()[i];
        (*ialltofewv_.rdispls)[src]    = importStarts[i];
        (*ialltofewv_.recvcounts)[src] = Teuchos::as<int>(importLengths[i]);
      }
    }

  }  // rootIdx

  // TODO: do we need to pass ExportValue{} here?
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<ExportValue>::getType(ExportValue{});
  // FIXME: is there a better way to do this?
  Teuchos::RCP<const Teuchos::MpiComm<int>> tMpiComm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm);
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> oMpiComm =
      tMpiComm->getRawMpiComm();
  MPI_Comm mpiComm = (*oMpiComm)();

  // don't care about send-side accessibility because it's not accessed through kokkos
  // rely on MPI to do the right thing
  constexpr bool recvDevAccess = Kokkos::SpaceAccessibility<
      Kokkos::DefaultExecutionSpace, typename ImpView::memory_space>::accessible;
  constexpr bool sendDevAccess = Kokkos::SpaceAccessibility<
      Kokkos::DefaultExecutionSpace, typename ExpView::memory_space>::accessible;
  static_assert(recvDevAccess == sendDevAccess, "sending across host/device");

  const int err = ialltofewv_.impl.post<recvDevAccess>(exports.data(), ialltofewv_.sendcounts->data(), ialltofewv_.sdispls->data(), rawType,
                                                       imports.data(), ialltofewv_.recvcounts->data(), ialltofewv_.rdispls->data(),
                                                       roots, nroots,
                                                       rawType,
                                                       mpiTag_, mpiComm, &(*ialltofewv_.req));

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "ialltofewv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");
}

template <class ExpView, class ImpView>
void DistributorActor::doPostsAllToAllImpl(const DistributorPlan& plan,
                                           const ExpView& exports,
                                           const SubViewLimits& exportSubViewLimits,
                                           const ImpView& imports,
                                           const SubViewLimits& importSubViewLimits) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");

  using size_type = Teuchos::Array<size_t>::size_type;

  auto comm = plan.getComm();
  std::vector<int> sendcounts(comm->getSize(), 0);
  std::vector<int> sdispls(comm->getSize(), 0);
  std::vector<int> recvcounts(comm->getSize(), 0);
  std::vector<int> rdispls(comm->getSize(), 0);

  auto& [importStarts, importLengths] = importSubViewLimits;
  auto& [exportStarts, exportLengths] = exportSubViewLimits;

  for (size_t pp = 0; pp < plan.getNumSends(); ++pp) {
    sdispls[plan.getProcsTo()[pp]] = exportStarts[pp];
    size_t numPackets              = exportLengths[pp];
    // numPackets is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(numPackets > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPostsAllToAll: "
                               "Send count for send "
                                   << pp << " (" << numPackets
                                   << ") is too large "
                                      "to be represented as int.");
    sendcounts[plan.getProcsTo()[pp]] = static_cast<int>(numPackets);
  }

  const size_type actualNumReceives =
      Teuchos::as<size_type>(plan.getNumReceives()) +
      Teuchos::as<size_type>(plan.hasSelfMessage() ? 1 : 0);

  for (size_type i = 0; i < actualNumReceives; ++i) {
    rdispls[plan.getProcsFrom()[i]] = importStarts[i];
    size_t totalPacketsFrom_i       = importLengths[i];
    // totalPacketsFrom_i is converted down to int, so make sure it can be
    // represented
    TEUCHOS_TEST_FOR_EXCEPTION(totalPacketsFrom_i > size_t(INT_MAX),
                               std::logic_error,
                               "Tpetra::Distributor::doPostsAllToAll: "
                               "Recv count for receive "
                                   << i << " (" << totalPacketsFrom_i
                                   << ") is too large "
                                      "to be represented as int.");
    recvcounts[plan.getProcsFrom()[i]] = static_cast<int>(totalPacketsFrom_i);
  }

  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiComm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm);
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> rawComm =
      mpiComm->getRawMpiComm();
  using T              = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  if (Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL == plan.getSendType()) {
    MPIX_Comm* mpixComm = *plan.getMPIXComm();
    TEUCHOS_TEST_FOR_EXCEPTION(!mpixComm, std::runtime_error,
                               "MPIX_Comm is null in doPostsAllToAll \""
                                   << __FILE__ << ":" << __LINE__);

    const int err = MPIX_Alltoallv(
        exports.data(), sendcounts.data(), sdispls.data(), rawType,
        imports.data(), recvcounts.data(), rdispls.data(), rawType, mpixComm);

    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "MPIX_Alltoallv failed with error \""
                                   << Teuchos::mpiErrorCodeToString(err)
                                   << "\".");

    return;
  }
#endif  // HAVE_TPETRACORE_MPI_ADVANCE

  const int err = MPI_Alltoallv(
      exports.data(), sendcounts.data(), sdispls.data(), rawType,
      imports.data(), recvcounts.data(), rdispls.data(), rawType, (*rawComm)());

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "MPI_Alltoallv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");
}

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
template <class ExpView, class ImpView>
void DistributorActor::doPostsNbrAllToAllVImpl(const DistributorPlan& plan,
                                               const ExpView& exports,
                                               const SubViewLimits& exportSubViewLimits,
                                               const ImpView& imports,
                                               const SubViewLimits& importSubViewLimits) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");

  const int myRank    = plan.getComm()->getRank();
  MPIX_Comm* mpixComm = *plan.getMPIXComm();
  using size_type     = Teuchos::Array<size_t>::size_type;

  const size_t numSends = plan.getNumSends() + plan.hasSelfMessage();
  const size_t numRecvs = plan.getNumReceives() + plan.hasSelfMessage();
  std::vector<int> sendcounts(numSends, 0);
  std::vector<int> sdispls(numSends, 0);
  std::vector<int> recvcounts(numRecvs, 0);
  std::vector<int> rdispls(numRecvs, 0);

  auto& [importStarts, importLengths] = importSubViewLimits;
  auto& [exportStarts, exportLengths] = exportSubViewLimits;

  for (size_t pp = 0; pp < numSends; ++pp) {
    sdispls[pp]       = exportStarts[pp];
    size_t numPackets = exportLengths[pp];
    // numPackets is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(numPackets > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPostsNbrAllToAllV: "
                               "Send count for send "
                                   << pp << " (" << numPackets
                                   << ") is too large "
                                      "to be represented as int.");
    sendcounts[pp] = static_cast<int>(numPackets);
  }

  for (size_type i = 0; i < numRecvs; ++i) {
    rdispls[i]                = importStarts[i];
    size_t totalPacketsFrom_i = importLengths[i];
    // totalPacketsFrom_i is converted down to int, so make sure it can be
    // represented
    TEUCHOS_TEST_FOR_EXCEPTION(totalPacketsFrom_i > size_t(INT_MAX),
                               std::logic_error,
                               "Tpetra::Distributor::doPostsNbrAllToAllV: "
                               "Recv count for receive "
                                   << i << " (" << totalPacketsFrom_i
                                   << ") is too large "
                                      "to be represented as int.");
    recvcounts[i] = static_cast<int>(totalPacketsFrom_i);
  }

  using T              = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());

  MPIX_Info* xinfo;
  MPIX_Topo* xtopo;
  MPIX_Info_init(&xinfo);
  MPIX_Topo_init(numRecvs, plan.getProcsFrom().data(), recvcounts.data(),
                 numSends, plan.getProcsTo().data(), sendcounts.data(), xinfo, &xtopo);
  const int err = MPIX_Neighbor_alltoallv_topo(
      exports.data(), sendcounts.data(), sdispls.data(), rawType,
      imports.data(), recvcounts.data(), rdispls.data(), rawType, xtopo, mpixComm);
  MPIX_Topo_free(&xtopo);
  MPIX_Info_free(&xinfo);

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "MPIX_Neighbor_alltoallv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");
}
#endif  // HAVE_TPETRACORE_MPI_ADVANCE
#endif  // HAVE_TPETRA_MPI

template <class ImpView>
void DistributorActor::doPostRecvs(const DistributorPlan& plan,
                                   size_t numPackets,
                                   const ImpView& imports) {
  auto importSubViewLimits = plan.getImportViewLimits(numPackets);
  doPostRecvsImpl(plan, imports, importSubViewLimits);
}

template <class ImpView>
void DistributorActor::doPostRecvs(const DistributorPlan& plan,
                                   const ImpView& imports,
                                   const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID) {
  auto importSubViewLimits = plan.getImportViewLimits(numImportPacketsPerLID);
  doPostRecvsImpl(plan, imports, importSubViewLimits);
}

template <class ImpView>
void DistributorActor::doPostRecvsImpl(const DistributorPlan& plan,
                                       const ImpView& imports,
                                       const SubViewLimits& importSubViewLimits) {
  static_assert(isKokkosView<ImpView>,
                "Data arrays for DistributorActor::doPostRecvs must be Kokkos::Views");
  using Kokkos::Compat::subview_offset;
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::ireceive;
  using size_type         = Array<size_t>::size_type;
  using imports_view_type = ImpView;

  // Set the MPI message tag for this round of communication.
  // The same tag will be used for recvs and sends. For every round of communication,
  // the tag gets incremented by one, until it eventually gets looped around back to a
  // small value. This logic is implemented in Teuchos.
  auto comm = plan.getComm();
  {
    auto non_const_comm = Teuchos::rcp_const_cast<Teuchos::Comm<int>>(comm);
    mpiTag_             = non_const_comm->incrementTag();
  }

#ifdef KOKKOS_ENABLE_CUDA
  static_assert(!std::is_same<typename ImpView::memory_space, Kokkos::CudaUVMSpace>::value,
                "Please do not use Tpetra::Distributor with UVM "
                "allocations.  See GitHub issue #1088.");
#endif  // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_SYCL
  static_assert(!std::is_same<typename ImpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value,
                "Please do not use Tpetra::Distributor with SharedUSM "
                "allocations.  See GitHub issue #1088 (corresponding to CUDA).");
#endif  // KOKKOS_ENABLE_SYCL

#if defined(HAVE_TPETRA_MPI)
  // All-to-all communication layout is quite different from
  // point-to-point, so we handle it separately.

  // These send options require no matching receives, so we just return.
  const Details::EDistributorSendType sendType = plan.getSendType();
  if ((sendType == Details::DISTRIBUTOR_ALLTOALL) || (sendType == Details::DISTRIBUTOR_IALLTOFEWV)
#ifdef HAVE_TPETRACORE_MPI_ADVANCE
      || (sendType == Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL) || (sendType == Details::DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV)
#endif
  ) {
    return;
  }
#endif  // HAVE_TPETRA_MPI

  ProfilingRegion pr("Tpetra::Distributor::doPostRecvs");

  const int myProcID = plan.getComm()->getRank();

  auto& [importStarts, importLengths] = importSubViewLimits;

  // Distributor uses requestsRecv_.size() and requestsSend_.size()
  // as the number of outstanding nonblocking message requests, so
  // we resize to zero to maintain this invariant.
  //
  // getNumReceives() does _not_ include the self message, if there is
  // one.  Here, we do actually send a message to ourselves, so we
  // include any self message in the "actual" number of receives to
  // post.
  //
  // NOTE (mfh 19 Mar 2012): Epetra_MpiDistributor::DoPosts()
  // doesn't (re)allocate its array of requests.  That happens in
  // CreateFromSends(), ComputeRecvs_(), DoReversePosts() (on
  // demand), or Resize_().
  const size_type actualNumReceives = as<size_type>(plan.getNumReceives()) +
                                      as<size_type>(plan.hasSelfMessage() ? 1 : 0);

#ifdef HAVE_TPETRA_DEBUG
  size_t totalNumImportPackets = 0;
  for (size_t i = 0; i < Teuchos::as<size_t>(actualNumReceives); ++i) {
    totalNumImportPackets += importLengths[i];
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
      imports.extent(0) < totalNumImportPackets, std::runtime_error,
      "Tpetra::Distributor::doPostRecvs: The 'imports' array must have "
      "enough entries to hold the expected number of import packets.  "
      "imports.extent(0) = "
          << imports.extent(0) << " < "
                                  "totalNumImportPackets = "
          << totalNumImportPackets << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(!requestsRecv_.empty(), std::logic_error,
                             "Tpetra::Distributor::"
                             "doPostRecvs: Process "
                                 << myProcID << ": requestsRecv_.size () = "
                                 << requestsRecv_.size() << " != 0.");
#endif  // HAVE_TPETRA_DEBUG

  requestsRecv_.resize(0);

  // Post the nonblocking receives.  It's common MPI wisdom to post
  // receives before sends.  In MPI terms, this means favoring
  // adding to the "posted queue" (of receive requests) over adding
  // to the "unexpected queue" (of arrived messages not yet matched
  // with a receive).
  {
    ProfilingRegion prr("Tpetra::Distributor::doPostRecvs MPI_Irecv");

    for (size_type i = 0; i < actualNumReceives; ++i) {
      size_t totalPacketsFrom_i = importLengths[Teuchos::as<size_t>(i)];
      TEUCHOS_TEST_FOR_EXCEPTION(totalPacketsFrom_i > size_t(INT_MAX),
                                 std::logic_error,
                                 "Tpetra::Distributor::doPostRecvs: "
                                 "Recv count for receive "
                                     << i << " (" << totalPacketsFrom_i << ") is too large "
                                                                           "to be represented as int.");
      if (plan.getProcsFrom()[i] != myProcID && totalPacketsFrom_i) {
        // If my process is receiving these packet(s) from another
        // process (not a self-receive), and if there is at least
        // one packet to receive:
        //
        // 1. Set up the persisting view (recvBuf) into the imports
        //    array, given the offset and size (total number of
        //    packets from process getProcsFrom()[i]).
        // 2. Start the Irecv and save the resulting request.
        imports_view_type recvBuf =
            subview_offset(imports, importStarts[i], totalPacketsFrom_i);
        requestsRecv_.push_back(ireceive<int>(recvBuf, plan.getProcsFrom()[i],
                                              mpiTag_, *plan.getComm()));
      }
    }
  }
}

template <class ExpView, class ImpView>
void DistributorActor::doPostSends(const DistributorPlan& plan,
                                   const ExpView& exports,
                                   size_t numPackets,
                                   const ImpView& imports) {
  auto exportSubViewLimits = plan.getExportViewLimits(numPackets);
  auto importSubViewLimits = plan.getImportViewLimits(numPackets);
  doPostSendsImpl(plan, exports, exportSubViewLimits, imports, importSubViewLimits);
}

template <class ExpView, class ImpView>
void DistributorActor::doPostSends(const DistributorPlan& plan,
                                   const ExpView& exports,
                                   const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                                   const ImpView& imports,
                                   const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID) {
  auto exportSubViewLimits = plan.getExportViewLimits(numExportPacketsPerLID);
  auto importSubViewLimits = plan.getImportViewLimits(numImportPacketsPerLID);
  doPostSendsImpl(plan, exports, exportSubViewLimits, imports, importSubViewLimits);
}

template <class ExpView, class ImpView>
void DistributorActor::doPostSendsImpl(const DistributorPlan& plan,
                                       const ExpView& exports,
                                       const SubViewLimits& exportSubViewLimits,
                                       const ImpView& imports,
                                       const SubViewLimits& importSubViewLimits) {
  static_assert(areKokkosViews<ExpView, ImpView>,
                "Data arrays for DistributorActor::doPostSends must be Kokkos::Views");
  using Kokkos::Compat::deep_copy_offset;
  using Kokkos::Compat::subview_offset;
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::isend;
  using Teuchos::send;
  using size_type         = Array<size_t>::size_type;
  using exports_view_type = ExpView;

#ifdef KOKKOS_ENABLE_CUDA
  static_assert(!std::is_same<typename ExpView::memory_space, Kokkos::CudaUVMSpace>::value &&
                    !std::is_same<typename ImpView::memory_space, Kokkos::CudaUVMSpace>::value,
                "Please do not use Tpetra::Distributor with UVM allocations.  "
                "See Trilinos GitHub issue #1088.");
#endif  // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_SYCL
  static_assert(!std::is_same<typename ExpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value &&
                    !std::is_same<typename ImpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value,
                "Please do not use Tpetra::Distributor with SharedUSM allocations.  "
                "See Trilinos GitHub issue #1088 (corresponding to CUDA).");
#endif  // KOKKOS_ENABLE_SYCL

  ProfilingRegion ps("Tpetra::Distributor::doPostSends");

  const int myRank = plan.getComm()->getRank();
  // Run-time configurable parameters that come from the input
  // ParameterList set by setParameterList().
  const Details::EDistributorSendType sendType = plan.getSendType();

  auto& [exportStarts, exportLengths] = exportSubViewLimits;
  auto& [importStarts, importLengths] = importSubViewLimits;

#if defined(HAVE_TPETRA_MPI)
  //  All-to-all communication layout is quite different from
  //  point-to-point, so we handle it separately.

  if (sendType == Details::DISTRIBUTOR_ALLTOALL) {
    doPostsAllToAllImpl(plan, exports, exportSubViewLimits, imports, importSubViewLimits);
    return;
  } else if (sendType == Details::DISTRIBUTOR_IALLTOFEWV) {
    doPostsIalltofewvImpl(plan, exports, exportSubViewLimits, imports, importSubViewLimits);
    return;
  }
#ifdef HAVE_TPETRACORE_MPI_ADVANCE
  else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL) {
    doPostsAllToAllImpl(plan, exports, exportSubViewLimits, imports, importSubViewLimits);
    return;
  } else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV) {
    doPostsNbrAllToAllVImpl(plan, exports, exportSubViewLimits, imports, importSubViewLimits);
    return;
  }
#endif  // defined(HAVE_TPETRACORE_MPI_ADVANCE)

#else   // HAVE_TPETRA_MPI
  if (plan.hasSelfMessage()) {
    // This is how we "send a message to ourself": we copy from
    // the export buffer to the import buffer.  That saves
    // Teuchos::Comm implementations other than MpiComm (in
    // particular, SerialComm) the trouble of implementing self
    // messages correctly.  (To do this right, SerialComm would
    // need internal buffer space for messages, keyed on the
    // message's tag.)
    size_t selfReceiveOffset = 0;
    deep_copy_offset(imports, exports, selfReceiveOffset,
                     exportStarts[0],
                     exportLengths[0]);
  }
  // should we just return here?
  // likely not as comm could be a serial comm
#endif  // HAVE_TPETRA_MPI

  size_t selfReceiveOffset = 0;

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(requestsSend_.size() != 0,
                             std::logic_error,
                             "Tpetra::Distributor::doPostSends: Process "
                                 << myRank << ": requestsSend_.size() = " << requestsSend_.size() << " != 0.");
#endif  // HAVE_TPETRA_DEBUG

  // Distributor uses requestsRecv_.size() and requestsSend_.size()
  // as the number of outstanding nonblocking message requests, so
  // we resize to zero to maintain this invariant.
  //
  // getNumReceives() does _not_ include the self message, if there is
  // one.  Here, we do actually send a message to ourselves, so we
  // include any self message in the "actual" number of receives to
  // post.
  //
  // NOTE (mfh 19 Mar 2012): Epetra_MpiDistributor::DoPosts()
  // doesn't (re)allocate its array of requests.  That happens in
  // CreateFromSends(), ComputeRecvs_(), DoReversePosts() (on
  // demand), or Resize_().
  const size_type actualNumReceives = as<size_type>(plan.getNumReceives()) +
                                      as<size_type>(plan.hasSelfMessage() ? 1 : 0);
  requestsSend_.resize(0);

  {
    for (size_type i = 0; i < actualNumReceives; ++i) {
      if (plan.getProcsFrom()[i] == myRank) {  // Receiving from myself
        selfReceiveOffset = importStarts[i];   // Remember the self-recv offset
      }
    }
  }

  ProfilingRegion pss("Tpetra::Distributor::doPostSends sends");

  // setup scan through getProcsTo() list starting with higher numbered procs
  // (should help balance message traffic)
  //
  // FIXME (mfh 20 Feb 2013) Why haven't we precomputed this?
  // It doesn't depend on the input at all.
  size_t numBlocks = plan.getNumSends() + plan.hasSelfMessage();
  size_t procIndex = 0;
  while ((procIndex < numBlocks) && (plan.getProcsTo()[procIndex] < myRank)) {
    ++procIndex;
  }
  if (procIndex == numBlocks) {
    procIndex = 0;
  }

  size_t selfNum   = 0;
  size_t selfIndex = 0;

  if (plan.getIndicesTo().is_null()) {
    const char isend_region[] = "Tpetra::Distributor::doPostSends MPI_Isend FAST";
    const char send_region[]  = "Tpetra::Distributor::doPostSends MPI_Send FAST";
    ProfilingRegion pssf((sendType == Details::DISTRIBUTOR_ISEND) ? isend_region : send_region);

    // Data are already blocked (laid out) by process, so we don't
    // need a separate send buffer (besides the exports array).
    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myRank) {
        if (exportLengths[p] == 0) {
          // Do not attempt to send messages of length 0.
          continue;
        }

        exports_view_type tmpSend = subview_offset(exports, exportStarts[p], exportLengths[p]);

        if (sendType == Details::DISTRIBUTOR_ISEND) {
          // NOTE: This looks very similar to the tmpSend above, but removing
          // tmpSendBuf and uses tmpSend leads to a performance hit on Arm
          // SerialNode builds
          exports_view_type tmpSendBuf =
              subview_offset(exports, exportStarts[p], exportLengths[p]);
          requestsSend_.push_back(isend<int>(tmpSendBuf, plan.getProcsTo()[p],
                                             mpiTag_, *plan.getComm()));
        } else {  // DISTRIBUTOR_SEND
          send<int>(tmpSend,
                    as<int>(tmpSend.size()),
                    plan.getProcsTo()[p], mpiTag_, *plan.getComm());
        }
      } else {  // "Sending" the message to myself
        selfNum = p;
      }
    }

    if (plan.hasSelfMessage()) {
      // This is how we "send a message to ourself": we copy from
      // the export buffer to the import buffer.  That saves
      // Teuchos::Comm implementations other than MpiComm (in
      // particular, SerialComm) the trouble of implementing self
      // messages correctly.  (To do this right, SerialComm would
      // need internal buffer space for messages, keyed on the
      // message's tag.)
      deep_copy_offset(imports, exports, selfReceiveOffset,
                       exportStarts[selfNum], exportLengths[selfNum]);
    }

  } else {  // data are not blocked by proc, use send buffer
    ProfilingRegion psss("Tpetra::Distributor::doPostSends: MPI_Send SLOW");

    using Packet = typename ExpView::non_const_value_type;
    using Layout = typename ExpView::array_layout;
    using Device = typename ExpView::device_type;
    using Mem    = typename ExpView::memory_traits;

    // This buffer is long enough for only one message at a time.
    // Thus, we use DISTRIBUTOR_SEND always in this case, regardless
    // of sendType requested by user.
    // This code path formerly errored out with message:
    //     Tpetra::Distributor::doPosts(3 args):
    //     The "send buffer" code path
    //     doesn't currently work with nonblocking sends.
    // Now, we opt to just do the communication in a way that works.
#ifdef HAVE_TPETRA_DEBUG
    if (sendType != Details::DISTRIBUTOR_SEND) {
      if (plan.getComm()->getRank() == 0)
        std::cout << "The requested Tpetra send type "
                  << DistributorSendTypeEnumToString(sendType)
                  << " requires Distributor data to be ordered by"
                  << " the receiving processor rank.  Since these"
                  << " data are not ordered, Tpetra will use Send"
                  << " instead." << std::endl;
    }
#endif

    size_t maxSendLength = 0;
    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      size_t sendArrayOffset = 0;
      size_t j               = plan.getStartsTo()[p];
      for (size_t k = 0; k < plan.getLengthsTo()[p]; ++k, ++j) {
        sendArrayOffset += exportLengths[j];
      }
      maxSendLength = std::max(maxSendLength, sendArrayOffset);
    }
    Kokkos::View<Packet*, Layout, Device, Mem> sendArray("sendArray", maxSendLength);

    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myRank) {
        size_t sendArrayOffset = 0;
        size_t j               = plan.getStartsTo()[p];
        for (size_t k = 0; k < plan.getLengthsTo()[p]; ++k, ++j) {
          packOffset(sendArray, exports, sendArrayOffset, exportStarts[j], exportLengths[j]);
          sendArrayOffset += exportLengths[j];
        }
        typename ExpView::execution_space().fence();

        ImpView tmpSend =
            subview_offset(sendArray, size_t(0), sendArrayOffset);

        send<int>(tmpSend,
                  as<int>(tmpSend.size()),
                  plan.getProcsTo()[p], mpiTag_, *plan.getComm());
      } else {  // "Sending" the message to myself
        selfNum   = p;
        selfIndex = plan.getStartsTo()[p];
      }
    }

    if (plan.hasSelfMessage()) {
      for (size_t k = 0; k < plan.getLengthsTo()[selfNum]; ++k) {
        packOffset(imports, exports, selfReceiveOffset, exportStarts[selfIndex], exportLengths[selfIndex]);
        selfReceiveOffset += exportLengths[selfIndex];
        ++selfIndex;
      }
    }
  }
}

}  // namespace Tpetra::Details

#endif
