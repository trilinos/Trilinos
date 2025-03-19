// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format off

#ifndef TPETRA_DETAILS_DISTRIBUTOR_ACTOR_HPP
#define TPETRA_DETAILS_DISTRIBUTOR_ACTOR_HPP
#define HAVE_TPETRA_DISTRIBUTOR_TIMINGS

#include "Tpetra_Details_DistributorPlan.hpp"
#include "Tpetra_Util.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include "Tpetra_Details_MpiTypeTraits.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Time.hpp"

#include "Kokkos_TeuchosCommAdapters.hpp"
#include "Kokkos_StdAlgorithms.hpp"

#ifdef HAVE_TPETRA_MPI
#include "mpi.h"
#endif

namespace Tpetra {
namespace Details {

template <class View1, class View2>
constexpr bool areKokkosViews = Kokkos::is_view<View1>::value && Kokkos::is_view<View2>::value;

class DistributorActor {
  static constexpr int DEFAULT_MPI_TAG = 1;

public:
  DistributorActor();
  DistributorActor(const DistributorActor& otherActor);

  template <class ExpView, class ImpView>
  void doPostsAndWaits(const DistributorPlan& plan,
                       const ExpView &exports,
                       size_t numPackets,
                       const ImpView &imports);
  
  template <class ExpView, class ImpView>
  void doPostsAndWaitsKokkos(const DistributorPlan& plan,
                             const ExpView &exports,
                             size_t numPackets,
                             const ImpView &imports);

  template <class ExpView, class ImpView>
  void doPostsAndWaits(const DistributorPlan& plan,
                       const ExpView &exports,
                       const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                       const ImpView &imports,
                       const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

  template <class ExpView, class ExpPacketsView, class ImpView, class ImpPacketsView>
  void doPostsAndWaitsKokkos(const DistributorPlan& plan,
                             const ExpView &exports,
                             const ExpPacketsView &numExportPacketsPerLID,
                             const ImpView &imports,
                             const ImpPacketsView &numImportPacketsPerLID);

  template <class ExpView, class ImpView>
  void doPosts(const DistributorPlan& plan,
               const ExpView& exports,
               size_t numPackets,
               const ImpView& imports);
  
  template <class ExpView, class ImpView>
  void doPostsKokkos(const DistributorPlan& plan,
                     const ExpView& exports,
                     size_t numPackets,
                     const ImpView& imports);

  template <class ExpView, class ImpView>
  void doPosts(const DistributorPlan& plan,
               const ExpView &exports,
               const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
               const ImpView &imports,
               const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID);

  template <class ExpView, class ExpPacketsView, class ImpView, class ImpPacketsView>
  void doPostsKokkos(const DistributorPlan& plan,
                     const ExpView &exports,
                     const ExpPacketsView &numExportPacketsPerLID,
                     const ImpView &imports,
                     const ImpPacketsView &numImportPacketsPerLID);
  
  template <class ExpView, class ImpView>
  void doPostsAllToAllKokkos(const DistributorPlan &plan, const ExpView &exports,
                             size_t numPackets, const ImpView &imports);

  template <class ExpView, class ExpPacketsView, class ImpView, class ImpPacketsView>
  void doPostsAllToAllKokkos(
      const DistributorPlan &plan, const ExpView &exports,
      const ExpPacketsView &numExportPacketsPerLID,
      const ImpView &imports,
      const ImpPacketsView &numImportPacketsPerLID);
  
  template <class ExpView, class ImpView>
  void doPostsNbrAllToAllVKokkos(const DistributorPlan &plan, const ExpView &exports,
                                 size_t numPackets, const ImpView &imports);
  
  template <class ExpView, class ExpPacketsView, class ImpView, class ImpPacketsView>
  void doPostsNbrAllToAllVKokkos(
      const DistributorPlan &plan, const ExpView &exports,
      const ExpPacketsView &numExportPacketsPerLID,
      const ImpView &imports,
      const ImpPacketsView &numImportPacketsPerLID);

  void doWaits(const DistributorPlan& plan);

  bool isReady() const;

private:
// clang-format on
#ifdef HAVE_TPETRA_MPI
  template <class ExpView, class ImpView>
  void doPostsAllToAll(const DistributorPlan &plan, const ExpView &exports,
                       size_t numPackets, const ImpView &imports);
  
  template <class ExpView, class ImpView>
  void doPostsAllToAll(
      const DistributorPlan &plan, const ExpView &exports,
      const Teuchos::ArrayView<const size_t> &numExportPacketsPerLID,
      const ImpView &imports,
      const Teuchos::ArrayView<const size_t> &numImportPacketsPerLID);
  
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  template <class ExpView, class ImpView>
  void doPostsNbrAllToAllV(const DistributorPlan &plan, const ExpView &exports,
                           size_t numPackets, const ImpView &imports);
  
  template <class ExpView, class ImpView>
  void doPostsNbrAllToAllV(
      const DistributorPlan &plan, const ExpView &exports,
      const Teuchos::ArrayView<const size_t> &numExportPacketsPerLID,
      const ImpView &imports,
      const Teuchos::ArrayView<const size_t> &numImportPacketsPerLID);
  
#endif // HAVE_TPETRACORE_MPI_ADVANCE
#endif // HAVE_TPETRA_CORE
  // clang-format off
  int mpiTag_;

  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest<int>>> requests_;

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::RCP<Teuchos::Time> timer_doPosts3KV_;
  Teuchos::RCP<Teuchos::Time> timer_doPosts4KV_;
  Teuchos::RCP<Teuchos::Time> timer_doWaits_;
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

template <class ExpView, class ImpView>
void DistributorActor::doPostsAndWaits(const DistributorPlan& plan,
                                       const ExpView& exports,
                                       size_t numPackets,
                                       const ImpView& imports)
{
  static_assert(areKokkosViews<ExpView, ImpView>,
      "Data arrays for DistributorActor::doPostsAndWaits must be Kokkos::Views");
  doPosts(plan, exports, numPackets, imports);
  doWaits(plan);
}

template <class ExpView, class ImpView>
void DistributorActor::doPostsAndWaitsKokkos(const DistributorPlan& plan,
                                             const ExpView& exports,
                                             size_t numPackets,
                                             const ImpView& imports)
{
  static_assert(areKokkosViews<ExpView, ImpView>,
      "Data arrays for DistributorActor::doPostsAndWaitsKokkos must be Kokkos::Views");
  doPostsKokkos(plan, exports, numPackets, imports);
  doWaits(plan);
}

template <class ExpView, class ImpView>
void DistributorActor::doPostsAndWaits(const DistributorPlan& plan,
                                       const ExpView& exports,
                                       const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                                       const ImpView& imports,
                                       const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
{
  static_assert(areKokkosViews<ExpView, ImpView>,
      "Data arrays for DistributorActor::doPostsAndWaits must be Kokkos::Views");
  doPosts(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
  doWaits(plan);
}


template <class ExpView, class ExpPacketsView, class ImpView, class ImpPacketsView>
void DistributorActor::doPostsAndWaitsKokkos(const DistributorPlan& plan,
                                             const ExpView &exports,
                                             const ExpPacketsView &numExportPacketsPerLID,
                                             const ImpView &imports,
                                             const ImpPacketsView &numImportPacketsPerLID)
{
  static_assert(areKokkosViews<ExpView, ImpView>,
      "Data arrays for DistributorActor::doPostsAndWaitsKokkos must be Kokkos::Views");
  static_assert(areKokkosViews<ExpPacketsView, ImpPacketsView>,
      "Num packets arrays for DistributorActor::doPostsAndWaitsKokkos must be Kokkos::Views");
  doPostsKokkos(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
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
           const size_t size)
{
  memcpy((void*) (dst.data()+dst_offset), src.data()+src_offset, size*sizeof(typename DstViewType::value_type));
}

template <typename DstViewType, typename SrcViewType>
enableIfNotHostAccessible<DstViewType, SrcViewType>
packOffset(const DstViewType& dst,
           const SrcViewType& src,
           const size_t dst_offset,
           const size_t src_offset,
           const size_t size)
{
  Kokkos::Compat::deep_copy_offset(dst, src, dst_offset, src_offset, size);
}

// clang-format on
#ifdef HAVE_TPETRA_MPI
template <class ExpView, class ImpView>
void DistributorActor::doPostsAllToAll(const DistributorPlan &plan,
                                       const ExpView &exports,
                                       size_t numPackets,
                                       const ImpView &imports) {
  using size_type = Teuchos::Array<size_t>::size_type;

  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");

  auto comm = plan.getComm();
  const int myRank = comm->getRank();
  std::vector<int> sendcounts(comm->getSize(), 0);
  std::vector<int> sdispls(comm->getSize(), 0);
  std::vector<int> recvcounts(comm->getSize(), 0);
  std::vector<int> rdispls(comm->getSize(), 0);

  size_t numBlocks = plan.getNumSends() + plan.hasSelfMessage();
  for (size_t p = 0; p < numBlocks; ++p) {
    sdispls[plan.getProcsTo()[p]] = plan.getStartsTo()[p] * numPackets;
    size_t sendcount = plan.getLengthsTo()[p] * numPackets;
    // sendcount is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(sendcount > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Send count for block "
                                   << p << " (" << sendcount
                                   << ") is too large "
                                      "to be represented as int.");
    sendcounts[plan.getProcsTo()[p]] = static_cast<int>(sendcount);
  }

  const size_type actualNumReceives =
      Teuchos::as<size_type>(plan.getNumReceives()) +
      Teuchos::as<size_type>(plan.hasSelfMessage() ? 1 : 0);
  size_t curBufferOffset = 0;
  for (size_type i = 0; i < actualNumReceives; ++i) {
    const size_t curBufLen = plan.getLengthsFrom()[i] * numPackets;
    TEUCHOS_TEST_FOR_EXCEPTION(
        curBufferOffset + curBufLen > static_cast<size_t>(imports.size()),
        std::logic_error,
        "Tpetra::Distributor::doPosts(3 args, Kokkos): "
        "Exceeded size of 'imports' array in packing loop on Process "
            << myRank << ".  imports.size() = " << imports.size()
            << " < "
               "curBufferOffset("
            << curBufferOffset << ") + curBufLen(" << curBufLen << ").");
    rdispls[plan.getProcsFrom()[i]] = curBufferOffset;
    // curBufLen is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(curBufLen > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Recv count for receive "
                                   << i << " (" << curBufLen
                                   << ") is too large "
                                      "to be represented as int.");
    recvcounts[plan.getProcsFrom()[i]] = static_cast<int>(curBufLen);
    curBufferOffset += curBufLen;
  }

  using T = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  if (Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL == plan.getSendType()) {
    MPIX_Comm *mpixComm = *plan.getMPIXComm();
    TEUCHOS_TEST_FOR_EXCEPTION(
        !mpixComm, std::runtime_error,
        "plan's MPIX_Comm null in doPostsAllToAll, but "
        "DISTRIBUTOR_MPIADVANCE_ALLTOALL set: plan.howInitialized()="
            << DistributorHowInitializedEnumToString(plan.howInitialized()));

    const int err = MPIX_Alltoallv(
        exports.data(), sendcounts.data(), sdispls.data(), rawType,
        imports.data(), recvcounts.data(), rdispls.data(), rawType, mpixComm);

    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "MPIX_Alltoallv failed with error \""
                                   << Teuchos::mpiErrorCodeToString(err)
                                   << "\".");

    return;
  }
#endif
  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiComm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm);
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> rawComm =
      mpiComm->getRawMpiComm();

  const int err = MPI_Alltoallv(
      exports.data(), sendcounts.data(), sdispls.data(), rawType,
      imports.data(), recvcounts.data(), rdispls.data(), rawType, (*rawComm)());

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "MPI_Alltoallv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");

  return;
}

template <class ExpView, class ImpView>
void DistributorActor::doPostsAllToAllKokkos(const DistributorPlan &plan,
                                             const ExpView &exports,
                                             size_t numPackets,
                                             const ImpView &imports) {
  using size_type = Teuchos::Array<size_t>::size_type;

  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");

  auto comm = plan.getComm();
  const int myRank = comm->getRank();
  Kokkos::View<int*, Kokkos::DefaultExecutionSpace> sendcounts_d("sendcounts_d", comm->getSize());
  Kokkos::View<int*, Kokkos::DefaultExecutionSpace> sdispls_d("sdispls_d", comm->getSize());
  Kokkos::View<int*, Kokkos::DefaultExecutionSpace> recvcounts_d("recvcounts_d", comm->getSize());
  Kokkos::View<int*, Kokkos::DefaultExecutionSpace> rdispls_d("rdispls_d", comm->getSize());

  auto getStartsTo = Kokkos::Compat::getKokkosViewDeepCopy<Kokkos::DefaultExecutionSpace> (plan.getStartsTo());
  auto getLengthsTo = Kokkos::Compat::getKokkosViewDeepCopy<Kokkos::DefaultExecutionSpace> (plan.getLengthsTo());
  auto getProcsTo = Kokkos::Compat::getKokkosViewDeepCopy<Kokkos::DefaultExecutionSpace> (plan.getProcsTo());

  size_t numBlocks = plan.getNumSends() + plan.hasSelfMessage();
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, numBlocks), KOKKOS_LAMBDA(const size_t p) {
    sdispls_d(getProcsTo(p)) = getStartsTo(p) * numPackets;
    size_t sendcount = getLengthsTo(p) * numPackets;
    /*
    // sendcount is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(sendcount > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Send count for block "
                                   << p << " (" << sendcount
                                   << ") is too large "
                                      "to be represented as int.");
    */
    sendcounts_d(getProcsTo(p)) = static_cast<int>(sendcount);
  });

  const size_type actualNumReceives =
      Teuchos::as<size_type>(plan.getNumReceives()) +
      Teuchos::as<size_type>(plan.hasSelfMessage() ? 1 : 0);
  auto getLengthsFrom = Kokkos::Compat::getKokkosViewDeepCopy<Kokkos::DefaultExecutionSpace> (plan.getLengthsFrom());
  auto getProcsFrom = Kokkos::Compat::getKokkosViewDeepCopy<Kokkos::DefaultExecutionSpace> (plan.getProcsFrom());
  
  Kokkos::parallel_scan(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, actualNumReceives), KOKKOS_LAMBDA(const size_type i, size_t& offset, bool is_final) {
    const size_t curBufLen = getLengthsFrom(i) * numPackets;
    /*
    TEUCHOS_TEST_FOR_EXCEPTION(
        curBufferOffset + curBufLen > static_cast<size_t>(imports.size()),
        std::logic_error,
        "Tpetra::Distributor::doPosts(3 args, Kokkos): "
        "Exceeded size of 'imports' array in packing loop on Process "
            << myRank << ".  imports.size() = " << imports.size()
            << " < "
               "curBufferOffset("
            << curBufferOffset << ") + curBufLen(" << curBufLen << ").");
    */
    if (is_final) rdispls_d(getProcsFrom(i)) = offset;
    /*
    // curBufLen is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(curBufLen > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Recv count for receive "
                                   << i << " (" << curBufLen
                                   << ") is too large "
                                      "to be represented as int.");
    */
    if (is_final) recvcounts_d(getProcsFrom(i)) = static_cast<int>(curBufLen);
    offset += curBufLen;
  });

  using T = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());
  
  auto sendcounts = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), sendcounts_d);
  auto sdispls = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), sdispls_d);
  auto recvcounts = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), recvcounts_d);
  auto rdispls = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), rdispls_d);

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  if (Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL == plan.getSendType()) {
    MPIX_Comm *mpixComm = *plan.getMPIXComm();
    TEUCHOS_TEST_FOR_EXCEPTION(
        !mpixComm, std::runtime_error,
        "plan's MPIX_Comm null in doPostsAllToAll, but "
        "DISTRIBUTOR_MPIADVANCE_ALLTOALL set: plan.howInitialized()="
            << DistributorHowInitializedEnumToString(plan.howInitialized()));

    const int err = MPIX_Alltoallv(
        exports.data(), sendcounts.data(), sdispls.data(), rawType,
        imports.data(), recvcounts.data(), rdispls.data(), rawType, mpixComm);

    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "MPIX_Alltoallv failed with error \""
                                   << Teuchos::mpiErrorCodeToString(err)
                                   << "\".");

    return;
  }
#endif
  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiComm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm);
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> rawComm =
      mpiComm->getRawMpiComm();

  const int err = MPI_Alltoallv(
      exports.data(), sendcounts.data(), sdispls.data(), rawType,
      imports.data(), recvcounts.data(), rdispls.data(), rawType, (*rawComm)());

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "MPI_Alltoallv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");

  return;
}

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
template <class ExpView, class ImpView>
void DistributorActor::doPostsNbrAllToAllV(const DistributorPlan &plan,
                                           const ExpView &exports,
                                           size_t numPackets,
                                           const ImpView &imports) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");

  const int myRank = plan.getComm()->getRank();
  MPIX_Comm *mpixComm = *plan.getMPIXComm();

  const size_t numSends = plan.getNumSends() + plan.hasSelfMessage();
  const size_t numRecvs = plan.getNumReceives() + plan.hasSelfMessage();
  std::vector<int> sendcounts(numSends, 0);
  std::vector<int> sdispls(numSends, 0);
  std::vector<int> recvcounts(numRecvs, 0);
  std::vector<int> rdispls(numRecvs, 0);

  for (size_t p = 0; p < numSends; ++p) {
    sdispls[p] = plan.getStartsTo()[p] * numPackets;
    const size_t sendcount = plan.getLengthsTo()[p] * numPackets;
    // sendcount is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(sendcount > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Send count for block "
                                   << p << " (" << sendcount
                                   << ") is too large "
                                      "to be represented as int.");
    sendcounts[p] = static_cast<int>(sendcount);
  }

  size_t curBufferOffset = 0;
  for (size_t i = 0; i < numRecvs; ++i) {
    const size_t curBufLen = plan.getLengthsFrom()[i] * numPackets;
    TEUCHOS_TEST_FOR_EXCEPTION(
        curBufferOffset + curBufLen > static_cast<size_t>(imports.size()),
        std::logic_error,
        "Tpetra::Distributor::doPosts(3 args, Kokkos): "
        "Exceeded size of 'imports' array in packing loop on Process "
            << myRank << ".  imports.size() = " << imports.size()
            << " < "
               "curBufferOffset("
            << curBufferOffset << ") + curBufLen(" << curBufLen << ").");
    rdispls[i] = curBufferOffset;
    // curBufLen is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(curBufLen > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Recv count for receive "
                                   << i << " (" << curBufLen
                                   << ") is too large "
                                      "to be represented as int.");
    recvcounts[i] = static_cast<int>(curBufLen);
    curBufferOffset += curBufLen;
  }

  using T = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());

  const int err = MPIX_Neighbor_alltoallv(
      exports.data(), sendcounts.data(), sdispls.data(), rawType,
      imports.data(), recvcounts.data(), rdispls.data(), rawType, mpixComm);

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "MPIX_Neighbor_alltoallv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");
}
template <class ExpView, class ImpView>
void DistributorActor::doPostsNbrAllToAllVKokkos(const DistributorPlan &plan,
                                                 const ExpView &exports,
                                                 size_t numPackets,
                                                 const ImpView &imports) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");

  const int myRank = plan.getComm()->getRank();
  MPIX_Comm *mpixComm = *plan.getMPIXComm();

  const size_t numSends = plan.getNumSends() + plan.hasSelfMessage();
  const size_t numRecvs = plan.getNumReceives() + plan.hasSelfMessage();
  std::vector<int> sendcounts(numSends, 0);
  std::vector<int> sdispls(numSends, 0);
  std::vector<int> recvcounts(numRecvs, 0);
  std::vector<int> rdispls(numRecvs, 0);

  for (size_t p = 0; p < numSends; ++p) {
    sdispls[p] = plan.getStartsTo()[p] * numPackets;
    const size_t sendcount = plan.getLengthsTo()[p] * numPackets;
    // sendcount is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(sendcount > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Send count for block "
                                   << p << " (" << sendcount
                                   << ") is too large "
                                      "to be represented as int.");
    sendcounts[p] = static_cast<int>(sendcount);
  }

  size_t curBufferOffset = 0;
  for (size_t i = 0; i < numRecvs; ++i) {
    const size_t curBufLen = plan.getLengthsFrom()[i] * numPackets;
    TEUCHOS_TEST_FOR_EXCEPTION(
        curBufferOffset + curBufLen > static_cast<size_t>(imports.size()),
        std::logic_error,
        "Tpetra::Distributor::doPosts(3 args, Kokkos): "
        "Exceeded size of 'imports' array in packing loop on Process "
            << myRank << ".  imports.size() = " << imports.size()
            << " < "
               "curBufferOffset("
            << curBufferOffset << ") + curBufLen(" << curBufLen << ").");
    rdispls[i] = curBufferOffset;
    // curBufLen is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(curBufLen > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Recv count for receive "
                                   << i << " (" << curBufLen
                                   << ") is too large "
                                      "to be represented as int.");
    recvcounts[i] = static_cast<int>(curBufLen);
    curBufferOffset += curBufLen;
  }

  using T = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());

  const int err = MPIX_Neighbor_alltoallv(
      exports.data(), sendcounts.data(), sdispls.data(), rawType,
      imports.data(), recvcounts.data(), rdispls.data(), rawType, mpixComm);

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "MPIX_Neighbor_alltoallv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");
}
#endif // HAVE_TPETRACORE_MPI_ADVANCE
#endif // HAVE_TPETRA_MPI
// clang-format off

template <class ExpView, class ImpView>
void DistributorActor::doPosts(const DistributorPlan& plan,
                               const ExpView& exports,
                               size_t numPackets,
                               const ImpView& imports)
{
  static_assert(areKokkosViews<ExpView, ImpView>,
      "Data arrays for DistributorActor::doPosts must be Kokkos::Views");
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::includesVerbLevel;
  using Teuchos::ireceive;
  using Teuchos::isend;
  using Teuchos::send;
  using Teuchos::TypeNameTraits;
  using Teuchos::typeName;
  using std::endl;
  using Kokkos::Compat::create_const_view;
  using Kokkos::Compat::create_view;
  using Kokkos::Compat::subview_offset;
  using Kokkos::Compat::deep_copy_offset;
  typedef Array<size_t>::size_type size_type;
  typedef ExpView exports_view_type;
  typedef ImpView imports_view_type;

#ifdef KOKKOS_ENABLE_CUDA
  static_assert
    (! std::is_same<typename ExpView::memory_space, Kokkos::CudaUVMSpace>::value &&
     ! std::is_same<typename ImpView::memory_space, Kokkos::CudaUVMSpace>::value,
     "Please do not use Tpetra::Distributor with UVM allocations.  "
     "See Trilinos GitHub issue #1088.");
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_SYCL
    static_assert
      (! std::is_same<typename ExpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value &&
       ! std::is_same<typename ImpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value,
       "Please do not use Tpetra::Distributor with SharedUSM allocations.  "
       "See Trilinos GitHub issue #1088 (corresponding to CUDA).");
#endif // KOKKOS_ENABLE_SYCL
  
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::TimeMonitor timeMon (*timer_doPosts3KV_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

  const int myRank = plan.getComm()->getRank ();
  // Run-time configurable parameters that come from the input
  // ParameterList set by setParameterList().
  const Details::EDistributorSendType sendType = plan.getSendType();

//clang-format on
#if defined(HAVE_TPETRA_MPI)
  //  All-to-all communication layout is quite different from
  //  point-to-point, so we handle it separately.

  if (sendType == Details::DISTRIBUTOR_ALLTOALL) {
    doPostsAllToAll(plan, exports,numPackets, imports);
    return;
  }
#ifdef HAVE_TPETRACORE_MPI_ADVANCE
  else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL) {
    doPostsAllToAll(plan, exports,numPackets, imports);
    return;
  } else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV) {
    doPostsNbrAllToAllV(plan, exports,numPackets, imports);
    return;
  }
#endif // defined(HAVE_TPETRACORE_MPI_ADVANCE)
// clang-format off
  
#else // HAVE_TPETRA_MPI
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
                       plan.getStartsTo()[0]*numPackets,
                       plan.getLengthsTo()[0]*numPackets);
    }
    // should we just return here?
    // likely not as comm could be a serial comm
#endif // HAVE_TPETRA_MPI

  size_t selfReceiveOffset = 0;

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION
    (requests_.size () != 0,
     std::logic_error,
     "Tpetra::Distributor::doPosts(3 args, Kokkos): Process "
     << myRank << ": requests_.size() = " << requests_.size () << " != 0.");
#endif // HAVE_TPETRA_DEBUG

  // Distributor uses requests_.size() as the number of outstanding
  // nonblocking message requests, so we resize to zero to maintain
  // this invariant.
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
  const size_type actualNumReceives = as<size_type> (plan.getNumReceives()) +
    as<size_type> (plan.hasSelfMessage() ? 1 : 0);
  requests_.resize (0);

  // Post the nonblocking receives.  It's common MPI wisdom to post
  // receives before sends.  In MPI terms, this means favoring
  // adding to the "posted queue" (of receive requests) over adding
  // to the "unexpected queue" (of arrived messages not yet matched
  // with a receive).
  {
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonRecvs (*timer_doPosts3KV_recvs_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    size_t curBufferOffset = 0;
    for (size_type i = 0; i < actualNumReceives; ++i) {
      const size_t curBufLen = plan.getLengthsFrom()[i] * numPackets;
      if (plan.getProcsFrom()[i] != myRank) {
        // If my process is receiving these packet(s) from another
        // process (not a self-receive):
        //
        // 1. Set up the persisting view (recvBuf) of the imports
        //    array, given the offset and size (total number of
        //    packets from process getProcsFrom()[i]).
        // 2. Start the Irecv and save the resulting request.
        TEUCHOS_TEST_FOR_EXCEPTION(
            curBufferOffset + curBufLen > static_cast<size_t> (imports.size ()),
            std::logic_error, "Tpetra::Distributor::doPosts(3 args, Kokkos): "
            "Exceeded size of 'imports' array in packing loop on Process " <<
            myRank << ".  imports.size() = " << imports.size () << " < "
            "curBufferOffset(" << curBufferOffset << ") + curBufLen(" <<
            curBufLen << ").");
        imports_view_type recvBuf =
          subview_offset (imports, curBufferOffset, curBufLen);
        requests_.push_back (ireceive<int> (recvBuf, plan.getProcsFrom()[i],
              mpiTag_, *plan.getComm()));
      }
      else { // Receiving from myself
        selfReceiveOffset = curBufferOffset; // Remember the self-recv offset
      }
      curBufferOffset += curBufLen;
    }
  }

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::TimeMonitor timeMonSends (*timer_doPosts3KV_sends_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

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

  size_t selfNum = 0;
  size_t selfIndex = 0;

  if (plan.getIndicesTo().is_null()) {

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonSends2 (*timer_doPosts3KV_sends_fast_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    // Data are already blocked (laid out) by process, so we don't
    // need a separate send buffer (besides the exports array).
    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myRank) {
        exports_view_type tmpSend = subview_offset(
            exports, plan.getStartsTo()[p]*numPackets, plan.getLengthsTo()[p]*numPackets);

        if (sendType == Details::DISTRIBUTOR_ISEND) {
          // NOTE: This looks very similar to the tmpSend above, but removing
          // tmpSendBuf and uses tmpSend leads to a performance hit on Arm
          // SerialNode builds
          exports_view_type tmpSendBuf =
            subview_offset (exports, plan.getStartsTo()[p] * numPackets,
                plan.getLengthsTo()[p] * numPackets);
          requests_.push_back (isend<int> (tmpSendBuf, plan.getProcsTo()[p],
                mpiTag_, *plan.getComm()));
        }
        else {  // DISTRIBUTOR_SEND
          send<int> (tmpSend,
              as<int> (tmpSend.size ()),
              plan.getProcsTo()[p], mpiTag_, *plan.getComm());
        }
      }
      else { // "Sending" the message to myself
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
          plan.getStartsTo()[selfNum]*numPackets,
          plan.getLengthsTo()[selfNum]*numPackets);
    }

  }
  else { // data are not blocked by proc, use send buffer

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonSends2 (*timer_doPosts3KV_sends_slow_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    typedef typename ExpView::non_const_value_type Packet;
    typedef typename ExpView::array_layout Layout;
    typedef typename ExpView::device_type Device;
    typedef typename ExpView::memory_traits Mem;

    // This buffer is long enough for only one message at a time.
    // Thus, we use DISTRIBUTOR_SEND always in this case, regardless
    // of sendType requested by user. 
    // This code path formerly errored out with message:
    //     Tpetra::Distributor::doPosts(3 args, Kokkos): 
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

    Kokkos::View<Packet*,Layout,Device,Mem> sendArray ("sendArray",
        plan.getMaxSendLength() * numPackets);

    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myRank) {
        size_t sendArrayOffset = 0;
        size_t j = plan.getStartsTo()[p];
        for (size_t k = 0; k < plan.getLengthsTo()[p]; ++k, ++j) {
          packOffset(sendArray, exports, sendArrayOffset, plan.getIndicesTo()[j]*numPackets, numPackets);
          sendArrayOffset += numPackets;
        }
        typename ExpView::execution_space().fence();

        ImpView tmpSend =
          subview_offset(sendArray, size_t(0), plan.getLengthsTo()[p]*numPackets);

        send<int> (tmpSend,
            as<int> (tmpSend.size ()),
            plan.getProcsTo()[p], mpiTag_, *plan.getComm());
      }
      else { // "Sending" the message to myself
        selfNum = p;
        selfIndex = plan.getStartsTo()[p];
      }
    }

    if (plan.hasSelfMessage()) {
      for (size_t k = 0; k < plan.getLengthsTo()[selfNum]; ++k) {
        packOffset(imports, exports, selfReceiveOffset, plan.getIndicesTo()[selfIndex]*numPackets, numPackets);
        ++selfIndex;
        selfReceiveOffset += numPackets;
      }
    }
  }
}

template <class ExpView, class ImpView>
void DistributorActor::doPostsKokkos(const DistributorPlan& plan,
                                     const ExpView& exports,
                                     size_t numPackets,
                                     const ImpView& imports)
{
  static_assert(areKokkosViews<ExpView, ImpView>,
      "Data arrays for DistributorActor::doPosts must be Kokkos::Views");
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::includesVerbLevel;
  using Teuchos::ireceive;
  using Teuchos::isend;
  using Teuchos::send;
  using Teuchos::TypeNameTraits;
  using Teuchos::typeName;
  using std::endl;
  using Kokkos::Compat::create_const_view;
  using Kokkos::Compat::create_view;
  using Kokkos::Compat::subview_offset;
  using Kokkos::Compat::deep_copy_offset;
  typedef Array<size_t>::size_type size_type;
  typedef ExpView exports_view_type;
  typedef ImpView imports_view_type;

#ifdef KOKKOS_ENABLE_CUDA
  static_assert
    (! std::is_same<typename ExpView::memory_space, Kokkos::CudaUVMSpace>::value &&
     ! std::is_same<typename ImpView::memory_space, Kokkos::CudaUVMSpace>::value,
     "Please do not use Tpetra::Distributor with UVM allocations.  "
     "See Trilinos GitHub issue #1088.");
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_SYCL
    static_assert
      (! std::is_same<typename ExpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value &&
       ! std::is_same<typename ImpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value,
       "Please do not use Tpetra::Distributor with SharedUSM allocations.  "
       "See Trilinos GitHub issue #1088 (corresponding to CUDA).");
#endif // KOKKOS_ENABLE_SYCL

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::TimeMonitor timeMon (*timer_doPosts3KV_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

  const int myRank = plan.getComm()->getRank ();
  // Run-time configurable parameters that come from the input
  // ParameterList set by setParameterList().
  const Details::EDistributorSendType sendType = plan.getSendType();

//clang-format on
#if defined(HAVE_TPETRA_MPI)
  //  All-to-all communication layout is quite different from
  //  point-to-point, so we handle it separately.

  if (sendType == Details::DISTRIBUTOR_ALLTOALL) {
    doPostsAllToAllKokkos(plan, exports,numPackets, imports);
    return;
  }
#ifdef HAVE_TPETRACORE_MPI_ADVANCE
  else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL) {
    doPostsAllToAllKokkos(plan, exports,numPackets, imports);
    return;
  } else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV) {
    doPostsNbrAllToAllVKokkos(plan, exports,numPackets, imports);
    return;
  }
#endif // defined(HAVE_TPETRACORE_MPI_ADVANCE)
// clang-format off
  
#else // HAVE_TPETRA_MPI
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
                       plan.getStartsTo()[0]*numPackets,
                       plan.getLengthsTo()[0]*numPackets);
    }
    // should we just return here?
    // likely not as comm could be a serial comm
#endif // HAVE_TPETRA_MPI

  size_t selfReceiveOffset = 0;

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION
    (requests_.size () != 0,
     std::logic_error,
     "Tpetra::Distributor::doPosts(3 args, Kokkos): Process "
     << myRank << ": requests_.size() = " << requests_.size () << " != 0.");
#endif // HAVE_TPETRA_DEBUG

  // Distributor uses requests_.size() as the number of outstanding
  // nonblocking message requests, so we resize to zero to maintain
  // this invariant.
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
  const size_type actualNumReceives = as<size_type> (plan.getNumReceives()) +
    as<size_type> (plan.hasSelfMessage() ? 1 : 0);
  requests_.resize (0);

  // Post the nonblocking receives.  It's common MPI wisdom to post
  // receives before sends.  In MPI terms, this means favoring
  // adding to the "posted queue" (of receive requests) over adding
  // to the "unexpected queue" (of arrived messages not yet matched
  // with a receive).
  {
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonRecvs (*timer_doPosts3KV_recvs_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    size_t curBufferOffset = 0;
    for (size_type i = 0; i < actualNumReceives; ++i) {
      const size_t curBufLen = plan.getLengthsFrom()[i] * numPackets;
      if (plan.getProcsFrom()[i] != myRank) {
        // If my process is receiving these packet(s) from another
        // process (not a self-receive):
        //
        // 1. Set up the persisting view (recvBuf) of the imports
        //    array, given the offset and size (total number of
        //    packets from process getProcsFrom()[i]).
        // 2. Start the Irecv and save the resulting request.
        TEUCHOS_TEST_FOR_EXCEPTION(
            curBufferOffset + curBufLen > static_cast<size_t> (imports.size ()),
            std::logic_error, "Tpetra::Distributor::doPosts(3 args, Kokkos): "
            "Exceeded size of 'imports' array in packing loop on Process " <<
            myRank << ".  imports.size() = " << imports.size () << " < "
            "curBufferOffset(" << curBufferOffset << ") + curBufLen(" <<
            curBufLen << ").");
        imports_view_type recvBuf =
          subview_offset (imports, curBufferOffset, curBufLen);
        requests_.push_back (ireceive<int> (recvBuf, plan.getProcsFrom()[i],
              mpiTag_, *plan.getComm()));
      }
      else { // Receiving from myself
        selfReceiveOffset = curBufferOffset; // Remember the self-recv offset
      }
      curBufferOffset += curBufLen;
    }
  }

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::TimeMonitor timeMonSends (*timer_doPosts3KV_sends_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

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

  size_t selfNum = 0;
  size_t selfIndex = 0;

  if (plan.getIndicesTo().is_null()) {

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonSends2 (*timer_doPosts3KV_sends_fast_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    // Data are already blocked (laid out) by process, so we don't
    // need a separate send buffer (besides the exports array).
    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myRank) {
        exports_view_type tmpSend = subview_offset(
            exports, plan.getStartsTo()[p]*numPackets, plan.getLengthsTo()[p]*numPackets);

        if (sendType == Details::DISTRIBUTOR_ISEND) {
          // NOTE: This looks very similar to the tmpSend above, but removing
          // tmpSendBuf and uses tmpSend leads to a performance hit on Arm
          // SerialNode builds
          exports_view_type tmpSendBuf =
            subview_offset (exports, plan.getStartsTo()[p] * numPackets,
                plan.getLengthsTo()[p] * numPackets);
          requests_.push_back (isend<int> (tmpSendBuf, plan.getProcsTo()[p],
                mpiTag_, *plan.getComm()));
        }
        else {  // DISTRIBUTOR_SEND
          send<int> (tmpSend,
              as<int> (tmpSend.size ()),
              plan.getProcsTo()[p], mpiTag_, *plan.getComm());
        }
      }
      else { // "Sending" the message to myself
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
          plan.getStartsTo()[selfNum]*numPackets,
          plan.getLengthsTo()[selfNum]*numPackets);
    }

  }
  else { // data are not blocked by proc, use send buffer

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonSends2 (*timer_doPosts3KV_sends_slow_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    typedef typename ExpView::non_const_value_type Packet;
    typedef typename ExpView::array_layout Layout;
    typedef typename ExpView::device_type Device;
    typedef typename ExpView::memory_traits Mem;

    // This buffer is long enough for only one message at a time.
    // Thus, we use DISTRIBUTOR_SEND always in this case, regardless
    // of sendType requested by user. 
    // This code path formerly errored out with message:
    //     Tpetra::Distributor::doPosts(3 args, Kokkos): 
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

    Kokkos::View<Packet*,Layout,Device,Mem> sendArray ("sendArray",
        plan.getMaxSendLength() * numPackets);

    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myRank) {
        size_t sendArrayOffset = 0;
        size_t j = plan.getStartsTo()[p];
        for (size_t k = 0; k < plan.getLengthsTo()[p]; ++k, ++j) {
          packOffset(sendArray, exports, sendArrayOffset, plan.getIndicesTo()[j]*numPackets, numPackets);
          sendArrayOffset += numPackets;
        }
        typename ExpView::execution_space().fence();

        ImpView tmpSend =
          subview_offset(sendArray, size_t(0), plan.getLengthsTo()[p]*numPackets);

        send<int> (tmpSend,
            as<int> (tmpSend.size ()),
            plan.getProcsTo()[p], mpiTag_, *plan.getComm());
      }
      else { // "Sending" the message to myself
        selfNum = p;
        selfIndex = plan.getStartsTo()[p];
      }
    }

    if (plan.hasSelfMessage()) {
      for (size_t k = 0; k < plan.getLengthsTo()[selfNum]; ++k) {
        packOffset(imports, exports, selfReceiveOffset, plan.getIndicesTo()[selfIndex]*numPackets, numPackets);
        ++selfIndex;
        selfReceiveOffset += numPackets;
      }
    }
  }
}

// clang-format on
#ifdef HAVE_TPETRA_MPI
template <class ExpView, class ImpView>
void DistributorActor::doPostsAllToAll(
    const DistributorPlan &plan, const ExpView &exports,
    const Teuchos::ArrayView<const size_t> &numExportPacketsPerLID,
    const ImpView &imports,
    const Teuchos::ArrayView<const size_t> &numImportPacketsPerLID) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");

  using size_type = Teuchos::Array<size_t>::size_type;

  auto comm = plan.getComm();
  std::vector<int> sendcounts(comm->getSize(), 0);
  std::vector<int> sdispls(comm->getSize(), 0);
  std::vector<int> recvcounts(comm->getSize(), 0);
  std::vector<int> rdispls(comm->getSize(), 0);

  size_t curPKToffset = 0;
  for (size_t pp = 0; pp < plan.getNumSends(); ++pp) {
    sdispls[plan.getProcsTo()[pp]] = curPKToffset;
    size_t numPackets = 0;
    for (size_t j = plan.getStartsTo()[pp];
         j < plan.getStartsTo()[pp] + plan.getLengthsTo()[pp]; ++j) {
      numPackets += numExportPacketsPerLID[j];
    }
    // numPackets is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(numPackets > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(4 args, Kokkos): "
                               "Send count for send "
                                   << pp << " (" << numPackets
                                   << ") is too large "
                                      "to be represented as int.");
    sendcounts[plan.getProcsTo()[pp]] = static_cast<int>(numPackets);
    curPKToffset += numPackets;
  }

  const size_type actualNumReceives =
      Teuchos::as<size_type>(plan.getNumReceives()) +
      Teuchos::as<size_type>(plan.hasSelfMessage() ? 1 : 0);

  size_t curBufferOffset = 0;
  size_t curLIDoffset = 0;
  for (size_type i = 0; i < actualNumReceives; ++i) {
    size_t totalPacketsFrom_i = 0;
    for (size_t j = 0; j < plan.getLengthsFrom()[i]; ++j) {
      totalPacketsFrom_i += numImportPacketsPerLID[curLIDoffset + j];
    }
    curLIDoffset += plan.getLengthsFrom()[i];

    rdispls[plan.getProcsFrom()[i]] = curBufferOffset;
    // totalPacketsFrom_i is converted down to int, so make sure it can be
    // represented
    TEUCHOS_TEST_FOR_EXCEPTION(totalPacketsFrom_i > size_t(INT_MAX),
                               std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Recv count for receive "
                                   << i << " (" << totalPacketsFrom_i
                                   << ") is too large "
                                      "to be represented as int.");
    recvcounts[plan.getProcsFrom()[i]] = static_cast<int>(totalPacketsFrom_i);
    curBufferOffset += totalPacketsFrom_i;
  }

  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiComm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm);
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> rawComm =
      mpiComm->getRawMpiComm();
  using T = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  if (Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL == plan.getSendType()) {
    MPIX_Comm *mpixComm = *plan.getMPIXComm();
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
#endif // HAVE_TPETRACORE_MPI_ADVANCE

  const int err = MPI_Alltoallv(
      exports.data(), sendcounts.data(), sdispls.data(), rawType,
      imports.data(), recvcounts.data(), rdispls.data(), rawType, (*rawComm)());

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "MPI_Alltoallv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");
}

template <class ExpView, class ExpPacketsView, class ImpView, class ImpPacketsView>
void DistributorActor::doPostsAllToAllKokkos(
    const DistributorPlan &plan, const ExpView &exports,
    const ExpPacketsView &numExportPacketsPerLID,
    const ImpView &imports,
    const ImpPacketsView &numImportPacketsPerLID) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");

  using size_type = Teuchos::Array<size_t>::size_type;
  using ExpExecSpace = typename ExpView::execution_space;
  using ImpExecSpace = typename ImpView::execution_space;

  auto comm = plan.getComm();
  Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> sendcounts("sendcounts", comm->getSize());
  Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> sdispls("sdispls", comm->getSize());
  Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> recvcounts("recvcounts", comm->getSize());
  Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> rdispls("rdispls", comm->getSize());

  auto sendcounts_d = Kokkos::create_mirror_view(ExpExecSpace(), sendcounts);
  auto sdispls_d = Kokkos::create_mirror_view(ExpExecSpace(), sdispls);
  auto recvcounts_d = Kokkos::create_mirror_view(ImpExecSpace(), recvcounts);
  auto rdispls_d = Kokkos::create_mirror_view(ImpExecSpace(), rdispls);

  auto getStartsTo = Kokkos::Compat::getKokkosViewDeepCopy<ExpExecSpace>(plan.getStartsTo());
  auto getLengthsTo = Kokkos::Compat::getKokkosViewDeepCopy<ExpExecSpace>(plan.getLengthsTo());
  auto getProcsTo = Kokkos::Compat::getKokkosViewDeepCopy<ExpExecSpace>(plan.getProcsTo());

  size_t curPKToffset = 0;
  Kokkos::parallel_scan("curPKToffsetAllToAll", Kokkos::RangePolicy<ExpExecSpace>(0, plan.getNumSends()), KOKKOS_LAMBDA(const size_t pp, size_t& offset, bool is_final) { 
    if (is_final) sdispls_d(getProcsTo(pp)) = offset;
    size_t numPackets = 0;
    for (size_t j = getStartsTo(pp); j < getStartsTo(pp) + getLengthsTo(pp); ++j) {
      numPackets += numExportPacketsPerLID(j);
    }
    if (is_final) sendcounts_d(getProcsTo(pp)) = static_cast<int>(numPackets);
    offset += numPackets;
  }, curPKToffset);

  int overflow;
  Kokkos::parallel_reduce("overflowAllToAll", Kokkos::RangePolicy<ExpExecSpace>(0, plan.getNumSends()), KOKKOS_LAMBDA(const size_t pp, int& index) {
    if(sendcounts_d(getProcsTo(pp)) < 0) {
      index = pp+1;
    }
  }, overflow);
  
  // numPackets is converted down to int, so make sure it can be represented
  TEUCHOS_TEST_FOR_EXCEPTION(overflow, std::logic_error,
                             "Tpetra::Distributor::doPostsKokkos(4 args, Kokkos): "
                             "Send count for send "
                                 << overflow-1 << " is too large "
                                    "to be represented as int.");

  const size_type actualNumReceives =
      Teuchos::as<size_type>(plan.getNumReceives()) +
      Teuchos::as<size_type>(plan.hasSelfMessage() ? 1 : 0);

  auto getLengthsFrom = Kokkos::Compat::getKokkosViewDeepCopy<ImpExecSpace>(plan.getLengthsFrom());
  auto getProcsFrom = Kokkos::Compat::getKokkosViewDeepCopy<ImpExecSpace>(plan.getProcsFrom());
  
  Kokkos::View<size_t*, ImpExecSpace> curLIDoffset("curLIDoffset", actualNumReceives);
  Kokkos::parallel_scan("curLIDoffsetAllToAll", Kokkos::RangePolicy<ImpExecSpace>(0, actualNumReceives), KOKKOS_LAMBDA(const size_type i, size_t& offset, bool is_final) {
    if(is_final) curLIDoffset(i) = offset;
    offset += getLengthsFrom(i); 
  });

  Kokkos::parallel_scan("receiveAllToAll", Kokkos::RangePolicy<ImpExecSpace>(0, actualNumReceives), KOKKOS_LAMBDA(const size_type i, size_t& curBufferOffset, bool is_final) {
    size_t totalPacketsFrom_i = 0;
    for(size_t j = 0; j < getLengthsFrom(i); j++) {
      totalPacketsFrom_i += numImportPacketsPerLID(curLIDoffset(i) + j);
    }

    if(is_final) rdispls_d(getProcsFrom(i)) = curBufferOffset;
    if(is_final) recvcounts_d(getProcsFrom(i)) = static_cast<int>(totalPacketsFrom_i);
    curBufferOffset += totalPacketsFrom_i;
  });
  
  Kokkos::parallel_reduce("overflowAllToAll2", Kokkos::RangePolicy<ExpExecSpace>(0, actualNumReceives), KOKKOS_LAMBDA(const size_type i, int& index) {
    if(recvcounts_d(getProcsFrom(i)) < 0) {
      index = i+1;
    }
  }, overflow);
  
  // totalPacketsFrom_i is converted down to int, so make sure it can be
  // represented
  TEUCHOS_TEST_FOR_EXCEPTION(overflow, std::logic_error,
                             "Tpetra::Distributor::doPostsKokkos(4 args, Kokkos): "
                             "Recv count for receive "
                                 << overflow-1 << " is too large "
                                    "to be represented as int.");
  
  Kokkos::deep_copy(sendcounts, sendcounts_d);
  Kokkos::deep_copy(sdispls, sdispls_d);
  Kokkos::deep_copy(recvcounts, recvcounts_d);
  Kokkos::deep_copy(rdispls, rdispls_d);
  
  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiComm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm);
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> rawComm =
      mpiComm->getRawMpiComm();
  using T = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  if (Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL == plan.getSendType()) {
    MPIX_Comm *mpixComm = *plan.getMPIXComm();
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
#endif // HAVE_TPETRACORE_MPI_ADVANCE

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
void DistributorActor::doPostsNbrAllToAllV(
    const DistributorPlan &plan, const ExpView &exports,
    const Teuchos::ArrayView<const size_t> &numExportPacketsPerLID,
    const ImpView &imports,
    const Teuchos::ArrayView<const size_t> &numImportPacketsPerLID) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");

  const Teuchos_Ordinal numSends = plan.getProcsTo().size();
  const Teuchos_Ordinal numRecvs = plan.getProcsFrom().size();

  auto comm = plan.getComm();
  std::vector<int> sendcounts(numSends, 0);
  std::vector<int> sdispls(numSends, 0);
  std::vector<int> recvcounts(numRecvs, 0);
  std::vector<int> rdispls(numRecvs, 0);

  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiComm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm);
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> rawComm =
      mpiComm->getRawMpiComm();
  using T = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());

  // unlike standard alltoall, entry `i` in sdispls and sendcounts
  // refer to the ith participating rank, rather than rank i
  size_t curPKToffset = 0;
  for (Teuchos_Ordinal pp = 0; pp < numSends; ++pp) {
    sdispls[pp] = curPKToffset;
    size_t numPackets = 0;
    for (size_t j = plan.getStartsTo()[pp];
         j < plan.getStartsTo()[pp] + plan.getLengthsTo()[pp]; ++j) {
      numPackets += numExportPacketsPerLID[j];
    }
    // numPackets is converted down to int, so make sure it can be represented
    TEUCHOS_TEST_FOR_EXCEPTION(numPackets > size_t(INT_MAX), std::logic_error,
                               "Tpetra::Distributor::doPosts(4 args, Kokkos): "
                               "Send count for send "
                                   << pp << " (" << numPackets
                                   << ") is too large "
                                      "to be represented as int.");
    sendcounts[pp] = static_cast<int>(numPackets);
    curPKToffset += numPackets;
  }
  size_t curBufferOffset = 0;
  size_t curLIDoffset = 0;
  for (Teuchos_Ordinal i = 0; i < numRecvs; ++i) {
    size_t totalPacketsFrom_i = 0;
    for (size_t j = 0; j < plan.getLengthsFrom()[i]; ++j) {
      totalPacketsFrom_i += numImportPacketsPerLID[curLIDoffset + j];
    }
    curLIDoffset += plan.getLengthsFrom()[i];

    rdispls[i] = curBufferOffset;
    // totalPacketsFrom_i is converted down to int, so make sure it can be
    // represented
    TEUCHOS_TEST_FOR_EXCEPTION(totalPacketsFrom_i > size_t(INT_MAX),
                               std::logic_error,
                               "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                               "Recv count for receive "
                                   << i << " (" << totalPacketsFrom_i
                                   << ") is too large "
                                      "to be represented as int.");
    recvcounts[i] = static_cast<int>(totalPacketsFrom_i);
    curBufferOffset += totalPacketsFrom_i;
  }

  MPIX_Comm *mpixComm = *plan.getMPIXComm();
  const int err = MPIX_Neighbor_alltoallv(
      exports.data(), sendcounts.data(), sdispls.data(), rawType,
      imports.data(), recvcounts.data(), rdispls.data(), rawType, mpixComm);

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "MPIX_Neighbor_alltoallv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");
}

template <class ExpView, class ExpPacketsView, class ImpView, class ImpPacketsView>
void DistributorActor::doPostsNbrAllToAllVKokkos(
    const DistributorPlan &plan, const ExpView &exports,
    const ExpPacketsView &numExportPacketsPerLID,
    const ImpView &imports,
    const ImpPacketsView &numImportPacketsPerLID) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      !plan.getIndicesTo().is_null(), std::runtime_error,
      "Send Type=\"Alltoall\" only works for fast-path communication.");
  
  using ExpExecSpace = typename ExpView::execution_space;
  using ImpExecSpace = typename ImpView::execution_space;

  const Teuchos_Ordinal numSends = plan.getProcsTo().size();
  const Teuchos_Ordinal numRecvs = plan.getProcsFrom().size();

  auto comm = plan.getComm();
  Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> sendcounts("sendcounts", numSends);
  Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> sdispls("sdispls", numSends);
  Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> recvcounts("recvcounts", numRecvs);
  Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> rdispls("rdispls", numRecvs);
  
  auto sendcounts_d = Kokkos::create_mirror_view(ExpExecSpace(), sendcounts);
  auto sdispls_d = Kokkos::create_mirror_view(ExpExecSpace(), sdispls);
  auto recvcounts_d = Kokkos::create_mirror_view(ImpExecSpace(), recvcounts);
  auto rdispls_d = Kokkos::create_mirror_view(ImpExecSpace(), rdispls);

  auto getStartsTo = Kokkos::Compat::getKokkosViewDeepCopy<ExpExecSpace>(plan.getStartsTo());
  auto getLengthsTo = Kokkos::Compat::getKokkosViewDeepCopy<ExpExecSpace>(plan.getLengthsTo());

  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiComm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm);
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> rawComm =
      mpiComm->getRawMpiComm();
  using T = typename ExpView::non_const_value_type;
  MPI_Datatype rawType = ::Tpetra::Details::MpiTypeTraits<T>::getType(T());
  
  // unlike standard alltoall, entry `i` in sdispls and sendcounts
  // refer to the ith participating rank, rather than rank i
  Kokkos::parallel_scan("sendNbrAllToAllV", Kokkos::RangePolicy<ExpExecSpace>(0, numSends), KOKKOS_LAMBDA(const Teuchos_Ordinal pp, size_t& curPKToffset, bool is_final) { 
    if (is_final) sdispls_d(pp) = curPKToffset;
    size_t numPackets = 0;
    for (size_t j = getStartsTo(pp); j < getStartsTo(pp) + getLengthsTo(pp); ++j) {
      numPackets += numExportPacketsPerLID(j);
    }
    if (is_final) sendcounts_d(pp) = static_cast<int>(numPackets);
    curPKToffset += numPackets;
  });
  
  int overflow;
  Kokkos::parallel_reduce("overflowNbrAllToAllV", Kokkos::RangePolicy<ExpExecSpace>(0, numSends), KOKKOS_LAMBDA(const Teuchos_Ordinal pp, int& index) {
    if(sendcounts_d(pp) < 0) {
      index = i+1;
    }
  }, overflow);
    
  // numPackets is converted down to int, so make sure it can be represented
  TEUCHOS_TEST_FOR_EXCEPTION(overflow, std::logic_error,
                             "Tpetra::Distributor::doPostsKokkos(4 args, Kokkos): "
                             "Send count for send "
                                 << overflow-1 << " is too large "
                                    "to be represented as int.");
  
  auto getLengthsFrom = Kokkos::Compat::getKokkosViewDeepCopy<ImpExecSpace>(plan.getLengthsFrom());
  
  Kokkos::View<size_t*, ImpExecSpace> curLIDoffset("curLIDoffset", numRecvs);
  Kokkos::parallel_scan("curLIDoffsetNbrAllToAllV", Kokkos::RangePolicy<ImpExecSpace>(0, numRecvs), KOKKOS_LAMBDA(const Teuchos_Ordinal i, size_t& offset, bool is_final) {
    if(is_final) curLIDoffset(i) = offset;
    offset += getLengthsFrom(i); 
  });
  
  Kokkos::parallel_scan("receiveNbrAllToAllV", Kokkos::RangePolicy<ImpExecSpace>(0, numRecvs), KOKKOS_LAMBDA(const Teuchos_Ordinal i, size_t& curBufferOffset, bool is_final) {
    if(is_final) rdispls_d(i) = curBufferOffset;
    size_t totalPacketsFrom_i = 0;
    for(size_t j = 0; j < getLengthsFrom(i); j++) {
      totalPacketsFrom_i += numImportPacketsPerLID(curLIDoffset(i) + j);
    }
  
    recvcounts_d(i) = static_cast<int>(totalPacketsFrom_i);
    curBufferOffset += totalPacketsFrom_i;
  });
  
  Kokkos::parallel_reduce("overflowNbrAllToAllV", Kokkos::RangePolicy<ImpExecSpace>(0, numRecvs), KOKKOS_LAMBDA(const Teuchos_Ordinal i, int& index) {
    if(recvcounts_d(pp) < 0) {
      index = i+1;
    }
  }, overflow);
  
  // totalPacketsFrom_i is converted down to int, so make sure it can be
  // represented
  TEUCHOS_TEST_FOR_EXCEPTION(overflow, std::logic_error,
                             "Tpetra::Distributor::doPostsKokkos(4 args, Kokkos): "
                             "Recv count for receive "
                                 << overflow-1 << ") is too large "
                                    "to be represented as int.");
  
  Kokkos::deep_copy(sendcounts, sendcounts_d);
  Kokkos::deep_copy(sdispls, sdispls_d);
  Kokkos::deep_copy(recvcounts, recvcounts_d);
  Kokkos::deep_copy(rdispls, rdispls_d);

  MPIX_Comm *mpixComm = *plan.getMPIXComm();
  const int err = MPIX_Neighbor_alltoallv(
      exports.data(), sendcounts.data(), sdispls.data(), rawType,
      imports.data(), recvcounts.data(), rdispls.data(), rawType, mpixComm);

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                             "MPIX_Neighbor_alltoallv failed with error \""
                                 << Teuchos::mpiErrorCodeToString(err)
                                 << "\".");
}
#endif // HAVE_TPETRACORE_MPI_ADVANCE
#endif // HAVE_TPETRA_MPI
       // clang-format off

template <class ExpView, class ImpView>
void DistributorActor::doPosts(const DistributorPlan& plan,
                               const ExpView &exports,
                               const Teuchos::ArrayView<const size_t>& numExportPacketsPerLID,
                               const ImpView &imports,
                               const Teuchos::ArrayView<const size_t>& numImportPacketsPerLID)
{
  static_assert(areKokkosViews<ExpView, ImpView>,
      "Data arrays for DistributorActor::doPosts must be Kokkos::Views");
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::ireceive;
  using Teuchos::isend;
  using Teuchos::send;
  using Teuchos::TypeNameTraits;
  using std::endl;
  using Kokkos::Compat::create_const_view;
  using Kokkos::Compat::create_view;
  using Kokkos::Compat::subview_offset;
  using Kokkos::Compat::deep_copy_offset;
  typedef Array<size_t>::size_type size_type;
  typedef ExpView exports_view_type;
  typedef ImpView imports_view_type;

#ifdef KOKKOS_ENABLE_CUDA
  static_assert (! std::is_same<typename ExpView::memory_space, Kokkos::CudaUVMSpace>::value &&
                 ! std::is_same<typename ImpView::memory_space, Kokkos::CudaUVMSpace>::value,
                 "Please do not use Tpetra::Distributor with UVM "
                 "allocations.  See GitHub issue #1088.");
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_SYCL
    static_assert (! std::is_same<typename ExpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value &&
                   ! std::is_same<typename ImpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value,
                   "Please do not use Tpetra::Distributor with SharedUSM "
                   "allocations.  See GitHub issue #1088 (corresponding to CUDA).");
#endif // KOKKOS_ENABLE_SYCL

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::TimeMonitor timeMon (*timer_doPosts4KV_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

  // Run-time configurable parameters that come from the input
  // ParameterList set by setParameterList().
  const Details::EDistributorSendType sendType = plan.getSendType();

#ifdef HAVE_TPETRA_MPI
  //  All-to-all communication layout is quite different from
  //  point-to-point, so we handle it separately.
  if (sendType == Details::DISTRIBUTOR_ALLTOALL) {
    doPostsAllToAll(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
    return;
  }
#ifdef HAVE_TPETRACORE_MPI_ADVANCE
  else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL)
  {
    doPostsAllToAll(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
    return;
  } else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV) {
    doPostsNbrAllToAllV(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
    return;
  }
#endif

#else // HAVE_TPETRA_MPI
    if (plan.hasSelfMessage()) {

      size_t selfReceiveOffset = 0;

      // setup arrays containing starting-offsets into exports for each send,
      // and num-packets-to-send for each send.
      Array<size_t> sendPacketOffsets(plan.getNumSends(),0), packetsPerSend(plan.getNumSends(),0);
      size_t maxNumPackets = 0;
      size_t curPKToffset = 0;
      for (size_t pp=0; pp<plan.getNumSends(); ++pp) {
        sendPacketOffsets[pp] = curPKToffset;
        size_t numPackets = 0;
        for (size_t j=plan.getStartsTo()[pp]; j<plan.getStartsTo()[pp]+plan.getLengthsTo()[pp]; ++j) {
          numPackets += numExportPacketsPerLID[j];
        }
        if (numPackets > maxNumPackets) maxNumPackets = numPackets;
        packetsPerSend[pp] = numPackets;
        curPKToffset += numPackets;
      }

      deep_copy_offset(imports, exports, selfReceiveOffset,
          sendPacketOffsets[0], packetsPerSend[0]);
    }
#endif // HAVE_TPETRA_MPI

  const int myProcID = plan.getComm()->getRank ();
  size_t selfReceiveOffset = 0;

#ifdef HAVE_TPETRA_DEBUG
  // Different messages may have different numbers of packets.
  size_t totalNumImportPackets = 0;
  for (size_type ii = 0; ii < numImportPacketsPerLID.size (); ++ii) {
    totalNumImportPackets += numImportPacketsPerLID[ii];
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
      imports.extent (0) < totalNumImportPackets, std::runtime_error,
      "Tpetra::Distributor::doPosts(4 args, Kokkos): The 'imports' array must have "
      "enough entries to hold the expected number of import packets.  "
      "imports.extent(0) = " << imports.extent (0) << " < "
      "totalNumImportPackets = " << totalNumImportPackets << ".");
  TEUCHOS_TEST_FOR_EXCEPTION
    (requests_.size () != 0, std::logic_error, "Tpetra::Distributor::"
     "doPosts(4 args, Kokkos): Process " << myProcID << ": requests_.size () = "
     << requests_.size () << " != 0.");
#endif // HAVE_TPETRA_DEBUG
  // Distributor uses requests_.size() as the number of outstanding
  // nonblocking message requests, so we resize to zero to maintain
  // this invariant.
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
  const size_type actualNumReceives = as<size_type> (plan.getNumReceives()) +
    as<size_type> (plan.hasSelfMessage() ? 1 : 0);
  requests_.resize (0);

  // Post the nonblocking receives.  It's common MPI wisdom to post
  // receives before sends.  In MPI terms, this means favoring
  // adding to the "posted queue" (of receive requests) over adding
  // to the "unexpected queue" (of arrived messages not yet matched
  // with a receive).
  {
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonRecvs (*timer_doPosts4KV_recvs_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    size_t curBufferOffset = 0;
    size_t curLIDoffset = 0;
    for (size_type i = 0; i < actualNumReceives; ++i) {
      size_t totalPacketsFrom_i = 0;
      for (size_t j = 0; j < plan.getLengthsFrom()[i]; ++j) {
        totalPacketsFrom_i += numImportPacketsPerLID[curLIDoffset+j];
      }
      // totalPacketsFrom_i is converted down to int, so make sure it can be represented
      TEUCHOS_TEST_FOR_EXCEPTION(totalPacketsFrom_i > size_t(INT_MAX),
                                 std::logic_error, "Tpetra::Distributor::doPosts(3 args, Kokkos): "
                                 "Recv count for receive " << i << " (" << totalPacketsFrom_i << ") is too large "
                                 "to be represented as int.");
      curLIDoffset += plan.getLengthsFrom()[i];
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
          subview_offset (imports, curBufferOffset, totalPacketsFrom_i);
        requests_.push_back (ireceive<int> (recvBuf, plan.getProcsFrom()[i],
              mpiTag_, *plan.getComm()));
      }
      else { // Receiving these packet(s) from myself
        selfReceiveOffset = curBufferOffset; // Remember the offset
      }
      curBufferOffset += totalPacketsFrom_i;
    }
  }

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::TimeMonitor timeMonSends (*timer_doPosts4KV_sends_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

  // setup arrays containing starting-offsets into exports for each send,
  // and num-packets-to-send for each send.
  Array<size_t> sendPacketOffsets(plan.getNumSends(),0), packetsPerSend(plan.getNumSends(),0);
  size_t maxNumPackets = 0;
  size_t curPKToffset = 0;
  for (size_t pp=0; pp<plan.getNumSends(); ++pp) {
    sendPacketOffsets[pp] = curPKToffset;
    size_t numPackets = 0;
    for (size_t j=plan.getStartsTo()[pp]; j<plan.getStartsTo()[pp]+plan.getLengthsTo()[pp]; ++j) {
      numPackets += numExportPacketsPerLID[j];
    }
    if (numPackets > maxNumPackets) maxNumPackets = numPackets;
    // numPackets will be used as a message length, so make sure it can be represented as int
    TEUCHOS_TEST_FOR_EXCEPTION(numPackets > size_t(INT_MAX),
                               std::logic_error, "Tpetra::Distributor::doPosts(4 args, Kokkos): "
                               "packetsPerSend[" << pp << "] = " << numPackets << " is too large "
                               "to be represented as int.");
    packetsPerSend[pp] = numPackets;
    curPKToffset += numPackets;
  }

  // setup scan through getProcsTo() list starting with higher numbered procs
  // (should help balance message traffic)
  size_t numBlocks = plan.getNumSends() + plan.hasSelfMessage();
  size_t procIndex = 0;
  while ((procIndex < numBlocks) && (plan.getProcsTo()[procIndex] < myProcID)) {
    ++procIndex;
  }
  if (procIndex == numBlocks) {
    procIndex = 0;
  }

  size_t selfNum = 0;
  size_t selfIndex = 0;
  if (plan.getIndicesTo().is_null()) {

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonSends2 (*timer_doPosts4KV_sends_fast_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    // Data are already blocked (laid out) by process, so we don't
    // need a separate send buffer (besides the exports array).
    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myProcID && packetsPerSend[p] > 0) {
        exports_view_type tmpSend =
          subview_offset(exports, sendPacketOffsets[p], packetsPerSend[p]);

        if (sendType == Details::DISTRIBUTOR_ISEND) {
          exports_view_type tmpSendBuf =
            subview_offset (exports, sendPacketOffsets[p], packetsPerSend[p]);
          requests_.push_back (isend<int> (tmpSendBuf, plan.getProcsTo()[p],
                mpiTag_, *plan.getComm()));
        }
        else { // DISTRIBUTOR_SEND
          send<int> (tmpSend,
              as<int> (tmpSend.size ()),
              plan.getProcsTo()[p], mpiTag_, *plan.getComm());
        }
      }
      else { // "Sending" the message to myself
        selfNum = p;
      }
    }

    if (plan.hasSelfMessage()) {
      deep_copy_offset(imports, exports, selfReceiveOffset,
          sendPacketOffsets[selfNum], packetsPerSend[selfNum]);
    }
  }
  else { // data are not blocked by proc, use send buffer

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonSends2 (*timer_doPosts4KV_sends_slow_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    // FIXME (mfh 05 Mar 2013) This may be broken for Isend.
    typedef typename ExpView::non_const_value_type Packet;
    typedef typename ExpView::array_layout Layout;
    typedef typename ExpView::device_type Device;
    typedef typename ExpView::memory_traits Mem;

    // This buffer is long enough for only one message at a time.
    // Thus, we use DISTRIBUTOR_SEND always in this case, regardless
    // of sendType requested by user.
    // This code path formerly errored out with message:
    //     Tpetra::Distributor::doPosts(4-arg, Kokkos):
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

    Kokkos::View<Packet*,Layout,Device,Mem> sendArray ("sendArray",
                                                        maxNumPackets);

    Array<size_t> indicesOffsets (numExportPacketsPerLID.size(), 0);
    size_t ioffset = 0;
    for (int j=0; j<numExportPacketsPerLID.size(); ++j) {
      indicesOffsets[j] = ioffset;
      ioffset += numExportPacketsPerLID[j];
    }

    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myProcID) {
        size_t sendArrayOffset = 0;
        size_t j = plan.getStartsTo()[p];
        size_t numPacketsTo_p = 0;
        for (size_t k = 0; k < plan.getLengthsTo()[p]; ++k, ++j) {
          numPacketsTo_p += numExportPacketsPerLID[j];
          deep_copy_offset(sendArray, exports, sendArrayOffset,
              indicesOffsets[j], numExportPacketsPerLID[j]);
          sendArrayOffset += numExportPacketsPerLID[j];
        }
        typename ExpView::execution_space().fence();

        if (numPacketsTo_p > 0) {
          ImpView tmpSend =
            subview_offset(sendArray, size_t(0), numPacketsTo_p);

          send<int> (tmpSend,
              as<int> (tmpSend.size ()),
              plan.getProcsTo()[p], mpiTag_, *plan.getComm());
        }
      }
      else { // "Sending" the message to myself
        selfNum = p;
        selfIndex = plan.getStartsTo()[p];
      }
    }

    if (plan.hasSelfMessage()) {
      for (size_t k = 0; k < plan.getLengthsTo()[selfNum]; ++k) {
        deep_copy_offset(imports, exports, selfReceiveOffset,
            indicesOffsets[selfIndex],
            numExportPacketsPerLID[selfIndex]);
        selfReceiveOffset += numExportPacketsPerLID[selfIndex];
        ++selfIndex;
      }
    }
  }
}

template <class ExpView, class ExpPacketsView, class ImpView, class ImpPacketsView>
void DistributorActor::doPostsKokkos(const DistributorPlan& plan,
                                     const ExpView &exports,
                                     const ExpPacketsView &numExportPacketsPerLID,
                                     const ImpView &imports,
                                     const ImpPacketsView &numImportPacketsPerLID)
{
  static_assert(areKokkosViews<ExpView, ImpView>,
      "Data arrays for DistributorActor::doPostsKokkos must be Kokkos::Views");
  static_assert(areKokkosViews<ExpPacketsView, ImpPacketsView>,
      "Num packets arrays for DistributorActor::doPostsKokkos must be Kokkos::Views");
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::ireceive;
  using Teuchos::isend;
  using Teuchos::send;
  using Teuchos::TypeNameTraits;
  using std::endl;
  using Kokkos::Compat::create_const_view;
  using Kokkos::Compat::create_view;
  using Kokkos::Compat::subview_offset;
  using Kokkos::Compat::deep_copy_offset;
  using ExpExecSpace = typename ExpView::execution_space;
  using ImpExecSpace = typename ImpView::execution_space;
  typedef Array<size_t>::size_type size_type;
  typedef ExpView exports_view_type;
  typedef ImpView imports_view_type;

#ifdef KOKKOS_ENABLE_CUDA
  static_assert (! std::is_same<typename ExpView::memory_space, Kokkos::CudaUVMSpace>::value &&
                 ! std::is_same<typename ImpView::memory_space, Kokkos::CudaUVMSpace>::value,
                 "Please do not use Tpetra::Distributor with UVM "
                 "allocations.  See GitHub issue #1088.");
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_SYCL
    static_assert (! std::is_same<typename ExpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value &&
                   ! std::is_same<typename ImpView::memory_space, Kokkos::Experimental::SYCLSharedUSMSpace>::value,
                   "Please do not use Tpetra::Distributor with SharedUSM "
                   "allocations.  See GitHub issue #1088 (corresponding to CUDA).");
#endif // KOKKOS_ENABLE_SYCL

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::TimeMonitor timeMon (*timer_doPosts4KV_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

  // Run-time configurable parameters that come from the input
  // ParameterList set by setParameterList().
  const Details::EDistributorSendType sendType = plan.getSendType();

#ifdef HAVE_TPETRA_MPI
  //  All-to-all communication layout is quite different from
  //  point-to-point, so we handle it separately.
  if (sendType == Details::DISTRIBUTOR_ALLTOALL) {
    doPostsAllToAllKokkos(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
    return;
  }
#ifdef HAVE_TPETRACORE_MPI_ADVANCE
  else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_ALLTOALL)
  {
    doPostsAllToAllKokkos(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
    return;
  } else if (sendType == Details::DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV) {
    doPostsNbrAllToAllVKokkos(plan, exports, numExportPacketsPerLID, imports, numImportPacketsPerLID);
    return;
  }
#endif

#else // HAVE_TPETRA_MPI
    if (plan.hasSelfMessage()) {
      size_t selfReceiveOffset = 0;

      // setup arrays containing starting-offsets into exports for each send,
      // and num-packets-to-send for each send.
      Array<size_t> sendPacketOffsets(plan.getNumSends(),0), packetsPerSend(plan.getNumSends(),0);
      size_t maxNumPackets = 0;
      size_t curPKToffset = 0;
      for (size_t pp=0; pp<plan.getNumSends(); ++pp) {
        sendPacketOffsets[pp] = curPKToffset;
        size_t numPackets = 0;
        for (size_t j=plan.getStartsTo()[pp]; j<plan.getStartsTo()[pp]+plan.getLengthsTo()[pp]; ++j) {
          numPackets += numExportPacketsPerLID[j];
        }
        if (numPackets > maxNumPackets) maxNumPackets = numPackets;
        packetsPerSend[pp] = numPackets;
        curPKToffset += numPackets;
      }

      deep_copy_offset(imports, exports, selfReceiveOffset,
          sendPacketOffsets[0], packetsPerSend[0]);
    }
#endif // HAVE_TPETRA_MPI

  const int myProcID = plan.getComm()->getRank ();
  size_t selfReceiveOffset = 0;

#ifdef HAVE_TPETRA_DEBUG
  // Different messages may have different numbers of packets.
  size_t totalNumImportPackets = Kokkos::Experimental::reduce(ImpExecSpace(), numImportPacketsPerLID);
  TEUCHOS_TEST_FOR_EXCEPTION(
      imports.extent (0) < totalNumImportPackets, std::runtime_error,
      "Tpetra::Distributor::doPostsKokkos(4 args, Kokkos): The 'imports' array must have "
      "enough entries to hold the expected number of import packets.  "
      "imports.extent(0) = " << imports.extent (0) << " < "
      "totalNumImportPackets = " << totalNumImportPackets << ".");
  TEUCHOS_TEST_FOR_EXCEPTION
    (requests_.size () != 0, std::logic_error, "Tpetra::Distributor::"
     "doPostsKokkos(4 args, Kokkos): Process " << myProcID << ": requests_.size () = "
     << requests_.size () << " != 0.");
#endif // HAVE_TPETRA_DEBUG
  // Distributor uses requests_.size() as the number of outstanding
  // nonblocking message requests, so we resize to zero to maintain
  // this invariant.
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
  const size_type actualNumReceives = as<size_type> (plan.getNumReceives()) +
    as<size_type> (plan.hasSelfMessage() ? 1 : 0);
  requests_.resize (0);

  // Post the nonblocking receives.  It's common MPI wisdom to post
  // receives before sends.  In MPI terms, this means favoring
  // adding to the "posted queue" (of receive requests) over adding
  // to the "unexpected queue" (of arrived messages not yet matched
  // with a receive).
  {
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonRecvs (*timer_doPosts4KV_recvs_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    size_t curBufferOffset = 0;
    size_t curLIDoffset = 0;

    for (size_type i = 0; i < actualNumReceives; ++i) {
      size_t totalPacketsFrom_i = 0;
      Kokkos::parallel_reduce("totalPacketsFrom_i", Kokkos::RangePolicy<ImpExecSpace>(0, plan.getLengthsFrom()[i]), KOKKOS_LAMBDA(const size_t j, size_t& total) {
        total += numImportPacketsPerLID(curLIDoffset+j);
      }, totalPacketsFrom_i);
      // totalPacketsFrom_i is converted down to int, so make sure it can be represented
      TEUCHOS_TEST_FOR_EXCEPTION(totalPacketsFrom_i > size_t(INT_MAX),
                                 std::logic_error, "Tpetra::Distributor::doPostsKokkos(3 args, Kokkos): "
                                 "Recv count for receive " << i << " (" << totalPacketsFrom_i << ") is too large "
                                 "to be represented as int.");
      curLIDoffset += plan.getLengthsFrom()[i];
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
          subview_offset (imports, curBufferOffset, totalPacketsFrom_i);
        requests_.push_back (ireceive<int> (recvBuf, plan.getProcsFrom()[i],
              mpiTag_, *plan.getComm()));
      }
      else { // Receiving these packet(s) from myself
        selfReceiveOffset = curBufferOffset; // Remember the offset
      }
      curBufferOffset += totalPacketsFrom_i;
    }
  }

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  Teuchos::TimeMonitor timeMonSends (*timer_doPosts4KV_sends_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

  // setup views containing starting-offsets into exports for each send,
  // and num-packets-to-send for each send.
  Kokkos::View<size_t*, Kokkos::DefaultHostExecutionSpace> sendPacketOffsets("sendPacketOffsets", plan.getNumSends());
  Kokkos::View<size_t*, Kokkos::DefaultHostExecutionSpace> packetsPerSend("packetsPerSend", plan.getNumSends());
  auto sendPacketOffsets_d = Kokkos::create_mirror_view(ExpExecSpace(), sendPacketOffsets);
  auto packetsPerSend_d = Kokkos::create_mirror_view(ExpExecSpace(), packetsPerSend);

  auto starts = Kokkos::Compat::getKokkosViewDeepCopy<ExpExecSpace>(plan.getStartsTo());
  auto lengths = Kokkos::Compat::getKokkosViewDeepCopy<ExpExecSpace>(plan.getLengthsTo());

  Kokkos::parallel_scan("sendPacketOffsets", Kokkos::RangePolicy<ExpExecSpace>(0, plan.getNumSends()), KOKKOS_LAMBDA(const size_t pp, size_t& curPKToffset, bool final_pass) {
    if(final_pass) sendPacketOffsets_d(pp) = curPKToffset;
    size_t numPackets = 0;
    for(size_t j = starts(pp); j < starts(pp) + lengths(pp); j++) {
      numPackets += numExportPacketsPerLID(j);
    }
    if(final_pass) packetsPerSend_d(pp) = numPackets;
    curPKToffset += numPackets;
  });

  size_t maxNumPackets;
  Kokkos::parallel_reduce("maxNumPackets", Kokkos::RangePolicy<ExpExecSpace>(0, plan.getNumSends()), KOKKOS_LAMBDA(const size_t pp, size_t& max) {
    if(packetsPerSend_d(pp) > max) {
      max = packetsPerSend_d(pp);
    }
  }, Kokkos::Max<size_t>(maxNumPackets));

  // numPackets will be used as a message length, so make sure it can be represented as int
  TEUCHOS_TEST_FOR_EXCEPTION(maxNumPackets > size_t(INT_MAX),
                             std::logic_error, "Tpetra::Distributor::doPostsKokkos(4 args, Kokkos): "
                             "numPackets = " << maxNumPackets << " is too large "
                             "to be represented as int.");

  Kokkos::deep_copy(sendPacketOffsets, sendPacketOffsets_d);
  Kokkos::deep_copy(packetsPerSend, packetsPerSend_d);

  // setup scan through getProcsTo() list starting with higher numbered procs
  // (should help balance message traffic)
  size_t numBlocks = plan.getNumSends() + plan.hasSelfMessage();
  size_t procIndex = 0;
  while ((procIndex < numBlocks) && (plan.getProcsTo()[procIndex] < myProcID)) {
    ++procIndex;
  }
  if (procIndex == numBlocks) {
    procIndex = 0;
  }

  size_t selfNum = 0;
  size_t selfIndex = 0;
  if (plan.getIndicesTo().is_null()) {

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonSends2 (*timer_doPosts4KV_sends_fast_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    // Data are already blocked (laid out) by process, so we don't
    // need a separate send buffer (besides the exports array).
    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myProcID && packetsPerSend[p] > 0) {
        exports_view_type tmpSend =
          subview_offset(exports, sendPacketOffsets[p], packetsPerSend[p]);

        if (sendType == Details::DISTRIBUTOR_ISEND) {
          exports_view_type tmpSendBuf =
            subview_offset (exports, sendPacketOffsets[p], packetsPerSend[p]);
          requests_.push_back (isend<int> (tmpSendBuf, plan.getProcsTo()[p],
                mpiTag_, *plan.getComm()));
        }
        else { // DISTRIBUTOR_SEND
          send<int> (tmpSend,
              as<int> (tmpSend.size ()),
              plan.getProcsTo()[p], mpiTag_, *plan.getComm());
        }
      }
      else { // "Sending" the message to myself
        selfNum = p;
      }
    }

    if (plan.hasSelfMessage()) {
      deep_copy_offset(imports, exports, selfReceiveOffset,
          sendPacketOffsets[selfNum], packetsPerSend[selfNum]);
    }
  }
  else { // data are not blocked by proc, use send buffer

#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor timeMonSends2 (*timer_doPosts4KV_sends_slow_);
#endif // HAVE_TPETRA_DISTRIBUTOR_TIMINGS

    // FIXME (mfh 05 Mar 2013) This may be broken for Isend.
    typedef typename ExpView::non_const_value_type Packet;
    typedef typename ExpView::array_layout Layout;
    typedef typename ExpView::device_type Device;
    typedef typename ExpView::memory_traits Mem;

    // This buffer is long enough for only one message at a time.
    // Thus, we use DISTRIBUTOR_SEND always in this case, regardless
    // of sendType requested by user.
    // This code path formerly errored out with message:
    //     Tpetra::Distributor::doPostsKokkos(4-arg, Kokkos):
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

    Kokkos::View<Packet*,Layout,Device,Mem> sendArray ("sendArray",
                                                        maxNumPackets);

    Kokkos::View<size_t*, ExpExecSpace> indicesOffsets ("indicesOffsets", numExportPacketsPerLID.extent(0));
    size_t ioffset = 0;
    Kokkos::parallel_scan("indicesOffsets", Kokkos::RangePolicy<ExpExecSpace>(0, numExportPacketsPerLID.extent(0)), KOKKOS_LAMBDA(const size_t j, size_t& offset, bool is_final) {
      if(is_final) indicesOffsets(j) = offset;
      offset += numExportPacketsPerLID(j);
    }, ioffset);

    for (size_t i = 0; i < numBlocks; ++i) {
      size_t p = i + procIndex;
      if (p > (numBlocks - 1)) {
        p -= numBlocks;
      }

      if (plan.getProcsTo()[p] != myProcID) {
        size_t j = plan.getStartsTo()[p];
        size_t numPacketsTo_p = 0;
        //mirror in case execspaces are different
        auto sendArrayMirror = Kokkos::create_mirror_view_and_copy(ExpExecSpace(), sendArray);
        auto exportsMirror = Kokkos::create_mirror_view_and_copy(ExpExecSpace(), exports);
        Kokkos::parallel_scan("deep_copy_offset", Kokkos::RangePolicy<ExpExecSpace>(0, plan.getLengthsTo()[p]), KOKKOS_LAMBDA(const size_t k, size_t& offset, bool is_final) {
          if(is_final) {
            const size_t dst_end = offset + numExportPacketsPerLID(j + k);
            const size_t src_end = indicesOffsets(j + k) + numExportPacketsPerLID(j + k);
            auto dst_sub = Kokkos::subview(sendArrayMirror, Kokkos::make_pair(offset, dst_end));
            auto src_sub = Kokkos::subview(exportsMirror, Kokkos::make_pair(indicesOffsets(j + k), src_end));
            Kokkos::Experimental::local_deep_copy(dst_sub, src_sub);
          }
          offset += numExportPacketsPerLID(j + k);
        }, numPacketsTo_p);
        Kokkos::deep_copy(sendArray, sendArrayMirror);
        typename ExpView::execution_space().fence();

        if (numPacketsTo_p > 0) {
          ImpView tmpSend =
            subview_offset(sendArray, size_t(0), numPacketsTo_p);

          send<int> (tmpSend,
              as<int> (tmpSend.size ()),
              plan.getProcsTo()[p], mpiTag_, *plan.getComm());
        }
      }
      else { // "Sending" the message to myself
        selfNum = p;
        selfIndex = plan.getStartsTo()[p];
      }
    }

    if (plan.hasSelfMessage()) {
      //mirror in case execspaces are different
      auto importsMirror = Kokkos::create_mirror_view_and_copy(ExpExecSpace(), imports);
      auto exportsMirror = Kokkos::create_mirror_view_and_copy(ExpExecSpace(), exports);
      size_t temp;
      Kokkos::parallel_scan("deep_copy_offset2", Kokkos::RangePolicy<ExpExecSpace>(0, plan.getLengthsTo()[selfNum]), KOKKOS_LAMBDA(const size_t k, size_t& offset, bool is_final) {
        if(is_final) {
          const size_t dst_end = selfReceiveOffset + offset + numExportPacketsPerLID(selfIndex + k);
          const size_t src_end = indicesOffsets(selfIndex + k) + numExportPacketsPerLID(selfIndex + k);
          auto dst_sub = Kokkos::subview(importsMirror, Kokkos::make_pair(selfReceiveOffset + offset, dst_end));
          auto src_sub = Kokkos::subview(exportsMirror, Kokkos::make_pair(indicesOffsets(selfIndex + k), src_end));
          Kokkos::Experimental::local_deep_copy(dst_sub, src_sub);
        }
        offset += numExportPacketsPerLID(selfIndex + k);
      }, temp);
      Kokkos::deep_copy(imports, importsMirror);
      selfIndex += plan.getLengthsTo()[selfNum];
      selfReceiveOffset += temp;
    }
  }
}

}
}

#endif
