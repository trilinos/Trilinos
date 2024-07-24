// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file 

   How we communicate (Send, ISend).
   How / whether a Distributor is initialized.
   Data structures tracking where (ranks) messages go and come from.

   Lengths are not in terms of bytes, but some kind of abstract object.

   Reverse plan: if X -> Y, then Y -> X in the reverse plan, e.g. in a halo exchange, we both want to receive ghost elements and send our own.

*/

#ifndef TPETRA_DETAILS_DISTRIBUTOR_PLAN_HPP
#define TPETRA_DETAILS_DISTRIBUTOR_PLAN_HPP

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include "TpetraCore_config.h"

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
#include <mpi_advance.h>
#endif

namespace Tpetra {
namespace Details {

/// \brief The type of MPI send that Distributor should use.
///
/// This is an implementation detail of Distributor.  Please do
/// not rely on these values in your code.
enum EDistributorSendType {
  DISTRIBUTOR_ISEND, // Use MPI_Isend (Teuchos::isend)
  DISTRIBUTOR_SEND,  // Use MPI_Send (Teuchos::send)
  DISTRIBUTOR_ALLTOALL // Use MPI_Alltoall
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  ,
  DISTRIBUTOR_MPIADVANCE_ALLTOALL,
  DISTRIBUTOR_MPIADVANCE_NBRALLTOALLV
#endif
};

/// \brief Convert an EDistributorSendType enum value to a string.
///
/// This is an implementation detail of Distributor.  Please do
/// not rely on this function in your code.
std::string
DistributorSendTypeEnumToString (EDistributorSendType sendType);

/// \brief Enum indicating how and whether a Distributor was initialized.
///
/// This is an implementation detail of Distributor.  Please do
/// not rely on these values in your code.
enum EDistributorHowInitialized {
  DISTRIBUTOR_NOT_INITIALIZED, // Not initialized yet
  DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS, // By createFromSends
  DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_RECVS, // By createFromRecvs
  DISTRIBUTOR_INITIALIZED_BY_CREATE_FROM_SENDS_N_RECVS, // By createFromSendsAndRecvs
  DISTRIBUTOR_INITIALIZED_BY_REVERSE, // By createReverseDistributor
  DISTRIBUTOR_INITIALIZED_BY_COPY, // By copy constructor
};

/// \brief Convert an EDistributorHowInitialized enum value to a string.
///
/// This is an implementation detail of Distributor.  Please do
/// not rely on this function in your code.
std::string
DistributorHowInitializedEnumToString (EDistributorHowInitialized how);

/// Instances of DistributorPlan take the following parameters that
/// control communication and debug output:
/// - "Send type" (<tt>std::string</tt>): When using MPI, the
///   variant of MPI_Send to use in do[Reverse]Posts().  Valid
///   values include "Isend",
///    and "Send".  The
///   default is "Send".  (The receive type is always MPI_Irecv, a
///   nonblocking receive.  Since we post receives first before
///   sends, this prevents deadlock, even if MPI_Send blocks and
///   does not buffer.)
class DistributorPlan : public Teuchos::ParameterListAcceptorDefaultBase {
  static constexpr int DEFAULT_MPI_TAG = 0;

public:
  DistributorPlan(Teuchos::RCP<const Teuchos::Comm<int>> comm);
  DistributorPlan(const DistributorPlan& otherPlan);

  size_t createFromSends(const Teuchos::ArrayView<const int>& exportProcIDs);
  void createFromRecvs(const Teuchos::ArrayView<const int>& remoteProcIDs);
  void createFromSendsAndRecvs(const Teuchos::ArrayView<const int>& exportProcIDs,
                               const Teuchos::ArrayView<const int>& remoteProcIDs);

  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& plist);

  Teuchos::RCP<DistributorPlan> getReversePlan() const;

  Teuchos::RCP<const Teuchos::Comm<int>> getComm() const { return comm_; }
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  Teuchos::RCP<MPIX_Comm*> getMPIXComm() const { return mpixComm_; }
#endif
  EDistributorSendType getSendType() const { return sendType_; }
  size_t getNumReceives() const { return numReceives_; }
  size_t getNumSends() const { return numSendsToOtherProcs_; }
  bool hasSelfMessage() const { return sendMessageToSelf_; }
  size_t getMaxSendLength() const { return maxSendLength_; }
  size_t getTotalReceiveLength() const { return totalReceiveLength_; }
  Teuchos::ArrayView<const int> getProcsFrom() const { return procsFrom_; }
  Teuchos::ArrayView<const int> getProcsTo() const { return procIdsToSendTo_; }
  Teuchos::ArrayView<const size_t> getLengthsFrom() const { return lengthsFrom_; }
  Teuchos::ArrayView<const size_t> getLengthsTo() const { return lengthsTo_; }
  Teuchos::ArrayView<const size_t> getStartsTo() const { return startsTo_; }
  Teuchos::ArrayView<const size_t> getIndicesTo() const { return indicesTo_; }
  Details::EDistributorHowInitialized howInitialized() const { return howInitialized_; }

private:

  // after the plan has been created we have the info we need to initialize the MPI advance communicator
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  void initializeMpiAdvance();
#endif

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  void createReversePlan() const;

  /// \brief Compute receive info from sends.
  ///
  /// This method computes numReceives_, lengthsFrom_, procsFrom_,
  /// totalReceiveLength_, indicesFrom_, and startsFrom_.
  ///
  /// \note This method currently ignores the sendType_ 
  ///   parameter, and always uses ireceive() /
  ///   send() for communication of the process IDs from which our
  ///   process is receiving and their corresponding receive packet
  ///   counts.
  void computeReceives();

  Teuchos::RCP<const Teuchos::Comm<int>> comm_;
#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
  Teuchos::RCP<MPIX_Comm*> mpixComm_;
#endif

  Details::EDistributorHowInitialized howInitialized_;
  mutable Teuchos::RCP<DistributorPlan> reversePlan_;

  //! @name Parameters read in from the Teuchos::ParameterList
  //@{
  EDistributorSendType sendType_;
  //@}

  bool sendMessageToSelf_;
  size_t numSendsToOtherProcs_;
  Teuchos::Array<int> procIdsToSendTo_;

  /// \brief Starting index of the block of Packets to send to each process.
  ///
  /// Given an export buffer that contains all of the data being
  /// sent by this process, the block of Packets to send to process
  /// p will start at position startsTo_[p].
  ///
  /// This array has length numSendsToOtherProcs_ + sendMessageToSelf_ (that is, it
  /// includes the self message, if there is one).
  Teuchos::Array<size_t> startsTo_;

  /// \brief Length (in number of Packets) of my process' send to each process.
  ///
  /// lengthsTo_[p] is the length of my process' send to process p.
  /// This array has length numSendsToOtherProcs_ + sendMessageToSelf_ (that is, it
  /// includes the self message, if there is one).
  Teuchos::Array<size_t> lengthsTo_;

  /// \brief The maximum send length (in number of Packets) to another process.
  ///
  /// maxSendLength_ = max(lengthsTo_[p]) for p != my process rank.
  size_t maxSendLength_;

  /// \brief Offset (by message, not by number of Packets) into exports array.
  ///
  /// This array is used by both versions of doPosts().  In that
  /// method, <tt>indicesTo_[j]*numPackets</tt> is the offset into
  /// the <tt>exports</tt> array, where <tt>j = startsTo_[p]</tt>
  /// and p is an index iterating through the sends in reverse order
  /// (starting with the process rank right before the self message,
  /// if there is a self message, else the largest process rank to
  /// which this process sends).
  ///
  /// This array is only used if export data are not blocked (laid
  /// out) by process rank, that is, if we need to use a send
  /// buffer.  Otherwise, this array has no entries.  (In fact,
  /// Distributor currently uses this in both overloads of doPosts()
  /// to test whether data are laid out by process.)
  Teuchos::Array<size_t> indicesTo_;

  /// \brief The number of messages received by my process from other processes.
  ///
  /// This does <i>not</i> count self receives.  If sendMessageToSelf_ is
  /// true, the actual number of receives is one more (we assume
  /// that we only receive zero or one messages from ourself).
  ///
  /// This value is computed by the \c computeReceives() method.
  /// That method first includes self receives in the count, but at
  /// the end subtracts one if sendMessageToSelf_ is true.
  size_t numReceives_;

  /// \brief sum(lengthsFrom_)
  ///
  /// This is computed by \c createFromSends() and is used to
  /// allocate the receive buffer.  The reverse communicator's total
  /// receive length is the total send length of the forward
  /// communicator.
  size_t totalReceiveLength_;

  /// \brief Array of lengths of incoming messages.
  ///
  /// This array has length numReceives_ + sendMessageToSelf_.  Incoming
  /// message i from process procsFrom_[i] has length
  /// lengthsFrom_[i].
  Teuchos::Array<size_t> lengthsFrom_;

  /// \brief Array of ranks of the process from which the calling
  ///   process will receive a message.
  ///
  /// This array has length numReceives_ + sendMessageToSelf_.  Incoming
  /// message i was sent by process procsFrom_[i].
  Teuchos::Array<int> procsFrom_;

  /// \brief Array of offsets of incoming messages.
  ///
  /// This array has length numReceives_ + sendMessageToSelf_.  It is an
  /// exclusive prefix sum of lengthsFrom_.  It is only used for
  /// constructing the reverse Distributor.
  Teuchos::Array<size_t> startsFrom_;

  /// \brief List that becomes the reverse communicator's indicesTo_.
  ///
  /// Array of length totalReceiveLength_.  When creating the
  /// reverse Distributor, this is assigned to the reverse
  /// Distributor's indicesTo_.
  Teuchos::Array<size_t> indicesFrom_;
};

}
}

#endif
