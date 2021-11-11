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

#ifndef TPETRA_DETAILS_DISTRIBUTOR_PLAN_HPP
#define TPETRA_DETAILS_DISTRIBUTOR_PLAN_HPP

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include "TpetraCore_config.h"

namespace Tpetra {
namespace Details {

namespace {
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  const bool barrierBetween_default = false;
#endif
  const bool useDistinctTags_default = true;
}

/// \brief The type of MPI send that Distributor should use.
///
/// This is an implementation detail of Distributor.  Please do
/// not rely on these values in your code.
enum EDistributorSendType {
  DISTRIBUTOR_ISEND, // Use MPI_Isend (Teuchos::isend)
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  DISTRIBUTOR_RSEND, // Use MPI_Rsend (Teuchos::readySend)
#endif
  DISTRIBUTOR_SEND   // Use MPI_Send (Teuchos::send)
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  , DISTRIBUTOR_SSEND  // Use MPI_Ssend (Teuchos::ssend)
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

/// Instances of Distributor take the following parameters that
/// control communication and debug output:
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
/// - "Barrier between receives and sends" (<tt>bool</tt>): (DEPRECATED)
///   Whether to execute a barrier between receives and sends in
///   do[Reverse]Posts().  
///   A barrier is required for correctness
///   when the "Send type" parameter is "Rsend".  Otherwise,
///   A barrier is correct and may be useful for debugging, but not
///   recommended, since it introduces useless synchronization.
#endif
/// - "Send type" (<tt>std::string</tt>): When using MPI, the
///   variant of MPI_Send to use in do[Reverse]Posts().  Valid
///   values include "Isend",
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
///   "Rsend (DEPRECATED)", "Ssend (DEPRECATED)",
#endif
///    and "Send".  The
///   default is "Send".  (The receive type is always MPI_Irecv, a
///   nonblocking receive.  Since we post receives first before
///   sends, this prevents deadlock, even if MPI_Send blocks and
///   does not buffer.)
class DistributorPlan : public Teuchos::ParameterListAcceptorDefaultBase {
public:
  DistributorPlan(Teuchos::RCP<const Teuchos::Comm<int>> comm);
  DistributorPlan(const DistributorPlan& otherPlan);

  //! Get the tag to use for receives and sends.
  ///
  /// See useDistinctTags_.  This is called in doPosts() (both
  /// variants) and computeReceives().
  int getTag(const int pathTag) const;

  size_t createFromSends(const Teuchos::ArrayView<const int>& exportProcIDs);
  void createFromRecvs(const Teuchos::ArrayView<const int>& remoteProcIDs);
  void createFromSendsAndRecvs(const Teuchos::ArrayView<const int>& exportProcIDs,
                               const Teuchos::ArrayView<const int>& remoteProcIDs);

  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& plist);

  Teuchos::RCP<DistributorPlan> getReversePlan() const;

  Teuchos::RCP<const Teuchos::Comm<int>> getComm() const { return comm_; }
  EDistributorSendType getSendType() const { return sendType_; }
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  bool barrierBetweenRecvSend() const { return barrierBetweenRecvSend_; }
#endif
  bool useDistinctTags() const { return useDistinctTags_; }
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
  Details::EDistributorHowInitialized howInitialized_;
  mutable Teuchos::RCP<DistributorPlan> reversePlan_;

  //! @name Parameters read in from the Teuchos::ParameterList
  //@{
  EDistributorSendType sendType_;
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  bool barrierBetweenRecvSend_;
#endif

  /// \brief Whether to use different tags for different code paths.
  ///
  /// There are currently three code paths in Distributor that post
  /// receives and sends:
  ///
  /// 1. Three-argument variant of doPosts()
  /// 2. Four-argument variant of doPosts()
  /// 3. computeReceives()
  ///
  /// If this option is true, Distributor will use a distinct
  /// message tag for each of these paths.
  bool useDistinctTags_;
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
