// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Tsqr_MessengerBase_hpp
#define __TSQR_Tsqr_MessengerBase_hpp

#include "Tsqr_ConfigDefs.hpp"

namespace TSQR {
  /// \class MessengerBase
  ///
  /// Interface for an object that performs collective communication.
  /// Each message contains some number of objects of scalar type
  /// Datum.  Datum must have a default constructor and a copy
  /// constructor, and taking its address must make sense (in terms of
  /// extracting the useful data).
  template<class Datum>
  class MessengerBase {
  public:
    virtual ~MessengerBase() = default;

    /// Send sendData[0:sendCount-1] to process destProc.
    ///
    /// \param sendData [in] Array of value-type elements to send
    /// \param sendCount [in] Number of elements in the array
    /// \param destProc [in] Rank of destination process
    /// \param tag [in] MPI tag (ignored)
    virtual void
    send (const Datum sendData[],
          const int sendCount,
          const int destProc,
          const int tag) = 0;

    /// Receive recvData[0:recvCount-1] from process srcProc.
    ///
    /// \param recvData [out] Array of value-type elements to receive
    /// \param recvCount [in] Number of elements to receive in the array
    /// \param srcProc [in] Rank of sending process
    /// \param tag [in] MPI tag (ignored)
    virtual void
    recv (Datum recvData[],
          const int recvCount,
          const int srcProc,
          const int tag) = 0;

    /// \brief Exchange data between processors.
    ///
    /// Exchange sencRecvCount elements of sendData with processor
    /// destProc, receiving the result into recvData.  Assume that
    /// sendData and recvData do not alias one another.
    ///
    /// \param sendData [in] Array of value-type elements to send
    /// \param recvData [out] Array of value-type elements to
    ///   receive.  Caller is responsible for making sure that
    ///   recvData does not alias sendData.
    /// \param sendRecvCount [in] Number of elements to send and
    ///   receive in the array
    /// \param destProc [in] The "other" process' rank (to which
    ///   this process is sending data, and from which this process is
    ///   receiving data)
    /// \param tag [in] MPI tag (ignored)
    virtual void
    swapData (const Datum sendData[],
              Datum recvData[],
              const int sendRecvCount,
              const int destProc,
              const int tag) = 0;

    //! Sum inDatum on all processors, and return the result.
    virtual Datum
    globalSum (const Datum& inDatum) = 0;

    //! Sum inData[0:count-1] over all processors into outData.
    virtual void
    globalVectorSum (const Datum inData[],
                     Datum outData[],
                     const int count) = 0;

    /// \brief Compute the global minimum over all processors.
    ///
    /// Assumes that Datum objects are less-than comparable.
    virtual Datum
    globalMin (const Datum& inDatum) = 0;

    /// \brief Compute the global maximum over all processors.
    ///
    /// Assumes that Datum objects are less-than comparable.
    virtual Datum
    globalMax (const Datum& inDatum) = 0;

    //! Broadcast data[0:count-1] from root to all processors.
    virtual void
    broadcast (Datum data[],
               const int count,
               const int root) = 0;

    //! Return this process' rank.
    virtual int rank () const = 0;

    //! Return the total number of ranks in the communicator.
    virtual int size () const = 0;

    //! Execute a barrier over the communicator.
    virtual void barrier () const = 0;
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_MessengerBase_hpp

