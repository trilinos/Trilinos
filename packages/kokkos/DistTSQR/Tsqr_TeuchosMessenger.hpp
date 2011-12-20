//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_TeuchosMessenger_hpp
#define __TSQR_TeuchosMessenger_hpp

#include <Teuchos_CommHelpers.hpp>
#include <Tsqr_MessengerBase.hpp>

#include <algorithm>


namespace TSQR { 

  /// \class TeuchosMessenger
  /// \brief Communication object for TSQR
  ///
  /// A thin wrapper around \c Teuchos::Comm<int>, for use by TSQR.
  /// The internode parallel part of TSQR communicates via a
  /// \c MessengerBase<Datum> interface.  \c TeuchosMessenger<Datum>
  /// implements that interface by wrapping \c Teuchos::Comm<int>.
  ///
  /// \tparam Datum A class with value-type semantics, whose instances
  ///   are less-than comparable.
  template<class Datum>
  class TeuchosMessenger : public MessengerBase<Datum> {
  public:
    typedef Teuchos::RCP<const Teuchos::Comm<int> > comm_ptr;

    //! Constructor, taking the communicator object to wrap.
    TeuchosMessenger (const comm_ptr& pComm) : pComm_ (pComm) {}

    //! Virtual destructor for memory safety of derived classes.
    virtual ~TeuchosMessenger() {}

    /// \brief Send sendData[0:sendCount-1] to process destProc.
    ///
    /// \param sendData [in] Array of value-type elements to send
    /// \param sendCount [in] Number of elements in the array
    /// \param destProc [in] Rank of destination process
    /// \param tag [in] MPI tag (ignored)
    void 
    send (const Datum sendData[], 
	  const int sendCount, 
	  const int destProc, 
	  const int tag) 
    {
      // NOTE (mfh 14 June 2010): Teuchos generates "tag" arguments to
      // MPI calls internally, so we ignore the tag here.  I don't use
      // tags for anything in TSQR, so it doesn't matter.
      Teuchos::send (*pComm_, sendCount, sendData, destProc);
    }

    /// \brief Receive recvData[0:recvCount-1] from process srcProc.
    ///
    /// \param recvData [out] Array of value-type elements to receive
    /// \param recvCount [in] Number of elements to receive in the array
    /// \param srcProc [in] Rank of sending process
    /// \param tag [in] MPI tag (ignored)
    void 
    recv (Datum recvData[], 
	  const int recvCount, 
	  const int srcProc, 
	  const int tag) 
    {
      // NOTE (mfh 14 June 2010): Teuchos generates "tag" arguments to
      // MPI calls internally, so we ignore the tag here.  I don't use
      // tags for anything in TSQR, so it doesn't matter.
      Teuchos::receive (*pComm_, srcProc, recvCount, recvData);
    }

    /// Exchange sencRecvCount elements of sendData with processor
    /// destProc, receiving the result into recvData.  Optimize for
    /// the case that sendData and recvData do not alias one another,
    /// which is the common case for TSQR.
    ///
    /// \warning This routine does only a rudimentary and incomplete
    ///   check for whether sendData and recvData alias one another.
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
    void 
    swapData (const Datum sendData[], 
	      Datum recvData[], 
	      const int sendRecvCount, 
	      const int destProc, 
	      const int tag)
    {
      if (destProc == rank())
	{
	  // If the sending and receiving processes are the same,
	  // then all we need to do is copy the data.  Hopefully in
	  // that case you aren't aliasing.  std::copy assumes that
	  // the third argument does not point to an element in the
	  // range of the first two arguments.
	  std::copy (sendData, sendData+sendRecvCount, recvData);
	}
      else
	{
	  using Teuchos::RCP;
	  using Teuchos::ArrayRCP;
	  using Teuchos::CommRequest;

	  const int srcProc = Teuchos::rank (*pComm_);

	  // If we can prove that sendData and recvData don't alias
	  // one another, use an isend and an ireceive to exchange
	  // them.  (Our test may not necessarily be safe in general,
	  // since we only check whether the pointers are equal and
	  // not whether the arrays overlap.  However, it is safe for
	  // the specific case of TSQR.)
	  //
	  // Otherwise, if the arrays do alias one another, safely
	  // perform a send and then a receive (or a receive and then
	  // a send, depending on whether this MPI process is the
	  // source or destination process).  
	  //
	  // (It would be nice if Teuchos had a sendRecv() routine, as
	  // of summer 2010 when this code was written.  As it stands,
	  // we have to do a send and then a receive.)
	  if (sendData == recvData)
	    {
	      // The smaller-rank process sends first, and the
	      // larger-rank process receives first.
	      //
	      // Teuchos::send() and Teuchos::recv() are blocking,
	      // so we may safely write to recvBuf even if it
	      // aliases sendBuf.
	      if (srcProc < destProc)
		{
		  Teuchos::send (*pComm_, sendRecvCount, sendData, destProc);
		  Teuchos::receive (*pComm_, destProc, sendRecvCount, recvData);
		}
	      else 
		{
		  Teuchos::receive (*pComm_, destProc, sendRecvCount, recvData);
		  Teuchos::send (*pComm_, sendRecvCount, sendData, destProc);
		}
	    }
	  else
	    {
	      ArrayRCP<const Datum> sendBuf (sendData, 0, sendRecvCount, false);
	      ArrayRCP<Datum> recvBuf (recvData, 0, sendRecvCount, false);

	      RCP<CommRequest> sendReq, recvReq;
	      if (srcProc < destProc)
		{
		  sendReq = Teuchos::isend (*pComm_, sendBuf, destProc);
		  recvReq = Teuchos::ireceive (*pComm_, recvBuf, destProc);
		}
	      else
		{
		  recvReq = Teuchos::ireceive (*pComm_, recvBuf, destProc);
		  sendReq = Teuchos::isend (*pComm_, sendBuf, destProc);
		}
	      // Wait on both the send and the receive to complete.  The
	      // two can happen independently, because sendBuf and recvBuf
	      // are different.  (We assert no aliasing of buffers here,
	      // and we've also checked above that destProc != rank().)
	      Teuchos::waitAll (*pComm_, Teuchos::tuple (sendReq, recvReq));
	    }
	}
    }

    //! Sum inDatum on all processors, and return the result.
    Datum 
    globalSum (const Datum& inDatum) 
    {
      Datum outDatum;
      Teuchos::reduceAll (*pComm_, Teuchos::REDUCE_SUM, inDatum, 
			  Teuchos::outArg(outDatum));
      return outDatum;
    }

    /// \brief Compute the global minimum over all processors.
    ///
    /// Assumes that Datum objects are less-than comparable.
    Datum 
    globalMin (const Datum& inDatum)
    {
      Datum outDatum;
      Teuchos::reduceAll (*pComm_, Teuchos::REDUCE_MIN, inDatum, 
			  Teuchos::outArg(outDatum));
      return outDatum;
    }

    /// \brief Compute the global maximum over all processors.
    ///
    /// Assumes that Datum objects are less-than comparable.
    Datum 
    globalMax (const Datum& inDatum)
    {
      Datum outDatum;
      Teuchos::reduceAll (*pComm_, Teuchos::REDUCE_MAX, inDatum, 
			  Teuchos::outArg(outDatum));
      return outDatum;
    }

    //! Sum inData[0:count-1] over all processors into outData.
    void
    globalVectorSum (const Datum inData[], 
		     Datum outData[], 
		     const int count) 
    {
      Teuchos::reduceAll (*pComm_, Teuchos::REDUCE_SUM, count, 
			  inData, outData);
    }

    //! Broadcast data[0:count-1] from root to all processors.
    void
    broadcast (Datum data[], 
	       const int count,
	       const int root)
    {
      Teuchos::broadcast (*pComm_, root, count, data);
    }

    //! Return this process' rank.
    int rank () const { return Teuchos::rank (*pComm_); }

    //! Return the total number of ranks in the communicator.
    int size () const { return Teuchos::size (*pComm_); }

    //! Execute a barrier over the communicator.
    void barrier () const { Teuchos::barrier (*pComm_); }

  private:

    //! Shared pointer to the the underlying communicator object.
    comm_ptr pComm_;
  };
} // namespace TSQR

#endif // __TSQR_TeuchosMessenger_hpp

