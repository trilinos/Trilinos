// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef __TSQR_Trilinos_TeuchosMessenger_hpp
#define __TSQR_Trilinos_TeuchosMessenger_hpp

#include <Teuchos_CommHelpers.hpp>
#include <Tsqr_MessengerBase.hpp>

#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR { 
  namespace Trilinos {
    /// \class TeuchosMessenger
    /// \brief Communication object for TSQR
    ///
    /// A thin wrapper around Teuchos::Comm, for use by TSQR.  The
    /// internode parallel part of TSQR communicates via a
    /// MessengerBase< Datum > interface.  TeuchosMessenger< Datum >
    /// implements that interface by wrapping Teuchos::Comm.
    ///
    /// \warning Datum should be a class with value-type semantics.
    template< class Datum >
    class TeuchosMessenger : public MessengerBase< Datum > {
    public:
      typedef Teuchos::RCP< const Teuchos::Comm<int> > comm_ptr;

      TeuchosMessenger (const comm_ptr& pComm) : pComm_ (pComm) {}

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
	    const int tag) 
      {
	/// \note (mfh 14 June 2010) Teuchos generates "tag" arguments to
	/// MPI calls internally, so we ignore the tag here.  I don't use
	/// tags for anything in TSQR, so it doesn't matter.
	Teuchos::send (*pComm_, sendCount, sendData, destProc);
      }

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
	    const int tag) 
      {
	/// \note (mfh 14 June 2010) Teuchos generates "tag" arguments to
	/// MPI calls internally, so we ignore the tag here.  I don't use
	/// tags for anything in TSQR, so it doesn't matter.
	Teuchos::receive (*pComm_, srcProc, recvCount, recvData);
      }

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

	    // FIXME (mfh 14 June 2010, 09 Jul 2010) It would be nice
	    // if Teuchos had a sendRecv() routine... as it is, we
	    // have to do a send and then a receive.  We could do an
	    // isend and an ireceive in order to exploit potential
	    // overlap of the two messages.  That works if sendData
	    // and recvData don't alias one another.  We only do a
	    // partial check for aliasing here (sendData == recvData).
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
		ArrayRCP< const Datum > sendBuf (sendData, 0, sendRecvCount, false);
		ArrayRCP< Datum > recvBuf (recvData, 0, sendRecvCount, false);

		RCP< CommRequest > sendReq, recvReq;
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

      /// Allreduce sum all processors' inDatum data, returning the
      /// result (on all processors).
      virtual Datum 
      globalSum (const Datum& inDatum) 
      {
	Datum outDatum;
	Teuchos::reduceAll (*pComm_, Teuchos::REDUCE_SUM, inDatum, 
			    Teuchos::outArg(outDatum));
	return outDatum;
      }

      ///
      /// Assumes that Datum objects are less-than comparable by the
      /// underlying communication protocol.
      ///
      virtual Datum 
      globalMin (const Datum& inDatum)
      {
	Datum outDatum;
	Teuchos::reduceAll (*pComm_, Teuchos::REDUCE_MIN, inDatum, 
			    Teuchos::outArg(outDatum));
	return outDatum;
      }

      ///
      /// Assumes that Datum objects are less-than comparable by the
      /// underlying communication protocol.
      ///
      virtual Datum 
      globalMax (const Datum& inDatum)
      {
	Datum outDatum;
	Teuchos::reduceAll (*pComm_, Teuchos::REDUCE_MAX, inDatum, 
			    Teuchos::outArg(outDatum));
	return outDatum;
      }

      /// Allreduce sum all processors' inData[0:count-1], storing the
      /// results (on all processors) in outData.
      virtual void
      globalVectorSum (const Datum inData[], 
		       Datum outData[], 
		       const int count) 
      {
	Teuchos::reduceAll (*pComm_, Teuchos::REDUCE_SUM, count, 
			    inData, outData);
      }

      virtual void
      broadcast (Datum data[], 
		 const int count,
		 const int root)
      {
	// Assumes that Datum has value semantics.
	Teuchos::broadcast (*pComm_, root, count, data);
      }

      /// 
      /// Return this process' rank.
      /// 
      virtual int rank () const { return Teuchos::rank (*pComm_); }

      /// 
      /// Return the total number of ranks in the Teuchos::Comm communicator.
      /// 
      virtual int size () const { return Teuchos::size (*pComm_); }

      /// 
      /// Execute a barrier over the communicator.
      /// 
      virtual void barrier () const { Teuchos::barrier (*pComm_); }

    private:
      /// 
      /// Shared pointer to the the underlying Teuchos::Comm object.
      ///
      comm_ptr pComm_;
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TeuchosMessenger_hpp

