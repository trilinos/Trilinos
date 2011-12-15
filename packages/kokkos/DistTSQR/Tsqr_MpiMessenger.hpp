/*
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
*/

#ifndef __TSQR_MpiMessenger_hpp
#define __TSQR_MpiMessenger_hpp

#include <Teuchos_ConfigDefs.hpp> // HAVE_MPI

#ifdef HAVE_MPI
#  include <Tsqr_MessengerBase.hpp>
#  include <Tsqr_MpiDatatype.hpp>
#  include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace MPI {

    /// \class MpiMessenger
    /// \brief Raw MPI-based communication object for TSQR
    ///
    /// A thin wrapper around an MPI_Comm, for use by TSQR.  I wrote
    /// this class specifically to experiment with the new MPI
    /// communicator objects, like MPI_COMM_NETWORK, MPI_COMM_NODE,
    /// and MPI_COMM_SOCKET.  These are not yet supported nicely in
    /// Teuchos.
    ///
    /// This class may also be used whenever it's preferred to work
    /// with a raw MPI_Comm object, rather than with the Teuchos comm
    /// wrappers.
    ///
    /// \warning Datum should be a class with value-type semantics.
    ///
    template< class Datum >
    class MpiMessenger : public TSQR::MessengerBase< Datum > {
    public:
      MpiMessenger (MPI_Comm comm) : comm_ (comm) {}
      virtual ~MpiMessenger () {}

      virtual void 
      send (const Datum sendData[], 
	    const int sendCount, 
	    const int destProc, 
	    const int tag)
      {
	const int err = 
	  MPI_Send (const_cast< Datum* const >(sendData), sendCount, 
		    mpiType_.get(), destProc, tag, comm_);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Send failed");
      }

      virtual void 
      recv (Datum recvData[], 
	    const int recvCount, 
	    const int srcProc, 
	    const int tag)
      {
	MPI_Status status;
	const int err = MPI_Recv (recvData, recvCount, mpiType_.get(), 
				  srcProc, tag, comm_, &status);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Recv failed");
      }

      virtual void 
      swapData (const Datum sendData[], 
		Datum recvData[], 
		const int sendRecvCount, 
		const int destProc, 
		const int tag)
      {
	MPI_Status status;
	const int err = 
	  MPI_Sendrecv (const_cast< Datum* const >(sendData), sendRecvCount, 
			mpiType_.get(), destProc, tag,
			recvData, sendRecvCount, mpiType_.get(), destProc, tag,
			comm_, &status);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Sendrecv failed");
      }

      virtual Datum 
      globalSum (const Datum& inDatum)
      {
	// Preserve const semantics of inDatum, by copying it and
	// using the copy in MPI_Allreduce().
	Datum input (inDatum);
	Datum output;

	int count = 1;
	const int err = MPI_Allreduce (&input, &output, count, 
				       mpiType_.get(), MPI_SUM, comm_);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Allreduce (MPI_SUM) failed");
	return output;
      }

      virtual Datum 
      globalMin (const Datum& inDatum)
      {
	// Preserve const semantics of inDatum, by copying it and
	// using the copy in MPI_Allreduce().
	Datum input (inDatum);
	Datum output;

	int count = 1;
	const int err = MPI_Allreduce (&input, &output, count, 
				       mpiType_.get(), MPI_MIN, comm_);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Allreduce (MPI_MIN) failed");
	return output;
      }

      virtual Datum 
      globalMax (const Datum& inDatum)
      {
	Datum input (inDatum);
	Datum output;

	int count = 1;
	const int err = MPI_Allreduce (&input, &output, count, 
				       mpiType_.get(), MPI_MAX, comm_);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Allreduce (MPI_MAX) failed");
	return output;
      }

      virtual void
      globalVectorSum (const Datum inData[], 
		       Datum outData[], 
		       const int count)
      {
	const int err = 
	  MPI_Allreduce (const_cast< Datum* const > (inData), outData, 
			 count, mpiType_.get(), MPI_SUM, comm_);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Allreduce failed");
      }

      virtual void
      broadcast (Datum data[], 
		 const int count,
		 const int root)
      {
	const int err = MPI_Bcast (data, count, mpiType_.get(), root, comm_);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Bcast failed");
      }

      virtual int 
      rank() const
      {
	int my_rank = 0;
	const int err = MPI_Comm_rank (comm_, &my_rank);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Comm_rank failed");
	return my_rank;
      }

      virtual int 
      size() const
      {
	int nprocs = 0;
	const int err = MPI_Comm_size (comm_, &nprocs);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Comm_size failed");
	else if (nprocs <= 0)
	  // We want to make sure that there is always at least one
	  // valid rank (at least rank() == 0).  The messenger can't
	  // do anything useful with MPI_COMM_NULL.
	  throw std::runtime_error ("MPI_Comm_size returned # processors <= 0");
	return nprocs;
      }

      virtual void 
      barrier () const 
      {
	const int err = MPI_Barrier (comm_);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Barrier failed");
      }

    private:
      /// \brief Handle to MPI communicator object
      ///
      /// Handle to the MPI communicator object used by this Messenger.
      /// This is a handle, so you can copy by value.  It's "mutable"
      /// because MPI functions don't take const arguments, even if
      /// they don't modify anything.  (C originally didn't have a
      /// notion of "const," and the MPI interface was designed for C.)
      mutable MPI_Comm comm_; 

      /// Wrapper around an MPI_Datatype corresponding to Datum.  The
      /// wrapper will handle allocation and deallocation of the
      /// MPI_Datatype object, if necessary.
      MpiDatatype< Datum > mpiType_;
    };
  } // namespace MPI
} // namespace TSQR

#endif // HAVE_MPI
#endif // __TSQR_MpiMessenger_hpp
