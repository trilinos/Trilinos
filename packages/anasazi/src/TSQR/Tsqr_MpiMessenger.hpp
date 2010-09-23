#ifndef __TSQR_MpiMessenger_hpp
#define __TSQR_MpiMessenger_hpp

#include <mpi.h>
#include <Tsqr_MessengerBase.hpp>
#include <Tsqr_MpiDatatype.hpp>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace MPI {

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
      size() const
      {
	int nprocs = 0;
	const int err = MPI_Comm_size (comm_, &nprocs);
	if (err != MPI_SUCCESS)
	  throw std::runtime_error ("MPI_Comm_size failed");
	else if (nprocs <= 0)
	  throw std::runtime_error ("MPI_Comm_size returned # processors <= 0");
	return nprocs;
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

      MpiDatatype< Datum > mpiType_;
    };
  } // namespace MPI
} // namespace TSQR

#endif // __TSQR_MpiMessenger_hpp
