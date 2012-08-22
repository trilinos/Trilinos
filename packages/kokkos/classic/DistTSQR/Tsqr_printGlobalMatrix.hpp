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

#ifndef __Tsqr_printGlobalMatrix_hpp
#define __Tsqr_printGlobalMatrix_hpp

#include <Tsqr_MessengerBase.hpp>
#include <Tsqr_Util.hpp>
#include <Tsqr_Matrix.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <limits>
#include <ostream>
#include <stdexcept>


namespace TSQR {

  /// \fn printGlobalMatrix
  ///
  /// Print a dense matrix distributed in block row fashion among all
  /// MPI processes in a participating communicator.  The given
  /// "MessengerBase" communicator wrapper objects should wrap the
  /// same underlying communicator.
  ///
  /// \param out [out] Output stream to which to write the matrix (on
  ///   MPI Proc 0 only, relative to the underlying communicator).
  /// \param A_local [in] Each MPI process' part of the matrix.
  /// \param scalarComm [in/out] Communicator wrapper for
  ///   ConstMatrixViewType::scalar_type objects.
  /// \param ordinalComm [in/out] Communicator wrapper for
  ///   ConstMatrixViewType::ordinal_type objects.
  template<class ConstMatrixViewType>
  void
  printGlobalMatrix (std::ostream& out,
		     const ConstMatrixViewType& A_local,
		     MessengerBase<typename ConstMatrixViewType::scalar_type>* const scalarComm,
		     MessengerBase<typename ConstMatrixViewType::ordinal_type>* const ordinalComm)
    {
      typedef typename ConstMatrixViewType::ordinal_type LocalOrdinal;
      typedef typename ConstMatrixViewType::scalar_type Scalar;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      using std::endl;

      const int myRank = scalarComm->rank ();
      const int nprocs = scalarComm->size ();
      const LocalOrdinal nrowsLocal = A_local.nrows();
      const LocalOrdinal ncols = A_local.ncols();
      const Scalar quiet_NaN = STS::nan();

      if (myRank == 0)
	{
	  // Print the remote matrix data
	  // out << "Processor " << my_rank << ":" << endl;
	  print_local_matrix (out, A_local.nrows(), A_local.ncols(), 
			      A_local.get(), A_local.lda());

	  // Space for remote matrix data.  Other processors are allowed
	  // to have different nrows_local values; we make space as
	  // necessary.
	  Matrix<LocalOrdinal, Scalar> A_remote (nrowsLocal, ncols, quiet_NaN);

	  // Loop through all the other processors in order.
	  // Fetch their matrix data and print it.
	  for (int srcProc = 1; srcProc < nprocs; ++srcProc)
	    {
	      // Get processor proc's local matrix dimensions
	      LocalOrdinal dims[2];
	      ordinalComm->recv (&dims[0], 2, srcProc, 0);

	      // Make space for the remote matrix data.
	      //
	      // mfh 13 Oct 2010: Teuchos::OrdinalTraits does not
	      // currently have this feature.  It's OK to use
	      // std::numeric_limits, since ordinal types in Trilinos
	      // are intended to be built-in types (like int or long
	      // long int).  std::numeric_limits only promises to work
	      // for built-in types, unless someone has defined an
	      // appropriate specialization.  Teuchos::ScalarTraits,
	      // in contrast, has to work for non-built-in Scalar
	      // types, like ARPREC or QD floating-point numbers.
	      if (std::numeric_limits<LocalOrdinal>::is_signed)
		{
		  if (dims[0] <= 0 || dims[1] <= 0)
		    throw std::runtime_error ("Invalid dimensions of remote matrix");
		}
	      else
		{
		  if (dims[0] == 0 || dims[1] == 0)
		    throw std::runtime_error ("Invalid dimensions of remote matrix");
		}
	      A_remote.reshape (dims[0], dims[1]);

	      // Receive the remote matrix data, which we assume is
	      // stored contiguously.
	      scalarComm->recv (A_remote.get(), dims[0]*dims[1], srcProc, 0);

	      // Print the remote matrix data
	      // out << "Processor " << proc << ":" << endl;
	      print_local_matrix (out, dims[0], dims[0], A_remote.get(), A_remote.lda());
	    }
	}
      else
	{
	  // Send my local matrix dimensions to proc 0.
	  int rootProc = 0;
	  LocalOrdinal dims[2];

	  dims[0] = nrowsLocal;
	  dims[1] = ncols;
	  ordinalComm->send (dims, 2, rootProc, 0);

	  // Create a (contiguous) buffer and copy the data into it.
	  Matrix< LocalOrdinal, Scalar > A_buf (nrowsLocal, ncols, quiet_NaN);
	  A_buf.copy (A_local);

	  // Send the actual data to proc 0.
	  scalarComm->send (A_buf.get(), nrowsLocal*ncols, rootProc, 0);
	}
      scalarComm->barrier ();
    }

} // namespace TSQR

#endif // __Tsqr_printGlobalMatrix_hpp
