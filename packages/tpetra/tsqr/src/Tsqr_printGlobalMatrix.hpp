// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Tsqr_printGlobalMatrix_hpp
#define __Tsqr_printGlobalMatrix_hpp

#include "Tsqr_MessengerBase.hpp"
#include "Tsqr_Util.hpp"
#include "Tsqr_Matrix.hpp"
#include "Teuchos_ScalarTraits.hpp"
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
  /// \param scalarComm [in/out] Communicator wrapper for scalar
  ///   objects.
  /// \param ordinalComm [in/out] Communicator wrapper for
  ///   ConstMatrixViewType::ordinal_type objects.
  template<class ConstMatrixViewType>
  void
  printGlobalMatrix (std::ostream& out,
                     const ConstMatrixViewType& A_local,
                     MessengerBase<typename ConstMatrixViewType::non_const_value_type>* const scalarComm,
                     MessengerBase<typename ConstMatrixViewType::ordinal_type>* const ordinalComm)
  {
    using LocalOrdinal = typename ConstMatrixViewType::ordinal_type;
    using Scalar = typename ConstMatrixViewType::non_const_value_type;
    using STS = Teuchos::ScalarTraits<Scalar>;
    using std::endl;

    const int myRank = scalarComm->rank ();
    const int nprocs = scalarComm->size ();
    const LocalOrdinal nrowsLocal = A_local.extent(0);
    const LocalOrdinal ncols = A_local.extent(1);
    const Scalar quiet_NaN = STS::nan();

    if (myRank == 0) {
      // Print the remote matrix data
      print_local_matrix (out, A_local.extent(0), A_local.extent(1),
                          A_local.data(), A_local.stride(1));

      // Space for remote matrix data.  Other processes are allowed to
      // have different nrows_local values; we make space as needed.
      Matrix<LocalOrdinal, Scalar> A_remote (nrowsLocal, ncols, quiet_NaN);

      // Loop through all the other processes in order.  Fetch their
      // matrix data and print it.
      for (int srcProc = 1; srcProc < nprocs; ++srcProc) {
        // Get local matrix dimensions
        LocalOrdinal dims[2];
        ordinalComm->recv (&dims[0], 2, srcProc, 0);

        // Make space for the remote matrix data.
        //
        // mfh 13 Oct 2010: Teuchos::OrdinalTraits does not currently
        // have this feature.  It's OK to use std::numeric_limits,
        // since ordinal types in Trilinos are intended to be built-in
        // types (like int or long long int).  std::numeric_limits
        // only promises to work for built-in types, unless someone
        // has defined an appropriate specialization.
        // Teuchos::ScalarTraits, in contrast, has to work for
        // non-built-in Scalar types, like ARPREC or QD floating-point
        // numbers.
        if (std::numeric_limits<LocalOrdinal>::is_signed) {
          if (dims[0] <= 0 || dims[1] <= 0) {
            throw std::runtime_error ("Invalid dimensions of remote matrix");
          }
        }
        else {
          if (dims[0] == 0 || dims[1] == 0) {
            throw std::runtime_error ("Invalid dimensions of remote matrix");
          }
        }
        A_remote.reshape (dims[0], dims[1]);

        // Receive the remote matrix data, which we assume is
        // stored contiguously.
        scalarComm->recv (A_remote.data(), dims[0]*dims[1], srcProc, 0);

        // Print the remote matrix data
        // out << "Processor " << proc << ":" << endl;
        print_local_matrix (out, dims[0], dims[0], A_remote.data(),
                            A_remote.stride(1));
      }
    }
    else {
      // Send my local matrix dimensions to proc 0.
      int rootProc = 0;
      LocalOrdinal dims[2];

      dims[0] = nrowsLocal;
      dims[1] = ncols;
      ordinalComm->send (dims, 2, rootProc, 0);

      // Create a (contiguous) buffer and copy the data into it.
      Matrix< LocalOrdinal, Scalar > A_buf (nrowsLocal, ncols, quiet_NaN);
      deep_copy (A_buf, A_local);

      // Send the actual data to proc 0.
      scalarComm->send (A_buf.data(), nrowsLocal*ncols, rootProc, 0);
    }
    scalarComm->barrier ();
  }
} // namespace TSQR

#endif // __Tsqr_printGlobalMatrix_hpp
