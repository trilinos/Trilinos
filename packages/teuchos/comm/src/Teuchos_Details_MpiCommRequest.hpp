// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_DETAILS_MPICOMMREQUEST_HPP
#define TEUCHOS_DETAILS_MPICOMMREQUEST_HPP

#include "Teuchos_DefaultMpiComm.hpp"

namespace Teuchos {
namespace Details {

/// \class MpiCommRequest
/// \brief MPI implementation of CommRequest<int>.
///
/// This class wraps MPI_Request, which is MPI's reification of a
/// nonblocking communication operation.
///
/// Users would not normally create an instance of this class.  Calls
/// to nonblocking communication operations (such as ireceive() or
/// isend()) return a pointer to a CommRequest.  If the Comm is an
/// MpiComm, then the returned CommRequest is an MpiCommRequest.
///
/// Users might wish to create an MpiCommRequest directly if they want
/// to encapsulate an MPI_Request returned by an external library or
/// by their own code.
///
/// The releaseRawMpiRequest() method (inherited from the base class)
/// does <i>not</i> free (decrement the reference count of) the
/// buffer, because relinquishing ownership does not actually force
/// the nonblocking operation to complete.  The destructor will free
/// the buffer automatically.
class MpiCommRequest : public MpiCommRequestBase<int> {
public:
  /// \brief Constructor.
  ///
  /// \param rawMpiRequest [in/out] Raw MPI_Request created by the
  ///   nonblocking communication operation.
  ///
  /// \param buffer [in] Buffer used by the nonblocking communication
  ///   operation.  We keep this only to delay its deallocation until
  ///   after the communication operation completes.  This may be null
  ///   if the message length is zero or if rawMpiRequest is
  ///   <tt>MPI_REQUEST_NULL</tt>.
  ///
  /// We do not template the buffer type on the \c Packet type (the
  /// type of the data being communicated).  You must do an
  /// rcp_reinterpret_cast on input of the buffer.  We require this
  /// because this class does not need to access the buffer; it merely
  /// needs to retain a reference to prevent premature deallocation.
  MpiCommRequest (MPI_Request rawMpiRequest,
		  const ArrayRCP<const char>& buffer);

  //! Destructor; cancels the request if it is still pending.
  virtual ~MpiCommRequest ();

private:
  //! The buffer for the nonblocking communication operation.
  ArrayRCP<const char> buffer_;

  MpiCommRequest (); // Not defined
  MpiCommRequest (const MpiCommRequest&); // Not defined
  MpiCommRequest& operator= (const MpiCommRequest&); // Not defined
};

/// \fn mpiCommRequest
/// \brief Nonmember constructor for MpiCommRequest.
/// \relates MpiCommRequest
///
/// \param rawMpiRequest [in] The raw MPI_Request opaque object.
/// \param buffer [in] Buffer used by the nonblocking communication
///   operation.
RCP<MpiCommRequest>
mpiCommRequest (MPI_Request rawMpiRequest,
		const ArrayRCP<const char>& buffer);

} // namespace Details
} // namespace Teuchos

#endif // TEUCHOS_DETAILS_MPICOMMREQUEST_HPP
