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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_SRCDISTOBJECTWITHPACK_DECL_HPP
#define TPETRA_SRCDISTOBJECTWITHPACK_DECL_HPP

/// \file Tpetra_SrcDistObjectWithPack.hpp
/// \brief Abstract base class for sources of an Import or Export,
///   that also know how to pack themselves.

#include <Tpetra_SrcDistObject.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Tpetra {

  // Forward declaration of Distributor.  We don't actually need to
  // include Distributor here, since we only refer to it by reference
  // in this file and don't call any of its methods or refer to any of
  // its fields.  Subclasses should be sure to include
  // Tpetra_Distributor.hpp.
  class Distributor;

  /// \class SrcDistObjectWithPack
  /// \brief Abstract base class for objects that can be the source of
  ///   an Import or Export operation, and that also know how to pack
  ///   their data to send to the target object.
  ///
  /// If an object implements SrcDistObjectWithPack, then that object
  /// acknowledges that it knows how to pack its data as the source
  /// object of an Import or Export operation.  The target object in
  /// general assumes responsibility for packing the source object's
  /// data.  However, the target object (in its packAndPrepare method)
  /// may ask the source object to pack its own data, if the source
  /// object implements SrcDistObjectWithPack.
  template<class Packet, class LocalOrdinal>
  class SrcDistObjectWithPack : public SrcDistObject {
  public:
    /// \brief Pack the object's data for an Import or Export.
    ///
    /// \param exportLIDs [in] List of the entries (as local IDs in
    ///   the source object) we will be sending to other processes.
    ///
    /// \param exports [out] On exit, the buffer for data to send.
    ///
    /// \param numPacketsPerLID [out] On exit, the implementation of
    ///   this method must do one of two things: set
    ///   numPacketsPerLID[i] to contain the number of packets to be
    ///   exported for exportLIDs[i] and set constantNumPackets to
    ///   zero, or set constantNumPackets to a nonzero value.  If the
    ///   latter, the implementation need not fill numPacketsPerLID.
    ///
    /// \param constantNumPackets [out] On exit, 0 if numPacketsPerLID
    ///   has variable contents (different size for each LID).  If
    ///   nonzero, then it is expected that the number of packets per
    ///   LID is constant, and that constantNumPackets is that value.
    ///
    /// \param distor [in] The Distributor object we are using.
    virtual void 
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
	  Teuchos::Array<Packet>& exports,
	  const Teuchos::ArrayView<size_t>& numPacketsPerLID,
	  size_t& constantNumPackets,
	  Distributor &distor) const = 0;
  };

} // namespace Tpetra

#endif /* TPETRA_SRCDISTOBJECTWITHPACK_DECL_HPP */
