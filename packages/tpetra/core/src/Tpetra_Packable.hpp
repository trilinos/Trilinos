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

#ifndef TPETRA_PACKABLE_DECL_HPP
#define TPETRA_PACKABLE_DECL_HPP

/// \file Tpetra_Packable.hpp
/// \brief Abstract base class for sources of an Import or Export,
///   that also know how to pack themselves.

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Tpetra {

  // Forward declaration of Distributor.  We don't actually need to
  // include Distributor here, since we only refer to it by reference
  // in this file and don't call any of its methods or refer to any of
  // its fields.  Subclasses should be sure to include
  // Tpetra_Distributor.hpp.
  class Distributor;

  /// \class Packable
  /// \brief Abstract base class for objects that can be the source of
  ///   an Import or Export operation, and that also know how to pack
  ///   their data to send to the target object.
  /// \tparam Packet The type of each entry of the array of packed
  ///   data to be sent in the Import or Export operation.  The type
  ///   of packed data may differ from the type of actual data stored
  ///   in the object.  For example, a sparse matrix might need to
  ///   pack both column indices (an integer type) and values
  ///   (typically, but not always, a floating-point type), and it
  ///   might choose any of various types to represent packed data.
  /// \tparam LocalOrdinal The type of local indices in the object.
  ///   This is a template parameter because the pack() method
  ///   includes as input a list of local indices to pack.  See the
  ///   documentation of Map for requirements.
  ///
  /// If an object implements Packable, then that object acknowledges
  /// that it knows how to pack its data as the source object of an
  /// Import or Export operation.  The target object in general
  /// assumes responsibility for packing the source object's data.
  /// However, the target object (in its packAndPrepare method) may
  /// ask the source object to pack its own data, if the source object
  /// implements Packable.
  ///
  /// It might make sense for Packable to inherit from SrcDistObject.
  /// However, that sets up the possibility of ambiguous multiple
  /// inheritance.  For example, RowGraph inherits from Packable, and
  /// CrsGraph inherits from both RowGraph and DistObject.
  /// Furthermore, it is not necessary for a source object of an
  /// Import or Export to know how to pack itself.  The ability to
  /// pack oneself is independent of the ability to be the source of
  /// an Import or Export.  Packable exists mainly for syntactic
  /// enforcement of the interface needed for an object to know how to
  /// pack itself for an Import or Export.
  template<class Packet, class LocalOrdinal>
  class Packable {
  public:
    /// \brief Pack the object's data for an Import or Export.
    ///
    /// \param exportLIDs [in] List of the local indices of the
    ///   entries which the source object will send out.
    ///
    /// \param exports [out] On exit, the buffer packed with data to
    ///   send.  This object may resize the array if necessary.
    ///
    /// \param numPacketsPerLID [out] On exit, the implementation of
    ///   this method must do one of two things: set
    ///   numPacketsPerLID[i] to contain the number of packets to be
    ///   exported for exportLIDs[i] and set constantNumPackets to
    ///   zero, or set constantNumPackets to a nonzero value.  If the
    ///   latter, the implementation need not fill numPacketsPerLID.
    ///
    /// \param constantNumPackets [out] On exit, zero if
    ///   <tt>numPacketsPerLID</tt> has variable contents (different
    ///   size for each local index).  If nonzero, then it is expected
    ///   that the number of packets per local index is constant, and
    ///   that <tt>constantNumPackets</tt> is that value.
    ///
    /// \param distor [in] The Distributor object we are using.
    ///   Implementations may ignore this object.  We provide it for
    ///   consistency with DistObject's packAndPrepare method.
    virtual void
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
          Teuchos::Array<Packet>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets,
          Distributor &distor) const = 0;

    //! Destructor (declared virtual for memory safety of derived classes).
    virtual ~Packable () {}
  };

} // namespace Tpetra

#endif /* TPETRA_PACKABLE_DECL_HPP */
