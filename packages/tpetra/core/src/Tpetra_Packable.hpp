// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_PACKABLE_DECL_HPP
#define TPETRA_PACKABLE_DECL_HPP

/// \file Tpetra_Packable.hpp
/// \brief Declaration of Tpetra::Packable

#include "Tpetra_Packable_fwd.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // Forward declarations
  template <class T> class Array;
  template <class T> class ArrayView;
} // namespace Teuchos

#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {
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
  template<class Packet,
           class LocalOrdinal>
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
    virtual void
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
          Teuchos::Array<Packet>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets) const = 0;

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~Packable () {}
  };
} // namespace Tpetra

#endif /* TPETRA_PACKABLE_DECL_HPP */
