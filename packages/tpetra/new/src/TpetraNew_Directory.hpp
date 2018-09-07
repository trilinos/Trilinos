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

#ifndef TPETRANEW_DIRECTORY_HPP
#define TPETRANEW_DIRECTORY_HPP

/// \file TpetraNew_Directory.hpp
/// \brief Declaration of TpetraNew::Directory

#include "TpetraCore_config.h"
#include "TpetraNew_Map.hpp"
#include "Tpetra_TieBreak.hpp"
#include "Teuchos_Describable.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace TpetraNew {
namespace Details {  
  class Directory; // forward declaration
} // namespace Details
} // namespace TpetraNew
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace TpetraNew {
  /// \class Directory
  /// \brief Implement mapping from global ID to process ID and local ID.
  ///
  /// This class is an implementation detail of Map.  It is mainly of
  /// interest to Tpetra developers and does not normally appear in
  /// users' code.
  ///
  /// Directory implements looking up the process ranks and local
  /// indices (IDs) corresponding to a given list of global indices
  /// (IDs).  Each Map owns a Directory object that does this.
  /// Map::getRemoteIndexList() calls the Map's directory's
  /// getDirectoryEntries() method directly.  Directory has four
  /// different ways to perform this lookup, depending on the kind of
  /// Map.
  ///
  /// 1. If the user's Map is not distributed (i.e., is serial or
  ///    locally replicated), then my process ID is the process ID for
  ///    all global IDs.  The Directory gets the local ID (if
  ///    requested) directly from the user's Map via its
  ///    getLocalIndex() method (which requires no communication).
  ///
  /// 2. If the user's Map is distributed, contiguous, and uniform,
  ///    the Directory computes a global ID's process ID and local ID
  ///    on that process using a simple mathematical formula.  This
  ///    requires \f$O(1)\f$ arithmetic operations per global ID, and
  ///    no additional storage (beyond what the Map already stores).
  ///
  /// 3. If the user's Map is distributed, contiguous, but not
  ///    uniform, then the Directory builds an array (replicated on
  ///    all processes) that maps from each process ID to the minimum
  ///    global ID that it owns.  This trades time and communication
  ///    for space (P+1 entries in the array if there are P processes
  ///    in the communicator), in that it allows lookups without
  ///    communication (once the array has been built).
  ///
  /// 4. If the user's Map is distributed and noncontiguous, then the
  ///    Directory must store a general mapping from global ID to
  ///    (process ID, local ID).  It can't afford to store the whole
  ///    mapping redundantly on all processes, so the Directory
  ///    distributes it using another Map (the "directory Map").  This
  ///    is a contiguous uniform Map whose keys are the global IDs.
  ///
  /// \note To Tpetra developers: As of fall 2013, Directory no longer
  ///   keeps a reference to the Map that created it.  This will
  ///   facilitate Tpetra's port to use (new) Kokkos data structures
  ///   and handle semantics, by removing this circular dependency.
  ///
  /// \note To Epetra developers: This class corresponds roughly to
  ///   Epetra_Directory or Epetra_BasicDirectory.  Epetra_BlockMap
  ///   creates its Epetra_Directory object on demand whenever the
  ///   map's RemoteIDList() method is called.  Map's
  ///   getRemoteIndexList() method assumes that the Map's directory
  ///   already exists.  Epetra_Directory is an abstract interface
  ///   with one implementation (Epetra_BasicDirectory); Directory is
  ///   a concrete implementation.
  class Directory : public Teuchos::Describable {
  public:
    //! Type of the Map specialization to give to the constructor.
    using map_type = Map;
    using local_ordinal_type = map_type::local_ordinal_type;
    using global_ordinal_type = map_type::global_ordinal_type;

    //! @name Constructors/Destructor.
    //@{

    /// \brief Default constructor: the only one you should use.
    ///
    /// To initialize the Directory, call the appropriate initialize()
    /// overload.
    Directory ();

    //! Destructor.
    ~Directory ();

    //! Initialize the Directory with its Map.
    void initialize (const map_type& map);

    //! Initialize the Directory, with its Map and a TieBreak object.
    void
    initialize (const map_type& map,
                const ::Tpetra::Details::TieBreak<local_ordinal_type, global_ordinal_type>& tieBreak);

    //! Whether the Directory is initialized.
    bool initialized () const;

    //@}
    //! @name Implementation of Teuchos::Describable.
    //@{

    //! A one-line human-readable description of this object.
    std::string description () const;

    //@}
    //! @name Query methods.
    //@{

    /// \brief Given a global ID list, return the list of their owning process IDs.
    ///
    /// Given a list \c globalIDs of global identifiers (GIDs), return
    /// the corresponding list \c processRanks of the process ranks which
    /// own those GIDs.  Tpetra uses this to figure out the locations
    /// of nonlocal Map entries.
    ///
    /// \param map [in] The Map given to this Directory instance's
    ///   constructor.
    ///
    /// \param globalIDs [in] List of global IDs to look up.
    ///
    /// \param processRanks [out] On input, an array view with the same
    ///   number of entries as \c globalIDs.  On output, processRanks[i] is
    ///   the ID of the process which owns globalIDs[i].  If
    ///   globalIDs[i] is not present in the directory (i.e., is not
    ///   owned by any process in the Directory's communicator),
    ///   processRanks[i] is -1.
    ///
    /// \return If at least one global ID was not present in the
    ///   directory, return IDNotPresent.  Otherwise, return
    ///   AllIDsPresent.
    ///
    /// \note If <tt>processRanks.size() != globalIDs.size()</tt>, then
    ///   this method throws std::runtime_error.
    ::Tpetra::LookupStatus
    getDirectoryEntries (const map_type& map,
                         const Teuchos::ArrayView<const global_ordinal_type>& globalIDs,
                         const Teuchos::ArrayView<int>& processRanks) const;

    /// \brief Given a global ID list, return a list of their owning
    ///   process IDs and their corresponding local IDs.
    ///
    /// Given a list \c globalIDs of global identifiers (GIDs), return
    /// the corresponding list \c processRanks of the process ranks which
    /// own those GIDs, as well as the list of the local identifiers
    /// (LIDs).  Tpetra uses this to figure out the locations of
    /// nonlocal Map entries.
    ///
    /// \param map [in] The Map given to this Directory instance's
    ///   constructor.
    ///
    /// \param globalIDs [in] List of global IDs to look up.
    ///
    /// \param processRanks [out] On input, an array view with the same
    ///   number of entries as \c globalIDs.  On output, processRanks[i] is
    ///   the ID of the process which owns globalIDs[i].  If
    ///   globalIDs[i] is not present in the directory (i.e., is not
    ///   owned by any process in the Directory's communicator), then
    ///   <tt>processRanks[i] == -1</tt>.
    ///
    /// \param localIDs [out] On input, an array view with the same
    ///   number of entries as \c globalIDs.  On output, \c
    ///   localIDs[i] is the local identifier corresponding to the
    ///   global identifier \c globalIDs[i].  If globalIDs[i] is not
    ///   present in the directory, then <tt>localIDs[i] ==
    ///   Teuchos::OrdinalTraits<local_ordinal_type>::invalid()</tt>.
    ///
    /// \return If at least one global ID was not present in the
    ///   directory, return IDNotPresent.  Otherwise, return
    ///   AllIDsPresent.
    ///
    /// \note If <tt>processRanks.size() != globalIDs.size()</tt> or
    ///   <tt>localIDs.size() != globalIDs.size()</tt>, then
    ///   this method throws std::runtime_error.
    ::Tpetra::LookupStatus
    getDirectoryEntries (const map_type& map,
                         const Teuchos::ArrayView<const global_ordinal_type>& globalIDs,
                         const Teuchos::ArrayView<int>& processRanks,
                         const Teuchos::ArrayView<local_ordinal_type>& localIDs) const;

    /// \brief Whether the Directory's input Map is (globally) one to one.
    ///
    /// This method should always be treated as a collective on all
    /// processes in the given Map's communicator, which must be the
    /// same as the input Map's communicator.  Not all implementations
    /// necessarily communicate.
    bool isOneToOne (const map_type& map) const;

    //@}
  private:
    /// \brief Type of the (base class) implementation of this object.
    ///
    /// \note To implementers: Directory has different
    ///   implementations, depending on characteristics of the input
    ///   Map (e.g., locally replicated or globally distributed,
    ///   contiguous or noncontiguous).
    using base_type = Details::Directory;

    /// \brief Implementation of this object.
    ///
    /// Directory creates this impl_ object either lazily, or on
    /// request (in initialize()), and caches it.
    const base_type* impl_ = nullptr;

    Directory (const Directory& directory) = delete;

    Directory&
    operator= (const Directory& source) = delete;
  }; // class Directory

} // namespace TpetraNew

#endif // TPETRANEW_DIRECTORY_HPP

