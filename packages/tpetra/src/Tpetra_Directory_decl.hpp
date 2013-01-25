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

#ifndef TPETRA_DIRECTORY_DECL_HPP
#define TPETRA_DIRECTORY_DECL_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_DirectoryImpl_decl.hpp"

namespace Tpetra {

  /// \class Directory
  /// \brief Implement mapping from global ID to process ID and local ID.
  ///
  /// This class is an implementation detail of \c Map.  It is mainly
  /// of interest to Tpetra developers and does not normally appear in
  /// users' code.  If using this with a Map, the template parameters
  /// of Directory should always be the same as the template
  /// parameters of Map.
  ///
  /// \tparam LocalOrdinal Same as Map's \c LocalOrdinal template
  ///   parameter.  The type of local IDs.
  ///
  /// \tparam GlobalOrdinal Same as Map's \c GlobalOrdinal template
  ///   parameter.  The type of global IDs.  Defaults to the same type
  ///   as LocalOrdinal.
  ///
  /// \tparam Node Same as Map's \c Node template parameter.  Defaults
  ///   to the default Kokkos Node type.
  ///
  /// Directory implements looking up the process IDs and local IDs
  /// corresponding to a given list of global IDs.  Each Map owns a
  /// Directory object that does this.  Map::getRemoteIndexList()
  /// calls the Map's directory's getDirectoryEntries() method
  /// directly.  Directory has three different ways to perform this
  /// lookup, depending on the kind of Map.
  ///
  /// 1. If the user's Map is not distributed (i.e., is serial or
  ///    locally replicated), then my process ID is the process ID for
  ///    all global IDs.  I can get the local ID (if requested)
  ///    directly from the user's Map via its \c getLocalElement()
  ///    method (which requires no communication).
  ///
  /// 2. If the user's Map is distributed but contiguous, then I build
  ///    an array (replicated on all processes) that maps from each
  ///    process ID to the minimum global ID that it owns.  This
  ///    trades time and communication for space (\f$P + 1\f$ entries
  ///    in the array if there are P processes in the communicator),
  ///    in that it allows lookups without communication (once the
  ///    array has been built).
  ///
  /// 3. If the user's Map is distributed and noncontiguous, then I
  ///    have to store a general mapping from global ID to (process
  ///    ID, local ID).  I can't afford to store the whole mapping
  ///    redundantly on all processes, so I distribute it using
  ///    another Map (the "directory Map").  This is a uniform
  ///    contiguous Map whose keys are the global IDs.
  ///
  /// This class is templated on the same \c LocalOrdinal and \c
  /// GlobalOrdinal types on which \c Map is templated.  Just as with
  /// Map, the \c LocalOrdinal type defaults to \c int if omitted.
  /// The \c GlobalOrdinal type defaults to the \c LocalOrdinal type,
  /// and the Node type defaults to Kokkos' default node type.
  ///
  /// \note (mfh 04 Jan 2012) To Epetra developers: This class
  ///   corresponds roughly to \c Epetra_Directory or \c
  ///   Epetra_BasicDirectory.  \c Epetra_BlockMap creates its \c
  ///   Epetra_Directory object on demand whenever the map's \c
  ///   RemoteIDList() method is called.  \c Tpetra::Map's
  ///   getRemoteIndexList() method assumes that the map's directory
  ///   already exists.  \c Epetra_Directory is an abstract interface
  ///   with one implementation (\c Epetra_BasicDirectory); \c
  ///   Tpetra::Directory is a concrete implementation.
  ///
  /// \note (mfh 10 May 2012) To Tpetra developers: The current design
  ///   of Directory is suboptimal.  Currently, the Directory has to
  ///   ask the Map what kind of Map it is (e.g., globally
  ///   distributed? contiguous?) in order to figure out how to
  ///   represent its data.  This is suboptimal, because the Map
  ///   already knows what kind of Map it is, based on which Map
  ///   constructor was invoked.  Thus, it should be the Map that
  ///   specifies which kind of Directory to create.  I've broken up
  ///   Directory into an interface and an implementation
  ///   (Details::Directory and its subclasses) in preparation for
  ///   making this change.
  template<class LocalOrdinal,
           class GlobalOrdinal = LocalOrdinal,
           class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Directory : public Teuchos::Describable {
  public:
    //! Type of the Map specialization to give to the constructor.
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

    //! @name Constructors/Destructor.
    //@{

    /// \brief Constructor.
    ///
    /// \param map [in] The Map object for which to create the Directory.
    ///
    /// \note This constructor is invoked by Map's constructor, using
    ///   the Map's <tt>this</tt> pointer as the input argument.
    explicit Directory (const Teuchos::RCP<const map_type>& map);

    //! Destructor.
    ~Directory ();

    /* Provide a shallow copy of this Directory for a different Node type, using a cloned Map.
       This is an advanced method that will not be called by end users. */
    template <class Node2>
    RCP<Directory<LocalOrdinal,GlobalOrdinal,Node2> >
    clone(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node2> > &clone_map) const
    {
      RCP<Directory<LocalOrdinal,GlobalOrdinal,Node2> > dir = rcp(new Directory<LocalOrdinal,GlobalOrdinal,Node2>());
      if (clone_map->isDistributed ()) {
        if (clone_map->isContiguous ()) {
          dir->impl_ = Teuchos::rcp_dynamic_cast<
                          const Details::DistributedContiguousDirectory<LocalOrdinal, GlobalOrdinal, Node> >
                              (impl_)->template clone<Node2>(clone_map);
        }
        else {
          dir->impl_ = Teuchos::rcp_dynamic_cast<
                          const Details::DistributedNoncontiguousDirectory<LocalOrdinal, GlobalOrdinal, Node> >
                              (impl_)->template clone<Node2>(clone_map);
        }
      }
      else {
        dir->impl_ = Teuchos::rcp_dynamic_cast<
                        const Details::ReplicatedDirectory<LocalOrdinal, GlobalOrdinal, Node> >
                            (impl_)->template clone<Node2>(clone_map);
      }
      return dir;
    }

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
    /// the corresponding list \c nodeIDs of the process ranks which
    /// own those GIDs.  Tpetra uses this to figure out the locations
    /// of nonlocal Map entries.
    ///
    /// \param globalIDs [in] List of global IDs to look up.
    ///
    /// \param nodeIDs [out] On input, an array view with the same
    ///   number of entries as \c globalIDs.  On output, nodeIDs[i] is
    ///   the ID of the process which owns globalIDs[i].  If
    ///   globalIDs[i] is not present in the directory (i.e., is not
    ///   owned by any process in the Directory's communicator),
    ///   nodeIDs[i] is -1.
    ///
    /// \return If at least one global ID was not present in the
    ///   directory, return IDNotPresent.  Otherwise, return
    ///   AllIDsPresent.
    ///
    /// \note If <tt>nodeIDs.size() != globalIDs.size()</tt>, then
    ///   this method throws a \c std::runtime_error exception.
    LookupStatus
    getDirectoryEntries (const Teuchos::ArrayView<const GlobalOrdinal>& globalIDs,
			 const Teuchos::ArrayView<int>& nodeIDs) const;

    /// \brief Given a global ID list, return a list of their owning
    ///   process IDs and their corresponding local IDs.
    ///
    /// Given a list \c globalIDs of global identifiers (GIDs), return
    /// the corresponding list \c nodeIDs of the process ranks which
    /// own those GIDs, as well as the list of the local identifiers
    /// (LIDs).  Tpetra uses this to figure out the locations of
    /// nonlocal Map entries.
    ///
    /// \param globalIDs [in] List of global IDs to look up.
    ///
    /// \param nodeIDs [out] On input, an array view with the same
    ///   number of entries as \c globalIDs.  On output, nodeIDs[i] is
    ///   the ID of the process which owns globalIDs[i].  If
    ///   globalIDs[i] is not present in the directory (i.e., is not
    ///   owned by any process in the Directory's communicator), then
    ///   <tt>nodeIDs[i] == -1</tt>.
    ///
    /// \param localIDs [out] On input, an array view with the same
    ///   number of entries as \c globalIDs.  On output, \c
    ///   localIDs[i] is the local identifier corresponding to the
    ///   global identifier \c globalIDs[i].  If globalIDs[i] is not
    ///   present in the directory, then <tt>localIDs[i] ==
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid()</tt>.
    ///
    /// \return If at least one global ID was not present in the
    ///   directory, return IDNotPresent.  Otherwise, return
    ///   AllIDsPresent.
    ///
    /// \note If <tt>nodeIDs.size() != globalIDs.size()</tt> or
    ///   <tt>localIDs.size() != globalIDs.size()</tt>, then
    ///   this method throws a \c std::runtime_error exception.
    LookupStatus
    getDirectoryEntries (const Teuchos::ArrayView<const GlobalOrdinal>& globalIDs,
			 const Teuchos::ArrayView<int>& nodeIDs,
			 const Teuchos::ArrayView<LocalOrdinal>& localIDs) const;
    //@}

  private:
    /// \brief Type of the (base class) implementation of this object.
    ///
    /// \note To implementers: Directory has different
    ///   implementations, depending on characteristics of the input
    ///   Map (e.g., locally replicated or globally distributed,
    ///   contiguous or noncontiguous).
    typedef Details::Directory<LocalOrdinal, GlobalOrdinal, Node> base_type;

    //! Implementation of this object.
    Teuchos::RCP<const base_type> impl_;

    //! Copy constructor: declared private but not defined on purpose.
    Directory (const Directory<LocalOrdinal, GlobalOrdinal, Node>& directory);

    template <class LO, class GO, class N> friend class Directory;

    //! Empty constructor, for delayed initialization by clone()
    Directory();

    //! Assignment operator: declared private but not defined on purpose.
    Directory<LocalOrdinal, GlobalOrdinal, Node>&
    operator= (const Directory<LocalOrdinal, GlobalOrdinal, Node>& source);

  }; // class Directory
} // namespace Tpetra

#endif // TPETRA_DIRECTORY_DECL_HPP

