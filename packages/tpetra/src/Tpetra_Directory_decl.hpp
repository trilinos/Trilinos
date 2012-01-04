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

namespace Tpetra {

  /// \class Directory
  /// \brief Implement mapping from global ID to process ID and local ID.
  ///
  /// This class is an implementation detail of \c Map.  It is mainly
  /// of interest to Tpetra developers and does not normally appear in
  /// users' code.
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
  ///    trades time and communication for space (\f$P + O(1)\f$
  ///    entries if there are P processes in the communicator), in
  ///    that it allows lookups without communication (once the array
  ///    has been built).
  ///
  /// 3. If the user's Map is distributed and noncontiguous, then I
  ///    have to store a general mapping from global ID to (process
  ///    ID, local ID).  I can't afford to store the whole mapping
  ///    redundantly on all processes, so I distribute it using
  ///    another Map (the "directory Map").  This is a uniform linear
  ///    Map whose keys are the global IDs.
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
  template<class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Directory : public Teuchos::Describable {
  public:
    //! @name Constructors/Destructor.
    //@{ 

    /// \brief Constructor.
    ///
    /// \param map [in] The Map object for which to create the Directory.
    ///
    /// \note This constructor is invoked by Map's constructor, using
    ///   the Map's <tt>this</tt> pointer as the input argument.
    explicit Directory(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map);

    //! Destructor.
    ~Directory();

    //@}

    //! @name Query methods.
    //@{ 

    /// \brief Given a GID list, return the list of their owning process IDs.
    ///
    /// Given a list \c globalIDs of global identifiers (GIDs), return
    /// the corresponding list \c nodeIDs of the process ranks which
    /// own those GIDs.  Tpetra uses this to figure out the locations
    /// of nonlocal Map entries.
    ///
    /// \param globalIDs [in] List of GIDs to look up.
    ///
    /// \param nodeIDs [out] On input, an array view with the same
    ///   number of entries as \c globalIDs.  On output, nodeIDs[i] is
    ///   the ID of the process which owns globalIDs[i].  If
    ///   globalIDs[i] is not present in the directory (i.e., is not
    ///   owned by any process in the Directory's communicator),
    ///   nodeIDs[i] is -1.
    ///
    /// \return If at least one GID was not present in the directory,
    ///   return IDNotPresent.  Otherwise, return AllIDsPresent.
    ///
    /// \note If <tt>nodeIDs.size() != globalIDs.size()</tt>, then
    ///   this method throws a \c std::runtime_error exception.
    LookupStatus getDirectoryEntries(const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
                                     const Teuchos::ArrayView<int> &nodeIDs) const;

    /// \brief Given a GID list, return a list of their owning process IDs and their corresponding LIDs.
    ///
    /// Given a list \c globalIDs of global identifiers (GIDs), return
    /// the corresponding list \c nodeIDs of the process ranks which
    /// own those GIDs, as well as the list of the local identifiers
    /// (LIDs).  Tpetra uses this to figure out the locations of
    /// nonlocal Map entries.
    ///
    /// \param globalIDs [in] List of GIDs to look up.
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
    /// \return If at least one GID was not present in the directory,
    ///   return IDNotPresent.  Otherwise, return AllIDsPresent.
    ///
    /// \note If <tt>nodeIDs.size() != globalIDs.size()</tt> or 
    ///   <tt>localIDs.size() != globalIDs.size()</tt>, then
    ///   this method throws a \c std::runtime_error exception.
    LookupStatus getDirectoryEntries(const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
                                     const Teuchos::ArrayView<int> &nodeIDs, 
                                     const Teuchos::ArrayView<LocalOrdinal> &localIDs) const;
    //@}
    
  private:
    //! The Map with which this Directory was created.
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > map_;

    /// \brief This Directory's Map which describes the distribution of its data.
    ///
    /// This Map is only instantiated if the user's Map (\c map_) is
    /// distributed and noncontiguous.  Otherwise, it remains null.
    /// It is instantiated in \c generateDirectory(), which is invoked
    /// by the constructor if necessary.
    ///
    /// We can't afford to store the whole directory redundantly on
    /// each process, so we distribute it.  This Map describes the
    /// distribution of the Directory.  It is a uniform contiguous map
    /// to prevent infinite recursion (since Map's constructor creates
    /// a Directory for the general case of a noncontiguous map).  The
    /// data which this Map distributes are nodeIDs_ and LIDs_: the
    /// process IDs and local IDs.  The "keys" or indices of this Map
    /// are the global IDs.  Thus, this Map has a range of elements
    /// from the minimum to the maximum GID of the user's Map, and its
    /// indexBase is the minimum GID over all processes in the user's
    /// Map.
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > directoryMap_;

    //! The communicator over which the Directory is distributed.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    /// \brief Minimum global ID for each node in the communicator.
    ///
    /// This array is only valid if the user's Map (\c map_) is
    /// distributed and contiguous.  Otherwise, it remains empty.  It
    /// is allocated in the constructor if necessary.
    ///
    /// This array has comm_->getSize()+1 entries.  Entry i contains
    /// the minimum global identifier (GID) of process i in the
    /// communicator comm_.  The last entry contains the maximum GID
    /// in the directory.  
    ///
    /// The directory uses this array to map from GID to process ID,
    /// when the GIDs are distributed contiguously in increasing order
    /// over the processes.  This array allows the directory to
    /// compute the mapping locally, without communication, for any
    /// given GID, whether or not it is owned by the local process.
    ///
    /// \note This is a potential memory bottleneck if the number of
    ///   processes P is large and the allowed memory usage per
    ///   process is small.
    std::vector<GlobalOrdinal> allMinGIDs_; 

    /// Array of the same length as the local number of entries in
    /// directoryMap_.  This array is only allocated and used if the
    /// user's map is distributed and noncontiguous.
    std::vector<int> nodeIDs_;
    /// Array of the same length as the local number of entries in
    /// directoryMap_.  This array is only allocated and used if the
    /// user's map is distributed and noncontiguous.
    std::vector<LocalOrdinal> LIDs_;

    //! Copy constructor: declared private but not defined on purpose.
    Directory (const Directory<LocalOrdinal,GlobalOrdinal,Node> &directory);

    //! Assignment operator: declared private but not defined on purpose.
    Directory<LocalOrdinal,GlobalOrdinal,Node> & operator = (const Directory<LocalOrdinal,GlobalOrdinal,Node> &source);

    /// \brief Common code for both versions of \c getDirectoryEntries().
    ///
    /// \param globalIDs [in] The global IDs to look up.
    /// \param nodeIDs [out] The process IDs corresponding to the
    ///   given global IDs.
    /// \param localIDs [out] If computeLIDs is true, fill with the
    ///   local IDs corresponding to the given global IDs.
    /// \param computeLIDs [in] Whether to fill in localIDs.
    ///
    /// \return If at least one GID was not present in the directory,
    ///   return IDNotPresent.  Otherwise, return AllIDsPresent.
    ///
    /// \note If <tt>nodeIDs.size() != globalIDs.size()</tt>, or if
    ///   computeLIDs is true and <tt>localIDs.size() !=
    ///   globalIDs.size()</tt>, then this method throws a \c
    ///   std::runtime_error exception.
    LookupStatus getEntries(const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
                            const Teuchos::ArrayView<int> &nodeIDs, 
                            const Teuchos::ArrayView<LocalOrdinal> &localIDs, 
                            bool computeLIDs) const;

    /// \brief Set up directory, if the user's Map is distributed and noncontiguous.
    ///
    /// This method is called in the constructor, and is only called
    /// if the user's Map is distributed and noncontiguous.  It sets
    /// up directoryMap_, nodeIDs_, and LIDs_.
    void generateDirectory();

  }; // class Directory

} // namespace Tpetra

#endif // TPETRA_DIRECTORY_DECL_HPP

