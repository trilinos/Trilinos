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

#ifndef __Tpetra_DirectoryImpl_decl_hpp
#define __Tpetra_DirectoryImpl_decl_hpp

#include <Tpetra_ConfigDefs.hpp>
#include "Tpetra_TieBreak.hpp"

//
// mfh 13-15 May 2013: HAVE_TPETRA_DIRECTORY_SPARSE_MAP_FIX governs
// the fix for Bug 5822.  The fix is enabled by default.  To disable
// the fix, uncomment out the three lines below that undefine
// HAVE_TPETRA_DIRECTORY_SPARSE_MAP_FIX, and comment out the three
// lines below them that define that macro.
//
// mfh 23 Mar 2014: I want Bug 5822 to stay fixed, so I am removing
// all references to HAVE_TPETRA_DIRECTORY_SPARSE_MAP_FIX.  I hope no
// downstream code is using that macro, but just in case, I will leave
// it defined.

#ifndef HAVE_TPETRA_DIRECTORY_SPARSE_MAP_FIX
#  define HAVE_TPETRA_DIRECTORY_SPARSE_MAP_FIX 1
#endif // HAVE_TPETRA_DIRECTORY_SPARSE_MAP_FIX

#include <Tpetra_Details_FixedHashTable_decl.hpp>


namespace Tpetra {
  // Forward declaration.
  template <class LocalOrdinal, class GlobalOrdinal, class Node> class Map;

  namespace Details {
    // Forward declaration.
    template <class LocalOrdinal, class GlobalOrdinal> class TieBreak;

    /// \class Directory
    /// \brief Computes the local ID and process ID corresponding to given global IDs.
    ///
    /// \note To implementers: This class and its subclasses implement
    ///   Tpetra::Directory.  We separate out the interface
    ///   (Tpetra::Directory) from the implementation in order to keep
    ///   backwards compatibility of the interface.
    template<class LocalOrdinal, class GlobalOrdinal, class NodeType>
    class Directory : public Teuchos::Describable {
    public:
      typedef LocalOrdinal local_ordinal_type;
      typedef GlobalOrdinal global_ordinal_type;
      typedef NodeType node_type;
      typedef ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, NodeType> map_type;

      /// \brief Constructor.
      ///
      /// Subclasses' constructors may only accept the Map to check
      /// its properties or to extract data from it in some way.  They
      /// may <i>not</i> keep a reference to the Map.  This prevents
      /// circular references, since the Map itself owns the
      /// Directory.
      Directory ();

      /// Find process IDs and (optionally) local IDs for the given global IDs.
      ///
      /// \pre nodeIDs.size() == globalIDs.size()
      /// \pre ! computeLIDs || localIDs.size() == globalIDs.size()
      ///
      /// \param map [in] The Directory's Map.  This must be the same
      ///   as given to the Directory's constructor.  Directory may
      ///   not keep a reference to the Map, in order to avoid
      ///   circular references between a Map and its Directory.
      ///
      /// \param globalIDs [in] The global IDs for which to find process
      ///   IDs (and optionally local IDs).
      ///
      /// \param nodeIDs [out] The process IDs corresponding to the
      ///   given global IDs.  If a global ID does not belong to any
      ///   process, the corresponding entry of nodeIDs will be -1.
      ///
      /// \param localIDs [out] If computeLIDs is true, we fill this
      ///   with the local IDs corresponding to the given global IDs.
      ///   If a given global ID does not correspond to a local ID, the
      ///   corresponding entry will be
      ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
      ///
      /// \param computeLIDs [in] Whether to fill in localIDs.
      ///
      /// \return If at least one GID was not present in the directory,
      ///   return IDNotPresent.  Otherwise, return AllIDsPresent.
      ///
      /// \note To implementers: The implementation of this method
      ///   first performs input validation, then invokes
      ///   getEntriesImpl() (implemented in the subclass) to do the
      ///   work.
      LookupStatus
      getEntries (const map_type& map,
                  const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                  const Teuchos::ArrayView<int> &nodeIDs,
                  const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                  const bool computeLIDs) const;

      /// \brief Whether the Directory's input Map is (globally) one to one.
      ///
      /// This method should always be treated as a collective on all
      /// processes in the given communicator, which must be the same
      /// as the input Map's communicator.  Not all implementations
      /// necessarily communicate.
      virtual bool isOneToOne (const Teuchos::Comm<int>& comm) const = 0;

    protected:
      //! Actually do the work of getEntries(), with no input validation.
      virtual LookupStatus
      getEntriesImpl (const map_type& map,
                      const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                      const Teuchos::ArrayView<int> &nodeIDs,
                      const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                      const bool computeLIDs) const = 0;
    };

    /// \class ReplicatedDirectory
    /// \brief Implementation of Directory for a locally replicated Map.
    template<class LocalOrdinal, class GlobalOrdinal, class NodeType>
    class ReplicatedDirectory :
      public Directory<LocalOrdinal, GlobalOrdinal, NodeType> {
    public:
      typedef Directory<LocalOrdinal, GlobalOrdinal, NodeType> base_type;
      typedef typename base_type::map_type map_type;

      //! Constructor (that takes a Map).
      ReplicatedDirectory (const map_type& map);

      //! Constructor (that takes no arguments).
      ReplicatedDirectory ();

      virtual bool isOneToOne (const Teuchos::Comm<int>& comm) const;

      template <class Node2>
      Directory<LocalOrdinal,GlobalOrdinal,Node2>*
      clone (const ::Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node2>& cloneMap) const
      {
        typedef ReplicatedDirectory<LocalOrdinal,GlobalOrdinal,Node2> Dir2;
        return new Dir2 (cloneMap);
      }

      //! @name Implementation of Teuchos::Describable.
      //@{

      //! A one-line human-readable description of this object.
      std::string description () const;
      //@}
    protected:
      //! Find process IDs and (optionally) local IDs for the given global IDs.
      LookupStatus
      getEntriesImpl (const map_type& map,
                      const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                      const Teuchos::ArrayView<int> &nodeIDs,
                      const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                      const bool computeLIDs) const;

    private:
      //! The number of process(es) in the input Map's communicator.
      const int numProcs_;
    };


    /// \class ContiguousUniformDirectory
    /// \brief Implementation of Directory for a contiguous, uniformly distributed Map.
    ///
    /// The Map may have any number of processes starting with one.
    /// Since the entries are uniformly distributed over the
    /// processes, this implementation of Directory can compute which
    /// process owns a GID (and the GID's corresponding LID) in
    /// \f$O(1)\f$ space and time.
    template<class LocalOrdinal, class GlobalOrdinal, class NodeType>
    class ContiguousUniformDirectory :
      public Directory<LocalOrdinal, GlobalOrdinal, NodeType> {
    private:
      // This friend declaration lets us implement clone().
      template <class LO, class GO, class N> friend class ContiguousUniformDirectory;

      //! Empty constructor for use by clone()
      ContiguousUniformDirectory () {}

    public:
      typedef Directory<LocalOrdinal, GlobalOrdinal, NodeType> base_type;
      typedef typename base_type::map_type map_type;

      //! Constructor.
      ContiguousUniformDirectory (const map_type& map);

      virtual bool isOneToOne (const Teuchos::Comm<int>&) const {
        return true;
      }

      template <class Node2>
      Directory<LocalOrdinal,GlobalOrdinal,Node2>*
      clone (const ::Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node2>& cloneMap) const
      {
        typedef ContiguousUniformDirectory<LocalOrdinal,GlobalOrdinal,Node2> Dir2;
        return new Dir2 (cloneMap);
      }

      //! @name Implementation of Teuchos::Describable.
      //@{

      //! A one-line human-readable description of this object.
      std::string description () const;
      //@}

    protected:
      //! Find process IDs and (optionally) local IDs for the given global IDs.
      LookupStatus
      getEntriesImpl (const map_type& map,
                      const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                      const Teuchos::ArrayView<int> &nodeIDs,
                      const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                      const bool computeLIDs) const;
    };


    /// \class DistributedContiguousDirectory
    /// \brief Implementation of Directory for a distributed contiguous Map.
    template<class LocalOrdinal, class GlobalOrdinal, class NodeType>
    class DistributedContiguousDirectory :
      public Directory<LocalOrdinal, GlobalOrdinal, NodeType> {
    private:
      template <class LO, class GO, class N> friend class DistributedContiguousDirectory;

      //! Empty constructor for use by clone()
      DistributedContiguousDirectory () {}

    public:
      typedef Directory<LocalOrdinal, GlobalOrdinal, NodeType> base_type;
      typedef typename base_type::map_type map_type;

      //! Constructor.
      DistributedContiguousDirectory (const map_type& map);

      virtual bool isOneToOne (const Teuchos::Comm<int>&) const {
        return true;
      }

      template <class Node2>
      Directory<LocalOrdinal,GlobalOrdinal,Node2>*
      clone (const ::Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node2>& cloneMap) const
      {
        typedef DistributedContiguousDirectory<LocalOrdinal,GlobalOrdinal,Node2> Dir2;
        Dir2* dir = new Dir2 (cloneMap);
        try {
          dir->allMinGIDs_ = allMinGIDs_;
        } catch (std::exception& e) {
          delete dir; // clean up just in case assignment throws (it shouldn't)
          throw;
        }
        return dir;
      }

      //! @name Implementation of Teuchos::Describable.
      //@{

      //! A one-line human-readable description of this object.
      std::string description () const;
      //@}

    protected:
      //! Find process IDs and (optionally) local IDs for the given global IDs.
      LookupStatus
      getEntriesImpl (const map_type& map,
                      const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                      const Teuchos::ArrayView<int> &nodeIDs,
                      const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                      const bool computeLIDs) const;

    private:
      /// \brief Minimum global ID for each process in the communicator.
      ///
      /// This array is only valid if the user's Map (\c map_) is
      /// distributed and contiguous.  Otherwise, it remains empty.  It
      /// is allocated in the constructor if necessary.
      ///
      /// This array has map_->getComm ()->getSize ()+1 entries.  Entry
      /// i contains the minimum global identifier (GID) of process i in
      /// map_'s communicator.  The last entry contains the maximum GID
      /// in the directory.
      ///
      /// The directory uses this array to map from GID to process ID,
      /// when the GIDs are distributed contiguously in increasing order
      /// over the processes.  This array allows the directory to
      /// compute the mapping locally, without communication, for any
      /// given GID, whether or not it is owned by the local process.
      ///
      /// \note To implementers: This is a potential memory bottleneck
      ///   if the number of processes P is large and the allowed memory
      ///   usage per process is small.  This should only be a problem
      ///   if \f$N/P \gg P\f$, where N is the global number of
      ///   elements.  In this case, it would be more memory scalable to
      ///   use reductions or scans to figure out who owns what.
      Teuchos::ArrayRCP<GlobalOrdinal> allMinGIDs_;
    };

    /// \class DistributedNoncontiguousDirectory
    /// \brief Implementation of Directory for a distributed noncontiguous Map.
    template<class LocalOrdinal, class GlobalOrdinal, class NodeType>
    class DistributedNoncontiguousDirectory :
      public Directory<LocalOrdinal, GlobalOrdinal, NodeType> {
    private:
      template <class LO, class GO, class N> friend class DistributedNoncontiguousDirectory;
      //! Private constructor for post-contruction initialization in clone()
      DistributedNoncontiguousDirectory () {}

    public:
      typedef Tpetra::Details::TieBreak<LocalOrdinal, GlobalOrdinal> tie_break_type;
      typedef Directory<LocalOrdinal, GlobalOrdinal, NodeType> base_type;
      typedef typename base_type::map_type map_type;

      //! Constructor.
      DistributedNoncontiguousDirectory (const map_type& map);

      //! Constructor.
      DistributedNoncontiguousDirectory (const map_type& map,
                                         const tie_break_type& tie_break);

      virtual bool isOneToOne (const Teuchos::Comm<int>& comm) const;

      template <class Node2>
      Directory<LocalOrdinal,GlobalOrdinal,Node2>*
      clone (const ::Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node2>& cloneMap) const
      {
        typedef DistributedNoncontiguousDirectory<LocalOrdinal,GlobalOrdinal,Node2> Dir2;
        typedef ::Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node2> output_map_type;
        Dir2* dir = new Dir2 (cloneMap);

        // This method returns a raw pointer.  Thus, take care to
        // check whether intermediate operations succeed, so that we
        // don't leak memory if they don't.
        RCP<const output_map_type> outDirMap;
        try {
          outDirMap = directoryMap_->template clone<Node2> (cloneMap.getNode ());
        }
        catch (...) {
          outDirMap = Teuchos::null; // deallocate
          throw;
        }

        dir->directoryMap_ = outDirMap;
        dir->PIDs_ = PIDs_;
        dir->LIDs_ = LIDs_;
        dir->lidToPidTable_ = lidToPidTable_;
        dir->lidToLidTable_ = lidToLidTable_;
        dir->useHashTables_ = useHashTables_;
        return dir;
      }

      //! @name Implementation of Teuchos::Describable.
      //@{

      //! A one-line human-readable description of this object.
      std::string description () const;
      //@}
    protected:
      //! Find process IDs and (optionally) local IDs for the given global IDs.
      LookupStatus
      getEntriesImpl (const map_type& map,
                      const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                      const Teuchos::ArrayView<int> &nodeIDs,
                      const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                      const bool computeLIDs) const;
    private:
      /// \brief Initialization routine that unifies the implementation of
      ///        the two constructors
      ///
      /// If the pointer to the TieBreak object is null this proceeds using
      /// a simple ordering to break any ownership ties. Otherwise the
      /// tie_break object is used to determine ownership.
      void
      initialize (const map_type& map,
                  Teuchos::Ptr<const tie_break_type> tie_break);

      /// \brief Whether the Directory is "locally" one to one.
      ///
      /// This means that the calling process' Directory does not own
      /// GIDs with multiple ownership on different processes.  If
      /// this method returns true on all processes in the Directory's
      /// communicator, then the Directory's input Map is one to one.
      /// If it returns false on at least one process in the
      /// Directory's communicator, then the Directory's input Map is
      /// <i>not</i> one to one.
      ///
      /// This method is protected because it is an implementation
      /// detail of isOneToOne().
      bool isLocallyOneToOne () const {
        return locallyOneToOne_;
      }

      /// \brief This Directory's Map which describes the distribution of its data.
      ///
      /// The Directory Map describes where to find the distributed
      /// global IDs (GIDs).  This is a different object from the Map
      /// with which this Directory was created.
      ///
      /// We can't afford to store the whole directory redundantly on
      /// each process, so we distribute it.  This Map describes the
      /// distribution of the Directory.  It is a uniform contiguous
      /// map to prevent infinite recursion (since Map's constructor
      /// creates a Directory for the general case of a noncontiguous
      /// map).  The data which this Map distributes are PIDs_ and
      /// LIDs_ (see below): the process IDs resp. local IDs.  The
      /// "keys" or indices of this Map are the global IDs.  Thus,
      /// this Map has a range of elements from the minimum to the
      /// maximum GID of the user's Map, and its indexBase is the
      /// minimum GID over all processes in the user's Map.
      Teuchos::RCP<const map_type> directoryMap_;

      //! \name First of two implementations of Directory storage
      //@{

      /// \brief Mapping from Directory Map LID to input Map PID.
      ///
      /// Array of the same length as the local number of entries in
      /// directoryMap_, containing the process IDs corresponding to the
      /// GIDs owned by the Directory Map on this process.
      Teuchos::ArrayRCP<int> PIDs_;

      /// \brief Mapping from Directory Map LID to input Map LID.
      ///
      /// Array of the same length as the local number of entries in
      /// directoryMap_, containing the LIDs corresponding to the GIDs
      /// owned by the Directory Map on this process.
      Teuchos::ArrayRCP<LocalOrdinal> LIDs_;

      //@}
      //! \name Second of two implementations of Directory storage
      //@{

      /// \brief Mapping from Directory Map LID to input Map PID.
      ///
      /// This hash table implements a mapping from an LID in the
      /// Directory Map (corresponding to a GID in the input Map) to
      /// the GID's owning PID in the input Map.
      Teuchos::RCP<Details::FixedHashTable<LocalOrdinal, int> > lidToPidTable_;

      /// \brief Mapping from Directory Map LID to input Map LID.
      ///
      /// This hash table implements a mapping from an LID in the
      /// Directory Map (corresponding to a GID in the input Map), to
      /// the GID's LID in the input Map on the GID's owning process.
      Teuchos::RCP<Details::FixedHashTable<LocalOrdinal, LocalOrdinal> > lidToLidTable_;
      //@}

      /// \brief The result of the first call to isOneToOne() on this object.
      ///
      /// If isOneToOne() has not yet been called on this object
      /// before, the value is ONE_TO_ONE_NOT_CALLED_YET.  Otherwise,
      /// if it returned false, the value is ONE_TO_ONE_FALSE; if it
      /// returned true, the value is ONE_TO_ONE_TRUE.
      mutable enum EOneToOneResult {
        ONE_TO_ONE_NOT_CALLED_YET,
        ONE_TO_ONE_FALSE,
        ONE_TO_ONE_TRUE
      } oneToOneResult_;

      /// \brief Whether this process is locally one-to-one.
      ///
      /// See documentation of isLocallyOneToOne() for a definition.
      bool locallyOneToOne_;

      /// \brief Whether this process is using hash tables for Directory storage.
      ///
      /// Directory may use either arrays or hash tables for Directory
      /// storage.  Lookups with arrays are faster, but hash tables
      /// use less memory if the input Map is sparse.  The choice of
      /// implementation is decided locally, and may differ from
      /// process to process.
      bool useHashTables_;
    };
  } // namespace Details
} // namespace Tpetra

#endif // __Tpetra_DirectoryImpl_decl_hpp
