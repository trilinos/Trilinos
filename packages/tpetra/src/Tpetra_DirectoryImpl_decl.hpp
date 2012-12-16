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

#include "Tpetra_ConfigDefs.hpp"

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class Node> class Map;

  namespace Details {
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
      typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, NodeType> map_type;

      //! Constructor.
      Directory (const Teuchos::RCP<const map_type>& map);

      /// Find process IDs and (optionally) local IDs for the given global IDs.
      ///
      /// \pre nodeIDs.size() == globalIDs.size()
      /// \pre ! computeLIDs || localIDs.size() == globalIDs.size()
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
      getEntries (const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                  const Teuchos::ArrayView<int> &nodeIDs,
                  const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                  const bool computeLIDs) const;

    protected:
      //! Actually do the work of getEntries(), with no input validation.
      virtual LookupStatus
      getEntriesImpl (const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                      const Teuchos::ArrayView<int> &nodeIDs,
                      const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                      const bool computeLIDs) const = 0;

      //! Get the Map with which this object was created.
      Teuchos::RCP<const map_type> getMap () const { return map_; }

      //! Set the Map for this object, for post-contructor initialization.
      void setMap (const Teuchos::RCP<const map_type> &map) { map_ = map; }

      //! Empty constructor for post-constructor initialization
      Directory() {}

    private:
      //! The Map with which this object was created.
      Teuchos::RCP<const map_type> map_;
    };

    /// \class ReplicatedDirectory
    /// \brief Implementation of Directory for a locally replicated Map.
    template<class LocalOrdinal, class GlobalOrdinal, class NodeType>
    class ReplicatedDirectory :
      public Directory<LocalOrdinal, GlobalOrdinal, NodeType> {
    public:
      typedef Directory<LocalOrdinal, GlobalOrdinal, NodeType> base_type;
      typedef typename base_type::map_type map_type;

      //! Constructor.
      ReplicatedDirectory (const Teuchos::RCP<const map_type>& map);

      template <class Node2>
      RCP<Directory<LocalOrdinal,GlobalOrdinal,Node2> >
      clone(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node2> > &clone_map) const
      {
        // this class has only the map as its data; cloning is trivial
        return rcp(new ReplicatedDirectory<LocalOrdinal,GlobalOrdinal,Node2>(clone_map));
      }

      //! @name Implementation of Teuchos::Describable.
      //@{

      //! A one-line human-readable description of this object.
      std::string description () const;
      //@}
    protected:
      //! Find process IDs and (optionally) local IDs for the given global IDs.
      LookupStatus
      getEntriesImpl (const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                      const Teuchos::ArrayView<int> &nodeIDs,
                      const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                      const bool computeLIDs) const;
    };

    /// \class ReplicatedDirectory
    /// \brief Implementation of Directory for a distributed contiguous Map.
    template<class LocalOrdinal, class GlobalOrdinal, class NodeType>
    class DistributedContiguousDirectory :
      public Directory<LocalOrdinal, GlobalOrdinal, NodeType> {
    private:
      template <class LO, class GO, class N> friend class DistributedContiguousDirectory;

      //! Empty constructor for used by clone()
      DistributedContiguousDirectory() {}

    public:
      typedef Directory<LocalOrdinal, GlobalOrdinal, NodeType> base_type;
      typedef typename base_type::map_type map_type;

      //! Constructor.
      DistributedContiguousDirectory (const Teuchos::RCP<const map_type>& map);

      template <class Node2>
      RCP<Directory<LocalOrdinal,GlobalOrdinal,Node2> >
      clone(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node2> > &clone_map) const
      {
        typedef DistributedContiguousDirectory<LocalOrdinal,GlobalOrdinal,Node2> Dir2;
        RCP<Dir2> dir = rcp(new Dir2());
        dir->setMap ( clone_map );
        dir->allMinGIDs_ = allMinGIDs_;
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
      getEntriesImpl (const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
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
      DistributedNoncontiguousDirectory() {}

    public:
      typedef Directory<LocalOrdinal, GlobalOrdinal, NodeType> base_type;
      typedef typename base_type::map_type map_type;

      //! Constructor.
      DistributedNoncontiguousDirectory (const Teuchos::RCP<const map_type>& map);

      template <class Node2>
      RCP<Directory<LocalOrdinal,GlobalOrdinal,Node2> >
      clone(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node2> > &clone_map) const
      {
        typedef DistributedNoncontiguousDirectory<LocalOrdinal,GlobalOrdinal,Node2> Dir2;
        RCP<Dir2> dir = rcp(new Dir2());
        dir->setMap ( clone_map );
        dir->directoryMap_ = directoryMap_->template clone<Node2>(clone_map->getNode());
        dir->nodeIDs_      = nodeIDs_;
        dir->LIDs_         = LIDs_;
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
      getEntriesImpl (const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs,
                      const Teuchos::ArrayView<int> &nodeIDs,
                      const Teuchos::ArrayView<LocalOrdinal> &localIDs,
                      const bool computeLIDs) const;
    private:
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
      /// map).  The data which this Map distributes are nodeIDs_ and
      /// LIDs_ (see below): the process IDs resp. local IDs.  The
      /// "keys" or indices of this Map are the global IDs.  Thus,
      /// this Map has a range of elements from the minimum to the
      /// maximum GID of the user's Map, and its indexBase is the
      /// minimum GID over all processes in the user's Map.
      Teuchos::RCP<const map_type> directoryMap_;

      /// Array of the same length as the local number of entries in
      /// directoryMap_, containing the process IDs corresponding to the
      /// GIDs owned by the Directory Map on this process.
      Teuchos::ArrayRCP<int> nodeIDs_;

      /// Array of the same length as the local number of entries in
      /// directoryMap_, containing the LIDs corresponding to the GIDs
      /// owned by the Directory Map on this process.
      Teuchos::ArrayRCP<LocalOrdinal> LIDs_;
    };
  } // namespace Details
} // namespace Tpetra

#endif // __Tpetra_DirectoryImpl_decl_hpp
