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

#ifndef TPETRANEW_DIRECTORYIMPL_HPP
#define TPETRANEW_DIRECTORYIMPL_HPP

/// \file TpetraNew_DirectoryImpl.hpp
/// \brief Declaration of implementation details of TpetraNew::Directory.

#include "TpetraNew_Map.hpp"
#include "Tpetra_TieBreak.hpp"
#include "Tpetra_Details_FixedHashTable_decl.hpp"
#include "Teuchos_Describable.hpp"

// #ifndef DOXYGEN_SHOULD_SKIP_THIS
// namespace Teuchos {
//   template<class OrdinalType> class Comm; // forward declaration
// } // namespace Teuchos

// namespace TpetraNew {
//   class Map; // forward declaration
// } // namespace TpetraNew
// #endif // DOXYGEN_SHOULD_SKIP_THIS

namespace TpetraNew {
namespace Details {

  /// \class Directory
  /// \brief Implementation of TpetraNew::Directory. Given global
  ///   indices, compute their corresponding process ranks and
  ///   (optionally) local indices.
  ///
  /// \note To implementers: This class and its subclasses implement
  ///   TpetraNew::Directory.  We separate out the interface
  ///   (TpetraNew::Directory) from the implementation in order to
  ///   keep backwards compatibility of the interface.
  class Directory : public ::Teuchos::Describable {
  public:
    using map_type = ::TpetraNew::Map;
    //! The type of local indices.
    using local_ordinal_type = map_type::local_ordinal_type;
    /// \brief The type of global indices.
    ///
    /// Depending on your configuration, this may be either 32 bits or
    /// 64 bits.  If 32 bits, then the <i>global</i> number of rows or
    /// columns in any of your data structures (e.g., sparse matrices)
    /// may be no more than \c INT_MAX, namely \f$2^{31} - 1\f$ (about
    /// two billion).  If you want to solve larger problems, you must
    /// reconfigure and rebuild Trilinos with the appropriate setting
    /// to make this 64 bits.
    using global_ordinal_type = map_type::global_ordinal_type;
      
    //! The Kokkos execution space.
    using execution_space = map_type::execution_space;
      
    //! The Kokkos memory space.
    using memory_space = map_type::memory_space;

    /// \brief The Kokkos device type over which to allocate Views and
    ///   perform work.
    ///
    /// A Kokkos::Device is an (execution_space, memory_space) pair.
    /// It defines where the Map's data live, and where Map might
    /// choose to execute parallel kernels.
    using device_type = map_type::device_type;

    /// \brief Constructor.
    ///
    /// Subclasses' constructors may only accept the Map to check
    /// its properties or to extract data from it in some way.  They
    /// may <i>not</i> keep a reference to the Map.  This prevents
    /// circular references, since the Map itself owns the
    /// Directory.
    Directory () = default;

    /// \brief Given global indices, find their corresponding
    ///   process ranks and (optionally) local indices.
    ///
    /// \pre <tt>processRanks.size() == globalIndices.size()</tt>
    /// \pre <tt>! computeLocalIndices || localIndices.size() == globalIndices.size()</tt>
    ///
    /// \param map [in] The Directory's Map.  This must be the same
    ///   as given to the Directory's constructor.  Directory may
    ///   not keep a reference to the Map, in order to avoid
    ///   circular references between a Map and its Directory.
    ///
    /// \param globalIndices [in] The global indices for which to
    ///   find process ranks (and optionally local indices).
    ///
    /// \param processRanks [out] The process ranks corresponding to
    ///   the given global indices.  If a global index does not
    ///   belong to any process, the corresponding entry of
    ///   processRanks will be -1.
    ///
    /// \param localIndices [out] If computeLocalIndices is true, we
    ///   fill this with the local indices corresponding to the
    ///   given global IDs.  If a given global index does not
    ///   correspond to a local local index, the corresponding entry
    ///   of localIndices will be
    ///   Teuchos::OrdinalTraits<local_ordinal_type>::invalid().
    ///
    /// \param computeLocalIndices [in] Whether to fill in localIndices.
    ///
    /// \return If at least one global index was not present in the
    ///   directory, return IDNotPresent.  Otherwise, return
    ///   AllIDsPresent.
    ///
    /// \note To implementers: The implementation of this method
    ///   first performs input validation, then invokes
    ///   getEntriesImpl() (implemented in the subclass) to do the
    ///   work.
    ::Tpetra::LookupStatus
    getEntries (const map_type& map,
		const Teuchos::ArrayView<const global_ordinal_type>& globalIndices,
		const Teuchos::ArrayView<int>& processRanks,
		const Teuchos::ArrayView<local_ordinal_type>& localIndices,
		const bool computeLocalIndices) const;

    /// \brief Whether the Directory's input Map is (globally) one to one.
    ///
    /// This method should always be treated as a collective on all
    /// processes in the given communicator, which must be the same
    /// as the input Map's communicator.  Not all implementations
    /// necessarily communicate.
    virtual bool isOneToOne (const Teuchos::Comm<int>& comm) const = 0;

  protected:
    //! Actually do the work of getEntries(), with no input validation.
    virtual ::Tpetra::LookupStatus
    getEntriesImpl (const map_type& map,
		    const Teuchos::ArrayView<const global_ordinal_type>& globalIndices,
		    const Teuchos::ArrayView<int>& processRanks,
		    const Teuchos::ArrayView<local_ordinal_type>& localIndices,
		    const bool computeLocalIndices) const = 0;
  };

  /// \class ReplicatedDirectory
  /// \brief Implementation of Directory for a locally replicated Map.
  class ReplicatedDirectory : public Directory {
  public:
    using base_type = Directory;
    using map_type = typename base_type::map_type;

    //! The type of local indices.
    using local_ordinal_type = base_type::local_ordinal_type;
    //! The type of global indices.
    using global_ordinal_type = base_type::global_ordinal_type;
    //! The Kokkos execution space.
    using execution_space = base_type::execution_space;
    //! The Kokkos memory space.
    using memory_space = base_type::memory_space;
    //! The Kokkos::Device specialization.
    using device_type = base_type::device_type;
      
    //! Constructor (that takes a Map).
    ReplicatedDirectory (const map_type& map);

    //! Constructor (that takes no arguments).
    ReplicatedDirectory ();

    virtual bool isOneToOne (const Teuchos::Comm<int>& comm) const;

    //! @name Implementation of Teuchos::Describable.
    //@{

    //! A one-line human-readable description of this object.
    std::string description () const;
    //@}
  protected:
    //! Find process IDs and (optionally) local IDs for the given global IDs.
    ::Tpetra::LookupStatus
    getEntriesImpl (const map_type& map,
		    const Teuchos::ArrayView<const global_ordinal_type>& globalIndices,
		    const Teuchos::ArrayView<int>& processRanks,
		    const Teuchos::ArrayView<local_ordinal_type>& localIndices,
		    const bool computeLocalIndices) const;

  private:
    //! The number of process(es) in the input Map's communicator.
    const int numProcs_;
  };

  /// \class ContiguousUniformDirectory
  /// \brief Implementation of Directory for a contiguous, uniformly
  ///   distributed Map.
  ///
  /// The Map may have any number of processes starting with one.
  /// Since the entries are uniformly distributed over the
  /// processes, this implementation of Directory can compute which
  /// process owns a GID (and the GID's corresponding LID) in
  /// \f$O(1)\f$ space and time.
  class ContiguousUniformDirectory : public Directory {
  public:
    using base_type = Directory;
    using map_type = typename base_type::map_type;
    //! The type of local indices.
    using local_ordinal_type = base_type::local_ordinal_type;
    //! The type of global indices.
    using global_ordinal_type = base_type::global_ordinal_type;
    //! The Kokkos execution space.
    using execution_space = base_type::execution_space;
    //! The Kokkos memory space.
    using memory_space = base_type::memory_space;
    //! The Kokkos::Device specialization.
    using device_type = base_type::device_type;

    //! Constructor.
    ContiguousUniformDirectory (const map_type& map);

    virtual bool isOneToOne (const Teuchos::Comm<int>&) const {
      return true;
    }

    //! @name Implementation of Teuchos::Describable.
    //@{

    //! A one-line human-readable description of this object.
    std::string description () const;
    //@}

  protected:
    //! Find process IDs and (optionally) local IDs for the given global IDs.
    ::Tpetra::LookupStatus
    getEntriesImpl (const map_type& map,
		    const Teuchos::ArrayView<const global_ordinal_type>& globalIndices,
		    const Teuchos::ArrayView<int>& processRanks,
		    const Teuchos::ArrayView<local_ordinal_type>& localIndices,
		    const bool computeLocalIndices) const;
  };

  /// \class DistributedContiguousDirectory
  /// \brief Implementation of Directory for a distributed contiguous Map.
  class DistributedContiguousDirectory : public Directory {
  public:
    using base_type = Directory;
    using map_type = typename base_type::map_type;
    //! The type of local indices.
    using local_ordinal_type = base_type::local_ordinal_type;
    //! The type of global indices.
    using global_ordinal_type = base_type::global_ordinal_type;
    //! The Kokkos execution space.
    using execution_space = base_type::execution_space;
    //! The Kokkos memory space.
    using memory_space = base_type::memory_space;
    //! The Kokkos::Device specialization.
    using device_type = base_type::device_type;

    //! Constructor.
    DistributedContiguousDirectory (const map_type& map);

    virtual bool isOneToOne (const Teuchos::Comm<int>&) const {
      return true;
    }

    //! @name Implementation of Teuchos::Describable.
    //@{

    //! A one-line human-readable description of this object.
    std::string description () const;
    //@}

  protected:
    /// \brief Given global indices, find their corresponding
    ///   process ranks and (optionally) local indices.
    ::Tpetra::LookupStatus
    getEntriesImpl (const map_type& map,
		    const Teuchos::ArrayView<const global_ordinal_type>& globalIndices,
		    const Teuchos::ArrayView<int>& processRanks,
		    const Teuchos::ArrayView<local_ordinal_type>& localIndices,
		    const bool computeLocalIndices) const;

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
    Teuchos::ArrayRCP<global_ordinal_type> allMinGIDs_;
  };

  /// \class DistributedNoncontiguousDirectory
  /// \brief Implementation of Directory for a distributed noncontiguous Map.
  class DistributedNoncontiguousDirectory : public Directory {
  public:
    using base_type = Directory;
    using map_type = typename base_type::map_type;
    //! The type of local indices.
    using local_ordinal_type = base_type::local_ordinal_type;
    //! The type of global indices.
    using global_ordinal_type = base_type::global_ordinal_type;
    //! The Kokkos execution space.
    using execution_space = base_type::execution_space;
    //! The Kokkos memory space.
    using memory_space = base_type::memory_space;
    //! The Kokkos::Device specialization.
    using device_type = base_type::device_type;
    using tie_break_type = ::Tpetra::Details::TieBreak<local_ordinal_type, global_ordinal_type>;

    //! Constructor that uses the default TieBreak.
    DistributedNoncontiguousDirectory (const map_type& map);

    //! Constructor that uses a custom TieBreak.
    DistributedNoncontiguousDirectory (const map_type& map,
				       const tie_break_type& tie_break);

    virtual bool isOneToOne (const Teuchos::Comm<int>& comm) const;

    //! @name Implementation of Teuchos::Describable.
    //@{

    //! A one-line human-readable description of this object.
    std::string description () const;
    //@}
  protected:
    //! Find process IDs and (optionally) local IDs for the given global IDs.
    ::Tpetra::LookupStatus
    getEntriesImpl (const map_type& map,
		    const Teuchos::ArrayView<const global_ordinal_type>& globalIndices,
		    const Teuchos::ArrayView<int>& processRanks,
		    const Teuchos::ArrayView<local_ordinal_type>& localIndices,
		    const bool computeLocalIndices) const;
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
    Teuchos::ArrayRCP<local_ordinal_type> LIDs_;

    //@}
    //! \name Second of two implementations of Directory storage
    //@{
    /// \brief Mapping from Directory Map LID to input Map PID.
    ///
    /// This hash table implements a mapping from an LID in the
    /// Directory Map (corresponding to a GID in the input Map) to
    /// the GID's owning PID in the input Map.
    Teuchos::RCP< ::Tpetra::Details::FixedHashTable<local_ordinal_type, int, device_type>> lidToPidTable_;

    /// \brief Mapping from Directory Map LID to input Map LID.
    ///
    /// This hash table implements a mapping from an LID in the
    /// Directory Map (corresponding to a GID in the input Map), to
    /// the GID's LID in the input Map on the GID's owning process.
    Teuchos::RCP< ::Tpetra::Details::FixedHashTable<local_ordinal_type, local_ordinal_type, device_type>> lidToLidTable_;
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
} // namespace TpetraNew

#endif // TPETRANEW_DIRECTORYIMPL_HPP
