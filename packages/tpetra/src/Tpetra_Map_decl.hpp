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

#ifndef TPETRA_MAP_DECL_HPP
#define TPETRA_MAP_DECL_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

#ifdef HAVE_TPETRA_UNORDERED_MAP
#  include <unordered_map>
#endif // HAVE_TPETRA_UNORDERED_MAP

// enums and defines
#include "Tpetra_ConfigDefs.hpp"

/** \file Tpetra_Map_decl.hpp 

    The declarations for the class Tpetra::Map and related non-member constructors.
 */

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // Forward declaration of Directory.
  template <class LO, class GO, class N> class Directory;
#endif

  /** \brief A class for partitioning distributed objects.

   This class is templated on \c LocalOrdinal and \c GlobalOrdinal. 
   The \c GlobalOrdinal type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Map : public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    /** \brief Map constructor with Tpetra-defined contiguous uniform distribution.
     *
     *   The elements are distributed among nodes so that the subsets
     *   of global elements are non-overlapping and contiguous and as
     *   evenly distributed across the nodes as possible.
     */
    Map (global_size_t numGlobalElements, 
	 GlobalOrdinal indexBase, 
	 const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
	 LocalGlobal lg=GloballyDistributed, 
	 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

    /** \brief Map constructor with a user-defined contiguous distribution.
     * 
     *  The elements are distributed among the nodes so that the
     *  subsets of global elements are non-overlapping and contiguous.
     *  
     *  If numGlobalElements ==
     *  Teuchos::OrdinalTraits<global_size_t>::invalid(), the number
     *  of global elements will be computed via a global
     *  communication.  Otherwise, it must be equal to the sum of the
     *  local elements across all nodes. This will only be verified if
     *  Trilinos' Teuchos package was built with debug support (CMake
     *  Boolean option TEUCHOS_ENABLE_DEBUG=ON).  If verification
     *  fails, a std::invalid_argument exception will be thrown.
     */
    Map (global_size_t numGlobalElements, 
	 size_t numLocalElements, 
	 GlobalOrdinal indexBase, 
	 const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
	 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

    /** \brief Map constructor with user-defined non-contiguous (arbitrary) distribution.
     *  
     *  If numGlobalElements ==
     *  Teuchos::OrdinalTraits<global_size_t>::invalid(), the number
     *  of global elements will be computed via a global
     *  communication.  Otherwise, it must be equal to the sum of the
     *  local elements across all nodes. This will only be verified if
     *  Trilinos' Teuchos package was built with debug support (CMake
     *  Boolean option TEUCHOS_ENABLE_DEBUG=ON).  If verification
     *  fails, a std::invalid_argument exception will be thrown.
     */
    Map (global_size_t numGlobalElements, 
	 const Teuchos::ArrayView<const GlobalOrdinal> &elementList, 
	 GlobalOrdinal indexBase, 
	 const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
	 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

    //! Map destructor. 
    ~Map();

    //@}


    //! @name Map Attribute Methods
    //@{ 

    //! The number of elements in this Map.
    inline global_size_t getGlobalNumElements() const { return numGlobalElements_; }

    //! The number of elements belonging to the calling node.
    inline size_t getNodeNumElements() const { return numLocalElements_; }

    //! The index base for this Map.
    inline GlobalOrdinal getIndexBase() const { return indexBase_; }

    //! The minimum local index.
    inline LocalOrdinal getMinLocalIndex() const { return Teuchos::OrdinalTraits<LocalOrdinal>::zero(); }

    //! The maximum local index.
    inline LocalOrdinal getMaxLocalIndex() const { return Teuchos::as<LocalOrdinal>(numLocalElements_-1); }

    //! The minimum global index owned by this node.
    inline GlobalOrdinal getMinGlobalIndex() const { return minMyGID_; }

    //! The maximum global index owned by this node
    inline GlobalOrdinal getMaxGlobalIndex() const { return maxMyGID_; }

    //! The minimum global index over all nodes.
    inline GlobalOrdinal getMinAllGlobalIndex() const { return minAllGID_; }

    //! The maximum global index over all nodes.
    inline GlobalOrdinal getMaxAllGlobalIndex() const { return maxAllGID_; }

    /// \brief The local index corresponding to the given global index.  
    ///
    /// If the given global index is not owned by this node, return
    /// <tt>Teuchos::OrdinalTraits<LocalOrdinal>::invalid()</tt>.
    LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const;

    /// \brief The global index corresponding to the given local index.
    ///
    /// If the given local index is not valid for this node, return
    /// <tt>Teuchos::OrdinalTraits<GlobalOrdinal>::invalid()</tt>.
    GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const;

    /// \brief Return the node IDs and corresponding local IDs for a given list of global IDs.
    ///
    /// \pre nodeIDList.size() == GIDList.size()
    /// \pre LIDList.size() == GIDList.size()
    ///
    /// \return IDNotPresent indicates that at least one global ID was
    ///   not present in the directory.  Otherwise, return
    ///   AllIDsPresent.
    ///
    /// \note For a distributed noncontiguous Map, this operation
    ///   requires communication.  This is crucial technology used in
    ///   \c Export, \c Import, \c CrsGraph, and \c CrsMatrix.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                    const Teuchos::ArrayView<                int> & nodeIDList, 
                                    const Teuchos::ArrayView<       LocalOrdinal> & LIDList) const;

    /// \brief Return the node IDs for a given list of global IDs.
    ///
    /// \pre nodeIDList.size() == GIDList.size()
    /// \pre nodeIDList.size() == GIDList.size()
    ///
    /// \return IDNotPresent indicates that at least one global ID was
    ///   not present in the directory.  Otherwise, return
    ///   AllIDsPresent.
    ///
    /// \note For a distributed noncontiguous Map, this operation
    ///   requires communication.  This is crucial technology used in
    ///   \c Export, \c Import, \c CrsGraph, and \c CrsMatrix.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                    const Teuchos::ArrayView<                int> & nodeIDList) const;

    //! Return a view of the global indices owned by this node.
    Teuchos::ArrayView<const GlobalOrdinal> getNodeElementList() const;

    //! True if the local index is valid for this Map on this node, else false.
    bool isNodeLocalElement(LocalOrdinal localIndex) const;

    //! True if the global index is found in this Map on this node, else false.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const;

    /// \brief True if this Map is distributed contiguously, else false.
    ///
    /// Currently, creating this Map using the constructor for a
    /// user-defined arbitrary distribution (that takes a list of
    /// global elements owned on each process) means that this method
    /// always returns false.  We currently make no effort to test
    /// whether the user-provided global indices are actually
    /// contiguous on all the processes.  Many operations may be
    /// faster for contiguous Maps.  Thus, if you know the indices are
    /// contiguous on all processes, you should consider using one of
    /// the constructors for contiguous elements.
    bool isContiguous() const;

    /// \brief Whether this Map is globally distributed or locally replicated.
    ///
    /// \return True if this Map is globally distributed, else false.
    ///
    /// "Globally distributed" means that <i>all</i> of the following
    /// are true:
    ///
    /// 1. The map's communicator has more than one process.
    ///
    /// 2. There is at least one process in the map's communicator,
    ///    whose local number of elements does not equal the number of
    ///    global elements.  (That is, not all the elements are
    ///    replicated over all the processes.)
    ///
    /// If at least one of the above are not true, then the map is
    /// "locally replicated."  (The two are mutually exclusive.)
    ///
    /// Calling this method requires no communication or computation,
    /// because the result is precomputed in Map's constructors.
    bool isDistributed() const;

    //@}

    //! @name Boolean Tests
    //@{ 

    /// \brief True if and only if \c map is compatible with this Map.
    ///
    /// Two Maps are "compatible" if all of the following are true:
    /// 1. They have the same global number of elements.
    /// 2. They have the same number of local elements on each process.
    ///
    /// Determining #2 requires communication (a reduction over this
    /// Map's communicator).  This method assumes that the input Map
    /// is valid on all processes in this Map's communicator.
    ///
    /// Compatibility is useful for determining correctness of certain
    /// operations, like assigning one MultiVector X to another Y.  If
    /// X and Y have the same number of columns, and if their Maps are
    /// compatible, then it is legal to assign X to Y or to assign Y
    /// to X.
    ///
    /// Notes for Tpetra developers:
    ///
    /// If the input Map and this Map have different communicators,
    /// the behavior of this method is currently undefined.  In
    /// general, Tpetra currently assumes that if users instigate
    /// interactions between Tpetra objects, then those Tpetra objects
    /// have the same communicator.  Also, defining semantics of
    /// interaction between Tpetra objects with different
    /// communicators may be tricky.  It seems like two Maps could be
    /// compatible even if they had different communicators, as long
    /// as their communicators have the same number of processes.
    /// Could two Maps with different communicators be the same (in
    /// the sense of \c isSameAs())?  It's not clear.
    /// 
    /// Checking whether two communicators are the same would require
    /// extending Teuchos::Comm to provide a comparison operator.
    /// This could be implemented for MPI communicators by returning
    /// \c true if \c MPI_Comm_compare() returns \c MPI_IDENT, and \c
    /// false otherwise.  (Presumably, \c MPI_Comm_compare() works
    /// even if the two communicators have different process counts;
    /// the MPI 2.2 standard doesn't say otherwise.)  All serial
    /// communicators have the same context and contain the same
    /// number of processes, so all serial communicators are equal.
    bool isCompatible (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const;

    /// \brief True if and only if \c map is identical to this Map.
    ///
    /// "Identical" is stronger than "compatible."  Two Maps are
    /// identical if all of the following are true:
    /// 1. They have the same min and max global indices.
    /// 2. They have the same global number of elements.
    /// 3. They are either both distributed, or both not distributed.
    /// 4. Their index bases are the same.
    /// 5. They have the same number of local elements on each process.
    /// 6. They have the same global indices on each process.
    ///
    /// #2 and #5 are exactly "compatibility" (see \c isCompatible()).
    /// Thus, "identical" includes, but is stronger than,
    /// "compatible."  
    ///
    /// A Map corresponds to a block permutation over process ranks
    /// and global element indices.  Two Maps with different numbers
    /// of processes in their communicators cannot be compatible, let
    /// alone identical.  Two identical Maps correspond to the same
    /// permutation.
    ///
    /// Notes for Tpetra developers:
    ///
    /// If the input Map and this Map have different communicators,
    /// the behavior of this method is currently undefined.  See
    /// further notes on \c isCompatible().
    bool isSameAs (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const;

    //@}

    //@{ Accessors for the \c Teuchos::Comm and Kokkos Node objects.

    //! Get the Comm object for this Map
    const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

    //! Get the Node object for this Map
    const Teuchos::RCP<Node> & getNode() const;

    //@}

    //@{ Implementation of \c Teuchos::Describable 

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print this object with the given verbosity level to the given \c FancyOStream object.
    void describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;

    //@}


  private:

    //! Create this Map's Directory, if it hasn't been created already.
    void setupDirectory();

    /// \brief Determine whether this map is globally distributed or locally replicated.
    /// 
    /// \return True if the map is globally distributed, else false.
    ///
    /// This operation requires communication (a single all-reduce).
    /// See the documentation of \c isDistributed() for definitions
    /// of "globally distributed" and "locally replicated."
    ///
    /// Map invokes this method in its constructors if necessary, to
    /// set the \c distributed_ flag (and thus the return value of \c
    /// isDistributed()).  Map doesn't need to call \c checkIsDist()
    /// when using the uniform contiguous constructor with
    /// lg=GloballyDistributed, since then checking the number of
    /// processes in the communicator suffices.
    bool checkIsDist() const;

    //! Copy constructor (declared but not defined; do not use).
    Map(const Map<LocalOrdinal,GlobalOrdinal,Node> & source);

    //! Assignment operator (declared but not defined; do not use).
    Map<LocalOrdinal,GlobalOrdinal,Node>& operator=(const Map<LocalOrdinal,GlobalOrdinal,Node> & source);

    //! The communicator over which this Map is distributed.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    /// \brief The Kokkos Node instance (for shared-memory parallelism).
    ///
    /// Map doesn't need node yet, but it likely will later. In the
    /// meantime, passing a Node to Map means that we don't have to
    /// pass a Node to downstream classes such as MultiVector, Vector,
    /// CrsGraph and CrsMatrix.
    Teuchos::RCP<Node> node_;

    //! The index base for global IDs in this Map.
    GlobalOrdinal indexBase_;
    //! The number of global IDs located in this Map across all nodes.
    global_size_t numGlobalElements_;
    //! The number of global IDs located in this Map on this node.
    size_t numLocalElements_;
    //! The minimum and maximum global IDs located in this Map on this node.
    GlobalOrdinal minMyGID_, maxMyGID_;
    //! The minimum and maximum global IDs located in this Map across all nodes.
    GlobalOrdinal minAllGID_, maxAllGID_;
    //! Whether the range of global indices are contiguous and ordered.
    bool contiguous_;
    //! Whether this map's global indices are non-identically distributed among different nodes.
    bool distributed_;

    /// \brief A mapping from local IDs to global IDs.
    ///
    /// By definition, this mapping is local; it only contains global
    /// IDs owned by this process.  This mapping is created in two
    /// cases:
    ///
    /// 1. It is always created for a noncontiguous Map, in the
    ///    noncontiguous version of the Map constructor.
    ///
    /// 2. In \c getNodeElementList(), on demand (if it wasn't created
    ///    before).  
    ///
    /// The potential for on-demand creation is why this member datum
    /// is declared "mutable".  Note that other methods, such as \c
    /// describe(), may invoke \c getNodeElementList().  
    mutable Teuchos::ArrayRCP<GlobalOrdinal> lgMap_;

    /// \typedef global_to_local_table_type
    /// \brief Type of the table that maps global IDs to local IDs.
    /// 
    /// The actual type depends on the Tpetra_ENABLE_UNORDERD_MAP
    /// configure-time option.  If the option was set, we use
    /// std::unordered_map (implemented as a hash table with linear
    /// chaining).  Otherwise, we use std::map (a sorted data
    /// structure which is typically implemented as a red-black tree).
#ifdef HAVE_TPETRA_UNORDERED_MAP
    typedef std::unordered_map<GlobalOrdinal, LocalOrdinal> global_to_local_table_type;
#else
    typedef std::map<GlobalOrdinal, LocalOrdinal> global_to_local_table_type;
#endif // HAVE_TPETRA_UNORDERED_MAP

    /// \brief A mapping from global IDs to local IDs.
    ///
    /// This is a local mapping.  \c Directory implements the global
    /// mapping for all global IDs (both remote and locally owned).
    /// This object corresponds roughly to Epetra_BlockMapData's
    /// LIDHash_ hash table (which also maps from global IDs to local
    /// IDs).
    ///
    /// This mapping is built only for a noncontiguous map, by the
    /// noncontiguous map constructor.  For noncontiguous maps, the \c
    /// getLocalElement() and \c isNodeGlobalElement() methods use
    /// this mapping.
    global_to_local_table_type glMap_;

    /// \brief A Directory for looking up nodes for this Map. 
    ///
    /// This directory is a nonowning RCP and is therefore not allowed
    /// to persist beyond the lifetime of this Map.  Never allow this
    /// pointer to escape the Map.  It must be a nonowning RCP since
    /// the directory in turn must hold an RCP to this Map; making
    /// this an owning RCP would cause a circular dependency which
    /// would break reference counting.
    Teuchos::RCP<Directory<LocalOrdinal,GlobalOrdinal,Node> > directory_;

  }; // Map class

  /** \brief Non-member constructor for a locally replicated Map with the default Kokkos Node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
  createLocalMap(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

  /** \brief Non-member constructor for a locally replicated Map with a specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node);

  /** \brief Non-member constructor for a uniformly distributed, contiguous Map with the default Kokkos Node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
  createUniformContigMap(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

  /** \brief Non-member constructor for a uniformly distributed, contiguous Map with a user-specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createUniformContigMapWithNode(global_size_t numElements,
                                 const Teuchos::RCP< const Teuchos::Comm< int > > &comm, 
				 const Teuchos::RCP< Node > &node);

  /** \brief Non-member constructor for a (potentially) non-uniformly distributed, contiguous Map with the default Kokkos Node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
  createContigMap (global_size_t numElements, 
		   size_t localNumElements, 
		   const Teuchos::RCP<const Teuchos::Comm<int> > &comm);

  /** \brief Non-member constructor for a (potentially) non-uniformly distributed, contiguous Map with a user-specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createContigMapWithNode (global_size_t numElements, 
			   size_t localNumElements, 
			   const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
			   const Teuchos::RCP<Node> &node);

  /** \brief Non-member constructor for a non-contiguous Map with the default Kokkos Node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
  createNonContigMap (const ArrayView<const GlobalOrdinal> &elementList,
		      const RCP<const Teuchos::Comm<int> > &comm);

  /** \brief Non-member constructor for a non-contiguous Map with a user-specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createNonContigMapWithNode (const ArrayView<const GlobalOrdinal> &elementList,
			      const RCP<const Teuchos::Comm<int> > &comm, 
			      const RCP<Node> &node);

  /** \brief Non-member constructor for a contiguous Map with user-defined weights and a user-specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createWeightedContigMapWithNode (int thisNodeWeight, 
				   global_size_t numElements, 
				   const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
				   const Teuchos::RCP<Node> &node);

} // Tpetra namespace

/** \brief  Returns true if \c map is identical to this map. Implemented in Tpetra::Map::isSameAs().
    \relatesalso Tpetra::Map */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator== (const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map1, 
		 const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map2)
{ return map1.isSameAs(map2); }

/** \brief Returns true if \c map is not identical to this map. Implemented in Tpetra::Map::isSameAs().
    \relatesalso Tpetra::Map */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator!= (const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map1, 
		 const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map2)
{ return !map1.isSameAs(map2); }

#endif // TPETRA_MAP_DECL_HPP

