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

/// \file Tpetra_Map_decl.hpp
/// \brief Declarations for the Tpetra::Map class and related nonmember constructors.
///
namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // Forward declaration of Directory.
  template <class LO, class GO, class N> class Directory;
#endif

  /// \class Map
  /// \brief Describes a parallel distribution of objects over processes.
  ///
  /// \tparam LocalOrdinal The type of local indices.  Should be an
  ///   integer, and generally should be signed.  A good model of \c
  ///   LocalOrdinal is \c int.  (In Epetra, this is always just \c
  ///   int.)
  ///
  /// \tparam GlobalOrdinal The type of global indices.  Should be an
  ///   integer, and generally should be signed.  Also,
  ///   <tt>sizeof(GlobalOrdinal)</tt> should be greater than to equal
  ///   to <tt>sizeof(LocalOrdinal)</tt>.  For example, if \c
  ///   LocalOrdinal is \c int, good models of \c GlobalOrdinal are \c
  ///   int, \c long, <tt>long long</tt> (if the configure-time option
  ///   Teuchos_ENABLE_LONG_LONG was set), or \c ptrdiff_t.
  ///
  /// \tparam Node A class implementing on-node shared-memory parallel
  ///   operations.  It must implement the
  ///   \ref kokkos_node_api "Kokkos Node API."
  ///   The default \c Node type should suffice for most users.
  ///   The actual default type depends on your Trilinos build options.
  ///
  /// This class describes a distribution of data elements over one or
  /// more processes in a communicator.  Each element has a global
  /// index (of type \c GlobalOrdinal) uniquely associated to it.
  /// Each global index in the Map is "owned" by one or more processes
  /// in the Map's communicator.  The user gets to decide what an
  /// "element" means; examples include a row or column of a sparse
  /// matrix (as in CrsMatrix), or a row of one or more vectors (as in
  /// MultiVector).
  ///
  /// \section Kokkos_Map_prereq Prerequisites
  ///
  /// Before reading the rest of this documentation, it helps to know
  /// something about the Teuchos memory management classes, in
  /// particular Teuchos::RCP, Teuchos::ArrayRCP, and
  /// Teuchos::ArrayView.  You should also know a little bit about MPI
  /// (the Message Passing Interface for distributed-memory
  /// programming).  You won't have to use MPI directly to use Map,
  /// but it helps to be familiar with the general idea of distributed
  /// storage of data over a communicator.
  ///
  /// \section Tpetra_Map_concepts Map concepts
  ///
  /// \subsection Tpetra_Map_local_vs_global Local and global indices
  ///
  /// The distinction between local and global indices and types might
  /// confuse new Tpetra users.  <i>Global</i> indices represent the
  /// elements of a distributed object (such as rows or columns of a
  /// CrsMatrix, or rows of a MultiVector) uniquely over the entire
  /// object, which may be distributed over multiple processes.
  /// <i>Local</i> indices are local to the process that owns them.
  /// If global index G is owned by process P, then there is a unique
  /// local index L on process P corresponding to G.  If the local
  /// index L is valid on process P, then there is a unique global
  /// index G owned by P corresponding to the pair (L, P).  However,
  /// multiple processes might own the same global index (an
  /// "overlapping Map"), so a global index G might correspond to
  /// multiple (L, P) pairs.  In summary, local indices on a process
  /// correspond to object elements (e.g., sparse matrix rows or
  /// columns) owned by that process.
  ///
  /// Tpetra differs from Epetra in that local and global indices may
  /// have different types.  In Epetra, local and global indices both
  /// have type \c int.  In Tpetra, you get to pick the type of each.
  /// For example, you can use a 64-bit integer \c GlobalOrdinal type
  /// to solve problems with more than \f$2^{31}\f$ unknowns, but a
  /// 32-bit integer \c LocalOrdinal type to save bandwidth in sparse
  /// matrix-vector multiply.
  ///
  /// \subsection Tpetra_Map_contig Contiguous or noncontiguous
  ///
  /// A Map is <i>contiguous</i> when each process' list of global IDs
  /// forms an interval and is strictly increasing, and the globally
  /// minimum global ID equals the index base.  Map optimizes for the
  /// contiguous case.  In particular, noncontiguous Maps require
  /// communication in order to figure out which process owns a
  /// particular global ID.  (This communication happens in
  /// getRemoteIndexList().)
  ///
  /// \subsection Tpetra_Map_dist_repl Globally distributed or locally replicated
  ///
  /// "Globally distributed" means that <i>all</i> of the following
  /// are true:
  ///
  /// 1. The map's communicator has more than one process.
  /// 2. There is at least one process in the map's communicator,
  ///    whose local number of elements does not equal the number of
  ///    global elements.  (That is, not all the elements are
  ///    replicated over all the processes.)
  ///
  /// If at least one of the above are not true, then the map is
  /// "locally replicated."  (The two are mutually exclusive.)
  ///
  /// Globally distributed objects are partitioned across multiple
  /// processes in a communicator.  Each process owns at least one
  /// element in the object's Map that is not owned by another
  /// process.  For locally replicated objects, each element in the
  /// object's Map is owned redundantly by all processes in the
  /// object's communicator.  Some algorithms use objects that are too
  /// small to be distributed across all processes.  The upper
  /// Hessenberg matrix in a GMRES iterative solve is a good example.
  /// In other cases, such as with block iterative methods, block dot
  /// product functions produce small dense matrices that are required
  /// by all images.  Replicated local objects handle these
  /// situations.
  template <class LocalOrdinal,
            class GlobalOrdinal = LocalOrdinal,
            class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Map : public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{
    
    //! The type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! The type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! The type of the Kokkos Node.
    typedef Node node_type;

    //! @name Constructors and destructor
    //@{

    /** \brief Constructor with Tpetra-defined contiguous uniform distribution.
     *
     * This constructor produces a Map with the given number of
     * elements distributed among processes in the given communicator,
     * so that the subsets of global elements are nonoverlapping,
     * contiguous, and as evenly distributed across the processes as
     * possible.
     *
     * \param numGlobalElements [in] Number of elements in the Map
     *   (over all processes)
     *
     * \param indexBase [in] The base of the global indices
     *   in the Map.
     *   This must be the same on every node in the comm.
     *    For this Map constructor, the index
     *   base will also be the smallest global ID in the Map.
     *   (If you don't know what this should be, use zero.)
     *
     * \param comm [in] Communicator over which to distribute the
     *   elements.
     *
     * \param node [in/out] Kokkos Node instance.  The type of this
     *   object must match the type of the Node template parameter of
     *   Map.  If Node is not the same as the default Node type
     *   Kokkos::DefaultNode::DefaultNodeType, you will need to
     *   provide a nondefault value.
     */
    Map (global_size_t numGlobalElements,
         GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
         LocalGlobal lg=GloballyDistributed,
         const Teuchos::RCP<Node> &node = Kokkos::Details::getNode<Node>());

    /** \brief Constructor with a user-defined contiguous distribution.
     *
     * If N is the sum of numLocalElements over all processes, then
     * this constructor produces a nonoverlapping Map distributed over
     * the processes in the given communicator, with numLocalElements
     * contiguous elements on the calling process.
     *
     * \param numGlobalElements [in] If numGlobalElements ==
     *   Teuchos::OrdinalTraits<global_size_t>::invalid(), the number
     *   of global elements will be computed (via a global
     *   communication) as the sum of the counts of local elements.
     *   Otherwise, it must equal the sum of the local elements over
     *   all processes.  This will only be checked if Trilinos'
     *   Teuchos package was built with debug support (CMake Boolean
     *   option TEUCHOS_ENABLE_DEBUG=ON).  If verification fails, the
     *   constructor will throw std::invalid_argument.
     *
     * \param numLocalElements [in] Number of elements that the
     *   calling process will own in the Map.
     *
     * \param indexBase [in] The base of the global indices
     *   in the Map.
     *   This must be the same on every node in the comm.
     *   For this Map constructor, the index
     *   base will also be the smallest global ID in the Map. If you
     *   don't know what this should be, use zero.
     *
     * \param comm [in] Communicator over which to distribute the
     *   elements.
     *
     * \param node [in/out] Kokkos Node instance.  The type of this
     *   object must match the type of the Node template parameter of
     *   Map.  If Node is not the same as the default Node type
     *   Kokkos::DefaultNode::DefaultNodeType, you will need to
     *   provide a nondefault value.
     */
    Map (global_size_t numGlobalElements,
         size_t numLocalElements,
         GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
         const Teuchos::RCP<Node> &node = Kokkos::Details::getNode<Node>());

    /** \brief Constructor with user-defined arbitrary (possibly noncontiguous) distribution.
     *
     * Call this constructor if you have an arbitrary list of global
     * IDs that you want each process in the given communicator to
     * own.  Those IDs need not be contiguous, and the sets of global
     * IDs on different processes may overlap.  This is the
     * constructor to use to make an overlapping distribution.
     *
     * \param numGlobalElements [in] If numGlobalElements ==
     *   Teuchos::OrdinalTraits<global_size_t>::invalid(), the number
     *   of global elements will be computed (via a global
     *   communication) as the sum of the counts of local elements.
     *   Otherwise, it must equal the sum of the local elements over
     *   all processes.  This will only be checked if Trilinos'
     *   Teuchos package was built with debug support (CMake Boolean
     *   option TEUCHOS_ENABLE_DEBUG=ON).  If verification fails, the
     *   constructor will throw std::invalid_argument.
     *
     * \param elementList [in] List of global IDs owned by the calling
     *   process.
     *
     * \param indexBase [in] The base of the global indices
     *   in the Map.
     *   This must be the same on every node in the comm.
     *   This must be less than all of the global IDs in \c elementList.
     *   (If you don't know what this should
     *   be, use zero.)
     *
     * \param comm [in] Communicator over which to distribute the
     *   elements.
     *
     * \param node [in/out] Kokkos Node instance.  The type of this
     *   object must match the type of the Node template parameter of
     *   Map.  If Node is not the same as the default Node type
     *   Kokkos::DefaultNode::DefaultNodeType, you will need to
     *   provide a nondefault value.
     */
    Map (global_size_t numGlobalElements,
         const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
         GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
         const Teuchos::RCP<Node> &node = Kokkos::Details::getNode<Node>());

    //! Destructor.
    ~Map();

    //@}
    //! @name Attributes
    //@{

    //! The number of elements in this Map.
    inline global_size_t getGlobalNumElements() const { return numGlobalElements_; }

    //! The number of elements belonging to the calling node.
    inline size_t getNodeNumElements() const { return numLocalElements_; }

    //! The index base for this Map.
    inline GlobalOrdinal getIndexBase() const { return indexBase_; }

    //! The minimum local index.
    inline LocalOrdinal getMinLocalIndex() const {
      return Teuchos::OrdinalTraits<LocalOrdinal>::zero();
    }

    /// \brief The maximum local index on the calling process.
    ///
    /// If getNodeNumElements() == 0, this returns
    /// Teuchos::OrdinalTraits<LO>::invalid().
    inline LocalOrdinal getMaxLocalIndex() const {
      if (getNodeNumElements () == 0) {
        return Teuchos::OrdinalTraits<LocalOrdinal>::invalid ();
      }
      else { // Local indices are always zero-based.
        return Teuchos::as<LocalOrdinal> (getNodeNumElements () - 1);
      }
    }

    //! The minimum global index owned by the calling process.
    inline GlobalOrdinal getMinGlobalIndex() const { return minMyGID_; }

    //! The maximum global index owned by the calling process.
    inline GlobalOrdinal getMaxGlobalIndex() const { return maxMyGID_; }

    //! The minimum global index over all processes in the communicator.
    inline GlobalOrdinal getMinAllGlobalIndex() const { return minAllGID_; }

    //! The maximum global index over all processes in the communicator.
    inline GlobalOrdinal getMaxAllGlobalIndex() const { return maxAllGID_; }

    /// \brief The local index corresponding to the given global index.
    ///
    /// If the given global index is not owned by this process, return
    /// Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    LocalOrdinal getLocalElement (GlobalOrdinal globalIndex) const;

    /// \brief The global index corresponding to the given local index.
    ///
    /// If the given local index is not valid on the calling process,
    /// return Teuchos::OrdinalTraits<GlobalOrdinal>::invalid().
    GlobalOrdinal getGlobalElement (LocalOrdinal localIndex) const;

    /// \brief Return the process IDs and corresponding local IDs for the given global IDs.
    ///
    /// This operation should always be called as a collective over
    /// all processes in the communicator.  For a distributed
    /// noncontiguous Map, this operation requires communication.
    ///
    /// \param GIDList [in] List of global IDs for which to find
    ///   process IDs and local IDs.  These global IDs need not be
    ///   owned by the calling process.  Indeed, they need not be
    ///   owned by any process.
    /// \param nodeIDList [out] List of process IDs corresponding to
    ///   the given global IDs.  If a global ID does not belong to any
    ///   process, the resulting process ID is -1.
    /// \param LIDList [out] List of local IDs (that is, the local ID
    ///   on the process that owns them) corresponding the given
    ///   global IDs.  If a global ID does not have a local ID, the
    ///   resulting local ID is
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    ///
    /// \pre nodeIDList.size() == GIDList.size()
    /// \pre LIDList.size() == GIDList.size()
    ///
    /// \return IDNotPresent indicates that for at least one global
    ///   ID, we could not find the corresponding process ID.
    ///   Otherwise, return AllIDsPresent.
    ///
    /// \note This is crucial technology used in \c Export, \c Import,
    ///   \c CrsGraph, and \c CrsMatrix.
    LookupStatus
    getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal>& GIDList,
                        const Teuchos::ArrayView<                int>& nodeIDList,
                        const Teuchos::ArrayView<       LocalOrdinal>& LIDList) const;

    /// \brief Return the process IDs for the given global IDs.
    ///
    /// This operation should always be called as a collective over
    /// all processes in the communicator.  For a distributed
    /// noncontiguous Map, this operation requires communication.
    ///
    /// \param GIDList [in] List of global IDs for which to find
    ///   process IDs and local IDs.  These global IDs need not be
    ///   owned by the calling process.  Indeed, they need not be
    ///   owned by any process.
    /// \param nodeIDList [out] List of process IDs corresponding to
    ///   the given global IDs.  If a global ID does not belong to any
    ///   process, the resulting process ID is -1.
    ///
    /// \pre nodeIDList.size() == GIDList.size()
    ///
    /// \return IDNotPresent indicates that for at least one global
    ///   ID, we could not find the corresponding process ID.
    ///   Otherwise, return AllIDsPresent.
    ///
    /// \note For a distributed noncontiguous Map, this operation
    ///   requires communication.  This is crucial technology used in
    ///   \c Export, \c Import, \c CrsGraph, and \c CrsMatrix.
    LookupStatus
    getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal> & GIDList,
                        const Teuchos::ArrayView<                int> & nodeIDList) const;

    //! Return a view of the global indices owned by this node.
    Teuchos::ArrayView<const GlobalOrdinal> getNodeElementList() const;

    //@}
    //! @name Boolean tests
    //@{

    //! True if the local index is valid for this Map on this node, else false.
    bool isNodeLocalElement (LocalOrdinal localIndex) const;

    //! True if the global index is found in this Map on this node, else false.
    bool isNodeGlobalElement (GlobalOrdinal globalIndex) const;

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
    bool isContiguous () const;

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
    bool isDistributed () const;

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
    //! Accessors for the \c Teuchos::Comm and Kokkos Node objects.
    //@{

    //! Get this Map's Comm object.
    const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

    //! Get this Map's Node object.
    const Teuchos::RCP<Node> & getNode() const;

    //@}
    //! Implementation of \c Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print this object with the given verbosity level to the given \c FancyOStream.
    void
    describe (Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

    //! Advanced methods
    //@{

    //! \brief Create a shallow copy of this map, templated on a different node type
    template <class Node2>
    RCP<const Map<LocalOrdinal, GlobalOrdinal, Node2> > clone(const RCP<Node2> &node2) const;

    //@}

  protected:

    template <class LO, class GO, class N> friend class Map;

    //! Empty constructor; used for post-construction initialization in clone()
    Map() {}

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
    Map<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const Map<LocalOrdinal,GlobalOrdinal,Node> & source);

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
    /// The actual type depends on the Tpetra_ENABLE_UNORDERED_MAP
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
    RCP<global_to_local_table_type> glMap_;

    /// \brief A Directory for looking up nodes for this Map.
    ///
    /// *** This directory is not allowed to persist beyond the lifetime of this Map ***
    ///
    /// Never allow this pointer to escape the Map.
    /// The directory must hold an RCP to this Map, which must be non-owning
    /// to prevent a circular dependency.
    /// Therefore, allowing the Directory to persist beyond this Map would result
    /// in a dangling RCP. We avoid this by not sharing the Directory.
    Teuchos::RCP<Directory<LocalOrdinal,GlobalOrdinal,Node> > directory_;

  }; // Map class

  /** \brief Non-member constructor for a locally replicated Map with the default Kokkos Node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal> >
  createLocalMap(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

  /** \brief Non-member constructor for a locally replicated Map with a specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node = Kokkos::Details::getNode<Node>());

  /** \brief Non-member constructor for a uniformly distributed, contiguous Map with the default Kokkos Node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal> >
  createUniformContigMap(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

  /** \brief Non-member constructor for a uniformly distributed, contiguous Map with a user-specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createUniformContigMapWithNode(global_size_t numElements,
                                 const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                                 const Teuchos::RCP< Node > &node = Kokkos::Details::getNode<Node>());

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

  /** \brief Creates a one-to-one version of the given Map where each GID is owned by only one process.

      The user must guarantee there are no duplicate GID on the same processor. Unexepected behavior may result.

      \relatesalso Map
   */
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createOneToOne(Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &M);

} // Tpetra namespace

#include "Tpetra_Directory_decl.hpp"

namespace Tpetra {
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  template <class Node2>
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node2> >
  Map<LocalOrdinal,GlobalOrdinal,Node>::clone(const RCP<Node2> &node2) const
  {
    typedef Map<LocalOrdinal,GlobalOrdinal,Node2> Map2;
    RCP<Map2> map = rcp(new Map2());
    // the same old stuff...
    map->comm_              = comm_;
    map->indexBase_         = indexBase_;
    map->numGlobalElements_ = numGlobalElements_;
    map->numLocalElements_  = numLocalElements_;
    map->minMyGID_          = minMyGID_;
    map->maxMyGID_          = maxMyGID_;
    map->minAllGID_         = minAllGID_;
    map->maxAllGID_         = maxAllGID_;
    map->contiguous_        = contiguous_;
    map->distributed_       = distributed_;
    map->lgMap_             = lgMap_;
    map->glMap_             = glMap_;
    // the hot new stuff!
    map->node_              = node2;
    if (directory_ != null) {
      map->directory_ = directory_->template clone<Node2>(map.create_weak());
    }
    return map;
  }
}

/// \brief True if map1 is the same as (in the sense of isSameAs()) map2, else false.
/// \relatesalso Tpetra::Map
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator== (const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map1,
                 const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map2)
{ return map1.isSameAs(map2); }

/// \brief True if map1 is not the same as (in the sense of isSameAs()) map2, else false.
/// \relatesalso Tpetra::Map
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator!= (const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map1,
                 const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map2)
{ return !map1.isSameAs(map2); }

#endif // TPETRA_MAP_DECL_HPP

