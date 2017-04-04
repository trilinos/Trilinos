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

/// \file Tpetra_Map_decl.hpp
/// \brief Declarations for the Tpetra::Map class and related
///   nonmember constructors.

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Details_LocalMap.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DualView.hpp"
#include "Teuchos_Describable.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // Forward declaration of Directory.
  template <class LO, class GO, class N> class Directory;
#endif // DOXYGEN_SHOULD_SKIP_THIS

  namespace Details {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    // Forward declaration of TieBreak
    template <class LO, class GO> class TieBreak;
#endif // DOXYGEN_SHOULD_SKIP_THIS

    /// \class MapCloner
    /// \brief Implementation detail of Map::clone().
    template<class OutMapType, class InMapType>
    struct MapCloner {
      typedef typename OutMapType::node_type out_node_type;
      typedef typename InMapType::node_type in_node_type;

      static OutMapType
      clone (const InMapType& mapIn,
             const Teuchos::RCP<out_node_type>& node2);
    };

  } // namespace Details

  template<class Node>
  Teuchos::RCP<Node> defaultArgNode() {
    // Workaround function for a deferred visual studio bug
    // http://connect.microsoft.com/VisualStudio/feedback/details/719847/erroneous-error-c2783-could-not-deduce-template-argument
    // Use this function for default arguments rather than calling
    // what is the return value below.  Also helps in reducing
    // duplication in various constructors.
    return KokkosClassic::Details::getNode<Node> ();
  }

  /// \class Map
  /// \brief Describes a parallel distribution of objects over processes.
  ///
  /// \tparam LocalOrdinal The type of local indices.  This
  ///   <i>must</i> be a built-in integer type, and <i>must</i> be
  ///   signed.  A good model of <tt>LocalOrdinal</tt> is
  ///   <tt>int</tt>, which is also the default.  (In Epetra, this is
  ///   always just <tt>int</tt>.)
  ///
  /// \tparam GlobalOrdinal The type of global indices.  This
  ///   <i>must</i> be a built-in integer type.  We allow either
  ///   signed or unsigned types here, but prefer signed types.  Also,
  ///   we require that <tt>GlobalOrdinal</tt> be no smaller than
  ///   <tt>LocalOrdinal</tt>, that is:
  ///   \code
  ///   sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal);
  ///   \endcode
  ///   If <tt>LocalOrdinal</tt> is <tt>int</tt>, good models of
  ///   <tt>GlobalOrdinal</tt> are
  ///   <ul>
  ///   <li> \c int </li>
  ///   <li> \c long </li>
  ///   <li> <tt>long long</tt> (if the configure-time option
  ///        <tt>Teuchos_ENABLE_LONG_LONG_INT</tt> was set) </li>
  ///   <li> \c ptrdiff_t (a signed integer type which is 64 bits on
  ///        64-bit machines) </li>
  ///   </ul>
  ///   If the default <tt>GlobalOrdinal</tt> is <tt>int</tt>, then
  ///   the <i>global</i> number of rows or columns in the matrix may
  ///   be no more than \c INT_MAX, which for typical 32-bit \c int is
  ///   \f$2^{31} - 1\f$ (about two billion).  If you want to solve
  ///   larger problems, you must use a 64-bit integer type here.
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
  /// \section Tpetra_Map_prereq Prerequisites
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
  /// A <i>contiguous</i> Map divides an interval of global indices
  /// over the processes in its communicator, such that each process
  /// gets a contiguous interval of zero or more of those global
  /// indices, with the indices owned by a process p strictly greater
  /// than those owned by process q if \f$p > q\f$.  Formally, we call
  /// a Map contiguous when all of the following hold:
  /// <ol>
  /// <li>the set of global indices (over all processes) forms an
  ///   interval, </li>
  /// <li>every global index in that interval is owned by exactly one
  ///   process in the Map's communicator, </li>
  /// <li>the (ordered) list of global indices on each process p in
  ///   the Map's communicator forms a contiguous interval, and </li>
  /// <li>if process p owns a global index \f$g_p\f$ and process q
  ///   owns a global index \f$g_q\f$, and if \f$p > q\f$, then
  ///   \f$g_p > g_q\f$. </li>
  /// </ol>
  /// Different processes may own different numbers of global indices.
  /// We call a Map <i>uniform</i> if it is contiguous, <i>and</i> if
  /// the user let the Map divide a global count of indices evenly
  /// over the Map's communicator's processes.  The latter happens by
  /// calling the version of Map's constructor that takes a global
  /// count of indices, rather than a local count or an arbitrary list
  /// of indices.
  ///
  /// Map optimizes for the contiguous case.  For example,
  /// noncontiguous Maps always require communication in order to
  /// figure out which process owns a particular global index.  (This
  /// communication happens in getRemoteIndexList().)  Contiguous but
  /// nonuniform Maps may also require communication in this case,
  /// though we may only need to perform that communication once (at
  /// Map setup time).  Contiguous Maps also can convert between
  /// global and local indices more efficiently.
  ///
  /// \subsection Tpetra_Map_dist_repl Globally distributed or locally replicated
  ///
  /// <i>Globally distributed</i> means that <i>all</i> of the
  /// following are true:
  /// <ol>
  /// <li> The map's communicator has more than one process. </li>
  /// <li>There is at least one process in the map's communicator,
  ///     whose local number of elements does not equal the number of
  ///     global elements.  (That is, not all the elements are
  ///     replicated over all the processes.) </li>
  /// </ol>
  /// If at least one of the above are not true, then the map is
  /// <i>locally replicated.</i> (The two are mutually exclusive.)
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
  template <class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
            class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
            class Node = ::Tpetra::Details::DefaultTypes::node_type>
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

    //! The Kokkos execution space.
    typedef typename Node::execution_space execution_space;
    //! The Kokkos memory space.
    typedef typename Node::memory_space memory_space;

    /// \brief The Kokkos device type over which to allocate Views and
    ///   perform work.
    ///
    /// A Kokkos::Device is an (execution_space, memory_space) pair.
    /// It defines where the Map's data live, and where Map might
    /// choose to execute parallel kernels.
    typedef typename Node::device_type device_type;

    /// \brief Type of the "local" Map.
    ///
    /// \warning This object's interface is not yet fixed.  We provide
    ///   this object currently only as a service to advanced users.
    ///
    /// The "local" Map is suitable for use in Kokkos parallel
    /// operations in the Map's native execution space, which is
    /// <tt>device_type::execution_space</tt>.
    ///
    /// By "local," we mean that the object performs no MPI
    /// communication, and can only access information that would
    /// never need MPI communication, no matter what kind of Map this
    /// is.
    typedef Details::LocalMap<LocalOrdinal, GlobalOrdinal, device_type>
      local_map_type;

    //@}
    //! @name Constructors and destructor
    //@{

    /** \brief Constructor with contiguous uniform distribution.
     *
     * Build a Map representing the following contiguous range of
     * <tt>numGlobalElements</tt> indices:
     * \code
     * [indexBase,
     *  indexBase + 1, ...,
     *  numGlobalElements + indexBase - 1]
     * \endcode
     * For example, if \c indexBase is 0 and \c numGlobalElements is
     * N and positive, the resulting contiguous range is [0, N-1].
     *
     * The \c lg argument determines whether the indices will be
     * distributed evenly over all the processes in the given
     * communicator \c comm, or replicated on all processes in the
     * communicator.  "Distributed evenly" (the default) means that
     * each process gets a contiguous range of either
     * numGlobalElements / P or (numGlobalElements / P) + 1 indices.
     * The resulting Map is nonoverlapping.  "Replicated" means that
     * every process shares the range [0, N-1]; the resulting Map is
     * an overlapping Map.
     *
     * This constructor must be called as a collective over the input
     * communicator.  It reserves the right to use MPI collectives to
     * check input values in a debug build.  If it does check and any
     * check fails, it will throw std::invalid_argument on all
     * processes in the given communicator.
     *
     * \param numGlobalElements [in] Global number of indices in the
     *   Map (over all processes).
     *
     * \param indexBase [in] The base of the global indices in the
     *   Map.  This must be the same on every process in the
     *   communicator.  The index base will also be the smallest
     *   global index in the Map.  (If you don't know what this should
     *   be, use zero.)
     *
     * \param lg [in] Either <tt>GloballyDistributed</tt> or
     *   <tt>LocallyReplicated</tt>.  If <tt>GloballyDistributed</tt>
     *   and the communicator contains P processes, then each process
     *   will own either <tt>numGlobalElements/P</tt> or
     *   <tt>numGlobalElements/P + 1</tt> nonoverlapping contiguous
     *   indices.  If <tt>LocallyReplicated</tt>, then all processes
     *   will get the same set of indices, namely <tt>indexBase,
     *   indexBase + 1, ..., numGlobalElements + indexBase - 1</tt>.
     *
     * \param comm [in] Communicator over which to distribute the
     *   indices.
     *
     * \param node [in/out] (OPTIONAL; default usually suffices)
     *   Kokkos Node instance.
     */
    Map (global_size_t numGlobalElements,
         GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
         LocalGlobal lg=GloballyDistributed,
         const Teuchos::RCP<Node> &node = defaultArgNode<Node>());

    /** \brief Constructor with contiguous, possibly nonuniform
     *    distribution.
     *
     * If N is the sum of \c numLocalElements over all processes, then
     * this constructor produces a nonoverlapping Map with N indices
     * in the global contiguous range [0, N-1], distributed over all
     * the processes in the given communicator \c comm, with a
     * contiguous range of \c numLocalElements indices on the calling
     * process.
     *
     * This constructor must be called as a collective over the input
     * communicator.  It reserves the right to use MPI collectives to
     * check input values in a debug build.  If it does check and any
     * check fails, it will throw std::invalid_argument on all
     * processes in the given communicator.
     *
     * \param numGlobalElements [in] If you want Tpetra to compute the
     *   global number of indices in the Map, set this to
     *   Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid().
     *   This costs a global all-reduce.  Otherwise, this must equal
     *   the sum of numLocalElements over all processes in the input
     *   communicator \c comm.
     *
     * \param numLocalElements [in] Number of indices that the calling
     *   process will own in the Map.
     *
     * \param indexBase [in] The base of the global indices in the
     *   Map.  This must be the same on every process in the given
     *   communicator.  For this Map constructor, the index base will
     *   also be the smallest global index in the Map.  If you don't
     *   know what this should be, use zero.
     *
     * \param comm [in] Communicator over which to distribute the
     *   elements.
     *
     * \param node [in/out] (OPTIONAL; default usually suffices)
     *   Kokkos Node instance.
     */
    Map (global_size_t numGlobalElements,
         size_t numLocalElements,
         GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
         const Teuchos::RCP<Node> &node = defaultArgNode<Node>());

    /** \brief Constructor with arbitrary (possibly noncontiguous
     *   and/or nonuniform and/or overlapping) distribution, taking
     *   the input global indices as a Kokkos::View.
     *
     * Call this constructor if you have an arbitrary list of global
     * indices for each process in the given communicator.  Those
     * indices need not be contiguous, and the sets of global indices
     * on different processes may overlap.  This is one of the
     * constructors to use to make a general, possibly overlapping
     * distribution.
     *
     * This constructor, like all Map constructors, must be called as
     * a collective over the input communicator.  It reserves the
     * right to use MPI collectives to check input values in a debug
     * build.  If it does check and any check fails, it will throw
     * std::invalid_argument on all processes in the given
     * communicator.
     *
     * \param numGlobalElements [in] If <tt>numGlobalElements ==
     *   Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()</tt>,
     *   the number of global elements will be computed (via a global
     *   communication) as the sum of the counts of local indices.
     *   Otherwise, it must equal the sum of the number of indices on
     *   each process, over all processes in the given communicator,
     *   and must be the same on all processes in the communicator.
     *
     * \param indexList [in] List of global indices owned by the
     *   calling process.  (This likely differs on different
     *   processes.)
     *
     * \param indexBase [in] The base of the global indices in the
     *   Map.  This must be the same on every process in the given
     *   communicator \c comm.  Currently, Map requires that this
     *   equal the global minimum index over all processes'
     *   <tt>entryList</tt> inputs.
     *
     * \param comm [in] Communicator over which to distribute the
     *   indices.  This constructor must be called as a collective
     *   over this communicator.
     */
    Map (const global_size_t numGlobalElements,
         const Kokkos::View<const GlobalOrdinal*, device_type>& indexList,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

    /** \brief Constructor with arbitrary (possibly noncontiguous
     *   and/or nonuniform and/or overlapping) distribution, taking
     *   the input global indices as a raw host pointer.
     *
     * Call this constructor if you have an arbitrary list of global
     * indices for each process in the given communicator.  Those
     * indices need not be contiguous, and the sets of global indices
     * on different processes may overlap.  This is one of the
     * constructors to use to make a general, possibly overlapping
     * distribution.
     *
     * This constructor, like all Map constructors, must be called as
     * a collective over the input communicator.  It reserves the
     * right to use MPI collectives to check input values in a debug
     * build.  If it does check and any check fails, it will throw
     * std::invalid_argument on all processes in the given
     * communicator.
     *
     * \param numGlobalElements [in] If <tt>numGlobalElements ==
     *   Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()</tt>,
     *   the number of global elements will be computed (via a global
     *   communication) as the sum of the counts of local elements.
     *   Otherwise, it must equal the sum of the local elements over
     *   all processes.  This value must be the same on all processes
     *   participating in the call.
     *
     * \param indexList [in] List of global indices owned by the
     *   calling process.
     *
     * \param indexListSize [in] Number of valid entries in indexList.
     *   This is a LocalOrdinal because the number of indices owned by
     *   each process must fit in LocalOrdinal.
     *
     * \param indexBase [in] The base of the global indices in the
     *   Map.  This must be the same on every process in the given
     *   communicator.  Currently, Map requires that this equal the
     *   global minimum index over all processes' \c indexList inputs.
     *
     * \param comm [in] Communicator over which to distribute the
     *   elements.
     */
    Map (const global_size_t numGlobalElements,
         const GlobalOrdinal indexList[],
         const LocalOrdinal indexListSize,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

    /** \brief Constructor with arbitrary (possibly noncontiguous
     *   and/or nonuniform and/or overlapping) distribution, taking
     *   the input global indices as a Teuchos::ArrayView (for
     *   backwards compatibility).
     *
     * Call this constructor if you have an arbitrary list of global
     * indices for each process in the given communicator.  Those
     * indices need not be contiguous, and the sets of global indices
     * on different processes may overlap.  This is one of the
     * constructors to use to make a general, possibly overlapping
     * distribution.
     *
     * This constructor, like all Map constructors, must be called as
     * a collective over the input communicator.  It reserves the
     * right to use MPI collectives to check input values in a debug
     * build.  If it does check and any check fails, it will throw
     * std::invalid_argument on all processes in the given
     * communicator.
     *
     * \param numGlobalElements [in] If <tt>numGlobalElements ==
     *   Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()</tt>,
     *   the number of global elements will be computed (via a global
     *   communication) as the sum of the counts of local elements.
     *   Otherwise, it must equal the sum of the number of indices on
     *   each process, over all processes in the given communicator,
     *   and must be the same on all processes in the communicator.
     *
     * \param indexList [in] List of global indices owned by the
     *   calling process.  (This likely differs on different
     *   processes.)
     *
     * \param indexBase [in] The base of the global indices in the
     *   Map.  This must be the same on every process in the given
     *   communicator.  Currently, Map requires that this equal the
     *   global minimum index over all processes' <tt>indexList</tt>
     *   inputs.
     *
     * \param comm [in] Communicator over which to distribute the
     *   indices.  This constructor must be called as a collective
     *   over this communicator.
     *
     * \param node [in/out] (OPTIONAL; default usually suffices)
     *   Kokkos Node instance.
     */
    Map (const global_size_t numGlobalElements,
         const Teuchos::ArrayView<const GlobalOrdinal>& indexList,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<Node>& node = defaultArgNode<Node>());

    /// \brief Default constructor (that does nothing).
    ///
    /// This creates an empty Map, with 0 (zero) indices total.  The
    /// Map's communicator only includes the calling process; in MPI
    /// terms, it behaves like MPI_COMM_SELF.
    ///
    /// This constructor exists mainly to support view semantics of
    /// Map.  That is, we can create an empty Map, and then assign a
    /// nonempty Map to it using operator=.  This constructor is also
    /// useful in methods like clone() and removeEmptyProcesses(),
    /// where we have the information to initialize the Map more
    /// efficiently ourselves, without going through one of the three
    /// usual Map construction paths.
    Map ();

    //! Destructor.
    ~Map ();

    //@}
    //! @name Attributes
    //@{

    /// \brief Whether the Map is one to one.
    ///
    /// This must be called collectively over all processes in the
    /// Map's communicator.
    bool isOneToOne () const;

    /// \brief The number of elements in this Map.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    global_size_t getGlobalNumElements () const {
      return numGlobalElements_;
    }

    /// \brief The number of elements belonging to the calling process.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    size_t getNodeNumElements () const {
      return numLocalElements_;
    }

    /// \brief The index base for this Map.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    GlobalOrdinal getIndexBase () const {
      return indexBase_;
    }

    /// \brief The minimum local index.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    LocalOrdinal getMinLocalIndex () const {
      return static_cast<LocalOrdinal> (0);
    }

    /// \brief The maximum local index on the calling process.
    ///
    /// If this process owns no elements, that is, if
    /// <tt>getNodeNumElements() == 0</tt>, then this method returns
    /// the same value as
    /// <tt>Teuchos::OrdinalTraits<LocalOrdinal>::invalid()</tt>.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    LocalOrdinal getMaxLocalIndex () const {
      if (this->getNodeNumElements () == 0) {
        return Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid ();
      } else { // Local indices are always zero-based.
        return static_cast<LocalOrdinal> (this->getNodeNumElements () - 1);
      }
    }

    /// \brief The minimum global index owned by the calling process.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    GlobalOrdinal getMinGlobalIndex () const {
      return minMyGID_;
    }

    /// \brief The maximum global index owned by the calling process.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    GlobalOrdinal getMaxGlobalIndex () const {
      return maxMyGID_;
    }

    /// \brief The minimum global index over all processes in the communicator.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    GlobalOrdinal getMinAllGlobalIndex () const {
      return minAllGID_;
    }

    /// \brief The maximum global index over all processes in the communicator.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    GlobalOrdinal getMaxAllGlobalIndex () const {
      return maxAllGID_;
    }

    /// \brief The local index corresponding to the given global index.
    ///
    /// \param globalIndex [in] The global index.
    ///
    /// \return If the given global index is owned by the calling
    ///   process, return the corresponding local index, else return
    ///   the same value as
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    LocalOrdinal getLocalElement (GlobalOrdinal globalIndex) const;

    /// \brief The global index corresponding to the given local index.
    ///
    /// \param localIndex [in] The local index.
    ///
    /// \return If the given local index is valid on the calling
    ///   process, return the corresponding global index, else return
    ///   the same value as
    ///   Teuchos::OrdinalTraits<GlobalOrdinal>::invalid().
    GlobalOrdinal getGlobalElement (LocalOrdinal localIndex) const;

    /// \brief Get the local Map for Kokkos kernels.
    ///
    /// \warning The interface of the local Map object is SUBJECT TO
    ///   CHANGE and is for EXPERT USERS ONLY.
    local_map_type getLocalMap () const;

    /// \brief Return the process ranks and corresponding local
    ///   indices for the given global indices.
    ///
    /// This operation must always be called as a collective over all
    /// processes in the Map's communicator.  For a distributed
    /// noncontiguous Map, this operation requires communication.
    ///
    /// \param GIDList [in] List of global indices for which to find
    ///   process ranks and local indices.  These global indices need
    ///   not be owned by the calling process.  Indeed, they need not
    ///   be owned by any process.
    /// \param nodeIDList [out] List of process rank corresponding to
    ///   the given global indices.  If a global index does not belong
    ///   to any process, the resulting process rank is -1.
    /// \param LIDList [out] List of local indices (that is, the local
    ///   index on the process that owns them) corresponding to the
    ///   given global indices.  If a global index does not have a
    ///   local index, the resulting local index has the same value as
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    ///
    /// \pre nodeIDList.size() == GIDList.size()
    /// \pre LIDList.size() == GIDList.size()
    ///
    /// \return IDNotPresent indicates that for at least one global
    ///   index, we could not find the corresponding process rank.
    ///   Otherwise, return AllIDsPresent.
    ///
    /// \note This is crucial technology used in Export, Import,
    ///   CrsGraph, and CrsMatrix.
    LookupStatus
    getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal>& GIDList,
                        const Teuchos::ArrayView<                int>& nodeIDList,
                        const Teuchos::ArrayView<       LocalOrdinal>& LIDList) const;

    /// \brief Return the process ranks for the given global indices.
    ///
    /// This method must always be called as a collective over all
    /// processes in the Map's communicator.  For a distributed
    /// noncontiguous Map, this operation requires communication.
    ///
    /// \param GIDList [in] List of global indices for which to find
    ///   process ranks and local indices.  These global indices need
    ///   not be owned by the calling process.  Indeed, they need not
    ///   be owned by any process.
    /// \param nodeIDList [out] List of process ranks corresponding to
    ///   the given global indices.  If a global index does not belong
    ///   to any process, the resulting process rank is -1.
    ///
    /// \pre nodeIDList.size() == GIDList.size()
    ///
    /// \return IDNotPresent indicates that for at least one global
    ///   index, we could not find the corresponding process rank.
    ///   Otherwise, return AllIDsPresent.
    ///
    /// \note For a distributed noncontiguous Map, this operation
    ///   requires communication.  This is crucial technology used in
    ///   Export, Import, CrsGraph, and CrsMatrix.
    LookupStatus
    getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal> & GIDList,
                        const Teuchos::ArrayView<                int> & nodeIDList) const;

  private:
    /// \brief Type of lgMap_ (see below); used to derive return type
    ///   of getMyGlobalIndices() (also below).
    ///
    /// \warning YOU ARE NOT ALLOWED TO REFER TO THIS TYPE BY NAME.
    ///   Use <tt>auto</tt> to refer to the type of the return value
    ///   of getMyGlobalIndices().
    ///
    /// I would have preferred not to have this typedef at all.  It
    /// exists only so that we could avoid needing to declare lgMap_
    /// before declaring the getMyGlobalIndices() method.  That would
    /// have made this class declaration harder to read.
    typedef Kokkos::View<const GlobalOrdinal*,
                         Kokkos::LayoutLeft,
                         device_type> global_indices_array_type;

  public:
    /// \brief Return a view of the global indices owned by this process.
    ///
    /// The returned "view" has some type that looks like
    /// <ul>
    /// <li> <tt> Kokkos::View<const GlobalOrdinal*, ...> </tt> or </li>
    /// <li> <tt> Teuchos::ArrayView<const GlobalOrdinal> </tt> </li>
    /// </ul>
    /// It implements operator[] and the size() method, and behaves as
    /// a one-dimensional array.  You may <i>not</i> modify its
    /// entries.
    ///
    /// \warning You are NOT allowed to refer to the return value's
    ///   type by name.  That name is private.  Use <tt>auto</tt>
    ///   instead.
    ///
    /// If you call this method on a contiguous Map, it will create
    /// and cache the list of global indices for later use.  Beware of
    /// calling this if the calling process owns a very large number
    /// of global indices.
    global_indices_array_type getMyGlobalIndices () const;

    /// \brief Return a NONOWNING view of the global indices owned by
    ///   this process.
    ///
    /// \warning This method may be deprecated at some point.  Please
    ///   consider using getMyGlobalIndices() (see above) instead.
    ///
    /// If you call this method on a contiguous Map, it will create
    /// and cache the list of global indices for later use.  Beware of
    /// calling this if the calling process owns a very large number
    /// of global indices.
    Teuchos::ArrayView<const GlobalOrdinal> getNodeElementList() const;

    //@}
    //! @name Boolean tests
    //@{

    /// \brief Whether the given local index is valid for this Map on
    ///   the calling process.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    bool isNodeLocalElement (LocalOrdinal localIndex) const;

    /// \brief Whether the given global index is owned by this Map on
    ///   the calling process.
    ///
    /// \note This function should be thread safe and thread scalable,
    ///   assuming that you refer to the Map by value or reference,
    ///   not by Teuchos::RCP.
    bool isNodeGlobalElement (GlobalOrdinal globalIndex) const;

    /// \brief Whether the range of global indices is uniform.
    ///
    /// This is a conservative quantity.  It need only be true if the
    /// Map was constructed using the first (uniform contiguous)
    /// constructor or a nonmember constructor that calls it.  We
    /// reserve the right to do more work to check this in the future.
    bool isUniform () const;

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

    /// \brief Whether this Map is globally distributed or locally
    ///   replicated.
    ///
    /// \return True if this Map is globally distributed, else false.
    ///
    /// "Globally distributed" means that <i>all</i> of the following
    /// are true:
    /// <ol>
    /// <li> The map's communicator has more than one process.</li>
    /// <li> There is at least one process in the map's communicator,
    ///    whose local number of elements does not equal the number of
    ///    global elements.  (That is, not all the elements are
    ///    replicated over all the processes.)</li>
    /// </ol>
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
    /// <ol>
    /// <li> Their communicators have the same numbers of processes.
    ///    (This is necessary even to call this method.)</li>
    /// <li> They have the same global number of elements.</li>
    /// <li> They have the same number of local elements on each process.</li>
    /// </ol>
    ///
    /// Determining #3 requires communication (a reduction over this
    /// Map's communicator).  This method assumes that the input Map
    /// is valid on all processes in this Map's communicator.
    ///
    /// Compatibility is useful for determining correctness of certain
    /// operations, like assigning one MultiVector X to another Y.  If
    /// X and Y have the same number of columns, and if their Maps are
    /// compatible, then it is legal to assign X to Y or to assign Y
    /// to X.
    ///
    /// The behavior of this method is undefined if the input Map's
    /// communicator and this Map's communicator have different
    /// numbers of processes.  This method must be called collectively
    /// over this Map's communicator.
    bool isCompatible (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const;

    /// \brief True if and only if \c map is identical to this Map.
    ///
    /// "Identical" is stronger than "compatible."  Two Maps are
    /// identical if all of the following are true:
    /// <ol>
    /// <li> Their communicators are <i>congruent</i> (have the same
    ///    number of processes, in the same order: this corresponds to
    ///    the \c MPI_IDENT or \c MPI_CONGRUENT return values of
    ///    MPI_Comm_compare).</li>
    /// <li> They have the same min and max global indices.</li>
    /// <li> They have the same global number of elements.</li>
    /// <li> They are either both distributed, or both not distributed.</li>
    /// <li> Their index bases are the same.</li>
    /// <li> They have the same number of local elements on each process.</li>
    /// <li> They have the same global indices on each process.</li>
    /// </ol>
    ///
    /// "Identical" (isSameAs()) includes and is stronger than
    /// "compatible" (isCompatible()).
    ///
    /// A Map corresponds to a block permutation over process ranks
    /// and global element indices.  Two Maps with different numbers
    /// of processes in their communicators cannot be compatible, let
    /// alone identical.  Two identical Maps correspond to the same
    /// permutation.
    ///
    /// The behavior of this method is undefined if the input Map's
    /// communicator and this Map's communicator have different
    /// numbers of processes.  This method must be called collectively
    /// over this Map's communicator.
    bool isSameAs (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const;

    /// \brief Is this Map locally the same as the input Map?
    ///
    /// "Locally the same" means that on the calling process, the two
    /// Maps' global indices are the same and occur in the same order.
    bool locallySameAs (const Map<LocalOrdinal, GlobalOrdinal, node_type>& map) const;

    //@}
    //! Accessors for the Teuchos::Comm and Kokkos Node objects.
    //@{

    //! Get this Map's communicator, as a Teuchos::Comm.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

    //! Get this Map's Node object.
    Teuchos::RCP<Node> getNode () const;

    //@}
    //! Implementation of \c Teuchos::Describable
    //@{

    //! Return a one-line description of this object.
    std::string description () const;

    /// \brief Describe this object in a human-readable way to the
    ///   given output stream.
    ///
    /// You must call this method as a collective over all processes
    /// in this object's communicator.
    ///
    /// \param out [out] Output stream to which to write.  Only
    ///   Process 0 in this object's communicator may write to the
    ///   output stream.
    ///
    /// \param verbLevel [in] Verbosity level.  This also controls
    ///   whether this method does any communication.  At verbosity
    ///   levels higher (greater) than Teuchos::VERB_LOW, this method
    ///   may behave as a collective over the object's communicator.
    ///
    /// Teuchos::FancyOStream wraps std::ostream.  It adds features
    /// like tab levels.  If you just want to wrap std::cout, try
    /// this:
    /// \code
    /// auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::out));
    /// \endcode
    void
    describe (Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel verbLevel =
                Teuchos::Describable::verbLevel_default) const;
    //@}
    //! Advanced methods
    //@{

    //! Create a shallow copy of this Map, with a different Node type.
    template <class NodeOut>
    Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, NodeOut> >
    clone (const Teuchos::RCP<NodeOut>& nodeOut) const;

    /// \brief Return a new Map with processes with zero elements removed.
    ///
    /// \warning This method is only for expert users.  Understanding
    ///   how to use this method correctly requires some familiarity
    ///   with semantics of MPI communicators.
    ///
    /// \warning We make no promises of backwards compatibility for
    ///   this method.  It may go away or change at any time.
    ///
    /// This method first computes a new communicator, which contains
    /// only those processes in this Map's communicator (the "original
    /// communicator") that have a nonzero number of elements in this
    /// Map (the "original Map").  It then returns a new Map
    /// distributed over the new communicator.  The new Map represents
    /// the same distribution as the original Map, except that
    /// processes containing zero elements are not included in the new
    /// Map or its communicator.  On processes not included in the new
    /// Map or communicator, this method returns
    /// <tt>Teuchos::null</tt>.
    ///
    /// The returned Map always has a distinct communicator from this
    /// Map's original communicator.  The new communicator contains a
    /// subset of processes from the original communicator.  Even if
    /// the number of processes in the new communicator equals the
    /// number of processes in the original communicator, the new
    /// communicator is distinct.  (In an MPI implementation, the new
    /// communicator is created using MPI_Comm_split.)  Furthermore, a
    /// process may have a different rank in the new communicator, so
    /// be wary of classes that like to store the rank rather than
    /// asking the communicator for it each time.
    ///
    /// This method must be called collectively on the original
    /// communicator.  It leaves the original Map and communicator
    /// unchanged.
    ///
    /// This method was intended for applications such as algebraic
    /// multigrid or other multilevel preconditioners.  Construction
    /// of each level of the multilevel preconditioner typically
    /// requires constructing sparse matrices, which in turn requires
    /// all-reduces over all participating processes at that level.
    /// Matrix sizes at successively coarser levels shrink
    /// geometrically.  At the coarsest levels, some processes might
    /// be left with zero rows of the matrix, or the multigrid
    /// implementation might "rebalance" (redistribute the matrix) and
    /// intentionally leave some processes with zero rows.  Removing
    /// processes with zero rows makes the all-reduces and other
    /// communication operations cheaper.
    Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
    removeEmptyProcesses () const;

    /// \brief Replace this Map's communicator with a subset communicator.
    ///
    /// \warning This method is only for expert users.  Understanding
    ///   how to use this method correctly requires some familiarity
    ///   with semantics of MPI communicators.
    ///
    /// \warning We make no promises of backwards compatibility for
    ///   this method.  It may go away or change at any time.
    ///
    /// \pre The input communicator's processes are a subset of this
    ///   Map's current communicator's processes.
    /// \pre On processes which are not included in the input
    ///   communicator, the input communicator is null.
    ///
    /// This method must be called collectively on the original
    /// communicator.  It leaves the original Map and communicator
    /// unchanged.
    ///
    /// \note This method differs from removeEmptyProcesses(), in that
    ///   it does not assume that excluded processes have zero
    ///   entries.  For example, one might wish to remove empty
    ///   processes from the row Map of a CrsGraph using
    ///   removeEmptyProcesses(), and then apply the resulting subset
    ///   communicator to the column, domain, and range Maps of the
    ///   same graph.  For the latter three Maps, one would in general
    ///   use this method instead of removeEmptyProcesses(), giving
    ///   the new row Map's communicator to this method.
    Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
    replaceCommWithSubset (const Teuchos::RCP<const Teuchos::Comm<int> >& newComm) const;
    //@}

  protected:
    // This lets other specializations of Map access all of this
    // specialization's internal methods and data, so that we can
    // implement clone() without exposing the details of Map to users.
    template <class LO, class GO, class N> friend class Map;

  private:
    template<class OutMapType, class InMapType>
    friend struct Details::MapCloner;

    /// \brief Print the calling process' verbose describe()
    ///   information to the returned string.
    ///
    /// This is an implementation detail of describe().
    std::string
    localDescribeToString (const Teuchos::EVerbosityLevel vl) const;

    /// \brief Create this Map's Directory, if it hasn't been created already.
    ///
    /// This method must be called collectively over all processes in
    /// the Map's communicator.
    ///
    /// This is declared "const" so that we can call it in
    /// getRemoteIndexList() to create the Directory on demand.
    void setupDirectory () const;

    /// \brief Determine whether this map is globally distributed or locally replicated.
    ///
    /// \return True if the map is globally distributed, else false.
    ///
    /// This operation requires communication (a single all-reduce).
    /// See the documentation of isDistributed() for definitions of
    /// "globally distributed" and "locally replicated."
    ///
    /// Map invokes this method in its constructors if necessary, to
    /// set the \c distributed_ flag (and thus the return value of
    /// isDistributed()).  Map doesn't need to call checkIsDist() when
    /// using the uniform contiguous constructor with
    /// <tt>lg=GloballyDistributed</tt>, since then checking the
    /// number of processes in the communicator suffices.
    bool checkIsDist() const;

    /// \brief Call at the beginning of the nonuniform constructors;
    ///   it does checks (with extra global communication) in a debug
    ///   build.  In a release build, it does nothing.
    ///
    /// \return In a debug build: The global sum of numLocalElements
    ///   over all processes in the given communicator.  In a release
    ///   build: 0 (zero).
    global_size_t
    initialNonuniformDebugCheck (const global_size_t numGlobalElements,
                                 const size_t numLocalElements,
                                 const GlobalOrdinal indexBase,
                                 const Teuchos::RCP<const Teuchos::Comm<int> >& comm) const;

    void
    initWithNonownedHostIndexList (const global_size_t numGlobalElements,
                                   const Kokkos::View<const GlobalOrdinal*,
                                     Kokkos::LayoutLeft,
                                     Kokkos::HostSpace,
                                     Kokkos::MemoryUnmanaged>& entryList,
                                   const GlobalOrdinal indexBase,
                                   const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

    //! The communicator over which this Map is distributed.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    /// \brief LEGACY Node instance.
    ///
    /// This object only exists for backwards compatibility.  Please
    /// do not rely upon it.  In particular, do not assume that this
    /// pointer is non-null.
    Teuchos::RCP<Node> node_;

    //! The index base for global indices in this Map.
    GlobalOrdinal indexBase_;

    /// \brief The total number of global indices in this Map over all
    ///   processes in its communicator \c comm (see above).
    global_size_t numGlobalElements_;

    //! The number of global indices owned by this process.
    size_t numLocalElements_;

    //! The min global index owned by this process.
    GlobalOrdinal minMyGID_;

    //! The max global index owned by this process.
    GlobalOrdinal maxMyGID_;

    /// \brief The min global index in this Map over all processes in
    ///   its communicator \c comm (see above).
    GlobalOrdinal minAllGID_;

    /// \brief The max global index in this Map over all processes in
    ///   its communicator \c comm (see above).
    GlobalOrdinal maxAllGID_;

    /// \brief First contiguous GID.
    ///
    /// This is only set if the Map was created using the
    /// noncontiguous constructor.  In that case, if the calling
    /// process owns at least one GID, this will always equal that
    /// first GID in the list of GIDs given to the constructor.
    GlobalOrdinal firstContiguousGID_;

    /// \brief Last contiguous GID.
    ///
    /// This is only set if the Map was created using the
    /// noncontiguous constructor.  In that case, if the calling
    /// process owns at least one GID, this will always equal the last
    /// GID (inclusive) that forms an initial sequence of contiguous
    /// GIDs, in the list of GIDs given to the constructor.
    ///
    /// For example, if the list is [42, 43, 44, 45, 100, 1001],
    /// firstContiguousGID_ will be 42 and lastContiguousGID_ will be
    /// 45.  If the list is [42, 100, 1001, 1002, 1003],
    /// firstContiguousGID_ will be 42 and lastContiguousGID_ will
    /// also be 42.
    GlobalOrdinal lastContiguousGID_;

    /// \brief Whether the range of global indices is uniform.
    ///
    /// This is only true if the Map was constructed using the first
    /// (uniform contiguous) constructor or a nonmember constructor
    /// that calls it.
    bool uniform_;

    //! Whether the range of global indices are contiguous and ordered.
    bool contiguous_;

    /// \brief Whether this map's global indices are distributed
    ///   (true), or locally replicated (false), over its communicator
    ///   \c comm (see above).
    ///
    /// This is true if the Map is globally distributed, and false
    /// otherwise (if the Map is locally replicated).  See the
    /// documentation of isDistributed() for a definition of these two
    /// mutually exclusive terms.
    bool distributed_;

    /// \brief A mapping from local IDs to global IDs.
    ///
    /// By definition, this mapping is local; it only contains global
    /// IDs owned by this process.  This mapping is created in two
    /// cases:
    ///
    /// <ol>
    /// <li> It is always created for a noncontiguous Map, in the
    ///    noncontiguous version of the Map constructor.</li>
    /// <li> In getNodeElementList(), on demand (if it wasn't created
    ///    before).</li>
    /// </ol>
    ///
    /// The potential for on-demand creation is why this member datum
    /// is declared "mutable".  Note that other methods, such as
    /// describe(), may invoke getNodeElementList().
    ///
    /// To clarify: If this is empty, then it could be either that the
    /// Map is contiguous (meaning that we don't need to store all the
    /// global indices explicitly), or that the Map really does
    /// contain zero indices on the calling process.
    ///
    /// NOTE: With CUDA, we assume UVM, in that host code can access
    /// the entries of this View.
    ///
    /// This has LayoutLeft so that we can call Kokkos::deep_copy to
    /// copy this between any two Kokkos Devices.  Otherwise, the
    /// Devices might have different default layouts, thus forbidding
    /// a deep_copy.  We use LayoutLeft instead of LayoutRight because
    /// LayoutRight is the default on non-CUDA Devices, and we want to
    /// make sure we catch assignment or copying from the default to
    /// the nondefault layout.
    mutable Kokkos::View<const GlobalOrdinal*,
                         Kokkos::LayoutLeft,
                         device_type> lgMap_;

    /// \brief Host View of lgMap_.
    ///
    /// This is allocated along with lgMap_, on demand (lazily), by
    /// getNodeElementList() (which see).  It is also used by
    /// getGlobalElement() (which is a host method, and therefore
    /// requires a host View) if necessary (only noncontiguous Maps
    /// need this).
    mutable typename decltype (lgMap_)::HostMirror lgMapHost_;

    //! Type of a mapping from global IDs to local IDs.
    typedef Details::FixedHashTable<GlobalOrdinal, LocalOrdinal, device_type>
      global_to_local_table_type;

    /// \brief A mapping from global IDs to local IDs.
    ///
    /// This is a local mapping.  Directory implements the global
    /// mapping for all global indices (both remote and locally
    /// owned).  This object corresponds roughly to
    /// Epetra_BlockMapData's LIDHash_ hash table (which also maps
    /// from global to local indices).
    ///
    /// This mapping is built only for a noncontiguous map, by the
    /// noncontiguous map constructor.  For noncontiguous maps, the
    /// getLocalElement() and isNodeGlobalElement() methods use
    /// this mapping.
    global_to_local_table_type glMap_;

    /// \brief Object that can find the process rank and local index
    ///   for any given global index.
    ///
    /// Initializing this object is a collective operation over all
    /// processes in the Map's communicator.  (Creating it is not.)
    /// getRemoteIndexList() is the only method that needs this
    /// object, and it also happens to be a collective.  Thus, we
    /// initialize the Directory on demand in getRemoteIndexList().
    /// This saves the communication cost of initializing the
    /// Directory, for some Maps which are never involved in an Import
    /// or Export operation.  For example, a nonsquare sparse matrix
    /// (CrsMatrix) with row and range Maps the same would never need
    /// to construct an Export object.  This is a common case for the
    /// prolongation or restriction operators in algebraic multigrid.
    ///
    /// \note This is declared "mutable" so that getRemoteIndexList()
    ///   can initialize the Directory on demand.
    ///
    /// \warning The Directory is an implementation detail of its Map.
    ///   It does not make sense to expose in the public interface of
    ///   Map.  Resist the temptation to do so.  There is no need,
    ///   because Map's public interface already exposes the one
    ///   useful feature of Directory, via getRemoteIndexList().
    ///
    /// \note We use Teuchos::RCP (one could also use std::shared_ptr)
    ///   here because different views of the same Map (remember that
    ///   Map implements view semantics) may share the same Directory.
    ///   Map's three creation constructors (not copy constructor!)
    ///   create the Directory, but creation is lightweight.  Since
    ///   the Directory then exists (directory_ is not null), multiple
    ///   views of the same object share the same Directory, and all
    ///   of them will see the result if the Directory is initialized
    ///   by one of them.  Otherwise, if directory_ were to start out
    ///   null, then previously existing views of a Map could not
    ///   benefit from lazy creation of the Directory.
    ///
    mutable Teuchos::RCP<Directory<LocalOrdinal,GlobalOrdinal,Node> > directory_;

  }; // Map class

  /// \brief Nonmember constructor for a locally replicated Map with
  ///   the default Kokkos Node.
  ///
  /// This method returns a Map instantiated on the default Kokkos
  /// Node type.  The Map is configured to use zero-based indexing.
  ///
  /// \param numElements [in] Number of elements on each process.
  ///   Each process gets the same set of elements, namely <tt>0, 1,
  ///   ..., numElements - 1</tt>.
  ///
  /// \param comm [in] The Map's communicator.
  ///
  /// \relatesalso Map
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal> >
  createLocalMap (const size_t numElements,
                  const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

  /// \brief Nonmember constructor for a locally replicated Map with
  ///   a specified Kokkos Node.
  ///
  /// This method returns a Map instantiated on the given Kokkos Node
  /// instance.  The Map is configured to use zero-based indexing.
  ///
  /// \param numElements [in] Number of elements on each process.
  ///   Each process gets the same set of elements, namely <tt>0, 1,
  ///   ..., numElements - 1</tt>.
  ///
  /// \param comm [in] The Map's communicator.
  ///
  /// \param node [in] The Kokkos Node instance.  If not provided, we
  ///   will construct an instance for you.
  ///
  /// \relatesalso Map
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createLocalMapWithNode (const size_t numElements,
                          const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                          const Teuchos::RCP<Node>& node = Teuchos::null);

  /// \brief Non-member constructor for a uniformly distributed,
  ///   contiguous Map with the default Kokkos Node.
  ///
  /// This method returns a Map instantiated on the Kokkos default
  /// Node type.  The resulting Map uses zero-based indexing.
  ///
  /// \relatesalso Map
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal> >
  createUniformContigMap (const global_size_t numElements,
                          const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

  /// \brief Non-member constructor for a uniformly distributed,
  ///   contiguous Map with a user-specified Kokkos Node.
  ///
  /// The resulting Map uses zero-based indexing.
  ///
  /// \relatesalso Map
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createUniformContigMapWithNode (const global_size_t numElements,
                                  const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                  const Teuchos::RCP<Node>& node = Teuchos::null);

  /** \brief Non-member constructor for a (potentially) non-uniformly distributed, contiguous Map with the default Kokkos Node.

      This method returns a Map instantiated on the Kokkos default node type.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal> >
  createContigMap (const global_size_t numElements,
                   const size_t localNumElements,
                   const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

  /** \brief Non-member constructor for a (potentially) non-uniformly distributed, contiguous Map with a user-specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createContigMapWithNode (const global_size_t numElements,
                           const size_t localNumElements,
                           const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                           const Teuchos::RCP<Node>& node =
                             defaultArgNode<Node> ());

  /** \brief Non-member constructor for a non-contiguous Map with the default Kokkos Node.

      This method returns a Map instantiated on the Kokkos default node type.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal> >
  createNonContigMap (const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
                      const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

  /** \brief Non-member constructor for a non-contiguous Map with a user-specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createNonContigMapWithNode (const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
                              const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                              const Teuchos::RCP<Node>& node = Teuchos::null);

  /** \brief Non-member constructor for a contiguous Map with user-defined weights and a user-specified Kokkos Node.

      The Map is configured to use zero-based indexing.

      \relatesalso Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createWeightedContigMapWithNode (const int thisNodeWeight,
                                   const global_size_t numElements,
                                   const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                                   const Teuchos::RCP<Node>& node = Teuchos::null);

  /** \brief Creates a one-to-one version of the given Map where each GID is owned by only one process.

      The user must guarantee there are no duplicate GID on the same processor. Unexepected behavior may result.

      \relatesalso Map
   */
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createOneToOne (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &M);

  /** \brief Creates a one-to-one version of the given Map where each GID is owned by only one process.
             The rule to break ties is specifed by the tie break object.

      The user must guarantee there are no duplicate GID on the same processor. Unexepected behavior may result.

      \relatesalso Map
   */
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createOneToOne(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &M,
                 const Tpetra::Details::TieBreak<LocalOrdinal,GlobalOrdinal> & tie_break);

} // Tpetra namespace

#include "Tpetra_Directory_decl.hpp"

namespace Tpetra {
  namespace Details {

    template<class OutMapType, class InMapType>
    OutMapType
    MapCloner<OutMapType, InMapType>::
    clone (const InMapType& mapIn,
           const Teuchos::RCP<out_node_type>& nodeOut)
    {
      static_assert (std::is_same<typename OutMapType::local_ordinal_type,
                                  typename InMapType::local_ordinal_type>::value,
                     "Tpetra::Map clone: The LocalOrdinal template parameter "
                     "of the input and output Map types must be the same.");
      static_assert (std::is_same<typename OutMapType::global_ordinal_type,
                                  typename InMapType::global_ordinal_type>::value,
                     "Tpetra::Map clone: The GlobalOrdinal template parameter "
                     "of the input and output Map types must be the same.");
      typedef typename OutMapType::local_ordinal_type LO;
      typedef typename OutMapType::global_ordinal_type GO;
      typedef ::Tpetra::Directory<LO, GO,
                                  typename OutMapType::node_type> out_dir_type;
      typedef typename OutMapType::global_to_local_table_type out_table_type;
      typedef typename OutMapType::device_type out_device_type;

      OutMapType mapOut; // Make an empty Map.

      // If out_node_type::execution_space hasn't been initialized
      // (because nodeOut is null) initialize it now, by letting
      // defaultArgNode create the out_node_type instance.  Do this
      // first, since we'll be creating Kokkos::View instances with
      // out_node_type::device_type below, and this requires that
      // out_node_type::execution_space be initialized.
      mapOut.node_ = nodeOut.is_null () ? defaultArgNode<out_node_type> () : nodeOut;

      // Fill the new Map with (possibly) shallow copies of all of the
      // original Map's data.  This is safe because Map is immutable,
      // so users can't change the original Map.
      mapOut.comm_              = mapIn.comm_;
      mapOut.indexBase_         = mapIn.indexBase_;
      mapOut.numGlobalElements_ = mapIn.numGlobalElements_;
      mapOut.numLocalElements_  = mapIn.numLocalElements_;
      mapOut.minMyGID_          = mapIn.minMyGID_;
      mapOut.maxMyGID_          = mapIn.maxMyGID_;
      mapOut.minAllGID_         = mapIn.minAllGID_;
      mapOut.maxAllGID_         = mapIn.maxAllGID_;
      mapOut.firstContiguousGID_= mapIn.firstContiguousGID_;
      mapOut.lastContiguousGID_ = mapIn.lastContiguousGID_;
      mapOut.uniform_           = mapIn.uniform_;
      mapOut.contiguous_        = mapIn.contiguous_;
      mapOut.distributed_       = mapIn.distributed_;
      {
        // mfh 25 Dec 2015, 11 Jan 2016: We really only need to make a
        // deep copy if the two Map types have different memory
        // spaces.  However, if you're calling clone(), it is likely
        // the case that the memory spaces differ, so it doesn't hurt
        // to make a deep copy here.
        Kokkos::View<GO*, Kokkos::LayoutLeft, out_device_type>
          lgMapOut ("lgMap", mapIn.lgMap_.dimension_0 ());
        Kokkos::deep_copy (lgMapOut, mapIn.lgMap_);
        mapOut.lgMap_ = lgMapOut; // cast to const

        // mfh 11 Apr 2016: We can't just assign mapIn.lgMapHost_ to
        // mapOut.lgMapHost_ either.  This is because the memory space
        // of the host mirror of a CudaUVMSpace View is also
        // CudaUVMSpace, but the memory space of the host mirror of a
        // HostSpace View is HostSpace.  We can't assign one View to
        // another View with a different memory space.
        //
        // What we _can_ do here, though, is avoid a deep_copy in case
        // we're not using CUDA, by exploiting host mirrors.

        static_assert (std::is_same<typename decltype (mapOut.lgMapHost_)::array_layout,
                         typename decltype (mapIn.lgMapHost_)::array_layout>::value,
          "mapOut.lgMapHost_ and MapIn.lgMapHost_ do not have the same "
          "array_layout.  Please report this bug to the Tpetra developers.");

        // lgMapOut is nonconst, so use it here instead of mapOut.lgMap_.
        auto lgMapHostOut = Kokkos::create_mirror_view (lgMapOut);
        Kokkos::deep_copy (lgMapHostOut, lgMapOut);
        mapOut.lgMapHost_ = lgMapHostOut;
      }
      // This makes a deep copy only if necessary.  We could have
      // defined operator= to do this, but that would violate
      // expectations.  (Kokkos::View::operator= only does a shallow
      // copy, EVER.)
      mapOut.glMap_ = out_table_type (mapIn.glMap_);

      // We could cleverly clone the Directory here if it is
      // initialized, but there is no harm in simply creating it
      // uninitialized.
      mapOut.directory_ = Teuchos::rcp (new out_dir_type ());

      return mapOut;
    }
  } // namespace Details


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  template <class NodeOut>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, NodeOut> >
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  clone (const Teuchos::RCP<NodeOut>& nodeOut) const
  {
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> in_map_type;
    typedef Map<LocalOrdinal, GlobalOrdinal, NodeOut> out_map_type;
    typedef Details::MapCloner<out_map_type, in_map_type> cloner_type;
    // Copy constructor does a shallow copy.
    return Teuchos::rcp (new out_map_type (cloner_type::clone (*this, nodeOut)));
  }

  namespace Details {
    /// \brief Is map1 locally fitted to map2?
    ///
    /// \param map1 [in] The first Map
    /// \param map2 [in] The second Map
    ///
    /// For Map instances map1 and map2, we say that map1 is
    /// <i>locally fitted</i> to map2 (on the calling process), when
    /// the initial indices of map1 (on the calling process) are the
    /// same and in the same order as those of map2.  "Fittedness" is
    /// entirely a local (per MPI process) property.
    ///
    /// The predicate "is map1 fitted to map2 ?" is <i>not</i>
    /// symmetric.  For example, map2 may have more entries than map1.
    ///
    /// Fittedness on a process can let Tpetra avoid deep copies in
    /// some Export or Import (communication) operations.  Tpetra
    /// could use this, for example, in optimizing its sparse
    /// matrix-vector multiply.
    template <class LocalOrdinal,class GlobalOrdinal, class Node>
    bool
    isLocallyFitted (const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map1,
                     const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map2);
  } // namespace Details

} // namespace Tpetra

/// \brief True if map1 is the same as (in the sense of isSameAs()) map2, else false.
/// \relatesalso Tpetra::Map
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator== (const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map1,
                 const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map2)
{ return map1.isSameAs (map2); }

/// \brief True if map1 is not the same as (in the sense of isSameAs()) map2, else false.
/// \relatesalso Tpetra::Map
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator!= (const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map1,
                 const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map2)
{ return ! map1.isSameAs (map2); }


#endif // TPETRA_MAP_DECL_HPP

