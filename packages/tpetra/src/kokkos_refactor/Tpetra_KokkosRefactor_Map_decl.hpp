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

#ifndef TPETRA_KOKKOSREFACTOR_MAP_DECL_HPP
#define TPETRA_KOKKOSREFACTOR_MAP_DECL_HPP

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include "Tpetra_KokkosRefactor_Details_Map.hpp"

namespace Tpetra {

  /// \brief Describes a parallel distribution of objects over
  ///   processes (Kokkos refactor specialization of Map).
  ///
  /// \tparam LocalOrdinal The type of local indices.  Should be an
  ///   integer, and generally should be signed.  A good model of
  ///   <tt>LocalOrdinal</tt> is \c int.  (In Epetra, this is always
  ///   just \c int.)
  ///
  /// \tparam GlobalOrdinal The type of global indices.  Should be an
  ///   integer, and generally should be signed.  Also,
  ///   <tt>sizeof(GlobalOrdinal)</tt> should be greater than to equal
  ///   to <tt>sizeof(LocalOrdinal)</tt>.  For example, if
  ///   <tt>LocalOrdinal</tt> is \c int, good models of
  ///   <tt>GlobalOrdinal</tt> are \c int, \c long, <tt>long long</tt>
  ///   (if the configure-time option
  ///   <tt>Teuchos_ENABLE_LONG_LONG_INT</tt> was set), or
  ///   <tt>ptrdiff_t</tt>.
  ///
  /// \tparam DeviceType The Kokkos Device type.
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
  /// \section Tpetra_KR_Map_prereq Prerequisites
  ///
  /// Before reading the rest of this documentation, it helps to know
  /// something about the Teuchos memory management classes, in
  /// particular Teuchos::RCP, Teuchos::ArrayRCP, and
  /// Teuchos::ArrayView, as well as Kokkos::View.  You should also
  /// know a little bit about MPI (the Message Passing Interface for
  /// distributed-memory programming).  You won't have to use MPI
  /// directly to use Map, but it helps to be familiar with the
  /// general idea of distributed storage of data over a communicator.
  ///
  /// \section Tpetra_KR_Map_concepts Map concepts
  ///
  /// \subsection Tpetra_KR_Map_local_vs_global Local and global indices
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
  /// \subsection Tpetra_KR_Map_contig Contiguous or noncontiguous
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
  ///   process in the Map's communicator, <li>
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
  /// \subsection Tpetra_KR_Map_dist_repl Globally distributed or locally replicated
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
  ///
  /// \section Tpetra_KR_Map_dev Host and device views
  ///
  /// Like all Tpetra objects that use the new Kokkos back-end, Map
  /// has "dual view" semantics.  This means that the data live in the
  /// preferred memory space of the \c DeviceType Kokkos device, but
  /// are also mirrored in host memory, if host memory is a different
  /// memory space.  Map is immutable, though, so it does not
  /// implement the synch() or modify() methods that let users
  /// transfer data between host and device memory, or mark either
  /// data set as modified.
  template <class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  class Map<LocalOrdinal,
            GlobalOrdinal,
            Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > :
    public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{

    //! The type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! The type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! The type of the Kokkos Node (less useful than device_type).
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;
    //! The type of the Kokkos Device (more useful than node_type).
    typedef DeviceType device_type;

  private:
    typedef Details::Map<LocalOrdinal, GlobalOrdinal, device_type> device_impl_type;
    typedef typename Kokkos::ViewTraits<LocalOrdinal,device_type>::host_mirror_space host_mirror_space;
    typedef Details::Map<LocalOrdinal, GlobalOrdinal, host_mirror_space> host_impl_type;

  public:

    /** \brief Constructor with Tpetra-defined contiguous uniform distribution.
     *
     * This constructor produces a Map with the following contiguous
     * range of <tt>numGlobalElements</tt> elements: <tt>indexBase,
     * indexBase + 1, ..., numGlobalElements + indexBase - 1</tt>.
     * Depending on the \c lg argument, the elements will either be
     * distributed evenly over all the processes in the given
     * communicator \c comm, or replicated on all processes in the
     * communicator.
     *
     * Preconditions on \c numGlobalElements and \c indexBase will
     * only be checked in a debug build (when Trilinos was configured
     * with CMake option <tt>Teuchos_ENABLE_DEBUG:BOOL=ON</tt>).  If
     * checks are enabled and any check fails, the constructor will
     * throw std::invalid_argument on all processes in the given
     * communicator.
     *
     * \param numGlobalElements [in] Number of elements in the Map
     *   (over all processes).
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
     *   elements.  If <tt>LocallyReplicated</tt>, then all processes
     *   will get the same set of elements, namely <tt>indexBase,
     *   indexBase + 1, ..., numGlobalElements + indexBase - 1</tt>.
     *
     * \param comm [in] Communicator over which to distribute the
     *   elements.
     *
     * \param node [in/out] Kokkos Node instance.  The type of this
     *   object must match the type of the Node template parameter.
     */
    Map (const global_size_t globalNumIndices,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const LocalGlobal lg = GloballyDistributed,
         const Teuchos::RCP<node_type> &node = KokkosClassic::Details::getNode<node_type> ());

    /** \brief Constructor with a user-defined contiguous distribution.
     *
     * If N is the sum of \c numLocalElements over all processes, then
     * this constructor produces a nonoverlapping Map with N elements
     * distributed over all the processes in the given communicator
     * <tt>comm</tt>, with either \c numLocalElements or
     * <tt>numLocalElements+1</tt> contiguous elements on the calling
     * process.
     *
     * Preconditions on \c numGlobalElements, \c numLocalElements, and
     * \c indexBase will only be checked in a debug build (when
     * Trilinos was configured with CMake option
     * <tt>TEUCHOS_ENABLE_DEBUG:BOOL=ON</tt>).  If checks are enabled
     * and any check fails, the constructor will throw
     * std::invalid_argument on all processes in the given
     * communicator.
     *
     * \param numGlobalElements [in] If <tt>numGlobalElements ==
     *   Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()</tt>,
     *   then the number of global elements will be computed (via a
     *   global communication) as the sum of numLocalElements over all
     *   processes.  Otherwise, it must equal the sum of
     *   numLocalElements over all processes.
     *
     * \param numLocalElements [in] Number of elements that the
     *   calling process will own in the Map.
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
     * \param node [in/out] Kokkos Node instance.  The type of this
     *   object must match the type of the Node template parameter.
     */
    Map (const global_size_t globalNumIndices,
         const size_t myNumIndices,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<node_type> &node = KokkosClassic::Details::getNode<node_type> ());

    /** \brief Constructor with user-defined arbitrary (possibly
     *    noncontiguous) distribution.
     *
     * Call this constructor if you have an arbitrary list of global
     * indices for each process in the given communicator.  Those
     * indices need not be contiguous, and the sets of global indices
     * on different processes may overlap.  This is the constructor to
     * use to make a general overlapping distribution.
     *
     * \param globalNumIndices [in] If <tt>globalNumIndices ==
     *   Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()</tt>,
     *   the number of global elements will be computed (via a global
     *   communication) as the sum of the counts of local elements.
     *   Otherwise, it must equal the sum of the local elements over
     *   all processes.  This will only be checked if Trilinos'
     *   Teuchos package was built with debug support (CMake Boolean
     *   option <tt>Teuchos_ENABLE_DEBUG:BOOL=ON</tt>).  If
     *   verification fails, the constructor will throw
     *   std::invalid_argument.
     *
     * \param myGlobalIndices [in] List of global indices owned by the
     *   calling process.
     *
     * \param indexBase [in] The base of the global indices in the
     *   Map.  This must be the same on every process in the given
     *   communicator.  Currently, Map requires that this equal the
     *   global minimum index over all processes' \c myGlobalIndices
     *   inputs.
     *
     * \param comm [in] Communicator over which to distribute the
     *   elements.
     *
     * \param node [in/out] Kokkos Node instance.  The type of this
     *   object must match the type of the Node template parameter.
     */
    Map (const global_size_t globalNumIndices,
         const Kokkos::View<const GlobalOrdinal*, device_type>& myGlobalIndices,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<node_type> &node = KokkosClassic::Details::getNode<node_type> ());

    /** \brief Constructor with user-defined arbitrary (possibly
     *    noncontiguous) distribution, as a Teuchos::ArrayView.
     *
     * This is like the previous noncontiguous constructor, except
     * that it takes the input list of global indices as a
     * Teuchos::ArrayView, rather than as a Kokkos::View.  We provide
     * this constructor for backwards compatibility.
     *
     * \param globalNumIndices [in] If <tt>globalNumIndices ==
     *   Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()</tt>,
     *   the number of global elements will be computed (via a global
     *   communication) as the sum of the counts of local elements.
     *   Otherwise, it must equal the sum of the local elements over
     *   all processes.  This will only be checked if Trilinos'
     *   Teuchos package was built with debug support (CMake Boolean
     *   option <tt>Teuchos_ENABLE_DEBUG:BOOL=ON</tt>).  If
     *   verification fails, the constructor will throw
     *   std::invalid_argument.
     *
     * \param myGlobalIndices [in] List of global indices owned by the
     *   calling process.
     *
     * \param indexBase [in] The base of the global indices in the
     *   Map.  This must be the same on every process in the given
     *   communicator.  Currently, Map requires that this equal the
     *   global minimum index over all processes' \c myGlobalIndices
     *   inputs.
     *
     * \param comm [in] Communicator over which to distribute the
     *   elements.
     *
     * \param node [in/out] Kokkos Node instance.  The type of this
     *   object must match the type of the Node template parameter.
     */
    Map (const global_size_t globalNumIndices,
         const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalIndices,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<node_type>& node = KokkosClassic::Details::getNode<node_type> ());

    /// \brief Default constructor (that does nothing).
    ///
    /// This only exists to support view semantics of Map.  That is,
    /// one can create an empty Map, and then assign a nonempty Map to
    /// it using operator=.
    ///
    /// This constructor is also useful in methods like clone() and
    /// removeEmptyProcesses(), where we have the information to
    /// initialize the Map more efficiently ourselves, without going
    /// through one of the three usual Map construction paths.
    Map ();

    //@}
    /// \name Methods safe to call in a Kokkos parallel kernel on the host.
    ///
    /// These methods are thread-safe and thread-scalable on the host.
    /// You may safely call them in a Kokkos parallel kernel on the
    /// host.  If DeviceType is not a host device (e.g., if it is
    /// Kokkos::Cuda), then you must get the underlying "device Map"
    /// object out in order to call these methods in a Kokkos parallel
    /// kernel on the device.  To do this, call getDeviceView().
    //@{

    /// \brief An invalid global index.
    ///
    /// Map's methods that return a global index use this index to
    /// indicate an error condition.  For example, if
    /// getGlobalElement() gets a local index that is not owned by the
    /// calling process, it returns an invalid global index, which
    /// equals the return value of this method.
    GlobalOrdinal invalidGlobalIndex () const;

    /// \brief An invalid local index.
    ///
    /// Map's methods that return a local index use this index to
    /// indicate an error condition.  For example, if
    /// getLocalElement() gets a global index that is not owned by the
    /// calling process, it returns an invalid local index, which
    /// equals the return value of this method.
    LocalOrdinal invalidLocalIndex () const;

    /// \brief Total number of indices over all processes in the Map's communicator.
    ///
    /// This is the sum, over all processes in the Map's communicator,
    /// of each process' number of owned indices.  Note that this
    /// counts global indices that are owned by multiple processes
    /// multiple times.
    global_size_t getGlobalNumElements() const;

    //! Number of indices owned by the calling process.
    size_t getNodeNumElements () const;

    //! Index base for this Map.
    GlobalOrdinal getIndexBase () const;

    /// \brief An invalid local index.
    ///
    /// Map's methods that return a local index use this index to
    /// indicate an error condition.  For example, if
    /// getLocalElement() gets a global index that is not owned by the
    /// calling process, it returns an invalid local index, which
    /// equals the return value of this method.
    LocalOrdinal getInvalidLocalIndex () const;

    //! The minimum local index (on the calling process).
    LocalOrdinal getMinLocalIndex () const;

    /// \brief The maximum local index on the calling process.
    ///
    /// If this process owns no indices, that is, if
    /// <tt>getNodeNumElements() == 0</tt>, then this method returns
    /// the same return value as getInvalidLocalIndex().
    LocalOrdinal getMaxLocalIndex () const;

    /// \brief An invalid global index.
    ///
    /// Map's methods that return a global index use this index to
    /// indicate an error condition.  For example, if
    /// getGlobalElement() gets a global index that is not owned by
    /// the calling process, it returns an invalid global index, which
    /// equals the return value of this method.
    GlobalOrdinal getInvalidGlobalIndex () const;

    //! The minimum global index owned by the calling process.
    GlobalOrdinal getMinGlobalIndex () const;

    //! The maximum global index owned by the calling process.
    GlobalOrdinal getMaxGlobalIndex () const;

    //! The minimum global index over all processes in the communicator.
    GlobalOrdinal getMinAllGlobalIndex () const;

    //! The maximum global index over all processes in the communicator.
    GlobalOrdinal getMaxAllGlobalIndex () const;

    /// \brief The local index corresponding to the given global index.
    ///
    /// If the given global index is not owned by this process, return
    /// the value returned by getInvalidLocalIndex();
    LocalOrdinal getLocalElement (const GlobalOrdinal globalIndex) const;

    /// \brief The global index corresponding to the given local index.
    ///
    /// If the given local index is not valid on the calling process,
    /// return the value returned by getInvalidGlobalIndex().
    GlobalOrdinal getGlobalElement (const LocalOrdinal localIndex) const;

    //! Whether the given local index is valid for this Map on this process.
    bool isNodeLocalElement (const LocalOrdinal localIndex) const;

    //! Whether the given global index is valid for this Map on this process.
    bool isNodeGlobalElement (const GlobalOrdinal globalIndex) const;

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

    /// \brief Whether this Map is globally distributed or locally replicated.
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

    //@}
    /// \name Methods <i>not</i> safe to call in a Kokkos parallel kernel.
    ///
    /// These methods are NOT thread-safe, on either host or device.
    /// You may <i>not</i> call them in a Kokkos parallel kernel.
    //@{

    /// \brief Whether the Map is one to one.
    ///
    /// This must be called collectively over all processes in the
    /// Map's communicator.
    bool isOneToOne () const;

    /// \fn getDeviceView
    /// \brief Get the version of the Map's implementation for the given Kokkos device.
    /// \tparam OutDeviceType The Kokkos device type of the returned object.
    ///
    /// This method will only compile if one of two conditions is met:
    /// <ol>
    /// <li> \c OutDeviceType is the same as \c device_type (this Map's
    ///      original Kokkos device type) </il>
    /// <li> \c OutDeviceType is compatible with
    ///      <tt>host_mirror_space</tt> </li>
    /// </ol>
    ///
    /// \return A Map implementation object of type
    ///   <tt>Details::Map<LocalOrdinal, GlobalOrdinal, OutDeviceType></tt>.
    ///   It has shallow-copy (handle or view) semantics.
    ///
    /// This method has a handy side effect: if it can do so cheaply,
    /// it initializes either the host or device Map implementation of
    /// this object, if it is not yet initialized but the other is.
    /// This is why the method is marked nonconst.
    template<class OutDeviceType>
    Details::Map<LocalOrdinal, GlobalOrdinal, OutDeviceType>
    getDeviceView () {
      // MapMirrorer finds which of mapDevice_ or mapHost_ (if either)
      // is initialized, and returns a mirror of that.  It's a shallow
      // copy if the device types are the same, else it's a deep copy.
      return Details::MapMirrorer<Details::Map<LocalOrdinal, GlobalOrdinal, OutDeviceType>,
        device_impl_type, host_impl_type>::mirror (mapDevice_, mapHost_);
    }

    /// \brief Return the process ranks and corresponding local
    ///   indices for the given global indices.
    ///
    /// This operation must always be called as a collective over all
    /// processes in the Map's communicator.  For a distributed
    /// noncontiguous Map, this operation requires communication.
    ///
    /// \param GIDs [in] List of global indices for which to find
    ///   process ranks and local indices.  These global indices need
    ///   not be owned by the calling process.  Indeed, they need not
    ///   be owned by any process.
    ///
    /// \param PIDs [out] List of process ranks corresponding to the
    ///   given global indices.  If a global index does not belong to
    ///   any process, the resulting process rank is -1.
    ///
    /// \param LIDs [out] List of local indices (that is, the local
    ///   index on the process that owns them) corresponding to the
    ///   given global indices.  If a global index does not have a
    ///   local index, the resulting local index has the value
    ///   returned by getInvalidLocalIndex().
    ///
    /// \pre PIDs.size() == GIDs.size()
    /// \pre LIDs.size() == GIDs.size()
    ///
    /// \return IDNotPresent indicates that for at least one global
    ///   index, we could not find the corresponding process rank.
    ///   Otherwise, return AllIDsPresent.
    ///
    /// \note This is crucial technology used in Export, Import,
    ///   CrsGraph, and CrsMatrix.
    LookupStatus
    getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal>& GIDs,
                        const Teuchos::ArrayView<                int>& PIDs,
                        const Teuchos::ArrayView<       LocalOrdinal>& LIDs) const;

    /// \brief Return the process ranks for the given global indices.
    ///
    /// This method must always be called as a collective over all
    /// processes in the Map's communicator.  For a distributed
    /// noncontiguous Map, this operation requires communication.
    ///
    /// \param GIDs [in] List of global indices for which to find
    ///   process ranks and local indices.  These global indices need
    ///   not be owned by the calling process.  Indeed, they need not
    ///   be owned by any process.
    ///
    /// \param PIDs [out] List of process ranks corresponding to the
    ///   given global indices.  If a global index does not belong to
    ///   any process, the resulting process rank is -1.
    ///
    /// \pre PIDs.size() == GIDs.size()
    ///
    /// \return IDNotPresent indicates that for at least one global
    ///   index, we could not find the corresponding process rank.
    ///   Otherwise, return AllIDsPresent.
    ///
    /// \note For a distributed noncontiguous Map, this operation
    ///   requires communication.  This is crucial technology used in
    ///   Export, Import, CrsGraph, and CrsMatrix.
    LookupStatus
    getRemoteIndexList (const Teuchos::ArrayView<const GlobalOrdinal> & GIDs,
                        const Teuchos::ArrayView<                int> & PIDs) const;

    /// \brief Return a view of the global indices owned by this process.
    ///
    /// If you call this method on a contiguous Map, it will create
    /// and cache the list of global indices for later use.  Beware of
    /// calling this if the calling process owns a very large number
    /// of global indices.
    Teuchos::ArrayView<const GlobalOrdinal> getNodeElementList () const;

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
    bool isCompatible (const Map<LocalOrdinal, GlobalOrdinal, node_type>& map) const;

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
    /// and global indices.  Two Maps with different numbers of
    /// processes in their communicators cannot be compatible, let
    /// alone identical.  Two identical Maps correspond to the same
    /// permutation.
    ///
    /// The behavior of this method is undefined if the input Map's
    /// communicator and this Map's communicator have different
    /// numbers of processes.  This method must be called collectively
    /// over this Map's communicator.
    bool isSameAs (const Map<LocalOrdinal, GlobalOrdinal, node_type>& map) const;

    //@}
    //! Accessors for the Teuchos::Comm and Kokkos Node objects.
    //@{

    //! Get this Map's Comm object.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

    //! Get this Map's Kokkos node instance.
    Teuchos::RCP<node_type> getNode () const;

    //@}
    //! Implementation of \c Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print this object with the given verbosity level to the given Teuchos::FancyOStream.
    void
    describe (Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}
    //! Advanced methods
    //@{

    //! Create a shallow copy of this Map, with a different Node type.
    template <class OutNodeType>
    Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, OutNodeType> >
    clone (const Teuchos::RCP<OutNodeType>& node2) const;

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
    RCP<const Map<LocalOrdinal, GlobalOrdinal, node_type> >
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
    RCP<const Map<LocalOrdinal, GlobalOrdinal, node_type> >
    replaceCommWithSubset (const Teuchos::RCP<const Teuchos::Comm<int> >& newComm) const;
    //@}


    /// \brief Is the given Map locally the same as the input Map?
    ///
    /// "Locally the same" means that on the calling process, the two
    /// Maps' global indices are the same and occur in the same order.
    bool locallySameAs (const Map<LocalOrdinal, GlobalOrdinal, node_type>& map) const;

  protected:
    // This lets other specializations of Map access all of this
    // specialization's internal methods and data, so that we can
    // implement clone() without exposing the details of Map to users.
    template <class LO, class GO, class N> friend class Map;

  private:
    template<class OutMapType, class InMapType>
    friend struct Details::MapCloner;

    /// \brief Create this Map's Directory, if it hasn't been created already.
    ///
    /// This method must be called collectively over all processes in
    /// the Map's communicator.
    ///
    /// This is declared "const" so that we can call it in
    /// getRemoteIndexList() to create the Directory on demand.
    void setupDirectory () const;

    //! The communicator over which this Map is distributed.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    /// \brief The Kokkos Node instance.
    ///
    /// Passing a Node to Map means that we don't have to pass a Node
    /// to downstream classes such as MultiVector, Vector, CrsGraph
    /// and CrsMatrix.
    Teuchos::RCP<node_type> node_;

    /// \brief Device-memory version of the Map's implementation.
    ///
    /// This has to be declared \c mutable, because Kokkos' idea of
    /// what methods should be const differs from current Tpetra's
    /// idea.
    mutable device_impl_type mapDevice_;

    /// \brief Host-memory version of the Map's implementation.
    ///
    /// This has to be declared \c mutable, because Kokkos' idea of
    /// what methods should be const differs from current Tpetra's
    /// idea.
    mutable host_impl_type mapHost_;

    /// \brief Object that can find the process rank and local index
    ///   for any given global index.
    ///
    /// Creating this object is a collective operation over all
    /// processes in the Map's communicator.  getRemoteIndexList() is
    /// the only method that needs this object, and it also happens to
    /// be a collective.  Thus, we create the Directory on demand in
    /// getRemoteIndexList().  This saves the communication cost of
    /// creating the Directory, for some Maps which are never involved
    /// in an Import or Export operation.  For example, a nonsquare
    /// sparse matrix (CrsMatrix) with row and range Maps the same
    /// would never need to construct an Export object.  This is a
    /// common case for the prolongation or restriction operators in
    /// algebraic multigrid.
    ///
    /// \note This is declared "mutable" so that getRemoteIndexList()
    ///   can create the Directory on demand.
    ///
    /// \warning The Directory is an implementation detail of its Map.
    ///   It does not make sense to expose in the public interface of
    ///   Map.  Resist the temptation to do so.  There is no need,
    ///   because Map's public interface already exposes the one
    ///   useful feature of Directory, via getRemoteIndexList().  We
    ///   only use Teuchos::RCP here because the more appropriate
    ///   std::unique_ptr is a C++11 feature and is therefore not
    ///   available to us.  Developers should not construe the use of
    ///   Teuchos::RCP as permission to share this object.
    mutable Teuchos::RCP<Directory<LocalOrdinal, GlobalOrdinal, node_type> > directory_;
  };



  namespace Details {
    template<class LO, class GO, class InDeviceType, class OutNodeType>
    struct MapCloner<Tpetra::Map<LO, GO, OutNodeType>,
                     Tpetra::Map<LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<InDeviceType> > > {
      typedef OutNodeType out_node_type;
      typedef Kokkos::Compat::KokkosDeviceWrapperNode<InDeviceType> in_node_type;
      typedef Tpetra::Map<LO, GO, OutNodeType> out_map_type;
      typedef Tpetra::Map<LO, GO, in_node_type> in_map_type;

      static out_map_type
      clone (const in_map_type& mapIn,
             const Teuchos::RCP<out_node_type>& nodeOut)
      {
        if (mapIn.isUniform ()) {
          const Tpetra::LocalGlobal lg = mapIn.isDistributed () ?
            Tpetra::GloballyDistributed : Tpetra::LocallyReplicated;
          return out_map_type (mapIn.getGlobalNumElements (),
                               mapIn.getIndexBase (),
                               mapIn.getComm (), lg, nodeOut);
        }
        else if (mapIn.isContiguous ()) {
          return out_map_type (mapIn.getGlobalNumElements (),
                               mapIn.getNodeNumElements (),
                               mapIn.getIndexBase (),
                               mapIn.getComm (), nodeOut);
        }
        else {
          return out_map_type (mapIn.getGlobalNumElements (),
                               mapIn.getNodeElementList (),
                               mapIn.getIndexBase (),
                               mapIn.getComm (), nodeOut);
        }
      }
    };

    template<class LO, class GO, class InNodeType, class OutDeviceType>
    struct MapCloner<Tpetra::Map<LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<OutDeviceType> >,
                     Tpetra::Map<LO, GO, InNodeType> > {
      typedef Kokkos::Compat::KokkosDeviceWrapperNode<OutDeviceType> out_node_type;
      typedef InNodeType in_node_type;
      typedef Tpetra::Map<LO, GO, out_node_type> out_map_type;
      typedef Tpetra::Map<LO, GO, in_node_type> in_map_type;

      static out_map_type
      clone (const in_map_type& mapIn,
             const Teuchos::RCP<out_node_type>& nodeOut)
      {
        if (mapIn.isUniform ()) {
          const Tpetra::LocalGlobal lg = mapIn.isDistributed () ?
            Tpetra::GloballyDistributed : Tpetra::LocallyReplicated;
          return out_map_type (mapIn.getGlobalNumElements (),
                               mapIn.getIndexBase (),
                               mapIn.getComm (), lg, nodeOut);
        }
        else if (mapIn.isContiguous ()) {
          return out_map_type (mapIn.getGlobalNumElements (),
                               mapIn.getNodeNumElements (),
                               mapIn.getIndexBase (),
                               mapIn.getComm (), nodeOut);
        }
        else {
          return out_map_type (mapIn.getGlobalNumElements (),
                               mapIn.getNodeElementList (),
                               mapIn.getIndexBase (),
                               mapIn.getComm (), nodeOut);
        }
      }
    };

    template<class LO, class GO, class InDeviceType, class OutDeviceType>
    struct MapCloner<Tpetra::Map<LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<OutDeviceType> >,
                     Tpetra::Map<LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<InDeviceType> > >
    {
      typedef Tpetra::Map<LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<OutDeviceType> > out_map_type;
      typedef Tpetra::Map<LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<InDeviceType> > in_map_type;
      typedef Kokkos::Compat::KokkosDeviceWrapperNode<OutDeviceType> out_node_type;

      static out_map_type
      clone (const in_map_type& mapIn,
             const Teuchos::RCP<out_node_type>& nodeOut)
      {
        typedef ::Tpetra::Directory<typename out_map_type::local_ordinal_type,
                                    typename out_map_type::global_ordinal_type,
                                    typename out_map_type::node_type> out_dir_type;

        out_map_type mapOut; // Make an empty Map.

        mapOut.comm_ = mapIn.comm_;
        mapOut.node_ = nodeOut;
        mapOut.mapDevice_.template create_copy_view<InDeviceType> (mapIn.mapDevice_);
        mapOut.mapHost_.template create_copy_view<typename out_map_type::host_mirror_space> (mapIn.mapHost_);

        // mfh 02 Apr 2013: While Map only needs to create the Directory
        // on demand in getRemoteIndexList, we have a Directory here that
        // we can clone inexpensively, so there is no harm in creating it
        // here.
        //
        // FIXME (mfh 09 May 2014) clone() doesn't quite work -- the
        // Directory has the wrong implementation type, for some
        // reason -- but it's always correct to let the output Map
        // initialize the Directory on its own.
        if (false && ! mapIn.directory_.is_null ()) {
          mapOut.directory_ = mapIn.directory_->template clone<out_node_type> (mapOut);
        } else {
          // It's created here, but not initialized yet.  The output
          // Map will initialize it on demand, if needed.
          mapOut.directory_ = Teuchos::rcp (new out_dir_type ());
        }
        return mapOut;
      }
    };

  } // namespace Details

  template <class LocalOrdinal, class GlobalOrdinal, class InDeviceType>
  template <class OutNodeType>
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, OutNodeType> >
  Map<LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<InDeviceType> >::
  clone (const Teuchos::RCP<OutNodeType>& outNode) const
  {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<InDeviceType> in_node_type;
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, in_node_type> in_map_type;
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, OutNodeType> out_map_type;
    typedef Tpetra::Details::MapCloner<out_map_type, in_map_type> cloner_type;
    // Copy constructor does a shallow copy.
    return Teuchos::rcp (new out_map_type (cloner_type::clone (*this, outNode)));
  }

} // namespace Tpetra

#endif // TPETRA_KOKKOSREFACTOR_MAP_DECL_HPP

