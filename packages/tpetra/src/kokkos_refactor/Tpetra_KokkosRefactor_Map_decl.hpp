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
#include <Kokkos_UnorderedMap.hpp> // includes Kokkos_View.hpp
#include "Tpetra_KokkosRefactor_Details_Map.hpp"

namespace Tpetra {

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

    Map (const global_size_t globalNumIndices,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const LocalGlobal lg = GloballyDistributed,
         const Teuchos::RCP<node_type> &node = KokkosClassic::Details::getNode<node_type> ());

    Map (const global_size_t globalNumIndices,
         const size_t myNumIndices,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<node_type> &node = KokkosClassic::Details::getNode<node_type> ());

    Map (const global_size_t globalNumIndices,
         const Kokkos::View<const GlobalOrdinal*, device_type>& myGlobalIndices,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<node_type> &node = KokkosClassic::Details::getNode<node_type> ());

    Map (const global_size_t globalNumIndices,
         const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalIndices,
         const GlobalOrdinal indexBase,
         const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<node_type>& node = KokkosClassic::Details::getNode<node_type> ());

    //@}
    /// \name Methods safe to call in a Kokkos parallel kernel on the host.
    ///
    /// These methods are thread-safe and thread-scalable on the host.
    /// You may safely call them in a Kokkos parallel kernel on the
    /// host.  If DeviceType is not a host device (e.g., if it is
    /// Kokkos::Cuda), then you must get the underlying "device Map"
    /// object out in order to call these methods in a Kokkos parallel
    /// kernel on the device.  To do this, call getDeviceMap().
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

    /// \brief Get the device-memory version of the Map's implementation.
    ///
    /// The returned object has shallow-copy (handle or view) semantics.
    Details::Map<LocalOrdinal, GlobalOrdinal, device_type> getDeviceMap () {
      return mapDevice_;
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

  protected:
    // This lets other specializations of Map access all of this
    // specialization's internal methods and data, so that we can
    // implement clone() without exposing the details of Map to users.
    template <class LO, class GO, class N> friend class Map;

    /// \brief Default constructor (that does nothing).
    ///
    /// We use this in clone() and removeEmptyProcesses(), where we
    /// have the information to initialize the Map more efficiently
    /// ourselves, without going through one of the three usual Map
    /// construction paths.
    Map () {}

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

    typedef Details::Map<LocalOrdinal, GlobalOrdinal, device_type> device_impl_type;
    typedef typename device_type::host_mirror_device_type host_mirror_device_type;
    typedef Details::Map<LocalOrdinal, GlobalOrdinal, host_mirror_device_type> host_impl_type;

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
    mutable device_impl_type mapHost_;

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

      static Teuchos::RCP<const out_map_type>
      clone (const in_map_type& mapIn,
             const Teuchos::RCP<out_node_type>& nodeOut)
      {
        using Teuchos::rcp;

        if (mapIn.isUniform ()) {
          const Tpetra::LocalGlobal lg = mapIn.isDistributed () ?
            Tpetra::GloballyDistributed : Tpetra::LocallyReplicated;
          return rcp (new out_map_type (mapIn.getGlobalNumElements (),
                                        mapIn.getIndexBase (),
                                        mapIn.getComm (), lg, nodeOut));
        }
        else if (mapIn.isContiguous ()) {
          return rcp (new out_map_type (mapIn.getGlobalNumElements (),
                                        mapIn.getNodeNumElements (),
                                        mapIn.getIndexBase (),
                                        mapIn.getComm (), nodeOut));
        }
        else {
          return rcp (new out_map_type (mapIn.getGlobalNumElements (),
                                        mapIn.getNodeElementList (),
                                        mapIn.getIndexBase (),
                                        mapIn.getComm (), nodeOut));
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

      static Teuchos::RCP<const out_map_type>
      clone (const in_map_type& mapIn,
             const Teuchos::RCP<out_node_type>& nodeOut)
      {
        using Teuchos::rcp;

        if (mapIn.isUniform ()) {
          const Tpetra::LocalGlobal lg = mapIn.isDistributed () ?
            Tpetra::GloballyDistributed : Tpetra::LocallyReplicated;
          return rcp (new out_map_type (mapIn.getGlobalNumElements (),
                                        mapIn.getIndexBase (),
                                        mapIn.getComm (), lg, nodeOut));
        }
        else if (mapIn.isContiguous ()) {
          return rcp (new out_map_type (mapIn.getGlobalNumElements (),
                                        mapIn.getNodeNumElements (),
                                        mapIn.getIndexBase (),
                                        mapIn.getComm (), nodeOut));
        }
        else {
          return rcp (new out_map_type (mapIn.getGlobalNumElements (),
                                        mapIn.getNodeElementList (),
                                        mapIn.getIndexBase (),
                                        mapIn.getComm (), nodeOut));
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

      static Teuchos::RCP<const out_map_type>
      clone (const in_map_type& mapIn,
             const Teuchos::RCP<out_node_type>& nodeOut)
      {
        using Teuchos::RCP;
        typedef typename OutDeviceType::host_mirror_device_type
          host_mirror_device_type;
        RCP<out_map_type> mapOut (new out_map_type ()); // Make an empty Map.

        mapOut->comm_ = mapIn.comm_;
        mapOut->node_ = nodeOut;
        mapOut->mapDevice_.template create_copy_view<OutDeviceType> (mapIn.mapDevice_);
        mapOut->mapHost_.template create_copy_view<host_mirror_device_type> (mapIn.mapHost_);

        // mfh 02 Apr 2013: While Map only needs to create the Directory
        // on demand in getRemoteIndexList, we have a Directory here that
        // we can clone inexpensively, so there is no harm in creating it
        // here.
        if (! mapIn.directory_.is_null ()) {
          mapOut->directory_ = mapIn.directory_->template clone<out_node_type> (*mapOut);
        } else {
          mapOut->directory_ = Teuchos::null; // created on demand
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
    return Details::MapCloner<out_map_type, in_map_type>::clone (*this, outNode);
  }

} // namespace Tpetra

#endif // TPETRA_KOKKOSREFACTOR_MAP_DECL_HPP

