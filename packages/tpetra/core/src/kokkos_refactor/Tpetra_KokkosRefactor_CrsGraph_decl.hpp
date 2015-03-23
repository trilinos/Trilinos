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

#ifndef TPETRA_KOKKOSREFACTOR_CRSGRAPH_DECL_HPP
#define TPETRA_KOKKOSREFACTOR_CRSGRAPH_DECL_HPP

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_StaticCrsGraph.hpp>

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration (needed for "friend" inside CrsGraph)
  template <class S, class LO, class GO, class N, const bool isClassic>
  class CrsMatrix;

  namespace Experimental {
    // forward declaration (needed for "friend" inside CrsGraph)
    template <class S, class LO, class GO, class N>
    class BlockCrsMatrix;
  } // namespace Experimental
#endif // DOXYGEN_SHOULD_SKIP_THIS

  namespace Details {
    /// \brief Status of the graph's or matrix's storage, when not in
    ///   a fill-complete state.
    ///
    /// When a CrsGraph or CrsMatrix is <i>not</i> fill complete, its
    /// data live in one of three storage formats:
    ///
    /// <ol>
    /// <li> "2-D storage": The graph stores column indices as "array
    ///   of arrays," and the matrix stores values as "array of
    ///   arrays."  The graph <i>must</i> have k_numRowEntries_
    ///   allocated.  This only ever exists if the graph was created
    ///   with DynamicProfile.  A matrix with 2-D storage must own its
    ///   graph, and the graph must have 2-D storage. </li>
    ///
    /// <li> "Unpacked 1-D storage": The graph uses a row offsets
    ///   array, and stores column indices in a single array.  The
    ///   matrix also stores values in a single array.  "Unpacked"
    ///   means that there may be extra space in each row: that is,
    ///   the row offsets array only says how much space there is in
    ///   each row.  The graph must use k_numRowEntries_ to find out
    ///   how many entries there actually are in the row.  A matrix
    ///   with unpacked 1-D storage must own its graph, and the graph
    ///   must have unpacked 1-D storage. </li>
    ///
    /// <li> "Packed 1-D storage": The matrix may or may not own the
    ///   graph.  "Packed" means that there is no extra space in each
    ///   row.  Thus, the k_numRowEntries_ array is not necessary and
    ///   may have been deallocated.  If the matrix was created with a
    ///   constant ("static") graph, this must be true. </li>
    /// </ol>
    ///
    /// With respect to the Kokkos refactor version of Tpetra, "2-D
    /// storage" should be considered a legacy option.
    ///
    /// The phrase "When not in a fill-complete state" is important.
    /// When the graph is fill complete, it <i>always</i> uses 1-D
    /// "packed" storage.  However, if storage is "not optimized," we
    /// retain the 1-D unpacked or 2-D format, and thus retain this
    /// enum value.
    enum EStorageStatus {
      STORAGE_2D, //<! 2-D storage
      STORAGE_1D_UNPACKED, //<! 1-D "unpacked" storage
      STORAGE_1D_PACKED, //<! 1-D "packed" storage
      STORAGE_UB //<! Invalid value; upper bound on enum values
    };
  } // namespace Details

  /// \brief Partial specialization of the "Kokkos refactor" (new)
  ///   version of CrsGraph.
  template <class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  class CrsGraph<LocalOrdinal, // the type of local indices
                 GlobalOrdinal, // the type of global indices
                 Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, // the Kokkos Node type
                 false> : // the Kokkos refactor version of CrsGraph, NOT the "classic" version
    public RowGraph<LocalOrdinal,
                    GlobalOrdinal,
                    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >,
    public DistObject<GlobalOrdinal,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >,
    public Teuchos::ParameterListAcceptorDefaultBase
  {
    template <class S, class LO, class GO, class N, const bool isClassic>
    friend class CrsMatrix;
    template <class LO2, class GO2, class N2, const bool isClassic>
    friend class CrsGraph;
    template <class S, class LO, class GO, class N>
    friend class ::Tpetra::Experimental::BlockCrsMatrix;

    //! The specialization of DistObject that is this class' parent class.
    typedef DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > dist_object_type;

  public:
    //! This class' first template parameter; the type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! This class' second template parameter; the type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! This class' Kokkos execution space.
    typedef typename DeviceType::execution_space execution_space;
    /// \brief The Kokkos Node type used by this class.
    ///
    /// This type depends on the DeviceType template parameter.  In
    /// this, the Kokkos refactor version of Tpetra, it exists only
    /// for backwards compatibility.
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;

    //! The type of the part of the sparse graph on each MPI process.
    typedef Kokkos::StaticCrsGraph<LocalOrdinal,
                                   Kokkos::LayoutLeft,
                                   execution_space> local_graph_type;
    //! DEPRECATED; use local_graph_type (above) instead.
    typedef local_graph_type LocalStaticCrsGraphType TPETRA_DEPRECATED;

    //! DEPRECATED; use <tt>local_graph_type::row_map_type</tt> instead.
    typedef typename local_graph_type::row_map_type t_RowPtrs TPETRA_DEPRECATED;
    //! DEPRECATED; use <tt>local_graph_type::row_map_type::non_const_type</tt> instead.
    typedef typename local_graph_type::row_map_type::non_const_type t_RowPtrsNC TPETRA_DEPRECATED;
    //! DEPRECATED; use <tt>local_graph_type::entries_type::non_const_type</tt> instead.
    typedef typename local_graph_type::entries_type::non_const_type t_LocalOrdinal_1D TPETRA_DEPRECATED;

    //! The Map specialization used by this class.
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, node_type> map_type;
    //! The Import specialization used by this class.
    typedef Tpetra::Import<LocalOrdinal, GlobalOrdinal, node_type> import_type;
    //! The Export specialization used by this class.
    typedef Tpetra::Export<LocalOrdinal, GlobalOrdinal, node_type> export_type;

    //! @name Constructor/Destructor Methods
    //@{

    /// \brief Constructor specifying a single upper bound for the
    ///   number of entries in all rows on the calling process.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              size_t maxNumEntriesPerRow,
              ProfileType pftype = DynamicProfile,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a (possibly different) upper
    ///   bound for the number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  If pftype==DynamicProfile, this is
    ///   only a hint.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed the allocated
    ///   number of entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
              const ProfileType pftype = DynamicProfile,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a (possibly different) upper
    ///   bound for the number of entries in each row (legacy
    ///   KokkosClassic version).
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  If pftype==DynamicProfile, this is
    ///   only a hint.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed the allocated
    ///   number of entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::ArrayRCP<const size_t>& numEntPerRow,
              const ProfileType pftype = DynamicProfile,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and a single upper
    ///   bound for the number of entries in all rows on the calling
    ///   process.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const size_t maxNumEntriesPerRow,
              const ProfileType pftype = DynamicProfile,
              const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying column Map and number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  If pftype==DynamicProfile, this is
    ///   only a hint.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed the allocated
    ///   number of entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
              ProfileType pftype = DynamicProfile,
              const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying column Map and number of entries
    ///   in each row (legacy KokkosClassic version).
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  If pftype==DynamicProfile, this is
    ///   only a hint.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed the allocated
    ///   number of entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const Teuchos::ArrayRCP<const size_t>& numEntPerRow,
              ProfileType pftype = DynamicProfile,
              const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying column Map and arrays containing the graph in sorted, local ids.
    ///
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param rowPointers [in] The beginning of each row in the graph,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the graph.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const typename local_graph_type::row_map_type& rowPointers,
              const typename local_graph_type::entries_type::non_const_type& columnIndices,
              const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying column Map and arrays containing the graph in sorted, local ids.
    ///
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param rowPointers [in] The beginning of each row in the graph,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the graph.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const Teuchos::ArrayRCP<size_t> & rowPointers,
              const Teuchos::ArrayRCP<LocalOrdinal> & columnIndices,
              const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying column Map and a local (sorted)
    ///   graph, which the resulting CrsGraph views.
    ///
    /// Unlike most other CrsGraph constructors, successful completion
    /// of this constructor will result in a fill-complete graph.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param lclGraph [in] A locally indexed Kokkos::StaticCrsGraph
    ///   whose local row indices come from the specified row Map, and
    ///   whose local column indices come from the specified column
    ///   Map.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const local_graph_type& lclGraph,
              const Teuchos::RCP<Teuchos::ParameterList>& params);

    /// \brief Create a cloned CrsGraph for a different Node type.
    ///
    /// This method creates a new CrsGraph on a specified Kokkos Node
    /// type, with all of the entries of this CrsGraph object.
    ///
    /// \param node2 [in] Kokkos Node instance for constructing the
    ///   clone CrsGraph and its constituent objects.
    ///
    /// \param params [in/out] Optional list of parameters. If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.  See the list below for valid options.
    ///
    /// Parameters accepted by this method:
    /// - "Static profile clone" [bool, default: true] If \c true,
    ///   creates the clone with a static allocation profile. If
    ///   false, a dynamic allocation profile is used.
    /// - "Locally indexed clone" [bool] If \c true, fills clone
    ///   using this graph's column map and local indices (requires
    ///   that this graph have a column map.) If false, fills clone
    ///   using global indices and does not provide a column map. By
    ///   default, will use local indices only if this graph is using
    ///   local indices.
    /// - "fillComplete clone" [boolean, default: true] If \c true,
    ///   calls fillComplete() on the cloned CrsGraph object, with
    ///   parameters from \c params sublist "CrsGraph". The domain map
    ///   and range maps passed to fillComplete() are those of the map
    ///   being cloned, if they exist. Otherwise, the row map is used.
    template<class Node2>
    Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node2> >
    clone (const Teuchos::RCP<Node2> &node2,
           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const
    {
      typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node2> output_crs_graph_type;
      typedef CrsGraph<LocalOrdinal, GlobalOrdinal, node_type> input_crs_graph_type;
      typedef Details::CrsGraphCopier<output_crs_graph_type, input_crs_graph_type> copier_type;
      return copier_type::clone (*this, node2, params);
    }

    //! Destructor.
    virtual ~CrsGraph();

    //@}
    //! @name Implementation of Teuchos::ParameterListAcceptor
    //@{

    //! Set the given list of parameters (must be nonnull).
    void setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Default parameter list suitable for validation.
    Teuchos::RCP<const ParameterList> getValidParameters () const;

    //@}
    //! @name Insertion/Removal Methods
    //@{

    /// \brief Insert global indices into the graph.
    ///
    /// \pre \c globalRow is a valid index in the row Map.  It need
    ///   not be owned by the calling process.
    /// \pre <tt>isLocallyIndexed() == false</tt>
    /// \pre <tt>isStorageOptimized() == false</tt>
    ///
    /// \post <tt>indicesAreAllocated() == true</tt>
    /// \post <tt>isGloballyIndexed() == true</tt>
    ///
    /// If \c globalRow does not belong to the graph on this process,
    /// then it will be communicated to the appropriate process when
    /// globalAssemble() is called.  (That method will be called
    /// automatically during the next call to fillComplete().)
    /// Otherwise, the entries will be inserted into the part of the
    /// graph owned by the calling process.
    ///
    /// If the graph row already contains entries at the indices
    /// corresponding to values in \c indices, then the redundant
    /// indices will be eliminated.  This may happen either at
    /// insertion or during the next call to fillComplete().
    void
    insertGlobalIndices (GlobalOrdinal globalRow,
                         const Teuchos::ArrayView<const GlobalOrdinal>& indices);

    //! Insert local indices into the graph.
    /**
       \pre \c localRow is a local row belonging to the graph on this process.
       \pre <tt>isGloballyIndexed() == false</tt>
       \pre <tt>isStorageOptimized() == false</tt>
       \pre <tt>hasColMap() == true</tt>

       \post <tt>indicesAreAllocated() == true</tt>
       \post <tt>isLocallyIndexed() == true</tt>

       \note If the graph row already contains entries at the indices
         corresponding to values in \c indices, then the redundant
         indices will be eliminated; this may happen at insertion or
         during the next call to fillComplete().
    */
    void
    insertLocalIndices (const LocalOrdinal localRow,
                        const Teuchos::ArrayView<const LocalOrdinal> &indices);

    //! Remove all graph indices from the specified local row.
    /**
       \pre \c localRow is a local row of this graph.
       \pre <tt>isGloballyIndexed() == false</tt>
       \pre <tt>isStorageOptimized() == false</tt>

       \post <tt>getNumEntriesInLocalRow(localRow) == 0</tt>
       \post <tt>indicesAreAllocated() == true</tt>
       \post <tt>isLocallyIndexed() == true</tt>
    */
    void removeLocalIndices (LocalOrdinal localRow);

    //@}
    //! @name Transformational Methods
    /**
       Each of the methods in this group is a global collective. It is
       necessary to call these mehtods on all nodes participating in the
       communicator associated with this graph.
    */
    //@{

    /// \brief Communicate non-local contributions to other processes.
    ///
    /// This method is called automatically by fillComplete().
    /// Most users do not need to call this themselves,
    /// though we do permit this.
    void globalAssemble ();

    /*! Resume fill operations.
      After calling fillComplete(), resumeFill() must be called before initiating any changes to the graph.

      resumeFill() may be called repeatedly.

      \post  <tt>isFillActive() == true<tt>
      \post  <tt>isFillComplete() == false<tt>
    */
    void resumeFill (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /*! \brief Signal that data entry is complete, specifying domain and range maps.

      Off-process indices are distributed (via globalAssemble()),
      indices are sorted, redundant indices are eliminated, and
      global indices are transformed to local indices.

      \pre <tt>isFillActive() == true<tt>
      \pre <tt>isFillComplete()() == false<tt>

      \post <tt>isFillActive() == false<tt>
      \post <tt>isFillComplete() == true<tt>

      Parameters:
      - "Optimize Storage" (\c bool): Default is false.  If true,
        then isStorageOptimized() returns true after fillComplete
        finishes.  See isStorageOptimized() for consequences.
    */
    void
    fillComplete (const Teuchos::RCP<const map_type> &domainMap,
                  const Teuchos::RCP<const map_type> &rangeMap,
                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /*! \brief Signal that data entry is complete.

      Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

      \note This method calls fillComplete( getRowMap(), getRowMap(), os ). See parameter options there.
    */
    void fillComplete (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Perform a fillComplete on a graph that already has
    ///   data, via setAllIndices().
    ///
    /// The graph must already have filled local 1-D storage.  If the
    /// graph has been constructed in any other way, this method will
    /// throw an exception.  This routine is needed to support other
    /// Trilinos packages and should not be called by ordinary users.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    expertStaticFillComplete (const Teuchos::RCP<const map_type> & domainMap,
                              const Teuchos::RCP<const map_type> & rangeMap,
                              const Teuchos::RCP<const import_type> &importer=Teuchos::null,
                              const Teuchos::RCP<const export_type> &exporter=Teuchos::null,
                              const Teuchos::RCP<Teuchos::ParameterList> &params=Teuchos::null);
    //@}
    //! @name Methods implementing RowGraph.
    //@{

    //! Returns the communicator.
    Teuchos::RCP<const Comm<int> > getComm() const;

    //! Returns the underlying node.
    Teuchos::RCP<node_type> getNode() const;

    //! Returns the Map that describes the row distribution in this graph.
    Teuchos::RCP<const map_type> getRowMap () const;

    //! \brief Returns the Map that describes the column distribution in this graph.
    Teuchos::RCP<const map_type> getColMap () const;

    //! Returns the Map associated with the domain of this graph.
    Teuchos::RCP<const map_type> getDomainMap () const;

    //! Returns the Map associated with the domain of this graph.
    Teuchos::RCP<const map_type> getRangeMap () const;

    //! Returns the importer associated with this graph.
    Teuchos::RCP<const import_type> getImporter () const;

    //! Returns the exporter associated with this graph.
    Teuchos::RCP<const export_type> getExporter () const;

    //! Returns the number of global rows in the graph.
    /** Undefined if isFillActive().
     */
    global_size_t getGlobalNumRows() const;

    //! \brief Returns the number of global columns in the graph.
    /** Returns the number of entries in the domain map of the matrix.
        Undefined if isFillActive().
    */
    global_size_t getGlobalNumCols() const;

    //! Returns the number of graph rows owned on the calling node.
    size_t getNodeNumRows() const;

    //! Returns the number of columns connected to the locally owned rows of this graph.
    /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
     */
    size_t getNodeNumCols() const;

    //! Returns the index base for global indices for this graph.
    GlobalOrdinal getIndexBase() const;

    //! Returns the global number of entries in the graph.
    /** Undefined if isFillActive().
     */
    global_size_t getGlobalNumEntries() const;

    //! Returns the local number of entries in the graph.
    size_t getNodeNumEntries() const;

    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
    size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
    size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;

    //! \brief Returns the total number of indices allocated for the graph, across all rows on this node.
    /*! This is the allocation available to the user. Actual allocation may be larger, for example, after
      calling fillComplete(), and thus this does not necessarily reflect the memory consumption of the
      this graph.

      This quantity is computed during the actual allocation. Therefore, if <tt>indicesAreAllocated() == false</tt>,
      this method returns <tt>OrdinalTraits<size_t>::invalid()</tt>.
    */
    size_t getNodeAllocationSize() const;

    //! \brief Returns the current number of allocated entries for this node in the specified global row .
    /** Throws exception std::runtime_error if the specified global row does not belong to this node. */
    size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const;

    //! Returns the current number of allocated entries on this node in the specified local row.
    /** Throws exception std::runtime_error if the specified local row is not valid for this node. */
    size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const;

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons.
    /** Undefined if isFillActive().
     */
    global_size_t getGlobalNumDiags() const;

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons.
    /** Undefined if isFillActive().
     */
    size_t getNodeNumDiags() const;

    /// \brief Maximum number of entries in all rows over all processes.
    ///
    /// \note Undefined if isFillActive().
    ///
    /// \note This is the same as the result of a global maximum of
    ///   getNodeMaxNumRowEntries() over all processes.  That may not
    ///   necessarily mean what you think it does if some rows of the
    ///   matrix are owned by multiple processes.  In particular, some
    ///   processes might only own some of the entries in a particular
    ///   row.  This method only counts the number of entries in each
    ///   row that a process owns, not the total number of entries in
    ///   the row over all processes.
    size_t getGlobalMaxNumRowEntries() const;

    //! \brief Maximum number of entries in all rows owned by the calling process.
    /** Undefined if isFillActive().
     */
    size_t getNodeMaxNumRowEntries() const;

    /// \brief Whether the graph has a column Map.
    ///
    /// A CrsGraph has a column Map either because it was given to its
    /// constructor, or because it was constructed in fillComplete().
    /// Calling fillComplete() always makes a column Map if the graph
    /// does not already have one.
    ///
    /// A column Map lets the graph
    ///
    ///   - use local indices for storing entries in each row, and
    ///   - compute an Import from the domain Map to the column Map.
    ///
    /// The latter is mainly useful for a graph associated with a
    /// CrsMatrix.
    bool hasColMap() const;

    /// \brief Whether the graph is locally lower triangular.
    ///
    /// \pre <tt>! isFillActive()</tt>.
    ///   If fill is active, this method's behavior is undefined.
    ///
    /// \note This is entirely a local property.  That means this
    ///   method may return different results on different processes.
    bool isLowerTriangular() const;

    /// \brief Whether the graph is locally upper triangular.
    ///
    /// \pre <tt>! isFillActive()</tt>.
    ///   If fill is active, this method's behavior is undefined.
    ///
    /// \note This is entirely a local property.  That means this
    ///   method may return different results on different processes.
    bool isUpperTriangular() const;

    //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
    bool isLocallyIndexed() const;

    //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
    bool isGloballyIndexed() const;

    //! Returns \c true if fillComplete() has been called and the graph is in compute mode.
    bool isFillComplete() const;

    //! Returns \c true if resumeFill() has been called and the graph is in edit mode.
    bool isFillActive() const;

    /// \brief Whether graph indices in all rows are known to be sorted.
    ///
    /// A fill-complete graph is always sorted, as is a newly
    /// constructed graph. A graph is sorted immediately after calling
    /// resumeFill(), but any changes to the graph may result in the
    /// sorting status becoming unknown (and therefore, presumed
    /// unsorted).
    bool isSorted() const;

    //! \brief Returns \c true if storage has been optimized.
    /**
       Optimized storage means that the allocation of each row is equal to the
       number of entries. The effect is that a pass through the matrix, i.e.,
       during a mat-vec, requires minimal memory traffic. One limitation of
       optimized storage is that no new indices can be added to the graph.
    */
    bool isStorageOptimized() const;

    //! Returns \c true if the graph was allocated with static data structures.
    ProfileType getProfileType() const;

    /// \brief Get a copy of the given row, using global indices.
    ///
    /// \param GlobalRow [in] Global index of the row.
    /// \param Indices [out] On output: Global column indices.
    /// \param NumIndices [out] Number of indices returned.
    void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      const Teuchos::ArrayView<GlobalOrdinal>& Indices,
                      size_t& NumIndices) const;

    /// \brief Get a copy of the given row, using local indices.
    ///
    /// \param LocalRow [in] Local index of the row.
    /// \param Indices [out] On output: Local column indices.
    /// \param NumIndices [out] Number of indices returned.
    ///
    /// \pre <tt>hasColMap()</tt>
    void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     const Teuchos::ArrayView<LocalOrdinal>& indices,
                     size_t& NumIndices) const;

    //! Extract a const, non-persisting view of global indices in a specified row of the graph.
    /*!
      \param GlobalRow - (In) Global row number for which indices are desired.
      \param Indices   - (Out) Global column indices corresponding to values.
      \pre <tt>isLocallyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

      Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
    */
    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      Teuchos::ArrayView<const GlobalOrdinal>& Indices) const;

    //! Extract a const, non-persisting view of local indices in a specified row of the graph.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices  - (Out) Global column indices corresponding to values.
      \pre <tt>isGloballyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

      Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
    */
    void
    getLocalRowView (LocalOrdinal LocalRow,
                     Teuchos::ArrayView<const LocalOrdinal>& indices) const;

    //@}
    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const;

    //! Print the object to the given output stream with given verbosity level.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;

    //@}
    //! \name Implementation of DistObject
    //@{

    virtual bool
    checkSizes (const SrcDistObject& source);

    virtual void
    copyAndPermute (const SrcDistObject& source,
                    size_t numSameIDs,
                    const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                    const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs);

    virtual void
    packAndPrepare (const SrcDistObject& source,
                    const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                    Teuchos::Array<GlobalOrdinal> &exports,
                    const Teuchos::ArrayView<size_t> & numPacketsPerLID,
                    size_t& constantNumPackets,
                    Distributor &distor);

    virtual void
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
          Teuchos::Array<GlobalOrdinal>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets,
          Distributor& distor) const;

    virtual void
    unpackAndCombine (const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                      const Teuchos::ArrayView<const GlobalOrdinal> &imports,
                      const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                      size_t constantNumPackets,
                      Distributor &distor,
                      CombineMode CM);
    //@}
    //! \name Advanced methods, at increased risk of deprecation.
    //@{

    /// \brief Get an upper bound on the number of entries that can be
    ///   stored in each row.
    ///
    /// When a CrsGraph is constructed, callers must give an upper
    /// bound on the number of entries in each local row.  They may
    /// either supply a single integer which is the upper bound for
    /// all local rows, or they may give an array with a possibly
    /// different upper bound for each local row.
    ///
    /// This method returns the upper bound for each row.  If
    /// numEntriesPerLocalRowBound is Teuchos::null on output and
    /// boundSameForAllLocalRows is true on output, then
    /// numEntriesAllLocalRowsBound is the upper bound for all local
    /// rows.  If boundSameForAllLocalRows is false on output, then
    /// numEntriesPerLocalRowBound has zero or more entries on output,
    /// and numEntriesPerLocalRowBound[i_local] is the upper bound for
    /// local row i_local.
    ///
    /// The output argument boundSameForAllLocalRows is conservative;
    /// it only tells us whether boundForAllLocalRows has a meaningful
    /// value on output.  We don't necessarily check whether all
    /// entries of boundPerLocalRow are the same.
    void
    getNumEntriesPerLocalRowUpperBound (Teuchos::ArrayRCP<const size_t>& boundPerLocalRow,
                                        size_t& boundForAllLocalRows,
                                        bool& boundSameForAllLocalRows) const;

    /// \brief Set the graph's data directly, using 1-D storage.
    ///
    /// \pre <tt>hasColMap() == true</tt>
    /// \pre <tt>rowPointers.size() != getNodeNumRows()+1</tt>
    /// \pre No insert routines have been called.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    setAllIndices (const typename local_graph_type::row_map_type& rowPointers,
                   const typename local_graph_type::entries_type::non_const_type& columnIndices);

    /// \brief Set the graph's data directly, using 1-D storage.
    ///
    /// \pre <tt>hasColMap() == true</tt>
    /// \pre <tt>rowPointers.size() != getNodeNumRows()+1</tt>
    /// \pre No insert routines have been called.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    setAllIndices (const Teuchos::ArrayRCP<size_t> & rowPointers,
                   const Teuchos::ArrayRCP<LocalOrdinal> & columnIndices);

    /// \brief Get a host view of the row offsets.
    ///
    /// \note Please prefer getLocalGraph() to get the row offsets.
    ///
    /// This may return either a copy or a view of the row offsets.
    /// In either case, it will <i>always</i> live in host memory,
    /// never in (CUDA) device memory.
    Teuchos::ArrayRCP<const size_t> getNodeRowPtrs () const;

    //! Get an Teuchos::ArrayRCP of the packed column-indices.
    /*!  The returned buffer exists in host-memory.
     */
    Teuchos::ArrayRCP<const LocalOrdinal> getNodePackedIndices() const;

    /// \brief Replace the graph's current column Map with the given Map.
    ///
    /// This <i>only</i> replaces the column Map.  It does <i>not</i>
    /// change the graph's current column indices, or otherwise apply
    /// a permutation.  For example, suppose that before calling this
    /// method, the calling process owns a row containing local column
    /// indices [0, 2, 4].  These indices do <i>not</i> change, nor
    /// does their order change, as a result of calling this method.
    ///
    /// \param newColMap [in] New column Map.  Must be nonnull.
    void replaceColMap (const Teuchos::RCP<const map_type>& newColMap);

    /// \brief Reindex the column indices in place, and replace the
    ///   column Map.  Optionally, replace the Import object as well.
    ///
    /// \pre On every calling process, every index owned by the
    ///   current column Map must also be owned by the new column Map.
    ///
    /// \pre If the new Import object is provided, the new Import
    ///   object's source Map must be the same as the current domain
    ///   Map, and the new Import's target Map must be the same as the
    ///   new column Map.
    ///
    /// \param newColMap [in] New column Map.  Must be nonnull.
    ///
    /// \param newImport [in] New Import object.  Optional; computed
    ///   if not provided or if null.  Computing an Import is
    ///   expensive, so it is worth providing this if you can.
    ///
    /// \param sortIndicesInEachRow [in] If true, sort the indices in
    ///   each row after reindexing.
    void
    reindexColumns (const Teuchos::RCP<const map_type>& newColMap,
                    const Teuchos::RCP<const import_type>& newImport = Teuchos::null,
                    const bool sortIndicesInEachRow = true);

    /// \brief Replace the current domain Map and Import with the given parameters.
    ///
    /// \warning This method is ONLY for use by experts.
    /// \warning We make NO promises of backwards compatibility.
    ///   This method may change or disappear at any time.
    ///
    /// \pre <tt>isFillComplete() == true<tt>
    /// \pre <tt>isFillActive() == false<tt>
    /// \pre Either the given Import object is null, or the target Map
    ///   of the given Import is the same as this graph's column Map.
    /// \pre Either the given Import object is null, or the source Map
    ///    of the given Import is the same as this graph's domain Map.
    void
    replaceDomainMapAndImporter (const Teuchos::RCP<const map_type>& newDomainMap,
                                 const Teuchos::RCP<const import_type>& newImporter);

    /// \brief Remove processes owning zero rows from the Maps and their communicator.
    ///
    /// \warning This method is ONLY for use by experts.  We highly
    ///   recommend using the nonmember function of the same name
    ///   defined in Tpetra_DistObject_decl.hpp.
    ///
    /// \warning We make NO promises of backwards compatibility.
    ///   This method may change or disappear at any time.
    ///
    /// \param newMap [in] This <i>must</i> be the result of calling
    ///   the removeEmptyProcesses() method on the row Map.  If it
    ///   is not, this method's behavior is undefined.  This pointer
    ///   will be null on excluded processes.
    ///
    /// This method satisfies the strong exception guarantee, as
    /// long the destructors of Export, Import, and Map do not throw
    /// exceptions.  This means that either the method returns
    /// normally (without throwing an exception), or there are no
    /// externally visible side effects.  However, this does not
    /// guarantee no deadlock when the graph's original communicator
    /// contains more than one process.  In order to prevent
    /// deadlock, you must still wrap this call in a try/catch block
    /// and do an all-reduce over all processes in the original
    /// communicator to test whether the call succeeded.  This
    /// safety measure should usually be unnecessary, since the
    /// method call should only fail on user error or failure to
    /// allocate memory.
    virtual void
    removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap);
    //@}

    template<class ViewType, class OffsetViewType >
    struct pack_functor {
      typedef typename ViewType::execution_space execution_space;
      ViewType src;
      ViewType dest;
      OffsetViewType src_offset;
      OffsetViewType dest_offset;
      typedef typename OffsetViewType::non_const_value_type ScalarIndx;

      pack_functor(ViewType dest_, ViewType src_, OffsetViewType dest_offset_, OffsetViewType src_offset_):
        src(src_),dest(dest_),src_offset(src_offset_),dest_offset(dest_offset_) {};

      KOKKOS_INLINE_FUNCTION
      void operator() (size_t row) const {
        ScalarIndx i = src_offset(row);
        ScalarIndx j = dest_offset(row);
        const ScalarIndx k = dest_offset(row+1);
        for(;j<k;j++,i++) {
          dest(j) = src(i);
        }
      }
    };

  protected:
    // these structs are conveniences, to cut down on the number of
    // arguments to some of the methods below.
    struct SLocalGlobalViews {
      Teuchos::ArrayView<const GlobalOrdinal> ginds;
      Teuchos::ArrayView<const LocalOrdinal>  linds;
    };
    struct SLocalGlobalNCViews {
      Teuchos::ArrayView<GlobalOrdinal>       ginds;
      Teuchos::ArrayView<LocalOrdinal>        linds;
    };

    bool indicesAreAllocated () const;
    void allocateIndices (const ELocalGlobal lg);

    template <class T>
    Teuchos::ArrayRCP<T> allocateValues1D () const {
      // FIXME (mfh 05 Aug 2014) This method only exists here to
      // placate the ETI macros for Kokkos classic.
      return Teuchos::null;
    }

    template <class T>
    Teuchos::ArrayRCP<Teuchos::Array<T> > allocateValues2D () const;

    template <class T>
    RowInfo updateLocalAllocAndValues (const RowInfo rowInfo,
                                       const size_t newAllocSize,
                                       Teuchos::Array<T>& rowVals)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( ! isLocallyIndexed () );
      TEUCHOS_TEST_FOR_EXCEPT( ! indicesAreAllocated() );
      TEUCHOS_TEST_FOR_EXCEPT( newAllocSize == 0 );
      TEUCHOS_TEST_FOR_EXCEPT( newAllocSize < rowInfo.allocSize );
      TEUCHOS_TEST_FOR_EXCEPT( ! rowMap_->isNodeLocalElement (rowInfo.localRow) );
#endif // HAVE_TPETRA_DEBUG

      // Teuchos::ArrayRCP::resize automatically copies over values on reallocation.
      lclInds2D_[rowInfo.localRow].resize (newAllocSize);
      rowVals.resize (newAllocSize);
      nodeNumAllocated_ += (newAllocSize - rowInfo.allocSize);

      RowInfo rowInfoOut = rowInfo;
      rowInfoOut.allocSize = newAllocSize;
      return rowInfoOut;
    }

    template <class T>
    RowInfo
    updateGlobalAllocAndValues (const RowInfo rowInfo,
                                const size_t newAllocSize,
                                Teuchos::Array<T>& rowVals)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( ! isGloballyIndexed () );
      TEUCHOS_TEST_FOR_EXCEPT( ! indicesAreAllocated () );
      TEUCHOS_TEST_FOR_EXCEPT( newAllocSize == 0 );
      TEUCHOS_TEST_FOR_EXCEPT( newAllocSize < rowInfo.allocSize );
      TEUCHOS_TEST_FOR_EXCEPT( ! rowMap_->isNodeLocalElement (rowInfo.localRow) );
#endif // HAVE_TPETRA_DEBUG

      // Teuchos::ArrayRCP::resize automatically copies over values on reallocation.
      gblInds2D_[rowInfo.localRow].resize (newAllocSize);
      rowVals.resize (newAllocSize);
      nodeNumAllocated_ += (newAllocSize - rowInfo.allocSize);

      RowInfo rowInfoOut = rowInfo;
      rowInfoOut.allocSize = newAllocSize;
      return rowInfoOut;
    }

    //! \name Methods governing changes between global and local indices
    //@{

    //! Make the graph's column Map, if it does not already have one.
    void makeColMap ();
    void makeIndicesLocal ();
    void makeImportExport ();

    //@}
    //! \name Methods for inserting indices or transforming values
    //@{

    template<ELocalGlobal lg>
    size_t filterIndices (const SLocalGlobalNCViews& inds) const
    {
      using Teuchos::ArrayView;
      Teuchos::CompileTimeAssert<lg != GlobalIndices && lg != LocalIndices> cta_lg;
      (void)cta_lg;

      const map_type& cmap = *colMap_;
      size_t numFiltered = 0;
#ifdef HAVE_TPETRA_DEBUG
      size_t numFiltered_debug = 0;
#endif
      if (lg == GlobalIndices) {
        ArrayView<GlobalOrdinal> ginds = inds.ginds;
        typename ArrayView<GlobalOrdinal>::iterator fend = ginds.begin();
        typename ArrayView<GlobalOrdinal>::iterator cptr = ginds.begin();
        while (cptr != ginds.end()) {
          if (cmap.isNodeGlobalElement(*cptr)) {
            *fend++ = *cptr;
#ifdef HAVE_TPETRA_DEBUG
            ++numFiltered_debug;
#endif
          }
          ++cptr;
        }
        numFiltered = fend - ginds.begin();
      }
      else if (lg == LocalIndices) {
        ArrayView<LocalOrdinal> linds = inds.linds;
        typename ArrayView<LocalOrdinal>::iterator fend = linds.begin();
        typename ArrayView<LocalOrdinal>::iterator cptr = linds.begin();
        while (cptr != linds.end()) {
          if (cmap.isNodeLocalElement(*cptr)) {
            *fend++ = *cptr;
#ifdef HAVE_TPETRA_DEBUG
            ++numFiltered_debug;
#endif
          }
          ++cptr;
        }
        numFiltered = fend - linds.begin();
      }
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( numFiltered != numFiltered_debug );
#endif
      return numFiltered;
    }


    template<class T>
    size_t
    filterGlobalIndicesAndValues (const Teuchos::ArrayView<GlobalOrdinal>& ginds,
                                  const Teuchos::ArrayView<T>& vals) const
    {
      using Teuchos::ArrayView;
      const map_type& cmap = *colMap_;
      size_t numFiltered = 0;
      typename ArrayView<T>::iterator fvalsend = vals.begin();
      typename ArrayView<T>::iterator valscptr = vals.begin();
#ifdef HAVE_TPETRA_DEBUG
      size_t numFiltered_debug = 0;
#endif
      typename ArrayView<GlobalOrdinal>::iterator fend = ginds.begin();
      typename ArrayView<GlobalOrdinal>::iterator cptr = ginds.begin();
      while (cptr != ginds.end()) {
        if (cmap.isNodeGlobalElement (*cptr)) {
          *fend++ = *cptr;
          *fvalsend++ = *valscptr;
#ifdef HAVE_TPETRA_DEBUG
          ++numFiltered_debug;
#endif
        }
        ++cptr;
        ++valscptr;
      }
      numFiltered = fend - ginds.begin();
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( numFiltered != numFiltered_debug );
      TEUCHOS_TEST_FOR_EXCEPT( valscptr != vals.end() );
      const size_t numFilteredActual =
        static_cast<size_t> (fvalsend - vals.begin ());
      TEUCHOS_TEST_FOR_EXCEPT( numFiltered != numFilteredActual );
#endif // HAVE_TPETRA_DEBUG
      return numFiltered;
    }

    template<class T>
    size_t
    filterLocalIndicesAndValues (const Teuchos::ArrayView<LocalOrdinal>& linds,
                                 const Teuchos::ArrayView<T>& vals) const
    {
      using Teuchos::ArrayView;
      const map_type& cmap = *colMap_;
      size_t numFiltered = 0;
      typename ArrayView<T>::iterator fvalsend = vals.begin();
      typename ArrayView<T>::iterator valscptr = vals.begin();
#ifdef HAVE_TPETRA_DEBUG
      size_t numFiltered_debug = 0;
#endif
      typename ArrayView<LocalOrdinal>::iterator fend = linds.begin();
      typename ArrayView<LocalOrdinal>::iterator cptr = linds.begin();
      while (cptr != linds.end()) {
        if (cmap.isNodeLocalElement (*cptr)) {
          *fend++ = *cptr;
          *fvalsend++ = *valscptr;
#ifdef HAVE_TPETRA_DEBUG
          ++numFiltered_debug;
#endif
        }
        ++cptr;
        ++valscptr;
      }
      numFiltered = fend - linds.begin();
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( numFiltered != numFiltered_debug );
      TEUCHOS_TEST_FOR_EXCEPT( valscptr != vals.end() );
      const size_t numFilteredActual =
        Teuchos::as<size_t> (fvalsend - vals.begin ());
      TEUCHOS_TEST_FOR_EXCEPT( numFiltered != numFilteredActual );
#endif
      return numFiltered;
    }

    /// \brief Insert indices into the given row.
    ///
    /// \pre <tt>! (lg == LocalIndices && I == GlobalIndices)</tt>.
    ///   It does not make sense to give this method local column
    ///   indices (meaning that the graph has a column Map), yet to
    ///   ask it to store global indices.
    ///
    /// \param rowInfo [in] Result of CrsGraph's getRowInfo() or
    ///   updateAllocAndValues() methods, for the locally owned row
    ///   (whose local index is <tt>rowInfo.localRow</tt>) for which
    ///   you want to insert indices.
    ///
    /// \param newInds [in] View of the column indices to insert.  If
    ///   <tt>lg == GlobalIndices</tt>, then newInds.ginds, a
    ///   <tt>Teuchos::ArrayView<const GlobalOrdinal></tt>, contains
    ///   the (global) column indices to insert.  Otherwise, if <tt>lg
    ///   == LocalIndices</tt>, then newInds.linds, a
    ///   <tt>Teuchos::ArrayView<const LocalOrdinal></tt>, contains
    ///   the (local) column indices to insert.
    ///
    /// \param lg If <tt>lg == GlobalIndices</tt>, then the input
    ///   indices (in \c newInds) are global indices.  Otherwise, if
    ///   <tt>lg == LocalIndices</tt>, the input indices are local
    ///   indices.
    ///
    /// \param I If <tt>lg == GlobalIndices</tt>, then this method
    ///   will store the input indices as global indices.  Otherwise,
    ///   if <tt>I == LocalIndices</tt>, this method will store the
    ///   input indices as local indices.
    size_t
    insertIndices (const RowInfo& rowInfo,
                   const SLocalGlobalViews& newInds,
                   const ELocalGlobal lg,
                   const ELocalGlobal I);

    /// \brief Insert indices and their values into the given row.
    ///
    /// \tparam Scalar The type of a single value.  When this method
    ///   is called by CrsMatrix, \c Scalar corresponds to the first
    ///   template parameter of CrsMatrix.
    ///
    /// \pre <tt>! (lg == LocalIndices && I == GlobalIndices)</tt>.
    ///   It does not make sense to give this method local column
    ///   indices (meaning that the graph has a column Map), yet to
    ///   ask it to store global indices.
    ///
    /// \param rowInfo [in] Result of CrsGraph's getRowInfo() or
    ///   updateAllocAndValues() methods, for the locally owned row
    ///   (whose local index is <tt>rowInfo.localRow</tt>) for which
    ///   you want to insert indices.
    ///
    /// \param newInds [in] View of the column indices to insert.  If
    ///   <tt>lg == GlobalIndices</tt>, then newInds.ginds, a
    ///   <tt>Teuchos::ArrayView<const GlobalOrdinal></tt>, contains
    ///   the (global) column indices to insert.  Otherwise, if <tt>lg
    ///   == LocalIndices</tt>, then newInds.linds, a
    ///   <tt>Teuchos::ArrayView<const LocalOrdinal></tt>, contains
    ///   the (local) column indices to insert.
    ///
    /// \param oldRowVals [out] View of the current values.  They will
    ///   be overwritten with the new values.
    ///
    /// \param newRowVals [in] View of the new values.  They will be
    ///   copied over the old values.
    ///
    /// \param lg If <tt>lg == GlobalIndices</tt>, then the input
    ///   indices (in \c newInds) are global indices.  Otherwise, if
    ///   <tt>lg == LocalIndices</tt>, the input indices are local
    ///   indices.
    ///
    /// \param I If <tt>lg == GlobalIndices</tt>, then this method
    ///   will store the input indices as global indices.  Otherwise,
    ///   if <tt>I == LocalIndices</tt>, this method will store the
    ///   input indices as local indices.
    template<class Scalar>
    void
    insertIndicesAndValues (const RowInfo& rowInfo,
                            const SLocalGlobalViews& newInds,
                            const Teuchos::ArrayView<Scalar>& oldRowVals,
                            const Teuchos::ArrayView<const Scalar>& newRowVals,
                            const ELocalGlobal lg,
                            const ELocalGlobal I);
    void
    insertGlobalIndicesImpl (const LocalOrdinal myRow,
                             const Teuchos::ArrayView<const GlobalOrdinal> &indices);
    void
    insertLocalIndicesImpl (const LocalOrdinal myRow,
                            const Teuchos::ArrayView<const LocalOrdinal> &indices);
    //! Like insertLocalIndices(), but with column Map filtering.
    void
    insertLocalIndicesFiltered (const LocalOrdinal localRow,
                                const Teuchos::ArrayView<const LocalOrdinal> &indices);

    //! Like insertGlobalIndices(), but with column Map filtering.
    void
    insertGlobalIndicesFiltered (const GlobalOrdinal localRow,
                                 const Teuchos::ArrayView<const GlobalOrdinal> &indices);

    /// \brief Transform the given values using local indices.
    ///
    /// \param rowInfo [in] Information about a given row of the graph.
    ///
    /// \param rowVals [in/out] The values to be transformed.  They
    ///   correspond to the row indicated by rowInfo.
    ///
    /// \param inds [in] The (local) indices in the row, for which
    ///   to transform the corresponding values in rowVals.
    ///
    /// \param newVals [in] Values to use for transforming rowVals.
    ///   It's probably OK for these to alias rowVals.
    ///
    /// \param f [in] A binary function used to transform rowVals.
    ///
    /// This method transforms the values using the expression
    /// \code
    /// newVals[k] = f( rowVals[k], newVals[j] );
    /// \endcode
    /// where k is the local index corresponding to
    /// <tt>inds[j]</tt>.  It ignores elements of \c inds that are
    /// not owned by the calling process.
    ///
    /// \return The number of valid local column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    template<class Scalar, class BinaryFunction>
    LocalOrdinal
    transformLocalValues (RowInfo rowInfo,
                          const Teuchos::ArrayView<Scalar>& rowVals,
                          const Teuchos::ArrayView<const LocalOrdinal>& inds,
                          const Teuchos::ArrayView<const Scalar>& newVals,
                          BinaryFunction f) const
    {
      typedef typename Teuchos::ArrayView<Scalar>::size_type size_type;
      const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
      const size_type numElts = inds.size ();
      size_t hint = 0; // Guess for the current index k into rowVals

      // Get a view of the column indices in the row.  This amortizes
      // the cost of getting the view over all the entries of inds.
      Teuchos::ArrayView<const LocalOrdinal> colInds = getLocalView (rowInfo);

      LocalOrdinal numValid = 0; // number of valid local column indices
      for (size_type j = 0; j < numElts; ++j) {
        const size_t k = findLocalIndex (rowInfo, inds[j], colInds, hint);
        if (k != STINV) {
          rowVals[k] = f (rowVals[k], newVals[j]); // use binary function f
          hint = k+1;
          ++numValid;
        }
      }
      return numValid;
    }

    /// \brief Transform the given values using global indices.
    ///
    /// \param rowInfo [in] Information about a given row of the graph.
    ///
    /// \param rowVals [in/out] The values to be transformed.  They
    ///   correspond to the row indicated by rowInfo.
    ///
    /// \param inds [in] The (global) indices in the row, for which
    ///   to transform the corresponding values in rowVals.
    ///
    /// \param newVals [in] Values to use for transforming rowVals.
    ///   It's probably OK for these to alias rowVals.
    ///
    /// \param f [in] A binary function used to transform rowVals.
    ///
    /// \return The number of valid local column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    template<class Scalar, class BinaryFunction>
    LocalOrdinal
    transformGlobalValues (RowInfo rowInfo,
                           const Teuchos::ArrayView<Scalar>& rowVals,
                           const Teuchos::ArrayView<const GlobalOrdinal>& inds,
                           const Teuchos::ArrayView<const Scalar>& newVals,
                           BinaryFunction f) const
    {
      typedef typename Teuchos::ArrayView<Scalar>::size_type size_type;
      const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
      const size_type numElts = inds.size ();
      size_t hint = 0; // guess at the index's relative offset in the row

      LocalOrdinal numValid = 0; // number of valid local column indices
      for (size_type j = 0; j < numElts; ++j) {
        const size_t k = findGlobalIndex (rowInfo, inds[j], hint);
        if (k != STINV) {
          rowVals[k] = f (rowVals[k], newVals[j]); // use binary function f
          hint = k+1;
          numValid++;
        }
      }
      return numValid;
    }

    //@}
    //! \name Methods for sorting and merging column indices.
    //@{

    //! Whether duplicate column indices in each row have been merged.
    bool isMerged () const;

    /// \brief Report that we made a local modification to its structure.
    ///
    /// Call this after making a local change to the graph's
    /// structure.  Changing the structure locally invalidates the "is
    /// sorted" and "is merged" states.
    void setLocallyModified ();

    //! Sort the column indices in all the rows.
    void sortAllIndices ();

    //! Sort the column indices in the given row.
    void sortRowIndices (const RowInfo rowinfo);

    /// \brief Sort the column indices and their values in the given row.
    ///
    /// \tparam Scalar The type of the values.  When calling this
    ///   method from CrsMatrix, this should be the same as the
    ///   <tt>Scalar</tt> template parameter of CrsMatrix.
    ///
    /// \param rowinfo [in] Result of getRowInfo() for the row.
    ///
    /// \param values [in/out] On input: values for the given row.  If
    ///   indices is an array of the column indices in the row, then
    ///   values and indices should have the same number of entries,
    ///   and indices[k] should be the column index corresponding to
    ///   values[k].  On output: the same values, but sorted in the
    ///   same order as the (now sorted) column indices in the row.
    template <class Scalar>
    void sortRowIndicesAndValues (const RowInfo rowinfo,
                                  const Teuchos::ArrayView<Scalar>& values);

    /// \brief Merge duplicate row indices in all of the rows.
    ///
    /// \pre The graph is locally indexed:
    ///   <tt>isGloballyIndexed() == false</tt>.
    ///
    /// \pre The graph has not already been merged: <tt>isMerged()
    ///   == false</tt>.  That is, this function would normally only
    ///   be called after calling sortIndices().
    void mergeAllIndices ();

    /// \brief Merge duplicate row indices in the given row.
    ///
    /// \pre The graph is not already storage optimized:
    ///   <tt>isStorageOptimized() == false</tt>
    void mergeRowIndices (RowInfo rowinfo);

    /// \brief Merge duplicate row indices in the given row, along
    ///   with their corresponding values.
    ///
    /// This method is only called by CrsMatrix, for a CrsMatrix whose
    /// graph is this CrsGraph instance.  It is only called when the
    /// matrix owns the graph, not when the matrix was constructed
    /// with a const graph.
    ///
    /// \pre The graph is not already storage optimized:
    ///   <tt>isStorageOptimized() == false</tt>
    template<class Scalar>
    void
    mergeRowIndicesAndValues (RowInfo rowinfo,
                              const Teuchos::ArrayView<Scalar>& rowValues);
    //@}

    /// Set the domain and range Maps, and invalidate the Import
    /// and/or Export objects if necessary.
    ///
    /// If the domain Map has changed, invalidate the Import object
    /// (if there is one).  Likewise, if the range Map has changed,
    /// invalidate the Export object (if there is one).
    ///
    /// \param domainMap [in] The new domain Map
    /// \param rangeMap [in] The new range Map
    void
    setDomainRangeMaps (const Teuchos::RCP<const map_type>& domainMap,
                        const Teuchos::RCP<const map_type>& rangeMap);

    void staticAssertions() const;
    void clearGlobalConstants();
    void computeGlobalConstants();

    /// \brief Get information about the locally owned row with local
    ///   index myRow.
    RowInfo getRowInfo (const size_t myRow) const;

    /// \brief Get a const, nonowned, locally indexed view of the
    ///   locally owned row myRow, such that rowinfo =
    ///   getRowInfo(myRow).
    Teuchos::ArrayView<const LocalOrdinal>
    getLocalView (const RowInfo rowinfo) const;

    /// \brief Get a nonconst, nonowned, locally indexed view of the
    ///   locally owned row myRow, such that rowinfo =
    ///   getRowInfo(myRow).
    Teuchos::ArrayView<LocalOrdinal>
    getLocalViewNonConst (const RowInfo rowinfo);

    /// \brief Get a const, nonowned, globally indexed view of the
    ///   locally owned row myRow, such that rowinfo =
    ///   getRowInfo(myRow).
    Teuchos::ArrayView<const GlobalOrdinal>
    getGlobalView (const RowInfo rowinfo) const;

    /// \brief Get a nonconst, nonowned, globally indexed view of the
    ///   locally owned row myRow, such that rowinfo =
    ///   getRowInfo(myRow).
    Teuchos::ArrayView<GlobalOrdinal>
    getGlobalViewNonConst (const RowInfo rowinfo);

    /// \brief Find the column offset corresponding to the given
    ///   (local) column index.
    ///
    /// The name of this method is a bit misleading.  It does not
    /// actually find the column index.  Instead, it takes a local
    /// column index \c ind, and returns the corresponding offset
    /// into the raw array of column indices (whether that be 1-D or
    /// 2-D storage).
    ///
    /// \param rowinfo [in] Result of getRowInfo() for the given row.
    /// \param ind [in] (Local) column index for which to find the offset.
    /// \param hint [in] Hint for where to find \c ind in the column
    ///   indices for the given row.  If colInds is the ArrayView of
    ///   the (local) column indices for the given row, and if
    ///   <tt>colInds[hint] == ind</tt>, then the hint is correct.
    ///   The hint is ignored if it is out of range (that is,
    ///   greater than or equal to the number of entries in the
    ///   given row).
    ///
    /// The hint optimizes for the case of calling this method
    /// several times with the same row (as it would be in
    /// transformLocalValues) when several index inputs occur in
    /// consecutive sequence.  This may occur (for example) when
    /// there are multiple degrees of freedom per mesh point, and
    /// users are handling the assignment of degrees of freedom to
    /// global indices manually (rather than letting BlockMap take
    /// care of it).  In that case, users might choose to assign the
    /// degrees of freedom for a mesh point to consecutive global
    /// indices.  Epetra implements the hint for this reason.
    ///
    /// The hint only costs two comparisons (one to check range, and
    /// the other to see if the hint was correct), and it can save
    /// searching for the indices (which may take a lot more than
    /// two comparisons).
    size_t
    findLocalIndex (RowInfo rowinfo,
                    LocalOrdinal ind,
                    size_t hint = 0) const;

    /// Find the column offset corresponding to the given (local)
    /// column index, given a view of the (local) column indices.
    ///
    /// The name of this method is a bit misleading.  It does not
    /// actually find the column index.  Instead, it takes a local
    /// column index \c ind, and returns the corresponding offset
    /// into the raw array of column indices (whether that be 1-D or
    /// 2-D storage).
    ///
    /// It is best to use this method if you plan to call it several
    /// times for the same row, like in transformLocalValues().  In
    /// that case, it amortizes the overhead of calling
    /// getLocalView().
    ///
    /// \param rowinfo [in] Result of getRowInfo() for the given row.
    /// \param ind [in] (Local) column index for which to find the offset.
    /// \param colInds [in] View of all the (local) column indices
    ///   for the given row.
    /// \param hint [in] Hint for where to find \c ind in the column
    ///   indices for the given row.  If colInds is the ArrayView of
    ///   the (local) column indices for the given row, and if
    ///   <tt>colInds[hint] == ind</tt>, then the hint is correct.
    ///   The hint is ignored if it is out of range (that is,
    ///   greater than or equal to the number of entries in the
    ///   given row).
    ///
    /// See the documentation of the three-argument version of this
    /// method for an explanation and justification of the hint.
    size_t
    findLocalIndex (RowInfo rowinfo,
                    LocalOrdinal ind,
                    Teuchos::ArrayView<const LocalOrdinal> colInds,
                    size_t hint = 0) const;

    /// \brief Find the column offset corresponding to the given
    ///   (global) column index.
    ///
    /// The name of this method is a bit misleading.  It does not
    /// actually find the column index.  Instead, it takes a global
    /// column index \c ind, and returns the corresponding offset
    /// into the raw array of column indices (whether that be 1-D or
    /// 2-D storage).
    size_t findGlobalIndex (RowInfo rowinfo, GlobalOrdinal ind, size_t hint = 0) const;

    /// \brief Get the local graph.
    ///
    /// This is only a valid representation of the local graph if the
    /// (global) graph is fill complete.
    local_graph_type getLocalGraph () const;

    //! Get the local graph (DEPRECATED: call getLocalGraph() instead).
    TPETRA_DEPRECATED local_graph_type getLocalGraph_Kokkos () const;

    void fillLocalGraph (const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Whether it is correct to call getRowInfo().
    bool hasRowInfo () const;

    //! Throw an exception if the internal state is not consistent.
    void checkInternalState () const;

    //! The Map describing the distribution of rows of the graph.
    Teuchos::RCP<const map_type> rowMap_;
    //! The Map describing the distribution of columns of the graph.
    Teuchos::RCP<const map_type> colMap_;
    //! The Map describing the range of the (matrix corresponding to the) graph.
    Teuchos::RCP<const map_type> rangeMap_;
    //! The Map describing the domain of the (matrix corresponding to the) graph.
    Teuchos::RCP<const map_type> domainMap_;

    /// \brief The Import from the domain Map to the column Map.
    ///
    /// This gets constructed by fillComplete.  It may be null if
    /// the domain Map and the column Map are the same, since no
    /// Import is necessary in that case for sparse matrix-vector
    /// multiply.
    Teuchos::RCP<const import_type> importer_;

    /// \brief The Export from the row Map to the range Map.
    ///
    /// This gets constructed by fillComplete.  It may be null if
    /// the row Map and the range Map are the same, since no Export
    /// is necessary in that case for sparse matrix-vector multiply.
    Teuchos::RCP<const export_type> exporter_;

    //! Local graph; only initialized after first fillComplete() call.
    local_graph_type lclGraph_;

    // Local and Global Counts
    // nodeNumEntries_ and nodeNumAllocated_ are required to be always consistent
    // nodeMaxNumEntries_, nodeNumDiags_ and the global quantities are computed during fillComplete() and only valid when isFillComplete()
    global_size_t globalNumEntries_, globalNumDiags_, globalMaxNumRowEntries_;
    size_t          nodeNumEntries_,   nodeNumDiags_,   nodeMaxNumRowEntries_, nodeNumAllocated_;

    //! Whether the graph was allocated with static or dynamic profile.
    ProfileType pftype_;

    /// \brief The maximum number of entries to allow in each locally
    ///   owned row, per row.
    ///
    /// This is an argument to some of the graph's constructors.
    /// Either this or numAllocForAllRows_ is used, but not both.
    /// allocateIndices, setAllIndices, and expertStaticFillComplete
    /// all deallocate this array once they are done with it.
    ///
    /// This array <i>only</i> exists on a process before the graph's
    /// indices are allocated on that process.  After that point, it
    /// is discarded, since the graph's allocation implicitly or
    /// explicitly represents the same information.
    ///
    /// FIXME (mfh 07 Aug 2014) We want graph's constructors to
    /// allocate, rather than doing lazy allocation at first insert.
    /// This will make both k_numAllocPerRow_ and numAllocForAllRows_
    /// obsolete.
    Kokkos::DualView<const size_t*, Kokkos::LayoutLeft, execution_space> k_numAllocPerRow_;

    /// \brief The maximum number of entries to allow in each locally owned row.
    ///
    /// This is an argument to some of the graph's constructors.
    /// Either this or k_numAllocPerRow_ is used, but not both.
    ///
    /// FIXME (mfh 07 Aug 2014) We want graph's constructors to
    /// allocate, rather than doing lazy allocation at first insert.
    /// This will make both k_numAllocPerRow_ and numAllocForAllRows_
    /// obsolete.
    size_t numAllocForAllRows_;

    //! \name 1-D storage (StaticProfile) data structures
    //@{

    /// \brief Local column indices for all rows.
    ///
    /// This is only allocated if
    ///
    ///   - The calling process has a nonzero number of entries
    ///   - The graph has StaticProfile (1-D storage)
    ///   - The graph is locally indexed
    typename local_graph_type::entries_type::non_const_type k_lclInds1D_;

    //! Type of the k_gblInds1D_ array of global column indices.
    typedef Kokkos::View<GlobalOrdinal*, execution_space> t_GlobalOrdinal_1D;

    /// \brief Global column indices for all rows.
    ///
    /// This is only allocated if
    ///
    ///   - The calling process has a nonzero number of entries
    ///   - The graph has StaticProfile (1-D storage)
    ///   - The graph is globally indexed
    t_GlobalOrdinal_1D k_gblInds1D_;

    /// \brief Legacy Kokkos classic version of k_gblInds1D_.
    ///
    /// This is just a view of k_gblInds1D_.  We create views using
    /// Kokkos::Compat::persistingView, so the Kokkos::View won't get
    /// deallocated until the ArrayRCP's reference count goes to zero.
    Teuchos::ArrayRCP<GlobalOrdinal> gblInds1D_;

    /// \brief Row offsets for "1-D" storage.
    ///
    /// This is only allocated if "1-D" (StaticProfile) storage is
    /// active.  In that case, if beg = k_rowPtrs_(i_lcl) and end =
    /// k_rowPtrs_(i_lcl+1) for local row index i_lcl, then
    ///
    ///   - if the graph is locally indexed, k_lclInds1D_(beg:end-1)
    ///     (inclusive range) is the space for any local column
    ///     indices in local row i_lcl, else
    ///   - if the graph is globally indexed, k_gblInds1D_(beg:end-1)
    ///     (inclusive range) is the space for any global column
    ///     indices in local row i_lcl.
    ///
    /// Only the first k_numRowEntries_(i_lcl) of these entries are
    /// actual valid column indices.  Any remaining entries are "extra
    /// space."  If the graph's storage is packed, then there is no
    /// extra space, and the k_numRowEntries_ array is invalid.
    ///
    /// Both the k_rowPtrs_ and k_numRowEntries_ arrays are not
    /// allocated if the graph has 2-D (DynamicProfile) storage.
    ///
    /// If it is allocated, k_rowPtrs_ has length getNodeNumRows()+1.
    /// The k_numRowEntries_ array has has length getNodeNumRows(),
    /// again if it is allocated.
    ///
    /// [we may delete this to save memory on fillComplete, if "Delete
    /// Row Pointers" is specified.]
    typename local_graph_type::row_map_type::const_type k_rowPtrs_;

    //@}
    /// \name 2-D storage (DynamicProfile) data structures
    ///
    /// 2-D storage exists only if the graph was allocated with
    /// DynamicProfile.  All of these data structures exist in host
    /// memory.  Currently, NONE of them are thread safe, let alone
    /// thread scalable.  These data structures only exist to support
    /// legacy use cases.  At some point, we may add a thread-scalable
    /// intermediate level of "dynamicity" between 2-D storage and 1-D
    /// storage (StaticProfile), which bounds the <i>total</i> number
    /// of entries allowed per process, but does <i>not</i> otherwise
    /// bound the number of entries per row.
    //@{

    /// \brief Local column indices for all rows.
    ///
    /// This is only allocated if
    ///
    ///   - The calling process has a nonzero number of entries
    ///   - The graph has DynamicProfile (2-D storage)
    ///   - The graph is locally indexed
    ///
    /// In that case, if i_lcl is the local index of a locally owned
    /// row, then lclInds2D_[i_lcl] stores the local column indices
    /// for that row.
    Teuchos::ArrayRCP<Teuchos::Array<LocalOrdinal> > lclInds2D_;

    /// \brief Global column indices for all rows.
    ///
    /// This is only allocated if
    ///
    ///   - The calling process has a nonzero number of entries
    ///   - The graph has DynamicProfile (2-D storage)
    ///   - The graph is globally indexed
    ///
    /// In that case, if i_gbl is the global index of a globally owned
    /// row, then gblInds2D_[i_gbl] stores the global column indices
    /// for that row.
    Teuchos::ArrayRCP<Teuchos::Array<GlobalOrdinal> > gblInds2D_;

    typedef Kokkos::DualView<size_t*, Kokkos::LayoutLeft, execution_space> t_numRowEntries_;

    /// \brief The number of local entries in each locally owned row.
    ///
    /// This is deallocated in fillComplete() if fillComplete()'s
    /// "Optimize Storage" parameter is set to \c true.
    ///
    /// This may also exist with 1-D storage, if storage is unpacked.
    t_numRowEntries_ k_numRowEntries_;

    //@}

    /// \brief Status of the graph's storage, when not in a
    ///   fill-complete state.
    ///
    /// The phrase "When not in a fill-complete state" is important.
    /// When the graph is fill complete, it <i>always</i> uses 1-D
    /// "packed" storage.  However, if the "Optimize Storage"
    /// parameter to fillComplete was false, the graph may keep
    /// unpacked 1-D or 2-D storage around and resume it on the next
    /// resumeFill call.
    Details::EStorageStatus storageStatus_;

    bool indicesAreAllocated_;
    bool indicesAreLocal_;
    bool indicesAreGlobal_;
    bool fillComplete_;

    //! Whether the graph is locally lower triangular.
    bool lowerTriangular_;
    //! Whether the graph is locally upper triangular.
    bool upperTriangular_;
    //! Whether the graph's indices are sorted in each row, on this process.
    bool indicesAreSorted_;
    /// \brief Whether the graph's indices are non-redundant (merged)
    ///   in each row, on this process.
    bool noRedundancies_;
    //! Whether this process has computed local constants.
    bool haveLocalConstants_;
    //! Whether all processes have computed global constants.
    bool haveGlobalConstants_;

    //! Nonlocal data given to insertGlobalValues or sumIntoGlobalValues.
    std::map<GlobalOrdinal, std::vector<GlobalOrdinal> > nonlocals_;

    /// \brief Whether to require makeColMap() (and therefore
    ///   fillComplete()) to order column Map GIDs associated with
    ///   each remote process in ascending order.
    ///
    /// makeColMap() always groups remote GIDs by process rank, so
    /// that all remote GIDs with the same owning rank occur
    /// contiguously.  By default, it always sorts remote GIDs in
    /// increasing order within those groups.  This behavior differs
    /// from Epetra, which does not sort remote GIDs with the same
    /// owning process.
    ///
    /// This is \c true by default, which means "sort remote GIDs."
    /// If you don't want to sort (for compatibility with Epetra),
    /// call sortGhostColumnGIDsWithinProcessBlock(false).
    bool sortGhostsAssociatedWithEachProcessor_;

  }; // class CrsGraph

} // namespace Tpetra

#endif
