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
// ************************************************************************
// @HEADER

#ifndef TPETRA_CRSGRAPH_DECL_HPP
#define TPETRA_CRSGRAPH_DECL_HPP

/// \file Tpetra_CrsGraph_decl.hpp
/// \brief Declaration of the Tpetra::CrsGraph class
///
/// If you want to use Tpetra::CrsGraph, include "Tpetra_CrsGraph.hpp"
/// (a file which CMake generates and installs for you).  If you only
/// want the declaration of Tpetra::CrsGraph, include this file
/// (Tpetra_CrsGraph_decl.hpp).

#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_BlockCrsMatrix_fwd.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Exceptions.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_Util.hpp" // need this here for sort2

#include "KokkosSparse_findRelOffset.hpp"
#include "Kokkos_DualView.hpp"
#include "Kokkos_StaticCrsGraph.hpp"

#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include <functional> // std::function
#include <memory>

namespace Tpetra {


  // Forward declaration for CrsGraph::swap() test
  template<class LocalOrdinal, class GlobalOrdinal, class Node> class crsGraph_Swap_Tester;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  namespace Details {
    template<class LocalOrdinal,
             class GlobalOrdinal>
    class CrsPadding;
  } // namespace Details

  namespace { // (anonymous)

    template<class ViewType>
    struct UnmanagedView {
      static_assert (Kokkos::is_view<ViewType>::value,
                     "ViewType must be a Kokkos::View specialization.");
      // FIXME (mfh 02 Dec 2015) Right now, this strips away other
      // memory traits.  Christian will add an "AllTraits" enum which is
      // the enum value of MemoryTraits<T>, that will help us fix this.
      typedef Kokkos::View<typename ViewType::data_type,
                           typename ViewType::array_layout,
                           typename ViewType::device_type,
                           Kokkos::MemoryUnmanaged> type;
    };

  } // namespace (anonymous)
#endif // DOXYGEN_SHOULD_SKIP_THIS

  /// \struct RowInfo
  /// \brief Allocation information for a locally owned row in a
  ///   CrsGraph or CrsMatrix
  ///
  /// A RowInfo instance identifies a locally owned row uniquely by
  /// its local index, and contains other information useful for
  /// inserting entries into the row.  It is the return value of
  /// CrsGraph's getRowInfo() or updateAllocAndValues() methods.
  struct RowInfo {
    size_t localRow;
    size_t allocSize;
    size_t numEntries;
    size_t offset1D;
  };

  enum ELocalGlobal {
    LocalIndices,
    GlobalIndices
  };

  namespace Details {
    /// \brief Status of the graph's or matrix's storage, when not in
    ///   a fill-complete state.
    ///
    /// When a CrsGraph or CrsMatrix is <i>not</i> fill complete and
    /// is allocated, then its data live in one of two storage
    /// formats:
    ///
    /// <ol>
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
    /// The phrase "When not in a fill-complete state" is important.
    /// When the graph is fill complete, it <i>always</i> uses 1-D
    /// "packed" storage.  However, if storage is "not optimized," we
    /// retain the 1-D unpacked format, and thus retain this enum
    /// value.
    enum EStorageStatus {
      STORAGE_1D_UNPACKED, //<! 1-D "unpacked" storage
      STORAGE_1D_PACKED, //<! 1-D "packed" storage
      STORAGE_UB //<! Invalid value; upper bound on enum values
    };

  } // namespace Details

  /// \class CrsGraph
  /// \brief A distributed graph accessed by rows (adjacency lists)
  ///   and stored sparsely.
  ///
  /// \tparam LocalOrdinal The type of local indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam GlobalOrdinal The type of global indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam Node The Kokkos Node type.  See the documentation of Map
  ///   for requirements.
  ///
  /// This class implements a distributed-memory parallel sparse
  /// graph.  It provides access by rows to the elements of the graph,
  /// as if the local data were stored in compressed sparse row format
  /// (adjacency lists, in graph terms).  (Implementations are
  /// <i>not</i> required to store the data in this way internally.)
  /// This class has an interface like that of Epetra_CrsGraph, but
  /// also allows insertion of data into nonowned rows, much like
  /// Epetra_FECrsGraph.
  ///
  /// \section Tpetra_CrsGraph_prereq Prerequisites
  ///
  /// Before reading the rest of this documentation, it helps to know
  /// something about the Teuchos memory management classes, in
  /// particular Teuchos::RCP, Teuchos::ArrayRCP, and
  /// Teuchos::ArrayView.  You should also know a little bit about MPI
  /// (the Message Passing Interface for distributed-memory
  /// programming).  You won't have to use MPI directly to use
  /// CrsGraph, but it helps to be familiar with the general idea of
  /// distributed storage of data over a communicator.  Finally, you
  /// should read the documentation of Map.
  ///
  /// \section Tpetra_CrsGraph_local_vs_global Local vs. global indices and nonlocal insertion
  ///
  /// Graph entries can be added using either local or global coordinates
  /// for the indices. The accessors isGloballyIndexed() and
  /// isLocallyIndexed() indicate whether the indices are currently
  /// stored as global or local indices. Many of the class methods are
  /// divided into global and local versions, which differ only in
  /// whether they accept/return indices in the global or local
  /// coordinate space. Some of these methods may only be used if the
  /// graph coordinates are in the appropriate coordinates.  For example,
  /// getGlobalRowView() returns a View to the indices in global
  /// coordinates; if the indices are not in global coordinates, then no
  /// such View can be created.
  ///
  /// The global/local distinction does distinguish between operation
  /// on the global/local graph. Almost all methods operate on the
  /// local graph, i.e., the rows of the graph associated with the
  /// local node, per the distribution specified by the row
  /// map. Access to non-local rows requires performing an explicit
  /// communication via the import/export capabilities of the CrsGraph
  /// object; see DistObject. However, the method
  /// insertGlobalIndices() is an exception to this rule, as non-local
  /// rows are allowed to be added via the local graph. These rows are
  /// stored in the local graph and communicated to the appropriate
  /// node on the next call to globalAssemble() or fillComplete() (the
  /// latter calls the former).
  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class CrsGraph :
    public RowGraph<LocalOrdinal, GlobalOrdinal, Node>,
    public DistObject<GlobalOrdinal,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node>,
    public Teuchos::ParameterListAcceptorDefaultBase
  {
    template <class S, class LO, class GO, class N>
    friend class CrsMatrix;
    template <class LO2, class GO2, class N2>
    friend class CrsGraph;
    template <class LO, class GO, class N>
    friend class FECrsGraph;

    //! The specialization of DistObject that is this class' parent class.
    using dist_object_type = DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;

  public:
    //! The type of the graph's local indices.
    using local_ordinal_type = LocalOrdinal;
    //! The type of the graph's global indices.
    using global_ordinal_type = GlobalOrdinal;
    //! This class' Kokkos device type.
    using device_type = typename Node::device_type;
    //! This class' Kokkos execution space.
    using execution_space = typename device_type::execution_space;

    /// \brief This class' Kokkos Node type.
    ///
    /// This is a leftover that will be deprecated and removed.
    /// See e.g., GitHub Issue #57.
    using node_type = Node;

    //! The type of the part of the sparse graph on each MPI process.
    using local_graph_type = Kokkos::StaticCrsGraph<local_ordinal_type,
                                                    Kokkos::LayoutLeft,
                                                    device_type,
                                                    void,
                                                    size_t>;

    //! The Map specialization used by this class.
    using map_type = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    //! The Import specialization used by this class.
    using import_type = ::Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
    //! The Export specialization used by this class.
    using export_type = ::Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>;

    //! @name Constructor/Destructor Methods
    //@{

    /// \brief Constructor specifying a single upper bound for the
    ///   number of entries in all rows on the calling process.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  This is a strict upper bound.
    ///
    /// \param pftype [in] If you specify this, then this must always
    ///   be StaticProfile.  No other values exist or are permitted.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const size_t maxNumEntriesPerRow,
              const ProfileType pftype = TPETRA_DEFAULT_PROFILE_TYPE,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a (possibly different) upper
    ///   bound for the number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  This is a strict upper bound.
    ///
    /// \param pftype [in] If you specify this, then this must always
    ///   be StaticProfile.  No other values exist or are permitted.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
              const ProfileType pftype = TPETRA_DEFAULT_PROFILE_TYPE,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a (possibly different) upper
    ///   bound for the number of entries in each row (legacy
    ///   KokkosClassic version).
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  This is a strict upper bound.
    ///
    /// \param pftype [in] If you specify this, then this must always
    ///   be StaticProfile.  No other values exist or are permitted.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::ArrayView<const size_t>& numEntPerRow,
              const ProfileType pftype = TPETRA_DEFAULT_PROFILE_TYPE,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


    /// \brief Constructor specifying column Map and a single upper
    ///   bound for the number of entries in all rows on the calling
    ///   process.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  This is a strict upper bound.
    ///
    /// \param pftype [in] If you specify this, then this must always
    ///   be StaticProfile.  No other values exist or are permitted.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const size_t maxNumEntriesPerRow,
              const ProfileType pftype = TPETRA_DEFAULT_PROFILE_TYPE,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  This is a strict upper bound.
    ///
    /// \param pftype [in] If you specify this, then this must always
    ///   be StaticProfile.  No other values exist or are permitted.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
              const ProfileType pftype = TPETRA_DEFAULT_PROFILE_TYPE,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and number of entries
    ///   in each row (legacy KokkosClassic version).
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  This is a strict upper bound.
    ///
    /// \param pftype [in] If you specify this, then this must always
    ///   be StaticProfile.  No other values exist or are permitted.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const Teuchos::ArrayView<const size_t>& numEntPerRow,
              const ProfileType pftype = TPETRA_DEFAULT_PROFILE_TYPE,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


    /// \brief Constructor specifying column Map and arrays containing
    ///   the graph. In almost all cases the indices must be sorted on input,
    ///   but if they aren't sorted, "sorted" must be set to false in params.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param rowPointers [in] The beginning of each row in the graph,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the graph.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the graph.
    ///   Entries in each row must be sorted (by local index).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const typename local_graph_type::row_map_type& rowPointers,
              const typename local_graph_type::entries_type::non_const_type& columnIndices,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and arrays containing
    ///   the graph. In almost all cases the indices must be sorted on input,
    ///   but if they aren't sorted, "sorted" must be set to false in params.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param rowPointers [in] The beginning of each row in the graph,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the graph.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the graph.
    ///   Entries in each row must be sorted (by local index).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const Teuchos::ArrayRCP<size_t>& rowPointers,
              const Teuchos::ArrayRCP<local_ordinal_type>& columnIndices,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and a local
    ///   graph, which the resulting CrsGraph views.
    ///   In almost all cases the local graph must be sorted on input,
    ///   but if it isn't sorted, "sorted" must be set to false in params.
    ///
    /// Unlike most other CrsGraph constructors, successful completion
    /// of this constructor will result in a fill-complete graph.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///   See replaceColMap() for the requirements.
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

    /// \brief Constructor specifying column, domain and range maps, and a
    ///   local graph, which the resulting CrsGraph views.
    ///   In almost all cases the local graph must be sorted on input,
    ///   but if it isn't sorted, "sorted" must be set to false in params.
    ///
    /// Unlike most other CrsGraph constructors, successful completion
    /// of this constructor will result in a fill-complete graph.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param domainMap [in] The graph's domain Map. MUST be one to
    ///   one!
    ///
    /// \param rangeMap [in] The graph's range Map.  MUST be one to
    ///   one!  May be, but need not be, the same as the domain Map.
    ///
    /// \param lclGraph [in] A locally indexed Kokkos::StaticCrsGraph
    ///   whose local row indices come from the specified row Map, and
    ///   whose local column indices come from the specified column
    ///   Map.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const local_graph_type& lclGraph,
              const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
              const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Create a fill-complete CrsGraph from all the things it needs.
    /// \param lclGraph [in] The local graph.  In almost all cases the
    ///   local graph must be sorted on input,
    ///   but if it isn't sorted, "sorted" must be set to false in params.
    CrsGraph (const local_graph_type& lclGraph,
              const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const Teuchos::RCP<const map_type>& domainMap,
              const Teuchos::RCP<const map_type>& rangeMap,
              const Teuchos::RCP<const import_type>& importer,
              const Teuchos::RCP<const export_type>& exporter,
              const Teuchos::RCP<Teuchos::ParameterList>& params =
                Teuchos::null);

    //! Copy constructor (default).
    CrsGraph (const CrsGraph<local_ordinal_type, global_ordinal_type, node_type>&) = default;

    //! Assignment operator (default).
    CrsGraph& operator= (const CrsGraph<local_ordinal_type, global_ordinal_type, node_type>&) = default;

    //! Move constructor (default).
    CrsGraph (CrsGraph<local_ordinal_type, global_ordinal_type, node_type>&&) = default;

    //! Move assignment (default).
    CrsGraph& operator= (CrsGraph<local_ordinal_type, global_ordinal_type, node_type>&&) = default;

    /// \brief Destructor (virtual for memory safety of derived classes).
    ///
    /// \note To Tpetra developers: See the C++ Core Guidelines C.21
    ///   ("If you define or <tt>=delete</tt> any default operation,
    ///   define or <tt>=delete</tt> them all"), in particular the
    ///   AbstractBase example, for why this destructor declaration
    ///   implies that we need the above four <tt>=default</tt>
    ///   declarations for copy construction, move construction, copy
    ///   assignment, and move assignment.
    virtual ~CrsGraph () = default;

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
    /// \brief True if and only if \c CrsGraph is identical to this CrsGraph
    ///
    /// \warning THIS METHOD IS FOR TPETRA DEVELOPERS ONLY.  DO NOT
    ///   RELY ON THIS METHOD.  WE MAKE NO PROMISES OF BACKWARDS
    ///   COMPATIBILITY.
    ///
    /// This performs _exact_ matches on objects with in the graphs. That is,
    /// internal data structures such as arrays must match exactly in both
    /// content and order. This is not performing any sort of isomorphic
    /// search.
    ///
    /// \param graph [in] a CrsGraph to compare against this one.
    ///
    /// \return True if the other CrsGraph's data structure is identical to this
    ///         CrsGraph.
    ///
    bool isIdenticalTo(const CrsGraph<LocalOrdinal, GlobalOrdinal, Node> &graph) const;

    //@}
    //! @name Implementation of Teuchos::ParameterListAcceptor
    //@{

    //! Set the given list of parameters (must be nonnull).
    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& params) override;

    //! Default parameter list suitable for validation.
    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const override;

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
    insertGlobalIndices (const global_ordinal_type globalRow,
                         const Teuchos::ArrayView<const global_ordinal_type>& indices);

    /// \brief Epetra compatibility version of insertGlobalIndices
    ///   (see above) that takes input as a raw pointer, rather than
    ///   Teuchos::ArrayView.
    ///
    /// Arguments are the same and in the same order as
    /// Epetra_CrsGraph::InsertGlobalIndices.
    void
    insertGlobalIndices (const global_ordinal_type globalRow,
                         const local_ordinal_type numEnt,
                         const global_ordinal_type inds[]);

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
    insertLocalIndices (const local_ordinal_type localRow,
                        const Teuchos::ArrayView<const local_ordinal_type>& indices);

    /// \brief Epetra compatibility version of insertLocalIndices
    ///   (see above) that takes input as a raw pointer, rather than
    ///   Teuchos::ArrayView.
    ///
    /// Arguments are the same and in the same order as
    /// Epetra_CrsGraph::InsertMyIndices.
    void
    insertLocalIndices (const local_ordinal_type localRow,
                        const local_ordinal_type numEnt,
                        const local_ordinal_type inds[]);

    //! Remove all graph indices from the specified local row.
    /**
       \pre \c localRow is a local row of this graph.
       \pre <tt>isGloballyIndexed() == false</tt>
       \pre <tt>isStorageOptimized() == false</tt>

       \post <tt>getNumEntriesInLocalRow(localRow) == 0</tt>
       \post <tt>indicesAreAllocated() == true</tt>
       \post <tt>isLocallyIndexed() == true</tt>
    */
    void removeLocalIndices (local_ordinal_type localRow);

    //@}
    //! @name Collective methods for changing the graph's global state
    //@{

    /// \brief Communicate nonlocal contributions to other processes.
    ///
    /// This method is called automatically by fillComplete().  Most
    /// users do not need to call this themselves.
    ///
    /// This method must be called collectively (that is, like any MPI
    /// collective) over all processes in the graph's communicator.
    void globalAssemble ();

    /// \brief Resume fill operations.
    ///
    /// After calling fillComplete(), resumeFill() must be called
    /// before initiating any changes to the graph.
    ///
    /// resumeFill() may be called repeatedly.
    ///
    /// \warning A CrsGraph instance does not currently (as of 23 Jul
    ///   2017) and never did support arbitrary structure changes
    ///   after the first fillComplete call on that instance.  The
    ///   safest thing to do is not to change structure at all after
    ///   first fillComplete.
    ///
    /// \post <tt>isFillActive() == true<tt>
    /// \post <tt>isFillComplete() == false<tt>
    ///
    /// This method must be called collectively (that is, like any MPI
    /// collective) over all processes in the graph's communicator.
    void
    resumeFill (const Teuchos::RCP<Teuchos::ParameterList>& params =
                  Teuchos::null);

    /// \brief Tell the graph that you are done changing its structure.
    ///
    /// This tells the graph to optimize its data structures for
    /// computational kernels, and to prepare (MPI) communication
    /// patterns.
    ///
    /// Off-process indices are distributed (via globalAssemble()),
    /// indices are sorted, redundant indices are eliminated, and
    /// global indices are transformed to local indices.
    ///
    /// This method must be called collectively (that is, like any MPI
    /// collective) over all processes in the graph's communicator.
    ///
    /// \warning The domain Map and row Map arguments to this method
    ///   MUST be one to one!  If you have Maps that are not one to
    ///   one, and you do not know how to make a Map that covers the
    ///   same global indices but <i>is</i> one to one, then you may
    ///   call Tpetra::createOneToOne() (see Map's header file) to
    ///   make a one-to-one version of your Map.
    ///
    /// \pre  <tt>   isFillActive() && ! isFillComplete() </tt>
    /// \post <tt> ! isFillActive() &&   isFillComplete() </tt>
    ///
    /// \param domainMap [in] The graph's domain Map.  MUST be one to
    ///   one!
    /// \param rangeMap [in] The graph's range Map.  MUST be one to
    ///   one!  May be, but need not be, the same as the domain Map.
    /// \param params [in/out] List of parameters controlling this
    ///   method's behavior.  See below for valid parameters.
    ///
    /// List of valid parameters in <tt>params</tt>:
    /// <ul>
    /// <li> "Optimize Storage" (\c bool): Default is false.  If true,
    ///      then isStorageOptimized() returns true after this method
    ///      finishes.  See isStorageOptimized() for consequences.
    /// </li>
    /// </ul>
    void
    fillComplete (const Teuchos::RCP<const map_type>& domainMap,
                  const Teuchos::RCP<const map_type>& rangeMap,
                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Tell the graph that you are done changing its
    ///   structure; set default domain and range Maps.
    ///
    /// See above three-argument version of fillComplete for full
    /// documentation.  If the graph does not yet have domain and
    /// range Maps (i.e., if fillComplete has not yet been called on
    /// this graph at least once), then this method uses the graph's
    /// row Map (result of this->getRowMap()) as both the domain Map
    /// and the range Map.  Otherwise, this method uses the graph's
    /// existing domain and range Maps.
    ///
    /// This method must be called collectively (that is, like any MPI
    /// collective) over all processes in the graph's communicator.
    ///
    /// \warning It is only valid to call this overload of
    ///   fillComplete if the row Map is one to one!  If the row Map
    ///   is NOT one to one, you must call the above three-argument
    ///   version of fillComplete, and supply one-to-one domain and
    ///   range Maps.  If you have Maps that are not one to one, and
    ///   you do not know how to make a Map that covers the same
    ///   global indices but <i>is</i> one to one, then you may call
    ///   Tpetra::createOneToOne() (see Map's header file) to make a
    ///   one-to-one version of your Map.
    ///
    /// \param params [in/out] List of parameters controlling this
    ///   method's behavior.  See documentation of the three-argument
    ///   version of fillComplete (above) for valid parameters.
    void
    fillComplete (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Perform a fillComplete on a graph that already has
    ///   data, via setAllIndices().
    ///
    /// The graph must already have filled local 1-D storage.  If the
    /// graph has been constructed in any other way, this method will
    /// throw an exception.  This routine is needed to support other
    /// Trilinos packages and should not be called by ordinary users.
    ///
    /// This method must be called collectively (that is, like any MPI
    /// collective) over all processes in the graph's communicator.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    ///
    /// \param domainMap [in] The graph's domain Map.  MUST be one to
    ///   one!
    /// \param rangeMap [in] The graph's range Map.  MUST be one to
    ///   one!  May be, but need not be, the same as the domain Map.
    /// \param importer [in] Import from the graph's domain Map to its
    ///   column Map.  If no Import is necessary (i.e., if the domain
    ///   and column Maps are the same, in the sense of
    ///   Tpetra::Map::isSameAs), then this may be Teuchos::null.
    /// \param exporter [in] Export from the graph's row Map to its
    ///   range Map.  If no Export is necessary (i.e., if the row and
    ///   range Maps are the same, in the sense of
    ///   Tpetra::Map::isSameAs), then this may be Teuchos::null.
    /// \param params [in/out] List of parameters controlling this
    ///   method's behavior.
    void
    expertStaticFillComplete (const Teuchos::RCP<const map_type>& domainMap,
                              const Teuchos::RCP<const map_type>& rangeMap,
                              const Teuchos::RCP<const import_type>& importer =
                                Teuchos::null,
                              const Teuchos::RCP<const export_type>& exporter =
                                Teuchos::null,
                              const Teuchos::RCP<Teuchos::ParameterList>& params =
                                Teuchos::null);
    //@}
    //! @name Methods implementing RowGraph.
    //@{

    //! Returns the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const override;


    //! Returns the Map that describes the row distribution in this graph.
    Teuchos::RCP<const map_type> getRowMap () const override;

    //! \brief Returns the Map that describes the column distribution in this graph.
    Teuchos::RCP<const map_type> getColMap () const override;

    //! Returns the Map associated with the domain of this graph.
    Teuchos::RCP<const map_type> getDomainMap () const override;

    //! Returns the Map associated with the domain of this graph.
    Teuchos::RCP<const map_type> getRangeMap () const override;

    //! Returns the importer associated with this graph.
    Teuchos::RCP<const import_type> getImporter () const override;

    //! Returns the exporter associated with this graph.
    Teuchos::RCP<const export_type> getExporter () const override;

    //! Returns the number of global rows in the graph.
    /** Undefined if isFillActive().
     */
    global_size_t getGlobalNumRows() const override;

    //! \brief Returns the number of global columns in the graph.
    /** Returns the number of entries in the domain map of the matrix.
        Undefined if isFillActive().
    */
    global_size_t getGlobalNumCols () const override;

    //! Returns the number of graph rows owned on the calling node.
    size_t getNodeNumRows () const override;

    //! Returns the number of columns connected to the locally owned rows of this graph.
    /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
     */
    size_t getNodeNumCols () const override;

    //! Returns the index base for global indices for this graph.
    global_ordinal_type getIndexBase () const override;

    //! Returns the global number of entries in the graph.
    /** Undefined if isFillActive ().
     */
    global_size_t getGlobalNumEntries () const override;

    /// \brief The local number of entries in the graph.
    ///
    /// "Local" means "local to the calling (MPI) process."
    ///
    /// \warning If the graph is not fill complete, this may launch a
    ///   thread-parallel computational kernel.  This is because we do
    ///   not store the number of entries as a separate integer field,
    ///   since doing so and keeping it updated would hinder
    ///   thread-parallel insertion of new entries.  See #1357.
    size_t getNodeNumEntries() const override;

    //! \brief Returns the current number of entries on this node in the specified global row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
    size_t
    getNumEntriesInGlobalRow (global_ordinal_type globalRow) const override;

    /// \brief Get the number of entries in the given row (local index).
    ///
    /// \return The number of entries in the given row, specified by
    ///   local index, on the calling MPI process.  If the specified
    ///   local row index is invalid on the calling process, return
    ///   <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>.
    size_t
    getNumEntriesInLocalRow (local_ordinal_type localRow) const override;

    /// \brief The local number of indices allocated for the graph,
    ///   over all rows on the calling (MPI) process.
    ///
    /// "Local" means "local to the calling (MPI) process."
    ///
    /// \warning If the graph is not fill complete, this may require
    ///   computation.  This is because we do not store the allocation
    ///   count as a separate integer field, since doing so and
    ///   keeping it updated would hinder thread-parallel insertion of
    ///   new entries.
    ///
    /// This is the allocation available to the user. Actual
    /// allocation may be larger, for example, after calling
    /// fillComplete().  Thus, this does not necessarily reflect the
    /// graph's memory consumption.
    ///
    /// \return If indicesAreAllocated() is true, the allocation size.
    ///   Otherwise,
    ///   <tt>Tpetra::Details::OrdinalTraits<size_t>::invalid()</tt>.
    size_t getNodeAllocationSize () const;

    /// \brief Current number of allocated entries in the given row on
    ///   the calling (MPI) process, using a global row index.
    ///
    /// \return If the given row index is in the row Map on the
    ///   calling process, then return this process' allocation size
    ///   for that row.  Otherwise, return
    ///   <tt>Tpetra::Details::OrdinalTraits<size_t>::invalid()</tt>.
    size_t getNumAllocatedEntriesInGlobalRow (global_ordinal_type globalRow) const;

    /// \brief Current number of allocated entries in the given row on
    ///   the calling (MPI) process, using a local row index.
    ///
    /// \return If the given row index is in the row Map on the
    ///   calling process, then return this process' allocation size
    ///   for that row.  Otherwise, return
    ///   <tt>Tpetra::Details::OrdinalTraits<size_t>::invalid()</tt>.
    size_t getNumAllocatedEntriesInLocalRow (local_ordinal_type localRow) const;

    /// \brief Maximum number of entries in any row of the graph,
    ///   over all processes in the graph's communicator.
    ///
    /// \pre <tt>! isFillActive()</tt>
    ///
    /// \note This is the same as the result of a global maximum of
    ///   getNodeMaxNumRowEntries() over all processes.  That may not
    ///   necessarily mean what you think it does if some rows of the
    ///   matrix are owned by multiple processes.  In particular, some
    ///   processes might only own some of the entries in a particular
    ///   row.  This method only counts the number of entries in each
    ///   row that a process owns, not the total number of entries in
    ///   the row over all processes.
    size_t getGlobalMaxNumRowEntries () const override;

    /// \brief Maximum number of entries in any row of the graph,
    ///   on this process.
    ///
    /// \pre <tt>! isFillActive()</tt>
    size_t getNodeMaxNumRowEntries () const override;

    /// \brief Whether the graph has a column Map.
    ///
    /// A CrsGraph has a column Map either because it was given to its
    /// constructor, or because it was constructed in fillComplete().
    /// Calling fillComplete() always makes a column Map if the graph
    /// does not already have one.
    ///
    /// A column Map lets the graph
    /// <ul>
    /// <li> use local indices for storing entries in each row, and </li>
    /// <li> compute an Import from the domain Map to the column Map. </li>
    /// </ul>
    ///
    /// The latter is mainly useful for a graph associated with a
    /// CrsMatrix.
    bool hasColMap () const override;


    /// \brief Whether the graph's column indices are stored as local indices.
    ///
    /// Weird quirk inherited from Epetra:
    /// <tt>! isLocallyIndexed() && ! isGloballyIndexed()</tt>
    /// means that there are no graph entries on the calling process.
    /// Please don't rely on this behavior, but note that it's
    /// possible.
    bool isLocallyIndexed () const override;

    /// \brief Whether the graph's column indices are stored as global indices.
    ///
    /// Weird quirk inherited from Epetra:
    /// <tt>! isLocallyIndexed() && ! isGloballyIndexed()</tt>
    /// means that there are no graph entries on the calling process.
    /// Please don't rely on this behavior, but note that it's
    /// possible.
    bool isGloballyIndexed () const override;

    //! Whether fillComplete() has been called and the graph is in compute mode.
    bool isFillComplete () const override;

    //! Whether resumeFill() has been called and the graph is in edit mode.
    bool isFillActive () const;

    /// \brief Whether graph indices in all rows are known to be sorted.
    ///
    /// A fill-complete graph is always sorted, as is a newly
    /// constructed graph. A graph is sorted immediately after calling
    /// resumeFill(), but any changes to the graph may result in the
    /// sorting status becoming unknown (and therefore, presumed
    /// unsorted).
    bool isSorted () const;

    //! \brief Returns \c true if storage has been optimized.
    /**
       Optimized storage means that the allocation of each row is equal to the
       number of entries. The effect is that a pass through the matrix, i.e.,
       during a mat-vec, requires minimal memory traffic. One limitation of
       optimized storage is that no new indices can be added to the graph.
    */
    bool isStorageOptimized () const;

    //! Returns \c true if the graph was allocated with static data structures.
    ProfileType getProfileType () const;

    /// \brief Get a copy of the given row, using global indices.
    ///
    /// \param gblRow [in] Global index of the row.
    /// \param gblColInds [out] On output: Global column indices.
    /// \param numColInds [out] Number of indices returned.
    void
    getGlobalRowCopy (global_ordinal_type gblRow,
                      const Teuchos::ArrayView<global_ordinal_type>& gblColInds,
                      size_t& numColInds) const override;

    /// \brief Get a copy of the given row, using local indices.
    ///
    /// \param lclRow [in] Local index of the row.
    /// \param lclColInds [out] On output: Local column indices.
    /// \param numColInds [out] Number of indices returned.
    ///
    /// \pre <tt>hasColMap()</tt>
    void
    getLocalRowCopy (local_ordinal_type lclRow,
                     const Teuchos::ArrayView<local_ordinal_type>& lclColInds,
                     size_t& numColInds) const override;

    /// \brief Get a const, non-persisting view of the given global
    ///   row's global column indices, as a Teuchos::ArrayView.
    ///
    /// \param gblRow [in] Global index of the row.
    /// \param gblColInds [out] Global column indices in the row.  If
    ///   the given row is not a valid row index on the calling
    ///   process, then the result has no entries (its size is zero).
    ///
    /// \pre <tt>! isLocallyIndexed()</tt>
    /// \post <tt>gblColInds.size() == getNumEntriesInGlobalRow(gblRow)</tt>
    void
    getGlobalRowView (const global_ordinal_type gblRow,
                      Teuchos::ArrayView<const global_ordinal_type>& gblColInds) const override;

    /// \brief Whether this class implements getLocalRowView() and
    ///   getGlobalRowView() (it does).
    bool supportsRowViews () const override;

    /// \brief Get a const, non-persisting view of the given local
    ///   row's local column indices, as a Teuchos::ArrayView.
    ///
    /// \param lclRow [in] Local index of the row.
    /// \param lclColInds [out] Local column indices in the row.  If
    ///   the given row is not a valid row index on the calling
    ///   process, then the result has no entries (its size is zero).
    ///
    /// \pre <tt>! isGloballyIndexed()</tt>
    /// \post <tt>lclColInds.size() == getNumEntriesInLocalRow(lclRow)</tt>
    void
    getLocalRowView (const local_ordinal_type lclRow,
                     Teuchos::ArrayView<const local_ordinal_type>& lclColInds) const override;

    //@}
    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a one-line human-readable description of this object.
    std::string description () const override;

    /// \brief Print this object to the given output stream with the
    ///   given verbosity level.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const override;

    //@}
    //! \name Implementation of DistObject
    //@{

    /// \typedef buffer_device_type
    /// \brief Kokkos::Device specialization for communication buffers.
    ///
    /// See #1088 for why this is not just <tt>device_type::device_type</tt>.
    using buffer_device_type = typename dist_object_type::buffer_device_type;
    //! Type of each entry of the DistObject communication buffer.
    using packet_type = global_ordinal_type;

    virtual bool
    checkSizes (const SrcDistObject& source) override;

    virtual void
    copyAndPermute
    (const SrcDistObject& source,
     const size_t numSameIDs,
     const Kokkos::DualView<const local_ordinal_type*,
       buffer_device_type>& permuteToLIDs,
     const Kokkos::DualView<const local_ordinal_type*,
       buffer_device_type>& permuteFromLIDs) override;

    using padding_type = Details::CrsPadding<
      local_ordinal_type, global_ordinal_type>;

    void
    applyCrsPadding(const padding_type& padding,
                    const bool verbose);

    std::unique_ptr<padding_type>
    computeCrsPadding(
      const RowGraph<local_ordinal_type, global_ordinal_type,
        node_type>& source,
      const size_t numSameIDs,
      const Kokkos::DualView<const local_ordinal_type*,
        buffer_device_type>& permuteToLIDs,
      const Kokkos::DualView<const local_ordinal_type*,
        buffer_device_type>& permuteFromLIDs,
      const bool verbose) const;

    // This actually modifies imports by sorting it.
    std::unique_ptr<padding_type>
    computeCrsPaddingForImports(
      const Kokkos::DualView<const local_ordinal_type*,
        buffer_device_type>& importLIDs,
      Kokkos::DualView<packet_type*, buffer_device_type> imports,
      Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
      const bool verbose) const;

    std::unique_ptr<padding_type>
    computePaddingForCrsMatrixUnpack(
      const Kokkos::DualView<const local_ordinal_type*,
        buffer_device_type>& importLIDs,
      Kokkos::DualView<char*, buffer_device_type> imports,
      Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
      const bool verbose) const;

    void
    computeCrsPaddingForSameIDs(
      padding_type& padding,
      const RowGraph<local_ordinal_type, global_ordinal_type,
        node_type>& source,
      const local_ordinal_type numSameIDs) const;

    void
    computeCrsPaddingForPermutedIDs(
      padding_type& padding,
      const RowGraph<local_ordinal_type, global_ordinal_type,
        node_type>& source,
      const Kokkos::DualView<const local_ordinal_type*,
        buffer_device_type>& permuteToLIDs,
      const Kokkos::DualView<const local_ordinal_type*,
      buffer_device_type>& permuteFromLIDs) const;

    virtual void
    packAndPrepare(
      const SrcDistObject& source,
      const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs,
      Kokkos::DualView<packet_type*, buffer_device_type>& exports,
      Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
      size_t& constantNumPackets,
      Distributor& distor) override;

    virtual void
    pack (const Teuchos::ArrayView<const local_ordinal_type>& exportLIDs,
          Teuchos::Array<global_ordinal_type>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets,
          Distributor& distor) const override;

    void
    packFillActive (const Teuchos::ArrayView<const local_ordinal_type>& exportLIDs,
                    Teuchos::Array<global_ordinal_type>& exports,
                    const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                    size_t& constantNumPackets,
                    Distributor& distor) const;

    void
    packFillActiveNew (const Kokkos::DualView<const local_ordinal_type*,
                         buffer_device_type>& exportLIDs,
                       Kokkos::DualView<packet_type*,
                         buffer_device_type>& exports,
                       Kokkos::DualView<size_t*,
                         buffer_device_type> numPacketsPerLID,
                       size_t& constantNumPackets,
                       Distributor& distor) const;

    virtual void
    unpackAndCombine
    (const Kokkos::DualView<const local_ordinal_type*,
       buffer_device_type>& importLIDs,
     Kokkos::DualView<packet_type*,
       buffer_device_type> imports,
     Kokkos::DualView<size_t*,
       buffer_device_type> numPacketsPerLID,
     const size_t constantNumPackets,
     Distributor& distor,
     const CombineMode combineMode) override;

    //@}
    //! \name Advanced methods, at increased risk of deprecation.
    //@{

    /// \brief Get offsets of the diagonal entries in the graph.
    ///
    /// \warning This method is only for expert users.
    /// \warning We make no promises about backwards compatibility
    ///   for this method.  It may disappear or change at any time.
    /// \warning This method must be called collectively.  We reserve
    ///   the right to do extra checking in a debug build that will
    ///   require collectives.
    ///
    /// This method helps users optimize
    /// Tpetra::CrsMatrix::getLocalDiagCopy and
    /// Tpetra::Experimental::BlockCrsMatrix::getLocalDiagCopy, for
    /// several calls when the graph's structure does not change.  The
    /// method fills an array of offsets of the local diagonal entries
    /// in the matrix.  getLocalDiagCopy uses the offsets to extract
    /// the diagonal entries directly, without needing to search for
    /// them using Map lookups and search in each row of the graph.
    ///
    /// The output array's contents are not defined in any other
    /// context other than for use in getLocalDiagCopy.  For example,
    /// you should not rely on \c offsets(i) being the index of the
    /// diagonal entry in the views returned by
    /// Tpetra::CrsMatrix::getLocalRowView.  This may be the case, but
    /// it need not be.  (For example, we may choose to optimize the
    /// lookups down to the optimized storage level, in which case the
    /// offsets will be computed with respect to the underlying
    /// storage format, rather than with respect to the views.)
    ///
    /// Changes to the graph's structure, or calling fillComplete on
    /// the graph (if its structure is not already fixed), may make
    /// the output array's contents invalid.  "Invalid" means that you
    /// must call this method again to recompute the offsets.
    ///
    /// \pre The graph must have a column Map.
    /// \pre All diagonal entries of the graph must be populated on
    ///   this process.  Results are undefined otherwise.
    /// \pre <tt>offsets.extent(0) >= this->getNodeNumRows()</tt>
    ///
    /// \param offsets [out] Output array of offsets.  This method
    ///   does NOT allocate the array; the caller must allocate.  Must
    ///   have getNodeNumRows() entries on the calling process.  (This
    ///   may be different on different processes.)
    void
    getLocalDiagOffsets (const Kokkos::View<size_t*, device_type, Kokkos::MemoryUnmanaged>& offsets) const;

    /// \brief Backwards compatibility overload of the above method.
    ///
    /// This method takes a Teuchos::ArrayRCP instead of a
    /// Kokkos::View.  It also reallocates the output array if it is
    /// not long enough.
    ///
    /// \param offsets [out] Output array of offsets.  This method
    ///   reallocates the array if it is not long enough.  This is why
    ///   the method takes the array by reference.
    void
    getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const;

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
    /// \pre columnIndices are sorted within rows
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
    /// \pre columnIndices are sorted within rows
    /// \pre <tt>hasColMap() == true</tt>
    /// \pre <tt>rowPointers.size() != getNodeNumRows()+1</tt>
    /// \pre No insert routines have been called.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    setAllIndices (const Teuchos::ArrayRCP<size_t> & rowPointers,
                   const Teuchos::ArrayRCP<local_ordinal_type> & columnIndices);

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
    Teuchos::ArrayRCP<const local_ordinal_type> getNodePackedIndices () const;

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
    ///   Within Tpetra, there are no particular restrictions on the column map.
    ///   However, if this graph will be used in Xpetra, Ifpack2, or MueLu,
    ///   the column map's list of global indices must follow "Aztec ordering":
    ///   locally owned GIDs (same order as in domain map), followed by remote GIDs
    ///   (in order of owning proc, and sorted within each proc).
    ///
    ///   It is strongly recommended to use \c Tpetra::Details::makeColMap()
    ///   to create the column map. makeColMap() follows Aztec ordering by default.
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

    /// \brief Replace the current domain Map and Import with the
    ///   given parameters.
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

    /// \brief Remove processes owning zero rows from the Maps and
    ///   their communicator.
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
    removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap) override;
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

  private:
    // Friend declaration for nonmember function.
    template<class CrsGraphType>
    friend Teuchos::RCP<CrsGraphType>
    importAndFillCompleteCrsGraph (const Teuchos::RCP<const CrsGraphType>& sourceGraph,
                                    const Import<typename CrsGraphType::local_ordinal_type,
                                                 typename CrsGraphType::global_ordinal_type,
                                                 typename CrsGraphType::node_type>& importer,
                                    const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                                 typename CrsGraphType::global_ordinal_type,
                                                                 typename CrsGraphType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                                 typename CrsGraphType::global_ordinal_type,
                                                                 typename CrsGraphType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

    // Friend declaration for nonmember function.
    template<class CrsGraphType>
    friend Teuchos::RCP<CrsGraphType>
    importAndFillCompleteCrsGraph (const Teuchos::RCP<const CrsGraphType>& sourceGraph,
                                    const Import<typename CrsGraphType::local_ordinal_type,
                                                 typename CrsGraphType::global_ordinal_type,
                                                 typename CrsGraphType::node_type>& rowImporter,
                                   const Import<typename CrsGraphType::local_ordinal_type,
                                                typename CrsGraphType::global_ordinal_type,
                                                typename CrsGraphType::node_type>& domainImporter,
                                    const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                                 typename CrsGraphType::global_ordinal_type,
                                                                 typename CrsGraphType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                                 typename CrsGraphType::global_ordinal_type,
                                                                 typename CrsGraphType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);


    // Friend declaration for nonmember function.
    template<class CrsGraphType>
    friend Teuchos::RCP<CrsGraphType>
    exportAndFillCompleteCrsGraph (const Teuchos::RCP<const CrsGraphType>& sourceGraph,
                                    const Export<typename CrsGraphType::local_ordinal_type,
                                                 typename CrsGraphType::global_ordinal_type,
                                                 typename CrsGraphType::node_type>& exporter,
                                    const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                                 typename CrsGraphType::global_ordinal_type,
                                                                 typename CrsGraphType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                                 typename CrsGraphType::global_ordinal_type,
                                                                 typename CrsGraphType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

    // Friend declaration for nonmember function.
    template<class CrsGraphType>
    friend Teuchos::RCP<CrsGraphType>
    exportAndFillCompleteCrsGraph (const Teuchos::RCP<const CrsGraphType>& sourceGraph,
                                    const Export<typename CrsGraphType::local_ordinal_type,
                                                 typename CrsGraphType::global_ordinal_type,
                                                 typename CrsGraphType::node_type>& rowExporter,
                                    const Export<typename CrsGraphType::local_ordinal_type,
                                                 typename CrsGraphType::global_ordinal_type,
                                                 typename CrsGraphType::node_type>& domainExporter,
                                    const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                                 typename CrsGraphType::global_ordinal_type,
                                                                 typename CrsGraphType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                                 typename CrsGraphType::global_ordinal_type,
                                                                 typename CrsGraphType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

  public:
    /// \brief Import from <tt>this</tt> to the given destination
    ///   graph, and make the result fill complete.
    ///
    /// If destGraph.is_null(), this creates a new graph as the
    /// destination.  (This is why destGraph is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// graph has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsGraph, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    importAndFillComplete (Teuchos::RCP<CrsGraph<local_ordinal_type, global_ordinal_type, Node> >& destGraph,
                           const import_type& importer,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

    /// \brief Import from <tt>this</tt> to the given destination
    ///   graph, and make the result fill complete.
    ///
    /// If destGraph.is_null(), this creates a new graph as the
    /// destination.  (This is why destGraph is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// graph has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsGraph, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    importAndFillComplete (Teuchos::RCP<CrsGraph<local_ordinal_type, global_ordinal_type, Node> >& destGraph,
                           const import_type& rowImporter,
                           const import_type& domainImporter,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params) const;


    /// \brief Export from <tt>this</tt> to the given destination
    ///   graph, and make the result fill complete.
    ///
    /// If destGraph.is_null(), this creates a new graph as the
    /// destination.  (This is why destGraph is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// graph has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsGraph, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    exportAndFillComplete (Teuchos::RCP<CrsGraph<local_ordinal_type, global_ordinal_type, Node> >& destGraph,
                           const export_type& exporter,
                           const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                           const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

    /// \brief Export from <tt>this</tt> to the given destination
    ///   graph, and make the result fill complete.
    ///
    /// If destGraph.is_null(), this creates a new graph as the
    /// destination.  (This is why destGraph is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// graph has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsGraph, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    exportAndFillComplete (Teuchos::RCP<CrsGraph<local_ordinal_type, global_ordinal_type, Node> >& destGraph,
                           const export_type& rowExporter,
                           const export_type& domainExporter,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params) const;


  private:
    /// \brief Transfer (e.g. Import/Export) from <tt>this</tt> to the
    ///   given destination graph, and make the result fill complete.
    ///
    /// If destGraph.is_null(), this creates a new graph, otherwise it
    /// checks for "pristine" status and throws if that is not the
    /// case.  This method implements importAndFillComplete and
    /// exportAndFillComplete, which in turn implemment the nonmember
    /// "constructors" importAndFillCompleteCrsGraph and
    /// exportAndFillCompleteCrsGraph.  It's convenient to put those
    /// nonmember constructors' implementations inside the CrsGraph
    /// class, so that we don't have to put much code in the _decl
    /// header file.
    ///
    /// The point of this method is to fuse three tasks:
    ///
    ///   1. Create a destination graph (CrsGraph constructor)
    ///   2. Import or Export this graph to the destination graph
    ///   3. Call fillComplete on the destination graph
    ///
    /// Fusing these tasks can avoid some communication and work.
    void
    transferAndFillComplete (Teuchos::RCP<CrsGraph<local_ordinal_type, global_ordinal_type, Node> >& destGraph,
                             const ::Tpetra::Details::Transfer<local_ordinal_type, global_ordinal_type, Node>& rowTransfer,
                             const Teuchos::RCP<const ::Tpetra::Details::Transfer<local_ordinal_type, global_ordinal_type, Node> > & domainTransfer,
                             const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                             const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
                             const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

  protected:
    // these structs are conveniences, to cut down on the number of
    // arguments to some of the methods below.
    struct SLocalGlobalViews {
      Teuchos::ArrayView<const global_ordinal_type> ginds;
      Teuchos::ArrayView<const local_ordinal_type>  linds;
    };
    struct SLocalGlobalNCViews {
      Teuchos::ArrayView<global_ordinal_type>       ginds;
      Teuchos::ArrayView<local_ordinal_type>        linds;
    };

    bool indicesAreAllocated () const;

    void
    allocateIndices(const ELocalGlobal lg, const bool verbose=false);

    //! \name Methods governing changes between global and local indices
    //@{

    /// \brief Make and set the graph's column Map.
    ///
    /// This method makes the column Map, even if the graph already
    /// has one.  It is the caller's responsibility not to call this
    /// method unnecessarily.
    ///
    /// \param remotePIDs [out] The process ranks corresponding to the
    ///   column Map's "remote" (not on the calling process in the
    ///   domain Map) indices.
    void makeColMap (Teuchos::Array<int>& remotePIDs);

    /// \brief Convert column indices from global to local.
    ///
    /// \pre The graph has a column Map.
    /// \post The graph is locally indexed.
    ///
    /// \param verbose [in] Whether to print verbose debugging output.
    ///   This exists because CrsMatrix may want to control output
    ///   independently of the CrsGraph that it owns.
    ///
    /// \return Error code and error string.  See below.
    ///
    /// First return value is the number of column indices on this
    /// process, counting duplicates, that could not be converted to
    /// local indices, because they were not in the column Map on the
    /// calling process.  If some error occurred before conversion
    /// happened, then this is
    /// <tt>Tpetra::Details::OrdinalTraits<size_t>::invalid()</tt>.
    ///
    /// Second return value is a human-readable error string.  If the
    /// first return value is zero, then the string may be empty.
    std::pair<size_t, std::string>
    makeIndicesLocal(const bool verbose=false);

    /// \brief Make the Import and Export objects, if needed.
    ///
    /// \param remotePIDs [in/out] On input: the output of
    ///   makeColMap().  May be modified on output.
    ///
    /// \param useRemotePIDs [in] Whether to use remotePIDs.  Use it
    ///   if we called makeColMap with this as the output argument,
    ///   else don't use it.
    void
    makeImportExport (Teuchos::Array<int>& remotePIDs,
                      const bool useRemotePIDs);

    //@}
    //! \name Methods for inserting indices or transforming values
    //@{

    /// \brief Insert indices into the given row.
    ///
    /// \pre <tt>! (lg == LocalIndices && I == GlobalIndices)</tt>.
    ///   It does not make sense to give this method local column
    ///   indices (meaning that the graph has a column Map), yet to
    ///   ask it to store global indices.
    ///
    /// This method does no allocation; it just inserts the indices.
    ///
    /// \param rowInfo [in/out] On input: Result of CrsGraph's
    ///   getRowInfo() or updateAllocAndValues() methods, for the
    ///   locally owned row (whose local index is
    ///   <tt>rowInfo.localRow</tt>) for which you want to insert
    ///   indices.  On output: numEntries field is updated.
    ///
    /// \param newInds [in] View of the column indices to insert.  If
    ///   <tt>lg == GlobalIndices</tt>, then newInds.ginds, a
    ///   <tt>Teuchos::ArrayView<const global_ordinal_type></tt>, contains
    ///   the (global) column indices to insert.  Otherwise, if <tt>lg
    ///   == LocalIndices</tt>, then newInds.linds, a
    ///   <tt>Teuchos::ArrayView<const local_ordinal_type></tt>, contains
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
    ///
    /// \return The number of indices inserted.
    size_t
    insertIndices (RowInfo& rowInfo,
                   const SLocalGlobalViews& newInds,
                   const ELocalGlobal lg,
                   const ELocalGlobal I);

    /// \brief Insert global indices, using an input <i>local</i> row index.
    ///
    /// \param rowInfo [in] Result of getRowInfo() on the row in which
    ///   to insert.
    /// \param inputGblColInds [in] Input global column indices.
    /// \param numInputInds [in] The number of input global column
    ///   indices.
    ///
    /// \return The number of indices inserted.
    size_t
    insertGlobalIndicesImpl (const local_ordinal_type lclRow,
                             const global_ordinal_type inputGblColInds[],
                             const size_t numInputInds);

    /// \brief Insert global indices, using an input RowInfo.
    ///
    /// \param rowInfo [in] Result of getRowInfo() on the row in which
    ///   to insert.
    /// \param inputGblColInds [in] Input global column indices.
    /// \param numInputInds [in] The number of input global column
    ///   indices.
    ///
    /// \return The number of indices inserted.
    size_t
    insertGlobalIndicesImpl (const RowInfo& rowInfo,
                             const global_ordinal_type inputGblColInds[],
                             const size_t numInputInds,
                             std::function<void(const size_t, const size_t, const size_t)> fun =
                                 std::function<void(const size_t, const size_t, const size_t)>());

    void
    insertLocalIndicesImpl (const local_ordinal_type lclRow,
                            const Teuchos::ArrayView<const local_ordinal_type>& gblColInds,
                            std::function<void(const size_t, const size_t, const size_t)> fun =
                                std::function<void(const size_t, const size_t, const size_t)>());

    /// \brief Finds indices in the given row.
    ///
    /// This method does no insertion; it just finds indices and calls
    /// a callback for each found index
    ///
    /// \param row [in] Row of interest
    ///
    /// \param indices [in] Column indices to find in row
    ///
    /// \param fun Call back function called at each found index.  Called as
    ///   fun(k, start, offset); where k is the index in to indices, start is
    ///   offset to the start of the row, and offset is the relative offset of
    ///   indices[k] in the graphs indices.
    ///
    /// \return The number of indices found.
    size_t
    findLocalIndices(const RowInfo& rowInfo,
                     const Teuchos::ArrayView<const local_ordinal_type>& indices,
                     std::function<void(const size_t, const size_t, const size_t)> fun) const;

    size_t
    findGlobalIndices(const RowInfo& rowInfo,
                      const Teuchos::ArrayView<const global_ordinal_type>& indices,
                      std::function<void(const size_t, const size_t, const size_t)> fun) const;

    /// \brief Like insertGlobalIndices(), but with column Map filtering.
    ///
    /// "Column Map filtering" means that any column indices not in
    /// the column Map on the calling process, get silently dropped.
    ///
    /// \param lclRow [in] Local index of the row in which to insert.
    ///   This row MUST be in the row Map on the calling process.
    /// \param gblColInds [in] The global column indices to insert
    ///   into that row.
    /// \param numGblColInds [in] The number of global column indices
    ///   to insert into that row.
    void
    insertGlobalIndicesFiltered (const local_ordinal_type lclRow,
                                 const global_ordinal_type gblColInds[],
                                 const local_ordinal_type numGblColInds);

    /// \brief Implementation of insertGlobalIndices for nonowned rows.
    ///
    /// A global row index is <i>nonowned</i> when it is not in the
    /// column Map on the calling process.
    ///
    /// \param gblRow [in] Global index of the row in which to insert.
    ///   This row must NOT be in the row Map on the calling process.
    /// \param gblColInds [in] The global column indices to insert
    ///   into that row.
    /// \param numGblColInds [in] The number of global column indices
    ///   to insert into that row.
    void
    insertGlobalIndicesIntoNonownedRows (const global_ordinal_type gblRow,
                                         const global_ordinal_type gblColInds[],
                                         const local_ordinal_type numGblColInds);

    /// \brief Whether transformLocalValues should use atomic updates
    ///   by default.
    ///
    /// \warning This is an implementation detail.
    static const bool useAtomicUpdatesByDefault =
#ifdef KOKKOS_ENABLE_SERIAL
      ! std::is_same<execution_space, Kokkos::Serial>::value;
#else
      true;
#endif // KOKKOS_ENABLE_SERIAL

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

  private:
    /// \brief Sort and merge the column indices in all the rows.
    ///
    /// \param sorted [in] Whether the indices are already sorted.
    /// \param merged [in] Whether the indices are already merged.
    void
    sortAndMergeAllIndices (const bool sorted, const bool merged);

    // mfh 08 May 2017: I only restore "protected" here for backwards
    // compatibility.
  protected:
    /// \brief Sort and merge duplicate column indices in the given row.
    ///
    /// \pre The graph is locally indexed:
    ///   <tt>isGloballyIndexed() == false</tt>.
    /// \pre The graph is not already storage optimized:
    ///   <tt>isStorageOptimized() == false</tt>
    ///
    /// \return The number of duplicate column indices eliminated from the row.
    size_t sortAndMergeRowIndices (const RowInfo& rowInfo,
                                   const bool sorted,
                                   const bool merged);
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

  public:
    //! Returns true if globalConstants have been computed; false otherwise
    bool haveGlobalConstants() const { return haveGlobalConstants_;}

    /// \brief Compute global constants, if they have not yet been computed.
    ///
    /// \warning This is an implementation detail of Tpetra.  It may
    ///   change or disappear at any time.  It is public only because
    ///   MueLu setup needs it to be public.
    ///
    /// Global constants include:
    /// <ul>
    /// <li> globalNumEntries_ </li>
    /// <li> globalMaxNumRowEntries_ </li>
    /// </ul>
    ///
    /// Always compute the following:
    /// <ul>
    /// <li> globalNumEntries_ </li>
    /// <li> globalMaxNumRowEntries_ </li>
    /// </ul>
    void computeGlobalConstants ();

  protected:
    /// \brief Compute local constants, if they have not yet been computed.
    ///
    /// \warning You MUST call fillLocalGraph (or
    ///   CrsMatrix::fillLocalGraphAndMatrix) before calling this
    ///   method!  This method depends on the Kokkos::StaticCrsGraph
    ///   (local_graph_type) object being ready.
    ///
    /// Local constants include:
    /// <ul>
    /// <li> nodeMaxNumRowEntries_ </li>
    /// </ul>
    ///
    /// Always compute the following:
    /// <ul>
    /// <li> nodeMaxNumRowEntries_ </li>
    /// </ul>
    ///
    /// computeGlobalConstants calls this method, if global constants
    /// have not yet been computed.
    void computeLocalConstants ();

    /// \brief Get information about the locally owned row with local
    ///   index myRow.
    RowInfo getRowInfo (const local_ordinal_type myRow) const;

    /// \brief Get information about the locally owned row with global
    ///   index gblRow.
    ///
    /// If \c gblRow is not locally owned, the \c localRow field of
    /// the returned struct is
    /// <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>.
    ///
    /// The point of this method is to fuse the global-to-local row
    /// index lookup for checking whether \c gblRow is locally owned,
    /// with other error checking that getRowInfo() does.  This avoids
    /// an extra global-to-local index lookup in methods like
    /// CrsMatrix::replaceGlobalValues().
    RowInfo getRowInfoFromGlobalRowIndex (const global_ordinal_type gblRow) const;

    /// \brief Get a const, nonowned, locally indexed view of the
    ///   locally owned row myRow, such that rowinfo =
    ///   getRowInfo(myRow).
    Teuchos::ArrayView<const local_ordinal_type>
    getLocalView (const RowInfo& rowinfo) const;

    /// \brief Get a nonconst, nonowned, locally indexed view of the
    ///   locally owned row myRow, such that rowinfo =
    ///   getRowInfo(myRow).
    Teuchos::ArrayView<local_ordinal_type>
    getLocalViewNonConst (const RowInfo& rowinfo);

    /// \brief Get a pointer to the local column indices of a locally
    ///   owned row, using the result of getRowInfo.
    ///
    /// \param lclInds [out] Pointer to the local column indices of
    ///   the given row.
    /// \param capacity [out] Capacity of (number of entries that can
    ///   fit in) the given row.
    /// \param rowInfo [in] Result of getRowInfo(lclRow) for the row
    ///   \c lclRow to view.
    ///
    /// \return 0 if successful, else a nonzero error code.
    local_ordinal_type
    getLocalViewRawConst (const local_ordinal_type*& lclInds,
                          local_ordinal_type& capacity,
                          const RowInfo& rowInfo) const;

  private:

    /// \brief Get a const nonowned view of the local column indices
    ///   indices of row rowinfo.localRow (only works if the matrix is
    ///   locally indexed on the calling process).
    ///
    /// \param rowInfo [in] Result of calling getRowInfo with the
    ///   index of the local row to view.
    Kokkos::View<const local_ordinal_type*, execution_space, Kokkos::MemoryUnmanaged>
    getLocalKokkosRowView (const RowInfo& rowInfo) const;

    /// \brief Get a nonconst nonowned view of the local column
    ///   indices of row rowinfo.localRow (only works if the matrix is
    ///   locally indexed on the calling process).
    ///
    /// \param rowInfo [in] Result of calling getRowInfo with the
    ///   index of the local row to view.
    Kokkos::View<local_ordinal_type*, execution_space, Kokkos::MemoryUnmanaged>
    getLocalKokkosRowViewNonConst (const RowInfo& rowInfo);

    /// \brief Get a const nonowned view of the global column indices
    ///   of row rowinfo.localRow (only works if the matrix is
    ///   globally indexed).
    ///
    /// \param rowInfo [in] Result of calling getRowInfo with the
    ///   index of the local row to view.
    Kokkos::View<const global_ordinal_type*, execution_space, Kokkos::MemoryUnmanaged>
    getGlobalKokkosRowView (const RowInfo& rowInfo) const;

  protected:

    /// \brief Get a const, nonowned, globally indexed view of the
    ///   locally owned row myRow, such that rowinfo =
    ///   getRowInfo(myRow).
    Teuchos::ArrayView<const global_ordinal_type>
    getGlobalView (const RowInfo& rowinfo) const;

    /// \brief Get a nonconst, nonowned, globally indexed view of the
    ///   locally owned row myRow, such that rowinfo =
    ///   getRowInfo(myRow).
    Teuchos::ArrayView<global_ordinal_type>
    getGlobalViewNonConst (const RowInfo& rowinfo);

    /// \brief Get a pointer to the global column indices of a locally
    ///   owned row, using the result of getRowInfoFromGlobalRowIndex.
    ///
    /// \param gblInds [out] Pointer to the global column indices of
    ///   the given row.
    /// \param capacity [out] Capacity of (number of entries that can
    ///   fit in) the given row.
    /// \param rowInfo [in] Result of
    ///   getRowInfoFromGlobalRowIndex(gblRow) for the row to view,
    ///   whose global row index is \c gblRow.
    ///
    /// \return 0 if successful, else a nonzero error code.
    local_ordinal_type
    getGlobalViewRawConst (const global_ordinal_type*& gblInds,
                           local_ordinal_type& capacity,
                           const RowInfo& rowInfo) const;


  public:


    /// \brief Get the local graph.
    ///
    /// \warning THIS IS AN EXPERT MODE FUNCTION.  THIS IS AN
    ///   IMPLEMENTATION DETAIL.  DO NOT CALL THIS FUNCTION!!!
    ///
    /// This is only a valid representation of the local graph if the
    /// (global) graph is fill complete.
    local_graph_type getLocalGraph () const;


  protected:


    void fillLocalGraph (const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Throw an exception if the internal state is not consistent.
    void checkInternalState () const;

    /// \brief Swaps the data from *this with the data and maps from graph.
    ///
    /// \param graph [in/out] a crsGraph
    void swap(CrsGraph<local_ordinal_type, global_ordinal_type, Node> & graph);

    // Friend the tester for CrsGraph::swap
    friend class Tpetra::crsGraph_Swap_Tester<local_ordinal_type, global_ordinal_type, Node>;


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

    /// \brief Local maximum of the number of entries in each row.
    ///
    /// Computed in computeLocalConstants; only valid when
    ///   isFillComplete() is true.
    size_t nodeMaxNumRowEntries_ =
      Teuchos::OrdinalTraits<size_t>::invalid();

    /// \brief Global number of entries in the graph.
    ///
    /// Only valid when isFillComplete() is true.
    global_size_t globalNumEntries_ =
      Teuchos::OrdinalTraits<global_size_t>::invalid();

    /// \brief Global maximum of the number of entries in each row.
    ///
    /// Computed in computeGlobalConstants(); only valid when
    ///   isFillComplete() is true.
    global_size_t globalMaxNumRowEntries_ =
      Teuchos::OrdinalTraits<global_size_t>::invalid();

    /// \brief The maximum number of entries to allow in each locally
    ///   owned row, per row.
    ///
    /// This comes in as an argument to some of the graph's
    /// constructors.  Either this or numAllocForAllRows_ is used, but
    /// not both.  allocateIndices(), setAllIndices(), and
    /// expertStaticFillComplete() all deallocate this array once they
    /// are done with it.
    ///
    /// This is a host View because it is only ever used on the host.
    /// It has the HostMirror type for backwards compatibility; this
    /// used to be a DualView with default layout, so making this a
    /// HostMirror ensures that we can still take it directly by
    /// assignment from the constructors that take DualView, without a
    /// deep copy.
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
    typename Kokkos::View<const size_t*, execution_space>::HostMirror
    k_numAllocPerRow_;

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

    //! \name Graph data structures (packed and unpacked storage).
    //@{

    /// \brief Local column indices for all rows.
    ///
    /// This is only allocated if
    ///
    ///   - The calling process has a nonzero number of entries
    ///   - The graph is locally indexed
    typename local_graph_type::entries_type::non_const_type k_lclInds1D_;

    //! Type of the k_gblInds1D_ array of global column indices.
    typedef Kokkos::View<global_ordinal_type*, execution_space> t_GlobalOrdinal_1D;

    /// \brief Global column indices for all rows.
    ///
    /// This is only allocated if
    ///
    ///   - The calling process has a nonzero number of entries
    ///   - The graph is globally indexed
    t_GlobalOrdinal_1D k_gblInds1D_;

    /// \brief Row offsets for "1-D" storage.
    ///
    /// This is only allocated if "1-D" storage is active.  In that
    /// case, if beg = k_rowPtrs_(i_lcl) and end = k_rowPtrs_(i_lcl+1)
    /// for local row index i_lcl, then
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
    /// If it is allocated, k_rowPtrs_ has length getNodeNumRows()+1.
    /// The k_numRowEntries_ array has has length getNodeNumRows(),
    /// again if it is allocated.
    typename local_graph_type::row_map_type::const_type k_rowPtrs_;

    /// \brief The type of k_numRowEntries_ (see below).
    ///
    /// This View gets used only on host.  However, making this
    /// literally a host View (of Kokkos::HostSpace) causes
    /// inexplicable test failures only on CUDA.  Thus, I left it as a
    /// HostMirror, which means (given Trilinos' current UVM
    /// requirement) that it will be a UVM allocation.
    typedef typename Kokkos::View<size_t*, Kokkos::LayoutLeft, device_type>::HostMirror num_row_entries_type;

    // typedef Kokkos::View<
    //   size_t*,
    //   Kokkos::LayoutLeft,
    //   Kokkos::Device<
    //     typename Kokkos::View<
    //       size_t*,
    //       Kokkos::LayoutLeft,
    //       device_type>::HostMirror::execution_space,
    //     Kokkos::HostSpace> > num_row_entries_type;

    /// \brief The number of local entries in each locally owned row.
    ///
    /// This is deallocated in fillComplete() if fillComplete()'s
    /// "Optimize Storage" parameter is set to \c true.
    ///
    /// This may also exist with 1-D storage, if storage is unpacked.
    num_row_entries_type k_numRowEntries_;

    //@}

    /// \brief Status of the graph's storage, when not in a
    ///   fill-complete state.
    ///
    /// The phrase "When not in a fill-complete state" is important.
    /// When the graph is fill complete, it <i>always</i> uses 1-D
    /// "packed" storage.  However, if the "Optimize Storage"
    /// parameter to fillComplete was false, the graph may keep
    /// unpacked 1-D storage around and resume it on the next
    /// resumeFill call.
    Details::EStorageStatus storageStatus_ =
      Details::STORAGE_1D_UNPACKED;

    bool indicesAreAllocated_ = false;
    bool indicesAreLocal_ = false;
    bool indicesAreGlobal_ = false;
    bool fillComplete_ = false;

    //! Whether the graph's indices are sorted in each row, on this process.
    bool indicesAreSorted_ = true;
    /// \brief Whether the graph's indices are non-redundant (merged)
    ///   in each row, on this process.
    bool noRedundancies_ = true;
    //! Whether this process has computed local constants.
    bool haveLocalConstants_ = false;
    //! Whether all processes have computed global constants.
    bool haveGlobalConstants_ = false;

    typedef typename std::map<global_ordinal_type, std::vector<global_ordinal_type> > nonlocals_type;

    //! Nonlocal data given to insertGlobalIndices.
    nonlocals_type nonlocals_;

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
    bool sortGhostsAssociatedWithEachProcessor_ = true;

  private:
    //! Get initial value of debug_ for this object.
    static bool getDebug();

    /// \brief Whether to do extra debug checks.
    /// This comes from Tpetra::Details::Behavior::debug("CrsGraph").
    bool debug_ = getDebug();

    //! Get initial value of verbose_ for this object.
    static bool getVerbose();

    /// \brief Whether to do extra debug checks.
    ///
    /// This comes from Tpetra::Details::Behavior::debug("CrsGraph").
    bool verbose_ = getVerbose();

  }; // class CrsGraph

  /// \brief Nonmember function to create an empty CrsGraph given a
  ///   row Map and the max number of entries allowed locally per row.
  ///
  /// \return A graph with at most the specified number of nonzeros
  ///   per row (defaults to zero).
  ///
  /// \relatesalso CrsGraph
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
  createCrsGraph(
    const Teuchos::RCP<
      const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
    size_t maxNumEntriesPerRow = 0,
    const Teuchos::RCP<Teuchos::ParameterList>& params =
      Teuchos::null)
  {
    using Teuchos::rcp;
    using graph_type = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
    const ProfileType pftype = TPETRA_DEFAULT_PROFILE_TYPE;
    return rcp(new graph_type(map, maxNumEntriesPerRow,
                              pftype, params));
  }

  /// \brief Nonmember CrsGraph constructor that fuses Import and
  ///   fillComplete().
  /// \relatesalso CrsGraph
  /// \tparam CrsGraphType A specialization of CrsGraph.
  ///
  /// A common use case is to create an empty destination CrsGraph,
  /// redistribute from a source CrsGraph (by an Import or Export
  /// operation), then call fillComplete() on the destination
  /// CrsGraph.  This constructor fuses these three cases, for an
  /// Import redistribution.
  ///
  /// Fusing redistribution and fillComplete() exposes potential
  /// optimizations.  For example, it may make constructing the column
  /// Map faster, and it may avoid intermediate unoptimized storage in
  /// the destination CrsGraph.
  ///
  /// The resulting graph is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source graph, and its range Map is the range Map of
  /// the source graph.
  ///
  /// \warning If the target Map of the Import is a subset of the
  ///   source Map of the Import, then you cannot use the default
  ///   range Map.  You should instead construct a nonoverlapping
  ///   version of the target Map and supply that as the nondefault
  ///   value of the range Map.
  ///
  /// \param sourceGraph [in] The source graph from which to
  ///   import.  The source of an Import must have a nonoverlapping
  ///   distribution.
  ///
  /// \param importer [in] The Import instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Import must be the same as the rowMap of sourceGraph unless
  ///   the "Reverse Mode" option on the params list, in which case
  ///   the targetMap of Import must match the rowMap of the sourceGraph
  ///
  /// \param domainMap [in] Domain Map of the returned graph.  If
  ///   null, we use the default, which is the domain Map of the
  ///   source graph.
  ///
  /// \param rangeMap [in] Range Map of the returned graph.  If
  ///   null, we use the default, which is the range Map of the
  ///   source graph.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsGraphType>
  Teuchos::RCP<CrsGraphType>
  importAndFillCompleteCrsGraph (const Teuchos::RCP<const CrsGraphType>& sourceGraph,
                                  const Import<typename CrsGraphType::local_ordinal_type,
                                               typename CrsGraphType::global_ordinal_type,
                                               typename CrsGraphType::node_type>& importer,
                                  const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                               typename CrsGraphType::global_ordinal_type,
                                                               typename CrsGraphType::node_type> >& domainMap = Teuchos::null,
                                  const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                               typename CrsGraphType::global_ordinal_type,
                                                               typename CrsGraphType::node_type> >& rangeMap = Teuchos::null,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
  {
    Teuchos::RCP<CrsGraphType> destGraph;
    sourceGraph->importAndFillComplete (destGraph,importer,domainMap, rangeMap, params);
    return destGraph;
  }

  /// \brief Nonmember CrsGraph constructor that fuses Import and fillComplete().
  /// \relatesalso CrsGraph
  /// \tparam CrsGraphType A specialization of CrsGraph.
  ///
  /// A common use case is to create an empty destination CrsGraph,
  /// redistribute from a source CrsGraph (by an Import or Export
  /// operation), then call fillComplete() on the destination
  /// CrsGraph.  This constructor fuses these three cases, for an
  /// Import redistribution.
  ///
  /// Fusing redistribution and fillComplete() exposes potential
  /// optimizations.  For example, it may make constructing the column
  /// Map faster, and it may avoid intermediate unoptimized storage in
  /// the destination CrsGraph.
  ///
  /// The resulting graph is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source graph, and its range Map is the range Map of
  /// the source graph.
  ///
  /// \warning If the target Map of the Import is a subset of the
  ///   source Map of the Import, then you cannot use the default
  ///   range Map.  You should instead construct a nonoverlapping
  ///   version of the target Map and supply that as the nondefault
  ///   value of the range Map.
  ///
  /// \param sourceGraph [in] The source graph from which to
  ///   import.  The source of an Import must have a nonoverlapping
  ///   distribution.
  ///
  /// \param rowImporter [in] The Import instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Import must be the same as the rowMap of sourceGraph unless
  ///   the "Reverse Mode" option on the params list, in which case
  ///   the targetMap of Import must match the rowMap of the sourceGraph
  ///
  /// \param domainImporter [in] The Import instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Import must be the same as the domainMap of sourceGraph unless
  ///   the "Reverse Mode" option on the params list, in which case
  ///   the targetMap of Import must match the domainMap of the sourceGraph
  ///
  /// \param domainMap [in] Domain Map of the returned graph.
  ///
  /// \param rangeMap [in] Range Map of the returned graph.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsGraphType>
  Teuchos::RCP<CrsGraphType>
  importAndFillCompleteCrsGraph (const Teuchos::RCP<const CrsGraphType>& sourceGraph,
                                  const Import<typename CrsGraphType::local_ordinal_type,
                                               typename CrsGraphType::global_ordinal_type,
                                               typename CrsGraphType::node_type>& rowImporter,
                                  const Import<typename CrsGraphType::local_ordinal_type,
                                              typename CrsGraphType::global_ordinal_type,
                                              typename CrsGraphType::node_type>& domainImporter,
                                  const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                               typename CrsGraphType::global_ordinal_type,
                                                               typename CrsGraphType::node_type> >& domainMap,
                                  const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                               typename CrsGraphType::global_ordinal_type,
                                                               typename CrsGraphType::node_type> >& rangeMap,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    Teuchos::RCP<CrsGraphType> destGraph;
    sourceGraph->importAndFillComplete (destGraph,rowImporter,domainImporter, domainMap, rangeMap, params);
    return destGraph;
  }

  /// \brief Nonmember CrsGraph constructor that fuses Export and fillComplete().
  /// \relatesalso CrsGraph
  /// \tparam CrsGraphType A specialization of CrsGraph.
  ///
  /// For justification, see the documentation of
  /// importAndFillCompleteCrsGraph() (which is the Import analog of
  /// this function).
  ///
  /// The resulting graph is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source graph, and its range Map is the range Map of
  /// the source graph.
  ///
  /// \param sourceGraph [in] The source graph from which to
  ///   export.  Its row Map may be overlapping, since the source of
  ///   an Export may be overlapping.
  ///
  /// \param exporter [in] The Export instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Export must be the same as the row Map of sourceGraph.
  ///
  /// \param domainMap [in] Domain Map of the returned graph.  If
  ///   null, we use the default, which is the domain Map of the
  ///   source graph.
  ///
  /// \param rangeMap [in] Range Map of the returned graph.  If
  ///   null, we use the default, which is the range Map of the
  ///   source graph.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsGraphType>
  Teuchos::RCP<CrsGraphType>
  exportAndFillCompleteCrsGraph (const Teuchos::RCP<const CrsGraphType>& sourceGraph,
                                  const Export<typename CrsGraphType::local_ordinal_type,
                                               typename CrsGraphType::global_ordinal_type,
                                               typename CrsGraphType::node_type>& exporter,
                                  const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                               typename CrsGraphType::global_ordinal_type,
                                                               typename CrsGraphType::node_type> >& domainMap = Teuchos::null,
                                  const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                               typename CrsGraphType::global_ordinal_type,
                                                               typename CrsGraphType::node_type> >& rangeMap = Teuchos::null,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
  {
    Teuchos::RCP<CrsGraphType> destGraph;
    sourceGraph->exportAndFillComplete (destGraph,exporter,domainMap, rangeMap, params);
    return destGraph;
  }

  /// \brief Nonmember CrsGraph constructor that fuses Export and fillComplete().
  /// \relatesalso CrsGraph
  /// \tparam CrsGraphType A specialization of CrsGraph.
  ///
  /// For justification, see the documentation of
  /// importAndFillCompleteCrsGraph() (which is the Import analog of
  /// this function).
  ///
  /// The resulting graph is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source graph, and its range Map is the range Map of
  /// the source graph.
  ///
  /// \param sourceGraph [in] The source graph from which to
  ///   export.  Its row Map may be overlapping, since the source of
  ///   an Export may be overlapping.
  ///
  /// \param rowExporter [in] The Export instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Export must be the same as the row Map of sourceGraph.
  ///
  /// \param domainExporter [in] The Export instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Export must be the same as the domain Map of sourceGraph.
  ///
  /// \param domainMap [in] Domain Map of the returned graph.
  ///
  /// \param rangeMap [in] Range Map of the returned graph.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsGraphType>
  Teuchos::RCP<CrsGraphType>
  exportAndFillCompleteCrsGraph (const Teuchos::RCP<const CrsGraphType>& sourceGraph,
                                  const Export<typename CrsGraphType::local_ordinal_type,
                                               typename CrsGraphType::global_ordinal_type,
                                               typename CrsGraphType::node_type>& rowExporter,
                                  const Export<typename CrsGraphType::local_ordinal_type,
                                               typename CrsGraphType::global_ordinal_type,
                                               typename CrsGraphType::node_type>& domainExporter,
                                  const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                               typename CrsGraphType::global_ordinal_type,
                                                               typename CrsGraphType::node_type> >& domainMap,
                                  const Teuchos::RCP<const Map<typename CrsGraphType::local_ordinal_type,
                                                               typename CrsGraphType::global_ordinal_type,
                                                               typename CrsGraphType::node_type> >& rangeMap,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    Teuchos::RCP<CrsGraphType> destGraph;
    sourceGraph->exportAndFillComplete (destGraph,rowExporter,domainExporter,domainMap, rangeMap, params);
    return destGraph;
  }


} // namespace Tpetra

#endif // TPETRA_CRSGRAPH_DECL_HPP
