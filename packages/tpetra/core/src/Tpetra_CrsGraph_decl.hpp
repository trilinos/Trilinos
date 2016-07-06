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

#ifndef TPETRA_CRSGRAPH_DECL_HPP
#define TPETRA_CRSGRAPH_DECL_HPP

/// \file Tpetra_CrsGraph_decl.hpp
/// \brief Declaration of the Tpetra::CrsGraph class
///
/// If you want to use Tpetra::CrsGraph, include "Tpetra_CrsGraph.hpp"
/// (a file which CMake generates and installs for you).  If you only
/// want the declaration of Tpetra::CrsGraph, include this file
/// (Tpetra_CrsGraph_decl.hpp).

#include "Tpetra_RowGraph.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Exceptions.hpp"
#include "Tpetra_Util.hpp" // need this here for sort2

#include "Kokkos_Sparse_findRelOffset.hpp"
#include "Kokkos_DualView.hpp"
#include "Kokkos_StaticCrsGraph.hpp"

#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  //
  // Dear users: These are just forward declarations.  Please skip
  // over them and go down to the CrsMatrix class declaration.  Thank
  // you.
  //
  template <class LO, class GO, class N, const bool isClassic>
  class CrsGraph;

  // forward declaration (needed for "friend" inside CrsGraph)
  template <class S, class LO, class GO, class N, const bool isClassic>
  class CrsMatrix;

  namespace Experimental {
    // forward declaration (needed for "friend" inside CrsGraph)
    template<class S, class LO, class GO, class N>
    class BlockCrsMatrix;
  } // namespace Experimental

  namespace Details {
    // Forward declaration of an implementation detail of CrsGraph::clone.
    template<class OutputCrsGraphType, class InputCrsGraphType>
    class CrsGraphCopier {
    public:
      static Teuchos::RCP<OutputCrsGraphType>
      clone (const InputCrsGraphType& graphIn,
             const Teuchos::RCP<typename OutputCrsGraphType::node_type> nodeOut,
             const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
    };
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

    template<class T, class BinaryFunction>
    T atomic_binary_function_update (volatile T* const dest, const T& inputVal, BinaryFunction f)
    {
      T oldVal = *dest;
      T assume;

      // NOTE (mfh 30 Nov 2015) I do NOT need a fence here for IBM
      // POWER architectures, because 'newval' depends on 'assume',
      // which depends on 'oldVal', which depends on '*dest'.  This
      // sets up a chain of read dependencies that should ensure
      // correct behavior given a sane memory model.
      do {
        assume = oldVal;
        T newVal = f (assume, inputVal);
        oldVal = Kokkos::atomic_compare_exchange (dest, assume, newVal);
      } while (assume != oldVal);

      return oldVal;
    }

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
  template <class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
            class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
            class Node = ::Tpetra::Details::DefaultTypes::node_type,
            const bool classic = Node::classic>
  class CrsGraph :
    public RowGraph<LocalOrdinal, GlobalOrdinal, Node>,
    public DistObject<GlobalOrdinal,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node>,
    public Teuchos::ParameterListAcceptorDefaultBase
  {
    static_assert (! classic, "The 'classic' version of Tpetra was deprecated long ago, and has been removed.");

    template <class S, class LO, class GO, class N, const bool isClassic>
    friend class CrsMatrix;
    template <class LO2, class GO2, class N2, const bool isClassic>
    friend class CrsGraph;
    template <class S, class LO, class GO, class N>
    friend class ::Tpetra::Experimental::BlockCrsMatrix;

    //! The specialization of DistObject that is this class' parent class.
    typedef DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> dist_object_type;

  public:
    //! This class' first template parameter; the type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! This class' second template parameter; the type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! This class' Kokkos Node type.
    typedef Node node_type;

    //! This class' Kokkos device type.
    typedef typename Node::device_type device_type;
    //! This class' Kokkos execution space.
    typedef typename device_type::execution_space execution_space;

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
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
    //! The Import specialization used by this class.
    typedef Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
    //! The Export specialization used by this class.
    typedef Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node> export_type;

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
    Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node2, Node2::classic> >
    clone (const Teuchos::RCP<Node2>& node2,
           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const
    {
      typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node2, Node2::classic> output_crs_graph_type;
      typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic> input_crs_graph_type;
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
    insertGlobalIndices (const GlobalOrdinal globalRow,
                         const Teuchos::ArrayView<const GlobalOrdinal>& indices);

    /// \brief Epetra compatibility version of insertGlobalIndices
    ///   (see above) that takes input as a raw pointer, rather than
    ///   Teuchos::ArrayView.
    ///
    /// Arguments are the same and in the same order as
    /// Epetra_CrsGraph::InsertGlobalIndices.
    void
    insertGlobalIndices (const GlobalOrdinal globalRow,
                         const LocalOrdinal numEnt,
                         const GlobalOrdinal inds[]);

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

    /// \brief Epetra compatibility version of insertLocalIndices
    ///   (see above) that takes input as a raw pointer, rather than
    ///   Teuchos::ArrayView.
    ///
    /// Arguments are the same and in the same order as
    /// Epetra_CrsGraph::InsertMyIndices.
    void
    insertLocalIndices (const LocalOrdinal localRow,
                        const LocalOrdinal numEnt,
                        const LocalOrdinal inds[]);

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

    /// \brief Get the number of entries in the given row (local index).
    ///
    /// \return The number of entries in the given row, specified by
    ///   local index, on the calling MPI process.  If the specified
    ///   local row index is invalid on the calling process, return
    ///   <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>.
    size_t getNumEntriesInLocalRow (LocalOrdinal localRow) const;

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

    /// \brief Get a const, non-persisting view of the given global
    ///   row's global column indices, as a Teuchos::ArrayView.
    ///
    /// \param GlobalRow [in] Global index of the row.
    /// \param Indices [out] Global column indices in the row.  If the
    ///   given row is not a valid row index on the calling process,
    ///   then the result has no entries (its size is zero).
    ///
    /// \pre <tt>! isLocallyIndexed()</tt>
    /// \post <tt>Indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>
    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      Teuchos::ArrayView<const GlobalOrdinal>& Indices) const;

    /// \brief Get a const, non-persisting view of the given local
    ///   row's local column indices, as a Teuchos::ArrayView.
    ///
    /// \param LocalRow [in] Local index of the row.
    /// \param indices [out] Local column indices in the row.  If the
    ///   given row is not a valid row index on the calling process,
    ///   then the result has no entries (its size is zero).
    ///
    /// \pre <tt>! isGloballyIndexed()</tt>
    /// \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>
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
    /// \pre <tt>offsets.dimension_0() >= this->getNodeNumRows()</tt>
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
    Teuchos::ArrayRCP<Teuchos::Array<T> > allocateValues2D () const
    {
      using Teuchos::arcp;
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      const char tfecfFuncName[] = "allocateValues2D: ";

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! indicesAreAllocated (), std::runtime_error,
         "Graph indices must be allocated before values.");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (getProfileType () != DynamicProfile, std::runtime_error,
         "Graph indices must be allocated in a dynamic profile.");

      ArrayRCP<Array<T> > values2D;
      values2D = arcp<Array<T> > (getNodeNumRows ());
      if (lclInds2D_ != null) {
        const size_t numRows = lclInds2D_.size ();
        for (size_t r = 0; r < numRows; ++r) {
          values2D[r].resize (lclInds2D_[r].size ());
        }
      }
      else if (gblInds2D_ != null) {
        const size_t numRows = gblInds2D_.size ();
        for (size_t r = 0; r < numRows; ++r) {
          values2D[r].resize (gblInds2D_[r].size ());
        }
      }
      return values2D;
    }

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
      static_assert (lg == GlobalIndices || lg == LocalIndices,
                     "Tpetra::CrsGraph::filterIndices: The template parameter "
                     "lg must be either GlobalIndices or LocalIndicies.");

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
                            const ELocalGlobal I)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char tfecfFuncName[] = "insertIndicesAndValues: ";
#endif // HAVE_TPETRA_DEBUG

#ifdef HAVE_TPETRA_DEBUG
      size_t numNewInds = 0;
      try {
        numNewInds = insertIndices (rowInfo, newInds, lg, I);
      } catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::runtime_error, "insertIndices threw an exception: "
           << e.what ());
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (numNewInds > static_cast<size_t> (oldRowVals.size ()),
         std::runtime_error, "numNewInds (" << numNewInds << ") > "
         "oldRowVals.size() (" << oldRowVals.size () << ".");
#else
      const size_t numNewInds = insertIndices (rowInfo, newInds, lg, I);
#endif // HAVE_TPETRA_DEBUG

      typedef typename Teuchos::ArrayView<Scalar>::size_type size_type;

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (rowInfo.numEntries + numNewInds > static_cast<size_t> (oldRowVals.size ()),
         std::runtime_error, "rowInfo.numEntries (" << rowInfo.numEntries << ")"
         " + numNewInds (" << numNewInds << ") > oldRowVals.size() ("
         << oldRowVals.size () << ").");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_type> (numNewInds) > newRowVals.size (),
         std::runtime_error, "numNewInds (" << numNewInds << ") > "
         "newRowVals.size() (" << newRowVals.size () << ").");
#endif // HAVE_TPETRA_DEBUG

      size_type oldInd = static_cast<size_type> (rowInfo.numEntries);

#ifdef HAVE_TPETRA_DEBUG
      try {
#endif // HAVE_TPETRA_DEBUG
        //NOTE: The code in the else branch fails on GCC 4.9 and newer in the assignement oldRowVals[oldInd] = newRowVals[newInd];
        //We supply a workaround n as well as other code variants which produce or not produce the error
#if defined(__GNUC__) && defined(__GNUC_MINOR__) && defined(__GNUC_PATCHLEVEL__)
#define GCC_VERSION __GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__
#if GCC_VERSION >= 490
#define GCC_WORKAROUND
#endif
#endif
#ifdef GCC_WORKAROUND
        size_type nNI = static_cast<size_type>(numNewInds);
        if (nNI > 0)
          memcpy(&oldRowVals[oldInd], &newRowVals[0], nNI*sizeof(Scalar));
        /*
        //Original Code Fails
        for (size_type newInd = 0; newInd < static_cast<size_type> (numNewInds);
        ++newInd, ++oldInd) {
        oldRowVals[oldInd] = newRowVals[newInd];
        }

        //char cast variant fails
        char* oldRowValPtr = (char*)&oldRowVals[oldInd];
        const char* newRowValPtr = (const char*) &newRowVals[0];

        for(size_type newInd = 0; newInd < (nNI * sizeof(Scalar)); newInd++) {
        oldRowValPtr[newInd] = newRowValPtr[newInd];
        }

        //Raw ptr variant fails
        Scalar* oldRowValPtr = &oldRowVals[oldInd];
        Scalar* newRowValPtr = const_cast<Scalar*>(&newRowVals[0]);

        for(size_type newInd = 0; newInd < nNI; newInd++) {
        oldRowValPtr[newInd] = newRowValPtr[newInd];
        }

        //memcpy works
        for (size_type newInd = 0; newInd < nNI; newInd++) {
        memcpy( &oldRowVals[oldInd+newInd], &newRowVals[newInd], sizeof(Scalar));
        }

        //just one loop index fails
        for (size_type newInd = 0; newInd < nNI; newInd++) {
        oldRowVals[oldInd+newInd] = newRowVals[newInd];
        }

        //inline increment fails
        for (size_type newInd = 0; newInd < numNewInds;) {
        oldRowVals[oldInd++] = newRowVals[newInd++];
        }

        */

#else // GCC Workaround above
        for (size_type newInd = 0; newInd < static_cast<size_type> (numNewInds);
             ++newInd, ++oldInd) {
          oldRowVals[oldInd] = newRowVals[newInd];
        }
#endif // GCC Workaround
#ifdef HAVE_TPETRA_DEBUG
      }
      catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::runtime_error, "for loop for copying values threw an "
           "exception: " << e.what ());
      }
#endif // HAVE_TPETRA_DEBUG
    }

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

    /// \brief Whether sumIntoLocalValues and transformLocalValues
    ///   should use atomic updates by default.
    ///
    /// \warning This is an implementation detail.
    static const bool useAtomicUpdatesByDefault =
#ifdef KOKKOS_HAVE_SERIAL
      ! std::is_same<execution_space, Kokkos::Serial>::value;
#else
      true;
#endif // KOKKOS_HAVE_SERIAL

    /// \brief Transform the given values using local indices.
    ///
    /// \tparam LocalIndicesViewType Kokkos::View specialization that
    ///   is a 1-D array of LocalOrdinal.
    /// \tparam OutputScalarViewType Kokkos::View specialization that is
    ///   a 1-D array of the type of values in the sparse matrix (that
    ///   type is CrsMatrix::impl_scalar_type).
    /// \tparam InputScalarViewType Kokkos::View specialization that
    ///   is a 1-D array of the type of values in the sparse matrix,
    ///   but with a possibly different memory space than how the
    ///   matrix stores its values.
    /// \tparam BinaryFunction The type of the binary function f to
    ///   use for updating the output value(s).  This should be
    ///   convertible to
    ///   std::function<impl_scalar_type (const impl_scalar_type&,
    ///                                   const impl_scalar_type&)>.
    ///
    /// \param rowInfo [in] Information about a given row of the graph.
    /// \param rowVals [in/out] The values to be transformed.  They
    ///   correspond to the row indicated by rowInfo.
    /// \param inds [in] The (local) indices in the row, for which
    ///   to transform the corresponding values in rowVals.
    /// \param newVals [in] Values to use for transforming rowVals.
    ///   It's probably OK for these to alias rowVals.
    /// \param f [in] A binary function used to transform rowVals.
    ///
    /// This method transforms the values using the expression
    /// \code
    /// newVals[k] = f( rowVals[k], newVals[j] );
    /// \endcode
    /// where k is the local index corresponding to <tt>inds[j]</tt>.
    /// It ignores invalid local column indices, but they are counted
    /// in the return value.
    ///
    /// \return The number of valid local column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    template<class OutputScalarViewType,
             class LocalIndicesViewType,
             class InputScalarViewType,
             class BinaryFunction>
    LocalOrdinal
    transformLocalValues (const RowInfo& rowInfo,
                          const typename UnmanagedView<OutputScalarViewType>::type& rowVals,
                          const typename UnmanagedView<LocalIndicesViewType>::type& inds,
                          const typename UnmanagedView<InputScalarViewType>::type& newVals,
                          BinaryFunction f,
                          const bool atomic = useAtomicUpdatesByDefault) const
    {
      // We use static_assert here to check the template parameters,
      // rather than std::enable_if (e.g., on the return value, to
      // enable compilation only if the template parameters match the
      // desired attributes).  This turns obscure link errors into
      // clear compilation errors.  It also makes the return value a
      // lot easier to see.
      static_assert (Kokkos::is_view<OutputScalarViewType>::value,
                     "Template parameter OutputScalarViewType must be "
                     "a Kokkos::View.");
      static_assert (Kokkos::is_view<LocalIndicesViewType>::value,
                     "Template parameter LocalIndicesViewType must be "
                     "a Kokkos::View.");
      static_assert (Kokkos::is_view<InputScalarViewType>::value,
                     "Template parameter InputScalarViewType must be a "
                     "Kokkos::View.");
      static_assert (static_cast<int> (OutputScalarViewType::rank) == 1,
                     "Template parameter OutputScalarViewType must "
                     "have rank 1.");
      static_assert (static_cast<int> (LocalIndicesViewType::rank) == 1,
                     "Template parameter LocalIndicesViewType must "
                     "have rank 1.");
      static_assert (static_cast<int> (InputScalarViewType::rank) == 1,
                     "Template parameter InputScalarViewType must have "
                     "rank 1.");
      static_assert (std::is_same<
                       typename OutputScalarViewType::non_const_value_type,
                       typename InputScalarViewType::non_const_value_type>::value,
                     "Template parameters OutputScalarViewType and "
                     "InputScalarViewType must contain values of the same "
                     "type.");
      static_assert (std::is_same<
                       typename LocalIndicesViewType::non_const_value_type,
                       local_ordinal_type>::value,
                     "Template parameter LocalIndicesViewType must "
                     "contain values of type local_ordinal_type.");

      typedef typename OutputScalarViewType::non_const_value_type ST;
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;

      if (newVals.dimension_0 () != inds.dimension_0 ()) {
        // The sizes of the input arrays must match.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }
      const LO numElts = static_cast<LO> (inds.dimension_0 ());
      const bool sorted = this->isSorted ();

      LO numValid = 0; // number of valid input column indices
      size_t hint = 0; // Guess for the current index k into rowVals

      if (isLocallyIndexed ()) {
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getLocalKokkosRowView (rowInfo);

        for (LO j = 0; j < numElts; ++j) {
          const LO lclColInd = inds(j);
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         lclColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            if (atomic) {
              // NOTE (mfh 30 Nov 2015) The commented-out code is
              // wrong because another thread may have changed
              // rowVals(offset) between those two lines of code.
              //
              //const ST newVal = f (rowVals(offset), newVals(j));
              //Kokkos::atomic_assign (&rowVals(offset), newVal);

              volatile ST* const dest = &rowVals(offset);
              (void) atomic_binary_function_update (dest, newVals(j), f);
            }
            else {
              // use binary function f
              rowVals(offset) = f (rowVals(offset), newVals(j));
            }
            hint = offset + 1;
            ++numValid;
          }
        }
      }
      else if (isGloballyIndexed ()) {
        // NOTE (mfh 26 Nov 2015) Dereferencing an RCP or reading its
        // pointer does NOT change its reference count.  Thus, this
        // code is still thread safe.
        if (colMap_.is_null ()) {
          // NO input column indices are valid in this case.  Either
          // the column Map hasn't been set yet (so local indices
          // don't exist yet), or the calling process owns no graph
          // entries.
          return numValid;
        }
        const map_type& colMap = *colMap_;
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getGlobalKokkosRowView (rowInfo);

        const GO GINV = Teuchos::OrdinalTraits<GO>::invalid ();
        for (LO j = 0; j < numElts; ++j) {
          const GO gblColInd = colMap.getGlobalElement (inds(j));
          if (gblColInd != GINV) {
            const size_t offset =
              KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                           gblColInd, hint, sorted);
            if (offset != rowInfo.numEntries) {
              if (atomic) {
                // NOTE (mfh 30 Nov 2015) The commented-out code is
                // wrong because another thread may have changed
                // rowVals(offset) between those two lines of code.
                //
                //const ST newVal = f (rowVals(offset), newVals(j));
                //Kokkos::atomic_assign (&rowVals(offset), newVal);

                volatile ST* const dest = &rowVals(offset);
                (void) atomic_binary_function_update (dest, newVals(j), f);
              }
              else {
                // use binary function f
                rowVals(offset) = f (rowVals(offset), newVals(j));
              }
              hint = offset + 1;
              numValid++;
            }
          }
        }
      }
      // If the graph is neither locally nor globally indexed on the
      // calling process, that means the calling process has no graph
      // entries.  Thus, none of the input column indices are valid.

      return numValid;
    }

    /// \brief Implementation detail of CrsMatrix::sumIntoLocalValues.
    ///
    /// \tparam LocalIndicesViewType Kokkos::View specialization that
    ///   is a 1-D array of LocalOrdinal.
    /// \tparam OutputScalarViewType Kokkos::View specialization that is
    ///   a 1-D array of the type of values in the sparse matrix (that
    ///   type is CrsMatrix::impl_scalar_type).
    /// \tparam InputScalarViewType Kokkos::View specialization that
    ///   is a 1-D array of the type of values in the sparse matrix,
    ///   but with a possibly different memory space than how the
    ///   matrix stores its values.
    ///
    /// \param rowInfo [in] Result of getRowInfo on the index of the
    ///   local row of the matrix to modify.
    /// \param rowVals [in/out] On input: Values of the row of the
    ///   sparse matrix to modify.  On output: The modified values.
    /// \param inds [in] Local column indices of that row to modify.
    /// \param newVals [in] For each k, increment the value in rowVals
    ///   corresponding to local column index inds[k] by newVals[k].
    /// \param atomic [in] Whether to use atomic updates (+=) when
    ///   incrementing values.
    template<class OutputScalarViewType,
             class LocalIndicesViewType,
             class InputScalarViewType>
    LocalOrdinal
    sumIntoLocalValues (const RowInfo& rowInfo,
                        const typename UnmanagedView<OutputScalarViewType>::type& rowVals,
                        const typename UnmanagedView<LocalIndicesViewType>::type& inds,
                        const typename UnmanagedView<InputScalarViewType>::type& newVals,
                        const bool atomic = useAtomicUpdatesByDefault) const
    {
      // We use static_assert here to check the template parameters,
      // rather than std::enable_if (e.g., on the return value, to
      // enable compilation only if the template parameters match the
      // desired attributes).  This turns obscure link errors into
      // clear compilation errors.  It also makes the return value a
      // lot easier to see.
      static_assert (Kokkos::is_view<OutputScalarViewType>::value,
                     "Template parameter OutputScalarViewType must be "
                     "a Kokkos::View.");
      static_assert (Kokkos::is_view<LocalIndicesViewType>::value,
                     "Template parameter LocalIndicesViewType must be "
                     "a Kokkos::View.");
      static_assert (Kokkos::is_view<InputScalarViewType>::value,
                     "Template parameter InputScalarViewType must be a "
                     "Kokkos::View.");
      static_assert (static_cast<int> (OutputScalarViewType::rank) == 1,
                     "Template parameter OutputScalarViewType must "
                     "have rank 1.");
      static_assert (static_cast<int> (LocalIndicesViewType::rank) == 1,
                     "Template parameter LocalIndicesViewType must "
                     "have rank 1.");
      static_assert (static_cast<int> (InputScalarViewType::rank) == 1,
                     "Template parameter InputScalarViewType must have "
                     "rank 1.");
      static_assert (std::is_same<
                       typename OutputScalarViewType::non_const_value_type,
                       typename InputScalarViewType::non_const_value_type>::value,
                     "Template parameters OutputScalarViewType and "
                     "InputScalarViewType must contain values of the same "
                     "type.");
      static_assert (std::is_same<
                       typename LocalIndicesViewType::non_const_value_type,
                       local_ordinal_type>::value,
                     "Template parameter LocalIndicesViewType must "
                     "contain values of type local_ordinal_type.");

      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;

      // Don't call this->hasColMap(), because that changes RCP's
      // reference count, which is not thread safe.  Just
      // dereferencing an RCP or calling RCP::is_null() does not
      // change its reference count.
      if (colMap_.is_null ()) {
        // No such thing as local column indices without a column Map.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }
      else if (newVals.dimension_0 () != inds.dimension_0 ()) {
        // The dimensions of the input arrays must match.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }

      const bool sorted = this->isSorted ();

      size_t hint = 0; // Guess for the current index k into rowVals
      LO numValid = 0; // number of valid local column indices

      // NOTE (mfh 11 Oct 2015) This method assumes UVM.  More
      // accurately, it assumes that the host execution space can
      // access data in both InputMemorySpace and ValsMemorySpace.

      if (isLocallyIndexed ()) {
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getLocalKokkosRowView (rowInfo);

        const LO numElts = static_cast<LO> (inds.dimension_0 ());
        for (LO j = 0; j < numElts; ++j) {
          const LO lclColInd = inds(j);
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         lclColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            if (atomic) {
              Kokkos::atomic_add (&rowVals(offset), newVals(j));
            }
            else {
              rowVals(offset) += newVals(j);
            }
            hint = offset + 1;
            ++numValid;
          }
        }
      }
      else if (isGloballyIndexed ()) {
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getGlobalKokkosRowView (rowInfo);

        const LO numElts = static_cast<LO> (inds.dimension_0 ());
        for (LO j = 0; j < numElts; ++j) {
          const GO gblColInd = this->colMap_->getGlobalElement (inds(j));
          if (gblColInd != Teuchos::OrdinalTraits<GO>::invalid ()) {
            const size_t offset =
              KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                           gblColInd, hint, sorted);
            if (offset != rowInfo.numEntries) {
              if (atomic) {
                Kokkos::atomic_add (&rowVals(offset), newVals(j));
              }
              else {
                rowVals(offset) += newVals(j);
              }
              hint = offset + 1;
              ++numValid;
            }
          }
        }
      }
      // NOTE (mfh 26 Jun 2014, 26 Nov 2015) In the current version of
      // CrsGraph and CrsMatrix, it's possible for a matrix (or graph)
      // to be neither locally nor globally indexed on a process.
      // This means that the graph or matrix has no entries on that
      // process.  Epetra also works like this.  It's related to lazy
      // allocation (on first insertion, not at graph / matrix
      // construction).  Lazy allocation will go away because it is
      // not thread scalable.

      return numValid;
    }

    /// \brief Implementation detail of CrsMatrix::replaceLocalValues.
    ///
    /// \tparam OutputScalarViewType Kokkos::View specialization that is
    ///   a 1-D array of the type of values in the sparse matrix (that
    ///   type is CrsMatrix::impl_scalar_type).
    /// \tparam LocalIndicesViewType Kokkos::View specialization that
    ///   is a 1-D array of LocalOrdinal.
    /// \tparam InputScalarViewType Kokkos::View specialization that
    ///   is a 1-D array of the type of values in the sparse matrix,
    ///   but with a possibly different memory space than how the
    ///   matrix stores its values.
    ///
    /// \param rowInfo [in] Result of getRowInfo on the index of the
    ///   local row of the matrix to modify.
    /// \param rowVals [in/out] On input: Values of the row of the
    ///   sparse matrix to modify.  On output: The modified values.
    /// \param inds [in] Local column indices of that row to modify.
    /// \param newVals [in] For each k, replace the value
    ///   corresponding to local column index inds[k] with newVals[k].
    ///
    /// \return The number of valid input column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    template<class OutputScalarViewType,
             class LocalIndicesViewType,
             class InputScalarViewType>
    LocalOrdinal
    replaceLocalValues (const RowInfo& rowInfo,
                        const typename UnmanagedView<OutputScalarViewType>::type& rowVals,
                        const typename UnmanagedView<LocalIndicesViewType>::type& inds,
                        const typename UnmanagedView<InputScalarViewType>::type& newVals) const
    {
      // We use static_assert here to check the template parameters,
      // rather than std::enable_if (e.g., on the return value, to
      // enable compilation only if the template parameters match the
      // desired attributes).  This turns obscure link errors into
      // clear compilation errors.  It also makes the return value a
      // lot easier to see.
      static_assert (Kokkos::is_view<OutputScalarViewType>::value,
                     "Template parameter OutputScalarViewType must be "
                     "a Kokkos::View.");
      static_assert (Kokkos::is_view<LocalIndicesViewType>::value,
                     "Template parameter LocalIndicesViewType must be "
                     "a Kokkos::View.");
      static_assert (Kokkos::is_view<InputScalarViewType>::value,
                     "Template parameter InputScalarViewType must be a "
                     "Kokkos::View.");
      static_assert (static_cast<int> (OutputScalarViewType::rank) == 1,
                     "Template parameter OutputScalarViewType must "
                     "have rank 1.");
      static_assert (static_cast<int> (LocalIndicesViewType::rank) == 1,
                     "Template parameter LocalIndicesViewType must "
                     "have rank 1.");
      static_assert (static_cast<int> (InputScalarViewType::rank) == 1,
                     "Template parameter InputScalarViewType must have "
                     "rank 1.");
      static_assert (std::is_same<
                       typename OutputScalarViewType::non_const_value_type,
                       typename InputScalarViewType::non_const_value_type>::value,
                     "Template parameters OutputScalarViewType and "
                     "InputScalarViewType must contain values of the same "
                     "type.");
      static_assert (std::is_same<
                       typename LocalIndicesViewType::non_const_value_type,
                       local_ordinal_type>::value,
                     "Template parameter LocalIndicesViewType must "
                     "contain values of type local_ordinal_type.");

      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;

      // Don't call this->hasColMap(), because that changes RCP's
      // reference count, which is not thread safe.  Just
      // dereferencing an RCP or calling RCP::is_null() does not
      // change its reference count.
      if (colMap_.is_null ()) {
        // No such thing as local column indices without a column Map.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }
      else if (newVals.dimension_0 () != inds.dimension_0 ()) {
        // The dimensions of the input arrays must match.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }

      const bool sorted = this->isSorted ();

      size_t hint = 0; // Guess for the current index k into rowVals
      LO numValid = 0; // number of valid local column indices

      // NOTE (mfh 11 Oct 2015) This method assumes UVM.  More
      // accurately, it assumes that the host execution space can
      // access data in all the Views.

      if (isLocallyIndexed ()) {
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getLocalKokkosRowView (rowInfo);

        const LO numElts = static_cast<LO> (inds.dimension_0 ());
        for (LO j = 0; j < numElts; ++j) {
          const LO lclColInd = inds(j);
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         lclColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            rowVals(offset) = newVals(j);
            hint = offset + 1;
            ++numValid;
          }
        }
      }
      else if (isGloballyIndexed ()) {
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getGlobalKokkosRowView (rowInfo);

        const LO numElts = static_cast<LO> (inds.dimension_0 ());
        for (LO j = 0; j < numElts; ++j) {
          const GO gblColInd = this->colMap_->getGlobalElement (inds(j));
          if (gblColInd != Teuchos::OrdinalTraits<GO>::invalid ()) {
            const size_t offset =
              KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                           gblColInd, hint, sorted);
            if (offset != rowInfo.numEntries) {
              rowVals(offset) = newVals(j);
              hint = offset + 1;
              ++numValid;
            }
          }
        }
      }
      // NOTE (mfh 26 Jun 2014, 26 Nov 2015) In the current version of
      // CrsGraph and CrsMatrix, it's possible for a matrix (or graph)
      // to be neither locally nor globally indexed on a process.
      // This means that the graph or matrix has no entries on that
      // process.  Epetra also works like this.  It's related to lazy
      // allocation (on first insertion, not at graph / matrix
      // construction).  Lazy allocation will go away because it is
      // not thread scalable.

      return numValid;
    }

    /// \brief Implementation detail of CrsMatrix::sumIntoGlobalValues.
    ///
    /// \tparam Scalar The type of each entry in the sparse matrix.
    /// \tparam InputMemorySpace Kokkos memory space / device in which
    ///   the input data live.  This may differ from the memory space
    ///   in which the current matrix values (rowVals) live.
    /// \tparam ValsMemorySpace Kokkos memory space / device in which
    ///   the matrix's current values live.  This may differ from the
    ///   memory space in which the input data (inds and newVals)
    ///   live.
    ///
    /// \param rowInfo [in] Result of getRowInfo on the index of the
    ///   local row of the matrix to modify.
    /// \param rowVals [in/out] On input: Values of the row of the
    ///   sparse matrix to modify.  On output: The modified values.
    /// \param inds [in] Global column indices of that row to modify.
    /// \param newVals [in] For each k, increment the value in rowVals
    ///   corresponding to global column index inds[k] by newVals[k].
    ///
    /// \return The number of valid input column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    template<class Scalar, class InputMemorySpace, class ValsMemorySpace>
    LocalOrdinal
    sumIntoGlobalValues (const RowInfo& rowInfo,
                         const Kokkos::View<Scalar*, ValsMemorySpace,
                           Kokkos::MemoryUnmanaged>& rowVals,
                         const Kokkos::View<const GlobalOrdinal*, InputMemorySpace,
                           Kokkos::MemoryUnmanaged>& inds,
                         const Kokkos::View<const Scalar*, InputMemorySpace,
                           Kokkos::MemoryUnmanaged>& newVals,
                         const bool atomic = useAtomicUpdatesByDefault) const
    {
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;

      if (newVals.dimension_0 () != inds.dimension_0 ()) {
        // The dimensions of the input arrays must match.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }

      const bool sorted = this->isSorted ();

      size_t hint = 0; // guess at the index's relative offset in the row
      LO numValid = 0; // number of valid input column indices

      // NOTE (mfh 11 Oct 2015) This method assumes UVM.  More
      // accurately, it assumes that the host execution space can
      // access data in both InputMemorySpace and ValsMemorySpace.

      if (isLocallyIndexed ()) {
        // NOTE (mfh 04 Nov 2015) Dereferencing an RCP or reading its
        // pointer does NOT change its reference count.  Thus, this
        // code is still thread safe.
        if (colMap_.is_null ()) {
          // NO input column indices are valid in this case, since if
          // the column Map is null on the calling process, then the
          // calling process owns no graph entries.
          return numValid;
        }
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getLocalKokkosRowView (rowInfo);
        const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();

        const LO numElts = static_cast<LO> (inds.dimension_0 ());
        for (LO j = 0; j < numElts; ++j) {
          const LO lclColInd = this->colMap_->getLocalElement (inds(j));
          if (lclColInd != LINV) {
            const size_t offset =
              KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                           lclColInd, hint, sorted);
            if (offset != rowInfo.numEntries) {
              if (atomic) {
                Kokkos::atomic_add (&rowVals(offset), newVals(j));
              }
              else {
                rowVals(offset) += newVals(j);
              }
              hint = offset + 1;
              numValid++;
            }
          }
        }
      }
      else if (isGloballyIndexed ()) {
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getGlobalKokkosRowView (rowInfo);

        const LO numElts = static_cast<LO> (inds.dimension_0 ());
        for (LO j = 0; j < numElts; ++j) {
          const GO gblColInd = inds(j);
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         gblColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            if (atomic) {
              Kokkos::atomic_add (&rowVals(offset), newVals(j));
            }
            else {
              rowVals(offset) += newVals(j);
            }
            hint = offset + 1;
            numValid++;
          }
        }
      }
      // If the graph is neither locally nor globally indexed on the
      // calling process, that means the calling process has no graph
      // entries.  Thus, none of the input column indices are valid.

      return numValid;
    }

    /// \brief Implementation detail of CrsMatrix::replaceGlobalValues.
    ///
    /// \tparam OutputScalarViewType Kokkos::View specialization that is
    ///   a 1-D array of the type of values in the sparse matrix (that
    ///   type is CrsMatrix::impl_scalar_type).
    /// \tparam GlobalIndicesViewType Kokkos::View specialization that
    ///   is a 1-D array of GlobalOrdinal.
    /// \tparam InputScalarViewType Kokkos::View specialization that
    ///   is a 1-D array of the type of values in the sparse matrix,
    ///   but with a possibly different memory space than how the
    ///   matrix stores its values.
    ///
    /// \param rowInfo [in] Result of getRowInfo on the index of the
    ///   local row of the matrix to modify.
    /// \param rowVals [in/out] On input: Values of the row of the
    ///   sparse matrix to modify.  On output: The modified values.
    /// \param inds [in] Global column indices of that row to modify.
    /// \param newVals [in] For each k, replace the value in rowVals
    ///   corresponding to global column index inds[k] with newVals[k].
    ///
    /// \return The number of valid input column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    template<class OutputScalarViewType,
             class GlobalIndicesViewType,
             class InputScalarViewType>
    LocalOrdinal
    replaceGlobalValues (const RowInfo& rowInfo,
                         const typename UnmanagedView<OutputScalarViewType>::type& rowVals,
                         const typename UnmanagedView<GlobalIndicesViewType>::type& inds,
                         const typename UnmanagedView<InputScalarViewType>::type& newVals) const
    {
      // We use static_assert here to check the template parameters,
      // rather than std::enable_if (e.g., on the return value, to
      // enable compilation only if the template parameters match the
      // desired attributes).  This turns obscure link errors into
      // clear compilation errors.  It also makes the return value a
      // lot easier to see.
      static_assert (Kokkos::is_view<OutputScalarViewType>::value,
                     "Template parameter OutputScalarViewType must be "
                     "a Kokkos::View.");
      static_assert (Kokkos::is_view<GlobalIndicesViewType>::value,
                     "Template parameter GlobalIndicesViewType must be "
                     "a Kokkos::View.");
      static_assert (Kokkos::is_view<InputScalarViewType>::value,
                     "Template parameter InputScalarViewType must be a "
                     "Kokkos::View.");
      static_assert (static_cast<int> (OutputScalarViewType::rank) == 1,
                     "Template parameter OutputScalarViewType must "
                     "have rank 1.");
      static_assert (static_cast<int> (GlobalIndicesViewType::rank) == 1,
                     "Template parameter GlobalIndicesViewType must "
                     "have rank 1.");
      static_assert (static_cast<int> (InputScalarViewType::rank) == 1,
                     "Template parameter InputScalarViewType must have "
                     "rank 1.");
      static_assert (std::is_same<
                       typename OutputScalarViewType::non_const_value_type,
                       typename InputScalarViewType::non_const_value_type>::value,
                     "Template parameters OutputScalarViewType and "
                     "InputScalarViewType must contain values of the same "
                     "type.");
      static_assert (std::is_same<
                       typename GlobalIndicesViewType::non_const_value_type,
                       global_ordinal_type>::value,
                     "Template parameter GlobalIndicesViewType must "
                     "contain values of type global_ordinal_type.");

      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;

      if (newVals.dimension_0 () != inds.dimension_0 ()) {
        // The dimensions of the input arrays must match.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }

      const bool sorted = this->isSorted ();

      size_t hint = 0; // guess at the index's relative offset in the row
      LO numValid = 0; // number of valid input column indices

      // NOTE (mfh 11 Oct 2015) This method assumes UVM.  More
      // accurately, it assumes that the host execution space can
      // access data in all the Views.

      if (isLocallyIndexed ()) {
        // NOTE (mfh 04 Nov 2015) Dereferencing an RCP or reading its
        // pointer does NOT change its reference count.  Thus, this
        // code is still thread safe.
        if (colMap_.is_null ()) {
          // NO input column indices are valid in this case, since if
          // the column Map is null on the calling process, then the
          // calling process owns no graph entries.
          return numValid;
        }
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getLocalKokkosRowView (rowInfo);
        const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();

        const LO numElts = static_cast<LO> (inds.dimension_0 ());
        for (LO j = 0; j < numElts; ++j) {
          const LO lclColInd = this->colMap_->getLocalElement (inds(j));
          if (lclColInd != LINV) {
            const size_t offset =
              KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                           lclColInd, hint, sorted);
            if (offset != rowInfo.numEntries) {
              rowVals(offset) = newVals(j);
              hint = offset + 1;
              numValid++;
            }
          }
        }
      }
      else if (isGloballyIndexed ()) {
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getGlobalKokkosRowView (rowInfo);

        const LO numElts = static_cast<LO> (inds.dimension_0 ());
        for (LO j = 0; j < numElts; ++j) {
          const GO gblColInd = inds(j);
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         gblColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            rowVals(offset) = newVals(j);
            hint = offset + 1;
            numValid++;
          }
        }
      }
      // If the graph is neither locally nor globally indexed on the
      // calling process, that means the calling process has no graph
      // entries.  Thus, none of the input column indices are valid.

      return numValid;
    }

    /// \brief Transform the given values using global indices.
    ///
    /// \tparam Scalar The type of each entry in the sparse matrix.
    /// \tparam BinaryFunction The type of the binary function f to
    ///   use for updating the sparse matrix's value(s).
    /// \tparam InputMemorySpace Kokkos memory space / device in which
    ///   the input data live.  This may differ from the memory space
    ///   in which the current matrix values (rowVals) live.
    /// \tparam ValsMemorySpace Kokkos memory space / device in which
    ///   the matrix's current values live.  This may differ from the
    ///   memory space in which the input data (inds and newVals)
    ///   live.
    ///
    /// \param rowInfo [in] Information about a given row of the graph.
    /// \param rowVals [in/out] The values to be transformed.  They
    ///   correspond to the row indicated by \c rowInfo.
    /// \param inds [in] The (global) indices in the row, for which
    ///   to transform the corresponding values in \c rowVals.
    /// \param newVals [in] Values to use for transforming \c rowVals.
    ///   These must NOT alias \c rowVals.
    /// \param f [in] A binary function used to transform \c rowVals.
    ///
    /// \return The number of valid local column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   <tt>Teuchos::OrdinalTraits<LocalOrdinal>::invalid()</tt>.
    template<class Scalar,
             class BinaryFunction,
             class InputMemorySpace,
             class ValsMemorySpace>
    LocalOrdinal
    transformGlobalValues (const RowInfo& rowInfo,
                           const Kokkos::View<Scalar*, ValsMemorySpace,
                             Kokkos::MemoryUnmanaged>& rowVals,
                           const Kokkos::View<const GlobalOrdinal*,
                             InputMemorySpace,
                             Kokkos::MemoryUnmanaged>& inds,
                           const Kokkos::View<const Scalar*,
                             InputMemorySpace,
                             Kokkos::MemoryUnmanaged>& newVals,
                           BinaryFunction f,
                           const bool atomic = useAtomicUpdatesByDefault) const
    {
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;

      if (newVals.dimension_0 () != inds.dimension_0 ()) {
        // The sizes of the input arrays must match.
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }

      const LO numElts = static_cast<LO> (inds.dimension_0 ());
      const bool sorted = this->isSorted ();

      LO numValid = 0; // number of valid input column indices
      size_t hint = 0; // guess at the index's relative offset in the row

      if (isLocallyIndexed ()) {
        // NOTE (mfh 04 Nov 2015) Dereferencing an RCP or reading its
        // pointer does NOT change its reference count.  Thus, this
        // code is still thread safe.
        if (colMap_.is_null ()) {
          // NO input column indices are valid in this case, since if
          // the column Map is null on the calling process, then the
          // calling process owns no graph entries.
          return numValid;
        }
        const map_type& colMap = *colMap_;
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getLocalKokkosRowView (rowInfo);

        const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();
        for (LO j = 0; j < numElts; ++j) {
          const LO lclColInd = colMap.getLocalElement (inds(j));
          if (lclColInd != LINV) {
            const size_t offset =
              KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                           lclColInd, hint, sorted);
            if (offset != rowInfo.numEntries) {
              if (atomic) {
                // NOTE (mfh 30 Nov 2015) The commented-out code is
                // wrong because another thread may have changed
                // rowVals(offset) between those two lines of code.
                //
                //const Scalar newVal = f (rowVals(offset), newVals(j));
                //Kokkos::atomic_assign (&rowVals(offset), newVal);

                volatile Scalar* const dest = &rowVals(offset);
                (void) atomic_binary_function_update (dest, newVals(j), f);
              }
              else {
                // use binary function f
                rowVals(offset) = f (rowVals(offset), newVals(j));
              }
              hint = offset + 1;
              numValid++;
            }
          }
        }
      }
      else if (isGloballyIndexed ()) {
        // Get a view of the column indices in the row.  This amortizes
        // the cost of getting the view over all the entries of inds.
        auto colInds = this->getGlobalKokkosRowView (rowInfo);

        for (LO j = 0; j < numElts; ++j) {
          const GO gblColInd = inds(j);
          const size_t offset =
            KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                         gblColInd, hint, sorted);
          if (offset != rowInfo.numEntries) {
            if (atomic) {
              // NOTE (mfh 30 Nov 2015) The commented-out code is
              // wrong because another thread may have changed
              // rowVals(offset) between those two lines of code.
              //
              //const Scalar newVal = f (rowVals(offset), newVals(j));
              //Kokkos::atomic_assign (&rowVals(offset), newVal);

              volatile Scalar* const dest = &rowVals(offset);
              (void) atomic_binary_function_update (dest, newVals(j), f);
            }
            else {
              // use binary function f
              rowVals(offset) = f (rowVals(offset), newVals(j));
            }
            hint = offset + 1;
            numValid++;
          }
        }
      }
      // If the graph is neither locally nor globally indexed on the
      // calling process, that means the calling process has no graph
      // entries.  Thus, none of the input column indices are valid.

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
    void
    sortRowIndicesAndValues (const RowInfo rowinfo,
                             const Teuchos::ArrayView<Scalar>& values)
    {
      if (rowinfo.numEntries > 0) {
        Teuchos::ArrayView<LocalOrdinal> inds_view =
          this->getLocalViewNonConst (rowinfo);
        sort2 (inds_view.begin (), inds_view.begin () + rowinfo.numEntries,
               values.begin ());
      }
    }

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
                              const Teuchos::ArrayView<Scalar>& rowValues)
    {
      using Teuchos::ArrayView;
      const char tfecfFuncName[] = "mergeRowIndicesAndValues: ";
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (isStorageOptimized(), std::logic_error, "It is invalid to call this "
         "method if the graph's storage has already been optimized.  Please "
         "report this bug to the Tpetra developers.");

      typedef typename ArrayView<Scalar>::iterator Iter;
      Iter rowValueIter = rowValues.begin ();
      ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst (rowinfo);
      typename ArrayView<LocalOrdinal>::iterator beg, end, newend;

      // beg,end define a half-exclusive interval over which to iterate.
      beg = inds_view.begin();
      end = inds_view.begin() + rowinfo.numEntries;
      newend = beg;
      if (beg != end) {
        typename ArrayView<LocalOrdinal>::iterator cur = beg + 1;
        Iter vcur = rowValueIter + 1;
        Iter vend = rowValueIter;
        cur = beg+1;
        while (cur != end) {
          if (*cur != *newend) {
            // new entry; save it
            ++newend;
            ++vend;
            (*newend) = (*cur);
            (*vend) = (*vcur);
          }
          else {
            // old entry; merge it
            //(*vend) = f (*vend, *vcur);
            (*vend) += *vcur;
          }
          ++cur;
          ++vcur;
        }
        ++newend; // one past the last entry, per typical [beg,end) semantics
      }
      const size_t mergedEntries = newend - beg;
#ifdef HAVE_TPETRA_DEBUG
      // merge should not have eliminated any entries; if so, the
      // assignment below will destroy the packed structure
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (isStorageOptimized() && mergedEntries != rowinfo.numEntries,
         std::logic_error,
         "Merge was incorrect; it eliminated entries from the graph.  "
         << "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

      k_numRowEntries_(rowinfo.localRow) = mergedEntries;
      nodeNumEntries_ -= (rowinfo.numEntries - mergedEntries);
    }

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
    RowInfo getRowInfo (const LocalOrdinal myRow) const;

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
    RowInfo getRowInfoFromGlobalRowIndex (const GlobalOrdinal gblRow) const;

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

    /// \brief Get a pointer to the local column indices of a locally
    ///   owned row, using the result of getRowInfo.
    ///
    /// \param lclInds [out] Pointer to the local column indices of
    ///   the given row.
    /// \param numEnt [out] Number of entries in the given row.
    /// \param rowinfo [in] Result of getRowInfo(lclRow) for the row
    ///   \c lclRow to view.
    ///
    /// \return 0 if successful, else a nonzero error code.
    LocalOrdinal
    getLocalViewRawConst (const LocalOrdinal*& lclInds,
                          LocalOrdinal& numEnt,
                          const RowInfo& rowinfo) const;

  private:

    /// \brief Get a const nonowned view of the local column indices
    ///   indices of row rowinfo.localRow (only works if the matrix is
    ///   locally indexed on the calling process).
    ///
    /// \param rowinfo [in] Result of calling getRowInfo with the
    ///   index of the local row to view.
    Kokkos::View<const LocalOrdinal*, execution_space, Kokkos::MemoryUnmanaged>
    getLocalKokkosRowView (const RowInfo& rowinfo) const;

    /// \brief Get a nonconst nonowned view of the local column
    ///   indices of row rowinfo.localRow (only works if the matrix is
    ///   locally indexed on the calling process).
    ///
    /// \param rowinfo [in] Result of calling getRowInfo with the
    ///   index of the local row to view.
    Kokkos::View<LocalOrdinal*, execution_space, Kokkos::MemoryUnmanaged>
    getLocalKokkosRowViewNonConst (const RowInfo& rowinfo);

    /// \brief Get a const nonowned view of the global column indices
    ///   of row rowinfo.localRow (only works if the matrix is
    ///   globally indexed).
    ///
    /// \param rowinfo [in] Result of calling getRowInfo with the
    ///   index of the local row to view.
    Kokkos::View<const GlobalOrdinal*, execution_space, Kokkos::MemoryUnmanaged>
    getGlobalKokkosRowView (const RowInfo& rowinfo) const;

  protected:

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

    /// \brief Get a pointer to the global column indices of a locally
    ///   owned row, using the result of getRowInfoFromGlobalRowIndex.
    ///
    /// \param gblInds [out] Pointer to the global column indices of
    ///   the given row.
    /// \param numEnt [out] Number of entries in the given row.
    /// \param rowinfo [in] Result of
    ///   getRowInfoFromGlobalRowIndex(gblRow) for the row to view,
    ///   whose global row index is \c gblRow.
    ///
    /// \return 0 if successful, else a nonzero error code.
    LocalOrdinal
    getGlobalViewRawConst (const GlobalOrdinal*& gblInds,
                           LocalOrdinal& numEnt,
                           const RowInfo& rowinfo) const;

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

    /// \brief The number of local entries in each locally owned row.
    ///
    /// This is deallocated in fillComplete() if fillComplete()'s
    /// "Optimize Storage" parameter is set to \c true.
    ///
    /// This may also exist with 1-D storage, if storage is unpacked.
    typename Kokkos::View<size_t*, Kokkos::LayoutLeft, device_type>::HostMirror
      k_numRowEntries_;
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

  /// \brief Nonmember function to create an empty CrsGraph given a
  ///   row Map and the max number of entries allowed locally per row.
  ///
  /// \return A dynamically allocated (DynamicProfile) graph with
  ///   specified number of nonzeros per row (defaults to zero).
  /// \relatesalso CrsGraph
  template <class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic = Node::classic>
  Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic> >
  createCrsGraph (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map,
                  size_t maxNumEntriesPerRow = 0,
                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
  {
    using Teuchos::rcp;
    typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node, classic> graph_type;
    return rcp (new graph_type (map, maxNumEntriesPerRow, DynamicProfile, params));
  }

  namespace Details {

    template<class LocalOrdinal,
             class GlobalOrdinal,
             class OutputNodeType,
             class InputNodeType>
    class CrsGraphCopier<CrsGraph<LocalOrdinal, GlobalOrdinal, OutputNodeType>,
                         CrsGraph<LocalOrdinal, GlobalOrdinal, InputNodeType> > {
    public:
      typedef CrsGraph<LocalOrdinal, GlobalOrdinal, InputNodeType> input_crs_graph_type;
      typedef CrsGraph<LocalOrdinal, GlobalOrdinal, OutputNodeType> output_crs_graph_type;

      static Teuchos::RCP<output_crs_graph_type>
      clone (const input_crs_graph_type& graphIn,
             const Teuchos::RCP<OutputNodeType> &nodeOut,
             const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
      {
        using Teuchos::arcp;
        using Teuchos::ArrayRCP;
        using Teuchos::ArrayView;
        using Teuchos::null;
        using Teuchos::outArg;
        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::REDUCE_MIN;
        using Teuchos::reduceAll;
        using Teuchos::sublist;
        using std::cerr;
        using std::endl;
        typedef LocalOrdinal LO;
        typedef GlobalOrdinal GO;
        typedef typename ArrayView<const GO>::size_type size_type;
        typedef ::Tpetra::Map<LO, GO, InputNodeType> input_map_type;
        typedef ::Tpetra::Map<LO, GO, OutputNodeType> output_map_type;
        const char prefix[] = "Tpetra::Details::CrsGraphCopier::clone: ";

        // Set parameters' default values.
        bool debug = false;
        bool fillCompleteClone = true;
        bool useLocalIndices = graphIn.hasColMap ();
        ProfileType pftype = StaticProfile;
        // If the user provided a ParameterList, get values from there.
        if (! params.is_null ()) {
          fillCompleteClone = params->get ("fillComplete clone", fillCompleteClone);
          useLocalIndices = params->get ("Locally indexed clone", useLocalIndices);
          if (params->get ("Static profile clone", true) == false) {
            pftype = DynamicProfile;
          }
          debug = params->get ("Debug", debug);
        }

        const Teuchos::Comm<int>& comm = * (graphIn.getRowMap ()->getComm ());
        const int myRank = comm.getRank ();

        TEUCHOS_TEST_FOR_EXCEPTION(
                                   ! graphIn.hasColMap () && useLocalIndices, std::runtime_error,
                                   prefix << "You asked clone() to use local indices (by setting the "
                                   "\"Locally indexed clone\" parameter to true), but the source graph "
                                   "does not yet have a column Map, so this is impossible.");

        if (debug) {
          std::ostringstream os;
          os << "Process " << myRank << ": Cloning row Map" << endl;
          cerr << os.str ();
        }

        RCP<const output_map_type> clonedRowMap =
          graphIn.getRowMap ()->template clone<OutputNodeType> (nodeOut);

        // Invoke the output graph's constructor, using the input graph's
        // upper bounds on the number of entries in each local row.
        RCP<output_crs_graph_type> clonedGraph; // returned by this function
        {
          ArrayRCP<const size_t> numEntriesPerRow;
          size_t numEntriesForAll = 0;
          bool boundSameForAllLocalRows = true;

          if (debug) {
            std::ostringstream os;
            os << "Process " << myRank << ": Getting per-row bounds" << endl;
            cerr << os.str ();
          }
          graphIn.getNumEntriesPerLocalRowUpperBound (numEntriesPerRow,
                                                      numEntriesForAll,
                                                      boundSameForAllLocalRows);
          if (debug) {
            std::ostringstream os;
            os << "Process " << myRank << ": numEntriesForAll = "
               << numEntriesForAll << endl;
            cerr << os.str ();
          }

          if (debug) {
            std::ostringstream os;
            os << "Process " << myRank << ": graphIn.getNodeMaxNumRowEntries() = "
               << graphIn.getNodeMaxNumRowEntries () << endl;
            cerr << os.str ();
          }

          RCP<ParameterList> graphparams;
          if (params.is_null ()) {
            graphparams = parameterList ("CrsGraph");
          } else {
            graphparams = sublist (params, "CrsGraph");
          }
          if (useLocalIndices) {
            RCP<const output_map_type> clonedColMap =
              graphIn.getColMap ()->template clone<OutputNodeType> (nodeOut);
            if (boundSameForAllLocalRows) {
              clonedGraph = rcp (new output_crs_graph_type (clonedRowMap, clonedColMap,
                                                            numEntriesForAll, pftype,
                                                            graphparams));
            } else {
              clonedGraph = rcp (new output_crs_graph_type (clonedRowMap, clonedColMap,
                                                            numEntriesPerRow, pftype,
                                                            graphparams));
            }
          } else {
            if (boundSameForAllLocalRows) {
              clonedGraph = rcp (new output_crs_graph_type (clonedRowMap,
                                                            numEntriesForAll, pftype,
                                                            graphparams));
            } else {
              clonedGraph = rcp (new output_crs_graph_type (clonedRowMap,
                                                            numEntriesPerRow,
                                                            pftype, graphparams));
            }
          }

          if (debug) {
            std::ostringstream os;
            os << "Process " << myRank << ": Invoked output graph's constructor" << endl;
            cerr << os.str ();
          }

          // done with these
          numEntriesPerRow = null;
          numEntriesForAll = 0;
        }

        const input_map_type& inputRowMap = * (graphIn.getRowMap ());
        const size_type numRows =
          static_cast<size_type> (inputRowMap.getNodeNumElements ());

        bool failed = false;

        if (useLocalIndices) {
          const LO localMinLID = inputRowMap.getMinLocalIndex ();
          const LO localMaxLID = inputRowMap.getMaxLocalIndex ();

          if (graphIn.isLocallyIndexed ()) {
            if (numRows != 0) {
              try {
                ArrayView<const LO> linds;
                for (LO lrow = localMinLID; lrow <= localMaxLID; ++lrow) {
                  graphIn.getLocalRowView (lrow, linds);
                  if (linds.size () != 0) {
                    clonedGraph->insertLocalIndices (lrow, linds);
                  }
                }
              }
              catch (std::exception& e) {
                std::ostringstream os;
                os << "Process " << myRank << ": copying (reading local by view, "
                  "writing local) indices into the output graph threw an "
                  "exception: " << e.what () << endl;
                cerr << os.str ();
                failed = true;
              }
            }
          }
          else { // graphIn.isGloballyIndexed()
            TEUCHOS_TEST_FOR_EXCEPTION(
                                       ! graphIn.hasColMap () && useLocalIndices, std::invalid_argument,
                                       prefix << "You asked clone() to use local indices (by setting the "
                                       "\"Locally indexed clone\" parameter to true), but the source graph "
                                       "does not yet have a column Map, so this is impossible.");

            // The input graph has a column Map, but is globally indexed.
            // That's a bit weird, but we'll run with it.  In this case,
            // getLocalRowView won't work, but getLocalRowCopy should
            // still work; it will just have to convert from global to
            // local indices internally.

            try {
              // Make space for getLocalRowCopy to put column indices.
              //
              // This is only a hint; we may have to resize in the loop
              // below.  getNodeMaxNumRowEntries() may return nonsense if
              // fill is active.  The key bool in CrsGraph is
              // haveLocalConstants_.
              size_t myMaxNumRowEntries =
                graphIn.isFillActive () ? static_cast<size_t> (0) :
                graphIn.getNodeMaxNumRowEntries ();

              Array<LO> linds (myMaxNumRowEntries);

              // Copy each row into the new graph, using local indices.
              for (LO lrow = localMinLID; lrow <= localMaxLID; ++lrow) {
                size_t theNumEntries = graphIn.getNumEntriesInLocalRow (lrow);
                if (theNumEntries > myMaxNumRowEntries) {
                  myMaxNumRowEntries = theNumEntries;
                  linds.resize (myMaxNumRowEntries);
                }
                graphIn.getLocalRowCopy (lrow, linds (), theNumEntries);
                if (theNumEntries != 0) {
                  clonedGraph->insertLocalIndices (lrow, linds (0, theNumEntries));
                }
              }
            }
            catch (std::exception& e) {
              std::ostringstream os;
              os << "Process " << myRank << ": copying (reading local by copy, "
                "writing local) indices into the output graph threw an exception: "
                 << e.what () << endl;
              cerr << os.str ();
              failed = true;
            }
          }
        }
        else { /* useGlobalIndices */
          if (numRows != 0) {
            const GlobalOrdinal localMinGID = inputRowMap.getMinGlobalIndex ();
            const GlobalOrdinal localMaxGID = inputRowMap.getMaxGlobalIndex ();
            const bool inputRowMapIsContiguous = inputRowMap.isContiguous ();

            if (graphIn.isGloballyIndexed ()) {
              ArrayView<const GlobalOrdinal> ginds;

              if (inputRowMapIsContiguous) {
                try {
                  for (GO grow = localMinGID; grow <= localMaxGID; ++grow) {
                    graphIn.getGlobalRowView (grow, ginds);
                    if (ginds.size () != 0) {
                      clonedGraph->insertGlobalIndices (grow, ginds);
                    }
                  }
                }
                catch (std::exception& e) {
                  std::ostringstream os;
                  os << "Process " << myRank << ": copying (reading global by view, "
                    "writing global) indices into the output graph threw an "
                    "exception: " << e.what () << endl;
                  cerr << os.str ();
                  failed = true;
                }
              }
              else { // input row Map is not contiguous
                try {
                  ArrayView<const GO> inputRowMapGIDs = inputRowMap.getNodeElementList ();
                  for (size_type k = 0; k < numRows; ++k) {
                    const GO grow = inputRowMapGIDs[k];
                    graphIn.getGlobalRowView (grow, ginds);
                    if (ginds.size () != 0) {
                      clonedGraph->insertGlobalIndices (grow, ginds);
                    }
                  }
                }
                catch (std::exception& e) {
                  std::ostringstream os;
                  os << "Process " << myRank << ": copying (reading global by view, "
                    "writing global) indices into the output graph threw an "
                    "exception: " << e.what () << endl;
                  cerr << os.str ();
                  failed = true;
                }
              }
            }
            else { // graphIn.isLocallyIndexed()
              // Make space for getGlobalRowCopy to put column indices.
              //
              // This is only a hint; we may have to resize in the loop
              // below.  getNodeMaxNumRowEntries() may return nonsense if
              // fill is active.  The key bool in CrsGraph is
              // haveLocalConstants_.
              size_t myMaxNumRowEntries =
                graphIn.isFillActive () ? static_cast<size_t> (0) :
                graphIn.getNodeMaxNumRowEntries ();

              Array<GO> ginds (myMaxNumRowEntries);

              if (inputRowMapIsContiguous) {
                try {
                  for (GO grow = localMinGID; grow <= localMaxGID; ++grow) {
                    size_t theNumEntries = graphIn.getNumEntriesInGlobalRow (grow);
                    if (theNumEntries > myMaxNumRowEntries) {
                      myMaxNumRowEntries = theNumEntries;
                      ginds.resize (myMaxNumRowEntries);
                    }
                    graphIn.getGlobalRowCopy (grow, ginds (), theNumEntries);
                    if (theNumEntries != 0) {
                      clonedGraph->insertGlobalIndices (grow, ginds (0, theNumEntries));
                    }
                  }
                }
                catch (std::exception& e) {
                  std::ostringstream os;
                  os << "Process " << myRank << ": copying (reading global by copy, "
                    "writing global) indices into the output graph threw an "
                    "exception: " << e.what () << endl;
                  cerr << os.str ();
                  failed = true;
                }
              }
              else { // input row Map is not contiguous
                try {
                  ArrayView<const GO> inputRowMapGIDs = inputRowMap.getNodeElementList ();
                  for (size_type k = 0; k < numRows; ++k) {
                    const GO grow = inputRowMapGIDs[k];

                    size_t theNumEntries = graphIn.getNumEntriesInGlobalRow (grow);
                    if (theNumEntries > myMaxNumRowEntries) {
                      myMaxNumRowEntries = theNumEntries;
                      ginds.resize (myMaxNumRowEntries);
                    }
                    graphIn.getGlobalRowCopy (grow, ginds (), theNumEntries);
                    if (theNumEntries != 0) {
                      clonedGraph->insertGlobalIndices (grow, ginds (0, theNumEntries));
                    }
                  }
                }
                catch (std::exception& e) {
                  std::ostringstream os;
                  os << "Process " << myRank << ": copying (reading global by copy, "
                    "writing global) indices into the output graph threw an "
                    "exception: " << e.what () << endl;
                  cerr << os.str ();
                  failed = true;
                }
              }
            }
          } // numRows != 0
        }

        if (debug) {
          std::ostringstream os;
          os << "Process " << myRank << ": copied entries" << endl;
          cerr << os.str ();
        }

        if (fillCompleteClone) {
          RCP<ParameterList> fillparams = params.is_null () ?
            parameterList ("fillComplete") :
            sublist (params, "fillComplete");
          try {
            RCP<const output_map_type> clonedRangeMap;
            RCP<const output_map_type> clonedDomainMap;
            if (! graphIn.getRangeMap ().is_null () &&
                graphIn.getRangeMap () != graphIn.getRowMap ()) {
              clonedRangeMap =
                graphIn.getRangeMap ()->template clone<OutputNodeType> (nodeOut);
            }
            else {
              clonedRangeMap = clonedRowMap;
            }
            if (! graphIn.getDomainMap ().is_null ()
                && graphIn.getDomainMap () != graphIn.getRowMap ()) {
              clonedDomainMap =
                graphIn.getDomainMap ()->template clone<OutputNodeType> (nodeOut);
            }
            else {
              clonedDomainMap = clonedRowMap;
            }

            if (debug) {
              std::ostringstream os;
              os << "Process " << myRank << ": About to call fillComplete on "
                "cloned graph" << endl;
              cerr << os.str ();
            }
            clonedGraph->fillComplete (clonedDomainMap, clonedRangeMap, fillparams);
          }
          catch (std::exception &e) {
            failed = true;
            std::ostringstream os;
            os << prefix << "Process " << myRank << ": Caught the following "
              "exception while calling fillComplete() on clone of type"
               << endl << Teuchos::typeName (*clonedGraph) << endl;
            cerr << os.str ();
          }
        }

        int lclSuccess = failed ? 0 : 1;
        int gblSuccess = 1;
        reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
        TEUCHOS_TEST_FOR_EXCEPTION(
                                   gblSuccess != 1, std::logic_error, prefix <<
                                   "Clone failed on at least one process.");

        if (debug) {
          std::ostringstream os;
          os << "Process " << myRank << ": Done with CrsGraph::clone" << endl;
          cerr << os.str ();
        }
        return clonedGraph;
      }
    };

  } // namespace Details
} // namespace Tpetra

#endif // TPETRA_CRSGRAPH_DECL_HPP
