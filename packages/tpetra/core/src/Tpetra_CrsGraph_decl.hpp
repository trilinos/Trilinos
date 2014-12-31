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

#include <Tpetra_ConfigDefs.hpp>

#include <Tpetra_RowGraph.hpp>
#include <Tpetra_DistObject.hpp>
#include <Tpetra_Exceptions.hpp>

#include <Kokkos_DefaultKernels.hpp>

#include <Teuchos_CompileTimeAssert.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>

#ifndef KOKKOS_CONFIGDEFS_HPP
#  error "KOKKOS_CONFIGDEFS_HPP not defined, so something is wrong with Kokkos_ConfigDefs.hpp"
#endif // KOKKOS_CONFIGDEFS_HPP

#ifndef KOKKOS_DEFAULT_KERNELS_HPP
#  error "KOKKOS_DEFAULT_KERNELS_HPP not defined, so something is wrong with Kokkos_ConfigDefs.hpp"
#endif // KOKKOS_DEFAULT_KERNELS_HPP

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  //
  // Dear users: These are just forward declarations.  Please skip
  // over them and go down to the CrsMatrix class declaration.  Thank
  // you.
  //
  template <class LO, class GO, class N, const bool isClassic>
  class CrsGraph;

  template <class S, class LO, class GO, class N>
  class CrsMatrix;

  namespace Experimental {
    template<class S, class LO, class GO, class N>
    class BlockCrsMatrix;
  }

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
  template <class LocalOrdinal = RowGraph<>::local_ordinal_type,
            class GlobalOrdinal = typename RowGraph<LocalOrdinal>::global_ordinal_type,
            class Node = typename RowGraph<LocalOrdinal, GlobalOrdinal>::node_type,
            const bool classic = Node::classic>
  class CrsGraph {
    // See partial specializations for documentation of methods.
  };

#if defined(HAVE_TPETRACLASSIC_SERIAL) || defined(HAVE_TPETRACLASSIC_TBB) || defined(HAVE_TPETRACLASSIC_THREADPOOL) || defined(HAVE_TPETRACLASSIC_OPENMP) || defined(HAVE_TPETRACLASSIC_THRUST)

  //! Partial specialization for the "classic" version of Tpetra.
  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class CrsGraph<LocalOrdinal, GlobalOrdinal, Node, true> :
    public RowGraph<LocalOrdinal,GlobalOrdinal,Node>,
    public DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>,
    public Teuchos::ParameterListAcceptorDefaultBase
  {
    template <class S, class LO, class GO, class N>
    friend class CrsMatrix;
    template <class S, class LO, class GO, class N>
    friend class Experimental::BlockCrsMatrix;
    template <class LO2, class GO2, class N2, const bool isClassic>
    friend class CrsGraph;
    template<class OutputCrsGraphType, class InputCrsGraphType>
    friend class Details::CrsGraphCopier;

  public:
    //! This class' first template parameter; the type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! This class' second template parameter; the type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! This class' third template parameter; the Kokkos Node type.
    typedef Node node_type;

    //! The Map specialization used by this class.
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, node_type> map_type;
    //! The Import specialization used by this class.
    typedef Tpetra::Import<LocalOrdinal, GlobalOrdinal, node_type> import_type;
    //! The Export specialization used by this class.
    typedef Tpetra::Export<LocalOrdinal, GlobalOrdinal, node_type> export_type;

    //! @name Constructor/Destructor Methods
    //@{

    /// \brief Constructor specifying fixed number of entries for each row.
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
              const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying (possibly different) number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of graph
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
              ProfileType pftype = DynamicProfile,
              const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying column Map and fixed number of entries for each row.
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
              size_t maxNumEntriesPerRow,
              ProfileType pftype = DynamicProfile,
              const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying column Map and number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the graph.
    ///
    /// \param colMap [in] Distribution of columns of the graph.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of graph
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
              const Teuchos::RCP<const map_type>& colMap,
              const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc,
              ProfileType pftype = DynamicProfile,
              const Teuchos::RCP<Teuchos::ParameterList>& params = null);

    /// \brief Constructor specifying column Map and arrays containing
    ///   the graph in sorted local indices.
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
           const Teuchos::RCP<Teuchos::ParameterList> &params = null) const
    {
      typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node2> output_crs_graph_type;
      typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node> input_crs_graph_type;
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
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters () const;

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
    void resumeFill (const Teuchos::RCP<Teuchos::ParameterList> &params = null);

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
                  const Teuchos::RCP<Teuchos::ParameterList> &params = null);

    /*! \brief Signal that data entry is complete.

      Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

      \note This method calls fillComplete( getRowMap(), getRowMap(), os ). See parameter options there.
    */
    void fillComplete (const Teuchos::RCP<Teuchos::ParameterList> &params = null);

    /// \brief Perform a fillComplete on a graph that already has data, via setAllIndices().
    ///
    /// The graph must already have filled local 1-D storage (lclInds1D_
    /// and rowPtrs_).  If the graph has been constructed in any other way,
    /// this method will throw an exception.  This routine is needed to
    /// support other Trilinos packages and should not be called by ordinary
    /// users.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    expertStaticFillComplete (const Teuchos::RCP<const map_type> & domainMap,
                              const Teuchos::RCP<const map_type> & rangeMap,
                              const Teuchos::RCP<const import_type>& importer = Teuchos::null,
                              const Teuchos::RCP<const export_type>& exporter = Teuchos::null,
                              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
    //@}
    //! @name Methods implementing RowGraph.
    //@{

    //! Returns the communicator.
    Teuchos::RCP<const Comm<int> > getComm() const;

    //! Returns the underlying node.
    Teuchos::RCP<Node> getNode() const;

    //! Returns the Map that describes the row distribution in this graph.
    Teuchos::RCP<const map_type> getRowMap() const;

    //! \brief Returns the Map that describes the column distribution in this graph.
    Teuchos::RCP<const map_type> getColMap() const;

    //! Returns the Map associated with the domain of this graph.
    Teuchos::RCP<const map_type> getDomainMap() const;

    //! Returns the Map associated with the domain of this graph.
    Teuchos::RCP<const map_type> getRangeMap() const;

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

    /// \brief Whether column indices are stored using local indices on the calling process.
    ///
    /// If the graph stores column indices as local indices on the
    /// calling process, this returns true.  Otherwise, it returns
    /// false.  Exactly one of the following is true:
    /// <ol>
    ///   <li> <tt>! isLocallyIndexed() && ! isGloballyIndexed()</tt> </li>
    ///   <li> <tt>! isLocallyIndexed() && isGloballyIndexed()</tt> </li>
    ///   <li> <tt>isLocallyIndexed() && ! isGloballyIndexed()</tt> </li>
    /// </ol>
    /// The first condition means that the calling process does not
    /// own any graph entries.  The second condition means that the
    /// graph does not (yet) have a column Map, so no process can be
    /// locally indexed.  The third condition means that the graph has
    /// a column Map.
    bool isLocallyIndexed() const;

    /// \brief Whether column indices are stored using global indices on the calling process.
    ///
    /// If the graph stores column indices as global indices on the
    /// calling process, this returns true.  Otherwise, it returns
    /// false.  Exactly one of the following is true:
    /// <ol>
    ///   <li> <tt>! isLocallyIndexed() && ! isGloballyIndexed()</tt> </li>
    ///   <li> <tt>! isLocallyIndexed() && isGloballyIndexed()</tt> </li>
    ///   <li> <tt>isLocallyIndexed() && ! isGloballyIndexed()</tt> </li>
    /// </ol>
    /// The first condition means that the calling process does not
    /// own any graph entries.  The second condition means that the
    /// graph does not (yet) have a column Map, so no process can be
    /// locally indexed.  The third condition means that the graph has
    /// a column Map.
    bool isGloballyIndexed() const;

    /// \brief Whether fillComplete() has been called and the graph is in compute mode.
    ///
    /// This state is the same on all processes in the graph's
    /// communicator.  Exactly one of the following is true:
    /// <ol>
    ///   <li> <tt>isFillActive() && ! isFillComplete()</tt> </li>
    ///   <li> <tt>! isFillActive() && isFillComplete()</tt> </li>
    /// </ol>
    bool isFillComplete() const;

    /// \brief Whether resumeFill() has been called and the graph is in edit mode.
    ///
    /// This method returns true under the following conditions:
    ///   - The graph has just been created, and fillComplete() has
    ///     not yet been called.
    ///   - resumeFill() has been called without an intervening call
    ///     to fillComplete().
    ///
    /// This state is the same on all processes in the graph's
    /// communicator.  Exactly one of the following is true:
    /// <ol>
    ///   <li> <tt>isFillActive() && ! isFillComplete()</tt> </li>
    ///   <li> <tt>! isFillActive() && isFillComplete()</tt> </li>
    /// </ol>
    bool isFillActive() const;

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
    bool isStorageOptimized() const;

    //! Returns \c true if the graph was allocated with static data structures.
    ProfileType getProfileType() const;

    /// \brief Get a copy of the column indices (as global indices) in the given row.
    ///
    /// \param GlobalRow [in] Global index of the row.
    ///
    /// \param Indices [out] Global column indices owned by the
    ///   calling process in row \c GlobalRow.  This is a view of an
    ///   array that you, the caller, must allocate in advance.  You
    ///   may call getNumEntriesInGlobalRow() to find the number of
    ///   entries in this row on the calling process, or
    ///   getNumAllocatedEntriesInGlobalRow() to get the maximum
    ///   number of entries in all rows on the calling process.  (The
    ///   latter may be useful when looping over all rows on the
    ///   calling process.)
    ///
    /// \param NumIndices [out] The number of column indices stored in
    ///   \c Indices.  If \c GlobalRow does not belong to (i.e., is
    ///   not in the row Map on) the calling process, then \c Indices
    ///   is unchanged and \c NumIndices is returned as
    ///   <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>.
    ///
    /// This method only returns the column indices owned by the
    /// calling process.  This is a <i>local</i> operation with
    /// respect to MPI interprocess communication.  It could be the
    /// case that multiple processes store different sets of column
    /// indices in the same row.  In that case, this method only
    /// returns the column indices stored on the calling process.  Use
    /// the row Map to determine which process(es) own which rows.
    ///
    /// This method throws std::runtime_error if the output array
    /// <tt>Indices</tt> is not large enough to hold the column
    /// indices in the given row.
    ///
    /// You may call this method at any time, whether or not the graph
    /// has a column Map or is fill complete.
    ///
    /// \note Stay tuned for the following possible future changes to
    ///   this method.  First, this method will return the number of
    ///   indices in the row, and the output argument \c NumIndices
    ///   will go away.  It will be the caller's responsibility to
    ///   test whether this is less than the size of \c Indices.
    ///   Second, instead of throwing an exception if \c GlobalRow is
    ///   not owned by the calling process, this method will simply
    ///   return zero.  (This indicates that the calling process does
    ///   not own any column indices in that row, which is true even
    ///   if the row index is not in the row Map on any process.)
    ///   Third, this method will take a Kokkos::View instead of a
    ///   Teuchos::ArrayView.
    void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      const Teuchos::ArrayView<GlobalOrdinal>& Indices,
                      size_t& NumIndices) const;

    /// \brief Get a copy of the column indices (as local indices) in the given row.
    ///
    /// \param LocalRow [in] Local index of the row.
    ///
    /// \param indices [out] Local column indices owned by the calling
    ///   process in row \c LocalRow.  This is a view of an array that
    ///   you, the caller, must allocate in advance.  You may call
    ///   getNumEntriesInLocalRow() to find the number of entries in
    ///   this row on the calling process, or
    ///   getNumAllocatedEntriesInLocalRow() to get the maximum number
    ///   of entries in all rows on the calling process.  (The latter
    ///   may be useful when looping over all rows on the calling
    ///   process.)
    ///
    /// \param NumIndices [out] The number of column indices stored in
    ///   \c indices.  If \c LocalRow does not belong to (i.e., is not
    ///   in the row Map on) the calling process, then \c indices is
    ///   unchanged and \c NumIndices is returned as
    ///   <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>.
    ///
    /// This method only returns the column indices owned by the
    /// calling process.  This is a <i>local</i> operation with
    /// respect to MPI interprocess communication.  It could be the
    /// case that multiple processes store different sets of column
    /// indices in the same row.  In that case, this method only
    /// returns the column indices stored on the calling process.  Use
    /// the row Map to determine which process(es) own which rows.
    ///
    /// This method throws std::runtime_error if the output array
    /// <tt>indices</tt> is not large enough to hold the column
    /// indices in the given row.
    ///
    /// \pre <tt>isLocallyIndexed() || hasColMap()</tt>.
    ///
    /// \note Stay tuned for the following possible future changes to
    ///   this method.  First, this method will return the number of
    ///   indices in the row, and the output argument \c NumIndices
    ///   will go away.  It will be the caller's responsibility to
    ///   test whether this is less than the size of \c indices.
    ///   Second, instead of throwing an exception if \c LocalRow is
    ///   not owned by the calling process, this method will simply
    ///   return zero.  (This indicates that the calling process does
    ///   not own any column indices in that row, which is true even
    ///   if the row index is not in the row Map on any process.)
    ///   Third, this method will take a Kokkos::View instead of a
    ///   Teuchos::ArrayView.
    void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     const Teuchos::ArrayView<LocalOrdinal>& indices,
                     size_t& NumIndices) const;

    /// \brief Return a const, nonpersisting view of global indices in the given row.
    ///
    /// \param GlobalRow [in] Global index of the row.
    /// \param Indices [out] View of the column indices (as global
    ///   indices) owned by the calling process in row \c GlobalRow.
    ///   If \c GlobalRow does not belong to the calling process, then
    ///   \c Indices is set to \c null.
    ///
    /// \pre <tt> ! isLocallyIndexed () </tt>
    /// \post <tt> indices.size () == getNumEntriesInGlobalRow (GlobalRow) </tt>
    ///
    /// \note Stay tuned for the following possible future changes to
    ///   this method.  First, instead of throwing an exception if
    ///   <tt>GlobalRow</tt> is not owned by the calling process, this
    ///   method will simply return an empty view.  (This indicates
    ///   that the calling process does not own any column indices in
    ///   that row, which is true even if the row index is not in the
    ///   row Map on any process.)  Second, this method will return a
    ///   Kokkos::View instead of a Teuchos::ArrayView.
    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      Teuchos::ArrayView<const GlobalOrdinal>& Indices) const;

    /// \brief Return a const, nonpersisting view of local indices in the given row.
    ///
    /// \param LocalRow [in] Local index of the row.
    /// \param Indices [out] View of the column indices (as local
    ///   indices) owned by the calling process in row \c LocalRow.
    ///   If \c LocalRow does not belong to the calling process, then
    ///   \c indices is set to \c null.
    ///
    /// \pre <tt> ! isGloballyIndexed () </tt>
    /// \post <tt> indices.size () == getNumEntriesInLocalRow (LocalRow) </tt>
    ///
    /// \note Stay tuned for the following possible future changes to
    ///   this method.  First, instead of throwing an exception if
    ///   <tt>LocalRow</tt> is not owned by the calling process, this
    ///   method will simply return an empty view.  (This indicates
    ///   that the calling process does not own any column indices in
    ///   that row, which is true even if the row index is not in the
    ///   row Map on any process.)  Second, this method will return a
    ///   Kokkos::View instead of a Teuchos::ArrayView.
    void
    getLocalRowView (LocalOrdinal LocalRow,
                     Teuchos::ArrayView<const LocalOrdinal>& indices) const;

    //@}
    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
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
    setAllIndices (const Teuchos::ArrayRCP<size_t> & rowPointers,
                   const Teuchos::ArrayRCP<LocalOrdinal> & columnIndices);

    //! Get an ArrayRCP of the row-offsets.
    /*!  The returned buffer exists in host-memory. This method may return Teuchos::null
      if "Delete Row Pointers" was \c true on fillComplete().
    */
    Teuchos::ArrayRCP<const size_t> getNodeRowPtrs() const;

    //! Get an ArrayRCP of the packed column-indices.
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

  protected:
    typedef typename KokkosClassic::DefaultKernels<void, LocalOrdinal, Node>::SparseOps
      sparse_ops_type;
    typedef typename sparse_ops_type::template graph<LocalOrdinal, Node>::graph_type
      local_graph_type;

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
    //
    // Allocation
    //
    bool indicesAreAllocated() const;

    void allocateIndices (ELocalGlobal lg);

    template <class T>
    Teuchos::ArrayRCP<T> allocateValues1D () const;
    template <class T>
    Teuchos::ArrayRCP<Array<T> > allocateValues2D () const;

    template <ELocalGlobal lg, class T>
    RowInfo updateAllocAndValues (RowInfo rowinfo, size_t newAllocSize, Array<T>& rowVals)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( ! rowMap_->isNodeLocalElement(rowinfo.localRow) );
      TEUCHOS_TEST_FOR_EXCEPT( newAllocSize < rowinfo.allocSize );
      TEUCHOS_TEST_FOR_EXCEPT( (lg == LocalIndices && ! isLocallyIndexed()) ||
                               (lg == GlobalIndices && ! isGloballyIndexed()) );
      TEUCHOS_TEST_FOR_EXCEPT( newAllocSize == 0 );
      TEUCHOS_TEST_FOR_EXCEPT( ! indicesAreAllocated() );
#endif
      // ArrayRCP::resize automatically copies over values on reallocation.
      if (lg == LocalIndices) {
        lclInds2D_[rowinfo.localRow].resize (newAllocSize);
      }
      else { // lg == GlobalIndices
        gblInds2D_[rowinfo.localRow].resize (newAllocSize);
      }
      rowVals.resize (newAllocSize);
      nodeNumAllocated_ += (newAllocSize - rowinfo.allocSize);
      rowinfo.allocSize = newAllocSize;
      return rowinfo;
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
    size_t filterIndices (const SLocalGlobalNCViews &inds) const;

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
#endif
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
        static_cast<size_t> (fvalsend - vals.begin ());
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
      const size_t numNewInds = insertIndices (rowInfo, newInds, lg, I);
      typename Teuchos::ArrayView<const Scalar>::const_iterator newRowValsBegin =
        newRowVals.begin ();
      std::copy (newRowValsBegin, newRowValsBegin + numNewInds,
                 oldRowVals.begin () + rowInfo.numEntries);
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
    void sortRowIndices (RowInfo rowinfo);

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
    setDomainRangeMaps (const Teuchos::RCP<const map_type> &domainMap,
                        const Teuchos::RCP<const map_type> &rangeMap);

    void staticAssertions() const;
    // global consts
    void clearGlobalConstants();
    void computeGlobalConstants();

    //! \name Methods to access the graph's data
    //@{

    /// \brief Return indexing information for the given (local) row index.
    ///
    /// \param myRow [in] Local row index, as a \c size_t.
    RowInfo getRowInfo (size_t myRow) const;

    /// \brief Get a const view of the local column indices in the given row.
    ///
    /// This is one of the methods that take a RowInfo object returned
    /// by getRowInfo().  Note that it returns a view of all space for
    /// column indices in the row.  The valid entries are the first
    /// <tt>rowinfo.numEntries</tt> entries of the returned view.
    Teuchos::ArrayView<const LocalOrdinal> getLocalView (RowInfo rowinfo) const;

    /// \brief Get a nonconst view of the local column indices in the given row.
    ///
    /// This is one of the methods that take a RowInfo object returned
    /// by getRowInfo().  Note that it returns a view of all space for
    /// column indices in the row.  The valid entries are the first
    /// <tt>rowinfo.numEntries</tt> entries of the returned view.
    Teuchos::ArrayView<LocalOrdinal> getLocalViewNonConst (RowInfo rowinfo);

    /// \brief Get a const view of the global column indices in the given row.
    ///
    /// This is one of the methods that take a RowInfo object returned
    /// by getRowInfo().  Note that it returns a view of all space for
    /// column indices in the row.  The valid entries are the first
    /// <tt>rowinfo.numEntries</tt> entries of the returned view.
    Teuchos::ArrayView<const GlobalOrdinal> getGlobalView (RowInfo rowinfo) const;

    /// \brief Get a nonconst view of the global column indices in the given row.
    ///
    /// This is one of the methods that take a RowInfo object returned
    /// by getRowInfo().  Note that it returns a view of all space for
    /// column indices in the row.  The valid entries are the first
    /// <tt>rowinfo.numEntries</tt> entries of the returned view.
    Teuchos::ArrayView<GlobalOrdinal> getGlobalViewNonConst (RowInfo rowinfo);

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

    /// \brief Find the column offset corresponding to the given (global) column index.
    ///
    /// The name of this method is a bit misleading.  It does not
    /// actually find the column index.  Instead, it takes a global
    /// column index \c ind, and returns the corresponding offset
    /// into the raw array of column indices (whether that be 1-D or
    /// 2-D storage).
    size_t findGlobalIndex (RowInfo rowinfo, GlobalOrdinal ind, size_t hint = 0) const;

    //@}
    //! \name Methods to set or access the local Kokkos graph.
    //@{

    void fillLocalGraph (const Teuchos::RCP<Teuchos::ParameterList> &params);
    const Teuchos::RCP<const local_graph_type> getLocalGraph() const;
    const Teuchos::RCP<local_graph_type> getLocalGraphNonConst();

    //@}

    /// \brief Check whether the graph's state is valid.
    ///
    /// This method is useful for debugging.  It gets called often in
    /// a debug build.
    ///
    /// \note This method is <i>not</i> a collective.  It does not
    ///   communicate (between MPI processes).  Developers: Please do
    ///   not invoke communication in this method!
    void checkInternalState() const;

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

    // local data, stored in a KokkosClassic::CrsGraph. only initialized after fillComplete()
    Teuchos::RCP<local_graph_type> lclGraph_;

    // Local and Global Counts
    // nodeNumEntries_ and nodeNumAllocated_ are required to be always consistent
    // nodeMaxNumEntries_, nodeNumDiags_ and the global quantities are computed during fillComplete() and only valid when isFillComplete()
    global_size_t globalNumEntries_, globalNumDiags_, globalMaxNumRowEntries_;
    size_t          nodeNumEntries_,   nodeNumDiags_,   nodeMaxNumRowEntries_, nodeNumAllocated_;

    //! Whether the graph was allocated with static or dynamic profile.
    ProfileType pftype_;

    /// \brief The maximum number of entries to allow in each locally owned row, per row.
    ///
    /// This is an argument to some of the graph's constructors.
    /// Either this or numAllocForAllRows_ is used, but not both.
    ///
    /// If this is not set in the constructor, it is allocated
    /// temporarily, if necessary, in allocateIndices().  In that same
    /// method, it is used to allocate the row offsets array, then
    /// discarded (set to null) unconditionally.
    Teuchos::ArrayRCP<const size_t> numAllocPerRow_;

    /// \brief The maximum number of entries to allow in each locally owned row.
    ///
    /// This is an argument to some of the graph's constructors.
    /// Either this or numAllocPerRow_ is used, but not both.
    size_t numAllocForAllRows_;

    // graph indices. before allocation, all are null.
    // after allocation, except during makeIndicesLocal(), one of local or global is null.
    // we will never have 1D and 2D structures being non-null
    // this is host memory
    // 1D == StaticAllocation, 2D == DynamicAllocation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // 1D/Static structures
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! lclInds1D_ are the indices for all rows
    Teuchos::ArrayRCP<LocalOrdinal> lclInds1D_;
    //! gblInds1D_ are the indices for all rows
    Teuchos::ArrayRCP<GlobalOrdinal> gblInds1D_;
    // offset to the beg entries of each row. only used for 1D (Static) allocation.
    // i.e., indices for row R are lclInds1D_[i] for i in [b,e) where b = rowPtrs_[R] and e = rowPtrs_[R+1]
    // only the first numRowEntries_[R] of these are valid
    // both of these are null for 2D (Dynamic) allocations
    // rowPtrs_ has length N+1, while numRowEntries_ has length N
    // we may delete this to save memory on fillComplete, if "Delete Row Pointers" is specified
    Teuchos::ArrayRCP<size_t> rowPtrs_;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // 2D/Dynamic structures.
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! <tt>lclInds2D_[r]</tt> are the indices for row \c r.
    Teuchos::ArrayRCP<Teuchos::Array<LocalOrdinal> > lclInds2D_;

    //! <tt>gblInds2D_[r]</tt> are the indices for row \c r.
    Teuchos::ArrayRCP<Teuchos::Array<GlobalOrdinal> > gblInds2D_;

    /// \brief The number of local entries in each locally owned row.
    ///
    /// This is deallocated in fillComplete() if fillComplete()'s
    /// "Optimize Storage" parameter is set to \c true.
    Teuchos::ArrayRCP<size_t> numRowEntries_;

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

    /// \brief Nonlocal data given to insertGlobalValues() or
    ///   sumIntoGlobalValues().
    ///
    /// "Nonlocal" means "corresponding to rows not in the row Map on
    /// the calling process."  These data are communicated and emptied
    /// out in globalAssemble().  That method is in turn called, if
    /// necessary, in fillComplete().
    std::map<GlobalOrdinal, std::vector<GlobalOrdinal> > nonlocals_;

    /// \brief Whether it is valid to call getRowInfo().
    ///
    /// FIXME (mfh 21 Oct 2013, 28 Sep 2014) As far as I can tell,
    /// this should <i>always</i> return true.  Why do we need it?  It
    /// looks like, historically, the graph (by default) "deleted row
    /// info" at fillComplete().  This is probably why many CrsGraph
    /// methods check whether hasRowInfo() returns true before doing
    /// anything.  However, I think that nonintuitive behavior was
    /// fixed later, such that hasRowInfo() should <i>always</i>
    /// return true.  If this is actually the case, then we should
    /// dispense with this method.
    bool hasRowInfo () const;

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

#endif // defined(HAVE_TPETRACLASSIC_SERIAL) || defined(HAVE_TPETRACLASSIC_TBB) || defined(HAVE_TPETRACLASSIC_THREADPOOL) || defined(HAVE_TPETRACLASSIC_OPENMP) || defined(HAVE_TPETRACLASSIC_THRUST)

  /// \brief Nonmember function to create an empty CrsGraph given a
  ///   row Map and the max number of entries allowed locally per row.
  ///
  /// \return A dynamically allocated (DynamicProfile) graph with
  ///   specified number of nonzeros per row (defaults to zero).
  /// \relatesalso CrsGraph
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
  createCrsGraph (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map,
                   size_t maxNumEntriesPerRow = 0,
                   const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
  {
    using Teuchos::rcp;
    typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node> graph_type;
    return rcp (new graph_type (map, maxNumEntriesPerRow, DynamicProfile, params));
  }

namespace Details {

template<class LocalOrdinal, class GlobalOrdinal, class OutputNodeType, class InputNodeType>
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

// Include KokkosRefactor partial specialisation if enabled
#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#include "Tpetra_KokkosRefactor_CrsGraph_decl.hpp"
#endif

#endif
