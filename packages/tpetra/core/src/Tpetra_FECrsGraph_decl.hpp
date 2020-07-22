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

#ifndef TPETRA_FECRSGRAPH_DECL_HPP
#define TPETRA_FECRSGRAPH_DECL_HPP

/// \file Tpetra_FECrsGraph_decl.hpp
/// \brief Declaration of the Tpetra::FECrsGraph class
///
/// If you want to use Tpetra::FECrsGraph, include "Tpetra_FECrsGraph.hpp"
/// (a file which CMake generates and installs for you).  If you only
/// want the declaration of Tpetra::FECrsGraph, include this file
/// (Tpetra_FECrsGraph_decl.hpp).

#include "Tpetra_FECrsGraph_fwd.hpp"
#include "Tpetra_CrsGraph_decl.hpp"

namespace Tpetra {
  /// \class FECrsGraph
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
  /// Note on some versions of the constructors:
  ///       consturctors appear in two fashions: one that accepts only one (optional) domainMap
  ///       (let's call it V1) and one that accepts two domainMaps (let's call it V2).
  ///       In V1, the input domainMap will be used for both owned and owned+shared graphs.
  ///       In this case, you *should* provide a map that works for the owned graph (since
  ///       that's the graph you will use in mat-vec and mat-mat produts, which will be checked).
  ///       The domain map for the owned+shared graph is *usually* only needed to build the colMap
  ///       of the owned+shared graph.
  ///
  ///       Whether you need V2 or are ok with V1 depends on what your "owned+shared" maps contain.
  ///
  ///       If you partition a mesh by node (N), then your owned+shared maps should contain all the
  ///       nodes that are in the patch of one of the owned nodes, which is precisely what
  ///       the colMap of the owned+shared graph should contain.
  ///       On the other hand, if you partition the mesh by cell (C), then your owned+shared maps 
  ///       should contain all nodes belonging to one of the owned elements. Such map does *NOT*
  ///       necessarily contain all the GID of all the nodes connected to one of the owned nodes
  ///       (which is what the colMap should contain). In this case, you are not likely to have
  ///       a colMap available, and you'll rely on Tpetra to build one for you once the graph
  ///       is filled.
  ///
  ///       If you do not provide a valid owned+shared colMap, then Tpetra creates one for you,
  ///       which is used for both the owned+shared and owned graphs. When doing this, Tpetra
  ///       guarantees that the domain map of the owned+shared graph is "locally fitted" to its
  ///       colMap (see Tpetra::Map documentation for the meaning of 'locally fitted').
  ///       If you use V1, then the domain map is the owned dofs only, which means that the
  ///       non-owned GIDs will appear in the colMap in an order that is out of your control.
  ///
  ///       Now, if you partition by Cell (C), then you don't have a colMap when you create
  ///       the graph, but only a owned+shared one. If you use version V1 of the consturctor,
  ///       the Local ID of the shared dofs in the graph's colMap is likely to NOT be the same
  ///       as in your owned+shared map, since the colMap creation only guarantees that the input
  ///       domainMap (which should contian owned dofs only) is going to be locally fitted in the
  ///       colMap. All other GIDs will appear in an order that you can't control.
  ///
  ///       On the other hand, if you use V2, then you are guaranteed that the owend+shared domain
  ///       map you provide will be locally fitted to the graph colMap.
  ///
  ///       This difference is important when you do local assembly with a FECrsMatrix: if the Local IDs
  ///       of the dofs in your owned+shared maps is THE SAME as in the graph's col map, then you can do
  ///       assembly using local ids. Otherwise, you need to have access to the GIDs of the dofs,
  ///       extract the colMap from the matrix, and convert GIDs into LIDs ACCORDINGLY TO THE COL MAP.
  ///
  ///       Notice that if your mesh is partitioned by node (N) OR you provide a valid colMap
  ///       to the graph's constructor, then using V1 is likely ok for you.
  ///
  ///       For more details, see GitHub issue #7455

  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class FECrsGraph :
    public CrsGraph<LocalOrdinal, GlobalOrdinal, Node>
  {
    //! The specialization of DistObject that is this class' parent class.
    typedef DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> dist_object_type;

    template <class S, class LO, class GO, class N>
    friend class FECrsMatrix;
  public:
    //! Parent class
    typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node> crs_graph_type;

    //! This class' first template parameter; the type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! This class' second template parameter; the type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! This class' Kokkos Node type.
    typedef Node node_type;

    //! The Kokkos device type.
    typedef typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::device_type device_type;
    //! The Kokkos execution space.
    typedef typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::execution_space execution_space;

    //! The type of the part of the sparse graph on each MPI process.
    typedef  typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type local_graph_type;

    //! The Map specialization used by this class.
    using map_type = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    //! The Import specialization used by this class.
    using import_type = ::Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
    //! The Export specialization used by this class.
    using export_type = ::Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>;

    //! @name Constructor/Destructor Methods
    //@{

    /// \brief Constructorfor globally-indexed assembly specifying a single upper bound for the
    ///   number of entries in all rows on the calling process.
    ///
    /// \param ownedRowMap [in] Distribution of rows of the owned graph.
    ///
    /// \param ownedPlusSharedRowMap [in] ownedMap plus the list of shared rows to which off-processor insertion is allowed
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  You cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param ownedPlusSharedToOwnedimporter [in] Optional importer between the ownedMap and ownedPlusSharedMap
    ///   This will be calculated by FECrsGraph if it is not provided
    ///
    /// \param domainMap [in] Optional domain map for both owned and owned+shared graphs. If this is null,
    ////  then ownedMap will be used as domain map for both graphs.
    ///
    /// \param ownedRangeMap [in] Optional range map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as range map for the owned graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not 
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    FECrsGraph(const Teuchos::RCP<const map_type> & ownedRowMap,
               const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
               const size_t maxNumEntriesPerRow,
               const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter = Teuchos::null,
               const Teuchos::RCP<const map_type> & domainMap = Teuchos::null,
               const Teuchos::RCP<const map_type> & ownedRangeMap = Teuchos::null,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructorfor globally-indexed assembly specifying a single upper bound for the
    ///   number of entries in all rows on the calling process.
    ///
    /// \param ownedRowMap [in] Distribution of rows of the owned graph.
    ///
    /// \param ownedPlusSharedRowMap [in] ownedMap plus the list of shared rows to which off-processor insertion is allowed
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  You cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param ownedPlusSharedDomainMap [in] Domain map for the owned+shared graph. If this is null,
    ///   then ownedPlusSharedRowMap will be used as domain map for the owned+shared graph.
    ///
    /// \param ownedPlusSharedToOwnedimporter [in] Optional importer between the ownedMap and ownedPlusSharedMap
    ///   This will be calculated by FECrsGraph if it is not provided
    ///
    /// \param ownedDomainMap [in] Optional domain map for the owned graph. If this is null,
    ////  then ownedMap will be used as domain map for the owned graph.
    ///
    /// \param ownedRangeMap [in] Optional range map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as range map for the owned graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
                const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
                const size_t maxNumEntriesPerRow,
                const Teuchos::RCP<const map_type> & ownedPlusSharedDomainMap,
                const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter = Teuchos::null,
                const Teuchos::RCP<const map_type> & ownedDomainMap = Teuchos::null,
                const Teuchos::RCP<const map_type> & ownedRangeMap = Teuchos::null,
                const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor for globally-indexed assembly specifying a (possibly different) upper
    ///   bound for the number of entries in each row.
    ///
    /// \param ownedRowMap [in] Distribution of rows of the owned graph.
    ///
    /// \param ownedPlusSharedRowMap [in] ownedMap plus the list of shared rows to which off-processor insertion is allowed
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  You cannot exceed the allocated
    ///   number of entries for any row.
    ///
    /// \param ownedPlusSharedToOwnedimporter [in] Optional importer between the ownedMap and ownedPlusSharedMap
    ///   This will be calculated by FECrsGraph if it is not provided
    ///
    /// \param domainMap [in] Optional domain map for both owned and owned+shared graphs. If this is null,
    ////  then ownedMap will be used as domain map for both graphs.
    ///
    /// \param ownedRangeMap [in] Optional range map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as range map for the owned graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
                const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
                const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
                const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter = Teuchos::null,
                const Teuchos::RCP<const map_type> & domainMap = Teuchos::null,
                const Teuchos::RCP<const map_type> & ownedRangeMap = Teuchos::null,
                const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor for globally-indexed assembly specifying a (possibly different) upper
    ///   bound for the number of entries in each row.
    ///
    /// \param ownedRowMap [in] Distribution of rows of the owned graph.
    ///
    /// \param ownedPlusSharedRowMap [in] ownedMap plus the list of shared rows to which off-processor insertion is allowed
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  You cannot exceed the allocated
    ///   number of entries for any row.
    ///
    /// \param ownedPlusSharedDomainMap [in] Domain map for the owned+shared graph. If this is null,
    ///   then ownedPlusSharedRowMap will be used in the call to endFill()
    ///
    /// \param ownedPlusSharedToOwnedimporter [in] Optional importer between the ownedMap and ownedPlusSharedMap
    ///   This will be calculated by FECrsGraph if it is not provided
    ///
    /// \param ownedDomainMap [in] Optional domain map for the owned graph. If this is not provided,
    ////  then ownedMap will be used as domain map for the owned graph.
    ///
    /// \param ownedRangeMap [in] Optional range map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as range map for the owned graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
                const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
                const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
                const Teuchos::RCP<const map_type> & ownedPlusSharedDomainMap,
                const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter = Teuchos::null,
                const Teuchos::RCP<const map_type> & ownedDomainMap = Teuchos::null,
                const Teuchos::RCP<const map_type> & ownedRangeMap = Teuchos::null,
                const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


    /// \brief Constructor for locally-indexed assembly specifying a single upper bound for the
    ///   number of entries in all rows on the calling process.
    ///
    /// \param ownedRowMap [in] Distribution of rows of the owned graph.
    ///
    /// \param ownedPlusSharedRowMap [in] ownedMap plus the list of shared rows to which off-processor insertion is allowed
    ///
    /// \param ownedPlusSharedColMap [in] list of owned and shared columns into which assertion is allowed.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  You cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param ownedPlusSharedToOwnedimporter [in] Optional importer between the ownedMap and ownedPlusSharedMap
    ///   This will be calculated by FECrsGraph if it is not provided
    ///
    /// \param domainMap [in] Optional domain map for the both owned and owned+shared graphs.  If this is not provided,
    ///   then ownedMap will be used as domain map for both graphs.
    ///
    /// \param ownedRangeMap [in] Optional range map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as range map for the owned graph
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    FECrsGraph(const Teuchos::RCP<const map_type> & ownedRowMap,
               const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
               const Teuchos::RCP<const map_type> & ownedPlusSharedColMap,
               const size_t maxNumEntriesPerRow,
               const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter = Teuchos::null,
               const Teuchos::RCP<const map_type> & domainMap = Teuchos::null,
               const Teuchos::RCP<const map_type> & ownedRangeMap = Teuchos::null,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor for locally-indexed assembly specifying a single upper bound for the
    ///   number of entries in all rows on the calling process.
    ///
    /// \param ownedRowMap [in] Distribution of rows of the owned graph.
    ///
    /// \param ownedPlusSharedRowMap [in] ownedMap plus the list of shared rows to which off-processor insertion is allowed
    ///
    /// \param ownedPlusSharedColMap [in] list of owned and shared columns into which assertion is allowed.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of graph
    ///   entries per row.  You cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param ownedPlusSharedDomainMap [in] Domain map for the owned+shared graph. If this is null, then ownedPlusSharedRowMap
    ///   will be used as domain map for the owned+shared graph.
    /// \param ownedPlusSharedToOwnedimporter [in] Optional importer between the ownedMap and ownedPlusSharedMap
    ///   This will be calculated by FECrsGraph if it is not provided
    ///
    /// \param ownedDomainMap [in] Optional domain map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as domain map for the owned graph.
    ///
    /// \param ownedRangeMap [in] Optional range map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as range map for the owned graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    FECrsGraph(const Teuchos::RCP<const map_type> & ownedRowMap,
               const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
               const Teuchos::RCP<const map_type> & ownedPlusSharedColMap,
               const size_t maxNumEntriesPerRow,
               const Teuchos::RCP<const map_type> & ownedPlusSharedDomainMap,
               const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter = Teuchos::null,
               const Teuchos::RCP<const map_type> & ownedDomainMap = Teuchos::null,
               const Teuchos::RCP<const map_type> & ownedRangeMap = Teuchos::null,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


    /// \brief Constructor for locally-indexed assembly specifying a (possibly different) upper
    ///   bound for the number of entries in each row.
    ///
    /// \param ownedRowMap [in] Distribution of rows of the owned graph.
    ///
    /// \param ownedPlusSharedRowMap [in] ownedMap plus the list of shared rows to which off-processor insertion is allowed
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  You cannot exceed the allocated
    ///   number of entries for any row.
    ///
    /// \param ownedPlusSharedToOwnedimporter [in] Optional importer between the ownedMap and ownedPlusSharedMap
    ///   This will be calculated by FECrsGraph if it is not provided
    ///
    /// \param domainMap [in] Optional domain map for both the owned and owned+shared graphs.  If this is not provided,
    ///   then ownedMap will be used as domain map for both graphs.
    ///
    /// \param ownedRangeMap [in] Optional range map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as range map for the owned graph
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
                const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
                const Teuchos::RCP<const map_type> & ownedPlusSharedColMap,
                const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
                const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter = Teuchos::null,
                const Teuchos::RCP<const map_type> & domainMap = Teuchos::null,
                const Teuchos::RCP<const map_type> & ownedRangeMap = Teuchos::null,
                const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor for locally-indexed assembly specifying a (possibly different) upper
    ///   bound for the number of entries in each row.
    ///
    /// \param ownedRowMap [in] Distribution of rows of the owned graph.
    ///
    /// \param ownedPlusSharedRowMap [in] ownedMap plus the list of shared rows to which off-processor insertion is allowed
    ///
    /// \param numEntPerRow [in] Maximum number of graph entries to
    ///   allocate for each row.  You cannot exceed the allocated
    ///   number of entries for any row.
    ///
    /// \param ownedPlusSharedDomainMap [in] Domain map for the owned+shared graph. If this is null, then ownedPlusSharedRowMap
    /// \param ownedPlusSharedToOwnedimporter [in] Optional importer between the ownedMap and ownedPlusSharedMap
    ///   This will be calculated by FECrsGraph if it is not provided
    ///
    /// \param ownedDomainMap [in] Optional domain map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as domain map for the owned graph.
    ///
    /// \param ownedRangeMap [in] Optional range map for the owned graph.  If this is not provided, then ownedMap
    ///   will be used as range map for the owned graph.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
                const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,
                const Teuchos::RCP<const map_type> & ownedPlusSharedColMap,
                const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
                const Teuchos::RCP<const map_type> & ownedPlusSharedDomainMap,
                const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter = Teuchos::null,
                const Teuchos::RCP<const map_type> & ownedDomainMap = Teuchos::null,
                const Teuchos::RCP<const map_type> & ownedRangeMap = Teuchos::null,
                const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Copy constructor (forbidden).
    FECrsGraph (const FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>&) = delete;

    //! Move constructor (forbidden).
    FECrsGraph (FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>&&) = delete;

    //! Copy assignment (forbidden).
    FECrsGraph&
    operator= (const FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>&) = delete;

    //! Move assignment (forbidden).
    FECrsGraph&
    operator= (FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>&&) = delete;

    /// \brief Destructor (virtual for memory safety of derived classes).
    ///
    /// \note To Tpetra developers: See the C++ Core Guidelines C.21
    ///   ("If you define or <tt>=delete</tt> any default operation,
    ///   define or <tt>=delete</tt> them all"), in particular the
    ///   AbstractBase example, for why this destructor declaration
    ///   implies that we need the above four <tt>=delete</tt>
    ///   declarations for copy construction, move construction, copy
    ///   assignment, and move assignment.
    virtual ~FECrsGraph () = default;

    //@}
    //! @name Collective methods for changing the graph's global state
    //@{

    //! Migrates data to the owned mode
    void endFill();

    //! Activates the owned+shared mode for assembly.  This can only be called once.
    void beginFill();

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
    /// \param ownedDomainMap [in] The graph's domain Map.  MUST be one to
    ///   one!
    /// \param ownedRangeMap [in] The graph's range Map.  MUST be one to
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
    fillComplete (const Teuchos::RCP<const map_type>& ownedDomainMap,
                  const Teuchos::RCP<const map_type>& ownedRangeMap,
                  const Teuchos::RCP<Teuchos::ParameterList>& /* params */ = Teuchos::null) {
      ownedDomainMap_ = ownedDomainMap;
      ownedRangeMap_ = ownedRangeMap;
      endFill();
    }

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
    fillComplete (const Teuchos::RCP<Teuchos::ParameterList>& /* params */ = Teuchos::null) {endFill();}


  private:
    /// \brief Migrate data from the owned+shared to the owned graph
    /// Since this is non-unique -> unique, we need a combine mode.
    /// Precondition: Must be FE_ACTIVE_OWNED_PLUS_SHARED mode
    void doOwnedPlusSharedToOwned(const CombineMode CM=Tpetra::ADD);

    /// \brief Migrate data from the owned to the owned+shared graph
    /// Precondition: Must be FE_ACTIVE_OWNED mode
    void doOwnedToOwnedPlusShared(const CombineMode CM=Tpetra::ADD);

  public:
    //! Switches which CrsGraph is active (without migrating data)
    void switchActiveCrsGraph();
    //@}

  private:

    // Common core guts of the constructor (the colMap argument is Teuchos::null if we're globally-indexed)
    void setup(const Teuchos::RCP<const map_type>  & ownedRowMap, const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,const Teuchos::RCP<const map_type> & ownedPlusSharedColMap, const Teuchos::RCP<Teuchos::ParameterList>& params);

    // Template to avoid long type names...lazy.
    // template<typename ViewType>
    // Teuchos::RCP<const map_type> makeOwnedColMap (ViewType ownedGraphIndices);

    // Enum for activity
    enum FEWhichActive
    {
      FE_ACTIVE_OWNED,
      FE_ACTIVE_OWNED_PLUS_SHARED
    };


    // This is whichever graph isn't currently active
    Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > inactiveCrsGraph_;

    // This is in RCP to make shallow copies of the FECrsGraph work correctly
    Teuchos::RCP<FEWhichActive> activeCrsGraph_;

    // The importer between the rowmaps of the two graphs
    Teuchos::RCP<const import_type> ownedRowsImporter_;

    // The domainMap to use in endFill() for the owned graph
    Teuchos::RCP<const map_type> ownedDomainMap_;

    // The rangeMap to use in endFill() for the owned graph
    Teuchos::RCP<const map_type> ownedRangeMap_;
  }; // class FECrsGraph


} // namespace Tpetra

#endif // TPETRA_FECRSGRAPH_DECL_HPP
