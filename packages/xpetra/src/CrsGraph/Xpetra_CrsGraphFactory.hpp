// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_CRSGRAPHFACTORY_HPP
#define XPETRA_CRSGRAPHFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_CrsGraph.hpp"

#include "Xpetra_TpetraCrsGraph.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsGraph.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class CrsGraphFactory {
  private:
    //! Private constructor. This is a static class.
    CrsGraphFactory() {}

  public:
    //! Constructor for empty graph (intended use is an import/export target - can't insert entries directly)
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build (const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(rowMap->lib() == UseEpetra, std::logic_error,
          "Can't create Xpetra::EpetraCrsMatrix with these scalar/LO/GO types");
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, 0) );

      XPETRA_FACTORY_END;
    }

    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t maxNumEntriesPerRow) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

      if (map->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> (map, maxNumEntriesPerRow) );

      XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

    //! Constructor specifying column Map and number of entries per row
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          size_t maxNumEntriesPerRow,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }


    //! Constructor specifying column Map and number of entries in each row.
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }


    //! Constructor using fused import
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const RCP<const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > >& sourceGraph,
          const Import< LocalOrdinal, GlobalOrdinal, Node > & importer,
          const RCP<const Map< LocalOrdinal, GlobalOrdinal, Node >>& domainMap = Teuchos::null,
          const RCP<const Map< LocalOrdinal, GlobalOrdinal, Node > >& rangeMap = Teuchos::null,
          const RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      if (sourceGraph->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(sourceGraph, importer, domainMap, rangeMap, params) );

      XPETRA_FACTORY_ERROR_IF_EPETRA(sourceGraph()->getRowMap()->lib());
      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }



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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::row_map_type& rowPointers,
          const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::entries_type::non_const_type& columnIndices,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap,
                                                                          colMap,
                                                                          rowPointers,
                                                                          columnIndices,
                                                                          plist));

      XPETRA_FACTORY_END;
    }

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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rowMap,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& colMap,
          const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type& lclGraph,
          const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap,
                                                                          colMap,
                                                                          lclGraph,
                                                                          params));

      XPETRA_FACTORY_END;
    }

    /// \brief Constructor specifying column, domain and range maps, and a
    ///   local (sorted) graph, which the resulting CrsGraph views.
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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
          Build(const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type& lclGraph,
                const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rowMap,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& colMap,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap = Teuchos::null,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap = Teuchos::null,
          const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(lclGraph,
                                                                          rowMap,
                                                                          colMap,
                                                                          domainMap,
                                                                          rangeMap,
                                                                          params));

      XPETRA_FACTORY_END;
    }

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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          const Teuchos::ArrayRCP<size_t>& rowPointers,
          const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap,
                                                                          colMap,
                                                                          rowPointers,
                                                                          columnIndices,
                                                                          plist));

      XPETRA_FACTORY_END;
    }

  };

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES))

  template <>
  class CrsGraphFactory<int, int, EpetraNode> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    CrsGraphFactory() {}

  public:
    //! Constructor for empty graph (intended use is an import/export target - can't insert entries directly)
    static RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build (const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap)
    {
      XPETRA_MONITOR("CrsMatrixFactory::Build");
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, 0) );
#ifdef HAVE_XPETRA_EPETRA
      if(rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<int,Node>(rowMap));
#endif
      XPETRA_FACTORY_END;
    }

    static RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t maxNumEntriesPerRow) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

      if (map->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> (map, maxNumEntriesPerRow) );

      if (map->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<int, Node>(map, maxNumEntriesPerRow) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<int, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

    //! Constructor specifying column Map and number of entries per row
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          size_t maxNumEntriesPerRow,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );
      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<int,Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }


    //! Constructor using fused import
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const RCP<const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > >& sourceGraph,
          const Import< LocalOrdinal, GlobalOrdinal, Node > & importer,
          const RCP<const Map< LocalOrdinal, GlobalOrdinal, Node >>& domainMap = Teuchos::null,
          const RCP<const Map< LocalOrdinal, GlobalOrdinal, Node > >& rangeMap = Teuchos::null,
          const RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      if (sourceGraph->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(sourceGraph, importer, domainMap, rangeMap, params) );
      if (sourceGraph->getRowMap()->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<int, Node>(sourceGraph, importer, domainMap, rangeMap, params) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::row_map_type& rowPointers,
          const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::entries_type::non_const_type& columnIndices,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap,
                                                                          colMap,
                                                                          rowPointers,
                                                                          columnIndices,
                                                                          plist));

      XPETRA_FACTORY_END;
    }

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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rowMap,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& colMap,
          const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type& lclGraph,
          const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap,
                                                                          colMap,
                                                                          lclGraph,
                                                                          params));

      XPETRA_FACTORY_END;
    }

    /// \brief Constructor specifying column, domain and range maps, and a
    ///   local (sorted) graph, which the resulting CrsGraph views.
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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
          Build(const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type& lclGraph,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rowMap,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& colMap,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap = Teuchos::null,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap = Teuchos::null,
          const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(lclGraph,
                                                                          rowMap,
                                                                          colMap,
                                                                          domainMap,
                                                                          rangeMap,
                                                                          params));

      XPETRA_FACTORY_END;
    }

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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          const Teuchos::ArrayRCP<size_t>& rowPointers,
          const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap,
                                                                          colMap,
                                                                          rowPointers,
                                                                          columnIndices,
                                                                          plist));

      XPETRA_FACTORY_END;
    }

  };
#endif

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES))

  template <>
  class CrsGraphFactory<int, long long, EpetraNode> {

    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    CrsGraphFactory() {}

  public:
    //! Constructor for empty graph (intended use is an import/export target - can't insert entries directly)
    static RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build (const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap)
    {
      XPETRA_MONITOR("CrsMatrixFactory::Build");
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, 0) );
#ifdef HAVE_XPETRA_EPETRA
      if(rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<long long,Node>(rowMap));
#endif
      XPETRA_FACTORY_END;
    }

    static RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t maxNumEntriesPerRow) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

      if (map->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> (map, maxNumEntriesPerRow) );

      if (map->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<long long, Node>(map, maxNumEntriesPerRow) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<long long, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

    //! Constructor specifying column Map and number of entries per row
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          size_t maxNumEntriesPerRow,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );
      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsGraphT<long long,Node>(rowMap, colMap, maxNumEntriesPerRow, plist) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }

    //! Constructor using fused import
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const RCP<const CrsGraph< LocalOrdinal, GlobalOrdinal, Node > >& sourceGraph,
          const Import< LocalOrdinal, GlobalOrdinal, Node > & importer,
          const RCP<const Map< LocalOrdinal, GlobalOrdinal, Node >>& domainMap = Teuchos::null,
          const RCP<const Map< LocalOrdinal, GlobalOrdinal, Node > >& rangeMap = Teuchos::null,
          const RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      if (sourceGraph->getRowMap()->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(sourceGraph, importer, domainMap, rangeMap, params) );
      if (sourceGraph->getRowMap()->lib() == UseTpetra)
        return rcp( new EpetraCrsGraphT<long long,Node><LocalOrdinal, GlobalOrdinal, Node>(sourceGraph, importer, domainMap, rangeMap, params) );

      XPETRA_FACTORY_END;
      TEUCHOS_UNREACHABLE_RETURN(null);
    }


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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::row_map_type& rowPointers,
          const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type::entries_type::non_const_type& columnIndices,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap,
                                                                          colMap,
                                                                          rowPointers,
                                                                          columnIndices,
                                                                          plist));

      XPETRA_FACTORY_END;
    }

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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rowMap,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& colMap,
          const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type& lclGraph,
          const Teuchos::RCP<Teuchos::ParameterList>& params) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap,
                                                                          colMap,
                                                                          lclGraph,
                                                                          params));

      XPETRA_FACTORY_END;
    }

    /// \brief Constructor specifying column, domain and range maps, and a
    ///   local (sorted) graph, which the resulting CrsGraph views.
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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type& lclGraph,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rowMap,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& colMap,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap = Teuchos::null,
          const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap = Teuchos::null,
          const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(lclGraph,
                                                                          rowMap,
                                                                          colMap,
                                                                          domainMap,
                                                                          rangeMap,
                                                                          params));

      XPETRA_FACTORY_END;
    }

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
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap,
          const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap,
          const Teuchos::ArrayRCP<size_t>& rowPointers,
          const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
          const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null) {
      XPETRA_MONITOR("CrsMatrixFactory::Build");

      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>(rowMap,
                                                                          colMap,
                                                                          rowPointers,
                                                                          columnIndices,
                                                                          plist));

      XPETRA_FACTORY_END;
    }

  };
#endif
}

#define XPETRA_CRSGRAPHFACTORY_SHORT
#endif
