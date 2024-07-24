// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CRSGRAPHTRANSPOSER_DEF_HPP
#define TPETRA_CRSGRAPHTRANSPOSER_DEF_HPP

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_shortSort.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "KokkosSparse_Utils.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosSparse_spadd.hpp"

namespace Tpetra {

  template<typename GO,
           typename LocalIndicesType,
           typename GlobalIndicesType,
           typename ColMapType>
  struct ConvertLocalToGlobalFunctor
  {
    ConvertLocalToGlobalFunctor(
                                const LocalIndicesType& colindsOrig_,
                                const GlobalIndicesType& colindsConverted_,
                                const ColMapType& colmap_) :
      colindsOrig (colindsOrig_),
      colindsConverted (colindsConverted_),
      colmap (colmap_)
    {}
    KOKKOS_INLINE_FUNCTION void
    operator() (const GO i) const
    {
      colindsConverted(i) = colmap.getGlobalElement(colindsOrig(i));
    }
    LocalIndicesType colindsOrig;
    GlobalIndicesType colindsConverted;
    ColMapType colmap;
  };

  template<class LO, class GO, class LOView, class GOView, class LocalMap>
  struct ConvertGlobalToLocalFunctor
  {
    ConvertGlobalToLocalFunctor(LOView& lids_, const GOView& gids_, const LocalMap localColMap_)
      : lids(lids_), gids(gids_), localColMap(localColMap_)
    {}

    KOKKOS_FUNCTION void operator() (const GO i) const
    {
      lids(i) = localColMap.getLocalElement(gids(i));
    }

    LOView lids;
    const GOView gids;
    const LocalMap localColMap;
  };


  template <typename size_type, typename ordinal_type,
            typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
            typename AcolindsT, typename BcolindsT, typename CcolindsT>
  struct SortedNumericIndicesOnlyFunctor {

    SortedNumericIndicesOnlyFunctor(const ArowptrsT& Arowptrs_,
                                    const BrowptrsT& Browptrs_,
                                    const CrowptrsT& Crowptrs_,
                                    const AcolindsT& Acolinds_,
                                    const BcolindsT& Bcolinds_,
                                    const CcolindsT& Ccolinds_)
      : Arowptrs(Arowptrs_),
        Browptrs(Browptrs_),
        Crowptrs(Crowptrs_),
        Acolinds(Acolinds_),
        Bcolinds(Bcolinds_),
        Ccolinds(Ccolinds_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const
    {
      const ordinal_type ORDINAL_MAX = Kokkos::ArithTraits<ordinal_type>::max();

      // count the union of nonzeros in Arow and Brow
      size_type ai         = 0;
      size_type bi         = 0;
      size_type Arowstart  = Arowptrs(i);
      size_type Arowlen    = Arowptrs(i + 1) - Arowstart;
      size_type Browstart  = Browptrs(i);
      size_type Browlen    = Browptrs(i + 1) - Browstart;
      ordinal_type Acol = (Arowlen == 0) ? ORDINAL_MAX : Acolinds(Arowstart);
      ordinal_type Bcol = (Browlen == 0) ? ORDINAL_MAX : Bcolinds(Browstart);
      size_type Coffset = Crowptrs(i);
      while (Acol != ORDINAL_MAX || Bcol != ORDINAL_MAX)
        {
          ordinal_type Ccol = (Acol < Bcol) ? Acol : Bcol;
          while(Acol == Ccol)
            {
              ai++;
              if(ai == Arowlen)
                Acol = ORDINAL_MAX;
              else
                Acol = Acolinds(Arowstart + ai);
            }
          while(Bcol == Ccol)
            {
              bi++;
              if(bi == Browlen)
                Bcol = ORDINAL_MAX;
              else
                Bcol = Bcolinds(Browstart + bi);
            }
          Ccolinds(Coffset) = Ccol;
          Coffset++;
        }
    }

    const ArowptrsT Arowptrs;
    const BrowptrsT Browptrs;
    const CrowptrsT Crowptrs;
    const AcolindsT Acolinds;
    const BcolindsT Bcolinds;
    CcolindsT Ccolinds;
  };

  template <typename size_type, typename ordinal_type,
            typename ArowptrsT, typename BrowptrsT, typename CrowptrsT,
            typename AcolindsT, typename BcolindsT, typename CcolindsT>
  struct UnsortedNumericIndicesOnlyFunctor {

    UnsortedNumericIndicesOnlyFunctor(
                                      const ArowptrsT Arowptrs_, const BrowptrsT Browptrs_, const CrowptrsT Crowptrs_,
                                      const AcolindsT Acolinds_, const BcolindsT Bcolinds_, CcolindsT Ccolinds_,
                                      const CcolindsT Apos_, const CcolindsT Bpos_)
      : Arowptrs(Arowptrs_),
        Browptrs(Browptrs_),
        Crowptrs(Crowptrs_),
        Acolinds(Acolinds_),
        Bcolinds(Bcolinds_),
        Ccolinds(Ccolinds_),
        Apos(Apos_),
        Bpos(Bpos_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const ordinal_type i) const {
      size_type CrowStart = Crowptrs(i);
      size_type ArowStart = Arowptrs(i);
      size_type ArowEnd   = Arowptrs(i + 1);
      size_type BrowStart = Browptrs(i);
      size_type BrowEnd   = Browptrs(i + 1);
      // add in A entries, while setting C colinds
      for (size_type j = ArowStart; j < ArowEnd; j++) {
        Ccolinds(CrowStart + Apos(j)) = Acolinds(j);
      }
      // add in B entries, while setting C colinds
      for (size_type j = BrowStart; j < BrowEnd; j++) {
        Ccolinds(CrowStart + Bpos(j)) = Bcolinds(j);
      }
    }
    const ArowptrsT Arowptrs;
    const BrowptrsT Browptrs;
    const CrowptrsT Crowptrs;
    const AcolindsT Acolinds;
    const BcolindsT Bcolinds;
    CcolindsT Ccolinds;
    const CcolindsT Apos;
    const CcolindsT Bpos;
  };


  template<class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
  CrsGraphTransposer<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraphTransposer (const Teuchos::RCP<const crs_graph_type>& origGraph,
                      const std::string& label)
    : origGraph_ (origGraph), label_ (label)
  {}

  template<class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
  Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
  CrsGraphTransposer<LocalOrdinal, GlobalOrdinal, Node>::
  symmetrize (const Teuchos::RCP<Teuchos::ParameterList> &params)
  {
    using Teuchos::RCP;
    using device_type = typename Node::device_type;
    using execution_space = typename device_type::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, size_t>;
    using local_graph_device_type = typename crs_graph_type::local_graph_device_type;
    using impl_scalar_type = ::Tpetra::Details::DefaultTypes::scalar_type;
    using row_ptrs_array = typename local_graph_device_type::row_map_type::non_const_type ;
    using col_inds_array = typename local_graph_device_type::entries_type::non_const_type;
    using local_map_type = typename map_type::local_map_type;
    using global_col_inds_array = typename Kokkos::View<GlobalOrdinal*, device_type>;

    auto graph = origGraph_;
    auto domain_map = graph->getDomainMap();
    auto range_map = graph->getRangeMap();
    auto row_map = graph->getRowMap();
    auto col_map = graph->getColMap();
    RCP<const map_type> col_map_sym;
    RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > importer;

    TEUCHOS_ASSERT(domain_map->isSameAs(*range_map));
    TEUCHOS_ASSERT(domain_map->isSameAs(*row_map));

    // Do the transpose
    RCP<crs_graph_type> graphT = createTranspose (params);

    auto col_map_T = graphT->getColMap();
    TEUCHOS_ASSERT(!col_map_T.is_null());
    TEUCHOS_ASSERT(domain_map->isSameAs(*graphT->getDomainMap()));

    bool graphSorted  = graph->isSorted();
    bool graphTSorted = graphT->isSorted();
    bool sorted = graphSorted && graphTSorted;
    bool matchingColMaps = col_map->isSameAs(*col_map_T);

    auto lclGraph  = graph->getLocalGraphDevice();
    auto lclGraphT = graphT->getLocalGraphDevice();

    using KKH_LO = KokkosKernels::Experimental::KokkosKernelsHandle<size_t, LocalOrdinal, impl_scalar_type,
                                                                    typename Node::execution_space, typename Node::memory_space, typename Node::memory_space>;
    using KKH_GO = KokkosKernels::Experimental::KokkosKernelsHandle<size_t, GlobalOrdinal, impl_scalar_type,
                                                                    typename Node::execution_space, typename Node::memory_space, typename Node::memory_space>;

    auto rowptrs  = lclGraph.row_map;
    auto rowptrsT = lclGraphT.row_map;
    auto colinds  = lclGraph.entries;
    auto colindsT = lclGraphT.entries;

    auto nrows = rowptrs.extent(0) - 1;
    auto rowptrsSym = row_ptrs_array(Kokkos::ViewAllocateWithoutInitializing("row ptrs sym"), nrows + 1);

    col_inds_array colindsSym;

    if(!matchingColMaps) {
      // convert indices of local graph to GlobalOrdinal
      auto lclColmap  = col_map->getLocalMap();
      global_col_inds_array colindsConverted(Kokkos::ViewAllocateWithoutInitializing("colinds (converted)"), colinds.extent(0));
      ConvertLocalToGlobalFunctor<GlobalOrdinal, col_inds_array, global_col_inds_array, local_map_type> convert(colinds, colindsConverted, lclColmap);
      Kokkos::parallel_for("colInds (converted)", range_type(0, colinds.extent(0)), convert);

      // convert indices of local graphT to GlobalOrdinal
      auto lclColmapT = col_map_T->getLocalMap();
      global_col_inds_array colindsTConverted(Kokkos::ViewAllocateWithoutInitializing("colindsT (converted)"), colindsT.extent(0));
      ConvertLocalToGlobalFunctor<GlobalOrdinal, col_inds_array, global_col_inds_array, local_map_type> convertT(colindsT, colindsTConverted, lclColmapT);
      Kokkos::parallel_for("colIndsT (converted)", range_type(0, colindsT.extent(0)), convertT);

      // sum graph and graphT in GlobalOrdinal
      KKH_GO handle;
      handle.create_spadd_handle(false);
      auto addHandle = handle.get_spadd_handle();

      global_col_inds_array globalColindsSym;

      KokkosSparse::Experimental::spadd_symbolic
        (&handle,
#if KOKKOSKERNELS_VERSION >= 40299
         nrows, graph->getGlobalNumCols(),
#endif
         rowptrs, colindsConverted, rowptrsT, colindsTConverted, rowptrsSym);
      globalColindsSym = global_col_inds_array(Kokkos::ViewAllocateWithoutInitializing("global colinds sym"), addHandle->get_c_nnz());

      UnsortedNumericIndicesOnlyFunctor<
        size_t, GlobalOrdinal,
        typename row_ptrs_array::const_type, typename row_ptrs_array::const_type, row_ptrs_array,
        typename global_col_inds_array::const_type, typename global_col_inds_array::const_type, global_col_inds_array>
        unsortedNumeric(rowptrs, rowptrsT, rowptrsSym,
                        colindsConverted, colindsTConverted, globalColindsSym,
                        addHandle->get_a_pos(), addHandle->get_b_pos());
      Kokkos::parallel_for("KokkosSparse::SpAdd:Numeric::InputNotSorted",
                           range_type(0, nrows), unsortedNumeric);

      // build column map for graphSym
      Tpetra::Details::makeColMap<LocalOrdinal, GlobalOrdinal, Node>
        (col_map_sym, domain_map, globalColindsSym);

      // convert indices of local graphSym to LocalOrdinal
      auto lclColmapSym = col_map_sym->getLocalMap();
      colindsSym = col_inds_array("colindsSym", globalColindsSym.extent(0));
      ConvertGlobalToLocalFunctor<LocalOrdinal, GlobalOrdinal, col_inds_array, global_col_inds_array, typename map_type::local_map_type> convertSym(colindsSym, globalColindsSym, lclColmapSym);
      Kokkos::parallel_for(range_type(0, globalColindsSym.extent(0)), convertSym);

    } else {

      // sum graph and graphT in LocalOrdinal
      KKH_LO handle;
      handle.create_spadd_handle(sorted);
      auto addHandle = handle.get_spadd_handle();

      KokkosSparse::Experimental::spadd_symbolic
        (&handle,
#if KOKKOSKERNELS_VERSION >= 40299
         nrows, graph->getGlobalNumCols(),
#endif
         rowptrs, colinds, rowptrsT, colindsT, rowptrsSym);
      colindsSym = col_inds_array(Kokkos::ViewAllocateWithoutInitializing("C colinds"), addHandle->get_c_nnz());

      if (sorted) {
        SortedNumericIndicesOnlyFunctor<
          size_t, LocalOrdinal,
          typename row_ptrs_array::const_type, typename row_ptrs_array::const_type, row_ptrs_array,
          typename col_inds_array::const_type, typename col_inds_array::const_type, col_inds_array>
          sortedNumeric(rowptrs, rowptrsT, rowptrsSym,
                        colinds, colindsT, colindsSym);
        Kokkos::parallel_for("KokkosSparse::SpAdd:Numeric::InputSorted",
                             range_type(0, nrows), sortedNumeric);

      } else {
        UnsortedNumericIndicesOnlyFunctor<
          size_t, LocalOrdinal,
          typename row_ptrs_array::const_type, typename row_ptrs_array::const_type, row_ptrs_array,
          typename col_inds_array::const_type, typename col_inds_array::const_type, col_inds_array>
          unsortedNumeric(rowptrs, rowptrsT, rowptrsSym,
                          colinds, colindsT, colindsSym,
                          addHandle->get_a_pos(), addHandle->get_b_pos());
        Kokkos::parallel_for("KokkosSparse::SpAdd:Numeric::InputNotSorted",
                             range_type(0, nrows), unsortedNumeric);
      }

      // column map for graphSym is graph's column map
      col_map_sym = col_map;
      importer = graph->getImporter();
    }

    bool sort = true;
    if (sort)
      KokkosSparse::sort_crs_graph<execution_space, row_ptrs_array, col_inds_array>(rowptrsSym, colindsSym);

    local_graph_device_type lclGraphSym = local_graph_device_type(colindsSym, rowptrsSym);

    RCP<Teuchos::ParameterList> graphParams = Teuchos::null;
    if(!sort) {
      graphParams = rcp(new Teuchos::ParameterList);
      graphParams->set("sorted", false);
    }

    return rcp (new crs_graph_type (lclGraphSym,
                                    row_map,
                                    col_map_sym,
                                    domain_map,
                                    range_map,
                                    importer,
                                    Teuchos::null,
                                    graphParams));
  }

  template<class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
  Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
  CrsGraphTransposer<LocalOrdinal, GlobalOrdinal, Node>::
  createTranspose (const Teuchos::RCP<Teuchos::ParameterList> &params)
  {
    using Teuchos::RCP;
    // Do the local transpose
    RCP<crs_graph_type> transGraphWithSharedRows = createTransposeLocal (params);

#ifdef HAVE_TPETRA_MMM_TIMINGS
    const std::string prefix = std::string ("Tpetra ") + label_ + ": ";
    using Teuchos::TimeMonitor;
    TimeMonitor MM (*TimeMonitor::getNewTimer (prefix + "Transpose TAFC"));
#endif

    // If transGraphWithSharedRows has an exporter, that's what we
    // want.  If it doesn't, the rows aren't actually shared, and we're
    // done!
    using export_type = Export<LocalOrdinal, GlobalOrdinal, Node>;
    RCP<const export_type> exporter =
      transGraphWithSharedRows->getExporter ();
    if (exporter.is_null ()) {
      return transGraphWithSharedRows;
    }
    else {
      Teuchos::ParameterList labelList;
#ifdef HAVE_TPETRA_MMM_TIMINGS
      labelList.set("Timer Label", label_);
#endif
      if(! params.is_null ()) {
        const char paramName[] = "compute global constants";
        labelList.set (paramName, params->get (paramName, true));
      }
      // Use the Export object to do a fused Export and fillComplete.
      // This always sorts the local graph after communication, so
      //   no need to set "sorted = false" in parameters.
      return exportAndFillCompleteCrsGraph<crs_graph_type>
        (transGraphWithSharedRows, *exporter, Teuchos::null,
         Teuchos::null, Teuchos::rcpFromRef (labelList));
    }
  }

  template<class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
  Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
  CrsGraphTransposer<LocalOrdinal, GlobalOrdinal, Node>::
  createTransposeLocal (const Teuchos::RCP<Teuchos::ParameterList> &params)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using import_type = Tpetra::Import<LO, GO, Node>;
    using export_type = Tpetra::Export<LO, GO, Node>;

#ifdef HAVE_TPETRA_MMM_TIMINGS
    std::string prefix = std::string("Tpetra ") + label_ + ": ";
    using Teuchos::TimeMonitor;
    TimeMonitor MM (*TimeMonitor::getNewTimer (prefix + "Transpose Local"));
#endif

    const bool sort = [&] () {
      constexpr bool sortDefault = true; // see #4607 discussion
      const char sortParamName[] = "sort";
      return params.get () == nullptr ? sortDefault :
        params->get (sortParamName, sortDefault);
    } ();

    using local_graph_device_type = typename crs_graph_type::local_graph_device_type;
    local_graph_device_type lclGraph = origGraph_->getLocalGraphDevice ();

    //Allocate views and call the other version of transpose_graph
    using c_rowmap_t = typename local_graph_device_type::row_map_type;
    using c_entries_t = typename local_graph_device_type::entries_type;
    using rowmap_t = typename local_graph_device_type::row_map_type::non_const_type;
    using entries_t = typename local_graph_device_type::entries_type::non_const_type;
    LocalOrdinal numCols = origGraph_->getColMap()->getLocalNumElements();
    rowmap_t lclGraphT_rowmap("Transpose rowmap", numCols + 1);
    entries_t lclGraphT_entries(
                                Kokkos::ViewAllocateWithoutInitializing("Transpose entries"), lclGraph.entries.extent(0));
    KokkosSparse::Impl::transpose_graph<
      c_rowmap_t, c_entries_t,
      rowmap_t, entries_t,
      rowmap_t, typename local_graph_device_type::execution_space>(
                                                                   lclGraph.numRows(), numCols,
                                                                   lclGraph.row_map, lclGraph.entries,
                                                                   lclGraphT_rowmap, lclGraphT_entries);

    if (sort)
      KokkosSparse::sort_crs_graph<
        typename local_graph_device_type::execution_space,
        rowmap_t, entries_t>(
                             lclGraphT_rowmap,
                             lclGraphT_entries);

    //And construct the transpose local_graph_device_type
    local_graph_device_type lclGraphT = local_graph_device_type(lclGraphT_entries, lclGraphT_rowmap);

    // Prebuild the importers and exporters the no-communication way,
    // flipping the importers and exporters around.
    const auto origExport = origGraph_->getExporter ();
    RCP<const import_type> myImport = origExport.is_null () ?
      Teuchos::null : rcp (new import_type (*origExport));
    const auto origImport = origGraph_->getImporter ();
    RCP<const export_type> myExport = origImport.is_null () ?
      Teuchos::null : rcp (new export_type (*origImport));

    RCP<Teuchos::ParameterList> graphParams = Teuchos::null;
    if(!sort) {
      graphParams = rcp(new Teuchos::ParameterList);
      graphParams->set("sorted", false);
    }

    return rcp (new crs_graph_type (lclGraphT,
                                    origGraph_->getColMap (),
                                    origGraph_->getRowMap (),
                                    origGraph_->getRangeMap (),
                                    origGraph_->getDomainMap (),
                                    myImport, myExport, graphParams));
  }

  //
  // Explicit instantiation macro
  //
  // Must be expanded from within the Tpetra namespace!
  //

#define TPETRA_CRSGRAPHTRANSPOSER_INSTANT(LO,GO,NODE)   \
  template class CrsGraphTransposer< LO , GO , NODE >;

} // namespace Tpetra

#endif
