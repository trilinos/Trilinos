// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LWGRAPHBASE_DECL_HPP
#define MUELU_LWGRAPHBASE_DECL_HPP

#include "Kokkos_Bitset.hpp"
#include "MueLu_ConfigDefs.hpp"

#include <Kokkos_StaticCrsGraph.hpp>
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <KokkosCompat_View.hpp>

#include <Xpetra_ConfigDefs.hpp>  // global_size_t
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <type_traits>

#include "MueLu_VerbosityLevel.hpp"
#include "MueLu_LWGraph_kokkos_fwd.hpp"

#include "MueLu_Exceptions.hpp"

namespace MueLu {

namespace {  // anonymous

template <class LocalOrdinal, class RowType>
class MaxNumRowEntriesFunctor {
 public:
  MaxNumRowEntriesFunctor(RowType rowPointers)
    : rowPointers_(rowPointers) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const LocalOrdinal i, size_t& maxLength) const {
    size_t d = rowPointers_(i + 1) - rowPointers_(i);

    maxLength = (d > maxLength ? d : maxLength);
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile size_t& dest, const volatile size_t& src) {
    dest = (dest > src ? dest : src);
  }

  KOKKOS_INLINE_FUNCTION
  void init(size_t& initValue) {
    initValue = 0;
  }

 private:
  RowType rowPointers_;
};

}  // namespace

/*!
  @class LWGraph_kokkos
  @brief Lightweight MueLu representation of a compressed row storage graph

  This class is lightweight in the sense that it holds to local graph
  information. These were built without using fillComplete.
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node, bool OnHost>
class LWGraphBase {
 public:
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using map_type            = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using crs_graph_type      = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using size_type           = size_t;

  using device_type     = typename std::conditional<OnHost,
                                                Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                                                typename Node::device_type>::type;
  using execution_space = typename device_type::execution_space;
  using memory_space    = typename device_type::memory_space;

  using local_graph_device_type = Kokkos::StaticCrsGraph<LocalOrdinal,
                                                         Kokkos::LayoutLeft,
                                                         typename Node::device_type,
                                                         void, size_t>;
  using local_graph_type        = typename std::conditional<OnHost, typename local_graph_device_type::HostMirror, local_graph_device_type>::type;
  using boundary_nodes_type     = Kokkos::View<bool*, memory_space>;
  using row_type                = typename local_graph_type::row_map_type;
  using entries_type            = typename local_graph_type::entries_type;
  using neighbor_vertices_type  = Kokkos::GraphRowViewConst<local_graph_type>;

#undef MUELU_LWGRAPHBASE_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  //! @name Constructors/Destructors.
  //@{

 private:
  void setup(const local_graph_type& graph,
             const RCP<const map_type>& domainMap,
             const RCP<const map_type>& importMap,
             const std::string& objectLabel) {
    using range_type = Kokkos::RangePolicy<local_ordinal_type, execution_space>;

    graph_       = graph;
    domainMap_   = domainMap;
    importMap_   = importMap;
    objectLabel_ = objectLabel;

    minLocalIndex_ = domainMap_->getMinLocalIndex();
    maxLocalIndex_ = domainMap_->getMaxLocalIndex();

    MaxNumRowEntriesFunctor<LO, typename local_graph_type::row_map_type> maxNumRowEntriesFunctor(graph_.row_map);
    Kokkos::parallel_reduce("MueLu:LocalLWGraph:LWGraph:maxnonzeros", range_type(0, graph_.numRows()), maxNumRowEntriesFunctor, maxNumRowEntries_);
  }

 public:
  //! LWGraph constructor
  //
  // @param[in] graph: local graph of type Kokkos::StaticCrsGraph containing CRS data
  // @param[in] domainMap: non-overlapping (domain) map for graph. Usually provided by AmalgamationFactory stored in UnAmalgamationInfo container
  // @param[in] importMap: overlapping map for graph. Usually provided by AmalgamationFactory stored in UnAmalgamationInfo container
  // @param[in] objectLabel: label string
  LWGraphBase(const local_graph_type& graph,
              const RCP<const map_type>& domainMap,
              const RCP<const map_type>& importMap,
              const std::string& objectLabel = "") {
    setup(graph, domainMap, importMap, objectLabel);
  }

  LWGraphBase(const RCP<const crs_graph_type>& graph,
              const std::string& objectLabel = "") {
    if constexpr (OnHost) {
      // We want the graph data to live on host
      if constexpr (std::is_same<local_graph_type,
                                 typename crs_graph_type::local_graph_type>::value)
        // The CrsGraph's data already lives on host.
        setup(graph->getLocalGraphHost(), graph->getRowMap(), graph->getColMap(), objectLabel);
      else {
        // We deep-copy the graph to host once during construction instead of using the host mirror.
        // This avoids issues with keeping a reference of the host mirror around.
        auto lclGraphDevice = graph->getLocalGraphDevice();
        auto rows           = typename local_graph_type::row_map_type::non_const_type("rows", lclGraphDevice.numRows() + 1);
        auto columns        = typename local_graph_type::entries_type::non_const_type("columns", lclGraphDevice.entries.extent(0));
        Kokkos::deep_copy(rows, lclGraphDevice.row_map);
        Kokkos::deep_copy(columns, lclGraphDevice.entries);
        local_graph_type lclGraph(columns, rows);
        setup(lclGraph, graph->getRowMap(), graph->getColMap(), objectLabel);
      }
    } else {
      // We want the graph data on device.
      setup(graph->getLocalGraphDevice(), graph->getRowMap(), graph->getColMap(), objectLabel);
    }
  }

  LWGraphBase(const row_type& rows,
              const entries_type& columns,
              const RCP<const map_type>& domainMap,
              const RCP<const map_type>& importMap,
              const std::string& objectLabel = "") {
    const local_graph_type graph(columns, rows);
    setup(graph, domainMap, importMap, objectLabel);
  }

  ~LWGraphBase() = default;
  //@}

  const RCP<const Teuchos::Comm<int>> GetComm() const {
    return domainMap_->getComm();
  }
  const RCP<const Map> GetDomainMap() const {
    return domainMap_;
  }
  //! Return overlapping import map (nodes).
  const RCP<const Map> GetImportMap() const {
    return importMap_;
  }

  //! Return number of graph vertices
  KOKKOS_INLINE_FUNCTION size_type GetNodeNumVertices() const {
    return graph_.numRows();
  }

  //! Return number of graph edges
  KOKKOS_INLINE_FUNCTION size_type GetNodeNumEdges() const {
    return graph_.row_map(GetNodeNumVertices());
  }

  //! Return global number of graph edges
  Xpetra::global_size_t GetGlobalNumEdges() const {
    Xpetra::global_size_t in = GetNodeNumEdges(), out;
    Teuchos::reduceAll(*domainMap_->getComm(), Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
    return out;
  }

  //! Returns the maximum number of entries across all rows/columns on this node
  KOKKOS_INLINE_FUNCTION size_type getLocalMaxNumRowEntries() const {
    return maxNumRowEntries_;
  }

  //! Return the list of vertices adjacent to the vertex 'v'.
  KOKKOS_INLINE_FUNCTION neighbor_vertices_type getNeighborVertices(LO i) const {
    auto rowView = graph_.rowConst(i);
    return rowView;
  }

  //! Return the list of vertices adjacent to the vertex 'v'.
  Teuchos::ArrayView<LO> getNeighborVertices_av(LO i) const {
    return Kokkos::Compat::getArrayView(Kokkos::subview(graph_.entries,
                                                        Kokkos::make_pair(graph_.row_map(i),
                                                                          graph_.row_map(i + 1))));
  }

  //! Return true if vertex with local id 'v' is on current process.
  KOKKOS_INLINE_FUNCTION bool isLocalNeighborVertex(LO i) const {
    return i >= minLocalIndex_ && i <= maxLocalIndex_;
  }

  //! Return the row pointers of the local graph
  KOKKOS_INLINE_FUNCTION row_type getRowPtrs() const {
    return graph_.row_map;
  }

  //! Return the list entries in the local graph
  KOKKOS_INLINE_FUNCTION entries_type getEntries() const {
    return graph_.entries;
  }

  //! Set boolean array indicating which rows correspond to Dirichlet boundaries.
  KOKKOS_INLINE_FUNCTION void SetBoundaryNodeMap(const boundary_nodes_type bndry) {
    dirichletBoundaries_ = bndry;
  }

  //! Returns map with global ids of boundary nodes.
  KOKKOS_INLINE_FUNCTION const boundary_nodes_type GetBoundaryNodeMap() const {
    return dirichletBoundaries_;
  }

  /// Return a simple one-line description of the Graph.
  std::string description() const {
    return "LWGraphBase (" + objectLabel_ + ")";
  }

  //! Print the Graph with some verbosity level to an FancyOStream object.
  void print(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const {
    if (verbLevel & Debug) {
      auto graph             = graph_;
      RCP<const Map> col_map = importMap_.is_null() ? domainMap_ : importMap_;
      int mypid              = col_map->getComm()->getRank();

      {
        std::ostringstream ss;
        ss << "[pid " << mypid << "] num entries=" << graph.entries.size();
        out << ss.str() << std::endl;
      }

      const size_t numRows = graph.numRows();
      auto rowPtrs         = graph.row_map;
      auto columns         = graph.entries;
      auto rowPtrs_h       = Kokkos::create_mirror_view(rowPtrs);
      auto columns_h       = Kokkos::create_mirror_view(columns);
      Kokkos::deep_copy(rowPtrs_h, rowPtrs);
      Kokkos::deep_copy(columns_h, columns);
      for (size_t i = 0; i < numRows; ++i) {
        std::ostringstream ss;
        ss << "[pid " << mypid << "] row " << domainMap_->getGlobalElement(i) << ":";
        ss << " (numEntries=" << rowPtrs_h(i + 1) - rowPtrs_h(i) << ")";

        for (size_t jj = rowPtrs_h(i); jj < rowPtrs_h(i + 1); jj++) {
          ss << " " << col_map->getGlobalElement(columns_h(jj));
        }
        out << ss.str() << std::endl;
      }
    }
  }

  local_graph_type& getGraph() const {
    return graph_;
  }

  RCP<crs_graph_type> GetCrsGraph() const {
    RCP<crs_graph_type> graph;
    if constexpr (std::is_same<local_graph_type,
                               typename crs_graph_type::local_graph_type>::value)
      graph = Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(GetDomainMap(), GetImportMap(), graph_, Teuchos::null);
    else {
      auto rows    = typename crs_graph_type::local_graph_type::row_map_type::non_const_type("rows", graph_.numRows() + 1);
      auto columns = typename crs_graph_type::local_graph_type::entries_type::non_const_type("columns", graph_.entries.extent(0));
      Kokkos::deep_copy(rows, graph_.row_map);
      Kokkos::deep_copy(columns, graph_.entries);
      typename crs_graph_type::local_graph_type lclGraph(columns, rows);
      graph = Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(GetDomainMap(), GetImportMap(), lclGraph, Teuchos::null);
    }
    return graph;
  }

  const std::string& getObjectLabel() const {
    return objectLabel_;
  }

 private:
  //! Underlying graph (with label)
  mutable local_graph_type graph_;

  //! Graph maps
  RCP<const map_type> domainMap_;
  RCP<const map_type> importMap_;

  //! Name of this graph.
  std::string objectLabel_;

  //! Boolean array marking Dirichlet rows.
  boundary_nodes_type dirichletBoundaries_;

  //! Local index boundaries (cached from domain map)
  LO minLocalIndex_, maxLocalIndex_;
  size_type maxNumRowEntries_;
};

}  // namespace MueLu

#define MUELU_LWGRAPHBASE_SHORT
#endif  // MUELU_LWGRAPHBASE_DECL_HPP
