// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_ProbingPreconditionerFactory.hpp"

#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestMesg.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include "KokkosGraph_Distance2ColorHandle.hpp"
#include "KokkosGraph_Distance2Color.hpp"
#include "KokkosKernels_Handle.hpp"

#include <Kokkos_Sort.hpp>
#include <cstdint>
#include <vector>

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace {

using ST = Teko::ST;
using LO = Teko::LO;
using GO = Teko::GO;
using NT = Teko::NT;

using exec_space = typename NT::device_type::execution_space;

template <class RowMapType, class EntriesType>
struct LocalCrsGraphViews {
  RowMapType row_map;
  EntriesType entries;
};

template <class MVType, class HostColorViewType>
void set_probe_by_color_host(const Teuchos::RCP<MVType>& probeVec,
                             const HostColorViewType& h_colors,
                             const std::vector<LO>& colLidToDomainLid, const LO color) {
  const LO invalid = Teuchos::OrdinalTraits<LO>::invalid();

  probeVec->putScalar(Teuchos::ScalarTraits<ST>::zero());

  auto hostView         = probeVec->getLocalViewHost(Tpetra::Access::ReadWrite);
  const LO numLocalCols = static_cast<LO>(colLidToDomainLid.size());

  for (LO lcol = 0; lcol < numLocalCols; ++lcol) {
    if (static_cast<LO>(h_colors(lcol)) == color) {
      const LO domainLid = colLidToDomainLid[static_cast<size_t>(lcol)];
      if (domainLid != invalid) {
        hostView(domainLid, 0) = Teuchos::ScalarTraits<ST>::one();
      }
    }
  }
}

template <class RowPtrViewType, class ColIndViewType, class ColorViewType, class ResponseViewType,
          class ValuesViewType>
void decode_probe_by_color(const RowPtrViewType& rowPtrs, const ColIndViewType& localColInds,
                           const ColorViewType& colors, const ResponseViewType& responseView,
                           const LO numRows, const LO color, const ValuesViewType& values) {
  Kokkos::parallel_for(
      "Teko::ProbingPreconditionerFactory::decode_probe",
      Kokkos::RangePolicy<exec_space>(0, numRows), KOKKOS_LAMBDA(const LO lrow) {
        const auto rowStart = rowPtrs(lrow);
        const auto rowEnd   = rowPtrs(lrow + 1);
        const ST rowValue   = responseView(lrow, 0);

        for (auto entry = rowStart; entry < rowEnd; ++entry) {
          const LO lcol = static_cast<LO>(localColInds(entry));
          if (colors(lcol) == color) {
            values(entry) = rowValue;
          }
        }
      });
}

template <class RowPtrViewType>
LO compute_max_row_length(const RowPtrViewType& rowPtrs, const LO numRows) {
  LO maxNumEntriesPerRow = 0;

  Kokkos::parallel_reduce(
      "Teko::ProbingPreconditionerFactory::max_row_length",
      Kokkos::RangePolicy<exec_space>(0, numRows),
      KOKKOS_LAMBDA(const LO r, LO& localMax) {
        const LO rowLen = static_cast<LO>(rowPtrs(r + 1) - rowPtrs(r));
        if (rowLen > localMax) localMax = rowLen;
      },
      Kokkos::Max<LO>(maxNumEntriesPerRow));

  return maxNumEntriesPerRow;
}

template <class RowMapViewType>
typename RowMapViewType::non_const_value_type count_candidate_col_graph_edges(
    const RowMapViewType& row_map_copy, const LO numRows) {
  using row_map_value_type = typename RowMapViewType::non_const_value_type;

  row_map_value_type candidateEdgeCount = 0;
  Kokkos::parallel_reduce(
      "Teko::ProbingPreconditionerFactory::count_candidate_col_graph_edges",
      Kokkos::RangePolicy<exec_space>(0, numRows),
      KOKKOS_LAMBDA(const LO row, row_map_value_type& update) {
        const row_map_value_type d = row_map_copy(row + 1) - row_map_copy(row);
        if (d > 1) update += d * (d - 1);
      },
      candidateEdgeCount);

  return candidateEdgeCount;
}

template <class RowMapDeviceView, class EntriesDeviceView>
LocalCrsGraphViews<Kokkos::View<typename RowMapDeviceView::non_const_value_type*,
                                typename RowMapDeviceView::device_type>,
                   Kokkos::View<typename EntriesDeviceView::non_const_value_type*,
                                typename EntriesDeviceView::device_type> >
build_column_intersection_graph_device(const RowMapDeviceView& row_map_in,
                                       const EntriesDeviceView& entries_in, const LO numRows,
                                       const LO numLocalCols) {
  using row_offset_type = typename RowMapDeviceView::non_const_value_type;
  using col_index_type  = typename EntriesDeviceView::non_const_value_type;
  using device_type     = typename RowMapDeviceView::device_type;
  using key_type        = std::uint64_t;

  using row_map_out_type = Kokkos::View<row_offset_type*, device_type>;
  using entries_out_type = Kokkos::View<col_index_type*, device_type>;

  const row_offset_type numCandidateKeys = count_candidate_col_graph_edges(row_map_in, numRows);

  Kokkos::View<row_offset_type*, device_type> candidateOffsets(
      Kokkos::ViewAllocateWithoutInitializing("teko_probe_candidate_offsets"),
      static_cast<size_t>(numRows) + 1);

  Kokkos::parallel_for(
      "Teko::ProbingPreconditionerFactory::fill_candidate_offsets",
      Kokkos::RangePolicy<exec_space>(0, numRows), KOKKOS_LAMBDA(const LO row) {
        const row_offset_type d = row_map_in(row + 1) - row_map_in(row);
        candidateOffsets(row)   = (d > 1 ? d * (d - 1) : 0);
      });

  Kokkos::parallel_scan(
      "Teko::ProbingPreconditionerFactory::scan_candidate_offsets",
      Kokkos::RangePolicy<exec_space>(0, numRows + 1),
      KOKKOS_LAMBDA(const LO i, row_offset_type& update, const bool final_pass) {
        const row_offset_type val = (i < numRows ? candidateOffsets(i) : row_offset_type(0));
        if (final_pass) candidateOffsets(i) = update;
        update += val;
      });

  Kokkos::View<key_type*, device_type> candidateKeys(
      Kokkos::ViewAllocateWithoutInitializing("teko_probe_candidate_keys"),
      static_cast<size_t>(numCandidateKeys));

  const key_type numColsKey = static_cast<key_type>(numLocalCols);

  Kokkos::parallel_for(
      "Teko::ProbingPreconditionerFactory::fill_candidate_keys",
      Kokkos::RangePolicy<exec_space>(0, numRows), KOKKOS_LAMBDA(const LO row) {
        const auto rowStart     = row_map_in(row);
        const auto rowEnd       = row_map_in(row + 1);
        const row_offset_type d = rowEnd - rowStart;
        row_offset_type pos     = candidateOffsets(row);

        for (row_offset_type i = 0; i < d; ++i) {
          const key_type src = static_cast<key_type>(entries_in(rowStart + i));
          for (row_offset_type j = 0; j < d; ++j) {
            if (i == j) continue;
            const key_type dst   = static_cast<key_type>(entries_in(rowStart + j));
            candidateKeys(pos++) = src * numColsKey + dst;
          }
        }
      });

  Kokkos::sort(candidateKeys);

  Kokkos::View<row_offset_type*, device_type> uniqueFlags(
      Kokkos::ViewAllocateWithoutInitializing("teko_probe_unique_flags"),
      static_cast<size_t>(numCandidateKeys));

  Kokkos::parallel_for(
      "Teko::ProbingPreconditionerFactory::mark_unique_candidate_keys",
      Kokkos::RangePolicy<exec_space>(0, numCandidateKeys), KOKKOS_LAMBDA(const row_offset_type i) {
        if (i == 0) {
          uniqueFlags(i) = 1;
        } else {
          uniqueFlags(i) = (candidateKeys(i) != candidateKeys(i - 1) ? 1 : 0);
        }
      });

  Kokkos::View<row_offset_type*, device_type> uniqueOffsets(
      Kokkos::ViewAllocateWithoutInitializing("teko_probe_unique_offsets"),
      static_cast<size_t>(numCandidateKeys) + 1);

  Kokkos::parallel_scan(
      "Teko::ProbingPreconditionerFactory::scan_unique_flags",
      Kokkos::RangePolicy<exec_space>(0, numCandidateKeys + 1),
      KOKKOS_LAMBDA(const row_offset_type i, row_offset_type& update, const bool final_pass) {
        const row_offset_type val = (i < numCandidateKeys ? uniqueFlags(i) : row_offset_type(0));
        if (final_pass) uniqueOffsets(i) = update;
        update += val;
      });

  auto h_uniqueCount = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), uniqueOffsets);
  const row_offset_type numUniqueKeys = h_uniqueCount(numCandidateKeys);

  Kokkos::View<key_type*, device_type> uniqueKeys(
      Kokkos::ViewAllocateWithoutInitializing("teko_probe_unique_keys"),
      static_cast<size_t>(numUniqueKeys));

  Kokkos::parallel_for(
      "Teko::ProbingPreconditionerFactory::pack_unique_keys",
      Kokkos::RangePolicy<exec_space>(0, numCandidateKeys), KOKKOS_LAMBDA(const row_offset_type i) {
        if (uniqueFlags(i)) {
          uniqueKeys(uniqueOffsets(i)) = candidateKeys(i);
        }
      });

  Kokkos::View<row_offset_type*, device_type> rowCounts(
      Kokkos::ViewAllocateWithoutInitializing("teko_probe_row_counts"),
      static_cast<size_t>(numLocalCols));
  Kokkos::deep_copy(rowCounts, row_offset_type(0));

  Kokkos::parallel_for(
      "Teko::ProbingPreconditionerFactory::count_unique_keys_per_row",
      Kokkos::RangePolicy<exec_space>(0, numUniqueKeys), KOKKOS_LAMBDA(const row_offset_type i) {
        const key_type key        = uniqueKeys(i);
        const row_offset_type src = static_cast<row_offset_type>(key / numColsKey);
        Kokkos::atomic_fetch_add(&rowCounts(src), row_offset_type(1));
      });

  row_map_out_type row_map_out(
      Kokkos::ViewAllocateWithoutInitializing("teko_probe_col_graph_row_map_dev"),
      static_cast<size_t>(numLocalCols) + 1);

  Kokkos::parallel_scan(
      "Teko::ProbingPreconditionerFactory::scan_col_graph_row_map",
      Kokkos::RangePolicy<exec_space>(0, numLocalCols + 1),
      KOKKOS_LAMBDA(const LO i, row_offset_type& update, const bool final_pass) {
        const row_offset_type val = (i < numLocalCols ? rowCounts(i) : row_offset_type(0));
        if (final_pass) row_map_out(i) = update;
        update += val;
      });

  entries_out_type entries_out(
      Kokkos::ViewAllocateWithoutInitializing("teko_probe_col_graph_entries_dev"),
      static_cast<size_t>(numUniqueKeys));

  Kokkos::View<row_offset_type*, device_type> fillOffsets(
      Kokkos::ViewAllocateWithoutInitializing("teko_probe_fill_offsets"),
      static_cast<size_t>(numLocalCols));

  Kokkos::parallel_for(
      "Teko::ProbingPreconditionerFactory::init_fill_offsets",
      Kokkos::RangePolicy<exec_space>(0, numLocalCols),
      KOKKOS_LAMBDA(const LO col) { fillOffsets(col) = row_map_out(col); });

  Kokkos::parallel_for(
      "Teko::ProbingPreconditionerFactory::fill_col_graph_entries",
      Kokkos::RangePolicy<exec_space>(0, numUniqueKeys), KOKKOS_LAMBDA(const row_offset_type i) {
        const key_type key        = uniqueKeys(i);
        const row_offset_type src = static_cast<row_offset_type>(key / numColsKey);
        const col_index_type dst  = static_cast<col_index_type>(key % numColsKey);
        const row_offset_type pos = Kokkos::atomic_fetch_add(&fillOffsets(src), row_offset_type(1));
        entries_out(pos)          = dst;
      });

  LocalCrsGraphViews<row_map_out_type, entries_out_type> out = {row_map_out, entries_out};
  return out;
}

}  // namespace

namespace Teko {

ProbingPreconditionerFactory::ProbingPreconditionerFactory() {}

LinearOp ProbingPreconditionerFactory::buildPreconditionerOperator(
    LinearOp& lo, PreconditionerState& state) const {
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > probedMat = probe(lo);

  LinearOp probedOp = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(probedMat->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(probedMat->getRangeMap()), probedMat);

  return Teko::buildInverse(*invFactory_, probedOp);
}

void ProbingPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList& pl) {
  RCP<const InverseLibrary> invLib = getInverseLibrary();

  const std::string inverse_type           = "Inverse Type";
  const std::string probing_graph_operator = "Probing Graph Operator";
  const std::string probing_graph          = "Probing Graph";
  const std::string user_graph             = "User Will Set Probing Graph";

  std::string invStr = "Ifpack2";
  if (pl.isParameter(inverse_type)) invStr = pl.get<std::string>(inverse_type);

  if (pl.isParameter(probing_graph_operator))
    setGraphOperator(pl.get<Teko::LinearOp>(probing_graph_operator));
  else if (pl.isParameter(probing_graph))
    setGraph(pl.get<RCP<const Tpetra::CrsGraph<LO, GO, NT> > >(probing_graph));
  else if (pl.isParameter(user_graph) && pl.get<bool>(user_graph)) {
    // noop
  } else {
    Teuchos::RCP<Teko::RequestHandler> rh = getRequestHandler();
    rh->preRequest<RCP<const Tpetra::CrsGraph<LO, GO, NT> > >(Teko::RequestMesg("Probing Graph"));
    setGraph(
        rh->request<RCP<const Tpetra::CrsGraph<LO, GO, NT> > >(Teko::RequestMesg("Probing Graph")));
  }

  setInverseFactory(invLib->getInverseFactory(invStr));
}

void ProbingPreconditionerFactory::setGraphOperator(const Teko::LinearOp& graphOp) {
  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp =
      rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(graphOp, true);
  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > crsMatrix =
      rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT> >(tOp->getConstTpetraOperator(),
                                                                 true);
  setGraph(crsMatrix->getCrsGraph());
}

void ProbingPreconditionerFactory::setGraph(
    const Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, NT> >& graph) {
  graph_ = graph;
}

RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ProbingPreconditionerFactory::probe(
    const LinearOp& lo) const {
  TEUCHOS_TEST_FOR_EXCEPTION(graph_ == Teuchos::null, std::runtime_error,
                             "ProbingPreconditionerFactory::probe: probing graph is null");

  RCP<const Thyra::LinearOpBase<ST> > thyraOp = lo;
  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp =
      rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(thyraOp, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > tpetraOp = tOp->getConstTpetraOperator();

  RCP<const Tpetra::Map<LO, GO, NT> > domainMap = tpetraOp->getDomainMap();
  RCP<const Tpetra::Map<LO, GO, NT> > rangeMap  = tpetraOp->getRangeMap();
  RCP<const Tpetra::Map<LO, GO, NT> > rowMap    = graph_->getRowMap();
  RCP<const Tpetra::Map<LO, GO, NT> > colMap    = graph_->getColMap();

  using local_graph_device_type = typename Tpetra::CrsGraph<LO, GO, NT>::local_graph_device_type;
  using row_map_type            = typename local_graph_device_type::row_map_type;
  using index_type              = typename local_graph_device_type::entries_type;
  using local_matrix_type = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::local_matrix_device_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using exec_space        = typename NT::device_type::execution_space;
  using memory_space      = typename NT::device_type::memory_space;
  using row_map_nonc_type = typename row_map_type::non_const_type;
  using entries_nonc_type = typename index_type::non_const_type;

  using kernel_handle_type =
      KokkosKernels::Experimental::KokkosKernelsHandle<typename row_map_type::value_type,
                                                       typename index_type::value_type, ST,
                                                       exec_space, memory_space, memory_space>;

  const LO numRows      = static_cast<LO>(graph_->getLocalNumRows());
  const LO numLocalCols = static_cast<LO>(colMap->getLocalNumElements());

  row_map_nonc_type row_map_copy;
  entries_nonc_type entries_copy;
  LO totalNumEntries     = 0;
  LO maxNumEntriesPerRow = 0;

  {
    auto localGraphDevice = graph_->getLocalGraphDevice();

    row_map_copy =
        row_map_nonc_type(Kokkos::ViewAllocateWithoutInitializing("teko_probe_row_map_copy"),
                          localGraphDevice.row_map.extent(0));
    entries_copy =
        entries_nonc_type(Kokkos::ViewAllocateWithoutInitializing("teko_probe_entries_copy"),
                          localGraphDevice.entries.extent(0));

    Kokkos::deep_copy(row_map_copy, localGraphDevice.row_map);
    Kokkos::deep_copy(entries_copy, localGraphDevice.entries);

    totalNumEntries     = static_cast<LO>(entries_copy.extent(0));
    maxNumEntriesPerRow = compute_max_row_length(row_map_copy, numRows);
  }

  auto colGraph =
      build_column_intersection_graph_device(row_map_copy, entries_copy, numRows, numLocalCols);

  kernel_handle_type kh;
  kh.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_DEFAULT);

  KokkosGraph::Experimental::graph_color_distance2(&kh, numLocalCols, colGraph.row_map,
                                                   colGraph.entries);

  auto coloringHandle = kh.get_distance2_graph_coloring_handle();
  auto colors         = coloringHandle->get_vertex_colors();
  const LO numColors  = static_cast<LO>(coloringHandle->get_num_colors());

  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<LO>(colors.extent(0)) < numLocalCols, std::runtime_error,
                             "ProbingPreconditionerFactory::probe: colors view extent ("
                                 << colors.extent(0) << ") is smaller than local column count ("
                                 << numLocalCols << ").");

  std::vector<LO> colLidToDomainLid(static_cast<size_t>(numLocalCols));
  {
    const LO invalid = Teuchos::OrdinalTraits<LO>::invalid();
    for (LO lcol = 0; lcol < numLocalCols; ++lcol) {
      const GO gid = colMap->getGlobalElement(lcol);
      const LO dlid =
          (gid == Teuchos::OrdinalTraits<GO>::invalid() ? invalid
                                                        : domainMap->getLocalElement(gid));
      colLidToDomainLid[static_cast<size_t>(lcol)] = dlid;
    }
  }

  values_type values(Kokkos::ViewAllocateWithoutInitializing("probed_values"), totalNumEntries);
  Kokkos::deep_copy(values, Teuchos::ScalarTraits<ST>::zero());

  auto h_colors = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), colors);

  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > probeVec =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(domainMap, 1));
  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > response =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(rangeMap, 1));

  for (LO color = 1; color <= numColors; ++color) {
    set_probe_by_color_host(probeVec, h_colors, colLidToDomainLid, color);

    tpetraOp->apply(*probeVec, *response, Teuchos::NO_TRANS, Teuchos::ScalarTraits<ST>::one(),
                    Teuchos::ScalarTraits<ST>::zero());

    auto responseView = response->getLocalViewDevice(Tpetra::Access::ReadOnly);
    decode_probe_by_color(row_map_copy, entries_copy, colors, responseView, numRows, color, values);
  }

  kh.destroy_graph_coloring_handle();

  auto lclMat = local_matrix_type("probed_local_matrix", numRows, maxNumEntriesPerRow,
                                  totalNumEntries, values, row_map_copy, entries_copy);

  return Teuchos::rcp(
      new Tpetra::CrsMatrix<ST, LO, GO, NT>(lclMat, rowMap, colMap, domainMap, rangeMap));
}

}  // namespace Teko