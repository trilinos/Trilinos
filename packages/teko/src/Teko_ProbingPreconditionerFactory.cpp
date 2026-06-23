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

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace {

using ST = Teko::ST;
using LO = Teko::LO;
using GO = Teko::GO;
using NT = Teko::NT;

using exec_space = typename NT::device_type::execution_space;

template <class MVType, class HostColorViewType>
void set_probe_by_color_host(const Teuchos::RCP<MVType>& probeVec,
                             const HostColorViewType& h_colors,
                             const LO numLocalCols,
                             const LO color) {
  probeVec->putScalar(Teuchos::ScalarTraits<ST>::zero());

  auto hostView = probeVec->getLocalViewHost(Tpetra::Access::ReadWrite);
  for (LO lcol = 0; lcol < numLocalCols; ++lcol) {
    if (static_cast<LO>(h_colors(lcol)) == color) {
      hostView(lcol, 0) = Teuchos::ScalarTraits<ST>::one();
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
      rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT> >(
          tOp->getConstTpetraOperator(), true);
  setGraph(crsMatrix->getCrsGraph());
}

void ProbingPreconditionerFactory::setGraph(
    const Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, NT> >& graph) {
  graph_ = graph;
}

RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> >
ProbingPreconditionerFactory::probe(const LinearOp& lo) const {
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

  using local_graph_device_type =
      typename Tpetra::CrsGraph<LO, GO, NT>::local_graph_device_type;
  using row_map_type      = typename local_graph_device_type::row_map_type;
  using index_type        = typename local_graph_device_type::entries_type;
  using local_matrix_type =
      typename Tpetra::CrsMatrix<ST, LO, GO, NT>::local_matrix_device_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using exec_space        = typename NT::device_type::execution_space;
  using memory_space      = typename NT::device_type::memory_space;
  using row_map_nonc_type = typename row_map_type::non_const_type;
  using entries_nonc_type = typename index_type::non_const_type;

  using kernel_handle_type =
      KokkosKernels::Experimental::KokkosKernelsHandle<
          typename row_map_type::value_type,
          typename index_type::value_type,
          ST,
          exec_space,
          memory_space,
          memory_space>;

  const LO numRows      = static_cast<LO>(graph_->getLocalNumRows());
  const LO numLocalCols = static_cast<LO>(colMap->getLocalNumElements());

  TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t>(numLocalCols) != domainMap->getLocalNumElements(),
      std::runtime_error,
      "ProbingPreconditionerFactory::probe: graph local column count ("
          << numLocalCols << ") does not match operator domain local size ("
          << domainMap->getLocalNumElements()
          << "). The probing code assumes local column indices index the domain vector.");

  row_map_nonc_type row_map_copy;
  entries_nonc_type entries_copy;
  LO totalNumEntries     = 0;
  LO maxNumEntriesPerRow = 0;

  {
    auto localGraphDevice = graph_->getLocalGraphDevice();

    row_map_copy = row_map_nonc_type(
        Kokkos::ViewAllocateWithoutInitializing("teko_probe_row_map_copy"),
        localGraphDevice.row_map.extent(0));
    entries_copy = entries_nonc_type(
        Kokkos::ViewAllocateWithoutInitializing("teko_probe_entries_copy"),
        localGraphDevice.entries.extent(0));

    Kokkos::deep_copy(row_map_copy, localGraphDevice.row_map);
    Kokkos::deep_copy(entries_copy, localGraphDevice.entries);

    totalNumEntries     = static_cast<LO>(entries_copy.extent(0));
    maxNumEntriesPerRow = compute_max_row_length(row_map_copy, numRows);
  }

  kernel_handle_type kh;
  kh.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_DEFAULT);

  KokkosGraph::Experimental::graph_color_distance2(
      &kh, numLocalCols, row_map_copy, entries_copy);

  auto coloringHandle = kh.get_distance2_graph_coloring_handle();
  auto colors         = coloringHandle->get_vertex_colors();
  const LO numColors  = static_cast<LO>(coloringHandle->get_num_colors());

  TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<LO>(colors.extent(0)) < numLocalCols,
      std::runtime_error,
      "ProbingPreconditionerFactory::probe: colors view extent ("
          << colors.extent(0) << ") is smaller than local column count ("
          << numLocalCols << ").");

  values_type values(Kokkos::ViewAllocateWithoutInitializing("probed_values"), totalNumEntries);
  Kokkos::deep_copy(values, Teuchos::ScalarTraits<ST>::zero());

  auto h_colors = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), colors);

  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > probeVec =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(domainMap, 1));
  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > response =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(rangeMap, 1));

  for (LO color = 1; color <= numColors; ++color) {
    set_probe_by_color_host(probeVec, h_colors, numLocalCols, color);

    tpetraOp->apply(*probeVec, *response, Teuchos::NO_TRANS,
                    Teuchos::ScalarTraits<ST>::one(),
                    Teuchos::ScalarTraits<ST>::zero());

    auto responseView = response->getLocalViewDevice(Tpetra::Access::ReadOnly);
    decode_probe_by_color(row_map_copy, entries_copy, colors,
                          responseView, numRows, color, values);
  }

  kh.destroy_distance2_graph_coloring_handle();

  auto lclMat = local_matrix_type("probed_local_matrix", numRows, maxNumEntriesPerRow,
                                  totalNumEntries, values,
                                  row_map_copy, entries_copy);

  return Teuchos::rcp(
      new Tpetra::CrsMatrix<ST, LO, GO, NT>(lclMat, rowMap, colMap, domainMap, rangeMap));
}

}  // namespace Teko