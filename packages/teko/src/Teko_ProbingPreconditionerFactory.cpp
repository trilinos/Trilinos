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

  using rowmap_t     = typename Tpetra::CrsGraph<LO, GO, NT>::local_graph_device_type::row_map_type;
  using colinds_t    = typename Tpetra::CrsGraph<LO, GO, NT>::local_graph_device_type::entries_type;
  using exec_space   = typename NT::device_type::execution_space;
  using memory_space = typename NT::device_type::memory_space;
  using kernel_handle_type =
      KokkosKernels::Experimental::KokkosKernelsHandle<typename rowmap_t::value_type,
                                                       typename colinds_t::value_type, ST,
                                                       exec_space, memory_space, memory_space>;

  kernel_handle_type kh;
  kh.create_distance2_graph_coloring_handle(KokkosGraph::COLORING_D2_DEFAULT);

  {
    auto localGraphDevice = graph_->getLocalGraphDevice();
    KokkosGraph::Experimental::graph_color_distance2(
        &kh, graph_->getLocalNumRows(), localGraphDevice.row_map, localGraphDevice.entries);
  }

  auto coloringHandle = kh.get_distance2_graph_coloring_handle();
  auto colors         = coloringHandle->get_vertex_colors();
  const LO numColors  = static_cast<LO>(coloringHandle->get_num_colors());
  auto colors_h       = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), colors);

  auto probedMat = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(graph_->getRowMap(), 0));

  auto localGraphHost = graph_->getLocalGraphHost();
  const auto colMap   = graph_->getColMap();

  Teuchos::Array<GO> colGids;
  Teuchos::Array<ST> vals;

  for (LO color = 1; color <= numColors; ++color) {
    RCP<Tpetra::MultiVector<ST, LO, GO, NT> > probeVec =
        rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(domainMap, 1));
    probeVec->putScalar(0.0);

    for (LO lrow = 0; lrow < static_cast<LO>(graph_->getLocalNumRows()); ++lrow) {
      if (colors_h(lrow) == color) {
        GO gid = domainMap->getGlobalElement(lrow);
        probeVec->replaceGlobalValue(gid, 0, Teuchos::ScalarTraits<ST>::one());
      }
    }

    RCP<Tpetra::MultiVector<ST, LO, GO, NT> > response =
        rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(rangeMap, 1));
    response->putScalar(0.0);

    Thyra::apply(*lo, Thyra::NOTRANS,
                 *Thyra::createConstMultiVector<ST, LO, GO, NT>(
                     probeVec, Thyra::tpetraVectorSpace<ST, LO, GO, NT>(domainMap)),
                 Thyra::createMultiVector<ST, LO, GO, NT>(
                     response, Thyra::tpetraVectorSpace<ST, LO, GO, NT>(rangeMap))
                     .ptr());

    auto responseView = response->getLocalViewHost(Tpetra::Access::ReadOnly);

    for (LO lrow = 0; lrow < static_cast<LO>(graph_->getLocalNumRows()); ++lrow) {
      auto row      = localGraphHost.rowConst(lrow);
      size_t rowLen = row.length;

      colGids.resize(rowLen);
      vals.resize(rowLen);

      GO rowGid = graph_->getRowMap()->getGlobalElement(lrow);

      size_t nnz = 0;
      for (size_t k = 0; k < rowLen; ++k) {
        LO lcol   = row.colidx(k);
        GO colGid = colMap->getGlobalElement(lcol);

        if (colors_h(lcol) == color) {
          colGids[nnz] = colGid;
          vals[nnz]    = responseView(lrow, 0);
          ++nnz;
        }
      }

      if (nnz > 0) {
        probedMat->sumIntoGlobalValues(rowGid,
                                       Teuchos::ArrayView<const GO>(colGids.getRawPtr(), nnz),
                                       Teuchos::ArrayView<const ST>(vals.getRawPtr(), nnz));
      }
    }
  }

  probedMat->fillComplete(domainMap, rangeMap);
  kh.destroy_distance2_graph_coloring_handle();

  return probedMat;
}

}  // namespace Teko
