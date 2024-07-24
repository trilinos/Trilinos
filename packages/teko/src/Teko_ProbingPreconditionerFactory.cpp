// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_ProbingPreconditionerFactory.hpp"

#ifdef Teko_ENABLE_Isorropia

#include "Teko_EpetraOperatorWrapper.hpp"

#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Epetra_CrsMatrix.h"

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

namespace Teko {

/*****************************************************/

ProbingPreconditionerFactory::ProbingPreconditionerFactory() {
  prober = rcp(new Isorropia::Epetra::Prober);
}

LinearOp ProbingPreconditionerFactory::buildPreconditionerOperator(
    LinearOp& lo, PreconditionerState& state) const {
  // make an epetra operator to be probed
  RCP<Epetra_Operator> epetraLo = rcp(new Teko::Epetra::EpetraOperatorWrapper(lo));

  // build color scheme
  prober->color();

  // probe operator: take me to your leader
  RCP<Epetra_CrsMatrix> retOp = prober->probe(*epetraLo);
  Teko::LinearOp probedOp     = Thyra::epetraLinearOp(retOp);

  return Teko::buildInverse(*invFactory_, probedOp);
}

void ProbingPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList& pl) {
  RCP<const InverseLibrary> invLib = getInverseLibrary();

  const std::string inverse_type           = "Inverse Type";
  const std::string probing_graph_operator = "Probing Graph Operator";
  const std::string probing_graph          = "Probing Graph";
  const std::string user_graph             = "User Will Set Probing Graph";

  // get string specifying default inverse
  std::string invStr = "Amesos";
  if (pl.isParameter(inverse_type)) invStr = pl.get<std::string>(inverse_type);

  if (pl.isParameter(probing_graph_operator))
    setGraphOperator(pl.get<Teko::LinearOp>(probing_graph_operator));
  else if (pl.isParameter(probing_graph))
    setGraph(pl.get<RCP<const Epetra_CrsGraph> >(probing_graph));
  else if (pl.isParameter(user_graph) && pl.get<bool>("User Will Set Probing Graph")) {
    // noop
  } else {
    Teuchos::RCP<Teko::RequestHandler> rh = getRequestHandler();
    rh->preRequest<RCP<const Epetra_CrsGraph> >(Teko::RequestMesg("Probing Graph"));
    setGraph(rh->request<RCP<const Epetra_CrsGraph> >(Teko::RequestMesg("Probing Graph")));
  }

  setInverseFactory(invLib->getInverseFactory(invStr));
}

void ProbingPreconditionerFactory::setGraphOperator(const Teko::LinearOp& graphOp) {
  RCP<const Epetra_CrsMatrix> crsMatrix =
      rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*graphOp));
  setGraph(Teuchos::rcpFromRef(crsMatrix->Graph()));
}

void ProbingPreconditionerFactory::setGraph(const Teuchos::RCP<const Epetra_CrsGraph>& graph) {
  prober->setGraph(graph);
}

void ProbingPreconditionerFactory::setProberList(const Teuchos::ParameterList& list) {
  prober->setList(list);
}

}  // end namespace Teko

#endif
