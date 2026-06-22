// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_ProbingPreconditionerFactory_hpp__
#define __Teko_ProbingPreconditionerFactory_hpp__

#include "Teko_PreconditionerState.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_ConfigDefs.hpp"

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Teko {

/** \brief Preconditioner factory that builds a probing approximation to an operator and then
 * "inverts" it as requested.
 *
 * NOTE: You need to set the probing Graph or GraphOperator for probing to know what to do.
 */
class ProbingPreconditionerFactory : public virtual Teko::PreconditionerFactory {
 public:
  //! @name Constructors.
  //@{

  /** Build an empty probing preconditioner factory
   */
  ProbingPreconditionerFactory();

  //@}

  /** Create the probed preconditioner operator.
   */
  LinearOp buildPreconditionerOperator(LinearOp& lo, PreconditionerState& state) const override;

  //! Initialize from a parameter list
  void initializeFromParameterList(const Teuchos::ParameterList& pl) override;

  void setGraphOperator(const Teko::LinearOp& graphOp);
  void setGraph(const Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, NT> >& graph);

  void setInverseFactory(const Teuchos::RCP<Teko::InverseFactory>& invFactory) {
    invFactory_ = invFactory;
  }

 protected:
  Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, NT> > graph_;
  Teuchos::RCP<Teko::InverseFactory> invFactory_;

 private:
  Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > probe(const LinearOp& lo) const;
};

}  // end namespace Teko

#endif
