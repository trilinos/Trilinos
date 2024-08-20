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

#include "Teko_Config.h"

#ifdef Teko_ENABLE_Isorropia

// Teko includes
#include "Teko_PreconditionerState.hpp"
#include "Teko_PreconditionerFactory.hpp"

// Isorropia includes
#include "Isorropia_EpetraProber.hpp"

namespace Teko {

/** \brief Preconditioner factory that builds a probing approximation to an operator and then
 * "inverts" it as requested.
 *
 * NOTE: You need to set the probing Graph or GraphOperator for probing to know what to do
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
  LinearOp buildPreconditionerOperator(LinearOp& lo, PreconditionerState& state) const;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList& pl);

  void setGraphOperator(const Teko::LinearOp& graphOp);
  void setGraph(const Teuchos::RCP<const Epetra_CrsGraph>& graph);

  void setProberList(const Teuchos::ParameterList& list);

  void setInverseFactory(const Teuchos::RCP<Teko::InverseFactory>& invFactory) {
    invFactory_ = invFactory;
  }

 protected:
  //! some members
  Teuchos::RCP<Isorropia::Epetra::Prober> prober;
  Teuchos::RCP<Teko::InverseFactory> invFactory_;
};

}  // end namespace Teko

#endif
#endif
