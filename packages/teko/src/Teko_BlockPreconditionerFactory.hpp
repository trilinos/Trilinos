// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_BlockPreconditionerFactory_hpp__
#define __Teko_BlockPreconditionerFactory_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

// Thyra includes
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_CloneFactory.hpp"
#include "Teko_PreconditionerState.hpp"
#include "Teko_PreconditionerFactory.hpp"

namespace Teko {

using Teuchos::ParameterList;
using Teuchos::rcp;
using Teuchos::RCP;
using Thyra::DefaultPreconditioner;
using Thyra::LinearOpBase;

/** \brief An implementation of a state object for block
 *        preconditioners.
 *
 * An implementation of the state object for block
 * preconditioners. T
 */
class BlockPreconditionerState : public PreconditionerState {
 public:
  //! Set the vector associated with this operator (think nonlinear system)
  virtual void setBlockSourceVector(const Teko::BlockedMultiVector &srcVec) {
    setSourceVector(srcVec);
  }

  //! Set the vector associated with this operator (think nonlinear system)
  virtual const Teko::BlockedMultiVector getBlockedSourceVector() const {
    return toBlockedMultiVector(getSourceVector());
  }
};

/** \brief Abstract class which block preconditioner factories in Teko
 *        should be based on.
 *
 * Abstract class which block preconditioner factories in Teko should
 * be based on. All that is needed is the implementation of
 * "buildPreconditionerOperator".
 */
class BlockPreconditionerFactory : public PreconditionerFactory {
 public:
  /** \brief Function that is called to build the preconditioner
   *        for the linear operator that is passed in.
   *
   * This function builds a preconditioner based on the passed
   * in BlockedLinearOp.
   *
   * \param[in] blo   Source linear operator that is to be preconditioned.
   * \param[in] state An object associated with this operator to store
   *                  the preconditioner state.
   *
   * \returns The preconditioner as a linear operator (i.e. to perform
   *           a matrix-vector operation simply call "apply").
   */
  virtual LinearOp buildPreconditionerOperator(BlockedLinearOp &blo,
                                               BlockPreconditionerState &state) const = 0;

  /** \brief Function that permits the construction of an arbitrary
   *        BlockPreconditionerState object.
   *
   * Function that permits the construction of an arbitrary
   * BlockPreconditionerState object. If the basic state object,
   * which takes a parameter list, is sufficient the default behavior
   * does precisely what is needed. Otherwise, an author of a
   * PreconditionerFactory would need to reimplement this method to
   * return a new state object.
   *
   * \returns A state object associated with this factory.
   */
  virtual RCP<PreconditionerState> buildPreconditionerState() const {
    return rcp(new BlockPreconditionerState());
  }

  /** \brief Function that constructs a BlockPreconditionerState object.
   *
   * This function is simply included for convenience.  Its implementation
   * just call "buildPreconditionerState" and returns a dynamically casted
   * <code>BlockPreconditionState</code> poniter.
   */
  RCP<BlockPreconditionerState> buildBlockPreconditionerState() const {
    return Teuchos::rcp_dynamic_cast<BlockPreconditionerState>(buildPreconditionerState());
  }

  virtual LinearOp buildPreconditionerOperator(LinearOp &blo, PreconditionerState &state) const;

  //! is this operator compatiable with the preconditioner factory?
  bool isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc) const;
};

}  // end namespace Teko

#endif
