/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

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
