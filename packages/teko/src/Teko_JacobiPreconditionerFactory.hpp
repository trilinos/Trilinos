// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_JacobiPreconditionerFactory_hpp__
#define __Teko_JacobiPreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_BlockInvDiagonalStrategy.hpp"

namespace Teko {

/** A block diagonal preconditioner.
  *
  * The use must specify an iverse for each diagonal.
  * If a specific integer is not specified, then the default
  * "Inverse Type" is used.
  \code
       <Parameter name="Type" type="string" value="Block Jacobi"/>
       <Parameter name="Inverse Type" type="string" value="<Some Inverse Factory>"/>
       <Parameter name="Inverse Type 1" type="string" value="<Some Inverse Factory>"/>
       <Parameter name="Inverse Type 2" type="string" value="<Some Inverse Factory>"/>
       <Parameter name="Inverse Type 3" type="string" value="<Some Inverse Factory>"/>
       <Parameter name="Inverse Type 4" type="string" value="<Some Inverse Factory>"/>
       <Parameter name="Inverse Type 5" type="string" value="<Some Inverse Factory>"/>
  \endcode
  */
class JacobiPreconditionerFactory : public BlockPreconditionerFactory {
 public:
  //! @name Constructors.
  //@{

  /** Construct a PreconditionerFactory assuming a specific block
   * \f$2\times2\f$ matrix. This case is a simple one.
   */
  JacobiPreconditionerFactory(const LinearOp& invD0, const LinearOp& invD1);

  /** The most flexible JacobiPreconditionerFactory constructor.
   * Pass in a generally defined BlockInvDiagonalStrategy to use the
   * full generality of this class.
   */
  JacobiPreconditionerFactory(const RCP<const BlockInvDiagonalStrategy>& strategy);

  /** Build an empty Jacobi preconditioner factory.
   */
  JacobiPreconditionerFactory();

  //@}

  /** \brief Create the Jacobi preconditioner operator.
   *
   * This method breaks apart the BlockLinearOp and builds a block
   * diagonal preconditioner. The inverse of the diagonals are specified
   * by the BlockInvDiagonalStrategy object.
   */
  LinearOp buildPreconditionerOperator(BlockedLinearOp& blo, BlockPreconditionerState& state) const;

  //! Get inv diagonal strategy
  Teuchos::RCP<const BlockInvDiagonalStrategy> getInvDiagStrategy() const {
    return invOpsStrategy_;
  }

 protected:
  using Teko::BlockPreconditionerFactory::buildPreconditionerOperator;

  //! some members
  Teuchos::RCP<const BlockInvDiagonalStrategy> invOpsStrategy_;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList& pl);
};

}  // end namespace Teko

#endif
