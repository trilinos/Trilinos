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

 public:
  /** \brief Builder function for creating strategies.
   *
   * Builder function for creating strategies.
   *
   * \param[in] name     String name of strategy to build
   * \param[in] settings Parameter list describing the parameters for the
   *                     strategy to build
   * \param[in] invLib   Inverse library for the strategy to use.
   *
   * \returns If the name is associated with a strategy
   *          a pointer is returned, otherwise Teuchos::null is returned.
   */
  static RCP<BlockInvDiagonalStrategy> buildStrategy(
      const std::string& name, const std::vector<Teuchos::RCP<InverseFactory> >& inverseFactories,
      const std::vector<Teuchos::RCP<InverseFactory> >& preconditionerFactories,
      const Teuchos::RCP<InverseFactory>& defaultInverseFact,
      const Teuchos::RCP<InverseFactory>& defaultPreconditionerFact);

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList& pl);

  /** \brief Add a strategy to the builder. This is done using the
   *        clone pattern.
   *
   * Add a strategy to the builder. This is done using the
   * clone pattern. If your class does not support the Cloneable interface then
   * you can use the AutoClone class to construct your object.
   *
   * \note If this method is called twice with the same string, the latter clone pointer
   *       will be used.
   *
   * \param[in] name String to associate with this object
   * \param[in] clone Pointer to Cloneable object
   */
  static void addStrategy(const std::string& name, const RCP<Cloneable>& clone);

 private:
  //! for creating the strategy objects
  static CloneFactory<BlockInvDiagonalStrategy> strategyBuilder_;

  //! This is where the default objects are put into the strategyBuilder_
  static void initializeStrategyBuilder();
};

}  // end namespace Teko

#endif
