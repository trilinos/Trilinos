// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_GaussSeidelPreconditionerFactory_hpp__
#define __Teko_GaussSeidelPreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_BlockInvDiagonalStrategy.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {

typedef enum { GS_UseLowerTriangle, GS_UseUpperTriangle } TriSolveType;

/** \brief A factory that creates a block Gauss Seidel preconditioner.
  *        The user must specify the solvers (or preconditioners) to use
  *        to approximately invert the diagonal operators.
  *
  * A factory that creates a block Gauss Seidel preconditioner.
  * The user must specify the solvers (or preconditioners) to use
  * to approximately invert the diagonal operators.
  *
  * To invoke this preconditioner using the XML file a diagonal inverse
  * needs to be specified. For example the following XML code creates
  * a Gauss-Seidel preconditioner called "GS-Outer" using Amesos
  * (a direct solver) to invert the diagonal blocks. This will invert the
  * lower triangular portion of the matrix.
  *
    \verbatim
    <ParameterList name="GS-Outer">
       <Parameter name="Type" type="string" value="Block Gauss-Seidel"/>
       <Parameter name="Use Upper Triangle" type="bool" value="false"/>
       <Parameter name="Inverse Type" type="string" value="Amesos"/>
    </ParameterList>
    \endverbatim
  *
  * Or if you want to specify a different inverse factory for a particular
  * diagonal you can use
  *
    \verbatim
    <ParameterList name="GS-Outer">
       <Parameter name="Type" type="string" value="Block Gauss-Seidel"/>
       <Parameter name="Use Upper Triangle" type="bool" value="false"/>
       <Parameter name="Inverse Type" type="string" value="DefaultInverse"/>
       <Parameter name="Inverse Type 1" type="string" value="InverseOfFirstDigonalEntry"/>
       <Parameter name="Inverse Type 3" type="string" value="InverseOfThirdDigonalEntry"/>
    </ParameterList>
    \endverbatim
  *
  * Notice that the "Inverse Type" parameter is now a default, and that you can
  * specify each diagonal inverse on its own. The diagonal entries run from 1...N where
  * N is the number of block rows. So the solver "InverseOfFirstDiagonalEntry" will
  * be used for the first diagonal block, for the second "DefaultInverse" will be used,
  * for the third "InverseOfThirdDigonalEntry" will be used, and for any further diagonal
  * blocks "DefaultInverse" will be used.
  */
class GaussSeidelPreconditionerFactory : public BlockPreconditionerFactory {
 public:
  //! @name Constructors.
  //@{

  /*! Construct a PreconditionerFactory assuming a specific block
      \f$2\times2\f$ matrix. This case is a simple one.
  */
  GaussSeidelPreconditionerFactory(TriSolveType solveType, const LinearOp& invD0,
                                   const LinearOp& invD1);

  /*! The most flexible JacobiPreconditionerFactory constructor.
      Pass in a generally defined BlockInvDiagonalStrategy to use the
      full generality of this class.
  */
  GaussSeidelPreconditionerFactory(TriSolveType solveType,
                                   const RCP<const BlockInvDiagonalStrategy>& strategy);

  /** Build an empty Gauss-Seidel preconditioner factory
   */
  GaussSeidelPreconditionerFactory();

  //@}

  /** \brief Create the Gauss-Seidel preconditioner operator.
   *
   * This method breaks apart the BlockLinearOp and builds a block
   * diagonal preconditioner. The inverse of the diagonals are specified
   * by the BlockInvDiagonalStrategy object.
   */
  LinearOp buildPreconditionerOperator(BlockedLinearOp& blo, BlockPreconditionerState& state) const;

 protected:
  using BlockPreconditionerFactory::buildPreconditionerOperator;

  //! some members
  Teuchos::RCP<const BlockInvDiagonalStrategy> invOpsStrategy_;
  TriSolveType solveType_;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList& pl);

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
        const std::string& name,
        const std::vector<Teuchos::RCP<InverseFactory> > &inverseFactories,
        const std::vector<Teuchos::RCP<InverseFactory> > &preconditionerFactories,
        const Teuchos::RCP<InverseFactory> &defaultInverseFact,
        const Teuchos::RCP<InverseFactory> &defaultPreconditionerFact);

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
