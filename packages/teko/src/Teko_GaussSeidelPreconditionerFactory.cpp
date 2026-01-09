// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_GaussSeidelPreconditionerFactory.hpp"

#include "Teko_BlockUpperTriInverseOp.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teuchos_CompilerCodeTweakMacros.hpp"

using Teuchos::rcp;
using Teuchos::RCP;

namespace Teko {

GaussSeidelPreconditionerFactory::GaussSeidelPreconditionerFactory(TriSolveType solveType,
                                                                   const LinearOp& invD0,
                                                                   const LinearOp& invD1)
    : invOpsStrategy_(rcp(new StaticInvDiagStrategy(invD0, invD1))), solveType_(solveType) {}

GaussSeidelPreconditionerFactory::GaussSeidelPreconditionerFactory(
    TriSolveType solveType, const RCP<const BlockInvDiagonalStrategy>& strategy)
    : invOpsStrategy_(strategy), solveType_(solveType) {}

GaussSeidelPreconditionerFactory::GaussSeidelPreconditionerFactory()
    : solveType_(GS_UseLowerTriangle) {}

LinearOp GaussSeidelPreconditionerFactory::buildPreconditionerOperator(
    BlockedLinearOp& blo, BlockPreconditionerState& state) const {
  int rows = blockRowCount(blo);
  int cols = blockColCount(blo);

  TEUCHOS_ASSERT(rows == cols);

  // get diagonal blocks
  std::vector<LinearOp> invDiag;
  invOpsStrategy_->getInvD(blo, state, invDiag);
  TEUCHOS_ASSERT(rows == (int)invDiag.size());

  if (solveType_ == GS_UseUpperTriangle) {
    // create a blocked linear operator
    BlockedLinearOp U = getUpperTriBlocks(blo);

    return createBlockUpperTriInverseOp(U, invDiag, "Gauss Seidel");
  } else if (solveType_ == GS_UseLowerTriangle) {
    // create a blocked linear operator
    BlockedLinearOp L = getLowerTriBlocks(blo);

    return createBlockLowerTriInverseOp(L, invDiag, "Gauss Seidel");
  }

  TEUCHOS_ASSERT(false);  // we should never make it here!
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

//! Initialize from a parameter list
void GaussSeidelPreconditionerFactory::initializeFromParameterList(
    const Teuchos::ParameterList& pl) {
  Teko_DEBUG_SCOPE("GaussSeidelPreconditionerFactory::initializeFromParameterList", 10);
  Teko_DEBUG_MSG_BEGIN(9);
  DEBUG_STREAM << "Parameter list: " << std::endl;
  pl.print(DEBUG_STREAM);
  Teko_DEBUG_MSG_END();

  const std::string inverse_type        = "Inverse Type";
  const std::string strategy_name       = "Strategy Name";
  const std::string preconditioner_type = "Preconditioner Type";
  std::vector<RCP<InverseFactory> > inverses;
  std::vector<RCP<InverseFactory> > preconditioners;

  RCP<const InverseLibrary> invLib = getInverseLibrary();

  // get string specifying default inverse
  std::string invStr = "";
#if defined(Teko_ENABLE_Amesos)
  invStr = "Amesos";
#elif defined(Teko_ENABLE_Amesos2)
  invStr = "Amesos2";
#endif
  std::string precStr = "None";
  if (pl.isParameter(inverse_type)) invStr = pl.get<std::string>(inverse_type);
  if (pl.isParameter(preconditioner_type)) precStr = pl.get<std::string>(preconditioner_type);
  if (pl.isParameter("Use Upper Triangle"))
    solveType_ = pl.get<bool>("Use Upper Triangle") ? GS_UseUpperTriangle : GS_UseLowerTriangle;

  std::string strategyStr = "Jacobi Strategy";
  if (pl.isParameter(strategy_name)) strategyStr = pl.get<std::string>(strategy_name);

  Teko_DEBUG_MSG("GSPrecFact: Building default inverse \"" << invStr << "\"", 5);
  RCP<InverseFactory> defaultInverse = invLib->getInverseFactory(invStr);
  RCP<InverseFactory> defaultPrec;
  if (precStr != "None") defaultPrec = invLib->getInverseFactory(precStr);

  // now check individual solvers
  Teuchos::ParameterList::ConstIterator itr;
  for (itr = pl.begin(); itr != pl.end(); ++itr) {
    std::string fieldName = itr->first;
    Teko_DEBUG_MSG("GSPrecFact: checking fieldName = \"" << fieldName << "\"", 9);

    // figure out what the integer is
    if (fieldName.compare(0, inverse_type.length(), inverse_type) == 0 &&
        fieldName != inverse_type) {
      int position = -1;
      std::string inverse, type;

      // figure out position
      std::stringstream ss(fieldName);
      ss >> inverse >> type >> position;

      if (position <= 0) {
        Teko_DEBUG_MSG("Gauss-Seidel \"Inverse Type\" must be a (strictly) positive integer", 1);
      }

      // inserting inverse factory into vector
      std::string invStr2 = pl.get<std::string>(fieldName);
      Teko_DEBUG_MSG("GSPrecFact: Building inverse " << position << " \"" << invStr2 << "\"", 5);
      if (position > (int)inverses.size()) {
        inverses.resize(position, defaultInverse);
        inverses[position - 1] = invLib->getInverseFactory(invStr2);
      } else
        inverses[position - 1] = invLib->getInverseFactory(invStr2);
    } else if (fieldName.compare(0, preconditioner_type.length(), preconditioner_type) == 0 &&
               fieldName != preconditioner_type) {
      int position = -1;
      std::string preconditioner, type;

      // figure out position
      std::stringstream ss(fieldName);
      ss >> preconditioner >> type >> position;

      if (position <= 0) {
        Teko_DEBUG_MSG("Gauss-Seidel \"Preconditioner Type\" must be a (strictly) positive integer",
                       1);
      }

      // inserting preconditioner factory into vector
      std::string precStr2 = pl.get<std::string>(fieldName);
      Teko_DEBUG_MSG(
          "GSPrecFact: Building preconditioner " << position << " \"" << precStr2 << "\"", 5);
      if (position > (int)preconditioners.size()) {
        preconditioners.resize(position, defaultPrec);
        preconditioners[position - 1] = invLib->getInverseFactory(precStr2);
      } else
        preconditioners[position - 1] = invLib->getInverseFactory(precStr2);
    }
  }

  // use default inverse
  if (inverses.size() == 0) inverses.push_back(defaultInverse);

  // based on parameter type build a strategy
  invOpsStrategy_ =
      buildStrategy(strategyStr, inverses, preconditioners, defaultInverse, defaultPrec);
}

//! for creating the preconditioner factories objects
CloneFactory<BlockInvDiagonalStrategy> GaussSeidelPreconditionerFactory::strategyBuilder_;

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
RCP<BlockInvDiagonalStrategy> GaussSeidelPreconditionerFactory::buildStrategy(
    const std::string& name, const std::vector<Teuchos::RCP<InverseFactory> >& inverseFactories,
    const std::vector<Teuchos::RCP<InverseFactory> >& preconditionerFactories,
    const Teuchos::RCP<InverseFactory>& defaultInverseFact,
    const Teuchos::RCP<InverseFactory>& defaultPreconditionerFact) {
  Teko_DEBUG_SCOPE("JacobiPreconditionerFactory::buildStrategy", 0);

  // initialize the defaults if necessary
  if (strategyBuilder_.cloneCount() == 0) initializeStrategyBuilder();

  Teko_DEBUG_MSG_BEGIN(1) std::vector<std::string> names;
  strategyBuilder_.getCloneNames(names);
  DEBUG_STREAM << "Strategy names = ";
  for (std::size_t i = 0; i < names.size(); i++) DEBUG_STREAM << names[i] << ", ";
  DEBUG_STREAM << std::endl;
  Teko_DEBUG_MSG_END()

      // request the preconditioner factory from the CloneFactory
      RCP<BlockInvDiagonalStrategy>
          strategy = strategyBuilder_.build(name);

  if (strategy == Teuchos::null) {
    Teko_DEBUG_MSG("Warning: Could not build BlockInvDiagonalStrategy named \""
                       << name << "\"...pressing on, failure expected",
                   0) return Teuchos::null;
  }

  // now that inverse library has been set, pass in the parameter list
  strategy->initialize(inverseFactories, preconditionerFactories, defaultInverseFact,
                       defaultPreconditionerFact);

  return strategy;
}

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
void GaussSeidelPreconditionerFactory::addStrategy(const std::string& name,
                                                   const RCP<Cloneable>& clone) {
  Teko_DEBUG_SCOPE("JacobiPreconditionerFactory::addStrategy", 10);

  // initialize the defaults if necessary
  if (strategyBuilder_.cloneCount() == 0) initializeStrategyBuilder();

  // add clone to builder
  strategyBuilder_.addClone(name, clone);
}

//! This is where the default objects are put into the strategyBuilder_
void GaussSeidelPreconditionerFactory::initializeStrategyBuilder() {
  Teko_DEBUG_SCOPE("JacobiPreconditionerFactory::initializeStrategyBuilder", 10);

  RCP<Cloneable> clone;

  // add various strategies to the factory
  clone = rcp(new AutoClone<InvFactoryDiagStrategy>());
  strategyBuilder_.addClone("Jacobi Strategy", clone);
}

}  // namespace Teko
