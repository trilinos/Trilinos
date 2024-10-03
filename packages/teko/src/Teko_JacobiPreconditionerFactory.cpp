// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_JacobiPreconditionerFactory.hpp"
#include "Teko_BlockDiagonalInverseOp.hpp"

using Teuchos::rcp;

namespace Teko {

JacobiPreconditionerFactory::JacobiPreconditionerFactory(const LinearOp& invD0,
                                                         const LinearOp& invD1)
    : invOpsStrategy_(rcp(new StaticInvDiagStrategy(invD0, invD1))) {}

JacobiPreconditionerFactory::JacobiPreconditionerFactory(
    const RCP<const BlockInvDiagonalStrategy>& strategy)
    : invOpsStrategy_(strategy) {}

/** Build a Jacobi preconditioner factory from a parameter list
 */
JacobiPreconditionerFactory::JacobiPreconditionerFactory() {}

LinearOp JacobiPreconditionerFactory::buildPreconditionerOperator(
    BlockedLinearOp& blo, BlockPreconditionerState& state) const {
  int rows = blo->productRange()->numBlocks();
  int cols = blo->productDomain()->numBlocks();

  TEUCHOS_ASSERT(rows == cols);

  // get diagonal blocks
  std::vector<LinearOp> invDiag;
  invOpsStrategy_->getInvD(blo, state, invDiag);
  TEUCHOS_ASSERT(rows == (int)invDiag.size());

  return createBlockDiagonalInverseOp(blo, invDiag, "Jacobi");
}

//! Initialize from a parameter list
void JacobiPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList& pl)
#if 0
{
   RCP<const InverseLibrary> invLib = getInverseLibrary();

   // get string specifying inverse
   std::string invStr = pl.get<std::string>("Inverse Type");
#if defined(Teko_ENABLE_Amesos)
   if(invStr=="") invStr = "Amesos";
#elif defined(Teko_ENABLE_Amesos2)
   if(invStr=="") invStr= "Amesos2";
#endif

   // based on parameter type build a strategy
   invOpsStrategy_ = rcp(new InvFactoryDiagStrategy(invLib->getInverseFactory(invStr)));
}
#endif
{
  Teko_DEBUG_SCOPE("JacobiPreconditionerFactory::initializeFromParameterList", 10);
  Teko_DEBUG_MSG_BEGIN(9);
  DEBUG_STREAM << "Parameter list: " << std::endl;
  pl.print(DEBUG_STREAM);
  Teko_DEBUG_MSG_END();

  const std::string inverse_type        = "Inverse Type";
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
  RCP<InverseFactory> defaultPrec;
  if (precStr != "None") defaultPrec = invLib->getInverseFactory(precStr);

  Teko_DEBUG_MSG("JacobiPrecFact: Building default inverse \"" << invStr << "\"", 5);
  RCP<InverseFactory> defaultInverse = invLib->getInverseFactory(invStr);

  // now check individual solvers
  Teuchos::ParameterList::ConstIterator itr;
  for (itr = pl.begin(); itr != pl.end(); ++itr) {
    std::string fieldName = itr->first;
    Teko_DEBUG_MSG("JacobiPrecFact: checking fieldName = \"" << fieldName << "\"", 9);

    // figure out what the integer is
    if (fieldName.compare(0, inverse_type.length(), inverse_type) == 0 &&
        fieldName != inverse_type) {
      int position = -1;
      std::string inverse, type;

      // figure out position
      std::stringstream ss(fieldName);
      ss >> inverse >> type >> position;

      if (position <= 0) {
        Teko_DEBUG_MSG("Jacobi \"Inverse Type\" must be a (strictly) positive integer", 1);
      }

      // inserting inverse factory into vector
      std::string invStr2 = pl.get<std::string>(fieldName);
      Teko_DEBUG_MSG("JacobiPrecFact: Building inverse " << position << " \"" << invStr2 << "\"",
                     5);
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
        Teko_DEBUG_MSG("Jacobi \"Preconditioner Type\" must be a (strictly) positive integer", 1);
      }

      // inserting preconditioner factory into vector
      std::string precStr2 = pl.get<std::string>(fieldName);
      Teko_DEBUG_MSG(
          "JacobiPrecFact: Building preconditioner " << position << " \"" << precStr2 << "\"", 5);
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
      rcp(new InvFactoryDiagStrategy(inverses, preconditioners, defaultInverse, defaultPrec));
}

}  // namespace Teko
