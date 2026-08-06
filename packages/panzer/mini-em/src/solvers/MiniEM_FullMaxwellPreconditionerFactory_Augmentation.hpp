// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_FullMaxwellPreconditionerFactory_Augmentation_hpp_
#define _MiniEM_FullMaxwellPreconditionerFactory_Augmentation_hpp_

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"

namespace mini_em {

/** \brief A Teko block preconditioner factory for the full (E,B)
  * Maxwell system, using the augmentation approach: an auxiliary
  * augmented (elliptic-like) operator built from the discrete
  * gradient is used to precondition the otherwise near-singular
  * curl-curl block.
  */
class FullMaxwellPreconditionerFactory_Augmentation : public Teko::BlockPreconditionerFactory {
public:

   FullMaxwellPreconditionerFactory_Augmentation() {}
   virtual ~FullMaxwellPreconditionerFactory_Augmentation() {}

   /** \brief Builds the block preconditioner operator for the given blocked Maxwell system, per the class description.
     *
     * \param[in] blo the blocked (E,B) system operator to precondition.
     * \param[in,out] state Teko's preconditioner state, used to cache sub-block operators/inverses across applications.
     */
   Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & state) const;

   //! Initialize the sub-solver inverse factories and options from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);

private:

   // Holds all inverse factories
   Teko::InverseLibrary invLib;

   bool use_discrete_gradient_;
   bool dump;
   bool doDebug;
   bool useAsPreconditioner;

   // parameters
   Teuchos::ParameterList params;
};

}

#endif
