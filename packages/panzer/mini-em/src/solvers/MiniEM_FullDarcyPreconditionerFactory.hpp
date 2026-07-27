// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_FullDarcyPreconditionerFactory_hpp_
#define _MiniEM_FullDarcyPreconditionerFactory_hpp_

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"

namespace mini_em {

/** \brief A Teko block preconditioner factory for the full (u,sigma)
  * mixed Darcy system, using a block-triangular factorization with a
  * Schur complement approximation for the sigma (flux) block (see
  * S_sigma_prec_type_).
  */
class FullDarcyPreconditionerFactory : public Teko::BlockPreconditionerFactory {
public:

   FullDarcyPreconditionerFactory() {}
   virtual ~FullDarcyPreconditionerFactory() {}

   /** \brief Builds the block preconditioner operator for the given blocked Darcy system, per the class description.
     *
     * \param[in] blo the blocked (u,sigma) system operator to precondition.
     * \param[in,out] state Teko's preconditioner state, used to cache sub-block operators/inverses across applications.
     */
   Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & state) const;

   //! Initialize the sub-solver inverse factories and options from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);

private:

   // Holds all inverse factories
   Teko::InverseLibrary invLib;

   bool use_discrete_div_;
   bool simplifyFaraday_;
   bool dump;
   bool doDebug;
   bool useAsPreconditioner;
   double dt;
   std::vector<int> pCoarsenSchedule_;

   // type of preconditioner for Schur complement
   std::string S_sigma_prec_type_;

   mutable Teko::InverseLinearOp S_sigma_prec_;

   // parameters
   Teuchos::ParameterList params;
};

}

#endif
