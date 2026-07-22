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

class FullDarcyPreconditionerFactory : public Teko::BlockPreconditionerFactory {
public:

   FullDarcyPreconditionerFactory() {}
   virtual ~FullDarcyPreconditionerFactory() {}

   Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & state) const;

   //! Initialize from a parameter list
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
