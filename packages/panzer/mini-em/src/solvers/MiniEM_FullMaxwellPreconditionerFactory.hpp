// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_FullMaxwellPreconditionerFactory_hpp_
#define _MiniEM_FullMaxwellPreconditionerFactory_hpp_

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"

namespace mini_em {

class FullMaxwellPreconditionerFactory : public Teko::BlockPreconditionerFactory {
public:
  
   FullMaxwellPreconditionerFactory() {}
   virtual ~FullMaxwellPreconditionerFactory() {}

   Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & state) const;

   //! Initialize from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);

private: 

   // Holds all inverse factories
   Teko::InverseLibrary invLib;

   bool use_discrete_curl_;
   bool simplifyFaraday_;
   bool dump;
   bool doDebug;
   bool useAsPreconditioner;
   double dt;

   // type of preconditioner for Schur complement
   std::string S_E_prec_type_;

   mutable Teko::InverseLinearOp S_E_prec_;

   // parameters
   Teuchos::ParameterList params;
};

}

#endif
