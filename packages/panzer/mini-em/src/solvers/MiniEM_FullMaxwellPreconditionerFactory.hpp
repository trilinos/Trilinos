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

/** \brief A Teko block preconditioner factory for the full (E,B)
  * Maxwell system, using a block-triangular factorization with a
  * Schur complement approximation for the E-field block (see
  * S_E_prec_type_).
  */
class FullMaxwellPreconditionerFactory : public Teko::BlockPreconditionerFactory {
public:

   FullMaxwellPreconditionerFactory() {}
   virtual ~FullMaxwellPreconditionerFactory() {}

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
