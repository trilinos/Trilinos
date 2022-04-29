#ifndef _MiniEM_HigherOrderMaxwellPreconditionerFactory_hpp_
#define _MiniEM_HigherOrderMaxwellPreconditionerFactory_hpp_

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"

namespace mini_em {

class HigherOrderMaxwellPreconditionerFactory : public Teko::BlockPreconditionerFactory {
public:

   HigherOrderMaxwellPreconditionerFactory() {}
   virtual ~HigherOrderMaxwellPreconditionerFactory() {}

   Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & state) const;

   //! Initialize from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);

private:

   // Holds all inverse factories
   Teko::InverseLibrary invLib;

   bool simplifyFaraday_;
   bool dump;
   bool doDebug;
   bool useAsPreconditioner;
   double dt;
   std::vector<int> pCoarsenSchedule_;

   // type of preconditioner for Schur complement
   std::string S_E_prec_type_;

   mutable Teko::InverseLinearOp S_E_prec_;

   // parameters
   Teuchos::ParameterList params;
};

}

#endif
