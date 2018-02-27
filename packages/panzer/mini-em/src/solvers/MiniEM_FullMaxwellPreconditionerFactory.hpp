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

   // Use refMaxwell for Schur complement solve
   bool use_refmaxwell;

   // Dump matrices? 
   bool dump;

   // parameters
   Teuchos::ParameterList params;
};

}

#endif
