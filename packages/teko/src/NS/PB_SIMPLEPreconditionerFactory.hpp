#ifndef __PB_SIMPLEPreconditionerFactory_hpp__
#define __PB_SIMPLEPreconditionerFactory_hpp__

#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_InverseFactory.hpp"

namespace PB {
namespace NS {

// Declaration of the preconditioner factory
class SIMPLEPreconditionerFactory : public BlockPreconditionerFactory {
public:
   // Constructor
   SIMPLEPreconditionerFactory(const Teuchos::RCP<const InverseFactory> & inverse,
                               double alpha);

   // Function inherited from BlockPreconditionerFactory
   LinearOp buildPreconditionerOperator(BlockedLinearOp & blo,
                                        BlockPreconditionerState & state) const;
    
protected:
   // class members
   Teuchos::RCP<const InverseFactory> inverse_;
   double alpha_;
   
};
 
} // end namespace NS
} // end namespace PB

#endif
