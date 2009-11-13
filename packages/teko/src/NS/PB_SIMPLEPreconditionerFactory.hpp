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
   SIMPLEPreconditionerFactory(const Teuchos::RCP<InverseFactory> & inverse,
                               double alpha);

   // Constructor
   SIMPLEPreconditionerFactory(const Teuchos::RCP<InverseFactory> & invVelFactory,
                               const Teuchos::RCP<InverseFactory> & invPrsFactory,
                               double alpha);

   //! Default constructor
   SIMPLEPreconditionerFactory();

   // Function inherited from BlockPreconditionerFactory
   LinearOp buildPreconditionerOperator(BlockedLinearOp & blo,
                                        BlockPreconditionerState & state) const;
    
protected:
   // class members
   Teuchos::RCP<InverseFactory> invVelFactory_;
   Teuchos::RCP<InverseFactory> invPrsFactory_;
   double alpha_;
   enum FInverseType {Diagonal,Lumped,AbsRowSum,Custom} fInverseType_;
   
   //! Initialize from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);
};
 
} // end namespace NS
} // end namespace PB

#endif
