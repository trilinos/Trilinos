#ifndef __PB_AddPreconditionerFactory_hpp__
#define __PB_AddPreconditionerFactory_hpp__

#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_Utilities.hpp"

namespace PB {

/** Preconditioning factories must supply a 'State' class which
  * is where data specific to the preconditioner construction 
  * is stored. The constructor will be invoked elsewhere.
  */
class AddPrecondState : public PB::BlockPreconditionerState {
public:
   AddPrecondState() {}

   Teuchos::RCP<BlockPreconditionerState> StateOne_;
   Teuchos::RCP<BlockPreconditionerState> StateTwo_;
};

/** Declaration of preconditioner factory that creates
  * a preconditioner which is the sum (additive) of two
  * other preconditioners.
  */
class AddPreconditionerFactory 
   : public PB::BlockPreconditionerFactory {
public:
   //! Constructor
   AddPreconditionerFactory(const Teuchos::RCP<const PB::BlockPreconditionerFactory> & FirstFactory,
                            const Teuchos::RCP<const PB::BlockPreconditionerFactory> & SecondFactory);

   AddPreconditionerFactory();

   //! Function inherited from PB::BlockPreconditionerFactory
   PB::LinearOp buildPreconditionerOperator(PB::BlockedLinearOp & blo,
                                            PB::BlockPreconditionerState & state) const;
    
   //! Build the AddPrecondState object
   virtual Teuchos::RCP<PB::BlockPreconditionerState> buildPreconditionerState() const;

protected:
   // class members
   Teuchos::RCP<const PB::BlockPreconditionerFactory> FirstFactory_;
   Teuchos::RCP<const PB::BlockPreconditionerFactory> SecondFactory_;

   //! Initialize from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);
};

} // end namespace PB

#endif
