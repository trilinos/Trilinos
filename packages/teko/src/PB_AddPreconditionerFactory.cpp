#include "PB_AddPreconditionerFactory.hpp"

namespace PB {

using Teuchos::RCP;

AddPreconditionerFactory::AddPreconditionerFactory(
                            const RCP<const BlockPreconditionerFactory> & FirstFactory,
                            const RCP<const BlockPreconditionerFactory> & SecondFactory)
   : FirstFactory_(FirstFactory), SecondFactory_(SecondFactory)
{}

//! Build the AddPrecondState object
RCP<BlockPreconditionerState> AddPreconditionerFactory::buildPreconditionerState() const
{ 
   AddPrecondState*   mystate = new AddPrecondState(); 
   mystate->StateOne_ = FirstFactory_->buildPreconditionerState();
   mystate->StateTwo_ = SecondFactory_->buildPreconditionerState();
   return rcp(mystate);
}

// Use the factory to build the preconditioner (this is where the work goes)
LinearOp AddPreconditionerFactory 
   ::buildPreconditionerOperator(BlockedLinearOp & blockOp,
                                 BlockPreconditionerState & state) const
{
   // The main tricky thing here is that we have to take the 'state' object
   // associated with AddPreconditionerFactory(), pull out the states for
   // the individual preconditioners, and pass these on to 
   // buildPreconditionerOperator() for each subpreconditioner.
   
   AddPrecondState *MyState = dynamic_cast<AddPrecondState *> (&state);
   TEUCHOS_ASSERT(MyState != 0);

   LinearOp M1 = FirstFactory_->buildPreconditionerOperator(blockOp, *MyState->StateOne_);
   LinearOp M2 = SecondFactory_->buildPreconditionerOperator(blockOp, *MyState->StateTwo_);

   LinearOp invA = add(M1, M2);

   // return fully constructed preconditioner
   return invA;
}

} // end namespace PB
