#include "PB_LU2x2DiagonalStrategy.hpp"

namespace PB {

LU2x2DiagonalStrategy::LU2x2DiagonalStrategy() 
{ }

//! Constructor to set the inverse factories.
LU2x2DiagonalStrategy::LU2x2DiagonalStrategy(const Teuchos::RCP<InverseFactory> & invFA,
                                             const Teuchos::RCP<InverseFactory> & invS)
   : invFactoryA00_(invFA), invFactoryS_(invS) 
{ }

/** returns the first (approximate) inverse of \f$A_{00}\f$ */
const PB::LinearOp
LU2x2DiagonalStrategy::getHatInvA00(const PB::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   initializeState(A,state);

   return state.getInverse("invA00");
}

/** returns the second (approximate) inverse of \f$A_{00}\f$ */
const PB::LinearOp
LU2x2DiagonalStrategy::getTildeInvA00(const PB::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   initializeState(A,state);

   return state.getInverse("invA00");
}

/** returns an (approximate) inverse of \f$S = -A_{11} + A_{10} \mbox{diag}(A_{00})^{-1} A_{01}\f$ */
const PB::LinearOp
LU2x2DiagonalStrategy::getInvS(const PB::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   initializeState(A,state);

   return state.getInverse("invS");
}

void LU2x2DiagonalStrategy::initializeState(const PB::BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   PB_DEBUG_SCOPE("LU2x2DiagonalStrategy::initializeState",10);

   // no work to be done
   if(state.isInitialized())
      return;

   // extract sub blocks
   LinearOp A00 = PB::getBlock(0,0,A);
   LinearOp A01 = PB::getBlock(0,1,A);
   LinearOp A10 = PB::getBlock(1,0,A);
   LinearOp A11 = PB::getBlock(1,1,A);

   // build the Schur complement
   /////////////////////////////////////////////
   PB_DEBUG_MSG("Building S",5);
   LinearOp diagA00 = getInvDiagonalOp(A00);

   // grab operators for building Schur complement 
   InverseLinearOp triple = state.getInverse("triple");
   InverseLinearOp S = state.getInverse("S");

   // build Schur-complement
   triple = explicitMultiply(A10,diagA00,A01,triple);
   S = explicitAdd(scale(-1.0,A11),triple,S);

   // add to state
   state.addInverse("trip",triple);
   state.addInverse("S",S);

   // build inverse S
   /////////////////////////////////////////////
   PB_DEBUG_MSG("Building inverse(S)",5);
   InverseLinearOp invS = state.getInverse("invS");
   if(invS==Teuchos::null)
      invS = buildInverse(*invFactoryS_,S);
   else
      rebuildInverse(*invFactoryS_,S,invS);
   state.addInverse("invS",invS);

   // build inverse A00
   /////////////////////////////////////////////
   PB_DEBUG_MSG("Building inverse(A00)",5);
   InverseLinearOp invA00 = state.getInverse("invA00");
   if(invA00==Teuchos::null)
      invA00 = buildInverse(*invFactoryA00_,A00);
   else
      rebuildInverse(*invFactoryA00_,A00,invA00);
   state.addInverse("invA00",invA00);

   // mark state as initialized
   state.setInitialized(true);
}

/** \brief This function builds the internals of the state from a parameter list.
  *        
  * This function builds the internals of the LU 2x2 state
  * from a parameter list. Furthermore, it allows a 
  * developer to easily add a factory to the build system.
  *
  * \param[in] settings Parameter list to use as the internal settings
  * \param[in] invLib Inverse library to use for building inverse factory objects
  *
  * \note The default implementation does nothing.
  */
void LU2x2DiagonalStrategy::initializeFromParameterList(const Teuchos::ParameterList & pl,
                                                        const InverseLibrary & invLib)
{
   PB_DEBUG_SCOPE("LU2x2DiagonalStrategy::initializeFromParameterList",10);

   std::string invStr="Amesos", invA00Str="", invSStr="";

   // "parse" the parameter list
   if(pl.isParameter("Inverse Type"))
      invStr = pl.get<std::string>("Inverse Type");
   if(pl.isParameter("Inverse A00 Type"))
      invA00Str = pl.get<std::string>("Inverse A00 Type");
   if(pl.isParameter("Inverse Schur Type"))
      invSStr = pl.get<std::string>("Inverse Schur Type");

   // set defaults as needed
   if(invA00Str=="") invA00Str = invStr;
   if(invSStr=="") invSStr = invStr;

   PB_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "LU2x2 Diagonal Strategy Parameters: " << std::endl;
      DEBUG_STREAM << "   inv type   = \"" << invStr  << "\"" << std::endl;
      DEBUG_STREAM << "   inv A00 type = \"" << invA00Str << "\"" << std::endl;
      DEBUG_STREAM << "   inv S type = \"" << invSStr << "\"" << std::endl;
      DEBUG_STREAM << "LU2x2 Diagonal Strategy Parameter list: " << std::endl;
      pl.print(DEBUG_STREAM);
   PB_DEBUG_MSG_END()

   // build velocity inverse factory
   invFactoryA00_ = invLib.getInverseFactory(invA00Str);
 
   if(invA00Str==invSStr)
      invFactoryS_ = invFactoryA00_;
   else
      invFactoryS_ = invLib.getInverseFactory(invSStr);
}
} // end namespace PB
