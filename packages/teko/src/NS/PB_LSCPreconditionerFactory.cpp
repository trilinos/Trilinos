#include "PB_LSCPreconditionerFactory.hpp"

#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "PB_LU2x2InverseOp.hpp"
#include "PB_Utilities.hpp"
#include "PB_BlockUpperTriInverseOp.hpp"
#include "PB_StaticLSCStrategy.hpp"
#include "PB_InvLSCStrategy.hpp"

#include "EpetraExt_RowMatrixOut.h"

#include "Teuchos_Time.hpp"

namespace PB {
namespace NS {

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;

using Thyra::multiply;
using Thyra::add;
using Thyra::identity;

// Stabilized constructor
LSCPreconditionerFactory::LSCPreconditionerFactory(const LinearOp & invF,const LinearOp & invBQBtmC,
                                                   const LinearOp & invD,const LinearOp & invMass)
      : invOpsStrategy_(rcp(new StaticLSCStrategy(invF,invBQBtmC,invD,invMass))), useMass_(false)
{ }

// Stable constructor
LSCPreconditionerFactory::LSCPreconditionerFactory(const LinearOp & invF, const LinearOp & invBQBtmC,
                                                   const LinearOp & invMass)
      : invOpsStrategy_(rcp(new StaticLSCStrategy(invF,invBQBtmC,invMass))), useMass_(false)
{ }

// fully generic constructor
LSCPreconditionerFactory::LSCPreconditionerFactory(const RCP<LSCStrategy> & strategy)
   : invOpsStrategy_(strategy), useMass_(false)
{ }

LSCPreconditionerFactory::LSCPreconditionerFactory() : useMass_(false)
{ }

// for PreconditionerFactoryBase
///////////////////////////////////////////////////////////////////////

// initialize a newly created preconditioner object
LinearOp LSCPreconditionerFactory::buildPreconditionerOperator(BlockedLinearOp & blockOp,BlockPreconditionerState & state) const
{
   PB_DEBUG_SCOPE("LSCPreconditionerFactory::buildPreconditionerOperator",10);
   PB_DEBUG_EXPR(Teuchos::Time timer(""));
   PB_DEBUG_EXPR(Teuchos::Time totalTimer(""));
   PB_DEBUG_EXPR(totalTimer.start());

   // extract sub-matrices from source operator 
   LinearOp F  = blockOp->getBlock(0,0);
   LinearOp B  = blockOp->getBlock(1,0);
   LinearOp Bt = blockOp->getBlock(0,1);

   // build what is neccessary for the state object
   PB_DEBUG_EXPR(timer.start(true));
   invOpsStrategy_->buildState(blockOp,state);
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("LSCPrecFact::buildPO BuildStateTime = " << timer.totalElapsedTime(),2);

   // extract operators from strategy
   PB_DEBUG_EXPR(timer.start(true));
   LinearOp invF      = invOpsStrategy_->getInvF(blockOp,state);
   LinearOp invBQBtmC = invOpsStrategy_->getInvBQBt(blockOp,state);
   LinearOp invBHBtmC = invOpsStrategy_->getInvBHBt(blockOp,state);
   LinearOp invAlphaD = invOpsStrategy_->getInvAlphaD(blockOp,state);

   // if necessary build an identity mass matrix
   LinearOp invMass   = invOpsStrategy_->getInvMass(blockOp,state);
   LinearOp HScaling  = invOpsStrategy_->getHScaling(blockOp,state);
   if(invMass==Teuchos::null)  invMass = identity<double>(F->range());
   if(HScaling==Teuchos::null) HScaling = identity<double>(F->range());
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("LSCPrecFact::buildPO GetInvTime = " << timer.totalElapsedTime(),2);

   // need to build Schur complement,  inv(P) = inv(B*Bt)*(B*F*Bt)*inv(B*Bt)

   // first construct middle operator: M = B * inv(Mass) * F * inv(Mass) * Bt
   LinearOp M = 
      //          (B * inv(Mass) ) * F * (inv(Mass) * Bt)
      multiply( multiply(B,invMass), F , multiply(HScaling,Bt));
      
   // now construct a linear operator schur complement
   LinearOp invPschur; 
   if(invAlphaD!=Teuchos::null)
      invPschur = add(multiply(invBQBtmC, M , invBHBtmC), invAlphaD);
   else
      invPschur = multiply(invBQBtmC, M , invBHBtmC);

   // build the preconditioner operator: Use LDU or upper triangular approximation
   if(invOpsStrategy_->useFullLDU()) { 
      PB_DEBUG_EXPR(totalTimer.stop());
      PB_DEBUG_MSG("LSCPrecFact::buildPO TotalTime = " << totalTimer.totalElapsedTime(),2);

      // solve using a full LDU decomposition
      return createLU2x2InverseOp(blockOp,invF,invPschur,"LSC-LDU");
   } else {
      // build diagonal operations
      std::vector<LinearOp> invDiag(2);
      invDiag[0] = invF;
      invDiag[1] = Thyra::scale(-1.0,invPschur);

      // get upper triangular matrix
      BlockedLinearOp U = getUpperTriBlocks(blockOp); 

      PB_DEBUG_EXPR(totalTimer.stop());
      PB_DEBUG_MSG("LSCPrecFact::buildPO TotalTime = " << totalTimer.totalElapsedTime(),2);

      // solve using only one inversion of F
      return createBlockUpperTriInverseOp(U,invDiag,"LSC-Upper");
   }
}

//! Initialize from a parameter list
void LSCPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   PB_DEBUG_MSG("Begin LSCPreconditionerFactory::initializeFromParameterList",10);

   RCP<const InverseLibrary> invLib = getInverseLibrary();

   std::string name = "Basic Inverse";
   if(pl.isParameter("Strategy Name"))
      name = pl.get<std::string>("Strategy Name");
   const Teuchos::ParameterEntry * pe = pl.getEntryPtr("Strategy Settings");

   // check for a mistake in input file
   if(name!="Basic Inverse" && pe==0) {
      RCP<Teuchos::FancyOStream> out = getOutputStream();
      *out << "LSC Construction failed: ";
      *out << "Strategy \"" << name << "\" requires a \"Strategy Settings\" sublist" << std::endl;
      throw std::runtime_error("LSC Construction failed: Strategy Settings not set");
   }

   // get the parameter list to construct the strategy
   Teuchos::RCP<const Teuchos::ParameterList> stratPL = Teuchos::rcpFromRef(pl);
   if(pe!=0)
      stratPL = Teuchos::rcpFromRef(pl.sublist("Strategy Settings"));

   // build the strategy object
   RCP<LSCStrategy> strategy = buildStrategy(name,*stratPL,invLib);
 
   // strategy could not be built
   if(strategy==Teuchos::null) {
      RCP<Teuchos::FancyOStream> out = getOutputStream();
      *out << "LSC Construction failed: ";
      *out << "Strategy \"" << name << "\" could not be constructed" << std::endl;
      throw std::runtime_error("LSC Construction failed: Strategy could not be constructed");
   }

   invOpsStrategy_ = strategy;

   PB_DEBUG_MSG("End LSCPreconditionerFactory::initializeFromParameterList",10);
}

//! For assiting in construction of the preconditioner
Teuchos::RCP<Teuchos::ParameterList> LSCPreconditionerFactory::getRequestedParameters() const
{
   return invOpsStrategy_->getRequestedParameters();
}

//! For assiting in construction of the preconditioner
bool LSCPreconditionerFactory::updateRequestedParameters(const Teuchos::ParameterList & pl)
{
   return invOpsStrategy_->updateRequestedParameters(pl);
}

/////////////////////////////////////////////////////
// Static members and methods
/////////////////////////////////////////////////////

//! for creating the preconditioner factories objects
CloneFactory<LSCStrategy> LSCPreconditionerFactory::strategyBuilder_;

/** \brief Builder function for creating strategies.
  *
  * Builder function for creating strategies.
  * 
  * \param[in] name     String name of strategy to build
  * \param[in] settings Parameter list describing the parameters for the
  *                     strategy to build
  * \param[in] invLib   Inverse library for the strategy to use.
  *
  * \returns If the name is associated with a strategy
  *          a pointer is returned, otherwise Teuchos::null is returned.
  */
RCP<LSCStrategy> LSCPreconditionerFactory::buildStrategy(const std::string & name, 
                                                             const Teuchos::ParameterList & settings,
                                                             const RCP<const InverseLibrary> & invLib)
{
   PB_DEBUG_MSG("Begin LSCPreconditionerFactory::buildStrategy",10);

   // initialize the defaults if necessary
   if(strategyBuilder_.cloneCount()==0) initializeStrategyBuilder();

   // request the preconditioner factory from the CloneFactory
   RCP<LSCStrategy> strategy = strategyBuilder_.build(name);

   if(strategy==Teuchos::null) return Teuchos::null;

   // now that inverse library has been set,
   // pass in the parameter list
   strategy->initializeFromParameterList(settings,*invLib);

   PB_DEBUG_MSG("End LSCPreconditionerFactory::buildStrategy",10);

   return strategy;
}

/** \brief Add a strategy to the builder. This is done using the
  *        clone pattern. 
  *
  * Add a strategy to the builder. This is done using the
  * clone pattern. If your class does not support the Cloneable interface then
  * you can use the AutoClone class to construct your object.
  *
  * \note If this method is called twice with the same string, the latter clone pointer
  *       will be used.
  *
  * \param[in] name String to associate with this object
  * \param[in] clone Pointer to Cloneable object
  */
void LSCPreconditionerFactory::addStrategy(const std::string & name,const RCP<Cloneable> & clone)
{
   // initialize the defaults if necessary
   if(strategyBuilder_.cloneCount()==0) initializeStrategyBuilder();

   // add clone to builder
   strategyBuilder_.addClone(name,clone); 
}

//! This is where the default objects are put into the strategyBuilder_
void LSCPreconditionerFactory::initializeStrategyBuilder()
{
   RCP<Cloneable> clone;

   // add various strategies to the factory
   clone = rcp(new AutoClone<InvLSCStrategy>());
   strategyBuilder_.addClone("Basic Inverse",clone);
}

} // end namespace NS
} // end namespace PB
