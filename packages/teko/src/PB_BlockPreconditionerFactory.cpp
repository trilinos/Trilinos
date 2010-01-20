#include "PB_BlockPreconditionerFactory.hpp"

#include "PB_InverseLibrary.hpp"
#include "PB_JacobiPreconditionerFactory.hpp"
#include "PB_GaussSeidelPreconditionerFactory.hpp"
#include "PB_AddPreconditionerFactory.hpp"
#include "PB_MultPreconditionerFactory.hpp"
#include "PB_LU2x2PreconditionerFactory.hpp"
#include "NS/Teko_LSCPreconditionerFactory.hpp"
#include "NS/Teko_SIMPLEPreconditionerFactory.hpp"

#include "Thyra_DefaultPreconditioner.hpp"

using namespace Thyra;

namespace Teko {
/////////////////////////////////////////////////////

//! Set parameters from a parameter list and return with default values.
void BlockPreconditionerState::setParameterList(const RCP<ParameterList> & paramList)
{
   paramList_ = paramList;
}

//! Get the parameter list that was set using setParameterList().
RCP<ParameterList> BlockPreconditionerState::getNonconstParameterList()
{
   if(paramList_==Teuchos::null)
      paramList_ = Teuchos::rcp(new Teuchos::ParameterList());

   return paramList_;
}

//! Unset the parameter list that was set using setParameterList(). 
RCP<ParameterList> BlockPreconditionerState::unsetParameterList()
{
   RCP<ParameterList> paramList = paramList_;
   paramList_ = Teuchos::null;
   return paramList;
}

/////////////////////////////////////////////////////

LinearOp BlockPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,BlockPreconditionerState & state) const
{
   // get the blocked linear operator
   RCP<LinearOpBase<double> > loA = Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(lo);
   BlockedLinearOp A = Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> >(loA);

   state.setInitialized(false);

   return buildPreconditionerOperator(A,state);
}

//! is this operator compatiable with the preconditioner factory?
bool BlockPreconditionerFactory::isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc) const
{
   RCP<const Thyra::PhysicallyBlockedLinearOpBase<double> > A 
         = Teuchos::rcp_dynamic_cast<const Thyra::PhysicallyBlockedLinearOpBase<double> >(fwdOpSrc.getOp());
   return A!=Teuchos::null;
}

//! create an instance of the preconditioner
RCP<Thyra::PreconditionerBase<double> > BlockPreconditionerFactory::createPrec() const
{
   // build a preconditioner, give it some inital state
   RCP<BlockPreconditioner> bp = rcp(new BlockPreconditioner());
   bp->setStateObject(buildPreconditionerState());
   bp->getStateObject()->setInitialized(false);

   return bp;
}

//! initialize a newly created preconditioner object
void BlockPreconditionerFactory::initializePrec(const RCP<const LinearOpSourceBase<double> > & ASrc,
                    PreconditionerBase<double> * prec,
                    const ESupportSolveUse supportSolveUse) const
{
   // get the blocked linear operator
   RCP<LinearOpBase<double> > loA = Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(ASrc->getOp());
   BlockedLinearOp A = Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> >(loA);

   BlockPreconditioner * blkPrec = dynamic_cast<BlockPreconditioner *>(prec);
   TEUCHOS_ASSERT(blkPrec!=0);
 
   // grab the state object
   RCP<BlockPreconditionerState> state = blkPrec->getStateObject();
   state->setInitialized(false);

   // build the preconditioner
   const RCP<const LinearOpBase<double> > M = buildPreconditionerOperator(A,*state);

   // must first cast that to be initialized
   DefaultPreconditioner<double> & dPrec = Teuchos::dyn_cast<DefaultPreconditioner<double> >(*prec);
   dPrec.initializeUnspecified(Teuchos::rcp_const_cast<LinearOpBase<double> >(M));
}

//! initialize a newly created preconditioner object
void BlockPreconditionerFactory::initializePrec(const RCP<const LinearOpSourceBase<double> > & ASrc,
                                                const RCP<const Thyra::MultiVectorBase<double> > & solnVec,
                                                PreconditionerBase<double> * prec,
                                                const ESupportSolveUse supportSolveUse) const
{
   BlockPreconditioner * blkPrec = dynamic_cast<BlockPreconditioner *>(prec);
   blkPrec->setSourceVector(Teuchos::rcp_const_cast<Thyra::MultiVectorBase<double> >(solnVec));

   initializePrec(ASrc,prec,supportSolveUse);
}

//! wipe clean a already initialized preconditioner object
void BlockPreconditionerFactory::uninitializePrec(PreconditionerBase<double> * prec, 
                      RCP<const LinearOpSourceBase<double> > * fwdOpSrc,
                      ESupportSolveUse *supportSolveUse) const
{
   // BlockPreconditioner * blkPrec = dynamic_cast<BlockPreconditioner *>(prec);

   // what do I do here?
   TEST_FOR_EXCEPT_MSG(true,"\"BlockPreconditionerFactory::uninitializePrec not implemented\"");
}

// for ParameterListAcceptor
///////////////////////////////////////////////////////////////////////

//! Set parameters from a parameter list and return with default values.
void BlockPreconditionerFactory::setParameterList(const RCP<ParameterList> & paramList)
{
   paramList_ = paramList;
}

//! Get the parameter list that was set using setParameterList().
RCP< ParameterList > BlockPreconditionerFactory::getNonconstParameterList()
{
   return paramList_;
}

//! Unset the parameter list that was set using setParameterList(). 
RCP< ParameterList > BlockPreconditionerFactory::unsetParameterList()
{
   RCP<ParameterList> _paramList = paramList_;
   paramList_ = Teuchos::null;
   return _paramList;
}

//! Set the inverse library used by this preconditioner factory
void BlockPreconditionerFactory::setInverseLibrary(const RCP<const InverseLibrary> & il)
{
   inverseLibrary_ = il;
}

//! Get the inverse library used by this preconditioner factory
RCP<const InverseLibrary> BlockPreconditionerFactory::getInverseLibrary() const
{
   // lazily build the inverse library only when needed
   if(inverseLibrary_==Teuchos::null)
      return InverseLibrary::buildFromStratimikos();

   return inverseLibrary_;
}

/////////////////////////////////////////////////////
// Static members and methods
/////////////////////////////////////////////////////

//! for creating the preconditioner factories objects
CloneFactory<BlockPreconditionerFactory> BlockPreconditionerFactory::precFactoryBuilder_;

/** \brief Builder function for creating preconditioner factories (yes
  *        this is a factory factory.
  *
  * Builder function for creating preconditioner factories (yes
  * this is a factory factory.
  * 
  * \param[in] name     String name of factory to build
  * \param[in] settings Parameter list describing the parameters for the
  *                     factory to build
  * \param[in] invLib   Inverse library for the factory to use.
  *
  * \returns If the name is associated with a preconditioner
  *          a pointer is returned, otherwise Teuchos::null is returned.
  */
RCP<BlockPreconditionerFactory> 
BlockPreconditionerFactory::buildPreconditionerFactory(const std::string & name,
                                                       const Teuchos::ParameterList & settings,
                                                       const RCP<const InverseLibrary> & invLib)
{
   Teko_DEBUG_MSG("Begin BlockPreconditionerFactory::buildPreconditionerFactory",10);

   // initialize the defaults if necessary
   if(precFactoryBuilder_.cloneCount()==0) initializePrecFactoryBuilder();

   // request the preconditioner factory from the CloneFactory
   RCP<BlockPreconditionerFactory> precFact = precFactoryBuilder_.build(name);

   Teko_DEBUG_MSG_BEGIN(5);
      DEBUG_STREAM << "Looked up \"" << name << "\"" << std::endl;
      DEBUG_STREAM << "Built " << precFact << std::endl;
   Teko_DEBUG_MSG_END();

   if(precFact==Teuchos::null)  
      return Teuchos::null;

   // add in the inverse library
   if(invLib!=Teuchos::null)
      precFact->setInverseLibrary(invLib);

   // now that inverse library has been set,
   // pass in the parameter list
   precFact->initializeFromParameterList(settings);

   Teko_DEBUG_MSG("End BlockPreconditionerFactory::buildPreconditionerFactory",10);

   return precFact;
}

/** \brief Add a preconditioner factory to the builder. This is done using the
  *        clone pattern. 
  *
  * Add a preconditioner factory to the builder. This is done using the
  * clone pattern. If your class does not support the Cloneable interface then
  * you can use the AutoClone class to construct your object.
  *
  * \note If this method is called twice with the same string, the latter clone pointer
  *       will be used.
  *
  * \param[in] name String to associate with this object
  * \param[in] clone Pointer to Cloneable object
  */
void BlockPreconditionerFactory::addPreconditionerFactory(const std::string & name,const RCP<Cloneable> & clone)
{
   // initialize the defaults if necessary
   if(precFactoryBuilder_.cloneCount()==0) initializePrecFactoryBuilder();

   // add clone to builder
   precFactoryBuilder_.addClone(name,clone); 
}

//! This is where the default objects are put into the precFactoryBuilder_
void BlockPreconditionerFactory::initializePrecFactoryBuilder()
{
   RCP<Cloneable> clone;

   // add various preconditioners to factory
   clone = rcp(new AutoClone<LU2x2PreconditionerFactory>());
   precFactoryBuilder_.addClone("Block LU2x2",clone);

   clone = rcp(new AutoClone<JacobiPreconditionerFactory>());
   precFactoryBuilder_.addClone("Block Jacobi",clone);

   clone = rcp(new AutoClone<GaussSeidelPreconditionerFactory>());
   precFactoryBuilder_.addClone("Block Gauss-Seidel",clone);

   clone = rcp(new AutoClone<AddPreconditionerFactory>());
   precFactoryBuilder_.addClone("Block Add",clone);

   clone = rcp(new AutoClone<MultPreconditionerFactory>());
   precFactoryBuilder_.addClone("Block Multiply",clone);

   clone = rcp(new AutoClone<NS::LSCPreconditionerFactory>());
   precFactoryBuilder_.addClone("NS LSC",clone);

   clone = rcp(new AutoClone<NS::SIMPLEPreconditionerFactory>());
   precFactoryBuilder_.addClone("NS SIMPLE",clone);
}

void BlockPreconditionerFactory::getPreconditionerFactoryNames(std::vector<std::string> & names)
{ 
   // initialize the defaults if necessary
   if(precFactoryBuilder_.cloneCount()==0) initializePrecFactoryBuilder();
   precFactoryBuilder_.getCloneNames(names); 
}

} // end namespace Teko
