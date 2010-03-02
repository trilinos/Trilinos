// Teko includes
#include "Teko_RepeatPreconditionerFactory.hpp"

namespace Teko {

//! Default constructor, for use with the AutoClone class.
RepeatPreconditionerFactory::RepeatPreconditionerFactory()
   : correctionNum_(0), precFactory_(Teuchos::null)
{ }

/** Construct a preconditioner factory that applies a specified
  * preconditioner, a fixed number of times.
  */
RepeatPreconditionerFactory::RepeatPreconditionerFactory(unsigned int correctionNum,
                            const Teuchos::RCP<Teko::InverseFactory> & precFactory)
   : correctionNum_(correctionNum), precFactory_(precFactory)
{ }

/** Construct a preconditioner factory that applies a specified
  * preconditioner, a fixed number of times.
  */
RepeatPreconditionerFactory::RepeatPreconditionerFactory(unsigned int correctionNum,
                            const Teuchos::RCP<Teko::PreconditionerFactory> & precFactory)
   : correctionNum_(correctionNum)
{
   precFactory_ = Teuchos::rcp(new Teko::PreconditionerInverseFactory(precFactory));
}


/** \brief Function that is called to build the preconditioner
  *        for the linear operator that is passed in.
  */
LinearOp RepeatPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const
{
   TEST_FOR_EXCEPTION(precFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::RepeatPreconditionerFactory::buildPreconditionerOperator requires that a "
                   << "preconditioner factory has been set. Currently it is null!");

   // build user specified preconditioner
   InverseLinearOp & invP = state.getModifiableOp("prec");
   if(invP==Teuchos::null)
      invP = Teko::buildInverse(precFactory_,lo);
   else
      Teko::rebuildInverse(precFactory_,lo,invP);

   // no repititions are required
   if(correctionNum_==0) return invP;

   // now build correction operator
   LinearOp I = Thyra::identity(lo->range());
   LinearOp AiP = multiply(lo,invP);
   LinearOp correction = add(I,scale(-1.0,AiP)); // I - A * iPA

   LinearOp resMap = I; // will map r_n to r_{n+1}
   for(int i=0;i<correctionNum_;i++) 
      resMap = multiply(correction,resMap); // resMap = (I-A*iP)*resMap
   
   // iP = (I-A*iP)^{correctionNum}
   return multiply(invP,resMap);
}

/** \brief This function builds the internals of the preconditioner factory
  *        from a parameter list.
  */
void RepeatPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & settings)
{
   correctionNum_ = settings.get<int>("Iteration Count",1);
 
   TEST_FOR_EXCEPTION(not settings.isParameter("Preconditioner Type"),std::runtime_error,
                      "Parameter \"Preconditioner Type\" is required by a Teko::RepeatPreconditionerFactory");
      
   // grab library and preconditioner name
   RCP<const InverseLibrary> il = getInverseLibrary();
   std::string precName = settings.get<std::string>("Preconditioner Type");

   // build preconditioner factory
   precFactory_ = il->getInverseFactory(precName);
}

/** \brief Request the additional parameters this preconditioner factory
  *        needs. 
  */
Teuchos::RCP<Teuchos::ParameterList> RepeatPreconditionerFactory::getRequestedParameters() const
{
   TEST_FOR_EXCEPTION(precFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::RepeatPreconditionerFactory::getRequestedParameters requires that a "
                   << "preconditioner factory has been set. Currently it is null!");

   return precFactory_->getRequestedParameters();
}

/** \brief Update this object with the fields from a parameter list.
  */
bool RepeatPreconditionerFactory::updateRequestedParameters(const Teuchos::ParameterList & pl)
{
   TEST_FOR_EXCEPTION(precFactory_==Teuchos::null,std::runtime_error,
                      "ERROR: Teko::RepeatPreconditionerFactory::updateRequestedParameters requires that a "
                   << "preconditioner factory has been set. Currently it is null!");

   return precFactory_->updatedRequestedParameters(pl);
}

} // end namespace Teko
