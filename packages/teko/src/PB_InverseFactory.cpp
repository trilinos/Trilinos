#include "PB_InverseFactory.hpp"

// Thyra includes
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// PB includes
#include "PB_Utilities.hpp"

using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;

namespace PB {

/** \brief Constructor that takes a Thyra solve factory and 
  *        makes it look like an InverseFactory
  *
  * Constructor that takes a Thyra solve factory and 
  * makes it look like an InverseFactory.
  * 
  * \param[in] lowsFactory Thyra LineaerOpWithSolveFactoryBase used for building 
  *                        the inverse.
  */
SolveInverseFactory::SolveInverseFactory(const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > & lowsFactory)
   : lowsFactory_(lowsFactory)
{ }

//! Copy constructor
SolveInverseFactory::SolveInverseFactory(const SolveInverseFactory & siFactory)
   : lowsFactory_(siFactory.lowsFactory_)
{ }

/** \brief Build an inverse operator
  *
  * Build the inverse operator using this factory.
  *
  * \param[in] linearOp Linear operator needing to be inverted.
  *
  * \returns New linear operator that functions as the inverse
  *          of <code>linearOp</code>.
  */
InverseLinearOp SolveInverseFactory::buildInverse(const LinearOp & linearOp) const
{
   // build and initialize inverse linear op with solve
   Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > invLOWS = lowsFactory_->createOp();
   lowsFactory_->initializeOp(Thyra::defaultLinearOpSource(linearOp),&*invLOWS,Thyra::SUPPORT_SOLVE_FORWARD_ONLY);
   
   return Thyra::nonconstInverse<double>(invLOWS);
}

/** \brief Pass in an already constructed inverse operator. Update
  *        the inverse operator based on the new source operator.
  *
  * Pass in an already constructed inverse operator. Update
  * the inverse operator based on the new source operator.
  *
  * \params[in]     source Source operator to be inverted.
  * \params[in,out] dest   Pre constructed inverse operator to be
  *                        rebuilt using the <code>source</code>
  *                        object.
  */
void SolveInverseFactory::rebuildInverse(const LinearOp & source,InverseLinearOp & dest) const
{
   // RCP<Thyra::LinearOpBase<double> > nonConstDest = rcp_const_cast<Thyra::LinearOpBase<double> >(dest);
   // RCP<Thyra::DefaultInverseLinearOp<double> > invDest = rcp_dynamic_cast<Thyra::DefaultInverseLinearOp<double> >(nonConstDest);
   RCP<Thyra::DefaultInverseLinearOp<double> > invDest = rcp_dynamic_cast<Thyra::DefaultInverseLinearOp<double> >(dest);
   RCP<Thyra::LinearOpWithSolveBase<double> > lows = invDest->getNonconstLows();

   lowsFactory_->initializeAndReuseOp(Thyra::defaultLinearOpSource(source),&*lows);
}

/** \brief A function that permits inspection of the parameters used to create
  *        this object.
  *
  * A function that permits inspection of the parameters used to create this
  * object. Useful for determining defaults and settings used.
  *
  * \returns A list used to parameterize this object.
  */
Teuchos::RCP<const Teuchos::ParameterList> SolveInverseFactory::getParameterList() const
{ 
   return lowsFactory_->getParameterList(); 
}

/** \brief Constructor that takes a Thyra solve factory and 
  *        makes it look like an InverseFactory
  *
  * Constructor that takes a Thyra solve factory and 
  * makes it look like an InverseFactory.
  * 
  * \param[in] precFactory Thyra PreconditionerFactoryBase used for building 
  *                        the inverse.
  */
PreconditionerInverseFactory::PreconditionerInverseFactory(const Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > & precFactory)
   : precFactory_(precFactory)
{ }

//! Copy constructor
PreconditionerInverseFactory::PreconditionerInverseFactory(const PreconditionerInverseFactory & pFactory)
   : precFactory_(pFactory.precFactory_)
{ }

/** \brief Build an inverse operator
  *
  * Build the inverse operator using this factory. 
  *
  * \param[in] linearOp Linear operator needing to be inverted.
  *
  * \returns New linear operator that functions as the inverse
  *          of <code>linearOp</code>.
  */
InverseLinearOp PreconditionerInverseFactory::buildInverse(const LinearOp & linearOp) const
{
   RCP<Thyra::PreconditionerBase<double> > prec = precFactory_->createPrec();
   precFactory_->initializePrec(Thyra::defaultLinearOpSource(linearOp),&*prec);

   RCP<Thyra::LinearOpBase<double> > precOp = prec->getNonconstUnspecifiedPrecOp();
   Teuchos::set_extra_data(prec,"prec",Teuchos::inOutArg(precOp));

   return precOp;
}

/** \brief Pass in an already constructed inverse operator. Update
  *        the inverse operator based on the new source operator.
  *
  * Pass in an already constructed inverse operator. Update
  * the inverse operator based on the new source operator.
  *
  * \params[in]     source Source operator to be inverted.
  * \params[in,out] dest   Pre constructed inverse operator to be
  *                        rebuilt using the <code>source</code>
  *                        object.
  */
void PreconditionerInverseFactory::rebuildInverse(const LinearOp & source,InverseLinearOp & dest) const
{
   RCP<Thyra::PreconditionerBase<double> > prec 
         = Teuchos::get_extra_data<RCP<Thyra::PreconditionerBase<double> > >(dest,"prec");

   precFactory_->initializePrec(Thyra::defaultLinearOpSource(source),&*prec);
}

/** \brief A function that permits inspection of the parameters used to create
  *        this object.
  *
  * A function that permits inspection of the parameters used to create this
  * object. Useful for determining defaults and settings used.
  *
  * \returns A list used to parameterize this object.
  */
Teuchos::RCP<const Teuchos::ParameterList> PreconditionerInverseFactory::getParameterList() const
{ 
   return precFactory_->getParameterList(); 
}

//! Build an inverse operator using a factory and a linear operator
InverseLinearOp buildInverse(const InverseFactory & factory,const LinearOp & A)
{
   return factory.buildInverse(A);
}

/** Using a prebuilt linear operator, use factory to build an inverse operator
  * given a new forward operator.
  */
void rebuildInverse(const InverseFactory & factory, const LinearOp & A, InverseLinearOp & invA)
{
   return factory.rebuildInverse(A,invA);
}

/** \brief Build an InverseFactory object from a ParameterList, as specified in Stratimikos.
  *
  * Build an InverseFactory object from a ParameterList, as specified in Stratimikos.
  * The specific inverse routine (either solver or preconditioner) to be chosen is specified
  * by a string.
  *
  * \param[in] list ParameterList that describes the available solvers/preconditioners.
  * \param[in] type String saying which solver/preconditioner to use.
  *
  * \returns An inverse factory using the specified inverse operation.
  */
RCP<const InverseFactory> invFactoryFromParamList(const Teuchos::ParameterList & list,const std::string & type)
{
   RCP<Teuchos::ParameterList> myList = rcp(new Teuchos::ParameterList(list));

   Stratimikos::DefaultLinearSolverBuilder strat;
   strat.setParameterList(myList);

   try {
      // try to build a preconditioner factory
      RCP<Thyra::PreconditionerFactoryBase<double> > precFact = strat.createPreconditioningStrategy(type);

      // string must map to a preconditioner
      return rcp(new PreconditionerInverseFactory(precFact));
   }
   catch(const Teuchos::Exceptions::InvalidParameterValue & exp) { }
 
   try {
      // try to build a solver factory
      RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solveFact = strat.createLinearSolveStrategy(type);

      // if its around, build a InverseFactory
      return rcp(new SolveInverseFactory(solveFact));
   }
   catch(const Teuchos::Exceptions::InvalidParameterValue & exp) { }

   return  Teuchos::null;;
}

/** \brief Get a valid parameter list for the inverse factory class.
  *
  * Get a valid parameter list for the inverse factory class. This will
  * specify the set of parameters for each possible "inverse".
  *
  * \returns A parameter list is returned that is suitable to be passed
  *          to <code>invFactoryFromParamList</code>.
  */
Teuchos::RCP<const Teuchos::ParameterList> invFactoryValidParameters()
{
   Stratimikos::DefaultLinearSolverBuilder strat;
 
   // extract valid parameter list from Stratimikos
   return strat.getValidParameters();
}

} // end namespace PB
