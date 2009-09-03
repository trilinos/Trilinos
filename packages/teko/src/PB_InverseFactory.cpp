#include "PB_InverseFactory.hpp"

// Thyra includes
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// PB includes
#include "PB_Utilities.hpp"
#include "PB_BlockPreconditionerFactory.hpp"

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

/** \brief Constructor that takes a Thyra solve factory and 
  *        makes it look like an InverseFactory. This constructor
  *        also permits the passing of an "Extra Parameters" parameter
  *        list.
  *
  * Constructor that takes a Thyra solve factory and 
  * makes it look like an InverseFactory.  This constructor
  * also permits the passing of an "Extra Parameters" parameter
  * list to be used and updated through the "RequestedParameters" function.
  * 
  * \param[in] precFactory Thyra PreconditionerFactoryBase used for building 
  *                        the inverse.
  * \param[in] xtraParam Parameter list containing extra parameters.
  */
PreconditionerInverseFactory::PreconditionerInverseFactory(
              const Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > & precFactory,
              const Teuchos::RCP<const Teuchos::ParameterList> & xtraParam)
   : precFactory_(precFactory), extraParams_(rcp(new Teuchos::ParameterList(*xtraParam)))
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
   PB_DEBUG_MSG("BEGIN PreconditionerInverseFactory::rebuildInverse",10);

   RCP<Thyra::PreconditionerBase<double> > prec 
         = Teuchos::get_extra_data<RCP<Thyra::PreconditionerBase<double> > >(dest,"prec");

   precFactory_->initializePrec(Thyra::defaultLinearOpSource(source),&*prec);

   PB_DEBUG_MSG("END PreconditionerInverseFactory::rebuildInverse",10);
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

/** \brief Request the additional parameters this preconditioner factory
  *        needs. 
  *
  * Request the additonal parameters needed by this preconditioner factory.
  * The parameter list will have a set of fields that can be filled with 
  * the requested values. These fields include all requirements, even those
  * of the sub-solvers if there are any.  Once correctly filled the object
  * can be updated by calling the updateRequestedParameters with the filled
  * parameter list.
  *
  * \returns A parameter list with the requested parameters.
  *
  * \note The default implementation returns Teuchos::null.
  */
Teuchos::RCP<Teuchos::ParameterList> PreconditionerInverseFactory::getRequestedParameters() const
{
   Teuchos::RCP<BlockPreconditionerFactory> bpf = rcp_dynamic_cast<BlockPreconditionerFactory>(precFactory_);

   // request the parameters from a BPF is required
   if(bpf!=Teuchos::null) 
      return bpf->getRequestedParameters();

   // for non block preconditioners see if there are user requested additional parameters
   return extraParams_;
}

/** \brief Update this object with the fields from a parameter list.
  *
  * Update the requested fields using a parameter list. This method is
  * expected to pair with the getRequestedParameters method (i.e. the fields
  * requested are going to be update using this method).
  *
  * \param[in] pl Parameter list containing the requested parameters.
  *
  * \returns If the method succeeded (found all its required parameters) this
  *          method returns true, otherwise it returns false.
  *
  * \note The default implementation returns true (it does nothing!).
  */
bool PreconditionerInverseFactory::updateRequestedParameters(const Teuchos::ParameterList & pl)
{
   Teuchos::RCP<BlockPreconditionerFactory> bpf = rcp_dynamic_cast<BlockPreconditionerFactory>(precFactory_);

   // update the parameters of a BPF is required
   if(bpf!=Teuchos::null)
      return bpf->updateRequestedParameters(pl);

   // for non block preconditioners see if there are user requested additional parameters
   if(extraParams_==Teuchos::null)
      return true;

   Teuchos::ParameterList::ConstIterator itr;
   RCP<Teuchos::ParameterList> srcPl = precFactory_->unsetParameterList();

   // find name of settings sublist
   std::string subName = "";
   for(itr=srcPl->begin();itr!=srcPl->end();++itr) {
      // search for string with "Settings" in name
      if(itr->first.find("Settings")!=string::npos) {
         subName = itr->first;
         continue;
      }
   }

   // update fails if no settings list was found
   if(subName=="") {
      precFactory_->setParameterList(srcPl);
      return false;
   }

   // add extra parameters to list
   Teuchos::ParameterList & settingsList = srcPl->sublist(subName);
   for(itr=pl.begin();itr!=pl.end();++itr)
      settingsList.setEntry(itr->first,itr->second);

   // set the parameter list
   precFactory_->setParameterList(srcPl);

   return true;
}

//! Build an inverse operator using a factory and a linear operator
InverseLinearOp buildInverse(const InverseFactory & factory,const LinearOp & A)
{
   InverseLinearOp inv;
   try {
      inv = factory.buildInverse(A);
   }
   catch(std::exception & e) {
      RCP<Teuchos::FancyOStream> out = PB::getOutputStream();

      *out << "PB: \"buildInverse\" could not construct the inverse operator using ";
      *out << "\"" << factory.toString() << "\"" << std::endl;
      *out << std::endl;
      *out << "*** THROWN EXCEPTION ***\n";
      *out << e.what() << std::endl;
      *out << "************************\n";
      
      throw e;
   }

   return inv;
}

/** Using a prebuilt linear operator, use factory to build an inverse operator
  * given a new forward operator.
  */
void rebuildInverse(const InverseFactory & factory, const LinearOp & A, InverseLinearOp & invA)
{
   InverseLinearOp inv;
   try {
      factory.rebuildInverse(A,invA);
   } 
   catch(std::exception & e) {
      RCP<Teuchos::FancyOStream> out = PB::getOutputStream();

      *out << "PB: \"rebuildInverse\" could not construct the inverse operator using ";
      *out << "\"" << factory.toString() << "\"" << std::endl;
      *out << std::endl;
      *out << "*** THROWN EXCEPTION ***\n";
      *out << e.what() << std::endl;
      *out << "************************\n";
      
      throw e;
   }
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
RCP<InverseFactory> invFactoryFromParamList(const Teuchos::ParameterList & list,const std::string & type)
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
