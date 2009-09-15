#include "PB_InverseLibrary.hpp"

#include "PB_BlockPreconditionerFactory.hpp"

#include <algorithm>

using Teuchos::RCP;
using Teuchos::rcp;

namespace PB {

InverseLibrary::InverseLibrary()
{
   // setup some valid Stratimikos parameters
   /////////////////////////////////////////////

   // set valid solve factory names
   stratValidSolver_.push_back("Belos"); 
   stratValidSolver_.push_back("Amesos"); 
   stratValidSolver_.push_back("AztecOO"); 

   // set valid preconditioner factory name
   stratValidPrecond_.push_back("ML"); 
   stratValidPrecond_.push_back("Ifpack"); 

   // set valid PB preconditioner factory names
   blockValidPrecond_.push_back("Block Jacobi"); 
   blockValidPrecond_.push_back("Block Gauss-Seidel"); 
   blockValidPrecond_.push_back("Block Add"); 
   blockValidPrecond_.push_back("Block Multiply"); 
   blockValidPrecond_.push_back("NS LSC");
   blockValidPrecond_.push_back("NS SIMPLE");
}

//! add an unspecified inverse to the library
void InverseLibrary::addInverse(const std::string & label,const Teuchos::ParameterList & pl)
{
   // strip out the label
   const std::string type = pl.get<std::string>("Type");

   // copy the parameter list so we can modify it
   Teuchos::ParameterList settingsList;
   settingsList.set(type,pl);
   settingsList.sublist(type).remove("Type");

   // is this a Stratimikos preconditioner or solver
   if(std::find(stratValidPrecond_.begin(),stratValidPrecond_.end(),type)!=stratValidPrecond_.end()) {
      // this is a Stratimikos preconditioner factory
      addStratPrecond(label,type,settingsList);
   }
   else if(std::find(stratValidSolver_.begin(),stratValidSolver_.end(),type)!=stratValidSolver_.end()) {
      // this is a Stratimikos preconditioner factory
      addStratSolver(label,type,settingsList);
   }
   else if(std::find(blockValidPrecond_.begin(),blockValidPrecond_.end(),type)!=blockValidPrecond_.end()) {
      // this is a PB preconditioner factory
      addBlockPrecond(label,type,settingsList);
   }
   else {
      Teuchos::FancyOStream & os = *PB::getOutputStream();
      os << "ERROR: Could not find inverse type \"" << type 
         << "\" required by inverse name \"" << label << "\"" << std::endl;
      TEUCHOS_ASSERT(false);
   }
}

//! Add a Stratimikos solver with a label to the library
void InverseLibrary::addStratSolver(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl)
{
   // add some additional parameters onto the list
   RCP<Teuchos::ParameterList> stratList = rcp(new Teuchos::ParameterList());
   stratList->set("Linear Solver Type",type);
   stratList->set("Linear Solver Types",pl);
   stratList->set("Preconditioner Type","None");

   stratSolver_[label] = stratList;
}

//! Add a Stratimikos preconditioner with a label to the library
void InverseLibrary::addStratPrecond(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl)
{
   // add some additional parameters onto the list
   RCP<Teuchos::ParameterList> stratList = rcp(new Teuchos::ParameterList());
   stratList->set("Preconditioner Type",type);
   stratList->set("Preconditioner Types",pl);

   stratPrecond_[label] = stratList;
}

//! Add a PB preconditioner to the library with a label
void InverseLibrary::addBlockPrecond(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl)
{
   // add some additional parameters onto the list
   RCP<Teuchos::ParameterList> blockList = rcp(new Teuchos::ParameterList());
   blockList->set("Preconditioner Type",type);
   blockList->set("Preconditioner Settings",pl.sublist(type));

   // add the PB preconditioner parameter list into the library
   blockPrecond_[label] = blockList;
}

/** Get the fully constructed parameter list for a particular label
  *
  * \param[in] label Name used for the desired solver.
  *
  * \returns If the label is found in the library the corresponding parameter list
  *          is returned, otherwise <code>Teuchos::null</code> is returned.
  */
Teuchos::RCP<const Teuchos::ParameterList> InverseLibrary::getParameterList(const std::string & label) const
{
   std::map<std::string,RCP<const Teuchos::ParameterList> >::const_iterator itr;
   
   // check preconditioners
   itr = stratPrecond_.find(label);
   if(itr!=stratPrecond_.end()) return itr->second;
    
   // check solvers
   itr = stratSolver_.find(label);
   if(itr!=stratSolver_.end()) return itr->second;
   
   // check solvers
   itr = blockPrecond_.find(label);
   if(itr!=blockPrecond_.end()) return itr->second;

   return Teuchos::null;
}

//! Get the inverse factory associated with a particular label
Teuchos::RCP<InverseFactory> InverseLibrary::getInverseFactory(const std::string & label) const
{
   PB_DEBUG_MSG("Begin InverseLibrary::getInverseFactory",10);

   std::map<std::string,RCP<const Teuchos::ParameterList> >::const_iterator itr;

   bool isStratSolver=false,isStratPrecond=false,isBlockPrecond=false;

   // is this a Stratimikos solver?
   itr = stratPrecond_.find(label);
   isStratPrecond = itr!=stratPrecond_.end();

   // is this a Stratimikos preconditioner?
   if(not isStratPrecond) {
      itr = stratSolver_.find(label);
      isStratSolver = itr!=stratSolver_.end();
   }

   // must be a "block" preconditioner
   if(not (isStratSolver || isStratPrecond)) {
      itr = blockPrecond_.find(label);
      isBlockPrecond = itr!=blockPrecond_.end();
   }

   PB_DEBUG_MSG("PB: Inverse \"" << label << "\" is of type " 
             << "strat prec = " << isStratPrecond << ", "
             << "strat solv = " << isStratSolver << ", " 
             << "block prec = " << isBlockPrecond,3);

   // Must be one of Strat solver, strat preconditioner, block preconditioner
   if(not (isStratSolver || isStratPrecond || isBlockPrecond)) {
      RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

      *out << "PB: getInverseFactory could not find \"" << label << "\" ... aborting\n";
      *out << std::endl;

      TEUCHOS_ASSERT(isStratSolver || isStratPrecond || isBlockPrecond);
   }
   
   RCP<const Teuchos::ParameterList> pl = itr->second;

   // build inverse factory
   if(isStratPrecond) {
      // remove required parameters
      RCP<Teuchos::ParameterList> plCopy = rcp(new Teuchos::ParameterList(*pl));
      std::string type = plCopy->get<std::string>("Preconditioner Type");
      RCP<Teuchos::ParameterList> xtraParams = rcp(new Teuchos::ParameterList(
            plCopy->sublist("Preconditioner Types").sublist(type).sublist("Required Parameters"))); 
      plCopy->sublist("Preconditioner Types").sublist(type).remove("Required Parameters"); 

      // print some debuggin info
      PB_DEBUG_MSG_BEGIN(10);
         DEBUG_STREAM << "Printing parameter list: " << std::endl; 
         PB_DEBUG_PUSHTAB(); plCopy->print(DEBUG_STREAM); PB_DEBUG_POPTAB();

         DEBUG_STREAM << "Printing extra parameters: " << std::endl; 
         PB_DEBUG_PUSHTAB(); xtraParams->print(DEBUG_STREAM); PB_DEBUG_POPTAB();
      PB_DEBUG_MSG_END();

      Stratimikos::DefaultLinearSolverBuilder strat;
      strat.setParameterList(plCopy);

      // try to build a preconditioner factory
      RCP<Thyra::PreconditionerFactoryBase<double> > precFact = strat.createPreconditioningStrategy(type);

      PB_DEBUG_MSG("End InverseLibrary::getInverseFactory (Stratimikos preconditioner)",10);

      // string must map to a preconditioner
      return rcp(new PreconditionerInverseFactory(precFact,xtraParams));
   }
   else if(isStratSolver) {
      Stratimikos::DefaultLinearSolverBuilder strat;
      strat.setParameterList(rcp(new Teuchos::ParameterList(*pl)));

      // try to build a solver factory
      std::string type = pl->get<std::string>("Linear Solver Type");
      RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solveFact = strat.createLinearSolveStrategy(type);

      PB_DEBUG_MSG("End InverseLibrary::getInverseFactory (Stratimikos solver)",10);

      // if its around, build a InverseFactory
      return rcp(new SolveInverseFactory(solveFact));
   }
   else if(isBlockPrecond) {
      try {
         std::string type = pl->get<std::string>("Preconditioner Type");
         const Teuchos::ParameterList & settings = pl->sublist("Preconditioner Settings");
   
         // build preconditioner factory from the string
         RCP<BlockPreconditionerFactory> precFact 
               = BlockPreconditionerFactory::buildPreconditionerFactory(type,settings,Teuchos::rcpFromRef(*this));
    
         TEUCHOS_ASSERT(precFact!=Teuchos::null);
   
         PB_DEBUG_MSG("End InverseLibrary::getInverseFactory (Block preconditioner)",10);

         // return the inverse factory object
         return rcp(new PreconditionerInverseFactory(precFact));   
      }
      catch(std::exception & e) {
         RCP<Teuchos::FancyOStream> out = PB::getOutputStream();
         
         *out << "PB: \"getInverseFactory\" failed, Parameter List =\n";
         pl->print(*out);

         *out << "*** THROWN EXCEPTION ***\n";
         *out << e.what() << std::endl;
         *out << "************************\n";
         
         throw e;
      }
   }

   PB_DEBUG_MSG("End InverseLibrary::getInverseFactory (FAILURE)",10);

   TEUCHOS_ASSERT(false);
}

/** \brief Build an inverse library from a parameter list.
  * 
  * Build an inverse library from a parameter list. This will
  * contain all the labeled inverses specified.
  *
  * \param[in] pl Parameter list to build the library from
  *
  * \returns A pointer to the inverse library created.
  */
RCP<InverseLibrary> InverseLibrary::buildFromParameterList(const Teuchos::ParameterList & pl,bool useStratDefaults)
{
   // build from Stratimikos or allocate a new inverse library
   RCP<InverseLibrary> invLib;
   if(useStratDefaults)
      invLib = InverseLibrary::buildFromStratimikos();
   else
      invLib = rcp(new InverseLibrary());

   // to convert the void* like entry
   Teuchos::ParameterList * temp;

   // loop over all entries in parameter list
   Teuchos::ParameterList::ConstIterator itr;
   for(itr=pl.begin();itr!=pl.end();++itr) {
      // get current entry
      std::string label             = itr->first;
      Teuchos::ParameterList & list = itr->second.getValue(temp);
      
      // add to library
      invLib->addInverse(label,list);
   }
   
   return invLib;
}

/** \brief Build an inverse library from Stratimikos
  * 
  * Build an inverse library from Stratimkos. The labels
  * will just be the names in Stratimikos.
  *
  * \param[in] strat Stratimikos object to use
  *
  * \returns A pointer to the inverse library created.
  */
Teuchos::RCP<InverseLibrary> InverseLibrary::buildFromStratimikos(const Stratimikos::DefaultLinearSolverBuilder & strat)
{
   RCP<InverseLibrary> invLib = rcp(new InverseLibrary());

   // get default inveres in Stratimikos
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList(*strat.getValidParameters()));
   Teuchos::ParameterList lst(pl->sublist("Linear Solver Types"));
   Teuchos::ParameterList pft(pl->sublist("Preconditioner Types"));

   Teuchos::ParameterList::ConstIterator itr;
   Teuchos::ParameterList * temp;

   // loop over all entries in solver list
   for(itr=lst.begin();itr!=lst.end();++itr) {
      // get current entry
      std::string label             = itr->first;
      Teuchos::ParameterList & list = itr->second.getValue(temp);
      list.set("Type",label);
      
      // add to library
      invLib->addInverse(label,list);
   }

   // loop over all entries in preconditioner list
   for(itr=pft.begin();itr!=pft.end();++itr) {
      // get current entry
      std::string label             = itr->first;
      Teuchos::ParameterList & list = itr->second.getValue(temp);
      list.set("Type",label);
      
      // add to library
      invLib->addInverse(label,list);
   }

   return invLib;
}

} // end namespace PB
