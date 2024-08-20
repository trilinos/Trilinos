// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_InverseLibrary.hpp"

#include "Teko_SolveInverseFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_BlockPreconditionerFactory.hpp"

#include "Teko_NeumannSeriesPreconditionerFactory.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teko_Utilities.hpp"

#include <algorithm>

using Teuchos::RCP;
using Teuchos::rcp;

namespace Teko {

/** This function adds some additional preconditioners to Stratimikos.
 * These are NOT block preconditioners.
 */
void addToStratimikosBuilder(const RCP<Stratimikos::DefaultLinearSolverBuilder>& builder) {
  typedef Thyra::PreconditionerFactoryBase<double> PrecFactory;

  RCP<const Teuchos::ParameterList> parameters = builder->getValidParameters();

  if (!parameters->sublist("Preconditioner Types").isSublist("Neumann Series")) {
    RCP<const Teuchos::AbstractFactory<Thyra::PreconditionerFactoryBase<double> > > factory;

    factory = Teuchos::abstractFactoryStd<PrecFactory,
                                          Teko::NeumannSeriesPreconditionerFactory<double> >();
    builder->setPreconditioningStrategyFactory(factory, "Neumann Series");
  }
}

InverseLibrary::InverseLibrary() {
  Teko_DEBUG_SCOPE("InverseLibrary::InverseLibrary", 10);

  defaultBuilder_ = Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder());
  addToStratimikosBuilder(defaultBuilder_);

  // setup some valid Stratimikos parameters
  /////////////////////////////////////////////

  // set valid solve factory names
  stratValidSolver_.push_back("Belos");
  stratValidSolver_.push_back("Amesos");
  stratValidSolver_.push_back("Amesos2");
  stratValidSolver_.push_back("AztecOO");

  // set valid preconditioner factory name
  stratValidPrecond_.push_back("ML");
  stratValidPrecond_.push_back("Ifpack");
  stratValidPrecond_.push_back("Neumann Series");
  stratValidPrecond_.push_back("MueLu");
  stratValidPrecond_.push_back("Ifpack2");

  // set valid Teko preconditioner factory names
  PreconditionerFactory::getPreconditionerFactoryNames(blockValidPrecond_);

  Teko_DEBUG_MSG_BEGIN(10) DEBUG_STREAM << "Loaded \"block\" preconditioners = ";
  for (std::size_t i = 0; i < blockValidPrecond_.size(); i++)
    DEBUG_STREAM << blockValidPrecond_[i] << ", ";
  DEBUG_STREAM << std::endl;
  Teko_DEBUG_MSG_END()
}

InverseLibrary::InverseLibrary(const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder>& strat)
    : defaultBuilder_(strat) {
  Teko_DEBUG_SCOPE("InverseLibrary::InverseLibrary", 10);

  RCP<Teuchos::ParameterList> pl =
      rcp(new Teuchos::ParameterList(*defaultBuilder_->getValidParameters()));
  Teuchos::ParameterList lst(pl->sublist("Linear Solver Types"));
  Teuchos::ParameterList pft(pl->sublist("Preconditioner Types"));

  Teuchos::ParameterList::ConstIterator itr;

  // set valid solve factory names
  for (itr = lst.begin(); itr != lst.end(); ++itr) stratValidSolver_.push_back(itr->first);

  Teko_DEBUG_MSG_BEGIN(10) DEBUG_STREAM << "Loaded \"Stratimikos\" solvers = ";
  for (std::size_t i = 0; i < stratValidSolver_.size(); i++)
    DEBUG_STREAM << stratValidSolver_[i] << ", ";
  DEBUG_STREAM << std::endl;
  Teko_DEBUG_MSG_END()

      // set valid prec factory names
      for (itr = pft.begin(); itr != pft.end(); ++itr) stratValidPrecond_.push_back(itr->first);

  Teko_DEBUG_MSG_BEGIN(10) DEBUG_STREAM << "Loaded \"Stratimikos\" preconditioners = ";
  for (std::size_t i = 0; i < stratValidPrecond_.size(); i++)
    DEBUG_STREAM << stratValidPrecond_[i] << ", ";
  DEBUG_STREAM << std::endl;
  Teko_DEBUG_MSG_END()

      // set valid Teko preconditioner factory names
      PreconditionerFactory::getPreconditionerFactoryNames(blockValidPrecond_);

  Teko_DEBUG_MSG_BEGIN(10) DEBUG_STREAM << "Loaded \"block\" preconditioners = ";
  for (std::size_t i = 0; i < blockValidPrecond_.size(); i++)
    DEBUG_STREAM << blockValidPrecond_[i] << ", ";
  DEBUG_STREAM << std::endl;
  Teko_DEBUG_MSG_END()
}

//! add an unspecified inverse to the library
void InverseLibrary::addInverse(const std::string& label, const Teuchos::ParameterList& pl) {
  // strip out the label
  const std::string type = pl.get<std::string>("Type");

  // copy the parameter list so we can modify it
  Teuchos::ParameterList settingsList;
  settingsList.set(type, pl);
  settingsList.sublist(type).remove("Type");

  // is this a Stratimikos preconditioner or solver
  if (std::find(stratValidPrecond_.begin(), stratValidPrecond_.end(), type) !=
      stratValidPrecond_.end()) {
    // this is a Stratimikos preconditioner factory
    addStratPrecond(label, type, settingsList);
  } else if (std::find(stratValidSolver_.begin(), stratValidSolver_.end(), type) !=
             stratValidSolver_.end()) {
    // this is a Stratimikos preconditioner factory
    addStratSolver(label, type, settingsList);
  } else if (std::find(blockValidPrecond_.begin(), blockValidPrecond_.end(), type) !=
             blockValidPrecond_.end()) {
    // this is a Teko preconditioner factory
    addBlockPrecond(label, type, settingsList);
  } else {
    Teuchos::FancyOStream& os = *Teko::getOutputStream();
    os << "ERROR: Could not find inverse type \"" << type << "\" required by inverse name \""
       << label << "\"" << std::endl;
    TEUCHOS_ASSERT(false);
  }
}

//! Add a Stratimikos solver with a label to the library
void InverseLibrary::addStratSolver(const std::string& label, const std::string& type,
                                    const Teuchos::ParameterList& pl) {
  // add some additional parameters onto the list
  RCP<Teuchos::ParameterList> stratList = rcp(new Teuchos::ParameterList());
  stratList->set("Linear Solver Type", type);
  stratList->set("Linear Solver Types", pl);
  stratList->set("Preconditioner Type", "None");

  stratSolver_[label] = stratList;
}

//! Add a Stratimikos preconditioner with a label to the library
void InverseLibrary::addStratPrecond(const std::string& label, const std::string& type,
                                     const Teuchos::ParameterList& pl) {
  // add some additional parameters onto the list
  RCP<Teuchos::ParameterList> stratList = rcp(new Teuchos::ParameterList());
  stratList->set("Preconditioner Type", type);
  stratList->set("Preconditioner Types", pl);

  stratPrecond_[label] = stratList;
}

//! Add a Teko preconditioner to the library with a label
void InverseLibrary::addBlockPrecond(const std::string& label, const std::string& type,
                                     const Teuchos::ParameterList& pl) {
  // add some additional parameters onto the list
  RCP<Teuchos::ParameterList> blockList = rcp(new Teuchos::ParameterList());
  blockList->set("Preconditioner Type", type);
  blockList->set("Preconditioner Settings", pl.sublist(type));

  // add the Teko preconditioner parameter list into the library
  blockPrecond_[label] = blockList;
}

/** Get the fully constructed parameter list for a particular label
 *
 * \param[in] label Name used for the desired solver.
 *
 * \returns If the label is found in the library the corresponding parameter list
 *          is returned, otherwise <code>Teuchos::null</code> is returned.
 */
Teuchos::RCP<const Teuchos::ParameterList> InverseLibrary::getParameterList(
    const std::string& label) const {
  std::map<std::string, RCP<const Teuchos::ParameterList> >::const_iterator itr;

  // check preconditioners
  itr = stratPrecond_.find(label);
  if (itr != stratPrecond_.end()) return itr->second;

  // check solvers
  itr = stratSolver_.find(label);
  if (itr != stratSolver_.end()) return itr->second;

  // check solvers
  itr = blockPrecond_.find(label);
  if (itr != blockPrecond_.end()) return itr->second;

  return Teuchos::null;
}

//! Get the inverse factory associated with a particular label
Teuchos::RCP<InverseFactory> InverseLibrary::getInverseFactory(const std::string& label) const {
  Teko_DEBUG_SCOPE("InverseLibrary::getInverseFactory", 10);

  std::map<std::string, RCP<const Teuchos::ParameterList> >::const_iterator itr;

  bool isStratSolver = false, isStratPrecond = false, isBlockPrecond = false;

  // is this a Stratimikos solver?
  itr            = stratPrecond_.find(label);
  isStratPrecond = itr != stratPrecond_.end();

  // is this a Stratimikos preconditioner?
  if (not isStratPrecond) {
    itr           = stratSolver_.find(label);
    isStratSolver = itr != stratSolver_.end();
  }

  // must be a "block" preconditioner
  if (not(isStratSolver || isStratPrecond)) {
    itr            = blockPrecond_.find(label);
    isBlockPrecond = itr != blockPrecond_.end();
  }

  Teko_DEBUG_MSG("Inverse \"" << label << "\" is of type "
                              << "strat prec = " << isStratPrecond << ", "
                              << "strat solv = " << isStratSolver << ", "
                              << "block prec = " << isBlockPrecond,
                 3);

  // Must be one of Strat solver, strat preconditioner, block preconditioner
  if (not(isStratSolver || isStratPrecond || isBlockPrecond)) {
    RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    *out << "InverseLibrary::getInverseFactory could not find \"" << label << "\" ... aborting\n";
    *out << "Choose one of: " << std::endl;

    *out << "   Stratimikos preconditioners = ";
    for (itr = stratPrecond_.begin(); itr != stratPrecond_.end(); ++itr)
      *out << "      \"" << itr->first << "\"\n";
    *out << std::endl;

    *out << "   Stratimikos solvers = ";
    for (itr = stratSolver_.begin(); itr != stratSolver_.end(); ++itr)
      *out << "      \"" << itr->first << "\"\n";
    *out << std::endl;

    *out << "   Block preconditioners = ";
    for (itr = blockPrecond_.begin(); itr != blockPrecond_.end(); ++itr)
      *out << "      \"" << itr->first << "\"\n";
    *out << std::endl;

    TEUCHOS_ASSERT(isStratSolver || isStratPrecond || isBlockPrecond);
  }

  RCP<const Teuchos::ParameterList> pl = itr->second;

  // build inverse factory
  if (isStratPrecond) {
    // remove required parameters
    RCP<Teuchos::ParameterList> plCopy = rcp(new Teuchos::ParameterList(*pl));
    std::string type                   = plCopy->get<std::string>("Preconditioner Type");
    RCP<Teuchos::ParameterList> xtraParams;
    if (plCopy->sublist("Preconditioner Types").sublist(type).isParameter("Required Parameters")) {
      xtraParams = rcp(new Teuchos::ParameterList(
          plCopy->sublist("Preconditioner Types").sublist(type).sublist("Required Parameters")));
      plCopy->sublist("Preconditioner Types").sublist(type).remove("Required Parameters");
    }

    // print some debuggin info
    Teko_DEBUG_MSG_BEGIN(10);
    DEBUG_STREAM << "Printing parameter list: " << std::endl;
    Teko_DEBUG_PUSHTAB();
    plCopy->print(DEBUG_STREAM);
    Teko_DEBUG_POPTAB();

    if (xtraParams != Teuchos::null) {
      DEBUG_STREAM << "Printing extra parameters: " << std::endl;
      Teko_DEBUG_PUSHTAB();
      xtraParams->print(DEBUG_STREAM);
      Teko_DEBUG_POPTAB();
    }
    Teko_DEBUG_MSG_END();

    // Stratimikos::DefaultLinearSolverBuilder strat;
    // addToStratimikosBuilder(strat);
    defaultBuilder_->setParameterList(plCopy);

    // try to build a preconditioner factory
    RCP<Thyra::PreconditionerFactoryBase<double> > precFact =
        defaultBuilder_->createPreconditioningStrategy(type);

    // string must map to a preconditioner
    RCP<Teko::PreconditionerInverseFactory> precInvFact =
        rcp(new PreconditionerInverseFactory(precFact, xtraParams, getRequestHandler()));
    precInvFact->setupParameterListFromRequestHandler();
    return precInvFact;
  } else if (isStratSolver) {
    RCP<Teuchos::ParameterList> solveList = rcp(new Teuchos::ParameterList(*pl));
    std::string type                      = solveList->get<std::string>("Linear Solver Type");

    // get preconditioner name, remove "Use Preconditioner" parameter
    Teuchos::ParameterList& solveSettings = solveList->sublist("Linear Solver Types").sublist(type);
    std::string precKeyWord               = "Use Preconditioner";
    std::string precName                  = "None";
    if (solveSettings.isParameter(precKeyWord)) {
      precName = solveSettings.get<std::string>(precKeyWord);
      solveSettings.remove(precKeyWord);
    }

    // build Thyra preconditioner factory
    RCP<Thyra::PreconditionerFactoryBase<double> > precFactory;
    if (precName != "None") {
      // we will manually set the preconditioner, so set this to null
      solveList->set<std::string>("Preconditioner Type", "None");

      // build inverse that preconditioner corresponds to
      RCP<PreconditionerInverseFactory> precInvFactory =
          Teuchos::rcp_dynamic_cast<PreconditionerInverseFactory>(getInverseFactory(precName));

      // extract preconditioner factory from preconditioner _inverse_ factory
      precFactory = precInvFactory->getPrecFactory();
    }

    // Stratimikos::DefaultLinearSolverBuilder strat;
    // addToStratimikosBuilder(strat);
    defaultBuilder_->setParameterList(solveList);

    // try to build a solver factory
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solveFact =
        defaultBuilder_->createLinearSolveStrategy(type);
    if (precFactory != Teuchos::null) solveFact->setPreconditionerFactory(precFactory, precName);

    // if its around, build a InverseFactory
    return rcp(new SolveInverseFactory(solveFact));
  } else if (isBlockPrecond) {
    try {
      std::string type                       = pl->get<std::string>("Preconditioner Type");
      const Teuchos::ParameterList& settings = pl->sublist("Preconditioner Settings");

      // build preconditioner factory from the string
      RCP<PreconditionerFactory> precFact = PreconditionerFactory::buildPreconditionerFactory(
          type, settings, Teuchos::rcpFromRef(*this));

      TEUCHOS_ASSERT(precFact != Teuchos::null);

      // return the inverse factory object
      return rcp(new PreconditionerInverseFactory(precFact, getRequestHandler()));
    } catch (std::exception& e) {
      RCP<Teuchos::FancyOStream> out = Teko::getOutputStream();

      *out << "Teko: \"getInverseFactory\" failed, Parameter List =\n";
      pl->print(*out);

      *out << "*** THROWN EXCEPTION ***\n";
      *out << e.what() << std::endl;
      *out << "************************\n";

      throw e;
    }
  }

  TEUCHOS_ASSERT(false);
}

//! Print the inverses and parameter lists available for use
void InverseLibrary::PrintAvailableInverses(std::ostream& os) const {
  std::map<std::string, Teuchos::RCP<const Teuchos::ParameterList> >::const_iterator itr;

  os << "Stratimikos Solvers: " << std::endl;
  os << "********************************" << std::endl;
  for (itr = stratSolver_.begin(); itr != stratSolver_.end(); ++itr) {
    os << "name = \"" << itr->first << "\"" << std::endl;
    itr->second->print(os);
    os << std::endl;
  }

  os << "Stratimikos Preconditioners: " << std::endl;
  os << "********************************" << std::endl;
  for (itr = stratPrecond_.begin(); itr != stratPrecond_.end(); ++itr) {
    os << "name = \"" << itr->first << "\"" << std::endl;
    itr->second->print(os);
    os << std::endl;
  }

  os << "Teko Preconditioners: " << std::endl;
  os << "********************************" << std::endl;
  for (itr = blockPrecond_.begin(); itr != blockPrecond_.end(); ++itr) {
    os << "name = \"" << itr->first << "\"" << std::endl;
    itr->second->print(os);
    os << std::endl;
  }
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
RCP<InverseLibrary> InverseLibrary::buildFromParameterList(const Teuchos::ParameterList& pl,
                                                           bool useStratDefaults) {
  // build from Stratimikos or allocate a new inverse library
  RCP<InverseLibrary> invLib;
  if (useStratDefaults)
    invLib = InverseLibrary::buildFromStratimikos();
  else
    invLib = rcp(new InverseLibrary());

  // to convert the void* like entry
  Teuchos::ParameterList* temp = 0;

  // loop over all entries in parameter list
  Teuchos::ParameterList::ConstIterator itr;
  for (itr = pl.begin(); itr != pl.end(); ++itr) {
    // get current entry
    std::string label            = itr->first;
    Teuchos::ParameterList& list = itr->second.getValue(temp);

    // add to library
    invLib->addInverse(label, list);
  }

  return invLib;
}

/** \brief Build an inverse library from a parameter list.
 *
 * Build an inverse library from a parameter list. This will
 * contain all the labeled inverses specified.
 *
 * \param[in] pl Parameter list to build the library from
 * \param[in] strat Stratimikos object to use
 *
 * \returns A pointer to the inverse library created.
 */
RCP<InverseLibrary> InverseLibrary::buildFromParameterList(
    const Teuchos::ParameterList& pl,
    const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder>& strat) {
  // if strat is set to null, use the defaults
  if (strat == Teuchos::null) return buildFromParameterList(pl, true);

  // build from Stratimikos or allocate a new inverse library
  RCP<InverseLibrary> invLib = InverseLibrary::buildFromStratimikos(strat);

  // to convert the void* like entry
  Teuchos::ParameterList* temp = 0;

  // loop over all entries in parameter list
  Teuchos::ParameterList::ConstIterator itr;
  for (itr = pl.begin(); itr != pl.end(); ++itr) {
    // get current entry
    std::string label            = itr->first;
    Teuchos::ParameterList& list = itr->second.getValue(temp);

    // add to library
    invLib->addInverse(label, list);
  }

  return invLib;
}

/** \brief Build an inverse library from Stratimikos
 *
 * Build an inverse library from Stratimkos. The labels
 * will just be the names in Stratimikos. Uses the
 * defaultBuilder with extra inverse added to it
 *
 * \returns A pointer to the inverse library created.
 */
Teuchos::RCP<InverseLibrary> InverseLibrary::buildFromStratimikos() {
  RCP<InverseLibrary> invLib = rcp(new InverseLibrary());

  // get default inveres in Stratimikos
  RCP<Stratimikos::DefaultLinearSolverBuilder> strat = invLib->defaultBuilder_;
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList(*strat->getValidParameters()));
  Teuchos::ParameterList lst(pl->sublist("Linear Solver Types"));
  Teuchos::ParameterList pft(pl->sublist("Preconditioner Types"));

  Teuchos::ParameterList::ConstIterator itr;
  Teuchos::ParameterList* temp = 0;

  // loop over all entries in solver list
  for (itr = lst.begin(); itr != lst.end(); ++itr) {
    // get current entry
    std::string label            = itr->first;
    Teuchos::ParameterList& list = itr->second.getValue(temp);
    list.set("Type", label);

    // add to library
    invLib->addInverse(label, list);
  }

  // loop over all entries in preconditioner list
  for (itr = pft.begin(); itr != pft.end(); ++itr) {
    // get current entry
    std::string label            = itr->first;
    Teuchos::ParameterList& list = itr->second.getValue(temp);
    list.set("Type", label);

    // add to library
    invLib->addInverse(label, list);
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
Teuchos::RCP<InverseLibrary> InverseLibrary::buildFromStratimikos(
    const Stratimikos::DefaultLinearSolverBuilder& strat) {
  RCP<InverseLibrary> invLib = rcp(new InverseLibrary());

  // get default inveres in Stratimikos
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList(*strat.getValidParameters()));
  Teuchos::ParameterList lst(pl->sublist("Linear Solver Types"));
  Teuchos::ParameterList pft(pl->sublist("Preconditioner Types"));

  Teuchos::ParameterList::ConstIterator itr;
  Teuchos::ParameterList* temp = 0;

  // loop over all entries in solver list
  for (itr = lst.begin(); itr != lst.end(); ++itr) {
    // get current entry
    std::string label            = itr->first;
    Teuchos::ParameterList& list = itr->second.getValue(temp);
    list.set("Type", label);

    // add to library
    invLib->addInverse(label, list);
  }

  // loop over all entries in preconditioner list
  for (itr = pft.begin(); itr != pft.end(); ++itr) {
    // get current entry
    std::string label            = itr->first;
    Teuchos::ParameterList& list = itr->second.getValue(temp);
    list.set("Type", label);

    // add to library
    invLib->addInverse(label, list);
  }

  return invLib;
}

/** \brief Build an inverse library from Stratimikos
 *
 * Build an inverse library from Stratimkos. The labels
 * will just be the names in Stratimikos.
 *
 * \param[in] strat Stratimikos pointer to use
 *
 * \returns A pointer to the inverse library created.
 */
Teuchos::RCP<InverseLibrary> InverseLibrary::buildFromStratimikos(
    const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder>& strat) {
  RCP<InverseLibrary> invLib = rcp(new InverseLibrary(strat));

  // get default inveres in Stratimikos
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList(*strat->getValidParameters()));
  Teuchos::ParameterList lst(pl->sublist("Linear Solver Types"));
  Teuchos::ParameterList pft(pl->sublist("Preconditioner Types"));

  Teuchos::ParameterList::ConstIterator itr;
  Teuchos::ParameterList* temp = 0;

  // loop over all entries in solver list
  for (itr = lst.begin(); itr != lst.end(); ++itr) {
    // get current entry
    std::string label            = itr->first;
    Teuchos::ParameterList& list = itr->second.getValue(temp);
    list.set("Type", label);

    // add to library
    invLib->addInverse(label, list);
  }

  // loop over all entries in preconditioner list
  for (itr = pft.begin(); itr != pft.end(); ++itr) {
    // get current entry
    std::string label            = itr->first;
    Teuchos::ParameterList& list = itr->second.getValue(temp);
    list.set("Type", label);

    // add to library
    invLib->addInverse(label, list);
  }

  return invLib;
}

}  // end namespace Teko
