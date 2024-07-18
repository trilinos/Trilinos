// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
// @header
//
// ***********************************************************************
//
//      teko: a package for block and physics based preconditioning
//                  copyright 2010 sandia corporation
//
// under the terms of contract de-ac04-94al85000 with sandia corporation,
// the u.s. government retains certain rights in this software.
//
// redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. neither the name of the corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// this software is provided by sandia corporation "as is" and any
// express or implied warranties, including, but not limited to, the
// implied warranties of merchantability and fitness for a particular
// purpose are disclaimed. in no event shall sandia corporation or the
// contributors be liable for any direct, indirect, incidental, special,
// exemplary, or consequential damages (including, but not limited to,
// procurement of substitute goods or services; loss of use, data, or
// profits; or business interruption) however caused and on any theory of
// liability, whether in contract, strict liability, or tort (including
// negligence or otherwise) arising in any way out of the use of this
// software, even if advised of the possibility of such damage.
//
// questions? contact eric c. cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @header
*/

#include "Teko_PreconditionerInverseFactory.hpp"

// Thyra includes
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Preconditioner.hpp"
#include "Teko_PreconditionerLinearOp.hpp"
#include "Teko_SolveInverseFactory.hpp"

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

namespace Teko {

/** \brief Constructor that takes a Thyra solve factory and
 *        makes it look like an InverseFactory
 *
 * Constructor that takes a Thyra solve factory and
 * makes it look like an InverseFactory.
 *
 * \param[in] precFactory Thyra PreconditionerFactoryBase used for building
 *                        the inverse.
 */
PreconditionerInverseFactory::PreconditionerInverseFactory(
    const Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> >& precFactory,
    const Teuchos::RCP<Teko::RequestHandler>& rh)
    : precFactory_(precFactory) {
  setRequestHandler(rh);
}

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
    const Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> >& precFactory,
    const Teuchos::RCP<const Teuchos::ParameterList>& xtraParam,
    const Teuchos::RCP<Teko::RequestHandler>& rh)
    : precFactory_(precFactory) {
  if (xtraParam != Teuchos::null)
    extraParams_ = rcp(new Teuchos::ParameterList(*xtraParam));
  else
    extraParams_ = Teuchos::null;  // make it explicit

  setRequestHandler(rh);
}

//! Copy constructor
PreconditionerInverseFactory::PreconditionerInverseFactory(
    const PreconditionerInverseFactory& pFactory)
    : precFactory_(pFactory.precFactory_) {
  setRequestHandler(pFactory.getRequestHandler());
}

/** \brief Build an inverse operator
 *
 * Build the inverse operator using this factory.
 *
 * \param[in] linearOp Linear operator needing to be inverted.
 *
 * \returns New linear operator that functions as the inverse
 *          of <code>linearOp</code>.
 */
InverseLinearOp PreconditionerInverseFactory::buildInverse(const LinearOp& linearOp) const {
  RCP<Thyra::PreconditionerBase<double> > prec = precFactory_->createPrec();
  precFactory_->initializePrec(Thyra::defaultLinearOpSource(linearOp), &*prec);

  RCP<Teko::PreconditionerLinearOp<double> > precOp =
      rcp(new Teko::PreconditionerLinearOp<double>(prec));

  return precOp;
}

/** \brief Build an inverse operator and make sure it aware of some parents state
 *        This functionality is only useful for Teko::PreconditionerFactory inverses.
 *
 * Build an inverse operator and make sure it aware of some parents state
 * This functionality is only useful for Teko::PreconditionerFactory inverses.
 *
 * \param[in] linearOp Linear operator needing to be inverted.
 * \param[in] parentState Current state object to be used. Only useful for preconditioners.
 *
 * \returns New linear operator that functions as the inverse
 *          of <code>linearOp</code>.
 */
InverseLinearOp PreconditionerInverseFactory::buildInverse(
    const LinearOp& linearOp, const PreconditionerState& parentState) const {
  Teko_DEBUG_SCOPE("PreconditionerInverseFactory::buildInverse(A,parentState)", 10);
  RCP<Thyra::PreconditionerBase<double> > prec = precFactory_->createPrec();

  {
    Teko_DEBUG_SCOPE("Casting to Teko::Preconditioner", 10);
    // pass state downward if a Teko::Preconditioner object is begin used
    RCP<Teko::Preconditioner> tekoPrec = Teuchos::rcp_dynamic_cast<Teko::Preconditioner>(prec);
    if (tekoPrec != Teuchos::null) {
      Teko_DEBUG_SCOPE("Merging states", 10);
      tekoPrec->mergeStateObject(parentState);
    }
  }

  precFactory_->initializePrec(Thyra::defaultLinearOpSource(linearOp), &*prec);

  RCP<Teko::PreconditionerLinearOp<double> > precOp =
      rcp(new Teko::PreconditionerLinearOp<double>(prec));

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
void PreconditionerInverseFactory::rebuildInverse(const LinearOp& source,
                                                  InverseLinearOp& dest) const {
  Teko_DEBUG_MSG("BEGIN PreconditionerInverseFactory::rebuildInverse", 10);

  RCP<Thyra::PreconditionerBase<double> > prec =
      Teuchos::rcp_dynamic_cast<Teko::PreconditionerLinearOp<double> >(dest)
          ->getNonconstPreconditioner();

  precFactory_->initializePrec(Thyra::defaultLinearOpSource(source), &*prec);

  Teko_DEBUG_MSG("END PreconditionerInverseFactory::rebuildInverse", 10);
}

/** \brief A function that permits inspection of the parameters used to create
 *        this object.
 *
 * A function that permits inspection of the parameters used to create this
 * object. Useful for determining defaults and settings used.
 *
 * \returns A list used to parameterize this object.
 */
Teuchos::RCP<const Teuchos::ParameterList> PreconditionerInverseFactory::getParameterList() const {
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
Teuchos::RCP<Teuchos::ParameterList> PreconditionerInverseFactory::getRequestedParameters() const {
  Teuchos::RCP<BlockPreconditionerFactory> bpf =
      rcp_dynamic_cast<BlockPreconditionerFactory>(precFactory_);

  // request the parameters from a BPF is required
  if (bpf != Teuchos::null) return bpf->getRequestedParameters();

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
bool PreconditionerInverseFactory::updateRequestedParameters(const Teuchos::ParameterList& pl) {
  Teuchos::RCP<BlockPreconditionerFactory> bpf =
      rcp_dynamic_cast<BlockPreconditionerFactory>(precFactory_);

  // update the parameters of a BPF is required
  if (bpf != Teuchos::null) return bpf->updateRequestedParameters(pl);

  // for non block preconditioners see if there are user requested additional parameters
  if (extraParams_ == Teuchos::null) return true;

  Teuchos::ParameterList::ConstIterator itr;
  RCP<Teuchos::ParameterList> srcPl = precFactory_->unsetParameterList();

  // find name of settings sublist
  std::string subName = "";
  for (itr = srcPl->begin(); itr != srcPl->end(); ++itr) {
    // search for std::string with "Settings" in name
    if (itr->first.find("Settings") != std::string::npos) {
      subName = itr->first;
      continue;
    }
  }

  // update fails if no settings list was found
  if (subName == "") {
    precFactory_->setParameterList(srcPl);
    return false;
  }

  // add extra parameters to list
  Teuchos::ParameterList& settingsList = srcPl->sublist(subName);
  for (itr = pl.begin(); itr != pl.end(); ++itr) {
    if (extraParams_->isParameter(itr->first)) settingsList.setEntry(itr->first, itr->second);
  }

  // set the parameter list
  precFactory_->setParameterList(srcPl);

  return true;
}

void PreconditionerInverseFactory::setupParameterListFromRequestHandler() {
  // for non block preconditioners see if there are user requested additional parameters
  if (extraParams_ == Teuchos::null) return;

  Teuchos::ParameterList::ConstIterator itr;
  RCP<Teuchos::ParameterList> srcPl = precFactory_->unsetParameterList();

  // find name of settings sublist
  std::string subName = "";
  for (itr = srcPl->begin(); itr != srcPl->end(); ++itr) {
    // search for std::string with "Settings" in name
    if (itr->first.find("Settings") != std::string::npos) {
      subName = itr->first;
      continue;
    }
  }

  // update fails if no settings list was found
  /*
  if(subName=="") {
     precFactory_->setParameterList(srcPl);
     return;
  }*/

  Teuchos::RCP<Teko::RequestHandler> rh = getRequestHandler();
  TEUCHOS_TEST_FOR_EXCEPTION(
      rh == Teuchos::null, std::runtime_error,
      "PreconditionerInverseFactory::setupParameterListFromRequestHandler: no request handler set");

  // add extra parameters to list
  rh->preRequest<Teuchos::RCP<Teuchos::ParameterList> >(RequestMesg(extraParams_));
  Teuchos::RCP<Teuchos::ParameterList> requestParams =
      rh->request<Teuchos::RCP<Teuchos::ParameterList> >(RequestMesg(extraParams_));

  TEUCHOS_TEST_FOR_EXCEPTION(requestParams == Teuchos::null, std::runtime_error,
                             "User specified request not satisfied!");

  // If there is no Settings sublist, assume that the list itself contains the settings
  if (subName == "") {
    for (itr = requestParams->begin(); itr != requestParams->end(); ++itr)
      srcPl->setEntry(itr->first, itr->second);
  } else {
    Teuchos::ParameterList& settingsList = srcPl->sublist(subName);
    for (itr = requestParams->begin(); itr != requestParams->end(); ++itr)
      settingsList.setEntry(itr->first, itr->second);
  }

  // reset with updated parameter list
  precFactory_->setParameterList(srcPl);
}

}  // namespace Teko
