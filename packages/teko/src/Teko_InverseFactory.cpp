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

#include "Teko_InverseFactory.hpp"

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
#include "Teko_PreconditionerInverseFactory.hpp"

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

namespace Teko {

//! Build an inverse operator using a factory and a linear operator
InverseLinearOp buildInverse(const InverseFactory& factory, const LinearOp& A) {
  InverseLinearOp inv;
  try {
    inv = factory.buildInverse(A);
  } catch (std::exception& e) {
    RCP<Teuchos::FancyOStream> out = Teko::getOutputStream();

    *out << "Teko: \"buildInverse\" could not construct the inverse operator using ";
    *out << "\"" << factory.toString() << "\"" << std::endl;
    *out << std::endl;
    *out << "*** THROWN EXCEPTION ***\n";
    *out << e.what() << std::endl;
    *out << "************************\n";

    throw e;
  }

  return inv;
}

/** Build an inverse operator using a factory and a linear operator
 *
 * \param[in] factory The inverse factory used to construct the inverse
 *                    operator
 * \param[in] precOp  Preconditioning operator
 * \param[in] A       Linear operator whose inverse is required
 *
 * \returns An (approximate) inverse operator is returned for the operator <code>A</code>.
 *
 * \relates InverseFactory
 */
InverseLinearOp buildInverse(const InverseFactory& factory, const LinearOp& A,
                             const LinearOp& precOp) {
  Teko_DEBUG_SCOPE("buildInverse(factory,A,precOp)", 10);
  InverseLinearOp inv;
  try {
    inv = factory.buildInverse(A, precOp);
  } catch (std::exception& e) {
    RCP<Teuchos::FancyOStream> out = Teko::getOutputStream();

    *out << "Teko: \"buildInverse\" could not construct the inverse operator using ";
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
void rebuildInverse(const InverseFactory& factory, const LinearOp& A, InverseLinearOp& invA) {
  InverseLinearOp inv;
  try {
    factory.rebuildInverse(A, invA);
  } catch (std::exception& e) {
    RCP<Teuchos::FancyOStream> out = Teko::getOutputStream();

    *out << "Teko: \"rebuildInverse\" could not construct the inverse operator using ";
    *out << "\"" << factory.toString() << "\"" << std::endl;
    *out << std::endl;
    *out << "*** THROWN EXCEPTION ***\n";
    *out << e.what() << std::endl;
    *out << "************************\n";

    throw e;
  }
}

/** Using a prebuilt linear operator, use factory to build an inverse operator
 * given a new forward operator.
 *
 * \note This function sometimes fails depending on the underlying type
 *       of the inverse factory.  Use with caution.
 *
 * \param[in] factory The inverse factory used to construct the inverse
 *                    operator
 * \param[in] A       Linear operator whose inverse is required
 * \param[in] precOp  Preconditioning operator
 * \param[in] invA    The inverse operator that is to be rebuilt using
 *                    the <code>A</code> operator.
 *
 * \relates InverseFactory
 */
void rebuildInverse(const InverseFactory& factory, const LinearOp& A, const LinearOp& precOp,
                    InverseLinearOp& invA) {
  InverseLinearOp inv;
  try {
    factory.rebuildInverse(A, precOp, invA);
  } catch (std::exception& e) {
    RCP<Teuchos::FancyOStream> out = Teko::getOutputStream();

    *out << "Teko: \"rebuildInverse\" could not construct the inverse operator using ";
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
RCP<InverseFactory> invFactoryFromParamList(const Teuchos::ParameterList& list,
                                            const std::string& type) {
  RCP<Teuchos::ParameterList> myList = rcp(new Teuchos::ParameterList(list));

  Stratimikos::DefaultLinearSolverBuilder strat;
  addToStratimikosBuilder(rcpFromRef(strat));
  strat.setParameterList(myList);

  try {
    // try to build a preconditioner factory
    RCP<Thyra::PreconditionerFactoryBase<double> > precFact =
        strat.createPreconditioningStrategy(type);

    // string must map to a preconditioner
    return rcp(new PreconditionerInverseFactory(precFact, Teuchos::null));
  } catch (const Teuchos::Exceptions::InvalidParameterValue& exp) {
  }

  try {
    // try to build a solver factory
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solveFact =
        strat.createLinearSolveStrategy(type);

    // if its around, build a InverseFactory
    return rcp(new SolveInverseFactory(solveFact));
  } catch (const Teuchos::Exceptions::InvalidParameterValue& exp) {
  }

  return Teuchos::null;
  ;
}

/** \brief Get a valid parameter list for the inverse factory class.
 *
 * Get a valid parameter list for the inverse factory class. This will
 * specify the set of parameters for each possible "inverse".
 *
 * \returns A parameter list is returned that is suitable to be passed
 *          to <code>invFactoryFromParamList</code>.
 */
Teuchos::RCP<const Teuchos::ParameterList> invFactoryValidParameters() {
  Stratimikos::DefaultLinearSolverBuilder strat;

  // extract valid parameter list from Stratimikos
  return strat.getValidParameters();
}

}  // end namespace Teko
