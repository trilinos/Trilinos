/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/
#ifndef THYRA_IFPACK2_PRECONDITIONERFACTORY_DEF_HPP
#define THYRA_IFPACK2_PRECONDITIONERFACTORY_DEF_HPP

#include "Thyra_Ifpack2PreconditionerFactory_decl.hpp"

#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include "Ifpack2_Factory.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Parameters.hpp"

#include "Tpetra_RowMatrix.hpp"

#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include <string>


namespace Thyra {


// Constructors/initializers/accessors


template <typename MatrixType>
Ifpack2PreconditionerFactory<MatrixType>::Ifpack2PreconditionerFactory()
{}


// Overridden from PreconditionerFactoryBase


template <typename MatrixType>
bool Ifpack2PreconditionerFactory<MatrixType>::isCompatible(
  const LinearOpSourceBase<scalar_type> &fwdOpSrc
  ) const
{
  const Teuchos::RCP<const LinearOpBase<scalar_type> > fwdOp = fwdOpSrc.getOp();

  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;

  typedef Thyra::TpetraLinearOp<scalar_type, local_ordinal_type, global_ordinal_type, node_type> ThyraTpetraLinOp;
  const Teuchos::RCP<const ThyraTpetraLinOp> thyraTpetraFwdOp = Teuchos::rcp_dynamic_cast<const ThyraTpetraLinOp>(fwdOp);

  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> TpetraLinOp;
  const Teuchos::RCP<const TpetraLinOp> tpetraFwdOp = Teuchos::nonnull(thyraTpetraFwdOp) ? thyraTpetraFwdOp->getConstTpetraOperator() : Teuchos::null;

  const Teuchos::RCP<const MatrixType> tpetraFwdMatrix = Teuchos::rcp_dynamic_cast<const MatrixType>(tpetraFwdOp);

  return Teuchos::nonnull(tpetraFwdMatrix);
}


template <typename MatrixType>
Teuchos::RCP<PreconditionerBase<typename Ifpack2PreconditionerFactory<MatrixType>::scalar_type> >
Ifpack2PreconditionerFactory<MatrixType>::createPrec() const
{
  return Teuchos::rcp(new DefaultPreconditioner<scalar_type>);
}


template <typename MatrixType>
void Ifpack2PreconditionerFactory<MatrixType>::initializePrec(
  const Teuchos::RCP<const LinearOpSourceBase<scalar_type> > &fwdOpSrc,
  PreconditionerBase<scalar_type> *prec,
  const ESupportSolveUse /* supportSolveUse */
  ) const
{
  // Check precondition

  TEUCHOS_ASSERT(Teuchos::nonnull(fwdOpSrc));
  TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
  TEUCHOS_ASSERT(prec);

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
    *out << "\nEntering Thyra::Ifpack2PreconditionerFactory::initializePrec(...) ...\n";
  }

  // Retrieve wrapped concrete Tpetra matrix from FwdOp

  const Teuchos::RCP<const LinearOpBase<scalar_type> > fwdOp = fwdOpSrc->getOp();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(fwdOp));

  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;

  typedef Thyra::TpetraLinearOp<scalar_type, local_ordinal_type, global_ordinal_type, node_type> ThyraTpetraLinOp;
  const Teuchos::RCP<const ThyraTpetraLinOp> thyraTpetraFwdOp = Teuchos::rcp_dynamic_cast<const ThyraTpetraLinOp>(fwdOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraTpetraFwdOp));

  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> TpetraLinOp;
  const Teuchos::RCP<const TpetraLinOp> tpetraFwdOp = thyraTpetraFwdOp->getConstTpetraOperator();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdOp));

  const Teuchos::RCP<const MatrixType> tpetraFwdMatrix = Teuchos::rcp_dynamic_cast<const MatrixType>(tpetraFwdOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdMatrix));

  // Retrieve concrete preconditioner object

  const Teuchos::Ptr<DefaultPreconditioner<scalar_type> > defaultPrec =
    Teuchos::ptr(dynamic_cast<DefaultPreconditioner<scalar_type> *>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  // Process parameter list

  Teuchos::RCP<const Teuchos::ParameterList> constParamList = paramList_;
  if (constParamList.is_null ()) {
    constParamList = getValidParameters ();
  }
  const std::string preconditionerType = Teuchos::getParameter<std::string>(*constParamList, "Prec Type");
  const Teuchos::RCP<const Teuchos::ParameterList> packageParamList = Teuchos::sublist(constParamList, "Ifpack2 Settings");

  // precTypeUpper is the upper-case version of preconditionerType.
  std::string precTypeUpper (preconditionerType);
  if (precTypeUpper.size () > 0) {
    std::locale locale;
    for (size_t k = 0; k < precTypeUpper.size (); ++k) {
      precTypeUpper[k] = std::toupper<char> (precTypeUpper[k], locale);
    }
  }
  
  // mfh 09 Nov 2013: If the Ifpack2 list doesn't already have the
  // "schwarz: overlap level" parameter, then override it with the
  // value of "Overlap".  This avoids use of the newly deprecated
  // three-argument version of Ifpack2::Factory::create() that takes
  // the overlap as an integer.
  if (constParamList->isType<int> ("Overlap") && ! packageParamList.is_null () && ! packageParamList->isType<int> ("schwarz: overlap level") &&
      precTypeUpper == "SCHWARZ") {
    const int overlap = constParamList->get<int> ("Overlap");
    Teuchos::RCP<Teuchos::ParameterList> nonconstPackageParamList =
      Teuchos::sublist (paramList_, "Ifpack2 Settings");
    nonconstPackageParamList->set ("schwarz: overlap level", overlap);
  }

  // Create the initial preconditioner

  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    *out << "\nCreating a new Ifpack2::Preconditioner object...\n";
  }
  timer.start(true);

  typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type> Ifpack2Prec;
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> row_matrix_type;
  const Teuchos::RCP<Ifpack2Prec> concretePrecOp =
    Ifpack2::Factory::create<row_matrix_type> (preconditionerType, tpetraFwdMatrix);

  timer.stop();
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    Teuchos::OSTab(out).o() << "> Creation time = " << timer.totalElapsedTime() << " sec\n";
  }

  // Initialize and compute the initial preconditioner

  concretePrecOp->setParameters(*packageParamList);
  concretePrecOp->initialize();
  concretePrecOp->compute();

  // Wrap concrete preconditioner

  const Teuchos::RCP<LinearOpBase<scalar_type> > thyraPrecOp = Thyra::createLinearOp(Teuchos::RCP<TpetraLinOp>(concretePrecOp));
  defaultPrec->initializeUnspecified(thyraPrecOp);

  totalTimer.stop();
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    *out << "\nTotal time in Thyra::Ifpack2PreconditionerFactory::initializePrec(...) = " << totalTimer.totalElapsedTime() << " sec\n";
  }

  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
    *out << "\nLeaving Thyra::Ifpack2PreconditionerFactory::initializePrec(...) ...\n";
  }
}


template <typename MatrixType>
void Ifpack2PreconditionerFactory<MatrixType>::uninitializePrec(
  PreconditionerBase<scalar_type> *prec,
  Teuchos::RCP<const LinearOpSourceBase<scalar_type> > *fwdOp,
  ESupportSolveUse *supportSolveUse
  ) const
{
  // Check precondition

  TEUCHOS_ASSERT(prec);

  // Retrieve concrete preconditioner object

  const Teuchos::Ptr<DefaultPreconditioner<scalar_type> > defaultPrec =
    Teuchos::ptr(dynamic_cast<DefaultPreconditioner<scalar_type> *>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  if (fwdOp) {
    // TODO: Implement properly instead of returning default value
    *fwdOp = Teuchos::null;
  }

  if (supportSolveUse) {
    // TODO: Implement properly instead of returning default value
    *supportSolveUse = Thyra::SUPPORT_SOLVE_UNSPECIFIED;
  }

  defaultPrec->uninitialize();
}


// Overridden from ParameterListAcceptor


template <typename MatrixType>
void Ifpack2PreconditionerFactory<MatrixType>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));

  const Teuchos::RCP<const Teuchos::ParameterList> validParamList = this->getValidParameters();
  paramList->validateParametersAndSetDefaults(*validParamList, 0);
  paramList->validateParameters(*validParamList, 1);

  paramList_ = paramList;
  Teuchos::readVerboseObjectSublist(paramList_.getRawPtr(), this);
}


template <typename MatrixType>
Teuchos::RCP<Teuchos::ParameterList>
Ifpack2PreconditionerFactory<MatrixType>::getNonconstParameterList()
{
  return paramList_;
}


template <typename MatrixType>
Teuchos::RCP<Teuchos::ParameterList>
Ifpack2PreconditionerFactory<MatrixType>::unsetParameterList()
{
  const Teuchos::RCP<Teuchos::ParameterList> savedParamList = paramList_;
  paramList_ = Teuchos::null;
  return savedParamList;
}


template <typename MatrixType>
Teuchos::RCP<const Teuchos::ParameterList>
Ifpack2PreconditionerFactory<MatrixType>::getParameterList() const
{
  return paramList_;
}

template <typename MatrixType>
Teuchos::RCP<const Teuchos::ParameterList>
Ifpack2PreconditionerFactory<MatrixType>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> validParamList;

  if (Teuchos::is_null(validParamList)) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("Ifpack2"));

    validParamList->set(
      "Prec Type", "ILUT",
      "Type of Ifpack2 preconditioner to use."
      );
    validParamList->set(
      "Overlap", 0,
      "Number of rows/columns overlapped between subdomains in different"
      "\nprocesses in the additive Schwarz-type domain-decomposition preconditioners."
      );
    Teuchos::ParameterList &packageParamList = validParamList->sublist(
      "Ifpack2 Settings", false,
      "Preconditioner settings that are passed onto the Ifpack preconditioners themselves."
      );
    Ifpack2::getValidParameters(packageParamList);

    Teuchos::setupVerboseObjectSublist(validParamList.getRawPtr());
  }

  return validParamList;
}


// Public functions overridden from Teuchos::Describable

template <typename MatrixType>
std::string Ifpack2PreconditionerFactory<MatrixType>::description() const
{
  return "Thyra::Ifpack2PreconditionerFactory";
}


} // namespace Thyra

#endif // THYRA_IFPACK2_PRECONDITIONERFACTORY_DEF_HPP
