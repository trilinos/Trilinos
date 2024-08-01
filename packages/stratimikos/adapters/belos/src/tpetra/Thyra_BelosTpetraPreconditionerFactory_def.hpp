// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DEF_HPP
#define THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DEF_HPP

#include "Thyra_BelosTpetraPreconditionerFactory_decl.hpp"

#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include "BelosTpetraOperator.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"

#include "Tpetra_MixedScalarMultiplyOp.hpp"

#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include <string>


// CAG: This is not entirely ideal, since it disables half-precision
//      altogether when we have e.g. double, float and
//      complex<double>, but not complex<float>, although a
//      double-float preconditioner would be possible.
#if (!defined(HAVE_TPETRA_INST_DOUBLE) || (defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT))) && \
    (!defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) || (defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)))
# define THYRA_BELOS_PREC_ENABLE_HALF_PRECISION
#endif

namespace Thyra {

// Constructors/initializers/accessors


template <typename MatrixType>
BelosTpetraPreconditionerFactory<MatrixType>::BelosTpetraPreconditionerFactory()
{}


// Overridden from PreconditionerFactoryBase


template <typename MatrixType>
bool BelosTpetraPreconditionerFactory<MatrixType>::isCompatible(
  const LinearOpSourceBase<scalar_type> &fwdOpSrc
  ) const
{
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;

  const Teuchos::RCP<const LinearOpBase<scalar_type> > fwdOp = fwdOpSrc.getOp();
  using TpetraExtractHelper = TpetraOperatorVectorExtraction<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  const auto tpetraFwdOp =  TpetraExtractHelper::getConstTpetraOperator(fwdOp);

  const Teuchos::RCP<const MatrixType> tpetraFwdMatrix = Teuchos::rcp_dynamic_cast<const MatrixType>(tpetraFwdOp,true);

  return Teuchos::nonnull(tpetraFwdMatrix);
}


template <typename MatrixType>
Teuchos::RCP<PreconditionerBase<typename BelosTpetraPreconditionerFactory<MatrixType>::scalar_type> >
BelosTpetraPreconditionerFactory<MatrixType>::createPrec() const
{
  return Teuchos::rcp(new DefaultPreconditioner<scalar_type>);
}


template <typename MatrixType>
void BelosTpetraPreconditionerFactory<MatrixType>::initializePrec(
  const Teuchos::RCP<const LinearOpSourceBase<scalar_type> > &fwdOpSrc,
  PreconditionerBase<scalar_type> *prec,
  const ESupportSolveUse /* supportSolveUse */
  ) const
{
  using Teuchos::rcp;
  using Teuchos::RCP;

  // Check precondition

  TEUCHOS_ASSERT(Teuchos::nonnull(fwdOpSrc));
  TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
  TEUCHOS_ASSERT(prec);

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
    *out << "\nEntering Thyra::BelosTpetraPreconditionerFactory::initializePrec(...) ...\n";
  }

  // Retrieve wrapped concrete Tpetra matrix from FwdOp

  const RCP<const LinearOpBase<scalar_type> > fwdOp = fwdOpSrc->getOp();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(fwdOp));

  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;

  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> TpetraLinOp;
  using TpetraExtractHelper = TpetraOperatorVectorExtraction<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  const auto tpetraFwdOp =  TpetraExtractHelper::getConstTpetraOperator(fwdOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdOp));

  // Belos-specific typedefs:
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> TpetraMV; 
  typedef Belos::TpetraOperator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> BelosTpOp;
  typedef Belos::LinearProblem<scalar_type, TpetraMV, TpetraLinOp> BelosTpLinProb;

#ifdef THYRA_BELOS_PREC_ENABLE_HALF_PRECISION
  // CAG: There is nothing special about the combination double-float,
  //      except that I feel somewhat confident that Trilinos builds
  //      with both scalar types.
  typedef typename Teuchos::ScalarTraits<scalar_type>::halfPrecision half_scalar_type;
  typedef Tpetra::Operator<half_scalar_type, local_ordinal_type, global_ordinal_type, node_type> TpetraLinOpHalf;
  typedef Tpetra::MultiVector<half_scalar_type, local_ordinal_type, global_ordinal_type, node_type> TpetraMVHalf; 
  typedef Belos::TpetraOperator<half_scalar_type, local_ordinal_type, global_ordinal_type, node_type> BelosTpOpHalf;
  typedef Belos::LinearProblem<half_scalar_type, TpetraMVHalf, TpetraLinOpHalf> BelosTpLinProbHalf;
#endif

  // Retrieve concrete preconditioner object

  const Teuchos::Ptr<DefaultPreconditioner<scalar_type> > defaultPrec =
    Teuchos::ptr(dynamic_cast<DefaultPreconditioner<scalar_type> *>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  // This check needed to address Issue #535. 
  RCP<Teuchos::ParameterList> innerParamList;
  if (paramList_.is_null ()) {
    innerParamList = rcp(new Teuchos::ParameterList(*getValidParameters()));
  }
  else {
    innerParamList = paramList_;
  }

  bool useHalfPrecision = false;
  if (innerParamList->isParameter("half precision"))
    useHalfPrecision = Teuchos::getParameter<bool>(*innerParamList, "half precision");

  const std::string solverType = Teuchos::getParameter<std::string>(*innerParamList, "BelosPrec Solver Type");
  const RCP<Teuchos::ParameterList> packageParamList = Teuchos::rcpFromRef(innerParamList->sublist("BelosPrec Solver Params"));

  // solverTypeUpper is the upper-case version of solverType.
  std::string solverTypeUpper (solverType);
  std::transform(solverTypeUpper.begin(), solverTypeUpper.end(),solverTypeUpper.begin(), ::toupper);

  // Create the initial preconditioner
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
    *out << "\nCreating a new BelosTpetra::Preconditioner object...\n";
  }
  RCP<LinearOpBase<scalar_type> > thyraPrecOp;

  if (useHalfPrecision) {
#ifdef THYRA_BELOS_PREC_ENABLE_HALF_PRECISION
    if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
      Teuchos::OSTab(out).o() << "> Creating half precision preconditioner\n";
    }
    const RCP<const MatrixType> tpetraFwdMatrix = Teuchos::rcp_dynamic_cast<const MatrixType>(tpetraFwdOp);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdMatrix));
    auto tpetraFwdMatrixHalf = tpetraFwdMatrix->template convert<half_scalar_type>();

    RCP<BelosTpLinProbHalf> belosLinProbHalf = rcp(new BelosTpLinProbHalf());
    belosLinProbHalf->setOperator(tpetraFwdMatrixHalf);
    RCP<TpetraLinOpHalf> belosOpRCPHalf = rcp(new BelosTpOpHalf(belosLinProbHalf, packageParamList, solverType, true));
    RCP<TpetraLinOp> wrappedOp = rcp(new Tpetra::MixedScalarMultiplyOp<scalar_type,half_scalar_type,local_ordinal_type,global_ordinal_type,node_type>(belosOpRCPHalf));

    thyraPrecOp = Thyra::createLinearOp(wrappedOp);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Solver does not have correct precisions enabled to use half precision.")
#endif
  } else {
    // Wrap concrete preconditioner
    RCP<BelosTpLinProb> belosLinProb = rcp(new BelosTpLinProb());
    belosLinProb->setOperator(tpetraFwdOp);
    RCP<TpetraLinOp> belosOpRCP = rcp(new BelosTpOp(belosLinProb, packageParamList, solverType, true));
    thyraPrecOp = Thyra::createLinearOp(belosOpRCP);
  }
  defaultPrec->initializeUnspecified(thyraPrecOp);

  totalTimer.stop();
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    *out << "\nTotal time in Thyra::BelosTpetraPreconditionerFactory::initializePrec(...) = " << totalTimer.totalElapsedTime() << " sec\n";
  }

  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
    *out << "\nLeaving Thyra::BelosTpetraPreconditionerFactory::initializePrec(...) ...\n";
  }
}


template <typename MatrixType>
void BelosTpetraPreconditionerFactory<MatrixType>::uninitializePrec(
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
void BelosTpetraPreconditionerFactory<MatrixType>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));

  const auto validParamList = this->getValidParameters();
  paramList->validateParametersAndSetDefaults(*validParamList, 0);
  paramList->validateParameters(*validParamList, 0);

  paramList_ = paramList;
  Teuchos::readVerboseObjectSublist(paramList_.getRawPtr(), this);
}


template <typename MatrixType>
Teuchos::RCP<Teuchos::ParameterList>
BelosTpetraPreconditionerFactory<MatrixType>::getNonconstParameterList()
{
  return paramList_;
}


template <typename MatrixType>
Teuchos::RCP<Teuchos::ParameterList>
BelosTpetraPreconditionerFactory<MatrixType>::unsetParameterList()
{
  const Teuchos::RCP<Teuchos::ParameterList> savedParamList = paramList_;
  paramList_ = Teuchos::null;
  return savedParamList;
}


template <typename MatrixType>
Teuchos::RCP<const Teuchos::ParameterList>
BelosTpetraPreconditionerFactory<MatrixType>::getParameterList() const
{
  return paramList_;
}

template <typename MatrixType>
Teuchos::RCP<const Teuchos::ParameterList>
BelosTpetraPreconditionerFactory<MatrixType>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> validParamList;

  if (Teuchos::is_null(validParamList)) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("BelosPrecTpetra"));
    validParamList->set(
      "BelosPrec Solver Type", "GMRES",
      "Name of Belos solver to be used as a preconditioner. (Use valid names for Belos::SolverFactory.)"
      );
    validParamList->set(
      "half precision", false,
      "Whether a half-of-standard-precision Belos-solver-as-preconditioner should be built."
      );
    validParamList->sublist(
      "BelosPrec Solver Params", false,
      "Belos solver settings that are passed onto the Belos solver itself."
      );
    Teuchos::setupVerboseObjectSublist(validParamList.getRawPtr());
  }

  return validParamList;
}


// Public functions overridden from Teuchos::Describable

template <typename MatrixType>
std::string BelosTpetraPreconditionerFactory<MatrixType>::description() const
{
  return "Thyra::BelosTpetraPreconditionerFactory";
}


} // namespace Thyra

#endif // THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DEF_HPP
