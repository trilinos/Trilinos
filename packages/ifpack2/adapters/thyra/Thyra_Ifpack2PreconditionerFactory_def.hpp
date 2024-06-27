// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
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
# define THYRA_IFPACK2_ENABLE_HALF_PRECISION
#endif

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
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;

  const Teuchos::RCP<const LinearOpBase<scalar_type> > fwdOp = fwdOpSrc.getOp();
  using TpetraExtractHelper = TpetraOperatorVectorExtraction<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  const auto tpetraFwdOp =  TpetraExtractHelper::getConstTpetraOperator(fwdOp);

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
    *out << "\nEntering Thyra::Ifpack2PreconditionerFactory::initializePrec(...) ...\n";
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

  const RCP<const MatrixType> tpetraFwdMatrix = Teuchos::rcp_dynamic_cast<const MatrixType>(tpetraFwdOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdMatrix));

  // Retrieve concrete preconditioner object

  const Teuchos::Ptr<DefaultPreconditioner<scalar_type> > defaultPrec =
    Teuchos::ptr(dynamic_cast<DefaultPreconditioner<scalar_type> *>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  // This check needed to address Issue #535. 
  // Corresponding fix exists in stratimikos/adapters/belos/tpetra/Thyra_BelosTpetraPreconditionerFactory_def.hpp
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

  const std::string preconditionerType = Teuchos::getParameter<std::string>(*innerParamList, "Prec Type");
  const RCP<Teuchos::ParameterList> packageParamList = Teuchos::rcpFromRef(innerParamList->sublist("Ifpack2 Settings"));

  // precTypeUpper is the upper-case version of preconditionerType.
  std::string precTypeUpper (preconditionerType);
  std::transform(precTypeUpper.begin(), precTypeUpper.end(),precTypeUpper.begin(), ::toupper);
  
  // mfh 09 Nov 2013: If the Ifpack2 list doesn't already have the
  // "schwarz: overlap level" parameter, then override it with the
  // value of "Overlap".  This avoids use of the newly deprecated
  // three-argument version of Ifpack2::Factory::create() that takes
  // the overlap as an integer.
  if (innerParamList->isType<int> ("Overlap") && ! packageParamList.is_null () && ! packageParamList->isType<int> ("schwarz: overlap level") &&
      precTypeUpper == "SCHWARZ") {
    const int overlap = innerParamList->get<int> ("Overlap");
    packageParamList->set ("schwarz: overlap level", overlap);
  }

  // Create the initial preconditioner

  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
    *out << "\nCreating a new Ifpack2::Preconditioner object...\n";
  }
  timer.start(true);

  RCP<LinearOpBase<scalar_type> > thyraPrecOp;

  typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type> Ifpack2Prec;
  RCP<Ifpack2Prec> concretePrecOp;

#ifdef THYRA_IFPACK2_ENABLE_HALF_PRECISION
  // CAG: There is nothing special about the combination double-float,
  //      except that I feel somewhat confident that Trilinos builds
  //      with both scalar types.
  typedef typename Teuchos::ScalarTraits<scalar_type>::halfPrecision half_scalar_type;
  typedef Ifpack2::Preconditioner<half_scalar_type, local_ordinal_type, global_ordinal_type, node_type> HalfIfpack2Prec;
  RCP<HalfIfpack2Prec> concretePrecOpHalf;
#endif

  if (useHalfPrecision) {
#ifdef THYRA_IFPACK2_ENABLE_HALF_PRECISION
    if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
      Teuchos::OSTab(out).o() << "> Creating half precision preconditioner\n";
    }


    typedef Tpetra::RowMatrix<half_scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> row_matrix_type;
    auto tpetraFwdMatrixHalf = tpetraFwdMatrix->template convert<half_scalar_type>();
    concretePrecOpHalf =
      Ifpack2::Factory::create<row_matrix_type> (preconditionerType, tpetraFwdMatrixHalf);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Ifpack2 does not have correct precisions enabled to use half precision.")
#endif
  } else {
    typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> row_matrix_type;
    concretePrecOp =
      Ifpack2::Factory::create<row_matrix_type> (preconditionerType, tpetraFwdMatrix);
  }

  timer.stop();
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    Teuchos::OSTab(out).o() << "> Creation time = " << timer.totalElapsedTime() << " sec\n";
  }

  // Initialize and compute the initial preconditioner

  if (useHalfPrecision) {
#ifdef THYRA_IFPACK2_ENABLE_HALF_PRECISION
    concretePrecOpHalf->setParameters(*packageParamList);
    concretePrecOpHalf->initialize();
    concretePrecOpHalf->compute();

    RCP<TpetraLinOp> wrappedOp = rcp(new Tpetra::MixedScalarMultiplyOp<scalar_type,half_scalar_type,local_ordinal_type,global_ordinal_type,node_type>(concretePrecOpHalf));

    thyraPrecOp = Thyra::createLinearOp(wrappedOp);
#endif
  } else {
    concretePrecOp->setParameters(*packageParamList);
    concretePrecOp->initialize();
    concretePrecOp->compute();

    // Wrap concrete preconditioner
    thyraPrecOp = Thyra::createLinearOp(RCP<TpetraLinOp>(concretePrecOp));
  }

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

  const auto validParamList = this->getValidParameters();
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
    validParamList->set(
      "half precision", false,
      "Whether a half-precision preconditioner should be built."
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
  std::ostringstream oss;
  oss << "Thyra::Ifpack2PreconditionerFactory";
  oss << "{";
  std::string *precTypePtr = nullptr;
  if (nonnull(paramList_) &&
    (precTypePtr = Teuchos::getParameterPtr<std::string>(*paramList_, "Prec Type")) )
  {
    oss << "\"Prec Type\" = \"" << *precTypePtr << "\"";
  }
  oss << "}";
  return oss.str();
}


} // namespace Thyra

#endif // THYRA_IFPACK2_PRECONDITIONERFACTORY_DEF_HPP
