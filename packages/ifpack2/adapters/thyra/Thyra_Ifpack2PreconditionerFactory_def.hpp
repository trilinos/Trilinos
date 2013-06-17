/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

#include "Teuchos_Ptr.hpp"
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
  const ESupportSolveUse supportSolveUse
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

  const Teuchos::RCP<const Teuchos::ParameterList> constParamList = paramList_;
  const std::string preconditionerType = Teuchos::getParameter<std::string>(*constParamList, "Prec Type");
  const int overlap = Teuchos::getParameter<int>(*constParamList, "Overlap");
  const Teuchos::RCP<const Teuchos::ParameterList> packageParamList = Teuchos::sublist(constParamList, "Ifpack2 Settings");

  // Create the initial preconditioner

  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    *out << "\nCreating a new Ifpack2::Preconditioner object...\n";
  }
  timer.start(true);

  typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type> Ifpack2Prec;
  const Teuchos::RCP<Ifpack2Prec> concretePrecOp = Ifpack2::Factory::create(preconditionerType, tpetraFwdMatrix, overlap);

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
