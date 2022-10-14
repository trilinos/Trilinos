/*@HEADER
// ***********************************************************************
//
//         Stratimikos: Thyra-based strategies for linear solvers
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
// Questions? Contact Jennifer A. Loe (jloe@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/
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
  const Teuchos::RCP<const LinearOpBase<scalar_type> > fwdOp = fwdOpSrc.getOp();

  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;

  typedef Thyra::TpetraLinearOp<scalar_type, local_ordinal_type, global_ordinal_type, node_type> ThyraTpetraLinOp;
  const Teuchos::RCP<const ThyraTpetraLinOp> thyraTpetraFwdOp = Teuchos::rcp_dynamic_cast<const ThyraTpetraLinOp>(fwdOp,true);

  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> TpetraLinOp;
  const Teuchos::RCP<const TpetraLinOp> tpetraFwdOp = Teuchos::nonnull(thyraTpetraFwdOp) ? thyraTpetraFwdOp->getConstTpetraOperator() : Teuchos::null;

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
    *out << "\nEntering Thyra::BelosTpetraPreconditionerFactory::initializePrec(...) ...\n";
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

  std::cout << "Retrieved wrapped Tpetra op... line 167" << std::endl;

  // Belos-specific typedefs:
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>     TpetraMV; 
  typedef Belos::TpetraOperator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> BelosTpOp;
  typedef Belos::LinearProblem<scalar_type, TpetraMV, TpetraLinOp> BelosTpLinProb;

#ifdef THYRA_BELOS_PREC_ENABLE_HALF_PRECISION
  // CAG: There is nothing special about the combination double-float,
  //      except that I feel somewhat confident that Trilinos builds
  //      with both scalar types.
  typedef typename Teuchos::ScalarTraits<scalar_type>::halfPrecision half_scalar_type;
  typedef Tpetra::Operator<half_scalar_type, local_ordinal_type, global_ordinal_type, node_type> TpetraLinOpHalf;
  typedef Tpetra::MultiVector<half_scalar_type, local_ordinal_type, global_ordinal_type, node_type>     TpetraMVHalf; 
  typedef Belos::TpetraOperator<half_scalar_type, local_ordinal_type, global_ordinal_type, node_type> BelosTpOpHalf;
  typedef Belos::LinearProblem<half_scalar_type, TpetraMVHalf, TpetraLinOpHalf> BelosTpLinProbHalf;
#endif

  // Retrieve concrete preconditioner object

  const Teuchos::Ptr<DefaultPreconditioner<scalar_type> > defaultPrec =
    Teuchos::ptr(dynamic_cast<DefaultPreconditioner<scalar_type> *>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  // Process parameter list

  bool useHalfPrecision = false;
  if (paramList_->isParameter("half precision"))
    useHalfPrecision = paramList_->get<bool>("half precision");

  Teuchos::RCP<const Teuchos::ParameterList> constParamList = paramList_;
  if (constParamList.is_null ()) {
    constParamList = getValidParameters ();
  }
  const std::string solverType = Teuchos::getParameter<std::string>(*constParamList, "BelosPrec Solver Type");
  const Teuchos::RCP<const Teuchos::ParameterList> packageParamList = Teuchos::sublist(constParamList, "BelosPrec Solver Params");
  Teuchos::RCP<Teuchos::ParameterList> nonconstPackageParamList = Teuchos::sublist(paramList_, "BelosPrec Solver Params");

  // solverTypeUpper is the upper-case version of solverType.
  std::string solverTypeUpper (solverType);
  if (solverTypeUpper.size () > 0) {
    for (size_t k = 0; k < solverTypeUpper.size (); ++k) {
      solverTypeUpper[k] = ::toupper(solverTypeUpper[k]);
    }
  }
   std::cout << "Processed parameter list... line 208" << std::endl;
  // Create the initial preconditioner

  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    *out << "\nCreating a new BelosTpetra::Preconditioner object...\n";
  }
  Teuchos::RCP<LinearOpBase<scalar_type> > thyraPrecOp;

  std::cout << "just set operator in belosLin Prob line 220." << std::endl; 
  
  //typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type> Ifpack2Prec;
  //Teuchos::RCP<Ifpack2Prec> concretePrecOp;
  std::cout << "Made thyraPrecOp line 229" << std::endl;

  if (useHalfPrecision) {
#ifdef THYRA_BELOS_PREC_ENABLE_HALF_PRECISION
    if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
      Teuchos::OSTab(out).o() << "> Creating half precision preconditioner\n";
    }
    std::cout << "Creating half precision prec" << std::endl;
    const Teuchos::RCP<const MatrixType> tpetraFwdMatrix = Teuchos::rcp_dynamic_cast<const MatrixType>(tpetraFwdOp);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdMatrix));
    auto tpetraFwdMatrixHalf = tpetraFwdMatrix->template convert<half_scalar_type>();

    Teuchos::RCP<BelosTpLinProbHalf> belosLinProbHalf = Teuchos::RCP<BelosTpLinProbHalf>(new BelosTpLinProbHalf());
    belosLinProbHalf->setOperator(tpetraFwdMatrixHalf);
    Teuchos::RCP<TpetraLinOpHalf> belosOpRCPHalf = Teuchos::RCP<TpetraLinOpHalf>(new BelosTpOpHalf(belosLinProbHalf, nonconstPackageParamList, solverType, true));
    RCP<TpetraLinOp> wrappedOp = rcp(new Tpetra::MixedScalarMultiplyOp<scalar_type,half_scalar_type,local_ordinal_type,global_ordinal_type,node_type>(belosOpRCPHalf));

    thyraPrecOp = Thyra::createLinearOp(wrappedOp);
#else
    if (Teuchos::nonnull(out)) {
      *out << "\nSolver does not have correct precisions enabled to use half precision.\n";
    }
    TEUCHOS_TEST_FOR_EXCEPT(true);
#endif
  } else {
    //concretePrecOp->setParameters(*packageParamList);

    // Wrap concrete preconditioner
    //thyraPrecOp = Thyra::createLinearOp(Teuchos::RCP<TpetraLinOp>(concretePrecOp));
    Teuchos::RCP<BelosTpLinProb> belosLinProb = Teuchos::RCP<BelosTpLinProb>(new BelosTpLinProb());
    belosLinProb->setOperator(tpetraFwdOp);
    Teuchos::RCP<TpetraLinOp> belosOpRCP = Teuchos::RCP<TpetraLinOp>(new BelosTpOp(belosLinProb, nonconstPackageParamList, solverType, true));
    thyraPrecOp = Thyra::createLinearOp(belosOpRCP);
    //TODO Update Belos TpetraOperator notes to make it clearer why this is true, not default false. 
  }
  std::cout << "Ready to init default Prec line 280." << std::endl;
  defaultPrec->initializeUnspecified(thyraPrecOp);

  totalTimer.stop();
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    *out << "\nTotal time in Thyra::Ifpack2PreconditionerFactory::initializePrec(...) = " << totalTimer.totalElapsedTime() << " sec\n";
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

  const Teuchos::RCP<const Teuchos::ParameterList> validParamList = this->getValidParameters();
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
    Teuchos::ParameterList &packageParamList = validParamList->sublist(
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
