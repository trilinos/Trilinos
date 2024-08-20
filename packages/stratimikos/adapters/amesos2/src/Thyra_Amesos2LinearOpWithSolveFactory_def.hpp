// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#include "Thyra_Amesos2LinearOpWithSolveFactory_decl.hpp"

#include "Thyra_Amesos2LinearOpWithSolve.hpp"
#include "Amesos2.hpp"
#include "Amesos2_Details_LinearSolverFactory.hpp"
#include "Amesos2_Version.hpp"
#include "Amesos2_Factory.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Thyra {


// Parameter names for Paramter List

template<typename Scalar>
const std::string Amesos2LinearOpWithSolveFactory<Scalar>::SolverType_name
  = "Solver Type";

template<typename Scalar>
const std::string Amesos2LinearOpWithSolveFactory<Scalar>::RefactorizationPolicy_name
  = "Refactorization Policy";

template<typename Scalar>
const std::string Amesos2LinearOpWithSolveFactory<Scalar>::ThrowOnPreconditionerInput_name
  = "Throw on Preconditioner Input";

template<typename Scalar>
const std::string Amesos2LinearOpWithSolveFactory<Scalar>::Amesos2_Settings_name
  = "Amesos2 Settings";

// Constructors/initializers/accessors

template<typename Scalar>
Amesos2LinearOpWithSolveFactory<Scalar>::~Amesos2LinearOpWithSolveFactory()
{
#ifdef TEUCHOS_DEBUG
  if(paramList_.get())
    paramList_->validateParameters(
      *this->getValidParameters(),0  // Only validate this level for now!
      );
#endif
}

template<typename Scalar>
Amesos2LinearOpWithSolveFactory<Scalar>::Amesos2LinearOpWithSolveFactory(
  const Amesos2::ESolverType solverType,
  const Amesos2::ERefactorizationPolicy refactorizationPolicy,
  const bool throwOnPrecInput
  )
  :solverType_(solverType)
  ,refactorizationPolicy_(refactorizationPolicy)
  ,throwOnPrecInput_(throwOnPrecInput)
{
}

// Overridden from LinearOpWithSolveFactoryBase

template<typename Scalar>
bool Amesos2LinearOpWithSolveFactory<Scalar>::isCompatible(
  const LinearOpSourceBase<Scalar> &fwdOpSrc
  ) const
{
  Teuchos::RCP<const LinearOpBase<Scalar> >
    fwdOp = fwdOpSrc.getOp();
  auto tpetraFwdOp = ConverterT::getConstTpetraOperator(fwdOp);
  if ( ! dynamic_cast<const MAT * >(&*tpetraFwdOp) )
    return false;
  return true;
}

template<typename Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
Amesos2LinearOpWithSolveFactory<Scalar>::createOp() const
{
  return Teuchos::rcp(new Amesos2LinearOpWithSolve<Scalar>());
}

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::initializeOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse /* supportSolveUse */
  ) const
{
  THYRA_FUNC_TIME_MONITOR("Stratimikos: Amesos2LOWSF");

  TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(fwdOpSrc.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(fwdOpSrc->getOp().get()==NULL);
  RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc->getOp();

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  //
  // Unwrap and get the forward Tpetra::Operator object
  //
  auto tpetraFwdOp = ConverterT::getConstTpetraOperator(fwdOp);
  auto tpetraCrsMat = Teuchos::rcp_dynamic_cast<const MAT>(tpetraFwdOp);
  // Get the Amesos2LinearOpWithSolve object
  Amesos2LinearOpWithSolve<Scalar>
    *amesos2Op = &Teuchos::dyn_cast<Amesos2LinearOpWithSolve<Scalar>>(*Op);

  //
  // Determine if we must start over or not
  //
  bool startOver = ( amesos2Op->get_amesos2Solver()==Teuchos::null );
  if (!startOver) {
    auto oldTpetraFwdOp = ConverterT::getConstTpetraOperator(amesos2Op->get_fwdOp());
    startOver =
      (
       tpetraFwdOp.get() != oldTpetraFwdOp.get()
       // Assuming that, like Amesos, Amesos2 must start over if the matrix changes
      );
  }
  //
  // Update the amesos2 solver
  //
  if (startOver) {
    //
    // This LOWS object has not be initialized yet or is not compatible with the existing
    //
    // so this is where we setup everything from the ground up.

    // Create the concrete solver
    Teuchos::RCP<Solver> amesos2Solver;
    {
      THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: Amesos2LOWSF:InitConstruct",
        InitConstruct);
      switch(solverType_) {
        case Thyra::Amesos2::KLU2:
          amesos2Solver = ::Amesos2::create<MAT,MV>("klu2", tpetraCrsMat);
          break;
#ifdef HAVE_AMESOS2_LAPACK
        case Thyra::Amesos2::LAPACK:
          amesos2Solver = ::Amesos2::create<MAT,MV>("lapack", tpetraCrsMat);
          break;
#endif
#ifdef HAVE_AMESOS2_SUPERLU
        case Thyra::Amesos2::SUPERLU:
          amesos2Solver = ::Amesos2::create<MAT,MV>("superlu", tpetraCrsMat);
          break;
#endif
#ifdef HAVE_AMESOS2_SUPERLUMT
        case Thyra::Amesos2::SUPERLUMT:
          amesos2Solver = ::Amesos2::create<MAT,MV>("superlumt", tpetraCrsMat);
          break;
#endif
#ifdef HAVE_AMESOS2_SUPERLUDIST
        case Thyra::Amesos2::SUPERLUDIST:
          amesos2Solver = ::Amesos2::create<MAT,MV>("superludist", tpetraCrsMat);
          break;
#  endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
        case Thyra::Amesos2::PARDISO_MKL:
          amesos2Solver = ::Amesos2::create<MAT,MV>("pardiso_mkl", tpetraCrsMat);
          break;
#endif
#ifdef HAVE_AMESOS2_CHOLMOD
        case Thyra::Amesos2::CHOLMOD:
          amesos2Solver = ::Amesos2::create<MAT,MV>("cholmod", tpetraCrsMat);
          break;
#endif
#ifdef HAVE_AMESOS2_BASKER
        case Thyra::Amesos2::BASKER:
          amesos2Solver = ::Amesos2::create<MAT,MV>("basker", tpetraCrsMat);
          break;
#endif
#ifdef HAVE_AMESOS2_MUMPS
        case Thyra::Amesos2::MUMPS:
          amesos2Solver = ::Amesos2::create<MAT,MV>("mumps", tpetraCrsMat);
          break;
#endif
          default:
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error
              ,"Error, the solver type ID = " << solverType_ << " is invalid!"
              );
      }
    }

    // Do the initial factorization
    {
      THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: Amesos2LOWSF:Symbolic", Symbolic);
      amesos2Solver->symbolicFactorization();
    }
    {
      THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: Amesos2LOWSF:Factor", Factor);
      amesos2Solver->numericFactorization();
    }

    // filter out the Stratimikos adapter parameters and hand
    // parameters down into the Solver
    const Teuchos::RCP<Teuchos::ParameterList> dup_list
      = Teuchos::rcp(new Teuchos::ParameterList(*paramList_));
    dup_list->remove(SolverType_name);
    dup_list->remove(RefactorizationPolicy_name);
    dup_list->remove(ThrowOnPreconditionerInput_name);
    dup_list->remove("VerboseObject");
    amesos2Solver->setParameters(dup_list);

    // Initialize the LOWS object and we are done!
    amesos2Op->initialize(fwdOp,fwdOpSrc,amesos2Solver);
  }
  else {
    //
    // This LOWS object has already be initialized once so we must just reset
    // the matrix and refactor it.
    auto amesos2Solver = amesos2Op->get_amesos2Solver();

    // set
    amesos2Solver->setA(tpetraCrsMat);

    // Do the initial factorization
    if(refactorizationPolicy_ == Amesos2::REPIVOT_ON_REFACTORIZATION) {
      THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: Amesos2LOWSF:Symbolic", Symbolic);
      amesos2Solver->symbolicFactorization();
    }
    {
      THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: Amesos2LOWSF::Factor", Factor);
      amesos2Solver->numericFactorization();
    }

    // Initialize the LOWS object and we are done!
    amesos2Op->initialize(fwdOp,fwdOpSrc,amesos2Solver);
  }
  amesos2Op->setOStream(this->getOStream());
  amesos2Op->setVerbLevel(this->getVerbLevel());
}

template<typename Scalar>
bool Amesos2LinearOpWithSolveFactory<Scalar>::supportsPreconditionerInputType(const EPreconditionerInputType /* precOpType */) const
{
  return false;
}

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const PreconditionerBase<Scalar> > &/* prec */,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->throwOnPrecInput_, std::logic_error,
    "Error, the concrete implementation described as \'"<<this->description()
    <<"\' does not support preconditioners"
    " and has been configured to throw this exception when the"
    " initializePreconditionedOp(...) function is called!"
    );
  this->initializeOp(fwdOpSrc,Op,supportSolveUse); // Ignore the preconditioner!
}

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const LinearOpSourceBase<Scalar> > &/* approxFwdOpSrc */,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->throwOnPrecInput_, std::logic_error,
    "Error, the concrete implementation described as \'"<<this->description()
    <<"\' does not support preconditioners"
    " and has been configured to throw this exception when the"
    " initializePreconditionedOp(...) function is called!"
    );
  this->initializeOp(fwdOpSrc,Op,supportSolveUse); // Ignore the preconditioner!
}

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::uninitializeOp(
  LinearOpWithSolveBase<Scalar> *Op,
  RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
  RCP<const PreconditionerBase<Scalar> > *prec,
  RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
  ESupportSolveUse * /* supportSolveUse */
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
#endif
  Amesos2LinearOpWithSolve<Scalar>
    *amesos2Op = &Teuchos::dyn_cast<Amesos2LinearOpWithSolve<Scalar>>(*Op);
  RCP<const LinearOpSourceBase<Scalar> >
    _fwdOpSrc = amesos2Op->extract_fwdOpSrc(); // Will be null if uninitialized!
  if(fwdOpSrc) *fwdOpSrc = _fwdOpSrc; // It is fine if the client does not want this object back!
  if(prec) *prec = Teuchos::null; // We never keep a preconditioner!
  if(approxFwdOpSrc) *approxFwdOpSrc = Teuchos::null; // never keep approx fwd op!
}

// Overridden from ParameterListAcceptor

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(paramList.get()==NULL);
  // Only validate this level for now here (expect Amesos2 to do its own
  // validation?)
  paramList->validateParameters(*this->getValidParameters(),0);
  paramList_ = paramList;
  solverType_ =
    Amesos2::solverTypeNameToEnumMap.get<Amesos2::ESolverType>(
      paramList_->get(
        SolverType_name
        ,Amesos2::toString(solverType_)
        )
      ,paramList_->name()+"->"+SolverType_name
      );
  refactorizationPolicy_ =
    Amesos2::refactorizationPolicyNameToEnumMap.get<Amesos2::ERefactorizationPolicy>(
      paramList_->get(
        RefactorizationPolicy_name
        ,Amesos2::toString(refactorizationPolicy_)
        )
      ,paramList_->name()+"->"+RefactorizationPolicy_name
      );
  throwOnPrecInput_ = paramList_->get(ThrowOnPreconditionerInput_name,throwOnPrecInput_);
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
}

template<typename Scalar>
RCP<Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::getNonconstParameterList()
{
  return paramList_;
}

template<typename Scalar>
RCP<Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

template<typename Scalar>
RCP<const Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return paramList_;
}

template<typename Scalar>
RCP<const Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  return generateAndGetValidParameters();
}

// Public functions overridden from Teuchos::Describable

template<typename Scalar>
std::string Amesos2LinearOpWithSolveFactory<Scalar>::description() const
{
  std::ostringstream oss;
  oss << "Thyra::Amesos2LinearOpWithSolveFactory{";
  oss << "solverType=" << toString(solverType_);
  oss << "}";
  return oss.str();
}

// private

template<typename Scalar>
RCP<const Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::generateAndGetValidParameters()
{
  static RCP<Teuchos::ParameterList> validParamList;
  if (validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("Amesos2"));
    validParamList->set(SolverType_name, Thyra::Amesos2::solverTypeNames[0]);
    validParamList->set(RefactorizationPolicy_name,
      Amesos2::toString(Amesos2::REPIVOT_ON_REFACTORIZATION));
    validParamList->set(ThrowOnPreconditionerInput_name,bool(true));
    Teuchos::setupVerboseObjectSublist(&*validParamList);
  }
  return validParamList;
}

} // namespace Thyra

#endif // THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
