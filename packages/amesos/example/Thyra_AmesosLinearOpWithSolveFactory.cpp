/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef __sun

#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"

#include "Thyra_AmesosLinearOpWithSolve.hpp"
#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Amesos.h"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#if defined(HAVE_AMESOS_PARDISO) || defined(HAVE_AMESOS_PARDISO_MKL)
#include "Amesos_Pardiso.h"
#endif
#ifdef HAVE_AMESOS_TAUCS
#include "Amesos_Taucs.h"
#endif
#ifdef HAVE_AMESOS_PARAKLETE
#include "Amesos_Paraklete.h"
#endif

namespace {

Teuchos::RCP<Teuchos::Time> overallTimer, constructTimer, symbolicTimer, factorTimer;

const std::string epetraFwdOp_str = "epetraFwdOp";

} // namespace

namespace Thyra {


// Parameter names for Paramter List

const std::string AmesosLinearOpWithSolveFactory::SolverType_name = "Solver Type";

const std::string AmesosLinearOpWithSolveFactory::RefactorizationPolicy_name = "Refactorization Policy";

const std::string AmesosLinearOpWithSolveFactory::ThrowOnPreconditionerInput_name = "Throw on Preconditioner Input";

const std::string AmesosLinearOpWithSolveFactory::Amesos_Settings_name = "Amesos Settings";

// Constructors/initializers/accessors

AmesosLinearOpWithSolveFactory::~AmesosLinearOpWithSolveFactory()
{
#ifdef TEUCHOS_DEBUG
  if(paramList_.get())
    paramList_->validateParameters(
      *this->getValidParameters(),0  // Only validate this level for now!
      );
#endif
}

AmesosLinearOpWithSolveFactory::AmesosLinearOpWithSolveFactory(
  const Amesos::ESolverType                            solverType
  ,const Amesos::ERefactorizationPolicy                refactorizationPolicy
  ,const bool                                          throwOnPrecInput
    )
  :epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
  ,solverType_(solverType)
  ,refactorizationPolicy_(refactorizationPolicy)
  ,throwOnPrecInput_(throwOnPrecInput)
{
  initializeTimers();
}

// Overridden from LinearOpWithSolveFactoryBase

bool AmesosLinearOpWithSolveFactory::isCompatible(
  const LinearOpSourceBase<double> &fwdOpSrc
  ) const
{
  Teuchos::RCP<const LinearOpBase<double> >
    fwdOp = fwdOpSrc.getOp();
  Teuchos::RCP<const Epetra_Operator> epetraFwdOp;
  ETransp                                     epetraFwdOpTransp;
  EApplyEpetraOpAs                            epetraFwdOpApplyAs;
  EAdjointEpetraOp                            epetraFwdOpAdjointSupport;
  double                                      epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOp
    ,&epetraFwdOp,&epetraFwdOpTransp,&epetraFwdOpApplyAs,&epetraFwdOpAdjointSupport,&epetraFwdOpScalar
    );
  if( !dynamic_cast<const Epetra_RowMatrix*>(&*epetraFwdOp) )
    return false;
  return true;
}

Teuchos::RCP<LinearOpWithSolveBase<double> >
AmesosLinearOpWithSolveFactory::createOp() const
{
  return Teuchos::rcp(new AmesosLinearOpWithSolve());
}

void AmesosLinearOpWithSolveFactory::initializeOp(
  const Teuchos::RCP<const LinearOpSourceBase<double> >    &fwdOpSrc
  ,LinearOpWithSolveBase<double>                                   *Op
  ,const ESupportSolveUse                                          supportSolveUse
  ) const
{
  Teuchos::TimeMonitor overallTimeMonitor(*overallTimer);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
#endif
  const Teuchos::RCP<const LinearOpBase<double> > 
    fwdOp = fwdOpSrc->getOp();
  //
  // Unwrap and get the forward Epetra_Operator object
  //
  Teuchos::RCP<const Epetra_Operator> epetraFwdOp;
  ETransp                                     epetraFwdOpTransp;
  EApplyEpetraOpAs                            epetraFwdOpApplyAs;
  EAdjointEpetraOp                            epetraFwdOpAdjointSupport;
  double                                      epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOp,&epetraFwdOp,&epetraFwdOpTransp,&epetraFwdOpApplyAs,&epetraFwdOpAdjointSupport,&epetraFwdOpScalar
    );
  // Get the AmesosLinearOpWithSolve object
  AmesosLinearOpWithSolve
    *amesosOp = &Teuchos::dyn_cast<AmesosLinearOpWithSolve>(*Op);
  //
  // Determine if we must start over or not
  //
  bool startOver = ( amesosOp->get_amesosSolver()==Teuchos::null );
  if(!startOver) {
    startOver =
      (
        epetraFwdOpTransp != amesosOp->get_amesosSolverTransp() ||
        epetraFwdOp.get() != amesosOp->get_epetraLP()->GetOperator()
        // We must start over if the matrix object changes.  This is a
        // weakness of Amesos but there is nothing I can do about this right
        // now!
        );
  }
  //
  // Update the amesos solver
  //
  if(startOver) {
    //
    // This LOWS object has not be initialized yet or is not compatible with the existing
    // 
    // so this is where we setup everything from the ground up.
    //
    // Create the linear problem and set the operator with memory of RCP to Epetra_Operator view!
    Teuchos::RCP<Epetra_LinearProblem>
      epetraLP = Teuchos::rcp(new Epetra_LinearProblem());
    epetraLP->SetOperator(const_cast<Epetra_Operator*>(&*epetraFwdOp));
    Teuchos::set_extra_data< Teuchos::RCP<const Epetra_Operator> >( epetraFwdOp,
      epetraFwdOp_str, Teuchos::inOutArg(epetraLP) );
    // Create the concrete solver
    Teuchos::RCP<Amesos_BaseSolver>
      amesosSolver;
    {
      Teuchos::TimeMonitor constructTimeMonitor(*constructTimer);
      switch(solverType_) {
        case Thyra::Amesos::LAPACK :
          amesosSolver = Teuchos::rcp(new Amesos_Lapack(*epetraLP));
          break;
#ifdef HAVE_AMESOS_KLU
        case Thyra::Amesos::KLU :
          amesosSolver = Teuchos::rcp(new Amesos_Klu(*epetraLP));
          break;
#endif
#ifdef HAVE_AMESOS_MUMPS
        case Thyra::Amesos::MUMPS :
          amesosSolver = Teuchos::rcp(new Amesos_Mumps(*epetraLP));
          break;
#endif
#ifdef HAVE_AMESOS_SCALAPACK
        case Thyra::Amesos::SCALAPACK :
          amesosSolver = Teuchos::rcp(new Amesos_Scalapack(*epetraLP));
          break;
#endif
#ifdef HAVE_AMESOS_UMFPACK
        case Thyra::Amesos::UMFPACK :
          amesosSolver = Teuchos::rcp(new Amesos_Umfpack(*epetraLP));
          break;
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
        case Thyra::Amesos::SUPERLUDIST :
          amesosSolver = Teuchos::rcp(new Amesos_Superludist(*epetraLP));
          break;
#endif
#ifdef HAVE_AMESOS_SUPERLU
        case Thyra::Amesos::SUPERLU :
          amesosSolver = Teuchos::rcp(new Amesos_Superlu(*epetraLP));
          break;
#endif
#ifdef HAVE_AMESOS_DSCPACK
        case Thyra::Amesos::DSCPACK :
          amesosSolver = Teuchos::rcp(new Amesos_Dscpack(*epetraLP));
          break;
#endif
#ifdef HAVE_AMESOS_PARDISO
        case Thyra::Amesos::PARDISO :
          amesosSolver = Teuchos::rcp(new Amesos_Pardiso(*epetraLP));
          break;
#endif
#ifdef HAVE_AMESOS_TAUCS
        case Thyra::Amesos::TAUCS :
          amesosSolver = Teuchos::rcp(new Amesos_Taucs(*epetraLP));
          break;
#endif
#ifdef HAVE_AMESOS_PARAKLETE
        case Thyra::Amesos::PARAKLETE :
          amesosSolver = Teuchos::rcp(new Amesos_Paraklete(*epetraLP));
          break;
#endif
        default:
          TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error
            ,"Error, the solver type ID = " << solverType_ << " is invalid!"
            );
      }
    }
    // Set the parameters
    if(paramList_.get()) amesosSolver->setParameterList(sublist(paramList_,"Amesos Settings"));
    // Do the initial factorization
    {
      Teuchos::TimeMonitor symbolicTimeMonitor(*symbolicTimer);
      amesosSolver->SymbolicFactorization();
    }
    {
      Teuchos::TimeMonitor factorTimeMonitor(*factorTimer);
      amesosSolver->NumericFactorization();
    }
    // Initialize the LOWS object and we are done!
    amesosOp->initialize(fwdOp,fwdOpSrc,epetraLP,amesosSolver,epetraFwdOpTransp,epetraFwdOpScalar);
  }
  else {
    //
    // This LOWS object has already be initialized once so we must just reset
    // the matrix and refactor it.
    //
    // Get non-const pointers to the linear problem and the amesos solver.
    // These const-casts are just fine since the amesosOp in non-const.
    Teuchos::RCP<Epetra_LinearProblem>
      epetraLP = Teuchos::rcp_const_cast<Epetra_LinearProblem>(amesosOp->get_epetraLP());
    Teuchos::RCP<Amesos_BaseSolver>
      amesosSolver = amesosOp->get_amesosSolver();
    // Reset the forward operator with memory of RCP to Epetra_Operator view!
    epetraLP->SetOperator(const_cast<Epetra_Operator*>(&*epetraFwdOp));
    Teuchos::get_extra_data< Teuchos::RCP<const Epetra_Operator> >(epetraLP,epetraFwdOp_str) = epetraFwdOp;
    // Reset the parameters
    if(paramList_.get()) amesosSolver->setParameterList(sublist(paramList_,Amesos_Settings_name));
    // Repivot if asked
    if(refactorizationPolicy_==Amesos::REPIVOT_ON_REFACTORIZATION) {
      Teuchos::TimeMonitor symbolicTimeMonitor(*symbolicTimer);
      amesosSolver->SymbolicFactorization();
    }
    {
      Teuchos::TimeMonitor factorTimeMonitor(*factorTimer);
      amesosSolver->NumericFactorization();
    }
    // Reinitialize the LOWS object and we are done! (we must do this to get the
    // possibly new transpose and scaling factors back in)
    amesosOp->initialize(fwdOp,fwdOpSrc,epetraLP,amesosSolver,epetraFwdOpTransp,epetraFwdOpScalar);
  }
  amesosOp->setOStream(this->getOStream());
  amesosOp->setVerbLevel(this->getVerbLevel());
}

bool AmesosLinearOpWithSolveFactory::supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const
{
  return false;
}

void AmesosLinearOpWithSolveFactory::initializePreconditionedOp(
  const Teuchos::RCP<const LinearOpSourceBase<double> >       &fwdOpSrc
  ,const Teuchos::RCP<const PreconditionerBase<double> >      &prec
  ,LinearOpWithSolveBase<double>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->throwOnPrecInput_, std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' does not support preconditioners "
    "and has been configured to throw this exception when the  initializePreconditionedOp(...) function is called!"
    );
  this->initializeOp(fwdOpSrc,Op,supportSolveUse); // Ignore the preconditioner!
}

void AmesosLinearOpWithSolveFactory::initializePreconditionedOp(
  const Teuchos::RCP<const LinearOpSourceBase<double> >       &fwdOpSrc
  ,const Teuchos::RCP<const LinearOpSourceBase<double> >      &approxFwdOpSrc
  ,LinearOpWithSolveBase<double>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->throwOnPrecInput_, std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' does not support preconditioners "
    "and has been configured to throw this exception when the  initializePreconditionedOp(...) function is called!"
    );
  this->initializeOp(fwdOpSrc,Op,supportSolveUse); // Ignore the preconditioner!
}

void AmesosLinearOpWithSolveFactory::uninitializeOp(
  LinearOpWithSolveBase<double>                               *Op
  ,Teuchos::RCP<const LinearOpSourceBase<double> >    *fwdOpSrc
  ,Teuchos::RCP<const PreconditionerBase<double> >    *prec
  ,Teuchos::RCP<const LinearOpSourceBase<double> >    *approxFwdOpSrc
  ,ESupportSolveUse                                           *supportSolveUse
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
#endif
  AmesosLinearOpWithSolve
    *amesosOp = &Teuchos::dyn_cast<AmesosLinearOpWithSolve>(*Op);
  Teuchos::RCP<const LinearOpSourceBase<double> >
    _fwdOpSrc = amesosOp->extract_fwdOpSrc(); // Will be null if uninitialized!
  if(_fwdOpSrc.get()) {
    // Erase the Epetra_Operator view of the forward operator!
    Teuchos::RCP<Epetra_LinearProblem> epetraLP = amesosOp->get_epetraLP();
    Teuchos::get_extra_data< Teuchos::RCP<const Epetra_Operator> >(
      epetraLP,epetraFwdOp_str
      )
      = Teuchos::null;
    // Note, we did not erase the address of the operator in
    // epetraLP->GetOperator() since it seems that the amesos solvers do not
    // recheck the value of GetProblem()->GetOperator() so you had better not
    // rest this!
  }
  if(fwdOpSrc) *fwdOpSrc = _fwdOpSrc; // It is fine if the client does not want this object back!
  if(prec) *prec = Teuchos::null; // We never keep a preconditioner!
  if(approxFwdOpSrc) *approxFwdOpSrc = Teuchos::null; // We never keep an approximate fwd operator!
}

// Overridden from ParameterListAcceptor

void AmesosLinearOpWithSolveFactory::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(paramList.get()==NULL);
  paramList->validateParameters(*this->getValidParameters(),0); // Only validate this level for now!
  paramList_ = paramList;
  solverType_ =
    Amesos::solverTypeNameToEnumMap.get<Amesos::ESolverType>(
      paramList_->get(
        SolverType_name
        ,Amesos::toString(solverType_)
        )
      ,paramList_->name()+"->"+SolverType_name
      );
  refactorizationPolicy_ = 
    Amesos::refactorizationPolicyNameToEnumMap.get<Amesos::ERefactorizationPolicy>(
      paramList_->get(
        RefactorizationPolicy_name
        ,Amesos::toString(refactorizationPolicy_)
        )
      ,paramList_->name()+"->"+RefactorizationPolicy_name
      );
  throwOnPrecInput_ = paramList_->get(ThrowOnPreconditionerInput_name,throwOnPrecInput_);
}

Teuchos::RCP<Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::getNonconstParameterList()
{
  return paramList_;
}

Teuchos::RCP<Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RCP<const Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::getParameterList() const
{
  return paramList_;
}

Teuchos::RCP<const Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::getValidParameters() const
{
  return generateAndGetValidParameters();
}

// Public functions overridden from Teuchos::Describable

std::string AmesosLinearOpWithSolveFactory::description() const
{
  std::ostringstream oss;
  oss << "Thyra::AmesosLinearOpWithSolveFactory{";
  oss << "solverType=" << toString(solverType_);
  oss << "}";
  return oss.str();
}

// private


void AmesosLinearOpWithSolveFactory::initializeTimers()
{
  if(!overallTimer.get()) {
    overallTimer    = Teuchos::TimeMonitor::getNewTimer("AmesosLOWSF");
    constructTimer  = Teuchos::TimeMonitor::getNewTimer("AmesosLOWSF:InitConstruct");
    symbolicTimer   = Teuchos::TimeMonitor::getNewTimer("AmesosLOWSF:Symbolic");
    factorTimer     = Teuchos::TimeMonitor::getNewTimer("AmesosLOWSF:Factor");
  }
}

Teuchos::RCP<const Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::generateAndGetValidParameters()
{
  static Teuchos::RCP<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("Amesos"));
    validParamList->set(
      SolverType_name
#ifdef HAVE_AMESOS_KLU
      ,Amesos::toString(Amesos::KLU)
#else
      ,Amesos::toString(Amesos::LAPACK)
#endif
      );
    validParamList->set(RefactorizationPolicy_name,Amesos::toString(Amesos::REPIVOT_ON_REFACTORIZATION));
    validParamList->set(ThrowOnPreconditionerInput_name,bool(true));
    validParamList->sublist(Amesos_Settings_name).setParameters(::Amesos::GetValidParameters());
  }
  return validParamList;
}

} // namespace Thyra

#endif // __sun
