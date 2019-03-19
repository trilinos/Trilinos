/*
// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP


#include "Thyra_BelosLinearOpWithSolveFactory_decl.hpp"
#include "Thyra_BelosLinearOpWithSolve.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"

#include "Thyra_BelosSolverFactory.hpp"
#include "BelosThyraAdapter.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StandardValidatorXMLConverters.hpp"


namespace Thyra {


// Parameter names for Parameter List

template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::SolverType_name = "Solver Type";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::SolverType_default = "PSEUDO BLOCK GMRES";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::SolverTypes_name = "Solver Types";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::ConvergenceTestFrequency_name = "Convergence Test Frequency";

namespace {
const std::string LeftPreconditionerIfUnspecified_name = "Left Preconditioner If Unspecified";
}

// Constructors/initializers/accessors


template<class Scalar>
BelosLinearOpWithSolveFactory<Scalar>::BelosLinearOpWithSolveFactory()
  : solverName_(SolverType_default),
    convergenceTestFrequency_(1)
{
  updateThisValidParamList();
}


template<class Scalar>
BelosLinearOpWithSolveFactory<Scalar>::BelosLinearOpWithSolveFactory(
  const RCP<PreconditionerFactoryBase<Scalar> > &precFactory
  )
  : solverName_(SolverType_default)
{
  this->setPreconditionerFactory(precFactory, "");
}


// Overridden from LinearOpWithSolveFactoryBase


template<class Scalar>
bool BelosLinearOpWithSolveFactory<Scalar>::acceptsPreconditionerFactory() const
{
  return true;
}


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::setPreconditionerFactory(
  const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
  const std::string &precFactoryName
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(!precFactory.get());
  RCP<const Teuchos::ParameterList>
    precFactoryValidPL = precFactory->getValidParameters();
  const std::string _precFactoryName =
    ( precFactoryName != ""
      ? precFactoryName
      : ( precFactoryValidPL.get() ? precFactoryValidPL->name() : "GENERIC PRECONDITIONER FACTORY" )
      );
  precFactory_ = precFactory;
  precFactoryName_ = _precFactoryName;
  updateThisValidParamList();
}


template<class Scalar>
RCP<PreconditionerFactoryBase<Scalar> >
BelosLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory() const
{
  return precFactory_;
}


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::unsetPreconditionerFactory(
  RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
  std::string *precFactoryName
  )
{
  if(precFactory) *precFactory = precFactory_;
  if(precFactoryName) *precFactoryName = precFactoryName_;
  precFactory_ = Teuchos::null;
  precFactoryName_ = "";
  updateThisValidParamList();
}


template<class Scalar>
bool BelosLinearOpWithSolveFactory<Scalar>::isCompatible(
  const LinearOpSourceBase<Scalar> &fwdOpSrc
  ) const
{
  if(precFactory_.get())
    return precFactory_->isCompatible(fwdOpSrc);
  return true; // Without a preconditioner, we are compatible with all linear operators!
}


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
BelosLinearOpWithSolveFactory<Scalar>::createOp() const
{
  return Teuchos::rcp(new BelosLinearOpWithSolve<Scalar>());
}


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  using Teuchos::null;
  initializeOpImpl(fwdOpSrc,null,null,false,Op,supportSolveUse);
}


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op
  ) const
{
  using Teuchos::null;
  initializeOpImpl(fwdOpSrc,null,null,true,Op,SUPPORT_SOLVE_UNSPECIFIED);
}


template<class Scalar>
bool BelosLinearOpWithSolveFactory<Scalar>::supportsPreconditionerInputType(
  const EPreconditionerInputType precOpType
  ) const
{
  if(precFactory_.get())
    return true;
  return (precOpType==PRECONDITIONER_INPUT_TYPE_AS_OPERATOR);
}


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const PreconditionerBase<Scalar> > &prec,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  using Teuchos::null;
  initializeOpImpl(fwdOpSrc,null,prec,false,Op,supportSolveUse);
}


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeApproxPreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  using Teuchos::null;
  initializeOpImpl(fwdOpSrc,approxFwdOpSrc,null,false,Op,supportSolveUse);
}


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::uninitializeOp(
  LinearOpWithSolveBase<Scalar> *Op,
  RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
  RCP<const PreconditionerBase<Scalar> > *prec,
  RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
  ESupportSolveUse *supportSolveUse
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
#endif
  BelosLinearOpWithSolve<Scalar>
    &belosOp = Teuchos::dyn_cast<BelosLinearOpWithSolve<Scalar> >(*Op);
  RCP<const LinearOpSourceBase<Scalar> > 
    _fwdOpSrc = belosOp.extract_fwdOpSrc();
  RCP<const PreconditionerBase<Scalar> >
    _prec = ( belosOp.isExternalPrec() ? belosOp.extract_prec() : Teuchos::null );
  // Note: above we only extract the preconditioner if it was passed in
  // externally.  Otherwise, we need to hold on to it so that we can reuse it
  // in the next initialization.
  RCP<const LinearOpSourceBase<Scalar> >
    _approxFwdOpSrc = belosOp.extract_approxFwdOpSrc();
  ESupportSolveUse
    _supportSolveUse = belosOp.supportSolveUse();
  if(fwdOpSrc) *fwdOpSrc = _fwdOpSrc;
  if(prec) *prec = _prec;
  if(approxFwdOpSrc) *approxFwdOpSrc = _approxFwdOpSrc;
  if(supportSolveUse) *supportSolveUse = _supportSolveUse;
}


// Overridden from ParameterListAcceptor


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(paramList.get()==NULL);
  std::string solverName, solverNameUC;
  solverName = paramList->get<std::string>(SolverType_name, SolverType_default);
  solverNameUC = Belos::Impl::upperCase (solverName);
  paramList->set(SolverType_name, solverNameUC);

  Teuchos::ParameterList SL;
  for (auto it = paramList->sublist(SolverTypes_name).begin(); it != paramList->sublist(SolverTypes_name).end(); it++) {
    solverName = it->first;
    solverNameUC = Belos::Impl::upperCase (solverName);
    SL.sublist(solverNameUC).setParameters(paramList->sublist(SolverTypes_name).sublist(solverName));
  }
  paramList->remove(SolverTypes_name);
  paramList->sublist(SolverTypes_name) = SL;

  paramList->validateParametersAndSetDefaults(*this->getValidParameters(), 0);
  paramList_ = paramList;
  solverName_ = paramList_->get<std::string>(SolverType_name);
  convergenceTestFrequency_ =
    Teuchos::getParameter<int>(*paramList_, ConvergenceTestFrequency_name);
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
}


template<class Scalar>
RCP<Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return paramList_;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  return thisValidParamList_;
}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string BelosLinearOpWithSolveFactory<Scalar>::description() const
{
  std::ostringstream oss;
  oss << "Thyra::BelosLinearOpWithSolveFactory";
  //oss << "{";
  // ToDo: Fill this in some!
  //oss << "}";
  return oss.str();
}


// private


template<class Scalar>
RCP<const Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::generateAndGetValidParameters()
{
  using Teuchos::as;
  using Teuchos::tuple;

  typedef MultiVectorBase<Scalar> MV_t;
  typedef LinearOpBase<Scalar> LO_t;
  static RCP<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("BelosLinearOpWithSolveFactory"));
    validParamList->set(ConvergenceTestFrequency_name, as<int>(1),
      "Number of linear solver iterations to skip between applying"
      " user-defined convergence test.");
    validParamList->set(
      LeftPreconditionerIfUnspecified_name, false,
      "If the preconditioner does not specify if it is left or right, and this\n"
      "option is set to true, put the preconditioner on the left side.\n"
      "Historically, preconditioning is on the right. Some solvers may not\n"
      "support left preconditioning.");
    Teuchos::ParameterList
      &solverTypesSL = validParamList->sublist(SolverTypes_name);

    validParamList->set(SolverType_name, SolverType_default);

    Belos::ThyraSolverFactory<Scalar> factory;
    Teuchos::Array<std::string> supportedSolvers = factory.supportedSolverNames();
    typedef Belos::SolverManager<Scalar,MV_t,LO_t> IterativeSolver_t;
    for (size_t i = 0; i < supportedSolvers.size(); i++) {
      try {
        RCP<IterativeSolver_t> iterativeSolver = factory.create ( supportedSolvers[i], Teuchos::null );
        solverTypesSL.sublist(supportedSolvers[i]) = *(iterativeSolver->getValidParameters());
      } catch (std::invalid_argument) {
        // pass
      }
    }
    validParamList->sublist(SolverTypes_name).setParameters(solverTypesSL);
  }
  return validParamList;
}


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::updateThisValidParamList()
{
  thisValidParamList_ = Teuchos::rcp(
    new Teuchos::ParameterList(*generateAndGetValidParameters())
    );
  Teuchos::setupVerboseObjectSublist(&*thisValidParamList_);
}


template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeOpImpl(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
  const RCP<const PreconditionerBase<Scalar> > &prec_in,
  const bool reusePrec,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{

  using Teuchos::rcp;
  using Teuchos::set_extra_data;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef MultiVectorBase<Scalar> MV_t;
  typedef LinearOpBase<Scalar> LO_t;

  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::BelosLinearOpWithSolveFactory<"<<ST::name()<<">::initializeOpImpl(...) ...\n";

  // These lines are changing the verbosity of the preconditioner, which has its own verbose object list,
  // so I am commenting these out, as it is not the job of the linear solver to dictate preconditioner verbosity.
  //typedef Teuchos::VerboseObjectTempState<PreconditionerFactoryBase<Scalar> > VOTSPF;
  //VOTSPF precFactoryOutputTempState(precFactory_,out,verbLevel);
  
  TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(fwdOpSrc.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(fwdOpSrc->getOp().get()==NULL);
  RCP<const LinearOpBase<Scalar> >
    fwdOp = fwdOpSrc->getOp(),
    approxFwdOp = ( approxFwdOpSrc.get() ? approxFwdOpSrc->getOp() : Teuchos::null );

  //
  // Get the BelosLinearOpWithSolve interface
  //

  BelosLinearOpWithSolve<Scalar>
    *belosOp = &Teuchos::dyn_cast<BelosLinearOpWithSolve<Scalar> >(*Op);

  //
  // Get/Create the preconditioner
  //

  RCP<PreconditionerBase<Scalar> > myPrec = Teuchos::null;
  RCP<const PreconditionerBase<Scalar> > prec = Teuchos::null;
  if(prec_in.get()) {
    // Use an externally defined preconditioner
    prec = prec_in;
  }
  else {
    // Try and generate a preconditioner on our own
    if(precFactory_.get()) {
      myPrec =
        ( !belosOp->isExternalPrec()
          ? Teuchos::rcp_const_cast<PreconditionerBase<Scalar> >(belosOp->extract_prec())
          : Teuchos::null
          );
      bool hasExistingPrec = false;
      if(myPrec.get()) {
        hasExistingPrec = true;
        // ToDo: Get the forward operator and validate that it is the same
        // operator that is used here!
      }
      else {
        hasExistingPrec = false;
        myPrec = precFactory_->createPrec();
      }
      if( hasExistingPrec && reusePrec ) {
        // Just reuse the existing preconditioner again!
      }
      else {
        // Update the preconditioner
        if(approxFwdOp.get())
          precFactory_->initializePrec(approxFwdOpSrc,&*myPrec);
        else
          precFactory_->initializePrec(fwdOpSrc,&*myPrec);
      }
      prec = myPrec;
    }
  }

  //
  // Uninitialize the current solver object
  //

  bool oldIsExternalPrec = false;
  RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> > oldLP = Teuchos::null;
  RCP<Belos::SolverManager<Scalar,MV_t,LO_t> > oldIterSolver = Teuchos::null;
  RCP<const LinearOpSourceBase<Scalar> > oldFwdOpSrc = Teuchos::null;
  RCP<const LinearOpSourceBase<Scalar> > oldApproxFwdOpSrc = Teuchos::null;   
  ESupportSolveUse oldSupportSolveUse = SUPPORT_SOLVE_UNSPECIFIED;

  belosOp->uninitialize( &oldLP, NULL, &oldIterSolver, &oldFwdOpSrc,
    NULL, &oldIsExternalPrec, &oldApproxFwdOpSrc, &oldSupportSolveUse );

  //
  // Create the Belos linear problem
  // NOTE:  If one exists already, reuse it.
  //

  typedef Belos::LinearProblem<Scalar,MV_t,LO_t> LP_t;
  RCP<LP_t> lp;
  if (oldLP != Teuchos::null) {
    lp = oldLP;
  }
  else {
    lp = rcp(new LP_t());
  }

  //
  // Set the operator
  //

  lp->setOperator(fwdOp);

  //
  // Set the preconditioner
  //

  if(prec.get()) {
    RCP<const LinearOpBase<Scalar> > unspecified = prec->getUnspecifiedPrecOp();
    RCP<const LinearOpBase<Scalar> > left = prec->getLeftPrecOp();
    RCP<const LinearOpBase<Scalar> > right = prec->getRightPrecOp();
    TEUCHOS_TEST_FOR_EXCEPTION(
      !( left.get() || right.get() || unspecified.get() ), std::logic_error
      ,"Error, at least one preconditoner linear operator objects must be set!"
      );
    if (nonnull(unspecified)) {
      if (paramList_->get<bool>(LeftPreconditionerIfUnspecified_name, false))
        lp->setLeftPrec(unspecified);
      else
        lp->setRightPrec(unspecified);
    }
    else if (nonnull(left)) {
      lp->setLeftPrec(left);
    }
    else if (nonnull(right)) {
      lp->setRightPrec(right);
    }
    else {
      // Set a left, right or split preconditioner
      TEUCHOS_TEST_FOR_EXCEPTION(
        nonnull(left) && nonnull(right),std::logic_error
        ,"Error, we can not currently handle split preconditioners!"
        );
    }
  }
  if(myPrec.get()) {
    set_extra_data<RCP<PreconditionerBase<Scalar> > >(myPrec,"Belos::InternalPrec",
      Teuchos::inOutArg(lp), Teuchos::POST_DESTROY, false);
  }
  else if(prec.get()) {
    set_extra_data<RCP<const PreconditionerBase<Scalar> > >(prec,"Belos::ExternalPrec",
      Teuchos::inOutArg(lp), Teuchos::POST_DESTROY, false);
  }

  //
  // Generate the parameter list.
  //

  Belos::ThyraSolverFactory<Scalar> factory;
  typedef Belos::SolverManager<Scalar,MV_t,LO_t> IterativeSolver_t;
  RCP<IterativeSolver_t> iterativeSolver = Teuchos::null;
  Teuchos::ParameterList &solverTypesPL = paramList_->sublist(SolverTypes_name);
  RCP<Teuchos::ParameterList> solverPL = Teuchos::rcpFromRef(solverTypesPL.sublist(solverName_));

  // Create the solver
  if (oldIterSolver != Teuchos::null) {
    iterativeSolver = oldIterSolver;
    iterativeSolver->setProblem( lp );
    iterativeSolver->setParameters( solverPL );
  }
  else {
    iterativeSolver = factory.create ( solverName_, solverPL );
    iterativeSolver->setProblem ( lp );
  }

  //
  // Initialize the LOWS object
  //

  belosOp->initialize(
    lp, solverPL, iterativeSolver,
    fwdOpSrc, prec, myPrec.get()==NULL, approxFwdOpSrc,
    supportSolveUse, convergenceTestFrequency_
    );
  belosOp->setOStream(out);
  belosOp->setVerbLevel(verbLevel);
#ifdef TEUCHOS_DEBUG
  if(paramList_.get()) {
    // Make sure we read the list correctly
    paramList_->validateParameters(*this->getValidParameters(),1); // Validate 0th and 1st level deep
  }
#endif
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nLeaving Thyra::BelosLinearOpWithSolveFactory<"<<ST::name()<<">::initializeOpImpl(...) ...\n";
  
}


} // namespace Thyra


#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
