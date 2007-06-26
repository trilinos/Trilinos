
#ifndef __sun

#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#include "Thyra_BelosLinearOpWithSolveFactoryDecl.hpp"
#include "Thyra_BelosLinearOpWithSolve.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosThyraAdapter.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace Thyra {

// Parameter names for Paramter List

template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::SolverType_name = "Solver Type";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::SolverType_default = "GMRES";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::MaxIters_name = "Maximum Iterations";
template<class Scalar>
const int         BelosLinearOpWithSolveFactory<Scalar>::MaxIters_default = 400;
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::MaxRestarts_name = "Maximum Restarts";
template<class Scalar>
const int         BelosLinearOpWithSolveFactory<Scalar>::MaxRestarts_default = 25;
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::BlockSize_name = "Block Size";
template<class Scalar>
const int         BelosLinearOpWithSolveFactory<Scalar>::BlockSize_default = 1; // ToDo: We need to make Belos robust when BlockSize > 1 !!!
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::AdaptiveBlockSize_name = "Adaptive Block Size";
template<class Scalar>
const bool        BelosLinearOpWithSolveFactory<Scalar>::AdaptiveBlockSize_default = true;
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::DefaultRelResNorm_name = "Default Rel Res Norm";
template<class Scalar>
const typename BelosLinearOpWithSolveFactory<Scalar>::MagnitudeType
                  BelosLinearOpWithSolveFactory<Scalar>::DefaultRelResNorm_default = 1e-6;
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::OrthoType_name = "Orthogonalization";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::OrthoType_default = "DGKS";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::Restart_Timers_name= "Restart Timers";
template<class Scalar>
const bool        BelosLinearOpWithSolveFactory<Scalar>::Restart_Timers_default = "true";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_name = "GMRES";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_MaxNumberOfKrylovVectors_name = "Max Number of Krylov Vectors";
template<class Scalar>
const int         BelosLinearOpWithSolveFactory<Scalar>::GMRES_MaxNumberOfKrylovVectors_default = 300; // Consistent with NOX::LinearSystemAztecOO
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_Variant_name = "Variant";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_Variant_default = "Standard";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::Outputter_name = "Outputter";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::Outputter_OutputFrequency_name = "Output Frequency";
template<class Scalar>
const int         BelosLinearOpWithSolveFactory<Scalar>::Outputter_OutputFrequency_default = 10;
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::Outputter_OutputMaxResOnly_name = "Show Maximum Residual Norm Only";
template<class Scalar>
const bool        BelosLinearOpWithSolveFactory<Scalar>::Outputter_OutputMaxResOnly_default = true;

// Constructors/initializers/accessors

template<class Scalar>
BelosLinearOpWithSolveFactory<Scalar>::BelosLinearOpWithSolveFactory()
  :useGmres_(true)
{
  updateThisValidParamList();
}

template<class Scalar>
BelosLinearOpWithSolveFactory<Scalar>::BelosLinearOpWithSolveFactory(
  const Teuchos::RCP<PreconditionerFactoryBase<Scalar> >  &precFactory
  )
  :useGmres_(true)
{
  this->setPreconditionerFactory(precFactory);
}

// Overridden from LinearOpWithSolveFactoryBase

template<class Scalar>
bool BelosLinearOpWithSolveFactory<Scalar>::acceptsPreconditionerFactory() const
{
  return true;
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::setPreconditionerFactory(
  const Teuchos::RCP<PreconditionerFactoryBase<Scalar> >  &precFactory
  ,const std::string                                              &precFactoryName
  )
{
  TEST_FOR_EXCEPT(!precFactory.get());
  Teuchos::RCP<const Teuchos::ParameterList>
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
Teuchos::RCP<PreconditionerFactoryBase<Scalar> >
BelosLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory() const
{
  return precFactory_;
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::unsetPreconditionerFactory(
  Teuchos::RCP<PreconditionerFactoryBase<Scalar> >  *precFactory
  ,std::string                                              *precFactoryName
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
Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
BelosLinearOpWithSolveFactory<Scalar>::createOp() const
{
  return Teuchos::rcp(new BelosLinearOpWithSolve<Scalar>());
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeOp(
  const Teuchos::RCP<const LinearOpSourceBase<Scalar> >    &fwdOpSrc
  ,LinearOpWithSolveBase<Scalar>                                   *Op
  ,const ESupportSolveUse                                          supportSolveUse
  ) const
{
  using Teuchos::null;
  initializeOpImpl(fwdOpSrc,null,null,false,Op,supportSolveUse);
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
  const Teuchos::RCP<const LinearOpSourceBase<Scalar> >    &fwdOpSrc
  ,LinearOpWithSolveBase<Scalar>                                   *Op
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
  const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
  ,const Teuchos::RCP<const PreconditionerBase<Scalar> >      &prec
  ,LinearOpWithSolveBase<Scalar>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  using Teuchos::null;
  initializeOpImpl(fwdOpSrc,null,prec,false,Op,supportSolveUse);
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeApproxPreconditionedOp(
  const Teuchos::RCP<const LinearOpSourceBase<Scalar> >      &fwdOpSrc
  ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >     &approxFwdOpSrc
  ,LinearOpWithSolveBase<Scalar>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  using Teuchos::null;
  initializeOpImpl(fwdOpSrc,approxFwdOpSrc,null,false,Op,supportSolveUse);
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::uninitializeOp(
  LinearOpWithSolveBase<Scalar>                               *Op
  ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >    *fwdOpSrc
  ,Teuchos::RCP<const PreconditionerBase<Scalar> >    *prec
  ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >    *approxFwdOpSrc
  ,ESupportSolveUse                                           *supportSolveUse
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(Op==NULL);
#endif
  BelosLinearOpWithSolve<Scalar>
    &belosOp = Teuchos::dyn_cast<BelosLinearOpWithSolve<Scalar> >(*Op);
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > 
    _fwdOpSrc = belosOp.extract_fwdOpSrc();
  Teuchos::RCP<const PreconditionerBase<Scalar> >
    _prec = ( belosOp.isExternalPrec() ? belosOp.extract_prec() : Teuchos::null );
  // Note: above we only extract the preconditioner if it was passed in
  // externally.  Otherwise, we need to hold on to it so that we can reuse it
  // in the next initialization.
  Teuchos::RCP<const LinearOpSourceBase<Scalar> >
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
void BelosLinearOpWithSolveFactory<Scalar>::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  paramList->validateParameters(*this->getValidParameters(),1); // Validate 0th and 1st level deep
  paramList_ = paramList;
  //
  if(precFactory_.get()) {
    // Only reset the PF's PL if the sublist exists or the PF does ot already
    // have a PL.  We don't want to overwrite an externally set PL for the PF
    // if we don't have a nested sublist defined here!
    const bool nestedPFSublistExists = paramList_->isSublist(precFactoryName_);
    const bool alreadyHasSublist = !is_null(precFactory_->getParameterList());
    if( nestedPFSublistExists || !alreadyHasSublist ) {
      precFactory_->setParameterList(Teuchos::sublist(paramList_,precFactoryName_));
    }
  }
  //
  const std::string &solverType = paramList_->get(SolverType_name,SolverType_default);
  if(solverType==GMRES_name)
    useGmres_ = true;
  else if(solverType=="CG")
    useGmres_ = false;
  else
    TEST_FOR_EXCEPTION(
      true,std::logic_error
      ,"Error, \"Solver Type\" = \"" << solverType << "\" is not recognized!"
      "  Valid values are \"GMRES\" and \"CG\""
      );
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
  // ToDo: Replace above with Teuchos::StringToIntMap ...
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::getParameterList()
{
  return paramList_;
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return paramList_;
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
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
Teuchos::RCP<const Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::generateAndGetValidParameters()
{
  static Teuchos::RCP<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("BelosLinearOpWithSolveFactory"));
    validParamList->set(SolverType_name,SolverType_default);
    validParamList->set(MaxIters_name,MaxIters_default);
    validParamList->set(MaxRestarts_name,MaxRestarts_default);
    validParamList->set(BlockSize_name,BlockSize_default);
    validParamList->set(AdaptiveBlockSize_name,AdaptiveBlockSize_default);
    validParamList->set(DefaultRelResNorm_name,DefaultRelResNorm_default);
    validParamList->set(OrthoType_name,OrthoType_default);
    Teuchos::ParameterList
      &gmresSL = validParamList->sublist(GMRES_name);
    gmresSL.set(GMRES_MaxNumberOfKrylovVectors_name,GMRES_MaxNumberOfKrylovVectors_default);
    gmresSL.set(GMRES_Variant_name,GMRES_Variant_default);
    gmresSL.set(Restart_Timers_name,Restart_Timers_default);
    Teuchos::ParameterList
      &outputterSL = validParamList->sublist(Outputter_name);
    outputterSL.set(Outputter_OutputFrequency_name,Outputter_OutputFrequency_default);
    outputterSL.set(Outputter_OutputMaxResOnly_name,Outputter_OutputMaxResOnly_default);
  }
  return validParamList;
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::updateThisValidParamList()
{
  thisValidParamList_ = Teuchos::rcp(
    new Teuchos::ParameterList(*generateAndGetValidParameters())
    );
  if(precFactory_.get()) {
    Teuchos::RCP<const Teuchos::ParameterList>
      precFactoryValidParamList = precFactory_->getValidParameters();
    if(precFactoryValidParamList.get()) {
      thisValidParamList_->sublist(precFactoryName_).setParameters(*precFactoryValidParamList);
    }
  }
  Teuchos::setupVerboseObjectSublist(&*thisValidParamList_);
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeOpImpl(
  const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
  ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >      &approxFwdOpSrc
  ,const Teuchos::RCP<const PreconditionerBase<Scalar> >      &prec_in
  ,const bool                                                         reusePrec
  ,LinearOpWithSolveBase<Scalar>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::set_extra_data;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef MultiVectorBase<Scalar>    MV_t;
  typedef LinearOpBase<Scalar>       LO_t;

  const Teuchos::RCP<Teuchos::FancyOStream> out       = this->getOStream();
  const Teuchos::EVerbosityLevel                    verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::BelosLinearOpWithSolveFactory<"<<ST::name()<<">::initializeOpImpl(...) ...\n";

  typedef Teuchos::VerboseObjectTempState<PreconditionerFactoryBase<Scalar> > VOTSPF;
  VOTSPF precFactoryOutputTempState(precFactory_,out,verbLevel);
  
  TEST_FOR_EXCEPT(Op==NULL);
  TEST_FOR_EXCEPT(fwdOpSrc.get()==NULL);
  TEST_FOR_EXCEPT(fwdOpSrc->getOp().get()==NULL);
  Teuchos::RCP<const LinearOpBase<Scalar> >
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
  RCP<PreconditionerBase<Scalar> >         myPrec = Teuchos::null;
  RCP<const PreconditionerBase<Scalar> >   prec = Teuchos::null;
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
  int oldMaxNumberOfKrylovVectors = 0;
  bool oldIsExternalPrec = false;
  RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> >     oldLP = Teuchos::null;
  RCP<Belos::SolverManager<Scalar,MV_t,LO_t> >     oldIterSolver = Teuchos::null;
  RCP<const LinearOpSourceBase<Scalar> >           oldFwdOpSrc = Teuchos::null;
  RCP<const LinearOpSourceBase<Scalar> >           oldApproxFwdOpSrc = Teuchos::null;   
  ESupportSolveUse                                 oldSupportSolveUse = SUPPORT_SOLVE_UNSPECIFIED;

  belosOp->uninitialize( &oldLP,
                         &oldMaxNumberOfKrylovVectors,
                         NULL,
                         &oldIterSolver,
                         &oldFwdOpSrc,
                         NULL,
                         &oldIsExternalPrec,
                         &oldApproxFwdOpSrc,
                         &oldSupportSolveUse );
  //
  // Create the Belos linear problem
  // NOTE:  If one exists already, reuse it.
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
    RCP<const LinearOpBase<Scalar> > left        = prec->getLeftPrecOp();
    RCP<const LinearOpBase<Scalar> > right       = prec->getRightPrecOp();
    TEST_FOR_EXCEPTION(
      !( left.get() || right.get() || unspecified.get() ), std::logic_error
      ,"Error, at least one preconditoner linear operator objects must be set!"
      );
    if(unspecified.get()) {
      lp->setRightPrec(unspecified);
      // ToDo: Allow user to determine whether this should be placed on the
      // left or on the right through a parameter in the parameter list!
    }
    else {
      // Set a left, right or split preconditioner
      TEST_FOR_EXCEPTION(
        left.get(),std::logic_error
        ,"Error, we can not currently handle a left preconditioner!"
        );
      lp->setRightPrec(right);
    }
  }
  if(myPrec.get()) {
    set_extra_data<RCP<PreconditionerBase<Scalar> > >(myPrec,"Belos::InternalPrec",
							      &lp, Teuchos::POST_DESTROY, false);
  }
  else if(prec.get()) {
    set_extra_data<RCP<const PreconditionerBase<Scalar> > >(prec,"Belos::ExternalPrec",
								    &lp, Teuchos::POST_DESTROY, false);
  }
  //
  // Generate the parameter list
  //
  typedef Belos::SolverManager<Scalar,MV_t,LO_t> IterativeSolver_t;
  RefCountPtr<IterativeSolver_t> iterativeSolver = Teuchos::null;
  RefCountPtr<Teuchos::ParameterList> solverPL = Teuchos::rcp( new Teuchos::ParameterList() );

  // Set the block size
  int blockSize = BlockSize_default;
  if(paramList_.get()) {
    blockSize = paramList_->get(BlockSize_name,blockSize);
  }
  solverPL->set(BlockSize_name, blockSize);

  const bool adaptiveBlockSize =
    ( paramList_.get()
      ? paramList_->get(AdaptiveBlockSize_name,AdaptiveBlockSize_default)
      : AdaptiveBlockSize_default );
  solverPL->set(AdaptiveBlockSize_name, adaptiveBlockSize);

  // Set the verbosity
  const int belosVerbLevel =
    (
      verbLevel == Teuchos::VERB_DEFAULT || static_cast<int>(verbLevel)>=static_cast<int>(Teuchos::VERB_LOW)
      ? Belos::Warnings | Belos::FinalSummary | Belos::IterationDetails
      : Belos::Errors 
      );
  solverPL->set("Verbosity", belosVerbLevel);

  // Set the status test parameters
  int         defaultMaxIterations = MaxIters_default;
  int         defaultMaxRestarts   = MaxRestarts_default;
  ScalarMag   defaultResNorm       = DefaultRelResNorm_default;
  int         outputFrequency      = Outputter_OutputFrequency_default;
  bool        outputMaxResOnly     = Outputter_OutputMaxResOnly_default;
  if(paramList_.get()) {
    defaultMaxIterations = paramList_->get(MaxIters_name,defaultMaxIterations);
    defaultMaxRestarts = paramList_->get(MaxRestarts_name,defaultMaxRestarts);
    defaultResNorm = paramList_->get(DefaultRelResNorm_name,defaultResNorm);
    Teuchos::ParameterList &outputterSL = paramList_->sublist(Outputter_name);
    outputFrequency = outputterSL.get(Outputter_OutputFrequency_name,outputFrequency);
    outputMaxResOnly = outputterSL.get(Outputter_OutputMaxResOnly_name,outputMaxResOnly);
  }
  solverPL->set(MaxIters_name,defaultMaxIterations);
  solverPL->set(MaxRestarts_name,defaultMaxRestarts);
  solverPL->set("Convergence Tolerance", defaultResNorm);
  solverPL->set(Outputter_OutputFrequency_name, outputFrequency);
  solverPL->set(Outputter_OutputMaxResOnly_name, outputMaxResOnly);
  
  // Set orthogonalization parameter
  std::string orthoType = OrthoType_default;
  if(paramList_.get()) {
    orthoType = paramList_->get(OrthoType_name,OrthoType_default);
  }
  solverPL->set(OrthoType_name, orthoType);

  int maxNumberOfKrylovVectors = -1; // Only gets used if getPL.get()!=NULL
  bool restartTimers = Restart_Timers_default;
  if(useGmres_) {
    // Set the PL
    solverPL->set("Num Blocks",1);
    // Note, the "Length" will be reset based on the number of RHS in the
    // BelosLOWS::solve(...) function!  This is needed to avoid memory
    // problems!  Above I just set it to 1 to avoid any memory allocation
    // problems!
    maxNumberOfKrylovVectors = GMRES_MaxNumberOfKrylovVectors_default;
    std::string GMRES_Variant = GMRES_Variant_default;
    if(paramList_.get()) {
      Teuchos::ParameterList &_gmresPL = paramList_->sublist(GMRES_name);
      maxNumberOfKrylovVectors = _gmresPL.get(GMRES_MaxNumberOfKrylovVectors_name,GMRES_MaxNumberOfKrylovVectors_default);
      GMRES_Variant = _gmresPL.get(GMRES_Variant_name,GMRES_Variant_default);
      restartTimers = _gmresPL.get(Restart_Timers_name,Restart_Timers_default);
    }
    solverPL->set(Restart_Timers_name, restartTimers);
    if (GMRES_Variant == "Flexible") {
      solverPL->set("Flexible Gmres", true);
    }
    // 
    // Create the solver
    // 
    if (oldIterSolver != Teuchos::null) {
      iterativeSolver = oldIterSolver;
      iterativeSolver->setProblem( lp );
      iterativeSolver->setParameters( solverPL );
    }
    else {    
      if (GMRES_Variant == "Pseudo") {
	iterativeSolver = rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
      }
      else {
	iterativeSolver = rcp(new Belos::BlockGmresSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
      }
    }
  }
  else {
    iterativeSolver = rcp(new Belos::BlockCGSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
  }
  
  //
  // Initialize the LOWS object
  //
  belosOp->initialize(
		      lp,maxNumberOfKrylovVectors,solverPL,iterativeSolver
		      ,fwdOpSrc,prec,myPrec.get()==NULL,approxFwdOpSrc,supportSolveUse
		      );
  belosOp->setOStream(out);
  belosOp->setVerbLevel(verbLevel);
#ifdef TEUCHOS_DEBUG
  if(paramList_.get()) {
    // Make sure we read the list correctly
    paramList_->validateParameters(*this->getValidParameters(),1); // Validate 0th and 1st level deep
  }
#endif
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nLeaving Thyra::BelosLinearOpWithSolveFactory<"<<ST::name()<<">::initializeOpImpl(...) ...\n";
  
}

} // namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#endif // __sun
