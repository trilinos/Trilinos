
#ifndef __sun

#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#include "Thyra_BelosLinearOpWithSolveFactoryDecl.hpp"
#include "Thyra_BelosLinearOpWithSolve.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "BelosBlockGmres.hpp"
#include "BelosBlockCG.hpp"
#include "BelosThyraAdapter.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestMaxRestarts.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestOutputter.hpp"
#include "BelosStatusTestCombo.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace Thyra {

// Parameter names for Paramter List

template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::SolverType_name = "Solver Type";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::MaxIters_name = "Max Iters";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::MaxRestarts_name = "Max Restarts";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::BlockSize_name = "Block Size";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::AdjustableBlockSize_name = "Adjustable Block Size";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::DefaultRelResNorm_name = "Default Rel Res Norm";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_name = "GMRES";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_Length_name = "Length";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_Variant_name = "Variant";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::Outputter_name = "Outputter";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::Outputter_OutputFrequency_name = "Output Frequency";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::Outputter_OutputMaxResOnly_name = "Output Max Res Only";
//template<class Scalar>
//const std::string BelosLinearOpWithSolveFactory<Scalar>::Preconditioner_name = "Preconditioner";

// Constructors/initializers/accessors

template<class Scalar>
BelosLinearOpWithSolveFactory<Scalar>::BelosLinearOpWithSolveFactory()
{
  updateThisValidParamList();
}

template<class Scalar>
BelosLinearOpWithSolveFactory<Scalar>::BelosLinearOpWithSolveFactory(
  const Teuchos::RefCountPtr<PreconditionerFactoryBase<Scalar> >  &precFactory
  )
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
  const Teuchos::RefCountPtr<PreconditionerFactoryBase<Scalar> >  &precFactory
  ,const std::string                                              &precFactoryName
  )
{
  TEST_FOR_EXCEPT(!precFactory.get());
  Teuchos::RefCountPtr<const Teuchos::ParameterList>
    precFactoryValidPL = precFactory->getValidParameters();
  const std::string _precFactoryName =
    ( precFactoryName.length()
      ? precFactoryName
      : ( precFactoryValidPL.get() ? precFactoryValidPL->name() : "GENERIC PRECONDITIONER FACTORY" )
      );
  precFactory_ = precFactory;
  precFactoryName_ = precFactoryName; // Remember what the user actually passed in!
  updateThisValidParamList();
}

template<class Scalar>
Teuchos::RefCountPtr<PreconditionerFactoryBase<Scalar> >
BelosLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory() const
{
  return precFactory_;
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::unsetPreconditionerFactory(
  Teuchos::RefCountPtr<PreconditionerFactoryBase<Scalar> >  *precFactory
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
  const LinearOpBase<Scalar> &fwdOp
  ) const
{
  return true; // Until we build a preconditioner internally
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
BelosLinearOpWithSolveFactory<Scalar>::createOp() const
{
  return Teuchos::rcp(new BelosLinearOpWithSolve<Scalar>());
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeOp(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &fwdOp
  ,LinearOpWithSolveBase<Scalar>                             *Op
  ,const ESupportSolveUse                                    supportSolveUse
  ) const
{

  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::set_extra_data;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef MultiVectorBase<Scalar>    MV_t;
  typedef LinearOpBase<Scalar>       LO_t;

  //
  // Unwrap the forward operator
  //

  // ToDo: Once we generate preconditioners internally

  //
  // Get the BelosLinearOpWithSolve interface
  //
  TEST_FOR_EXCEPT(Op==NULL);
  BelosLinearOpWithSolve<Scalar>
    *belosOp = &Teuchos::dyn_cast<BelosLinearOpWithSolve<Scalar> >(*Op);
  //
  // Get/Create the preconditioner
  //
  RefCountPtr<PreconditionerBase<Scalar> >         myPrec = Teuchos::null;
  RefCountPtr<const PreconditionerBase<Scalar> >   prec = Teuchos::null;
  if(precFactory_.get()) {
    if(paramList_.get()) {
      precFactory_->setParameterList(Teuchos::sublist(paramList_,precFactoryName_));
    }
    myPrec = precFactory_->createPrec();
    precFactory_->initializePrec(fwdOp,&*myPrec);
    prec = myPrec;
  }
  //
  // Create the Belos linear problem
  //
  typedef Belos::LinearProblem<Scalar,MV_t,LO_t> LP_t;
  RefCountPtr<LP_t>
    lp = rcp(new LP_t());
  //
  // Set the operator
  //
  lp->SetOperator(fwdOp);
  //
  // Set the preconditioner
  //
  if(prec.get()) {
    RefCountPtr<const LinearOpBase<Scalar> > unspecified = prec->getUnspecifiedPrecOp();
    RefCountPtr<const LinearOpBase<Scalar> > left        = prec->getLeftPrecOp();
    RefCountPtr<const LinearOpBase<Scalar> > right       = prec->getRightPrecOp();
    TEST_FOR_EXCEPTION(
      !( left.get() || right.get() || unspecified.get() ), std::logic_error
      ,"Error, at least one preconditoner linear operator objects must be set!"
      );
    if(unspecified.get()) {
      lp->SetRightPrec(unspecified);
      // ToDo: Allow user to determine whether this should be placed on the
      // left or on the right through a parameter in the parameter list!
    }
    else {
      // Set a left, right or split preconditioner
      TEST_FOR_EXCEPTION(
        left.get(),std::logic_error
        ,"Error, we can not currently handle a left preconditioner!"
        );
      lp->SetRightPrec(right);
    }
  }
  if(myPrec.get()) {
    set_extra_data<RefCountPtr<PreconditionerBase<Scalar> > >(myPrec,"Belos::InternalPrec",&lp);
  }
  else if(prec.get()) {
    set_extra_data<RefCountPtr<const PreconditionerBase<Scalar> > >(prec,"Belos::ExternalPrec",&lp);
  }
  //
  // Set the block size
  //
  int blockSize = 1;
  if(paramList_.get()) {
    blockSize = paramList_->get(BlockSize_name,blockSize);
  }
  lp->SetBlockSize(blockSize);
  //
  // Create the output manager 
  //
  typedef Belos::OutputManager<Scalar> OutputManager_t;
  //int belosVerbLevel = Belos::Warnings | Belos::FinalSummary;
  int belosVerbLevel = Belos::Warnings | Belos::FinalSummary | Belos::IterationDetails;
  // ToDo: Set these from the verbosity level
  RefCountPtr<Teuchos::FancyOStream>
    ostream = this->getOStream();
  RefCountPtr<OutputManager_t>
    outputManager = rcp(new OutputManager_t(0,belosVerbLevel));
  // Note: The stream itself will be set in the BelosLinearOpWithSolve object!
  //
  // Create the default status test
  //
  int         defaultMaxIterations = 400;
  int         defaultMaxRestarts = 25;
  ScalarMag   defaultResNorm = 1e-6;
  int         outputFrequency = 1; // ToDo: Set from the parameter list!
  bool        outputMaxResOnly = false; // ToDo: Set from the parameter list!
  if(paramList_.get()) {
    defaultMaxIterations = paramList_->get(MaxIters_name,defaultMaxIterations);
    defaultMaxRestarts = paramList_->get(MaxRestarts_name,defaultMaxRestarts);
    defaultResNorm = paramList_->get(DefaultRelResNorm_name,defaultResNorm);
    Teuchos::ParameterList &outputterSL = paramList_->sublist(Outputter_name);
    outputFrequency = outputterSL.get(Outputter_OutputFrequency_name,outputFrequency);
    outputMaxResOnly = outputterSL.get(Outputter_OutputMaxResOnly_name,outputMaxResOnly);
  }
  //
  typedef Belos::StatusTestResNorm<Scalar,MV_t,LO_t>      StatusTestResNorm_t;
  typedef Belos::StatusTestMaxIters<Scalar,MV_t,LO_t>     StatusTestMaxIters_t;
  typedef Belos::StatusTestMaxRestarts<Scalar,MV_t,LO_t>  StatusTestMaxRestarts_t;
  typedef Belos::StatusTestOutputter<Scalar,MV_t,LO_t>    StatusTestOutputter_t;
  typedef Belos::StatusTestCombo<Scalar,MV_t,LO_t>        StatusTestCombo_t;
  RefCountPtr<StatusTestMaxIters_t>
    maxItersST = rcp(new StatusTestMaxIters_t(defaultMaxIterations));
  RefCountPtr<StatusTestMaxRestarts_t>
    maxRestartsST = rcp(new StatusTestMaxRestarts_t(defaultMaxRestarts));
  RefCountPtr<StatusTestResNorm_t>
    resNormST = rcp(new StatusTestResNorm_t(defaultResNorm));
  RefCountPtr<StatusTestOutputter_t>
    outputterResNormST = rcp(new StatusTestOutputter_t());
  outputterResNormST->outputFrequency(outputFrequency);
  outputterResNormST->outputMaxResOnly(outputMaxResOnly);
  outputterResNormST->resString("||A*x-b||/||b||");
  outputterResNormST->set_resNormStatusTest(resNormST);
  outputterResNormST->set_outputManager(outputManager);
  RefCountPtr<StatusTestCombo_t>
    maxItersOrRestartsST = rcp(new StatusTestCombo_t(StatusTestCombo_t::OR,*maxItersST,*maxRestartsST));
  set_extra_data(maxItersST,"maxItersST",&maxItersOrRestartsST);
  set_extra_data(maxRestartsST,"maxRestartsST",&maxItersOrRestartsST);
  RefCountPtr<StatusTestCombo_t>
    comboST = rcp(new StatusTestCombo_t(StatusTestCombo_t::OR,*maxItersOrRestartsST,*outputterResNormST));
  set_extra_data(maxItersOrRestartsST,"maxItersOrRestartsST",&comboST);
  set_extra_data(outputterResNormST,"resNormST",&comboST);
  //
  // Generate the solver
  //
  bool useGmres = true;
  if(paramList_.get()) {
    const std::string &solverType = paramList_->get(SolverType_name,GMRES_name);
    if(solverType==GMRES_name)
      useGmres = true;
    else if(solverType=="CG")
      useGmres = false;
    else
      TEST_FOR_EXCEPTION(
        true,std::logic_error
        ,"Error, \"Solver Type\" = \"" << solverType << "\" is not recognized!"
        "  Valid values are \"GMRES\" and \"CG\""
        );
  }
  typedef Belos::IterativeSolver<Scalar,MV_t,LO_t> IterativeSolver_t;
  RefCountPtr<IterativeSolver_t> iterativeSolver = Teuchos::null;
  if(useGmres) {
    RefCountPtr<Teuchos::ParameterList>
      gmresPL;
    if(paramList_.get()) {
      gmresPL = Teuchos::sublist(paramList_,GMRES_name);
    }
    else {
      gmresPL = Teuchos::rcp(new Teuchos::ParameterList()); // BlockGmres requires a PL!
    }
    iterativeSolver = rcp(new Belos::BlockGmres<Scalar,MV_t,LO_t>(lp,comboST,outputManager,gmresPL));
  }
  else {
    iterativeSolver = rcp(new Belos::BlockCG<Scalar,MV_t,LO_t>(lp,comboST,outputManager));
  }
  //
  // Initialize the LOWS object
  //
  const bool adjustableBlockSize = ( paramList_.get() ? paramList_->get(AdjustableBlockSize_name,bool(true)) : true );
  belosOp->initialize(lp,adjustableBlockSize,resNormST,iterativeSolver,outputManager,prec,false,Teuchos::null,supportSolveUse);
#ifdef _DEBUG
  if(paramList_.get()) {
    // Make sure we read the list correctly
    paramList_->validateParameters(*this->getValidParameters(),1); // Validate 0th and 1st level deep
  }
#endif
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &fwdOp
  ,LinearOpWithSolveBase<Scalar>                             *Op
  ) const
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
bool BelosLinearOpWithSolveFactory<Scalar>::supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const
{
  TEST_FOR_EXCEPT(true);
  return false;
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >             &fwdOp
  ,const Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >      &prec
  ,LinearOpWithSolveBase<Scalar>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::initializeApproxPreconditionedOp(
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >             &fwdOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >            &approxFwdOp
  ,LinearOpWithSolveBase<Scalar>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::uninitializeOp(
  LinearOpWithSolveBase<Scalar>                               *Op
  ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >          *fwdOp
  ,Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >    *prec
  ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >          *approxFwdOp
  ,ESupportSolveUse                                           *supportSolveUse
  ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(Op==NULL);
#endif
  BelosLinearOpWithSolve<Scalar>
    &belosOp = Teuchos::dyn_cast<BelosLinearOpWithSolve<Scalar> >(*Op);
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > 
    _fwdOp = belosOp.extract_fwdOp();
  Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >
    _prec = ( belosOp.isExternalPrec() ? belosOp.extract_prec() : Teuchos::null );
  // Note: above we only extract the preconditioner if it was passed in
  // externally.  Otherwise, we need to hold on to it so that we can reuse it
  // in the next initialization.
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
    _approxFwdOp = belosOp.extract_approxFwdOp();
  ESupportSolveUse
    _supportSolveUse = belosOp.supportSolveUse();
  if(fwdOp) *fwdOp = _fwdOp;
  if(prec) *prec = _prec;
  if(approxFwdOp) *approxFwdOp = _approxFwdOp;
  if(supportSolveUse) *supportSolveUse = _supportSolveUse;
}

// Overridden from ParameterListAcceptor

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  paramList->validateParameters(*this->getValidParameters(),1); // Validate 0th and 1st level deep
  paramList_ = paramList;
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::getParameterList()
{
  return paramList_;
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return paramList_;
}

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList>
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
  oss << "{";
  // ToDo: Fill this in some!
  oss << "}";
  return oss.str();
}

// private

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::generateAndGetValidParameters()
{
  static Teuchos::RefCountPtr<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("BelosLinearOpWithSolveFactory"));
    validParamList->set(SolverType_name,GMRES_name);
    validParamList->set(MaxIters_name,int(400));
    validParamList->set(MaxRestarts_name,int(25));
    validParamList->set(BlockSize_name,int(1));
    validParamList->set(AdjustableBlockSize_name,bool(true));
    validParamList->set(DefaultRelResNorm_name,double(1e-6));
    Teuchos::ParameterList
      &gmresSL = validParamList->sublist(GMRES_name);
    gmresSL.set(GMRES_Length_name,int(25));
    gmresSL.set(GMRES_Variant_name,"Standard");
    Teuchos::ParameterList
      &outputterSL = validParamList->sublist(Outputter_name);
    outputterSL.set(Outputter_OutputFrequency_name,int(0));
    outputterSL.set(Outputter_OutputMaxResOnly_name,bool(true));
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
    Teuchos::RefCountPtr<const Teuchos::ParameterList>
      precFactoryValidParamList = precFactory_->getValidParameters();
    if(precFactoryValidParamList.get()) {
      //thisValidParamList_->set(Preconditioner_name,precFactoryName_);
      thisValidParamList_->sublist(precFactoryName_).setParameters(*precFactoryValidParamList);
    }
  }
}

} // namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#endif // __sun
