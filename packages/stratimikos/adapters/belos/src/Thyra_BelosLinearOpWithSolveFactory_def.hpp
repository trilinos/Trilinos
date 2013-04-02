
#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP


#include "Thyra_BelosLinearOpWithSolveFactory_decl.hpp"
#include "Thyra_BelosLinearOpWithSolve.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosMinresSolMgr.hpp"
#include "BelosThyraAdapter.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StandardValidatorXMLConverters.hpp"


namespace Thyra {


// Parameter names for Paramter List

template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::SolverType_name = "Solver Type";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::SolverType_default = "Pseudo Block GMRES";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::SolverTypes_name = "Solver Types";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::BlockGMRES_name = "Block GMRES";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::PseudoBlockGMRES_name = "Pseudo Block GMRES";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::BlockCG_name = "Block CG";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::PseudoBlockCG_name = "Pseudo Block CG";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GCRODR_name = "GCRODR";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::RCG_name = "RCG";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::MINRES_name = "MINRES";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::ConvergenceTestFrequency_name = "Convergence Test Frequency";

// Constructors/initializers/accessors


template<class Scalar>
BelosLinearOpWithSolveFactory<Scalar>::BelosLinearOpWithSolveFactory()
  :solverType_(SOLVER_TYPE_PSEUDO_BLOCK_GMRES),
   convergenceTestFrequency_(1)
{
  updateThisValidParamList();
}


template<class Scalar>
BelosLinearOpWithSolveFactory<Scalar>::BelosLinearOpWithSolveFactory(
  const RCP<PreconditionerFactoryBase<Scalar> > &precFactory
  )
  :solverType_(SOLVER_TYPE_PSEUDO_BLOCK_GMRES)
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
  paramList->validateParametersAndSetDefaults(*this->getValidParameters(), 1);
  paramList_ = paramList;
  solverType_ =
    Teuchos::getIntegralValue<EBelosSolverType>(*paramList_, SolverType_name);
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
  using Teuchos::setStringToIntegralParameter;
Teuchos::ValidatorXMLConverterDB::addConverter(
  Teuchos::DummyObjectGetter<
    Teuchos::StringToIntegralParameterEntryValidator<EBelosSolverType> 
  >::getDummyObject(),
  Teuchos::DummyObjectGetter<Teuchos::StringToIntegralValidatorXMLConverter<
    EBelosSolverType> >::getDummyObject());

  typedef MultiVectorBase<Scalar> MV_t;
  typedef LinearOpBase<Scalar> LO_t;
  static RCP<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("BelosLinearOpWithSolveFactory"));
    setStringToIntegralParameter<EBelosSolverType>(
      SolverType_name, SolverType_default,
      "Type of linear solver algorithm to use.",
      tuple<std::string>(
        "Block GMRES",
        "Pseudo Block GMRES",
        "Block CG",
        "Pseudo Block CG",
        "GCRODR",
	"RCG",
        "MINRES"
        ),
      tuple<std::string>(
        "Block GMRES solver for nonsymmetric linear systems.  It can also solve "
	"single right-hand side systems, and can also perform Flexible GMRES "
	"(where the preconditioner may change at every iteration, for example "
	"for inner-outer iterations), by setting options in the \"Block GMRES\" "
	"sublist.",

        "GMRES solver for nonsymmetric linear systems, that performs single "
	"right-hand side solves on multiple right-hand sides at once.  It "
	"exploits operator multivector multiplication in order to amortize "
        "global communication costs.  Individual linear systems are deflated "
	"out as they are solved.",

        "Block CG solver for symmetric (Hermitian in complex arithmetic) "
	"positive definite linear systems.  It can also solve single "
	"right-hand-side systems.",

        "CG solver that performs single right-hand side CG on multiple right-hand "
	"sides at once.  It exploits operator multivector multiplication in order "
	"to amortize global communication costs.  Individual linear systems are "
	"deflated out as they are solved.",

        "GMRES solver for nonsymmetric linear systems, that performs subspace "
	"recycling to accelerate convergence for sequences of related linear "
	"systems.",

	"CG solver for symmetric (Hermitian in complex arithmetic) positive "
	"definite linear systems, that performs subspace recycling to "
	"accelerate convergence for sequences of related linear systems.",

        "MINRES solver for symmetric indefinite linear systems.  It performs "
	"single-right-hand-side solves on multiple right-hand sides sequentially."
        ),
      tuple<EBelosSolverType>(
        SOLVER_TYPE_BLOCK_GMRES,
        SOLVER_TYPE_PSEUDO_BLOCK_GMRES,
        SOLVER_TYPE_BLOCK_CG,
        SOLVER_TYPE_PSEUDO_BLOCK_CG,
        SOLVER_TYPE_GCRODR,
	SOLVER_TYPE_RCG,
        SOLVER_TYPE_MINRES
        ),
      &*validParamList
      );
    validParamList->set(ConvergenceTestFrequency_name, as<int>(1),
      "Number of linear solver iterations to skip between applying"
      " user-defined convergence test.");
    Teuchos::ParameterList
      &solverTypesSL = validParamList->sublist(SolverTypes_name);
    {
      Belos::BlockGmresSolMgr<Scalar,MV_t,LO_t> mgr;
      solverTypesSL.sublist(BlockGMRES_name).setParameters(
        *mgr.getValidParameters()
        );
    }
    {
      Belos::PseudoBlockGmresSolMgr<Scalar,MV_t,LO_t> mgr;
      solverTypesSL.sublist(PseudoBlockGMRES_name).setParameters(
        *mgr.getValidParameters()
        );
    }
    {
      Belos::BlockCGSolMgr<Scalar,MV_t,LO_t> mgr;
      solverTypesSL.sublist(BlockCG_name).setParameters(
        *mgr.getValidParameters()
        );
    }
    {
      Belos::PseudoBlockCGSolMgr<Scalar,MV_t,LO_t> mgr;
      solverTypesSL.sublist(PseudoBlockCG_name).setParameters(
        *mgr.getValidParameters()
        );
    }
    {
      Belos::GCRODRSolMgr<Scalar,MV_t,LO_t> mgr;
      solverTypesSL.sublist(GCRODR_name).setParameters(
        *mgr.getValidParameters()
        );
    }
    {
      Belos::RCGSolMgr<Scalar,MV_t,LO_t> mgr;
      solverTypesSL.sublist(RCG_name).setParameters(
        *mgr.getValidParameters()
        );
    }
    {
      Belos::MinresSolMgr<Scalar,MV_t,LO_t> mgr;
      solverTypesSL.sublist(MINRES_name).setParameters(
        *mgr.getValidParameters()
        );
    }
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
    if(unspecified.get()) {
      lp->setRightPrec(unspecified);
      // ToDo: Allow user to determine whether this should be placed on the
      // left or on the right through a parameter in the parameter list!
    }
    else {
      // Set a left, right or split preconditioner
      TEUCHOS_TEST_FOR_EXCEPTION(
        left.get(),std::logic_error
        ,"Error, we can not currently handle a left preconditioner!"
        );
      lp->setRightPrec(right);
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

  typedef Belos::SolverManager<Scalar,MV_t,LO_t> IterativeSolver_t;
  RCP<IterativeSolver_t> iterativeSolver = Teuchos::null;
  RCP<Teuchos::ParameterList> solverPL = Teuchos::rcp( new Teuchos::ParameterList() );
  
  switch(solverType_) {
    case SOLVER_TYPE_BLOCK_GMRES: 
    {
      // Set the PL
      if(paramList_.get()) {
        Teuchos::ParameterList &solverTypesPL = paramList_->sublist(SolverTypes_name);
        Teuchos::ParameterList &gmresPL = solverTypesPL.sublist(BlockGMRES_name);
        solverPL = Teuchos::rcp( &gmresPL, false );
      }
      // Create the solver
      if (oldIterSolver != Teuchos::null) {
        iterativeSolver = oldIterSolver;
        iterativeSolver->setProblem( lp );
        iterativeSolver->setParameters( solverPL );
      } 
      else {
        iterativeSolver = rcp(new Belos::BlockGmresSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
      }
      break;
    }
    case SOLVER_TYPE_PSEUDO_BLOCK_GMRES:
    {
      // Set the PL
      if(paramList_.get()) {
        Teuchos::ParameterList &solverTypesPL = paramList_->sublist(SolverTypes_name);
        Teuchos::ParameterList &pbgmresPL = solverTypesPL.sublist(PseudoBlockGMRES_name);
        solverPL = Teuchos::rcp( &pbgmresPL, false );
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
        iterativeSolver = rcp(new Belos::PseudoBlockGmresSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
      }
      break;
    }
    case SOLVER_TYPE_BLOCK_CG:
    {
      // Set the PL
      if(paramList_.get()) {
        Teuchos::ParameterList &solverTypesPL = paramList_->sublist(SolverTypes_name);
        Teuchos::ParameterList &cgPL = solverTypesPL.sublist(BlockCG_name);
        solverPL = Teuchos::rcp( &cgPL, false );
      }
      // Create the solver
      if (oldIterSolver != Teuchos::null) {
        iterativeSolver = oldIterSolver;
        iterativeSolver->setProblem( lp );
        iterativeSolver->setParameters( solverPL );
      }
      else {
        iterativeSolver = rcp(new Belos::BlockCGSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
      }
      break;
    }
    case SOLVER_TYPE_PSEUDO_BLOCK_CG:
    {
      // Set the PL
      if(paramList_.get()) {
        Teuchos::ParameterList &solverTypesPL = paramList_->sublist(SolverTypes_name);
        Teuchos::ParameterList &pbcgPL = solverTypesPL.sublist(PseudoBlockCG_name);
        solverPL = Teuchos::rcp( &pbcgPL, false );
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
        iterativeSolver = rcp(new Belos::PseudoBlockCGSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
      }
      break;
    }
    case SOLVER_TYPE_GCRODR:
    {
      // Set the PL
      if(paramList_.get()) {
        Teuchos::ParameterList &solverTypesPL = paramList_->sublist(SolverTypes_name);
        Teuchos::ParameterList &gcrodrPL = solverTypesPL.sublist(GCRODR_name);
        solverPL = Teuchos::rcp( &gcrodrPL, false );
      }
      // Create the solver
      if (oldIterSolver != Teuchos::null) {
        iterativeSolver = oldIterSolver;
        iterativeSolver->setProblem( lp );
        iterativeSolver->setParameters( solverPL );
      } 
      else {
        iterativeSolver = rcp(new Belos::GCRODRSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
      }
      break;
    }
    case SOLVER_TYPE_RCG:
    {
      // Set the PL
      if(paramList_.get()) {
        Teuchos::ParameterList &solverTypesPL = paramList_->sublist(SolverTypes_name);
        Teuchos::ParameterList &rcgPL = solverTypesPL.sublist(RCG_name);
        solverPL = Teuchos::rcp( &rcgPL, false );
      }
      // Create the solver
      if (oldIterSolver != Teuchos::null) {
        iterativeSolver = oldIterSolver;
        iterativeSolver->setProblem( lp );
        iterativeSolver->setParameters( solverPL );
      } 
      else {
        iterativeSolver = rcp(new Belos::RCGSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
      }
      break;
    }
    case SOLVER_TYPE_MINRES:
    {
      // Set the PL
      if(paramList_.get()) {
        Teuchos::ParameterList &solverTypesPL = paramList_->sublist(SolverTypes_name);
        Teuchos::ParameterList &minresPL = solverTypesPL.sublist(MINRES_name);
        solverPL = Teuchos::rcp( &minresPL, false );
      }
      // Create the solver
      if (oldIterSolver != Teuchos::null) {
        iterativeSolver = oldIterSolver;
        iterativeSolver->setProblem( lp );
        iterativeSolver->setParameters( solverPL );
      }
      else {
        iterativeSolver = rcp(new Belos::MinresSolMgr<Scalar,MV_t,LO_t>(lp,solverPL));
      }
      break;
    }

    default:
    {
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
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
