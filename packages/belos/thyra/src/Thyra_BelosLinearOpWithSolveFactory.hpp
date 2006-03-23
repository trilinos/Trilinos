
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
const std::string BelosLinearOpWithSolveFactory<Scalar>::DefaultRelResNorm_name = "Default Rel Res Norm";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_name = "GMRES";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_Length_name = "Length";
template<class Scalar>
const std::string BelosLinearOpWithSolveFactory<Scalar>::GMRES_Variant_name = "Variant";

// Constructors/initializers/accessors

// Overridden from LinearOpWithSolveFactoryBase

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
  ,const ESupportSolveUse                                     supportSolveUse
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
  // Set the block size
  //
  int blockSize = 1;
  if(paramList_.get()) {
    blockSize = paramList_->get(BlockSize_name,blockSize);
  }
  lp->SetBlockSize(blockSize);
  //
  // Create the default status test
  //
  int         defaultMaxIterations = 400;
  int         defaultMaxRestarts = 25;
  ScalarMag   defaultResNorm = 1e-6;
  if(paramList_.get()) {
    defaultMaxIterations = paramList_->get(MaxIters_name,defaultMaxIterations);
    defaultMaxRestarts = paramList_->get(MaxRestarts_name,defaultMaxRestarts);
    defaultResNorm = paramList_->get(DefaultRelResNorm_name,defaultResNorm);
  }
  //
  typedef Belos::StatusTestResNorm<Scalar,MV_t,LO_t>      StatusTestResNorm_t;
  typedef Belos::StatusTestMaxIters<Scalar,MV_t,LO_t>     StatusTestMaxIters_t;
  typedef Belos::StatusTestMaxRestarts<Scalar,MV_t,LO_t>  StatusTestMaxRestarts_t;
  typedef Belos::StatusTestCombo<Scalar,MV_t,LO_t>        StatusTestCombo_t;
  RefCountPtr<StatusTestMaxIters_t>
    maxItersST = rcp(new StatusTestMaxIters_t(defaultMaxIterations));
  RefCountPtr<StatusTestMaxRestarts_t>
    maxRestartsST = rcp(new StatusTestMaxRestarts_t(defaultMaxRestarts));
  RefCountPtr<StatusTestResNorm_t>
    resNormST = rcp(new StatusTestResNorm_t(defaultResNorm));
  RefCountPtr<StatusTestCombo_t>
    maxItersOrRestartsST = rcp(new StatusTestCombo_t(StatusTestCombo_t::OR,*maxItersST,*maxRestartsST));
  set_extra_data(maxItersST,"maxItersST",&maxItersOrRestartsST);
  set_extra_data(maxRestartsST,"maxRestartsST",&maxItersOrRestartsST);
  RefCountPtr<StatusTestCombo_t>
    comboST = rcp(new StatusTestCombo_t(StatusTestCombo_t::OR,*maxItersOrRestartsST,*resNormST));
  set_extra_data(maxItersOrRestartsST,"maxItersOrRestartsST",&comboST);
  set_extra_data(resNormST,"resNormST",&comboST);
  //
  // Create the output manager 
  //
  typedef Belos::OutputManager<Scalar> OutputManager_t;
  int belosVerbLevel = Belos::Warnings | Belos::FinalSummary;
  // ToDo: Set these from the verbosity level
  RefCountPtr<Teuchos::FancyOStream>
    ostream = this->getOStream();
  RefCountPtr<OutputManager_t>
    outputManager = rcp(new OutputManager_t(0,belosVerbLevel));
  // Note: The stream itself will be set in the BelosLinearOpWithSolve object!
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
      gmresPL = Teuchos::rcp(&paramList_->sublist(GMRES_name),false);
      set_extra_data(paramList_,"topList",&gmresPL);
    }
    iterativeSolver = rcp(new Belos::BlockGmres<Scalar,MV_t,LO_t>(lp,comboST,outputManager,gmresPL));
  }
  else {
    iterativeSolver = rcp(new Belos::BlockCG<Scalar,MV_t,LO_t>(lp,comboST,outputManager));
  }
  //
  // Initialize the LOWS object
  //
  belosOp->initialize(lp,resNormST,iterativeSolver,outputManager);
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
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >     &fwdOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &precOp
  ,const EPreconditionerInputType                             precOpType
  ,LinearOpWithSolveBase<Scalar>                              *Op
  ,const ESupportSolveUse                                     supportSolveUse
  ) const
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void BelosLinearOpWithSolveFactory<Scalar>::uninitializeOp(
  LinearOpWithSolveBase<Scalar>                       *Op
  ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >  *fwdOp
  ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >  *precOp
  ,EPreconditionerInputType                           *precOpType
  ,ESupportSolveUse                                   *supportSolveUse
  ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(Op==NULL);
#endif
  TEST_FOR_EXCEPT(true);
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
  return generateAndGetValidParameters();
}

// Public functions overridden from Teuchos::Describable

template<class Scalar>
std::string BelosLinearOpWithSolveFactory<Scalar>::description() const
{
  std::ostringstream oss;
  TEST_FOR_EXCEPT(true);
  return oss.str();
}

// private

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList>
BelosLinearOpWithSolveFactory<Scalar>::generateAndGetValidParameters()
{
  static Teuchos::RefCountPtr<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("Thyra::BelosLinearOpWithSolveFactory"));
    validParamList->set(SolverType_name,GMRES_name);
    validParamList->set(MaxIters_name,int(400));
    validParamList->set(MaxRestarts_name,int(25));
    validParamList->set(BlockSize_name,1);
    validParamList->set(DefaultRelResNorm_name,double(1e-6));
    Teuchos::ParameterList
      &gmresSL = validParamList->sublist(GMRES_name);
    gmresSL.set(GMRES_Length_name,int(25));
    gmresSL.set(GMRES_Variant_name,"Standard");
  }
  return validParamList;
}

} // namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#endif // __sun
