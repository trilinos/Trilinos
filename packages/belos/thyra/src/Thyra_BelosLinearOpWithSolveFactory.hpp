
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

// Constructors/initializers/accessors

// Overridden from LinearOpWithSolveFactoryBase

template<class Scalar>
bool BelosLinearOpWithSolveFactory<Scalar>::isCompatible(
  const LinearOpBase<Scalar> &fwdOp
  ) const
{
  //TEST_FOR_EXCEPT(true);
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
  // Create the default status test
  //
  // ToDo: Read these from a parameter list!
  const int         defaultMaxIterations  = 1000;
  const int         defaultMaxRestarts = 25;
  const ScalarMag   defaultResNorm = 1e-6;
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
  RefCountPtr<OutputManager_t>
    outputManager = rcp(new OutputManager_t(0)); // ToDo: Use VerboseObject::getOStream()
  //
  // Get the parameter list 
  //
  RefCountPtr<Teuchos::ParameterList>
    paramList = rcp(new Teuchos::ParameterList()); // ToDo: Grab this from the Teuchos::ParameterListAcceptor interface ...
  //
  // Generate the solver
  //
  const bool useGmres = true; // Determined form parameter list
  typedef Belos::IterativeSolver<Scalar,MV_t,LO_t> IterativeSolver_t;
  RefCountPtr<IterativeSolver_t> iterativeSolver = Teuchos::null;
  if(useGmres) {
    iterativeSolver = rcp(new Belos::BlockGmres<Scalar,MV_t,LO_t>(lp,comboST,outputManager,paramList));
  }
  else {
    iterativeSolver = rcp(new Belos::BlockCG<Scalar,MV_t,LO_t>(lp,comboST,outputManager));
  }
  //
  // Initialize the LOWS object
  //
  belosOp->initialize(lp,resNormST,iterativeSolver);
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

} // namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#endif // __sun
