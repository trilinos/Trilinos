
#ifndef __sun

#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_BelosLinearOpWithSolveDecl.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
BelosLinearOpWithSolve<Scalar>::BelosLinearOpWithSolve()
{}

template<class Scalar>
BelosLinearOpWithSolve<Scalar>::BelosLinearOpWithSolve(
  const Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
  ,const Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    &resNormST
  ,const Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      &iterativeSolver
  ,const Teuchos::RefCountPtr<Belos::OutputManager<Scalar> >                  &outputManager
  )
{
  initialize(lp,resNormST,iterativeSolver,outputManager);
}

template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::initialize(
  const Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
  ,const Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    &resNormST
  ,const Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      &iterativeSolver
  ,const Teuchos::RefCountPtr<Belos::OutputManager<Scalar> >                  &outputManager
  )
{
  this->setLinePrefix("BELOS/T");
  // ToDo: Validate input
  lp_ = lp;
  resNormST_ = resNormST;
  iterativeSolver_ = iterativeSolver;
  outputManager_ = outputManager;
  defaultTol_ = resNormST_->GetTolerance(); // We need to remember this!
}

// Overridden from LinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::range() const
{
  return lp_->GetOperator()->range();
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::domain() const
{
  return lp_->GetOperator()->domain();
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::string BelosLinearOpWithSolve<Scalar>::description() const
{
  std::ostringstream oss;
  oss << "Thyra::BelosLinearOpWithSolve<"<<Teuchos::ScalarTraits<Scalar>::name()<<">";
  if(lp_->GetOperator().get()) {
    oss << "(";
    oss << "fwdOp=\'"<<lp_->GetOperator()->description()<<"\'";
    oss << ",iterativeSolver=\'"<<iterativeSolver_->description()<<"\'";
    oss << ")";
  }
  // ToDo: Make Belos::IterativeSolver derive from Teuchos::Describable so
  // that we can get better information.
  return oss.str();
}

// ToDo: Add more detailed describe() function override to show all of the good stuff!

// protected

// Overridden from SingleScalarLinearOpBase

template<class Scalar>
bool BelosLinearOpWithSolve<Scalar>::opSupported(ETransp M_trans) const
{
  return ::Thyra::opSupported(*lp_->GetOperator(),M_trans);
}

template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  ::Thyra::apply(*lp_->GetOperator(),M_trans,X,Y,alpha,beta);
}

// Overridden from SingleScalarLinearOpWithSolveBase

template<class Scalar>
bool BelosLinearOpWithSolve<Scalar>::solveSupportsTrans(ETransp M_trans) const
{
  if(real_trans(M_trans)==NOTRANS) return true;
  return false; // ToDo: Support adjoint solves!
}

template<class Scalar>
bool BelosLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureType(ETransp M_trans, const SolveMeasureType& solveMeasureType) const
{
  if(real_trans(M_trans)==NOTRANS) {
    return (
      solveMeasureType.useDefault()
      ||
      (solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS) && resNormST_.get())
      );
  }
  // TRANS
  return false; // ToDo: Support adjoint solves!
}

template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::solve(
  const ETransp                         M_trans
  ,const MultiVectorBase<Scalar>        &B
  ,MultiVectorBase<Scalar>              *X
  ,const int                            numBlocks
  ,const BlockSolveCriteria<Scalar>     blockSolveCriteria[]
  ,SolveStatus<Scalar>                  blockSolveStatus[]
  ) const
{
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  
  TEST_FOR_EXCEPT(numBlocks > 1); // ToDo: Deal with multiple solve criteria later if needed

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = this->getOStream();
  Teuchos::EVerbosityLevel
    verbLevel = this->getVerbLevel();

  OSTab tab = this->getOSTab();

  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
    *out << "\nStarting iterations with Belos solver of type \""<<iterativeSolver_->description()<<"\" ...\n\n";
  
  //
  // Set RHS and LHS
  //
  lp_->Reset(Teuchos::rcp(X,false),Teuchos::rcp(&B,false));
  //
  // Set the solution criteria
  //
  SolveMeasureType solveMeasureType;
  if(numBlocks==1) {
    solveMeasureType = blockSolveCriteria[0].solveCriteria.solveMeasureType;
    const ScalarMag requestedTol = blockSolveCriteria[0].solveCriteria.requestedTol;
    TEST_FOR_EXCEPT( !( solveMeasureType.useDefault() || solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS) ) );
    if( solveMeasureType.useDefault() ) {
        resNormST_->ResetTolerance(defaultTol_);
    }
    else if( solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS) ) {
      if( requestedTol != SolveCriteria<Scalar>::unspecifiedTolerance() )
        resNormST_->ResetTolerance(requestedTol);
      else
        resNormST_->ResetTolerance(defaultTol_);
      resNormST_->DefineScaleForm(StatusTestResNorm_t::NormOfRHS,Belos::TwoNorm);
    }
    else {
      TEST_FOR_EXCEPT(true); // Should never get there.
    }
  }
  else {
    resNormST_->ResetTolerance(defaultTol_);
  }
  //
  // Solve the linear system
  //
  iterativeSolver_->GetStatusTest()->Reset(); 
  iterativeSolver_->Reset();
  if(1){
    outputManager_->SetOStream(out);
    Teuchos::OSTab tab(out,1,"BELOS");
    iterativeSolver_->Solve();
  }
  //
  // Report the solve status
  //
  const Belos::StatusType belosSolveStatus = resNormST_->GetStatus();
  const std::vector<ScalarMag>* resNormValues = resNormST_->GetTestValue();
  TEST_FOR_EXCEPT(resNormValues==NULL);
  const ScalarMag belosAchievedTol = *std::max_element(resNormValues->begin(),resNormValues->end());
  TEST_FOR_EXCEPTION(
    belosSolveStatus==Belos::Failed || belosSolveStatus==Belos::NaN, CatastrophicSolveFailure
    ,"Error, something really bad happened in the Belos solver!"
    ); // ToDo: Get Belos to return a more informative error mesage to embed here?
  //
  ESolveStatus solveStatus = SOLVE_STATUS_UNKNOWN;
  ScalarMag    achievedTol = SolveStatus<Scalar>::unknownTolerance();
  if(solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS)) {
    switch(belosSolveStatus) {
      case Belos::Unchecked:
        solveStatus = SOLVE_STATUS_UNKNOWN;
        break;
      case Belos::Unconverged:
        solveStatus = SOLVE_STATUS_UNCONVERGED;
        break;
      case Belos::Converged:
        solveStatus = SOLVE_STATUS_CONVERGED;
        break;
      default:
        TEST_FOR_EXCEPT(true); // Should never get here!
    }
    achievedTol = belosAchievedTol;
  }
  const int iterations = iterativeSolver_->GetNumIters();
  //
  std::ostringstream ossmessage;
  ossmessage
    << "The Belos solver of type \""<<iterativeSolver_->description()<<"\" returned a solve status of \""
    << toString(belosSolveStatus) << "\" in " << iterations << " iterations and achieved an approximate tolerance of "
    << SolveStatus<Scalar>::achievedTolToString(achievedTol);
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
    *out << "\n" << ossmessage.str() << "\n";
  //
  if( numBlocks && blockSolveStatus ) {
    blockSolveStatus[0].solveStatus = solveStatus;
    blockSolveStatus[0].achievedTol = achievedTol;
    blockSolveStatus[0].message     = ossmessage.str();
  }
}

}	// end namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP

#endif // __sun
