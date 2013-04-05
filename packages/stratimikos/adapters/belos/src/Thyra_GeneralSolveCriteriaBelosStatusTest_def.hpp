
#ifndef THYRA_GENERAL_SOLVE_CRITERIA_BELOS_STATUS_TEST_DEF_HPP
#define THYRA_GENERAL_SOLVE_CRITERIA_BELOS_STATUS_TEST_DEF_HPP

#include "Thyra_GeneralSolveCriteriaBelosStatusTest.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "BelosThyraAdapter.hpp"
#include "Teuchos_DebugDefaultAsserts.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_as.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
GeneralSolveCriteriaBelosStatusTest<Scalar>::GeneralSolveCriteriaBelosStatusTest()
  :convergenceTestFrequency_(-1),
   compute_x_(false),
   compute_r_(false),
   lastCurrIter_(-1),
   lastRtnStatus_(Belos::Undefined)
{
  GeneralSolveCriteriaBelosStatusTest<Scalar>::reset();
}


template<class Scalar>
void GeneralSolveCriteriaBelosStatusTest<Scalar>::setSolveCriteria(
  const SolveCriteria<Scalar> &solveCriteria,
  const int convergenceTestFrequency
  )
{

  using Teuchos::as;
  typedef ScalarTraits<ScalarMag> SMT;

  // Make sure the solve criteria is valid

  TEUCHOS_ASSERT_INEQUALITY(solveCriteria.requestedTol, >=, SMT::zero());
  TEUCHOS_ASSERT_INEQUALITY(solveCriteria.solveMeasureType.numerator, !=, SOLVE_MEASURE_ONE);
  TEUCHOS_ASSERT(nonnull(solveCriteria.numeratorReductionFunc) ||
    nonnull(solveCriteria.denominatorReductionFunc) );

  // Remember the solve criteria sett

  solveCriteria_ = solveCriteria;
  convergenceTestFrequency_ = convergenceTestFrequency;

  // Determine what needs to be computed on each check

  compute_r_ = solveCriteria.solveMeasureType.contains(SOLVE_MEASURE_NORM_RESIDUAL);

  compute_x_ = (compute_r_ ||
    solveCriteria.solveMeasureType.contains(SOLVE_MEASURE_NORM_SOLUTION));

}


template<class Scalar>
ArrayView<const typename ScalarTraits<Scalar>::magnitudeType>
GeneralSolveCriteriaBelosStatusTest<Scalar>::achievedTol() const
{
  return lastAchievedTol_;
}


// Overridden public functions from Belos::StatusTest


template <class Scalar>
Belos::StatusType
GeneralSolveCriteriaBelosStatusTest<Scalar>::checkStatus(
  Belos::Iteration<Scalar,MV,OP> *iSolver
  )
{

  using Teuchos::null;
  
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(iSolver);
#endif

  const int currIter = iSolver->getNumIters();

  if (currIter == 0 || currIter % convergenceTestFrequency_ != 0) {
    // We will skip this check!
    return Belos::Undefined;
  }

  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  const Belos::LinearProblem<Scalar,MV,OP>& lp = iSolver->getProblem();
  const int numRhs = lp.getRHS()->domain()->dim();
  
  // Compute the rhs norm if requested and not already computed
  if (solveCriteria_.solveMeasureType.contains(SOLVE_MEASURE_NORM_RHS)) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "ToDo: Handle ||b||");
  }

  // Compute the initial residual norm if requested and not already computed
  if (solveCriteria_.solveMeasureType.contains(SOLVE_MEASURE_NORM_INIT_RESIDUAL)) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "ToDo: Handle ||r0||");
  }

  // Compute X if requested
  RCP<MultiVectorBase<Scalar> > X;
  if (compute_x_) {
    RCP<MV> X_update = iSolver->getCurrentUpdate();
    X = lp.updateSolution(X_update);
  }

  // Compute R if requested
  RCP<MultiVectorBase<Scalar> > R;
  if (compute_r_) {
    R = createMembers(lp.getOperator()->range(), X->domain());
    lp.computeCurrResVec(&*R, &*X);
  }

  // Form numerators and denominaotors gN(vN)/gD(vD) for each system

  lastNumerator_.resize(numRhs);
  lastDenominator_.resize(numRhs);

  for (int j = 0; j < numRhs; ++j) {
    const RCP<const VectorBase<Scalar> > x_j = (nonnull(X) ? X->col(j) : null);
    const RCP<const VectorBase<Scalar> > r_j = (nonnull(R) ? R->col(j) : null);
    lastNumerator_[j] = computeReductionFunctional(
      solveCriteria_.solveMeasureType.numerator,
      solveCriteria_.numeratorReductionFunc.ptr(),
      x_j.ptr(), r_j.ptr() );
    lastDenominator_[j] = computeReductionFunctional(
      solveCriteria_.solveMeasureType.denominator,
      solveCriteria_.denominatorReductionFunc.ptr(),
      x_j.ptr(), r_j.ptr() );
  }

  // Determine if convRatio <= requestedTol

  bool systemsAreConverged = true;
  lastAchievedTol_.resize(numRhs);

  for (int j = 0; j < numRhs; ++j) {
    const ScalarMag convRatio = lastNumerator_[j] / lastDenominator_[j]; 
    lastAchievedTol_[j] = convRatio;
    const bool sys_converged_j = (convRatio <= solveCriteria_.requestedTol);
    if (includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
      printRhsStatus(currIter, j, *out);
    }
    if (!sys_converged_j) {
      systemsAreConverged = false;
    }
  }
  
  lastRtnStatus_ = (systemsAreConverged ? Belos::Passed : Belos::Failed);
  lastCurrIter_ = currIter;

  return lastRtnStatus_;

}

template <class Scalar>
Belos::StatusType
GeneralSolveCriteriaBelosStatusTest<Scalar>::getStatus() const
{
  return lastRtnStatus_;
}


template <class Scalar>
void GeneralSolveCriteriaBelosStatusTest<Scalar>::reset()
{
  r0_nrm_.clear();
  b_nrm_.clear();
  lastNumerator_.clear();
  lastDenominator_.clear();
  lastAchievedTol_.clear();
  lastCurrIter_ = -1;
  lastRtnStatus_ = Belos::Undefined;
}


template <class Scalar>
void GeneralSolveCriteriaBelosStatusTest<Scalar>::print(
  std::ostream& os, int indent
  ) const
{
  const int numRhs = lastNumerator_.size();
  for (int j = 0; j < numRhs; ++j) {
    printRhsStatus(lastCurrIter_, j, os, indent);
  }
}


// private


template <class Scalar>
typename GeneralSolveCriteriaBelosStatusTest<Scalar>::ScalarMag
GeneralSolveCriteriaBelosStatusTest<Scalar>::computeReductionFunctional(
  ESolveMeasureNormType measureType,
  const Ptr<const ReductionFunctional<Scalar> > &reductFunc,
  const Ptr<const VectorBase<Scalar> > &x,
  const Ptr<const VectorBase<Scalar> > &r
  ) const
{
  typedef ScalarTraits<ScalarMag> SMT;
  ScalarMag rtn = -SMT::one();
  Ptr<const VectorBase<Scalar> > v;
  switch(measureType) {
    case SOLVE_MEASURE_ONE:
      rtn = SMT::one();
      break;
    case SOLVE_MEASURE_NORM_RESIDUAL:
      v = r;
      break;
    case SOLVE_MEASURE_NORM_SOLUTION:
      v = x;
      break;
    case SOLVE_MEASURE_NORM_INIT_RESIDUAL:
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "ToDo: Handle ||r0||!)")
      break;
    case SOLVE_MEASURE_NORM_RHS:
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "ToDo: Handle ||b||!)");
      break;
    TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT();
   }
  if (rtn >= SMT::zero()) {
    // We already have what we need!
  }
  else if (nonnull(v) && rtn < SMT::zero()) {
    if (nonnull(reductFunc)) {
      rtn = reductFunc->reduce(*v);
    }
    else {
      rtn = norm(*v);
    }
  }
  TEUCHOS_IF_ELSE_DEBUG_ASSERT();
  return rtn;
}


template <class Scalar>
void
GeneralSolveCriteriaBelosStatusTest<Scalar>::printRhsStatus(
  const int currIter, const int j, std::ostream &out,
  int indent
  ) const
{
  const ScalarMag convRatio = lastNumerator_[j] / lastDenominator_[j]; 
  const bool sys_converged_j = (convRatio <= solveCriteria_.requestedTol);
  for (int i = 0; i < indent; ++i) { out << " "; }
  out
    << "["<<currIter<<"] "
    << "gN(vN("<<j<<"))/gD(vD("<<j<<")) = "
    << lastNumerator_[j] << "/" << lastDenominator_[j] << " = "
    << convRatio << " <= " << solveCriteria_.requestedTol << " : "
    << (sys_converged_j ? " true" : "false")
    << "\n";
}


} // namespace Thyra


#endif	// THYRA_GENERAL_SOLVE_CRITERIA_BELOS_STATUS_TEST_DEF_HPP
