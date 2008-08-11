
#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_BelosLinearOpWithSolveDecl.hpp"
#include "Thyra_LinearOpWithSolveHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
BelosLinearOpWithSolve<Scalar>::BelosLinearOpWithSolve()
  :isExternalPrec_(false)
  ,supportSolveUse_(SUPPORT_SOLVE_UNSPECIFIED)
  ,defaultTol_(-1.0)
{}

template<class Scalar>
BelosLinearOpWithSolve<Scalar>::BelosLinearOpWithSolve(
  const Teuchos::RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
  ,const Teuchos::RCP<Teuchos::ParameterList>                         &solverPL
  ,const Teuchos::RCP<Belos::SolverManager<Scalar,MV_t,LO_t> >        &iterativeSolver
  ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >              &fwdOpSrc
  ,const Teuchos::RCP<const PreconditionerBase<Scalar> >              &prec
  ,const bool                                                         isExternalPrec
  ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >              &approxFwdOpSrc
  ,const ESupportSolveUse                                             &supportSolveUse
  )
{
  initialize(
    lp,solverPL,iterativeSolver
    ,fwdOpSrc,prec,isExternalPrec,approxFwdOpSrc,supportSolveUse
    );
}

template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::initialize(
  const Teuchos::RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
  ,const Teuchos::RCP<Teuchos::ParameterList>                         &solverPL
  ,const Teuchos::RCP<Belos::SolverManager<Scalar,MV_t,LO_t> >        &iterativeSolver
  ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >              &fwdOpSrc
  ,const Teuchos::RCP<const PreconditionerBase<Scalar> >              &prec
  ,const bool                                                         isExternalPrec
  ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >              &approxFwdOpSrc
  ,const ESupportSolveUse                                             &supportSolveUse
  )
{
  this->setLinePrefix("BELOS/T");
  // ToDo: Validate input
  lp_ = lp;
  solverPL_ = solverPL;
  iterativeSolver_ = iterativeSolver;
  fwdOpSrc_ = fwdOpSrc;
  prec_ = prec;
  isExternalPrec_ = isExternalPrec;
  approxFwdOpSrc_ = approxFwdOpSrc;
  supportSolveUse_ = supportSolveUse;
  // Check if "Convergence Tolerance" is in the solver parameter list.  If not, use the default from the solver.
  if ( !is_null(solverPL_) ) {
    if (solverPL_->isParameter("Convergence Tolerance")) {
      defaultTol_ = solverPL_->get<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>("Convergence Tolerance");
    }
  }
  else {
    Teuchos::RCP<const Teuchos::ParameterList> defaultPL = iterativeSolver->getValidParameters();
    defaultTol_ = defaultPL->get<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>("Convergence Tolerance");
  }
}

template<class Scalar>
Teuchos::RCP<const LinearOpSourceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::extract_fwdOpSrc()
{
  Teuchos::RCP<const LinearOpSourceBase<Scalar> >
    _fwdOpSrc = fwdOpSrc_;
  fwdOpSrc_ = Teuchos::null;
  return _fwdOpSrc;
}

template<class Scalar>
Teuchos::RCP<const PreconditionerBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::extract_prec()
{
  Teuchos::RCP<const PreconditionerBase<Scalar> >
    _prec = prec_;
  prec_ = Teuchos::null;
  return _prec;
}

template<class Scalar>
bool BelosLinearOpWithSolve<Scalar>::isExternalPrec() const
{
  return isExternalPrec_;
}

template<class Scalar>
Teuchos::RCP<const LinearOpSourceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::extract_approxFwdOpSrc()
{
  Teuchos::RCP<const LinearOpSourceBase<Scalar> >
    _approxFwdOpSrc = approxFwdOpSrc_;
  approxFwdOpSrc_ = Teuchos::null;
  return _approxFwdOpSrc;
}

template<class Scalar>
ESupportSolveUse BelosLinearOpWithSolve<Scalar>::supportSolveUse() const
{
  return supportSolveUse_;
}

template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::uninitialize(
  Teuchos::RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> >         *lp
  ,Teuchos::RCP<Teuchos::ParameterList>                         *solverPL
  ,Teuchos::RCP<Belos::SolverManager<Scalar,MV_t,LO_t> >        *iterativeSolver
  ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >              *fwdOpSrc
  ,Teuchos::RCP<const PreconditionerBase<Scalar> >              *prec
  ,bool                                                         *isExternalPrec
  ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >              *approxFwdOpSrc
  ,ESupportSolveUse                                             *supportSolveUse
  )
{
  if(lp) *lp = lp_;
  if(solverPL) *solverPL = solverPL_;
  if(iterativeSolver) *iterativeSolver = iterativeSolver_;
  if(fwdOpSrc) *fwdOpSrc = fwdOpSrc_;
  if(prec) *prec = prec_;
  if(isExternalPrec) *isExternalPrec = isExternalPrec_;
  if(approxFwdOpSrc) *approxFwdOpSrc = approxFwdOpSrc_;
  if(supportSolveUse) *supportSolveUse = supportSolveUse_;

  lp_ = Teuchos::null;
  solverPL_ = Teuchos::null;
  iterativeSolver_ = Teuchos::null;
  fwdOpSrc_ = Teuchos::null;
  prec_ = Teuchos::null;
  isExternalPrec_ = false;
  approxFwdOpSrc_ = Teuchos::null;
  supportSolveUse_ = SUPPORT_SOLVE_UNSPECIFIED;
}

// Overridden from LinearOpBase

template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::range() const
{
  if (!is_null(lp_))
    return lp_->getOperator()->range();
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::domain() const
{
  if (!is_null(lp_))
    return lp_->getOperator()->domain();
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::string BelosLinearOpWithSolve<Scalar>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if( !is_null(lp_) && !is_null(lp_->getOperator()) ) {
    oss << "{";
    oss << "iterativeSolver=\'"<<iterativeSolver_->description()<<"\'";
    oss << ",fwdOp=\'"<<lp_->getOperator()->description()<<"\'";
    if(lp_->getLeftPrec().get())
      oss << ",leftPrecOp=\'"<<lp_->getLeftPrec()->description()<<"\'";
    if(lp_->getRightPrec().get())
      oss << ",rightPrecOp=\'"<<lp_->getRightPrec()->description()<<"\'";
    oss << "}";
  }
  // ToDo: Make Belos::SolverManager derive from Teuchos::Describable so
  // that we can get better information.
  return oss.str();
}

template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RCP;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << Teuchos::Describable::description()<< "{"
        << "rangeDim=" << this->range()->dim()
        << ",domainDim=" << this->domain()->dim() << "}\n";
      if(lp_->getOperator().get()) {
        OSTab tab(out);
        *out
          << "iterativeSolver = "<<describe(*iterativeSolver_,verbLevel)
          << "fwdOp = " << describe(*lp_->getOperator(),verbLevel);
        if(lp_->getLeftPrec().get())
          *out << "leftPrecOp = "<<describe(*lp_->getLeftPrec(),verbLevel);
        if(lp_->getRightPrec().get())
          *out << "rightPrecOp = "<<describe(*lp_->getRightPrec(),verbLevel);
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}

// protected

// Overridden from SingleScalarLinearOpBase

template<class Scalar>
bool BelosLinearOpWithSolve<Scalar>::opSupported(EOpTransp M_trans) const
{
  return ::Thyra::opSupported(*lp_->getOperator(),M_trans);
}

template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::apply(
  const EOpTransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  ::Thyra::apply(*lp_->getOperator(),M_trans,X,Y,alpha,beta);
}

// Overridden from SingleScalarLinearOpWithSolveBase

template<class Scalar>
bool BelosLinearOpWithSolve<Scalar>::solveSupportsTrans(EOpTransp M_trans) const
{
  if(real_trans(M_trans)==NOTRANS) return true;
  return false; // ToDo: Support adjoint solves!
}

template<class Scalar>
bool BelosLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureType(EOpTransp M_trans, const SolveMeasureType& solveMeasureType) const
{
  if(real_trans(M_trans)==NOTRANS) {
    return (
      solveMeasureType.useDefault()
      ||
      solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS)
      ||
      solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_INIT_RESIDUAL)
      );
  }
  // TRANS
  return false; // ToDo: Support adjoint solves!
}

template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::solve(
  const EOpTransp                         M_trans
  ,const MultiVectorBase<Scalar>        &B
  ,MultiVectorBase<Scalar>              *X
  ,const int                            numBlocks
  ,const BlockSolveCriteria<Scalar>     blockSolveCriteria[]
  ,SolveStatus<Scalar>                  blockSolveStatus[]
  ) const
{

  TEUCHOS_FUNC_TIME_MONITOR("BelosLOWS");

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);
  
  TEST_FOR_EXCEPT(numBlocks > 1); // ToDo: Deal with multiple solve criteria later if needed
  TEST_FOR_EXCEPT(!this->solveSupportsTrans(M_trans)); // ToDo: Support adjoint solves later!

  const int numRhs = B.domain()->dim();
  const int numEquations = B.range()->dim();

  Teuchos::RCP<Teuchos::FancyOStream>  out = this->getOStream();
  Teuchos::EVerbosityLevel                     verbLevel = this->getVerbLevel();
  OSTab tab = this->getOSTab();
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE)) {
    *out << "\nStarting iterations with Belos:\n";
    OSTab tab(out);
    *out << "Using forward operator = " << describe(*fwdOpSrc_->getOp(),verbLevel);
    *out << "Using iterative solver = " << describe(*iterativeSolver_,verbLevel);
    *out << "With #Eqns="<<numEquations<<", #RHSs="<<numRhs<<" ...\n";
  }
  //
  // Set RHS and LHS
  //
  bool ret = lp_->setProblem(Teuchos::rcp(X,false),Teuchos::rcp(&B,false));
  TEST_FOR_EXCEPTION(
    ret == false, CatastrophicSolveFailure
    ,"Error, the Belos::LinearProblem could not be set for the current solve!"
    );
  //
  // Set the solution criteria
  //
  Teuchos::RCP<Teuchos::ParameterList> tmpPL = Teuchos::rcp( new Teuchos::ParameterList() );

  SolveMeasureType solveMeasureType;
  if(numBlocks==1) {
    solveMeasureType = blockSolveCriteria[0].solveCriteria.solveMeasureType;
    const ScalarMag requestedTol = blockSolveCriteria[0].solveCriteria.requestedTol;
    assertSupportsSolveMeasureType(*this,M_trans,solveMeasureType);
    if( solveMeasureType.useDefault() ) {
      tmpPL->set("Convergence Tolerance", defaultTol_);
    }
    else if( solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS) ) {
      if( requestedTol != SolveCriteria<Scalar>::unspecifiedTolerance() )
        tmpPL->set("Convergence Tolerance", requestedTol);
      else
        tmpPL->set("Convergence Tolerance", defaultTol_);
      tmpPL->set("Explicit Residual Scaling", "Norm of RHS");
    }
    else if( solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_INIT_RESIDUAL) ) {
      if( requestedTol != SolveCriteria<Scalar>::unspecifiedTolerance() )
        tmpPL->set("Convergence Tolerance", requestedTol);
      else
        tmpPL->set("Convergence Tolerance", defaultTol_);
      tmpPL->set("Explicit Residual Scaling", "Norm of Initial Residual");
    }
    else {
      TEST_FOR_EXCEPT(true); // Should never get there.
    }
  }
  else {
    tmpPL->set("Convergence Tolerance", defaultTol_);
  }
  //
  // Reset the blocksize if we adding more vectors than half the number of equations,
  // orthogonalization will fail on the first iteration!
  //
  Teuchos::RCP<const Teuchos::ParameterList> solverParams = iterativeSolver_->getCurrentParameters();
  const int currBlockSize = Teuchos::getParameter<int>(*solverParams, "Block Size");
  bool isNumBlocks = false;
  int currNumBlocks = 0;
  if (Teuchos::isParameterType<int>(*solverParams, "Num Blocks")) {
    currNumBlocks = Teuchos::getParameter<int>(*solverParams, "Num Blocks");
    isNumBlocks = true;
  }
  const int newBlockSize = TEUCHOS_MIN(currBlockSize,numEquations/2);
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE) && newBlockSize != currBlockSize)
    *out << "\nAdjusted block size = " << newBlockSize << "\n";
  //
  tmpPL->set("Block Size",newBlockSize);
  //
  // Set the number of Krylov blocks if we are using a GMRES solver, or a solver
  // that recognizes "Num Blocks". Otherwise the solver will throw an error!
  if (isNumBlocks) {
    const int Krylov_length = (currNumBlocks*currBlockSize)/newBlockSize;
    tmpPL->set("Num Blocks",Krylov_length);
  
    if(newBlockSize != currBlockSize) {
      if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
        *out
          << "\nAdjusted max number of Krylov basis blocks = " << Krylov_length << "\n";
    }
  }
  //
  // Solve the linear system
  //
  Belos::ReturnType belosSolveStatus;
  {
    RCP<std::ostream>
      outUsed =
      ( static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE)
        ? out
        : rcp(new FancyOStream(rcp(new Teuchos::oblackholestream())))
        );
    Teuchos::OSTab tab(outUsed,1,"BELOS");
    tmpPL->set("Output Stream", outUsed);
    iterativeSolver_->setParameters( tmpPL );
    belosSolveStatus = iterativeSolver_->solve();
  }
  //
  // Report the solve status
  //
  totalTimer.stop();
  //
  ESolveStatus solveStatus = SOLVE_STATUS_UNKNOWN;
  ScalarMag    achievedTol = SolveStatus<Scalar>::unknownTolerance();
  switch(belosSolveStatus) {
    case Belos::Unconverged:
      solveStatus = SOLVE_STATUS_UNCONVERGED;
      break;
    case Belos::Converged:
      solveStatus = SOLVE_STATUS_CONVERGED;
      achievedTol = tmpPL->get("Convergence Tolerance", defaultTol_);
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
  //const int iterations = iterativeSolver_->getNumIters();

  std::ostringstream ossmessage;
  ossmessage
    << "The Belos solver of type \""<<iterativeSolver_->description()<<"\" returned a solve status of \""
    << toString(solveStatus) << 
    /* "\" in " << iterations << " iterations and achieved an approximate max tolerance of "
       << belosAchievedTol << 
    */
    " with total CPU time of " << totalTimer.totalElapsedTime() << " sec" ;
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
    *out << "\n" << ossmessage.str() << "\n";
  //
  if( numBlocks && blockSolveStatus ) {
    blockSolveStatus[0].solveStatus = solveStatus;
    blockSolveStatus[0].achievedTol = achievedTol;
    blockSolveStatus[0].message     = ossmessage.str();
  }

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal solve time = "<<totalTimer.totalElapsedTime()<<" sec\n";
}

}	// end namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP
