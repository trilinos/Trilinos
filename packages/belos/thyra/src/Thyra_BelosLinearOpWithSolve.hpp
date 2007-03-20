
#ifndef __sun

#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_BelosLinearOpWithSolveDecl.hpp"
#include "Teuchos_Time.hpp"

namespace Thyra {

namespace PrivateUtilities {

template<class LP>
class BelosLinearProblemTempBlockSizeSetter {
public:
  BelosLinearProblemTempBlockSizeSetter( const Teuchos::RefCountPtr<LP> &lp, const int newBlockSize )
    :lp_(lp), oldBlockSize_(lp->GetBlockSize())
    {
      lp_->SetBlockSize(newBlockSize);
    }
  ~BelosLinearProblemTempBlockSizeSetter()
    {
      lp_->SetBlockSize(oldBlockSize_);
    }
private:
  Teuchos::RefCountPtr<LP> lp_;
  int                      oldBlockSize_;
  // Not defined and not to be called
  BelosLinearProblemTempBlockSizeSetter();
  BelosLinearProblemTempBlockSizeSetter(const BelosLinearProblemTempBlockSizeSetter&);
  BelosLinearProblemTempBlockSizeSetter& operator=(const BelosLinearProblemTempBlockSizeSetter&);
};

}

// Constructors/initializers/accessors

template<class Scalar>
BelosLinearOpWithSolve<Scalar>::BelosLinearOpWithSolve()
  :isExternalPrec_(false)
  ,supportSolveUse_(SUPPORT_SOLVE_UNSPECIFIED)
  ,defaultTol_(-1.0)
{}

template<class Scalar>
BelosLinearOpWithSolve<Scalar>::BelosLinearOpWithSolve(
  const Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
  ,const bool                                                                 adjustableBlockSize
  ,const int                                                                  maxNumberOfKrylovVectors
  ,const Teuchos::RefCountPtr<Teuchos::ParameterList>                         &gmresPL
  ,const Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    &resNormST
  ,const Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      &iterativeSolver
  ,const Teuchos::RefCountPtr<Belos::OutputManager<Scalar> >                  &outputManager
  ,const Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              &fwdOpSrc
  ,const Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >              &prec
  ,const bool                                                                 isExternalPrec
  ,const Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              &approxFwdOpSrc
  ,const ESupportSolveUse                                                     &supportSolveUse
  )
{
  initialize(
    lp,adjustableBlockSize,maxNumberOfKrylovVectors,gmresPL,resNormST,iterativeSolver
    ,outputManager,fwdOpSrc,prec,isExternalPrec,approxFwdOpSrc,supportSolveUse
    );
}

template<class Scalar>
void BelosLinearOpWithSolve<Scalar>::initialize(
  const Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
  ,const bool                                                                 adjustableBlockSize
  ,const int                                                                  maxNumberOfKrylovVectors
  ,const Teuchos::RefCountPtr<Teuchos::ParameterList>                         &gmresPL
  ,const Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    &resNormST
  ,const Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      &iterativeSolver
  ,const Teuchos::RefCountPtr<Belos::OutputManager<Scalar> >                  &outputManager
  ,const Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              &fwdOpSrc
  ,const Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >              &prec
  ,const bool                                                                 isExternalPrec
  ,const Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              &approxFwdOpSrc
  ,const ESupportSolveUse                                                     &supportSolveUse
  )
{
  this->setLinePrefix("BELOS/T");
  // ToDo: Validate input
  lp_ = lp;
  adjustableBlockSize_ = adjustableBlockSize;
  maxNumberOfKrylovVectors_ = maxNumberOfKrylovVectors;
  gmresPL_ = gmresPL;
  resNormST_ = resNormST;
  iterativeSolver_ = iterativeSolver;
  outputManager_ = outputManager;
  fwdOpSrc_ = fwdOpSrc;
  prec_ = prec;
  isExternalPrec_ = isExternalPrec;
  approxFwdOpSrc_ = approxFwdOpSrc;
  supportSolveUse_ = supportSolveUse;
  defaultTol_ = resNormST_->GetTolerance(); // We need to remember this!
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::extract_fwdOpSrc()
{
  Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >
    _fwdOpSrc = fwdOpSrc_;
  fwdOpSrc_ = Teuchos::null;
  return _fwdOpSrc;
}

template<class Scalar>
Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::extract_prec()
{
  Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >
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
Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >
BelosLinearOpWithSolve<Scalar>::extract_approxFwdOpSrc()
{
  Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >
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
  Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >         *lp
  ,bool                                                                 *adjustableBlockSize
  ,int                                                                  *maxNumberOfKrylovVectors
  ,Teuchos::RefCountPtr<Teuchos::ParameterList>                         *gmresPL
  ,Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    *resNormST
  ,Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      *iterativeSolver
  ,Teuchos::RefCountPtr<Belos::OutputManager<Scalar> >                  *outputManager
  ,Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              *fwdOpSrc
  ,Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >              *prec
  ,bool                                                                 *isExternalPrec
  ,Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              *approxFwdOpSrc
  ,ESupportSolveUse                                                     *supportSolveUse
  )
{
  TEST_FOR_EXCEPT(true); // ToDo: Implement when needed!
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
  oss << Teuchos::Describable::description();
  if(lp_->GetOperator().get()) {
    oss << "{";
    oss << "iterativeSolver=\'"<<iterativeSolver_->description()<<"\'";
    oss << ",fwdOp=\'"<<lp_->GetOperator()->description()<<"\'";
    if(lp_->GetLeftPrec().get())
      oss << ",leftPrecOp=\'"<<lp_->GetLeftPrec()->description()<<"\'";
    if(lp_->GetRightPrec().get())
      oss << ",rightPrecOp=\'"<<lp_->GetRightPrec()->description()<<"\'";
    oss << "}";
  }
  // ToDo: Make Belos::IterativeSolver derive from Teuchos::Describable so
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
  using Teuchos::RefCountPtr;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  RefCountPtr<FancyOStream> out = rcp(&out_arg,false);
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
      if(lp_->GetOperator().get()) {
        OSTab tab(out);
        *out
          << "iterativeSolver = "<<describe(*iterativeSolver_,verbLevel)
          << "fwdOp = " << describe(*lp_->GetOperator(),verbLevel);
        if(lp_->GetLeftPrec().get())
          *out << "leftPrecOp = "<<describe(*lp_->GetLeftPrec(),verbLevel);
        if(lp_->GetRightPrec().get())
          *out << "rightPrecOp = "<<describe(*lp_->GetRightPrec(),verbLevel);
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
  using Teuchos::RefCountPtr;
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

  Teuchos::RefCountPtr<Teuchos::FancyOStream>  out = this->getOStream();
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
  // Temporarily reset the blocksize if we are allowed to due so
  //
  const int currBlockSize = lp_->GetBlockSize();
  const int newBlockSize =
    ( adjustableBlockSize_
      ? TEUCHOS_MIN(numRhs,TEUCHOS_MIN(currBlockSize,numEquations/2))
      // Above: Don't add more vectors than are in the RHS and don't add more
      // vectors than half the number of equations or orthogonalization will
      // fail on the first iteration!
      : currBlockSize
      // Above: Just use the user-set block size if we are not allowed to
      // ajust this!
      );
  // Note: Above, we may also need to adjust the block size to be an even
  // multiple of the number of RHS.  Also, it would be good if
  // Belos::LinearProblem itself allowed for dynamic resizing of the block
  // size between restarts.
  PrivateUtilities::BelosLinearProblemTempBlockSizeSetter<Belos::LinearProblem<Scalar,MV_t,LO_t> >
    tempBlockSizeSetter(lp_,newBlockSize);
  if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
    *out << "\nAdjusted block size = " << newBlockSize << "\n";
  //
  // Set the number of Krylov blocks if we are using GMRES!
  //
  if(gmresPL_.get()) {
    const int BlockGmres_length = maxNumberOfKrylovVectors_/newBlockSize; 
    gmresPL_->set("Length",BlockGmres_length);
    if(out.get() && static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE))
      *out
        << "\nAdjusted max number of GMRES Krylov basis blocks (maxNumberOfKrylovVectors/blockSize) = ("
        << maxNumberOfKrylovVectors_ << ")/(" << newBlockSize << ") = " << BlockGmres_length << "\n";
  }
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
  {
    RefCountPtr<FancyOStream>
      outUsed =
      ( static_cast<int>(verbLevel) > static_cast<int>(Teuchos::VERB_NONE)
        ? out
        : rcp(new FancyOStream(rcp(new Teuchos::oblackholestream())))
        );
    outputManager_->SetOStream(outUsed);
    Teuchos::OSTab tab(outUsed,1,"BELOS");
    iterativeSolver_->Solve();
  }
  //
  // Report the solve status
  //
  totalTimer.stop();
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
    << toString(belosSolveStatus) << "\" in " << iterations << " iterations and achieved an approximate max tolerance of "
    << belosAchievedTol << " with total CPU time of " << totalTimer.totalElapsedTime() << " sec" ;
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

#endif // __sun
