
/// The type of the error tolerance
enum EErrTolType { REL_RESIDUAL_NORM, REL_SOLUTION_ERR_NORM };

/// Exception type throw on blow up
class CatastrophicSolveFalire;

/// Solve tolerance
template <class Scalar>
struct SolveTolerance {
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  /** \brief . */
  static const ScalarMag  DEFAULT_SOLVE_TOLERANCE = -1;
  /** \brief . */
  ScalarMag      requestedTol;  ///< DEFAULT_SOLVE_TOLERANCE means use default
  EErrTolType    errTolType;
  /** \brief . */
  SolveTolerance()
    : requestedTol(DEFAULT_SOLVE_TOLERANCE), errTolType(REL_RESIDUAL_ERR) {}
  SolveTolerance(ScalarMag _requestedTol, EErrorTolType _errTolType)
    : requestedTol(_requestedTol), errTolType(_errTolType) {}
};

template <class Scalar>
struct BlockSolveTolerance {
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  /** \brief . */
  ...
  /** \brief . */
  SolveTolerance<Scalar>   solveTolerance;
  Range1D                  blockRange;
};

/// The type of the solution
enum ESolveReturnStatus { SOLVE_STATUS_CONVERGED, SOLVE_STATUS_UNCONVERGED };

/// Solve return
template <class Scalar>
struct SolveReturn {
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  /** \brief . */
  ScalarMag             achievedTol;         // < -1 means I don't know
  ESolveReturnStatus    solveReturnStatus;
  int                   numIterations;
};

/// A linear operator with a solve function
template <class Scalar>
class LinearOpWithSolve : virtual public LinearOp<Scalar> {
public:

  /** \brief Solve (or try to solve) a block system with one set of tolerances.
   *
   * With throw CatastrophicSolveFalire on complete failure!
   */
  virtual SolveReturn<Scalar> solve(
    const MultiVectorBase<Scalar>   &B
    ,MultiVectorBase<Scalar>        *X                                          // [inout]
    ,const SolveTolerance<Scalar>   &solveTolerance = SolveTolerance<Scalar>()  // Default is to use built-in tol
    ) const = 0;

  /** \brief Solve (or try to solve) a block system with different targeted tolerances.
   *
   * With throw CatastrophicSolveFalire on complete failure!
   */
  virtual SolveReturn<Scalar> solve(
    const MultiVectorBase<Scalar>        &B
    ,MultiVectorBase<Scalar>             *X
    ,const int                           numBlocks
    ,const BlockSolveTolerance<Scalar>   blockSolveTolerances[]    // ==NULL means use built-in tol
    ) const = 0;

};

//
// Use cases
//

template <class Scalar>
void solveUsingDefaultTol(
  const LinearOpWithSolve<Scalar>    &A
  ,const MultiVectorBase<Scalar>     &B
  ,MultiVectorBase<Scalar>           *X
  )
{
  SolveReturn<Scalar>
    solveReturn = A.solve(B,X);
  // Check return
  ...
}

template <class Scalar>
void solveUsingResidualTol(
  const LinearOpWithSolve<Scalar>                      &A
  ,const MultiVectorBase<Scalar>                       &B
  ,const Teuchos::ScalarTraits<Scalar>::magnitudeType  tol
  ,MultiVectorBase<Scalar>                             *X
  )
{
  SolveReturn<Scalar>
    solveReturn = A.solve(B,X,SolveTolerance<Scalar>(tol,REL_RESIDUAL_NORM));
  // Check return
  ...
}

template <class Scalar>
void solveUsingSolutionTol(
  const LinearOpWithSolve<Scalar>                      &A
  ,const MultiVectorBase<Scalar>                       &B
  ,const Teuchos::ScalarTraits<Scalar>::magnitudeType  tol
  ,MultiVectorBase<Scalar>                             *X
  )
{
  SolveReturn<Scalar>
    solveReturn = A.solve(B,X,SolveTolerance<Scalar>(tol,REL_SOLUTION_ERR_NORM));
  // Check return
  ...
}
