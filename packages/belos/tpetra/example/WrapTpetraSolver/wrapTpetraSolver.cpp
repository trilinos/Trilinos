#include "BelosSolverManager.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "KokkosBlas1_mult.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_ParameterList.hpp"
#include <cstdlib> // EXIT_SUCCESS
#include <iostream>

/// \file wrapTpetraSolver.cpp
/// \brief Example of how to wrap a "native" solver as a Belos solver
///
/// By a "native" solver, I mean a solver written for a particular
/// linear algebra implementation.  This corresponds to a particular
/// (ScalarType, MV, and OP) Belos template parameter combination.
/// "Wrap as a Belos solver" means to make available as a
/// Belos::SolverManager subclass.
///
/// This example includes a "stub" generic implementation of the
/// Belos::SolverManager subclass, as well as the actual non-generic
/// implementation specifically for our linear algebra of choice
/// (Tpetra in this case).  In order to make this example shorter,
/// I've chosen a linear algebra implementation for which Belos has
/// MultiVecTraits and OperatorTraits specializations already
/// implemented.

template<class SC, class LO, class GO, class NT>
class SolverInput {
private:
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType mag_type;

public:
  mag_type r_norm_orig;
  mag_type tol;
  int maxNumIters;
  bool needToScale;
};


/// \brief Result of a linear solve.
///
/// "A linear solve" refers either to a solve Ax=b with a single
/// right-hand side b, or to an aggregation of results of two such
/// solves with the same matrix, but different right-hand sides.
/// "Aggregation" here just means reporting a single group of metrics.
/// For example, for two solves, we report the max of the two
/// iteration counts, and the max of the two residual norms.
template<class SC, class LO, class GO, class NT>
class SolverOutput {
public:
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType mag_type;

private:
  typedef Teuchos::ScalarTraits<mag_type> STM;

public:
  //! Absolute residual norm.
  mag_type absResid;
  //! Relative residual norm (if applicable).
  mag_type relResid;
  //! Number of iterations executed.
  int numIters;
  //! Whether the solve converged.
  bool converged;

  /// \brief Default constructor.
  ///
  /// The default constructor creates output corresponding to "solving
  /// a linear system with zero right-hand sides."  This means that
  /// the solve succeeded trivially (converged == true), in zero
  /// iterations, with zero residual norm.
  SolverOutput () :
    absResid (STM::zero ()),
    relResid (STM::zero ()),
    numIters (0),
    converged (true)
  {}

  /// \brief Combine two solver outputs.
  ///
  /// "Combining" two solver outputs means aggregating the results of
  /// solving A x_1 = b_1 and A x_2 = b_2, that is, two solves with
  /// the same matrix, but different right-hand sides.  Combining is
  /// associative and commutative.
  void
  combine (const SolverOutput<SC, LO, GO, NT>& src)
  {
    // max of the residuals and iteration counts
    relResid = relResid > src.relResid ? relResid : src.relResid;
    absResid = absResid > src.absResid ? absResid : src.absResid;
    numIters = numIters > src.numIters ? numIters : src.numIters;
    // "converged" if all converged
    converged = converged && src.converged;
  }
};

template<class SC, class LO, class GO, class NT>
std::ostream&
operator<< (std::ostream& out,
            const SolverOutput<SC, LO, GO, NT>& so)
{
  using std::endl;

  out << "Solver output:" << endl
      << "  Absolute residual norm: " << so.absResid << endl
      << "  Relative residual norm: " << so.relResid << endl
      << "  Number of iterations: " << so.numIters << endl
      << "  Converged: " << (so.converged ? "true" : "false");
  return out;
}

/// \brief Tpetra implementation of CG
///
/// This CG implementation can solve linear systems with multiple
/// right-hand sides, but it solves them one right-hand side at a
/// time.  The reported convergence results in the case of multiple
/// right-hand sides is the max of the residuals and iteration counts,
/// and the AND of the "did it converge" Boolean values.
template<class SC = Tpetra::Operator<>::scalar_type,
         class LO = typename Tpetra::Operator<SC>::local_ordinal_type,
         class GO = typename Tpetra::Operator<SC, LO>::global_ordinal_type,
         class NT = typename Tpetra::Operator<SC, LO, GO>::node_type>
class CG {
public:
  typedef Tpetra::Operator<SC, LO, GO, NT> op_type;
  typedef Tpetra::MultiVector<SC, LO, GO, NT> mv_type;

private:
  typedef Tpetra::Vector<SC, LO, GO, NT> vec_type;
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType mag_type;
  typedef Teuchos::ScalarTraits<mag_type> STM;
  typedef typename NT::device_type device_type;

public:
  CG () :
    tol_ (STM::squareroot (STS::eps ())),
    maxNumIters_ (100),
    verbosity_ (0),
    needToScale_ (true),
    relResid_ (STM::zero ()),
    absResid_ (STM::zero ()),
    numIters_ (0),
    converged_ (false)
  {}

  CG (const Teuchos::RCP<const op_type>& A) :
    A_ (A),
    tol_ (STM::squareroot (STS::eps ())),
    maxNumIters_ (100),
    verbosity_ (0),
    needToScale_ (true),
    relResid_ (STM::zero ()),
    absResid_ (STM::zero ()),
    numIters_ (0),
    converged_ (false)
  {}

  //! Set the matrix A in the linear system to solve.
  void setMatrix (const Teuchos::RCP<const op_type>& A) {
    if (A_.getRawPtr () != A.getRawPtr ()) {
      A_ = A;
    }
  }

  Teuchos::RCP<const op_type> getMatrix () const {
    return A_;
  }

  /// \brief Fill \c params with all parameters this solver
  ///   understands, and either their current values, or their default
  ///   values.
  ///
  /// \param params [out] To be filled with this solver's parameters.
  /// \param defaultValues [in] Whether to use default values (true)
  ///   or current values (false).
  void
  getParameters (Teuchos::ParameterList& params,
                 const bool defaultValues) const
  {
    // Yes, the inner STS is supposed to be STS.  STS::eps() returns
    // mag_type.  It's SC's machine epsilon.
    const mag_type tol =
      defaultValues ? STM::squareroot (STS::eps ()) : tol_;
    const int maxNumIters = defaultValues ? 100 : maxNumIters_;
    const int verbosity = defaultValues ? 0 : verbosity_;
    const bool needToScale = defaultValues ? true : needToScale_;
    const std::string implResScal = needToScale ?
      "Norm of Preconditioned Initial Residual" :
      "None";

    params.set ("Convergence Tolerance", tol);
    params.set ("Implicit Residual Scaling", implResScal);
    params.set ("Maximum Iterations", maxNumIters);
    params.set ("Verbosity", verbosity);
  }

  /// \brief Set the solver's parameters.
  ///
  /// This solver takes a subset of the parameters that
  /// Belos::PseudoBlockCGSolMgr (Belos' CG implementation) takes, and
  /// ignores the rest.  If it takes a parameter but doesn't implement
  /// all that parameter's options, it throws an exception if it
  /// encounters an option that it does not implement.  The point is
  /// compatibility with Belos.
  void setParameters (Teuchos::ParameterList& params) {
    mag_type tol = tol_;
    bool needToScale = needToScale_;
    if (params.isParameter ("Convergence Tolerance")) {
      tol = params.get<mag_type> ("Convergence Tolerance");
      TEUCHOS_TEST_FOR_EXCEPTION
        (tol < STM::zero (), std::invalid_argument,
         "\"Convergence tolerance\" = " << tol << " < 0.");
    }
    if (params.isParameter ("Implicit Residual Scaling")) {
      const std::string implScal =
        params.get<std::string> ("Implicit Residual Scaling");
      if (implScal == "Norm of Initial Residual") {
        // FIXME (mfh 26 Oct 2016) Once we implement left
        // preconditioning, we'll have to keep separate preconditioned
        // and unpreconditioned absolute residuals.
        needToScale = true;
      }
      else if (implScal == "Norm of Preconditioned Initial Residual") {
        needToScale = true;
      }
      else if (implScal == "Norm of RHS") {
        // FIXME (mfh 26 Oct 2016) If we want to implement this, it
        // would make sense to combine that all-reduce with the
        // all-reduce for computing the initial residual norms.  We
        // could modify computeResiduals to have an option to do this.
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error,
           "\"Norm of RHS\" scaling option not implemented");
      }
      else if (implScal == "None") {
        needToScale = false;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument, "\"Implicit Residual Scaling\""
           " has an invalid value \"" << implScal << "\".");
      }
    }

    int maxNumIters = maxNumIters_;
    if (params.isParameter ("Maximum Iterations")) {
      maxNumIters = params.get<int> ("Maximum Iterations");
      TEUCHOS_TEST_FOR_EXCEPTION
        (maxNumIters < 0, std::invalid_argument,
         "\"Maximum Iterations\" = " << maxNumIters << " < 0.");
    }

    int verbosity = verbosity_;
    if (params.isType<int> ("Verbosity")) {
      verbosity = params.get<int> ("Verbosity");
    }
    else if (params.isType<bool> ("Verbosity")) {
      const bool verbBool = params.get<bool> ("Verbosity");
      verbosity = verbBool ? 1 : 0;
    }

    tol_ = tol;
    maxNumIters_ = maxNumIters;
    verbosity_ = verbosity;
    needToScale_ = needToScale;
  }

private:
  static void
  computeResiduals (Kokkos::DualView<mag_type*, device_type>& norms,
                    mv_type& R,
                    const op_type& A,
                    const mv_type& X,
                    const mv_type& B)
  {
    typedef typename device_type::memory_space dev_mem_space;

    const SC ONE = STS::one ();
    A.apply (X, R);
    R.update (ONE, B, -ONE); // R := B - A*X

    norms.template modify<dev_mem_space> ();
    Kokkos::View<mag_type*, device_type> norms_d =
      norms.template view<dev_mem_space> ();
    R.norm2 (norms_d);
  }

  static mag_type
  computeResidual (vec_type& R,
                   const op_type& A,
                   const vec_type& X,
                   const vec_type& B)
  {
    const SC ONE = STS::one ();
    A.apply (X, R);
    R.update (ONE, B, -ONE); // R := B - A*X
    const mag_type r_norm = R.norm2 ();
    return r_norm;
  }

public:
  //! Solve the linear system(s) AX=B.
  SolverOutput<SC, LO, GO, NT>
  solve (mv_type& X, const mv_type& B)
  {
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::oblackholestream;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;

    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.is_null (), std::runtime_error, "Matrix A is null.  Please call "
       "setMatrix() with a nonnull argument before calling solve().");
    TEUCHOS_TEST_FOR_EXCEPTION
      (X.getNumVectors () != B.getNumVectors (), std::runtime_error,
       "X.getNumVectors() = " << X.getNumVectors () <<
       " != B.getNumVectors() = " << B.getNumVectors () << ".");

    RCP<FancyOStream> outPtr;
    if (verbosity_) {
      const int myRank = A_->getDomainMap ()->getComm ()->getRank ();
      if (myRank == 0) {
        outPtr = getFancyOStream (rcpFromRef (std::cout));
      }
      else {
        outPtr = getFancyOStream (rcp (new oblackholestream ()));
      }
    }
    return solveImpl (outPtr.getRawPtr (), X, B);
  }

private:
  //! Solve the linear system(s) AX=B.
  SolverOutput<SC, LO, GO, NT>
  solveImpl (Teuchos::FancyOStream* outPtr, mv_type& X, const mv_type& B)
  {
    using Teuchos::RCP;
    using std::endl;

    const size_t numVecs = B.getNumVectors ();
    Kokkos::DualView<mag_type*, device_type> norms ("norms", numVecs);
    mv_type R (B.getMap (), numVecs);

    computeResiduals (norms, R, *A_, X, B);
    norms.template sync<Kokkos::HostSpace> ();
    auto norms_h = norms.template view<Kokkos::HostSpace> ();

    SolverInput<SC, LO, GO, NT> input;
    input.tol = tol_;
    input.maxNumIters = maxNumIters_;
    input.needToScale = needToScale_;
    SolverOutput<SC, LO, GO, NT> allOutput;

    for (size_t j = 0; j < numVecs; ++j) {
      if (outPtr != NULL) {
        *outPtr << "Solve for column " << (j+1) << " of " << numVecs << ":" << endl;
        outPtr->pushTab ();
      }
      RCP<vec_type> R_j = R.getVectorNonConst (j);
      RCP<vec_type> X_j = X.getVectorNonConst (j);
      input.r_norm_orig = norms_h(j);
      SolverOutput<SC, LO, GO, NT> curOutput;
      solveOneVec (outPtr, curOutput, *X_j, *R_j, *A_, input);
      allOutput.combine (curOutput);
      if (outPtr != NULL) {
        outPtr->popTab ();
      }
    }
    return allOutput;
  }

  static mag_type
  getConvergenceMetric (const mag_type r_norm_new,
                        const SolverInput<SC, LO, GO, NT>& input)
  {
    if (input.needToScale) {
      return input.r_norm_orig == STM::zero () ?
        r_norm_new :
        (r_norm_new / input.r_norm_orig);
    }
    else {
      return r_norm_new;
    }
  }

  static void
  solveOneVec (Teuchos::FancyOStream* outPtr,
               SolverOutput<SC, LO, GO, NT>& output,
               vec_type& X, // in/out
               vec_type& R, // in/out
               const op_type& A,
               const SolverInput<SC, LO, GO, NT>& input)
  {
    using std::endl;

    const SC ONE = STS::one ();

    vec_type P (R, Teuchos::Copy);
    vec_type AP (R.getMap ());
    mag_type r_norm_old = input.r_norm_orig; // R.norm2 ()

    if (r_norm_old == STM::zero () || r_norm_old <= input.tol) {
      if (outPtr != NULL) {
        *outPtr << "Initial guess' residual norm " << r_norm_old
                << " meets tolerance " << input.tol << endl;
      }
      output.absResid = r_norm_old;
      output.relResid = r_norm_old;
      output.numIters = 0;
      output.converged = true;
      return;
    }

    mag_type r_norm_new = STM::zero ();
    for (int iter = 0; iter < input.maxNumIters; ++iter) {
      if (outPtr != NULL) {
        *outPtr << "Iteration " << (iter+1) << " of " << input.maxNumIters << ":" << endl;
        outPtr->pushTab ();
        *outPtr << "r_norm_old: " << r_norm_old << endl;
      }

      A.apply (P, AP);
      const mag_type PAP = P.dot (AP);
      if (outPtr != NULL) {
        *outPtr << "PAP: " << PAP << endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION
        (PAP <= STM::zero (), std::runtime_error, "At iteration " << (iter+1)
         << " out of " << input.maxNumIters << ", P.dot(AP) = " << PAP <<
         " <= 0.  This usually means that the matrix A is not symmetric "
         "(Hermitian) positive definite.");

      const mag_type alpha = (r_norm_old * r_norm_old) / PAP;
      if (outPtr != NULL) {
        *outPtr << "alpha: " << alpha << endl;
      }
      X.update (static_cast<SC> (alpha), P, ONE);
      R.update (static_cast<SC> (-alpha), AP, ONE);

      r_norm_new = R.norm2 ();
      if (outPtr != NULL) {
        *outPtr << "r_norm_new: " << r_norm_new << endl;
      }
      const mag_type metric = getConvergenceMetric (r_norm_new, input);
      if (outPtr != NULL) {
        *outPtr << "metric: " << metric << endl;
      }
      if (metric <= input.tol) {
        output.absResid = r_norm_new;
        output.relResid = input.r_norm_orig == STM::zero () ?
          r_norm_new :
          (r_norm_new / input.r_norm_orig);
        output.numIters = iter + 1;
        output.converged = true;
        return;
      }
      else if (iter + 1 < input.maxNumIters) { // not last iteration
        const mag_type beta = (r_norm_new * r_norm_new) /
          (r_norm_old * r_norm_old);
        P.update (ONE, R, static_cast<SC> (beta));
        r_norm_old = r_norm_new;
      }

      if (outPtr != NULL) {
        outPtr->popTab ();
      }
    }

    // Reached max iteration count without converging
    output.absResid = r_norm_new;
    output.relResid = input.r_norm_orig == STM::zero () ?
      r_norm_new :
      (r_norm_new / input.r_norm_orig);
    output.numIters = input.maxNumIters;
    output.converged = false;
  }

private:
  Teuchos::RCP<const op_type> A_;

  mag_type tol_;
  int maxNumIters_;
  int verbosity_;
  // FIXME (mfh 26 Oct 2016) Once we implement left preconditioning,
  // we'll have to keep separate preconditioned and unpreconditioned
  // absolute residuals.
  bool needToScale_;

  mag_type relResid_;
  mag_type absResid_;
  int numIters_;
  bool converged_;
};


template<class ScalarType, class MV, class OP>
class CgWrapper :
  public Belos::SolverManager<ScalarType, MV, OP>
{
public:
  virtual ~CgWrapper () {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  const Belos::LinearProblem<ScalarType,MV,OP>& getProblem () const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters () const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters () const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  int getNumIters () const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  bool isLOADetected () const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  void setProblem (const Teuchos::RCP<Belos::LinearProblem<ScalarType,MV,OP> >& problem) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  void reset (const Belos::ResetType type) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  Belos::ReturnType solve () {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  virtual Teuchos::RCP<Belos::SolverManager<ScalarType, MV, OP> > clone () const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "clone() not implemented");
  }
};

template<class SC, class LO, class GO, class NT>
class CgWrapper<SC,
                Tpetra::MultiVector<SC, LO, GO, NT>,
                Tpetra::Operator<SC, LO, GO, NT> > :
  public Belos::SolverManager<SC,
                              Tpetra::MultiVector<SC, LO, GO, NT>,
                              Tpetra::Operator<SC, LO, GO, NT> >
{
public:
  typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
  typedef Tpetra::Operator<SC, LO, GO, NT> OP;
  typedef Belos::LinearProblem<SC, MV, OP> belos_problem_type;

  virtual ~CgWrapper () {}

  const Belos::LinearProblem<SC,MV,OP>& getProblem () const {
    TEUCHOS_TEST_FOR_EXCEPTION
      (problem_.is_null (), std::runtime_error, "The linear problem has not "
       "yet been set.  Please call setProblem with a nonnull argument before "
       "calling this method.");
    return *problem_;
  }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters () const {
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    RCP<ParameterList> params (new ParameterList ("CG"));
    const bool defaultValues = true;
    solver_.getParameters (*params, defaultValues);
    return params;
  }

  Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters () const {
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    RCP<ParameterList> params (new ParameterList ("CG"));
    const bool defaultValues = false;
    solver_.getParameters (*params, defaultValues);
    return params;
  }

  int getNumIters () const {
    return lastSolverOutput_.numIters;
  }

  bool isLOADetected () const {
    return false; // this solver doesn't attempt to detect loss of accuracy
  }

  void
  setProblem (const Teuchos::RCP<belos_problem_type>& problem)
  {
    if (problem.is_null ()) {
      solver_.setMatrix (Teuchos::null);
    }
    else if (solver_.getMatrix ().getRawPtr () !=
             problem->getOperator ().getRawPtr ()) {
      // setMatrix resets state, so only call if necessary.
      solver_.setMatrix (problem->getOperator ());
    }
    problem_ = problem;
  }

  void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params) {
    if (! params.is_null ()) {
      solver_.setParameters (*params);
    }
  }

  void reset (const Belos::ResetType /* type */ ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  Belos::ReturnType solve () {
    using Teuchos::RCP;
    TEUCHOS_TEST_FOR_EXCEPTION
      (problem_.is_null (), std::runtime_error, "The linear problem has not "
       "yet been set.  Please call setProblem with a nonnull argument before "
       "calling this method.");
    RCP<const MV> B = problem_->getRHS ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (B.is_null (), std::runtime_error, "The linear problem's right-hand "
       "side(s) B has/have not yet been set.  Please call setProblem with "
       "a nonnull argument before calling this method.");
    RCP<MV> X = problem_->getLHS ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (X.is_null (), std::runtime_error, "The linear problem's left-hand "
       "side(s) X has/have not yet been set.  Please call setProblem with "
       "a nonnull argument before calling this method.");
    SolverOutput<SC, LO, GO, NT> result = solver_.solve (*X, *B);
    lastSolverOutput_ = result;
    return result.converged ? Belos::Converged : Belos::Unconverged;
  }

  virtual Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > clone () const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "clone() not implemented");
  }

private:
  typedef SolverOutput<SC, LO, GO, NT> output_type;

  CG<SC, LO, GO, NT> solver_;
  //! Output of the last solve.
  ///
  /// This does not include the solution (multi)vector, just things
  /// like the residual nor, the iteration count, and whether the
  /// solve converged.  See SolverOutput documentation above for
  /// details.
  output_type lastSolverOutput_;
  Teuchos::RCP<belos_problem_type> problem_;
};


template<class SC = Tpetra::Operator<>::scalar_type,
         class LO = typename Tpetra::Operator<SC>::local_ordinal_type,
         class GO = typename Tpetra::Operator<SC, LO>::global_ordinal_type,
         class NT = typename Tpetra::Operator<SC, LO, GO>::node_type>
class TestDiagonalOperator : public Tpetra::Operator<SC, LO, GO, NT> {
private:
  typedef typename NT::device_type device_type;

public:
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType mag_type;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  TestDiagonalOperator (const Teuchos::RCP<const map_type>& map,
                        const mag_type minSingularValue) :
    map_ (map)
  {
    typedef Kokkos::View<mag_type*, device_type> dev_view_type;
    typedef typename dev_view_type::HostMirror host_view_type;
    typedef Teuchos::ScalarTraits<SC> STS;

    const LO lclNumRows = map_.is_null () ?
      static_cast<LO> (0) :
      static_cast<LO> (map_->getNodeNumElements ());
    dev_view_type diag_d ("diag", lclNumRows);
    host_view_type diag_h = Kokkos::create_mirror_view (diag_d);

    if (lclNumRows != 0) {
      const SC ONE = STS::one ();
      const SC exponent = ONE / static_cast<mag_type> (lclNumRows - 1);
      const SC scalingFactor = ::pow (minSingularValue, exponent);
      diag_h(0) = ONE;
      for (LO lclRow = 1; lclRow < lclNumRows; ++lclRow) {
        diag_h(lclRow) = diag_h(lclRow-1) * scalingFactor;
      }
    }
    Kokkos::deep_copy (diag_d, diag_h);
    diag_ = diag_d;
  }

  Teuchos::RCP<const map_type> getDomainMap () const {
    return map_;
  }

  Teuchos::RCP<const map_type> getRangeMap () const {
    return map_;
  }

  void
  apply (const Tpetra::MultiVector<SC, LO, GO, NT>& X,
         Tpetra::MultiVector<SC, LO, GO, NT>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         SC alpha = Teuchos::ScalarTraits<SC>::one (),
         SC beta = Teuchos::ScalarTraits<SC>::zero ()) const
  {
    using Teuchos::RCP;
    typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
    typedef typename MV::impl_scalar_type ISC;

    TEUCHOS_TEST_FOR_EXCEPTION
      (mode == Teuchos::CONJ_TRANS && Teuchos::ScalarTraits<SC>::isComplex,
       std::logic_error, "Conjugate transpose case not implemented.");

    const size_t numVecs = X.getNumVectors ();
    for (size_t j = 0; j < numVecs; ++j) {
      RCP<const MV> X_j = X.getVector (j);
      RCP<MV> Y_j = Y.getVectorNonConst (j);

      auto X_j_lcl_2d = X_j->template getLocalView<device_type> ();
      auto X_j_lcl = Kokkos::subview (X_j_lcl_2d, Kokkos::ALL (), 0);
      auto Y_j_lcl_2d = Y_j->template getLocalView<device_type> ();
      auto Y_j_lcl = Kokkos::subview (Y_j_lcl_2d, Kokkos::ALL (), 0);

      KokkosBlas::mult (static_cast<ISC> (beta), Y_j_lcl,
                        static_cast<ISC> (alpha), X_j_lcl, diag_);
    }
  }

private:
  Teuchos::RCP<const map_type> map_;
  Kokkos::View<const mag_type*, device_type> diag_;
};


int
main (int argc, char* argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::MultiVector<> mv_type;
  typedef Tpetra::Operator<> op_type;
  typedef mv_type::scalar_type SC;
  typedef map_type::global_ordinal_type GO;
  typedef mv_type::mag_type mag_type;
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef Teuchos::ScalarTraits<mag_type> STM;
  const SC ZERO = STS::zero ();
  const SC ONE = STS::one ();

  Tpetra::initialize (&argc, &argv);
  auto comm = Tpetra::getDefaultComm ();

  const int myRank = comm->getRank ();

  const GO gblNumRows = 10000;
  const GO indexBase = 0;
  RCP<const map_type> map (new map_type (gblNumRows, indexBase, comm));
  const mag_type minSingVal = 0.1;
  RCP<TestDiagonalOperator<> > A =
    rcp (new TestDiagonalOperator<> (map, minSingVal));

  mv_type X (A->getDomainMap (), 1);
  mv_type B (A->getRangeMap (), 1);

  B.randomize ();

  CG<> solver (A);
  auto out = solver.solve (X, B);
  if (myRank == 0) {
    cout << out << endl;
  }

  {
    // Check the residual.
    mv_type X_copy (X, Teuchos::Copy);
    mv_type R (B.getMap (), B.getNumVectors ());
    A->apply (X_copy, R);
    R.update (ONE, B, -ONE);
    Teuchos::Array<mag_type> R_norms (R.getNumVectors ());
    R.norm2 (R_norms ());
    Teuchos::Array<mag_type> B_norms (B.getNumVectors ());
    B.norm2 (B_norms ());
    if (myRank == 0) {
      for (size_t j = 0; j < R.getNumVectors (); ++j) {
        const mag_type relResNorm = (B_norms[j] == STM::zero ()) ?
          R_norms[j] :
          R_norms[j] / B_norms[j];
        cout << "Column " << (j+1) << " of " << R.getNumVectors ()
             << ": Absolute residual norm: " << R_norms[j]
             << ", Relative residual norm: " << relResNorm
             << endl;
      }
      cout << endl;
    }
  }

  // Get ready for next solve by resetting initial guess to zero.
  X.putScalar (ZERO);

  typedef Belos::LinearProblem<SC, mv_type, op_type> belos_problem_type;
  {
    RCP<Belos::LinearProblem<SC, mv_type, op_type> > lp =
      rcp (new belos_problem_type (A, rcpFromRef (X), rcpFromRef (B)));
    // Our CG implementation bypasses this, but we call it anyway, just
    // for interface consistency with Belos' other solvers.
    lp->setProblem ();

    CgWrapper<SC, mv_type, op_type> solverWrapper;
    solverWrapper.setProblem (lp);
    const Belos::ReturnType belosResult = solverWrapper.solve ();
    if (myRank == 0) {
      cout << "Belos solver wrapper result: "
           << (belosResult == Belos::Converged ? "Converged" : "Unconverged")
           << endl
           << "Number of iterations: " << solverWrapper.getNumIters ()
           << endl;
    }
  }

  {
    // Check the residual.
    mv_type X_copy (X, Teuchos::Copy);
    mv_type R (B.getMap (), B.getNumVectors ());
    A->apply (X_copy, R);
    R.update (ONE, B, -ONE);
    Teuchos::Array<mag_type> R_norms (R.getNumVectors ());
    R.norm2 (R_norms ());
    Teuchos::Array<mag_type> B_norms (B.getNumVectors ());
    B.norm2 (B_norms ());
    if (myRank == 0) {
      for (size_t j = 0; j < R.getNumVectors (); ++j) {
        const mag_type relResNorm = (B_norms[j] == STM::zero ()) ?
          R_norms[j] :
          R_norms[j] / B_norms[j];
        cout << "Column " << (j+1) << " of " << R.getNumVectors ()
             << ": Absolute residual norm: " << R_norms[j]
             << ", Relative residual norm: " << relResNorm
             << endl;
      }
      cout << endl;
    }
  }

  // Prepare for next linear solve by resetting initial guess to zero.
  X.putScalar (ZERO);

  Belos::PseudoBlockCGSolMgr<SC, mv_type, op_type> belosSolver;
  {
    RCP<Belos::LinearProblem<SC, mv_type, op_type> > lp =
      rcp (new belos_problem_type (A, rcpFromRef (X), rcpFromRef (B)));
    // Our CG implementation bypasses this, but we call it anyway, just
    // for interface consistency with Belos' other solvers.
    lp->setProblem ();

    belosSolver.setProblem (lp);
    const Belos::ReturnType belosResult = belosSolver.solve ();
    if (myRank == 0) {
      cout << "Belos solver (PseudoBlockCGSolMgr) result: "
           << (belosResult == Belos::Converged ? "Converged" : "Unconverged")
           << endl
           << "Number of iterations: " << belosSolver.getNumIters ()
           << endl;
    }

    // Check the residual.
    mv_type X_copy (X, Teuchos::Copy);
    mv_type R (B.getMap (), B.getNumVectors ());
    A->apply (X_copy, R);
    R.update (ONE, B, -ONE);
    Teuchos::Array<mag_type> R_norms (R.getNumVectors ());
    R.norm2 (R_norms ());
    Teuchos::Array<mag_type> B_norms (B.getNumVectors ());
    B.norm2 (B_norms ());
    if (myRank == 0) {
      for (size_t j = 0; j < R.getNumVectors (); ++j) {
        const mag_type relResNorm = (B_norms[j] == STM::zero ()) ?
          R_norms[j] :
          R_norms[j] / B_norms[j];
        cout << "Column " << (j+1) << " of " << R.getNumVectors ()
             << ": Absolute residual norm: " << R_norms[j]
             << ", Relative residual norm: " << relResNorm
             << endl;
      }
    }
  }

  Tpetra::finalize ();
  return EXIT_SUCCESS;
}
