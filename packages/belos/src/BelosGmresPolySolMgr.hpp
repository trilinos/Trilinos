// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//

#ifndef BELOS_GMRES_POLY_SOLMGR_HPP
#define BELOS_GMRES_POLY_SOLMGR_HPP

/// \file BelosGmresPolySolMgr.hpp
/// \brief Declaration and definition of Belos::GmresPolySolMgr
///   (hybrid block GMRES linear solver).

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"
#include "BelosGmresPolyOp.hpp"
#include "BelosSolverFactory_Generic.hpp"
#include "BelosOrthoManagerFactory.hpp"
#include "Teuchos_as.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

/** \example epetra/example/BlockGmres/BlockGmresPolyEpetraExFile.cpp
    This is an example of how to use the Belos::GmresPolySolMgr with Epetra.
*/

namespace Belos {

//! @name GmresPolySolMgr Exceptions
//@{

/** \brief GmresPolySolMgrLinearProblemFailure is thrown when the linear problem is
 * not setup (i.e. setProblem() was not called) when solve() is called.
 *
 * This std::exception is thrown from the GmresPolySolMgr::solve() method.
 *
 */
class GmresPolySolMgrLinearProblemFailure : public BelosError {public:
  GmresPolySolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

/** \brief GmresPolySolMgrPolynomialFailure is thrown when their is a problem generating
 * the GMRES polynomial for this linear problem.
 *
 * This std::exception is thrown from the GmresPolySolMgr::solve() method.
 *
 */
class GmresPolySolMgrPolynomialFailure : public BelosError {public:
  GmresPolySolMgrPolynomialFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

/// \class Belos::GmresPolySolMgr
/// \brief Hybrid block GMRES iterative linear solver.
/// \author Heidi Thornquist and Jennifer Loe
/// \ingroup belos_solver_framework
///
/// The GMRES polynomial solver manager can perform two types of linear
/// solves. First the solver runs block GMRES, storing the resulting 
/// coefficients (or roots), which can be used to form a matrix polynomial.
/// It then can reuse this polynomial as either a surrogate operator, or
/// as a preconditioner for an outer solver.
/// By applying the GMRES polynomial as an operator or preconditioner, one 
/// avoids the cost of the inner products and norms in GMRES, thus reducing
/// communication costs.  
//
/// The GMRES polynomial can be created in conjunction with any standard preconditioner.
/// Simply pass the preconditioner to the LinearProblem before calling the GmresPolySolMgr
/// and your preconditioner will be combined with the polynomial automatically.
///
/// Here is a list of all the parameters that this solver accepts:
///   - "Polynomial Type" (\c std::string): The desired polynomial type: 
///      Roots, Arnoldi, or Gmres.  Default: "Roots"
///   - "Polynomial Tolerance" (\c MagnitudeType): The level that
///     residual norms must reach to decide convergence. Default:
///     1e-8.
///   - "Maximum Degree" (\c int): Requested maximum degree for the polynomial. 
///      The preconditioned problem Ap(A) will have this degree, while the polynomial
///       p(A) will have degree deg-1.  Default: 25
///   - "Random RHS" (\c bool): to generate the polynomial using a random vector. Default: true
///   - "Add Roots" (\c bool): to add roots to the polynomial as needed for stability. Default: true
///   - "Damp Poly" (\c bool): to damp polynomial. Default: false
///   - "Orthogonalization" (\c std::string): The desired
///     orthogonalization method to create polynomial.  Default: "ICGS".
///   - "Verbosity" (Belos::MsgType): A sum of Belos::MsgType values
///     specifying what kinds of messages to print.  Default:
///     Belos::Errors.
///   - "Outer Solver" (\c std::string): Name of the outer solver in Belos solver factory.
///   - "Outer Solver Params" (\c Teuchos::parameterList): List of parameters for the outer solver
///   - "Timer Label" (\c std::string): Label for timers with polynomial solve.
///
/// This solver manager provides three different implementations of the same polynomial preconditioner.
/// The polynomial is the minimum residual polynomial from GMRES.
/// The "Roots" version is default.  It is the only implementation which provides the option of added 
/// roots for stability.  These added roots can allow for high-degree polynomials.  
/// The "Arnoldi" version typically gives similar results to the "Roots" version
/// but is slightly more expensive to apply. Both of these polynomials can be "damped", which is 
/// sometimes useful for indefinite or other ill-conditioned problems. 
/// The "Gmres" version is based on a power-basis implementation
/// and is only stable for well-conditioned problems and low-degree polynomials. 
//
/// For more information on the implementation and formulas, see the following references:
/// "Roots" version: https://arxiv.org/abs/1806.08020 (Includes explanation of root-adding and damping.)
/// "Arnoldi" version: https://scholarship.rice.edu/handle/1911/17630
/// "Gmres" version: https://epubs.siam.org/doi/pdf/10.1137/140968276
///
///
/// Like all Belos solvers, parameters have relative or "delta"
/// semantics.  This means the following:
///   - Any parameter that was <i>never</i> set has its default value
///   - Any parameter not explicitly set in the input ParameterList
///     retains its current value

template<class ScalarType, class MV, class OP>
class GmresPolySolMgr : public SolverManager<ScalarType,MV,OP> {
private:

  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

  typedef Teuchos::ScalarTraits<MagnitudeType> MTS;
  typedef Belos::GmresPolyOp<ScalarType,MV,OP> gmres_poly_t;
  typedef Belos::GmresPolyMv<ScalarType,MV>    gmres_poly_mv_t;

public:

  //! @name Constructors/Destructor
  //@{

  /*! \brief Empty constructor for GmresPolySolMgr.
   * This constructor takes no arguments and sets the default values for the solver.
   * The linear problem must be passed in using setProblem() before solve() is called on this object.
   * The solver values can be changed using setParameters().
   */
  GmresPolySolMgr();

  /*! \brief Basic constructor for GmresPolySolMgr.
   *
   * This constructor accepts the LinearProblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - "Polynomial Type" -a \c std::string specifying the type of polynomial: Roots, Arnoldi, or Gmres.  Default: "Roots"
   *   - "Maximum Degree" - a \c int specifying the maximum degree of the polynomial. Default: 25
   *   - "Random RHS" - a \c bool indicates whether to generate polynomial using a random vector.  Default: true
   *   - "Add Roots" - a \c bool to add roots to the polynomial as needed for stability. Default: true
   *   - "Damp Poly" - a \c bool to damp polynomial. Default: false
   *   - "Orthogonalization" - a \c std::string specifying the desired orthogonalization to create the 
   *                            polynomial:  DGKS, ICGS, and IMGS. Default: "ICGS"
   *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Belos::Errors
   *   - "Polynomial Tolerance" - a \c MagnitudeType specifying the polynomial tolerance (sometimes) used to 
   *                            generate polynomial. Default: 1e-8
   *   - "Outer Solver" -a \c std::string specifying name of outer solver in Belos solver factory.  Default: ""
   *   - "Outer Solver Params" -a \c Teuchos::ParameterList giving parameters for the outer solver.
   *   - "Timer Label" -a \c std::string specifying the label on polynomial solve timers.
   */
  GmresPolySolMgr( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
    const Teuchos::RCP<Teuchos::ParameterList> &pl );

  //! Destructor.
  virtual ~GmresPolySolMgr() {};

  //! clone for Inverted Injection (DII)
  Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const override {
    return Teuchos::rcp(new GmresPolySolMgr<ScalarType,MV,OP>);
  }
  //@}

  //! @name Accessor methods
  //@{

  /*! \brief Get current linear problem being solved for in this object.
   */
  const LinearProblem<ScalarType,MV,OP>& getProblem() const override {
    return *problem_;
  }

  /*! \brief Get a parameter list containing the valid parameters for this object.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /*! \brief Get a parameter list containing the current parameters for this object.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const override { return params_; }

  /*! \brief Tolerance achieved by the last \c solve() invocation.
    
     This is the maximum over all right-hand sides' achieved
     convergence tolerances, and is set whether or not the solve
     actually managed to achieve the desired convergence tolerance.
    
     \warning This solver manager may be use as either a polynomial
       preconditioned iterative method or a polynomial preconditioner.
       In the later case, where a static polynomial is being applied
       through each call to solve(), there is not an outer solve that 
       can provide the achieved tolerance.
     \warning This result may not be meaningful if there was a loss
       of accuracy during the outer solve.  You should first call \c
       isLOADetected() to check for a loss of accuracy during the
       last solve.
  */
    MagnitudeType achievedTol() const override {
      return achievedTol_;
    }

  /*! \brief Return the timers for this object.
   *
   * The timers are ordered as follows:
   *   - time spent in solve() routine
   */
  Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
    return Teuchos::tuple(timerPoly_);
  }

  //! Get the iteration count for the most recent call to \c solve().
  int getNumIters() const override {
    return numIters_;
  }

  /*! \brief Return whether a loss of accuracy was detected by this solver during the most current solve.
      \note This flag will be reset the next time solve() is called.
   */
  bool isLOADetected() const override { return loaDetected_; }

  //@}

  //! @name Set methods
  //@{

  //! Set the linear problem that needs to be solved.
  void setProblem( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem ) override { problem_ = problem; }

  //! Set the parameters the solver manager should use to solve the linear problem.
  void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params ) override;

  //@}
  //! @name Reset methods
  //@{

  /// \brief Reset the solver.
  ///
  /// \param type [in] How to reset the solver.
  ///
  /// If type includes Belos::Problem, then reset the solver's state.
  /// This clears out the stored coefficients, so that the next call
  /// to solve() actually computes a full block GMRES solve, instead
  /// of just reusing the coefficients from the first solve.
  void reset( const ResetType type ) override {
    if ((type & Belos::Problem) && ! problem_.is_null ()) {
      problem_->setProblem ();
      poly_Op_ = Teuchos::null;
      poly_dim_ = 0;  // Rebuild the GMRES polynomial
    }
  }

  //@}
  //! @name Solver application methods
  //@{

  /*! \brief This method performs possibly repeated calls to the underlying linear solver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to
   * quit.
   *
   * This method calls BlockGmresIter::iterate(), which will return either because a specially constructed status test evaluates to
   * ::Passed or an std::exception is thrown.
   *
   * A return from BlockGmresIter::iterate() signifies one of the following scenarios:
   *    - the maximum number of restarts has been exceeded. In this scenario, the current solutions to the linear system
   *      will be placed in the linear problem and return ::Unconverged.
   *    - global convergence has been met. In this case, the current solutions to the linear system will be placed in the linear
   *      problem and the solver manager will return ::Converged
   *
   * \returns ::ReturnType specifying:
   *     - ::Converged: the linear problem was solved to the specification required by the solver manager.
   *     - ::Unconverged: the linear problem was not solved to the specification desired by the solver manager.
   */
  ReturnType solve() override;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief Method to return description of the hybrid block GMRES solver manager */
  std::string description() const override;

  //@}

private:

  // Linear problem.
  Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;

  // Output manager.
  Teuchos::RCP<std::ostream> outputStream_;

  // Current parameter list.
  Teuchos::RCP<Teuchos::ParameterList> params_;
  Teuchos::RCP<Teuchos::ParameterList> outerParams_;

  // Default solver values.
  static constexpr int maxDegree_default_ = 25;
  static constexpr int verbosity_default_ = Belos::Errors;
  static constexpr const char * label_default_ = "Belos";
  static constexpr const char * outerSolverType_default_ = "";
  static constexpr const char * polyType_default_ = "Arnoldi";
  static constexpr const char * orthoType_default_ = "ICGS";
  static constexpr bool addRoots_default_ = true;
  static constexpr bool dampPoly_default_ = false;
  static constexpr bool randomRHS_default_ = true; 

  // Current solver values.
  MagnitudeType polyTol_, achievedTol_;
  int maxDegree_, numIters_;
  int verbosity_;
  bool hasOuterSolver_;
  bool randomRHS_;
  bool damp_;
  bool addRoots_;
  std::string polyType_;
  std::string outerSolverType_;
  std::string orthoType_;

  // Polynomial storage
  int poly_dim_;
  Teuchos::RCP<gmres_poly_t> poly_Op_;

  // Timers.
  std::string label_;
  Teuchos::RCP<Teuchos::Time> timerPoly_;

  // Internal state variables.
  bool isSet_;
  bool loaDetected_;

  //! Cached default (valid) parameters.
  mutable Teuchos::RCP<const Teuchos::ParameterList> validPL_;
};


template<class ScalarType, class MV, class OP>
GmresPolySolMgr<ScalarType,MV,OP>::GmresPolySolMgr () :
  outputStream_ (Teuchos::rcpFromRef(std::cout)),
  polyTol_ (DefaultSolverParameters::polyTol),
  achievedTol_(MTS::zero()),
  maxDegree_ (maxDegree_default_),
  numIters_ (0),
  verbosity_ (verbosity_default_),
  hasOuterSolver_ (false),
  randomRHS_ (randomRHS_default_),
  damp_ (dampPoly_default_),
  addRoots_ (addRoots_default_),
  polyType_ (polyType_default_),
  outerSolverType_ (outerSolverType_default_),
  orthoType_ (orthoType_default_),
  poly_dim_ (0),
  label_ (label_default_),
  isSet_ (false),
  loaDetected_ (false)
{}


template<class ScalarType, class MV, class OP>
GmresPolySolMgr<ScalarType,MV,OP>::
GmresPolySolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                 const Teuchos::RCP<Teuchos::ParameterList> &pl) :
  problem_ (problem),
  outputStream_ (Teuchos::rcpFromRef(std::cout)),
  polyTol_ (DefaultSolverParameters::polyTol),
  maxDegree_ (maxDegree_default_),
  numIters_ (0),
  verbosity_ (verbosity_default_),
  hasOuterSolver_ (false),
  randomRHS_ (randomRHS_default_),
  damp_ (dampPoly_default_),
  addRoots_ (addRoots_default_),
  polyType_ (polyType_default_),
  outerSolverType_ (outerSolverType_default_),
  orthoType_ (orthoType_default_),
  poly_dim_ (0),
  label_ (label_default_),
  isSet_ (false),
  loaDetected_ (false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    problem_.is_null (), std::invalid_argument,
    "Belos::GmresPolySolMgr: The given linear problem is null.  "
    "Please call this constructor with a nonnull LinearProblem argument, "
    "or call the constructor that does not take a LinearProblem.");

  // If the input parameter list is null, then the parameters take
  // default values.
  if (! pl.is_null ()) {
    setParameters (pl);
  }
}


template<class ScalarType, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList>
GmresPolySolMgr<ScalarType,MV,OP>::getValidParameters() const
{
  if (validPL_.is_null ()) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList ();

    // The static_cast is to resolve an issue with older clang versions which
    // would cause the constexpr to link fail. With c++17 the problem is resolved.
    pl->set("Polynomial Type", static_cast<const char *>(polyType_default_),
      "The type of GMRES polynomial that is used as a preconditioner: Roots, Arnoldi, or Gmres.");
    pl->set("Polynomial Tolerance", static_cast<MagnitudeType>(DefaultSolverParameters::polyTol),
      "The relative residual tolerance that used to construct the GMRES polynomial.");
    pl->set("Maximum Degree", static_cast<int>(maxDegree_default_),
      "The maximum degree allowed for any GMRES polynomial.");
    pl->set("Outer Solver", static_cast<const char *>(outerSolverType_default_),
      "The outer solver that this polynomial is used to precondition.");
    pl->set("Outer Solver Params", Teuchos::ParameterList(),
      "Parameter list for the outer solver.");
    pl->set("Verbosity", static_cast<int>(verbosity_default_),
      "What type(s) of solver information should be outputted\n"
      "to the output stream.");
    pl->set("Output Stream", Teuchos::rcpFromRef(std::cout),
      "A reference-counted pointer to the output stream where all\n"
      "solver output is sent.");
    pl->set("Timer Label", static_cast<const char *>(label_default_),
      "The string to use as a prefix for the timer labels.");
    pl->set("Orthogonalization", static_cast<const char *>(orthoType_default_),
      "The type of orthogonalization to use to generate polynomial: DGKS, ICGS, or IMGS.");
    pl->set("Random RHS", static_cast<bool>(randomRHS_default_),
      "Add roots to polynomial for stability.");
    pl->set("Add Roots", static_cast<bool>(addRoots_default_),
      "Add roots to polynomial for stability.");
    pl->set("Damp Poly", static_cast<bool>(dampPoly_default_),
      "Damp polynomial for ill-conditioned problems.");
    validPL_ = pl;
  }
  return validPL_;
}


template<class ScalarType, class MV, class OP>
void GmresPolySolMgr<ScalarType,MV,OP>::
setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  // Create the internal parameter list if ones doesn't already exist.
  if (params_.is_null ()) {
    params_ = Teuchos::parameterList (*getValidParameters ());
  }
  else {
    params->validateParameters (*getValidParameters (),0);
  }

  // Check which Gmres polynomial to use
  if (params->isParameter("Polynomial Type")) {
    polyType_ = params->get("Polynomial Type", polyType_default_);
  }

  // Update the outer solver in our list.
  params_->set("Polynomial Type", polyType_);

  // Check if there is an outer solver for this Gmres Polynomial
  if (params->isParameter("Outer Solver")) {
    outerSolverType_ = params->get("Outer Solver", outerSolverType_default_);
  }

  // Update the outer solver in our list.
  params_->set("Outer Solver", outerSolverType_);

  // Check if there is a parameter list for the outer solver
  if (params->isSublist("Outer Solver Params")) {
    outerParams_ = Teuchos::parameterList( params->get<Teuchos::ParameterList>("Outer Solver Params") );
  }   

  // Check for maximum polynomial degree
  if (params->isParameter("Maximum Degree")) {
    maxDegree_ = params->get("Maximum Degree",maxDegree_default_);
  }

  // Update parameter in our list.
  params_->set("Maximum Degree", maxDegree_);

  // Check to see if the timer label changed.
  if (params->isParameter("Timer Label")) {
    std::string tempLabel = params->get("Timer Label", label_default_);

    // Update parameter in our list and solver timer
    if (tempLabel != label_) {
      label_ = tempLabel;
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      std::string polyLabel = label_ + ": GmresPolyOp creation time";
      timerPoly_ = Teuchos::TimeMonitor::getNewCounter(polyLabel);
#endif
    }
  }

  // Update timer label
  params_->set("Timer Label", label_);

  // Check if the orthogonalization changed.
  if (params->isParameter("Orthogonalization")) {
    std::string tempOrthoType = params->get("Orthogonalization",orthoType_default_);
    OrthoManagerFactory<ScalarType, MV, OP> factory;
    // Ensure that the specified orthogonalization type is valid.
    if (! factory.isValidName (tempOrthoType)) {
      std::ostringstream os;
      os << "Belos::GCRODRSolMgr: Invalid orthogonalization name \""
         << tempOrthoType << "\".  The following are valid options "
         << "for the \"Orthogonalization\" name parameter: ";
      factory.printValidNames (os);
      throw std::invalid_argument (os.str());
    }
    if (tempOrthoType != orthoType_) {
      orthoType_ = tempOrthoType;
    }
  }

  params_->set("Orthogonalization", orthoType_);

  // Check for a change in verbosity level
  if (params->isParameter("Verbosity")) {
    if (Teuchos::isParameterType<int>(*params,"Verbosity")) {
      verbosity_ = params->get("Verbosity", verbosity_default_);
    } else {
      verbosity_ = (int)Teuchos::getParameter<Belos::MsgType>(*params,"Verbosity");
    }
  }

  // Update parameter in our list.
  params_->set("Verbosity", verbosity_);

  // output stream
  if (params->isParameter("Output Stream")) {
    outputStream_ = Teuchos::getParameter<Teuchos::RCP<std::ostream> >(*params,"Output Stream");
  }

  // Update parameter in our list.
  params_->set("Output Stream", outputStream_);

  // Convergence
  // Check for polynomial convergence tolerance
  if (params->isParameter("Polynomial Tolerance")) {
    if (params->isType<MagnitudeType> ("Polynomial Tolerance")) {
      polyTol_ = params->get ("Polynomial Tolerance",
                              static_cast<MagnitudeType> (DefaultSolverParameters::polyTol));
    }
    else {
      polyTol_ = params->get ("Polynomial Tolerance", DefaultSolverParameters::polyTol);
    }
  }

  // Update parameter in our list and residual tests.
  params_->set("Polynomial Tolerance", polyTol_);

  // Check for maximum polynomial degree
  if (params->isParameter("Random RHS")) {
    randomRHS_ = params->get("Random RHS",randomRHS_default_);
  }

  // Update parameter in our list.
  params_->set("Random RHS", randomRHS_);
  
  
  // Check for polynomial damping
  if (params->isParameter("Damped Poly")) {
    damp_ = params->get("Damped Poly",dampPoly_default_);
  }
  // Update parameter in our list.
  params_->set("Damped Poly", damp_);

  // Check: Should we add roots for stability if needed?
  if (params->isParameter("Add Roots")) {
    addRoots_ = params->get("Add Roots",addRoots_default_);
  }

  // Update parameter in our list.
  params_->set("Add Roots", addRoots_);

  // Create the timers if we need to.
#ifdef BELOS_TEUCHOS_TIME_MONITOR
  if (timerPoly_ == Teuchos::null) {
    std::string polyLabel = label_ + ": GmresPolyOp creation time";
    timerPoly_ = Teuchos::TimeMonitor::getNewCounter(polyLabel);
  }
#endif

  // Check if we are going to perform an outer solve.
  if (outerSolverType_ != "") { 
    hasOuterSolver_ = true;
  }

  // Inform the solver manager that the current parameters were set.
  isSet_ = true;
}


template<class ScalarType, class MV, class OP>
ReturnType GmresPolySolMgr<ScalarType,MV,OP>::solve ()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;

  // Assume convergence is achieved if user does not require strict convergence.
  ReturnType ret = Belos::Converged;

  // Set the current parameters if they were not set before.  NOTE:
  // This may occur if the user generated the solver manager with the
  // default constructor and then didn't set any parameters using
  // setParameters().
  if (! isSet_) {
    setParameters (Teuchos::parameterList (*getValidParameters ()));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    problem_.is_null (), GmresPolySolMgrLinearProblemFailure,
    "Belos::GmresPolySolMgr::solve: The linear problem has not been set yet, "
    "or was set to null.  Please call setProblem() with a nonnull input before "
    "calling solve().");

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! problem_->isProblemSet (), GmresPolySolMgrLinearProblemFailure,
    "Belos::GmresPolySolMgr::solve: The linear problem is not ready.  Please "
    "call setProblem() on the LinearProblem object before calling solve().");

  // If the GMRES polynomial has not been constructed for this
  // (nmatrix, preconditioner) pair, generate it.
  if (!poly_dim_ && maxDegree_) {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor slvtimer(*timerPoly_);
#endif
    poly_Op_ = Teuchos::rcp( new gmres_poly_t( problem_, params_ ) );
    poly_dim_ = poly_Op_->polyDegree();

    TEUCHOS_TEST_FOR_EXCEPTION( !poly_dim_, GmresPolySolMgrPolynomialFailure,
      "Belos::GmresPolyOp: Failed to generate polynomial that satisfied requirements.");
  }


  // Solve the linear system using the polynomial
  if (hasOuterSolver_  && maxDegree_) {

    // Then the polynomial will be used as an operator for an outer solver.
    // Use outer solver parameter list passed in a sublist.
    Belos::GenericSolverFactory<ScalarType, MultiVec<ScalarType>, Operator<ScalarType> > factory;
    RCP<SolverManager<ScalarType, MultiVec<ScalarType>, Operator<ScalarType> > > solver = factory.create( outerSolverType_, outerParams_ );
    TEUCHOS_TEST_FOR_EXCEPTION( solver == Teuchos::null, std::invalid_argument,
      "Belos::GmresPolySolMgr::solve(): Selected solver is not valid.");

    // Create a copy of the linear problem that uses the polynomial as a preconditioner.
    // The original initial solution and right-hand side are thinly wrapped in the gmres_poly_mv_t
    RCP<gmres_poly_mv_t> new_lhs = rcp( new gmres_poly_mv_t( problem_->getLHS() ) );
    RCP<gmres_poly_mv_t> new_rhs = rcp( new gmres_poly_mv_t( rcp_const_cast<MV>( problem_->getRHS() ) ) );
    RCP<gmres_poly_t> A = rcp( new gmres_poly_t( problem_ ) );  // This just performs problem_->applyOp
    RCP<LinearProblem<ScalarType,MultiVec<ScalarType>,Operator<ScalarType> > > newProblem =
      rcp( new LinearProblem<ScalarType,MultiVec<ScalarType>,Operator<ScalarType> >( A, new_lhs, new_rhs ) );
    std::string solverLabel = label_ + ": Hybrid Gmres";
    newProblem->setLabel(solverLabel); 

    // If the preconditioner is left preconditioner, use Gmres poly as a left preconditioner.
    if (problem_->getLeftPrec() != Teuchos::null)
      newProblem->setLeftPrec( poly_Op_ );
    else 
      newProblem->setRightPrec( poly_Op_ );
    // Set the initial residual vector, if it has already been set in the original problem.
    // Don't set the preconditioned residual vector, because it is not the GmresPoly preconditioned residual vector.
    if (problem_->getInitResVec() != Teuchos::null)
      newProblem->setInitResVec( rcp( new gmres_poly_mv_t( rcp_const_cast<MV>( problem_->getInitResVec() ) ) ) );
    newProblem->setProblem();

    solver->setProblem( newProblem );
    
    ret = solver->solve();
    numIters_ = solver->getNumIters();
    loaDetected_ = solver->isLOADetected();
    achievedTol_ = solver->achievedTol();

  } // if (hasOuterSolver_ && maxDegree_)
  else if (hasOuterSolver_) {

    // There is no polynomial, just create the outer solver with the outerSolverType_ and outerParams_.
    Belos::SolverFactory<ScalarType, MV, OP> factory;
    RCP<SolverManager<ScalarType, MV, OP> > solver = factory.create( outerSolverType_, outerParams_ );
    TEUCHOS_TEST_FOR_EXCEPTION( solver == Teuchos::null, std::invalid_argument,
      "Belos::GmresPolySolMgr::solve(): Selected solver is not valid.");

    solver->setProblem( problem_ );
    
    ret = solver->solve();
    numIters_ = solver->getNumIters();
    loaDetected_ = solver->isLOADetected();
    achievedTol_ = solver->achievedTol();

  }
  else if (maxDegree_) {

    // Apply the polynomial to the current linear system
    poly_Op_->ApplyPoly( *problem_->getRHS(), *problem_->getLHS() );
    achievedTol_ = MTS::one();

  }

  return ret;
}


template<class ScalarType, class MV, class OP>
std::string GmresPolySolMgr<ScalarType,MV,OP>::description () const
{
  std::ostringstream out;

  out << "\"Belos::GmresPolySolMgr\": {"
      << "ScalarType: " << Teuchos::TypeNameTraits<ScalarType>::name ()
      << ", Poly Degree: " << poly_dim_
      << ", Poly Max Degree: " << maxDegree_
      << ", Poly Tol: " << polyTol_;
  out << "}";
  return out.str ();
}

} // namespace Belos

#endif // BELOS_GMRES_POLY_SOLMGR_HPP
