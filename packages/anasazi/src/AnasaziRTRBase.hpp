// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziRTRBase.hpp
  \brief Base class for Implicit Riemannian Trust-Region solvers
*/

#ifndef ANASAZI_RTRBASE_HPP
#define ANASAZI_RTRBASE_HPP

#include "AnasaziTypes.hpp"

#include "AnasaziEigensolver.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "AnasaziGenOrthoManager.hpp"
#include "AnasaziSolverUtils.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!     \class Anasazi::RTRBase

        \brief This class is an abstract base class for Implicit Riemannian Trust-Region
        based eigensolvers. The class provides the interfaces shared by the %IRTR
        solvers (e.g., getState() and initialize()) as well as the shared
        implementations (e.g., inner products). 

        %IRTR eigensolvers are capable of solving symmetric/Hermitian
        %eigenvalue problems. These solvers may be used to compute
        either the leftmost (smallest real, "SR") or rightmost (largest real,
        "LR") eigenvalues.  For more information, see the publications at the
        <a href="https://people.scs.fsu.edu/~cbaker/RTRESGEV/">RTR eigensolvers
        page</a>.

        This class is abstract and objects cannot be instantiated. Instead,
        instantiate one of the concrete derived classes: IRTR and SIRTR, 
        the caching and non-caching implementations of this solver. The main difference between
        these solver is the memory allocated by the solvers in support of the %IRTR iteration.

        The reduction in memory usage is effected by eliminating the caching of
        operator applications. This also results in a reduction in vector
        arithmetic required to maintain these caches. The cost is an increase
        in the number of operator applications. For inexpensive operator
        applications, SIRTR should provide better performance over IRTR. As the
        operator applications becomes more expensive, the performance scale
        tips towards the IRTR solver. <b>Note</b>, the trajectory of both
        solvers is identical in exact arithmetic. However, the effects of
        round-off error in the cached results mean that some difference between
        the solvers may exist. This effect is seen when a large number of
        iterations are required to solve the trust-region subproblem in
        solveTRSubproblem(). <b>Also note</b>, the inclusion of auxiliary
        vectors increases the memory requirements of these solvers linearly
        with the number of auxiliary vectors. The required storage is listed in
        the following table:

        <center>
        <table>
        <tr><td align=center colspan=5>Number of vectors (bS == blockSize())</td></tr>
        <tr><td>Solver</td><td>Base requirement</td><td>Generalized/B != null</td><td>Preconditioned</td><td>Generalized and Preconditioned</td></tr>
        <tr><td>IRTR</td><td>10*bS</td><td>11*bS</td><td>12*bS</td><td>13*bS</td></tr>
        <tr><td>SIRTR</td><td>6*bS</td><td>7*bS</td><td>7*bS</td><td>8*bS</td></tr>
        </table>
        </center>

        \ingroup anasazi_solver_framework

        \author Chris Baker
*/

namespace Anasazi {

  //! @name RTRBase Structures
  //@{ 

  /** \brief Structure to contain pointers to RTR state variables.
   *
   * This struct is utilized by RTRBase::initialize() and RTRBase::getState().
   */
  template <class ScalarType, class MV>
  struct RTRState {
    //! The current eigenvectors.
    Teuchos::RCP<const MV> X; 
    //! The image of the current eigenvectors under A, or Teuchos::null is we implement a skinny solver.
    Teuchos::RCP<const MV> AX; 
    //! The image of the current eigenvectors under B, or Teuchos::null if B was not specified.
    Teuchos::RCP<const MV> BX;
    //! The current residual vectors.
    Teuchos::RCP<const MV> R;
    //! The current Ritz values.
    Teuchos::RCP<const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > T;
    /*! \brief The current rho value.
     *  This is only valid if the debugging level of verbosity is enabled.
     */ 
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType rho;
    RTRState() : X(Teuchos::null),AX(Teuchos::null),BX(Teuchos::null),
                 R(Teuchos::null),T(Teuchos::null),rho(0) {};
  };

  //@}

  //! @name RTR Exceptions
  //@{ 

  /** \brief RTRRitzFailure is thrown when the RTR solver is unable to
   *  continue a call to RTRBase::iterate() due to a failure of the algorithm.
   *
   *  This signals that the Rayleigh-Ritz analysis of <tt>X + Eta</tt>
   *  detected ill-conditioning of the projected mass matrix
   *  and the inability to generate a set of orthogonal eigenvectors for 
   *  the projected problem (if thrown from iterate()) or that the analysis of 
   *  the initial iterate failed in RTRBase::initialize().
   *
   *  After catching this exception, the user can recover the subspace via
   *  RTRBase::getState(). This information can be used to restart the solver.
   *
   */
  class RTRRitzFailure : public AnasaziError {public:
    RTRRitzFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  /** \brief RTRInitFailure is thrown when the RTR solver is unable to
   * generate an initial iterate in the RTRBase::initialize() routine. 
   *
   * This exception is thrown from the RTRBase::initialize() method, which is
   * called by the user or from the RTRBase::iterate() method when isInitialized()
   * == \c false.
   *
   */
  class RTRInitFailure : public AnasaziError {public:
    RTRInitFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  /** \brief RTROrthoFailure is thrown when an orthogonalization attempt 
   * fails.
   *
   * This is thrown in one of two scenarios. After preconditioning the residual,
   * the orthogonalization manager is asked to orthogonalize the preconditioned
   * residual (H) against the auxiliary vectors. If full orthogonalization
   * is enabled, H is also orthogonalized against X and P and normalized.
   *
   * The second scenario involves the generation of new X and P from the
   * basis [X H P]. When full orthogonalization is enabled, an attempt is
   * made to select coefficients for X and P so that they will be
   * mutually orthogonal and orthonormal.
   *
   * If either of these attempts fail, the solver throws an RTROrthoFailure
   * exception.
   */
  class RTROrthoFailure : public AnasaziError {public:
    RTROrthoFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};


  //@}


  template <class ScalarType, class MV, class OP>
  class RTRBase : public Eigensolver<ScalarType,MV,OP> {
  public:

    //! @name Constructor/Destructor
    //@{ 

    /*! \brief %RTRBase constructor with eigenproblem, solver utilities, and parameter list of solver options.
     *
     * The RTRBase class is abstract and cannot be instantiated; this constructor is called by derived classes
     * IRTR and RTR.
     */
    RTRBase(const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem, 
            const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
            const Teuchos::RCP<OutputManager<ScalarType> > &printer,
            const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
            const Teuchos::RCP<GenOrthoManager<ScalarType,MV,OP> > &ortho,
                  Teuchos::ParameterList &params,
            const std::string &solverLabel, bool skinnySolver
           );

    //! %RTRBase destructor
    virtual ~RTRBase() {};

    //@}

    //! @name Solver methods
    //@{

    /*! \brief This method performs %RTR iterations until the status test
     * indicates the need to stop or an error occurs (in which case, an
     * exception is thrown).
     *
     * iterate() will first determine whether the solver is initialized; if
     * not, it will call initialize() using default arguments.  After
     * initialization, the solver performs %RTR iterations until the status
     * test evaluates as ::Passed, at which point the method returns to the
     * caller.
     *
     * The %RTR iteration proceeds as follows:
     * -# the trust-region subproblem at \c X is solved for update \c Eta via a call to solveTRSubproblem()
     * -# the new iterate is the Ritz vectors with respect to <tt>X+Eta</tt>
     * -# the eigenproblem residuals are formed with respect to the new iterate
     *
     * The status test is queried at the beginning of the iteration.
     *
     * Possible exceptions thrown include std::logic_error, std::invalid_argument or
     * one of the RTR-specific exceptions.
     *
    */
    virtual void iterate() = 0;

    /*! \brief Initialize the solver to an iterate, optionally providing the
     * Ritz values and residual.
     *
     * The %RTR eigensolver contains a certain amount of state relating to
     * the current iterate.
     *
     * initialize() gives the user the opportunity to manually set these,
     * although this must be done with caution, abiding by the rules
     * given below. All notions of orthogonality and orthonormality are derived
     * from the inner product specified by the orthogonalization manager.
     *
     * \post 
     *   - isInitialized() == true (see post-conditions of isInitialize())
     *
     * The user has the option of specifying any component of the state using
     * initialize(). However, these arguments are assumed to match the
     * post-conditions specified under isInitialized(). Any component of the
     * state (i.e., AX) not given to initialize() will be generated.
     *
     * If the Ritz values relative to <tt>newstate.X</tt> are passed in <tt>newstate.T</tt>,
     * then <tt>newstate.X</tt> is assume to contain Ritz vectors, i.e., <tt>newstate.T</tt> 
     * must be B-orthonormal and it must partially diagonalize A.
     *
     */
    void initialize(RTRState<ScalarType,MV>& newstate);

    /*! \brief Initialize the solver with the initial vectors from the eigenproblem
     *  or random data.
     */
    void initialize();

    /*! \brief Indicates whether the solver has been initialized or not.
     *
     * \return bool indicating the state of the solver.
     * \post
     * If isInitialized() == \c true:
     *   - X is orthogonal to auxiliary vectors and has orthonormal columns
     *   - AX == A*X
     *   - BX == B*X if B != Teuchos::null\n
     *     Otherwise, BX == Teuchos::null
     *   - getRitzValues() returns the sorted Ritz values with respect to X
     *   - getResidualVecs() returns the residual vectors with respect to X
     */
    bool isInitialized() const;

    /*! \brief Get the current state of the eigensolver.
     * 
     * The data is only valid if isInitialized() == \c true. 
     *
     * \returns An RTRState object containing const pointers to the current
     * solver state.
     */
    RTRState<ScalarType,MV> getState() const;

    //@}

    //! @name Status methods
    //@{

    //! \brief Get the current iteration count.
    int getNumIters() const;

    //! \brief Reset the iteration count.
    void resetNumIters();

    /*! \brief Get the Ritz vectors from the previous iteration.
      
        \return A multivector with getBlockSize() vectors containing 
        the sorted Ritz vectors corresponding to the most significant Ritz values.
        The i-th vector of the return corresponds to the i-th Ritz vector; there is no need to use
        getRitzIndex().
     */
    Teuchos::RCP<const MV> getRitzVectors();

    /*! \brief Get the Ritz values from the previous iteration.
     *
     *  \return A vector of length getCurSubspaceDim() containing the Ritz values from the
     *  previous projected eigensolve.
     */
    std::vector<Value<ScalarType> > getRitzValues();

    /*! \brief Get the index used for extracting Ritz vectors from getRitzVectors().
     *
     * Because BlockDavidson is a Hermitian solver, all Ritz values are real and all Ritz vectors can be represented in a 
     * single column of a multivector. Therefore, getRitzIndex() is not needed when using the output from getRitzVectors().
     *
     * \return An \c int vector of size getCurSubspaceDim() composed of zeros.
     */
    std::vector<int> getRitzIndex();

    /*! \brief Get the current residual norms
     *
     *  \return A vector of length getCurSubspaceDim() containing the norms of the
     *  residuals, with respect to the orthogonalization manager norm() method.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getResNorms();


    /*! \brief Get the current residual 2-norms
     *
     *  \return A vector of length getCurSubspaceDim() containing the 2-norms of the
     *  residuals. 
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRes2Norms();


    /*! \brief Get the 2-norms of the Ritz residuals.
     *
     *  \return A vector of length getCurSubspaceDim() containing the 2-norms of the Ritz residuals.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRitzRes2Norms();


    /*! \brief Get the dimension of the search subspace used to generate the current eigenvectors and eigenvalues.
     *
     *  %RTR employs a sequential subspace iteration, maintaining a fixed-rank basis, as opposed to an expanding subspace
     *  mechanism employed by Krylov-subspace solvers like BlockKrylovSchur and BlockDavidson.
     *  
     *  \return An integer specifying the rank of the subspace generated by the eigensolver. If isInitialized() == \c false, 
     *  the return is 0. Otherwise, the return will be getBlockSize().
     */
    int getCurSubspaceDim() const;

    /*! \brief Get the maximum dimension allocated for the search subspace. For %RTR, this always returns getBlockSize().
     */
    int getMaxSubspaceDim() const;

    //@}

    //!  @name Accessor routines from Eigensolver
    //@{

    //! Set a new StatusTest for the solver.
    void setStatusTest(Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test);

    //! Get the current StatusTest used by the solver.
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > getStatusTest() const;

    //! Get a constant reference to the eigenvalue problem.
    const Eigenproblem<ScalarType,MV,OP>& getProblem() const;


    /*! \brief Set the blocksize to be used by the iterative solver in solving
     * this eigenproblem.
     *  
     *  If the block size is reduced, then the new iterate (and residual and
     *  search direction) are chosen as the subset of the current iterate
     *  preferred by the sort manager.  Otherwise, the solver state is set to
     *  uninitialized.
     */
    void setBlockSize(int blockSize);


    //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
    int getBlockSize() const;


    /*! \brief Set the auxiliary vectors for the solver.
     *
     *  Because the current iterate X cannot be assumed orthogonal to the new
     *  auxiliary vectors, a call to setAuxVecs() with a non-empty argument
     *  will reset the solver to the uninitialized state.
     *
     *  In order to preserve the current state, the user will need to extract
     *  it from the solver using getState(), orthogonalize it against the new
     *  auxiliary vectors, and manually reinitialize the solver using
     *  initialize().
     *
     *  <b>NOTE:</b> The requirements of the %IRTR solvers is such that the 
     *  auxiliary vectors must be moved into contiguous storage with the 
     *  current iterate. As a result, the multivector data in \c auxvecs will be copied,
     *  and the multivectors in \c auxvecs will no longer be referenced. The (unchanged) 
     *  internal copies of the auxilliary vectors will be made available to the caller
     *  by the getAuxVecs() routine. This allows the caller to delete the caller's copies and
     *  instead use the copies owned by the solver, avoiding the duplication of data. This 
     *  is not necessary, however. The partitioning of the auxiliary vectors passed to setAuxVecs() will be preserved.
     */
    void setAuxVecs(const Teuchos::Array<Teuchos::RCP<const MV> > &auxvecs);

    //! Get the current auxiliary vectors.
    Teuchos::Array<Teuchos::RCP<const MV> > getAuxVecs() const;

    //@}

    //!  @name Output methods
    //@{

    //! This method requests that the solver print out its current status to screen.
    virtual void currentStatus(std::ostream &os);

    //@}

  protected:
    //
    // Convenience typedefs
    //
    typedef SolverUtils<ScalarType,MV,OP> Utils;
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MAT;
    const MagnitudeType ONE;  
    const MagnitudeType ZERO; 
    const MagnitudeType NANVAL;
    typedef typename std::vector<MagnitudeType>::iterator vecMTiter;
    typedef typename std::vector<ScalarType>::iterator    vecSTiter;
    //
    // Internal structs
    //
    struct CheckList {
      bool checkX, checkAX, checkBX;
      bool checkEta, checkAEta, checkBEta;
      bool checkR, checkQ, checkBR;
      bool checkZ, checkPBX;
      CheckList() : checkX(false),checkAX(false),checkBX(false),
                    checkEta(false),checkAEta(false),checkBEta(false),
                    checkR(false),checkQ(false),checkBR(false),
                    checkZ(false), checkPBX(false) {};
    };
    //
    // Internal methods
    //
    std::string accuracyCheck(const CheckList &chk, const std::string &where) const;
    // solves the model minimization
    virtual void solveTRSubproblem() = 0;
    // Riemannian metric
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType ginner(const MV &xi) const;
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType ginner(const MV &xi, const MV &zeta) const;
    void ginnersep(const MV &xi, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > &d) const;
    void ginnersep(const MV &xi, const MV &zeta, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > &d) const;
    //
    // Classes input through constructor that define the eigenproblem to be solved.
    //
    const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> >     problem_;
    const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> >           sm_;
    const Teuchos::RCP<OutputManager<ScalarType> >          om_;
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> >             tester_;
    const Teuchos::RCP<GenOrthoManager<ScalarType,MV,OP> >  orthman_;
    //
    // Information obtained from the eigenproblem
    //
    Teuchos::RCP<const OP> AOp_;
    Teuchos::RCP<const OP> BOp_;
    Teuchos::RCP<const OP> Prec_;
    bool hasBOp_, hasPrec_, olsenPrec_;
    //
    // Internal timers
    //
    Teuchos::RCP<Teuchos::Time> timerAOp_, timerBOp_, timerPrec_,
                                timerSort_, 
                                timerLocalProj_, timerDS_,
                                timerLocalUpdate_, timerCompRes_,
                                timerOrtho_, timerInit_;
    //
    // Counters
    //
    // Number of operator applications
    int counterAOp_, counterBOp_, counterPrec_;

    //
    // Algorithmic parameters.
    //
    // blockSize_ is the solver block size
    int blockSize_;
    //
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;
    //
    // nevLocal_ reflects how much of the current basis is valid (0 <= nevLocal_ <= blockSize_)
    // this tells us how many of the values in theta_ are valid Ritz values
    int nevLocal_;
    // 
    // are we implementing a skinny solver? (SkinnyIRTR)
    bool isSkinny_;
    // 
    // are we computing leftmost or rightmost eigenvalue?
    bool leftMost_;
    //
    // State Multivecs
    //
    // if we are implementing a skinny solver (SkinnyIRTR), 
    // then some of these will never be allocated
    // 
    // In order to handle auxiliary vectors, we need to handle the projector
    //   P_{[BQ BX],[BQ BX]}
    // Using an orthomanager with B-inner product, this requires calling with multivectors
    // [BQ,BX] and [Q,X].
    // These multivectors must be combined because <[BQ,BX],[Q,X]>_B != I
    // Hence, we will create two multivectors V and BV, which store
    //   V = [Q,X]
    //  BV = [BQ,BX]
    // 
    //  In the context of preconditioning, we may need to apply the projector
    //   P_{prec*[BQ,BX],[BQ,BX]}
    //  Because [BQ,BX] do not change during the supproblem solver, we will apply 
    //  the preconditioner to [BQ,BX] only once. This result is stored in PBV.
    // 
    // X,BX are views into V,BV
    // We don't need views for Q,BQ
    // Inside the subproblem solver, X,BX are constant, so we can allow these
    // views to exist alongside the full view of V,BV without violating
    // view semantics.
    // 
    // Skinny solver allocates 6/7/8 multivectors:
    //    V_, BV_ (if hasB)
    //    PBV_ (if hasPrec and olsenPrec)
    //    R_, Z_  (regardless of hasPrec)
    //    eta_, delta_, Hdelta_
    //
    // Hefty solver allocates 10/11/12/13 multivectors:
    //    V_, AX_, BV_ (if hasB)
    //    PBV_ (if hasPrec and olsenPrec)
    //    R_, Z_ (if hasPrec)
    //    eta_, Aeta_, Beta_
    //    delta_, Adelta_, Bdelta_, Hdelta_
    //
    Teuchos::RCP<MV> V_, BV_, PBV_,                     // V = [Q,X]; B*V; Prec*B*V
                     AX_, R_,                           // A*X_; residual, gradient, and residual of model minimization
                     eta_, Aeta_, Beta_,                // update vector from model minimization
                     delta_, Adelta_, Bdelta_, Hdelta_, // search direction in model minimization
                     Z_;                                // preconditioned residual
    Teuchos::RCP<const MV> X_, BX_;
    // 
    // auxiliary vectors
    Teuchos::Array<Teuchos::RCP<const MV> > auxVecs_;
    int numAuxVecs_;
    //
    // Number of iterations that have been performed.
    int iter_;
    // 
    // Current eigenvalues, residual norms
    std::vector<MagnitudeType> theta_, Rnorms_, R2norms_, ritz2norms_;
    // 
    // are the residual norms current with the residual?
    bool Rnorms_current_, R2norms_current_;
    // 
    // parameters solver and inner solver
    MagnitudeType conv_kappa_, conv_theta_;
    MagnitudeType rho_;
    // 
    // current objective function value
    MagnitudeType fx_;
  };




  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  template <class ScalarType, class MV, class OP>
  RTRBase<ScalarType,MV,OP>::RTRBase(
        const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> >    &problem, 
        const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
        const Teuchos::RCP<OutputManager<ScalarType> >         &printer,
        const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &tester,
        const Teuchos::RCP<GenOrthoManager<ScalarType,MV,OP> > &ortho,
        Teuchos::ParameterList                                 &params,
        const std::string &solverLabel, bool skinnySolver
        ) :
    ONE(Teuchos::ScalarTraits<MagnitudeType>::one()),
    ZERO(Teuchos::ScalarTraits<MagnitudeType>::zero()),
    NANVAL(Teuchos::ScalarTraits<MagnitudeType>::nan()),
    // problem, tools
    problem_(problem), 
    sm_(sorter),
    om_(printer),
    tester_(tester),
    orthman_(ortho),
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    timerAOp_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Operation A*x")),
    timerBOp_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Operation B*x")),
    timerPrec_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Operation Prec*x")),
    timerSort_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Sorting eigenvalues")),
    timerLocalProj_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Local projection")),
    timerDS_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Direct solve")),
    timerLocalUpdate_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Local update")),
    timerCompRes_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Computing residuals")),
    timerOrtho_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Orthogonalization")),
    timerInit_(Teuchos::TimeMonitor::getNewTimer("Anasazi: "+solverLabel+"::Initialization")),
#endif
    counterAOp_(0),
    counterBOp_(0),
    counterPrec_(0),
    // internal data
    blockSize_(0),
    initialized_(false),
    nevLocal_(0),
    isSkinny_(skinnySolver),
    leftMost_(true),
    auxVecs_( Teuchos::Array<Teuchos::RCP<const MV> >(0) ), 
    numAuxVecs_(0),
    iter_(0),
    Rnorms_current_(false),
    R2norms_current_(false),
    conv_kappa_(.1), 
    conv_theta_(1),
    rho_( MAT::nan() ),
    fx_( MAT::nan() )
  {
    TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,std::invalid_argument,
        "Anasazi::RTRBase::constructor: user passed null problem pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(sm_ == Teuchos::null,std::invalid_argument,
        "Anasazi::RTRBase::constructor: user passed null sort manager pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(om_ == Teuchos::null,std::invalid_argument,
        "Anasazi::RTRBase::constructor: user passed null output manager pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(tester_ == Teuchos::null,std::invalid_argument,
        "Anasazi::RTRBase::constructor: user passed null status test pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(orthman_ == Teuchos::null,std::invalid_argument,
        "Anasazi::RTRBase::constructor: user passed null orthogonalization manager pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(problem_->isProblemSet() == false, std::invalid_argument,
        "Anasazi::RTRBase::constructor: problem is not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(problem_->isHermitian() == false, std::invalid_argument,
        "Anasazi::RTRBase::constructor: problem is not Hermitian.");

    // get the problem operators
    AOp_   = problem_->getOperator();
    TEUCHOS_TEST_FOR_EXCEPTION(AOp_ == Teuchos::null, std::invalid_argument,
                       "Anasazi::RTRBase::constructor: problem provides no A matrix.");
    BOp_  = problem_->getM();
    Prec_ = problem_->getPrec();
    hasBOp_ = (BOp_ != Teuchos::null);
    hasPrec_ = (Prec_ != Teuchos::null);
    olsenPrec_ = params.get<bool>("Olsen Prec", true);

    TEUCHOS_TEST_FOR_EXCEPTION(orthman_->getOp() != BOp_,std::invalid_argument,
        "Anasazi::RTRBase::constructor: orthogonalization manager must use mass matrix for inner product.");

    // set the block size and allocate data
    int bs = params.get("Block Size", problem_->getNEV());
    setBlockSize(bs);

    // leftmost or rightmost?
    leftMost_ = params.get("Leftmost",leftMost_);

    conv_kappa_ = params.get("Kappa Convergence",conv_kappa_);
    TEUCHOS_TEST_FOR_EXCEPTION(conv_kappa_ <= 0 || conv_kappa_ >= 1,std::invalid_argument,
                       "Anasazi::RTRBase::constructor: kappa must be in (0,1).");
    conv_theta_ = params.get("Theta Convergence",conv_theta_);
    TEUCHOS_TEST_FOR_EXCEPTION(conv_theta_ <= 0,std::invalid_argument,
                       "Anasazi::RTRBase::constructor: theta must be strictly postitive.");
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void RTRBase<ScalarType,MV,OP>::setBlockSize (int blockSize) 
  {
    // time spent here counts towards timerInit_
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor lcltimer( *timerInit_ );
#endif

    // This routine only allocates space; it doesn't not perform any computation
    // if solver is initialized and size is to be decreased, take the first blockSize vectors of all to preserve state
    // otherwise, shrink/grow/allocate space and set solver to unitialized

    Teuchos::RCP<const MV> tmp;
    // grab some Multivector to Clone
    // in practice, getInitVec() should always provide this, but it is possible to use a 
    // Eigenproblem with nothing in getInitVec() by manually initializing with initialize(); 
    // in case of that strange scenario, we will try to Clone from R_
    // we like R_ for this, because it has minimal size (blockSize_), as opposed to V_ (blockSize_+numAuxVecs_)
    if (blockSize_ > 0) {
      tmp = R_;
    }
    else {
      tmp = problem_->getInitVec();
      TEUCHOS_TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::logic_error,
          "Anasazi::RTRBase::setBlockSize(): Eigenproblem did not specify initial vectors to clone from");
    }

    TEUCHOS_TEST_FOR_EXCEPTION(blockSize <= 0 || blockSize > MVT::GetGlobalLength(*tmp), std::invalid_argument, 
        "Anasazi::RTRBase::setBlockSize was passed a non-positive block size");

    // last chance to quit before causing side-effects
    if (blockSize == blockSize_) {
      // do nothing
      return;
    }

    // clear views
    X_  = Teuchos::null;
    BX_ = Teuchos::null;

    // regardless of whether we preserve any data, any content in V, BV and PBV corresponding to the
    // auxiliary vectors must be retained
    // go ahead and do these first
    // 
    // two cases here: 
    // * we are growing (possibly, from empty) 
    //   any aux data must be copied over, nothing from X need copying
    // * we are shrinking
    //   any aux data must be copied over, go ahead and copy X material (initialized or not)
    //
    if (blockSize > blockSize_)
    {
      // GROWING 
      // get a pointer for Q-related material, and an index for extracting/setting it
      Teuchos::RCP<const MV> Q;
      std::vector<int> indQ(numAuxVecs_);
      for (int i=0; i<numAuxVecs_; i++) indQ[i] = i;
      // if numAuxVecs_ > 0, then necessarily blockSize_ > 0 (we have already been allocated once)
      TEUCHOS_TEST_FOR_EXCEPTION(numAuxVecs_ > 0 && blockSize_ == 0, std::logic_error,
          "Anasazi::RTRBase::setSize(): logic error. Please contact Anasazi team.");
      // V
      if (numAuxVecs_ > 0) Q = MVT::CloneView(*V_,indQ);
      V_ = MVT::Clone(*tmp,numAuxVecs_ + blockSize);
      if (numAuxVecs_ > 0) MVT::SetBlock(*Q,indQ,*V_);
      // BV
      if (hasBOp_) {
        if (numAuxVecs_ > 0) Q = MVT::CloneView(*BV_,indQ);
        BV_ = MVT::Clone(*tmp,numAuxVecs_ + blockSize);
        if (numAuxVecs_ > 0) MVT::SetBlock(*Q,indQ,*BV_);
      }
      else {
        BV_ = V_;
      }
      // PBV
      if (hasPrec_ && olsenPrec_) {
        if (numAuxVecs_ > 0) Q = MVT::CloneView(*PBV_,indQ);
        PBV_ = MVT::Clone(*tmp,numAuxVecs_ + blockSize);
        if (numAuxVecs_ > 0) MVT::SetBlock(*Q,indQ,*PBV_);
      }
    }
    else 
    {
      // SHRINKING
      std::vector<int> indV(numAuxVecs_+blockSize);
      for (int i=0; i<numAuxVecs_+blockSize; i++) indV[i] = i;
      // V
      V_ = MVT::CloneCopy(*V_,indV);
      // BV
      if (hasBOp_) {
        BV_ = MVT::CloneCopy(*BV_,indV);
      }
      else {
        BV_ = V_;
      }
      // PBV
      if (hasPrec_ && olsenPrec_) {
        PBV_ = MVT::CloneCopy(*PBV_,indV);
      }
    }

    if (blockSize < blockSize_) {
      // shrink vectors
      blockSize_ = blockSize;

      theta_.resize(blockSize_);
      ritz2norms_.resize(blockSize_);
      Rnorms_.resize(blockSize_);
      R2norms_.resize(blockSize_);

      if (initialized_) {
        // shrink multivectors with copy
        std::vector<int> ind(blockSize_);
        for (int i=0; i<blockSize_; i++) ind[i] = i;

        // Z can be deleted, no important info there
        Z_ = Teuchos::null;
        
        // we will not use tmp for cloning, clear it and free some space
        tmp = Teuchos::null;

        R_      = MVT::CloneCopy(*R_     ,ind);
        eta_    = MVT::CloneCopy(*eta_   ,ind);
        delta_  = MVT::CloneCopy(*delta_ ,ind);
        Hdelta_ = MVT::CloneCopy(*Hdelta_,ind);
        if (!isSkinny_) {
          AX_     = MVT::CloneCopy(*AX_    ,ind);
          Aeta_   = MVT::CloneCopy(*Aeta_  ,ind);
          Adelta_ = MVT::CloneCopy(*Adelta_,ind);
        }
        else {
          AX_     = Teuchos::null;
          Aeta_   = Teuchos::null;
          Adelta_ = Teuchos::null;
        }

        if (hasBOp_) {
          if (!isSkinny_) {
            Beta_   = MVT::CloneCopy(*Beta_,ind);
            Bdelta_ = MVT::CloneCopy(*Bdelta_,ind);
          }
          else {
            Beta_   = Teuchos::null;
            Bdelta_ = Teuchos::null;
          }
        }
        else {
          Beta_   = eta_;
          Bdelta_ = delta_;
        }
        
        // we need Z if we have a preconditioner
        // we also use Z for temp storage in the skinny solvers
        if (hasPrec_ || isSkinny_) {
          Z_ = MVT::Clone(*V_,blockSize_);
        }
        else {
          Z_ = R_;
        }

      }
      else {
        // release previous multivectors
        // shrink multivectors without copying
        AX_     = Teuchos::null;
        R_      = Teuchos::null;
        eta_    = Teuchos::null;
        Aeta_   = Teuchos::null;
        delta_  = Teuchos::null;
        Adelta_ = Teuchos::null;
        Hdelta_ = Teuchos::null;
        Beta_   = Teuchos::null;
        Bdelta_ = Teuchos::null;
        Z_      = Teuchos::null;

        R_      = MVT::Clone(*tmp,blockSize_);
        eta_    = MVT::Clone(*tmp,blockSize_);
        delta_  = MVT::Clone(*tmp,blockSize_);
        Hdelta_ = MVT::Clone(*tmp,blockSize_);
        if (!isSkinny_) {
          AX_     = MVT::Clone(*tmp,blockSize_);
          Aeta_   = MVT::Clone(*tmp,blockSize_);
          Adelta_ = MVT::Clone(*tmp,blockSize_);
        }

        if (hasBOp_) {
          if (!isSkinny_) {
            Beta_   = MVT::Clone(*tmp,blockSize_);
            Bdelta_ = MVT::Clone(*tmp,blockSize_);
          }
        }
        else {
          Beta_   = eta_;
          Bdelta_ = delta_;
        }

        // we need Z if we have a preconditioner
        // we also use Z for temp storage in the skinny solvers
        if (hasPrec_ || isSkinny_) {
          Z_ = MVT::Clone(*tmp,blockSize_);
        }
        else {
          Z_ = R_;
        }
      } // if initialized_
    } // if blockSize is shrinking
    else {  // if blockSize > blockSize_
      // this is also the scenario for our initial call to setBlockSize(), in the constructor
      initialized_ = false;

      // grow/allocate vectors
      theta_.resize(blockSize,NANVAL);
      ritz2norms_.resize(blockSize,NANVAL);
      Rnorms_.resize(blockSize,NANVAL);
      R2norms_.resize(blockSize,NANVAL);

      // deallocate old multivectors
      AX_     = Teuchos::null;
      R_      = Teuchos::null;
      eta_    = Teuchos::null;
      Aeta_   = Teuchos::null;
      delta_  = Teuchos::null;
      Adelta_ = Teuchos::null;
      Hdelta_ = Teuchos::null;
      Beta_   = Teuchos::null;
      Bdelta_ = Teuchos::null;
      Z_      = Teuchos::null;

      // clone multivectors off of tmp
      R_      = MVT::Clone(*tmp,blockSize);
      eta_    = MVT::Clone(*tmp,blockSize);
      delta_  = MVT::Clone(*tmp,blockSize);
      Hdelta_ = MVT::Clone(*tmp,blockSize);
      if (!isSkinny_) {
        AX_     = MVT::Clone(*tmp,blockSize);
        Aeta_   = MVT::Clone(*tmp,blockSize);
        Adelta_ = MVT::Clone(*tmp,blockSize);
      }

      if (hasBOp_) {
        if (!isSkinny_) {
          Beta_   = MVT::Clone(*tmp,blockSize);
          Bdelta_ = MVT::Clone(*tmp,blockSize);
        }
      }
      else {
        Beta_   = eta_;
        Bdelta_ = delta_;
      }
      if (hasPrec_ || isSkinny_) {
        Z_ = MVT::Clone(*tmp,blockSize);
      }
      else {
        Z_ = R_;
      }
      blockSize_ = blockSize;
    }

    // get view of X from V, BX from BV
    // these are located after the first numAuxVecs columns
    {
      std::vector<int> indX(blockSize_);
      for (int i=0; i<blockSize_; i++) indX[i] = numAuxVecs_+i;
      X_ = MVT::CloneView(*V_,indX);
      if (hasBOp_) {
        BX_ = MVT::CloneView(*BV_,indX);
      }
      else {
        BX_ = X_;
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set a new StatusTest for the solver.
  template <class ScalarType, class MV, class OP>
  void RTRBase<ScalarType,MV,OP>::setStatusTest(Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test) {
    TEUCHOS_TEST_FOR_EXCEPTION(test == Teuchos::null,std::invalid_argument,
        "Anasazi::RTRBase::setStatusTest() was passed a null StatusTest.");
    tester_ = test;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the current StatusTest used by the solver.
  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > RTRBase<ScalarType,MV,OP>::getStatusTest() const {
    return tester_;
  }
  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the auxiliary vectors
  template <class ScalarType, class MV, class OP>
  void RTRBase<ScalarType,MV,OP>::setAuxVecs(const Teuchos::Array<Teuchos::RCP<const MV> > &auxvecs) {
    typedef typename Teuchos::Array<Teuchos::RCP<const MV> >::const_iterator tarcpmv;

    // set new auxiliary vectors
    auxVecs_.resize(0);
    auxVecs_.reserve(auxvecs.size());

    numAuxVecs_ = 0;
    for (tarcpmv v=auxvecs.begin(); v != auxvecs.end(); ++v) {
      numAuxVecs_ += MVT::GetNumberVecs(**v);
    }

    // If the solver has been initialized, X is not necessarily orthogonal to new auxiliary vectors
    if (numAuxVecs_ > 0 && initialized_) {
      initialized_ = false;
    }

    // clear X,BX views
    X_   = Teuchos::null;
    BX_  = Teuchos::null;
    // deallocate, we'll clone off R_ below
    V_   = Teuchos::null;
    BV_  = Teuchos::null;
    PBV_ = Teuchos::null;

    // put auxvecs into V, update BV and PBV if necessary
    if (numAuxVecs_ > 0) {
      V_ = MVT::Clone(*R_,numAuxVecs_ + blockSize_);
      int numsofar = 0;
      for (tarcpmv v=auxvecs.begin(); v != auxvecs.end(); ++v) {
        std::vector<int> ind(MVT::GetNumberVecs(**v));
        for (size_t j=0; j<ind.size(); j++) ind[j] = numsofar++;
        MVT::SetBlock(**v,ind,*V_);
        auxVecs_.push_back(MVT::CloneView(*Teuchos::rcp_static_cast<const MV>(V_),ind));
      }
      TEUCHOS_TEST_FOR_EXCEPTION(numsofar != numAuxVecs_, std::logic_error,
          "Anasazi::RTRBase::setAuxVecs(): Logic error. Please contact Anasazi team.");
      // compute B*V, Prec*B*V
      if (hasBOp_) {
        BV_ = MVT::Clone(*R_,numAuxVecs_ + blockSize_);
        OPT::Apply(*BOp_,*V_,*BV_);
      }
      else {
        BV_ = V_;
      }
      if (hasPrec_ && olsenPrec_) {
        PBV_ = MVT::Clone(*R_,numAuxVecs_ + blockSize_);
        OPT::Apply(*Prec_,*BV_,*V_);
      }
    }
    // 

    if (om_->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkQ = true;
      om_->print( Debug, accuracyCheck(chk, "in setAuxVecs()") );
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  /* Initialize the state of the solver
   * 
   * POST-CONDITIONS:
   *
   * initialized_ == true
   * X is orthonormal, orthogonal to auxVecs_
   * AX = A*X if not skinnny
   * BX = B*X if hasBOp_
   * theta_ contains Ritz values of X
   * R = AX - BX*diag(theta_)
   */
  template <class ScalarType, class MV, class OP>
  void RTRBase<ScalarType,MV,OP>::initialize(RTRState<ScalarType,MV>& newstate)
  {
    // NOTE: memory has been allocated by setBlockSize(). Use SetBlock below; do not Clone
    // NOTE: Overall time spent in this routine is counted to timerInit_; portions will also be counted towards other primitives

    // clear const views to X,BX in V,BV
    // work with temporary non-const views
    X_  = Teuchos::null;
    BX_ = Teuchos::null;
    Teuchos::RCP<MV> X, BX;
    {
      std::vector<int> ind(blockSize_);
      for (int i=0; i<blockSize_; ++i) ind[i] = numAuxVecs_+i;
      X = MVT::CloneViewNonConst(*V_,ind);
      if (hasBOp_) {
        BX = MVT::CloneViewNonConst(*BV_,ind);
      }
      else {
        BX = X;
      }
    }

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor inittimer( *timerInit_ );
#endif

    std::vector<int> bsind(blockSize_);
    for (int i=0; i<blockSize_; i++) bsind[i] = i;

    // in RTR, X (the subspace iterate) is primary
    // the order of dependence follows like so.
    // --init->                 X
    //    --op apply->          AX,BX
    //       --ritz analysis->  theta
    // 
    // if the user specifies all data for a level, we will accept it.
    // otherwise, we will generate the whole level, and all subsequent levels.
    //
    // the data members are ordered based on dependence, and the levels are
    // partitioned according to the amount of work required to produce the
    // items in a level.
    //
    // inconsitent multivectors widths and lengths will not be tolerated, and
    // will be treated with exceptions.

    // set up X, AX, BX: get them from "state" if user specified them
    if (newstate.X != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.X) != MVT::GetGlobalLength(*X),
                          std::invalid_argument, "Anasazi::RTRBase::initialize(newstate): vector length of newstate.X not correct." );
      // newstate.X must have blockSize_ vectors; any more will be ignored
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.X) < blockSize_,
                          std::invalid_argument, "Anasazi::RTRBase::initialize(newstate): newstate.X must have at least block size vectors.");

      // put data in X
      MVT::SetBlock(*newstate.X,bsind,*X);

      // put data in AX
      // if we are implementing a skinny solver, then we don't have memory allocated for AX
      // in this case, point AX at Z (skinny solvers allocate Z) and use it for temporary storage
      // we will clear the pointer later
      if (isSkinny_) {
        AX_ = Z_;
      }
      if (newstate.AX != Teuchos::null) {
        TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.AX) != MVT::GetGlobalLength(*X),
                            std::invalid_argument, "Anasazi::RTRBase::initialize(newstate): vector length of newstate.AX not correct." );
        // newstate.AX must have blockSize_ vectors; any more will be ignored
        TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.AX) < blockSize_,
                            std::invalid_argument, "Anasazi::RTRBase::initialize(newstate): newstate.AX must have at least block size vectors.");
        MVT::SetBlock(*newstate.AX,bsind,*AX_);
      }
      else {
        {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor lcltimer( *timerAOp_ );
#endif
          OPT::Apply(*AOp_,*X,*AX_);
          counterAOp_ += blockSize_;
        }
        // we generated AX; we will generate R as well
        newstate.R = Teuchos::null;
      }

      // put data in BX
      // skinny solvers always allocate BX if hasB, so this is unconditionally appropriate
      if (hasBOp_) {
        if (newstate.BX != Teuchos::null) {
          TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.BX) != MVT::GetGlobalLength(*X),
                              std::invalid_argument, "Anasazi::RTRBase::initialize(newstate): vector length of newstate.BX not correct." );
          // newstate.BX must have blockSize_ vectors; any more will be ignored
          TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.BX) < blockSize_,
                              std::invalid_argument, "Anasazi::RTRBase::initialize(newstate): newstate.BX must have at least block size vectors.");
          MVT::SetBlock(*newstate.BX,bsind,*BX);
        }
        else {
          {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor lcltimer( *timerBOp_ );
#endif
            OPT::Apply(*BOp_,*X,*BX);
            counterBOp_ += blockSize_;
          }
          // we generated BX; we will generate R as well
          newstate.R = Teuchos::null;
        }
      }
      else {
        // the assignment BX_==X_ would be redundant; take advantage of this opportunity to debug a little
        TEUCHOS_TEST_FOR_EXCEPTION(BX != X, std::logic_error, "Anasazi::RTRBase::initialize(): solver invariant not satisfied (BX==X).");
      }

    }
    else {
      // user did not specify X

      // clear state so we won't use any data from it below
      newstate.R  = Teuchos::null;
      newstate.T  = Teuchos::null;

      // generate something and projectAndNormalize
      Teuchos::RCP<const MV> ivec = problem_->getInitVec();
      TEUCHOS_TEST_FOR_EXCEPTION(ivec == Teuchos::null,std::logic_error,
                         "Anasazi::RTRBase::initialize(): Eigenproblem did not specify initial vectors to clone from.");

      int initSize = MVT::GetNumberVecs(*ivec);
      if (initSize > blockSize_) {
        // we need only the first blockSize_ vectors from ivec; get a view of them
        initSize = blockSize_;
        std::vector<int> ind(blockSize_);
        for (int i=0; i<blockSize_; i++) ind[i] = i;
        ivec = MVT::CloneView(*ivec,ind);
      }

      // assign ivec to first part of X
      if (initSize > 0) {
        std::vector<int> ind(initSize);
        for (int i=0; i<initSize; i++) ind[i] = i;
        MVT::SetBlock(*ivec,ind,*X);
      }
      // fill the rest of X with random
      if (blockSize_ > initSize) {
        std::vector<int> ind(blockSize_ - initSize);
        for (int i=0; i<blockSize_ - initSize; i++) ind[i] = initSize + i;
        Teuchos::RCP<MV> rX = MVT::CloneViewNonConst(*X,ind);
        MVT::MvRandom(*rX);
        rX = Teuchos::null;
      }

      // put data in BX
      if (hasBOp_) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerBOp_ );
#endif
        OPT::Apply(*BOp_,*X,*BX);
        counterBOp_ += blockSize_;
      }
      else {
        // the assignment BX==X would be redundant; take advantage of this opportunity to debug a little
        TEUCHOS_TEST_FOR_EXCEPTION(BX != X, std::logic_error, "Anasazi::RTRBase::initialize(): solver invariant not satisfied (BX==X).");
      }
  
      // remove auxVecs from X and normalize it
      if (numAuxVecs_ > 0) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
#endif
        Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummyC;
        int rank = orthman_->projectAndNormalizeMat(*X,auxVecs_,dummyC,Teuchos::null,BX);
        TEUCHOS_TEST_FOR_EXCEPTION(rank != blockSize_, RTRInitFailure,
                           "Anasazi::RTRBase::initialize(): Couldn't generate initial basis of full rank.");
      }
      else {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
#endif
        int rank = orthman_->normalizeMat(*X,Teuchos::null,BX);
        TEUCHOS_TEST_FOR_EXCEPTION(rank != blockSize_, RTRInitFailure,
                           "Anasazi::RTRBase::initialize(): Couldn't generate initial basis of full rank.");
      }

      // put data in AX
      if (isSkinny_) {
        AX_ = Z_;
      }
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerAOp_ );
#endif
        OPT::Apply(*AOp_,*X,*AX_);
        counterAOp_ += blockSize_;
      }

    } // end if (newstate.X != Teuchos::null)


    // set up Ritz values
    if (newstate.T != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION( (signed int)(newstate.T->size()) < blockSize_,
                          std::invalid_argument, "Anasazi::RTRBase::initialize(newstate): newstate.T must contain at least block size Ritz values.");
      for (int i=0; i<blockSize_; i++) {
        theta_[i] = (*newstate.T)[i];
      }
    }
    else {
      // get ritz vecs/vals
      Teuchos::SerialDenseMatrix<int,ScalarType> AA(blockSize_,blockSize_),
                                                 BB(blockSize_,blockSize_),
                                                  S(blockSize_,blockSize_);
      // project A
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerLocalProj_ );
#endif
        MVT::MvTransMv(ONE,*X,*AX_,AA);
        if (hasBOp_) {
          MVT::MvTransMv(ONE,*X,*BX,BB);
        }
      }
      nevLocal_ = blockSize_;

      // solve the projected problem
      // project B
      int ret;
      if (hasBOp_) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerDS_ );
#endif
        ret = Utils::directSolver(blockSize_,AA,Teuchos::rcpFromRef(BB),S,theta_,nevLocal_,1);
      }
      else {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerDS_ );
#endif
        ret = Utils::directSolver(blockSize_,AA,Teuchos::null,S,theta_,nevLocal_,10);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(ret != 0,RTRInitFailure,
          "Anasazi::RTRBase::initialize(): failure solving projected eigenproblem after retraction. LAPACK returns " << ret);
      TEUCHOS_TEST_FOR_EXCEPTION(nevLocal_ != blockSize_,RTRInitFailure,"Anasazi::RTRBase::initialize(): retracted iterate failed in Ritz analysis.");

      // We only have blockSize_ ritz pairs, ergo we do not need to select.
      // However, we still require them to be ordered correctly
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerSort_ );
#endif
        
        std::vector<int> order(blockSize_);
        // 
        // sort the first blockSize_ values in theta_
        sm_->sort(theta_, Teuchos::rcpFromRef(order), blockSize_);   // don't catch exception
        //
        // apply the same ordering to the primitive ritz vectors
        Utils::permuteVectors(order,S);
      }

      // compute ritz residual norms
      {
        Teuchos::BLAS<int,ScalarType> blas;
        Teuchos::SerialDenseMatrix<int,ScalarType> RR(blockSize_,blockSize_);
        // RR = AA*S - BB*S*diag(theta)
        int info;
        if (hasBOp_) {
          info = RR.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,BB,S,ZERO);
          TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, "Anasazi::RTRBase::initialize(): Logic error calling SerialDenseMatrix::multiply.");
        }
        else {
          RR.assign(S);
        }
        for (int i=0; i<blockSize_; i++) {
          blas.SCAL(blockSize_,theta_[i],RR[i],1);
        }
        info = RR.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,AA,S,-ONE);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, "Anasazi::RTRBase::initialize(): Logic error calling SerialDenseMatrix::multiply.");
        for (int i=0; i<blockSize_; i++) {
          ritz2norms_[i] = blas.NRM2(blockSize_,RR[i],1);
        }
      }

      // update the solution
      // use R as local work space
      // Z may already be in use as work space (holding AX if isSkinny)
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerLocalUpdate_ );
#endif
        // X <- X*S
        MVT::MvAddMv( ONE, *X, ZERO, *X, *R_ );        
        MVT::MvTimesMatAddMv( ONE, *R_, S, ZERO, *X );
        // AX <- AX*S
        MVT::MvAddMv( ONE, *AX_, ZERO, *AX_, *R_ );        
        MVT::MvTimesMatAddMv( ONE, *R_, S, ZERO, *AX_ );
        if (hasBOp_) {
          // BX <- BX*S
          MVT::MvAddMv( ONE, *BX, ZERO, *BX, *R_ );        
          MVT::MvTimesMatAddMv( ONE, *R_, S, ZERO, *BX );
        }
      }
    }
    
    // done modifying X,BX; clear temp views and setup const views
    X  = Teuchos::null;
    BX = Teuchos::null;
    {
      std::vector<int> ind(blockSize_);
      for (int i=0; i<blockSize_; ++i) ind[i] = numAuxVecs_+i;
      this->X_ = MVT::CloneView(static_cast<const MV&>(*V_),ind);
      if (hasBOp_) {
        this->BX_ = MVT::CloneView(static_cast<const MV&>(*BV_),ind);
      }
      else {
        this->BX_ = this->X_;
      }
    }


    // get objective function value
    fx_ = std::accumulate(theta_.begin(),theta_.end(),ZERO);

    // set up R
    if (newstate.R != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.R) < blockSize_,
                          std::invalid_argument, "Anasazi::RTRBase::initialize(newstate): newstate.R must have blockSize number of vectors." );
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*R_),
                          std::invalid_argument, "Anasazi::RTRBase::initialize(newstate): vector length of newstate.R not correct." );
      MVT::SetBlock(*newstate.R,bsind,*R_);
    }
    else {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerCompRes_ );
#endif
      // form R <- AX - BX*T
      MVT::MvAddMv(ZERO,*AX_,ONE,*AX_,*R_);
      Teuchos::SerialDenseMatrix<int,ScalarType> T(blockSize_,blockSize_);
      T.putScalar(ZERO);
      for (int i=0; i<blockSize_; i++) T(i,i) = theta_[i];
      MVT::MvTimesMatAddMv(-ONE,*BX_,T,ONE,*R_);
    }

    // R has been updated; mark the norms as out-of-date
    Rnorms_current_ = false;
    R2norms_current_ = false;

    // if isSkinny, then AX currently points to Z for temp storage
    // set AX back to null
    if (isSkinny_) {
      AX_ = Teuchos::null;
    }

    // finally, we are initialized
    initialized_ = true;

    if (om_->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkX = true;
      chk.checkAX = true;
      chk.checkBX = true;
      chk.checkR = true;
      chk.checkQ = true;
      om_->print( Debug, accuracyCheck(chk, "after initialize()") );
    }
  }

  template <class ScalarType, class MV, class OP>
  void RTRBase<ScalarType,MV,OP>::initialize()
  {
    RTRState<ScalarType,MV> empty;
    initialize(empty);
  }




  //////////////////////////////////////////////////////////////////////////////////////////////////
  // compute/return residual B-norms
  template <class ScalarType, class MV, class OP>
  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> 
  RTRBase<ScalarType,MV,OP>::getResNorms() {
    if (Rnorms_current_ == false) {
      // Update the residual norms
      orthman_->norm(*R_,Rnorms_);
      Rnorms_current_ = true;
    }
    return Rnorms_;
  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // compute/return residual 2-norms
  template <class ScalarType, class MV, class OP>
  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> 
  RTRBase<ScalarType,MV,OP>::getRes2Norms() {
    if (R2norms_current_ == false) {
      // Update the residual 2-norms 
      MVT::MvNorm(*R_,R2norms_);
      R2norms_current_ = true;
    }
    return R2norms_;
  }




  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check accuracy, orthogonality, and other debugging stuff
  // 
  // bools specify which tests we want to run (instead of running more than we actually care about)
  //
  // we don't bother checking the following because they are computed explicitly:
  //   AH == A*H
  //
  // 
  // checkX    : X orthonormal
  //             orthogonal to auxvecs
  // checkAX   : check AX == A*X
  // checkBX   : check BX == B*X
  // checkEta  : check that Eta'*B*X == 0
  //             orthogonal to auxvecs
  // checkAEta : check that AEta = A*Eta
  // checkBEta : check that BEta = B*Eta
  // checkR    : check R orthogonal to X
  // checkBR   : check R B-orthogonal to X, valid in and immediately after solveTRSubproblem
  // checkQ    : check that auxiliary vectors are actually orthonormal
  //
  // TODO: 
  //  add checkTheta 
  //
  template <class ScalarType, class MV, class OP>
  std::string RTRBase<ScalarType,MV,OP>::accuracyCheck( const CheckList &chk, const std::string &where ) const 
  {
    using std::setprecision;
    using std::scientific;
    using std::setw;
    using std::endl;
    std::stringstream os;
    MagnitudeType tmp;

    os << " Debugging checks: " << where << endl;

    // X and friends
    if (chk.checkX && initialized_) {
      tmp = orthman_->orthonormError(*X_);
      os << " >> Error in X^H B X == I :    " << scientific << setprecision(10) << tmp << endl;
      for (Array_size_type i=0; i<auxVecs_.size(); i++) {
        tmp = orthman_->orthogError(*X_,*auxVecs_[i]);
        os << " >> Error in X^H B Q[" << i << "] == 0 : " << scientific << setprecision(10) << tmp << endl;
      }
    }
    if (chk.checkAX && !isSkinny_ && initialized_) {
      tmp = Utils::errorEquality(*X_, *AX_, AOp_);
      os << " >> Error in AX == A*X    :    " << scientific << setprecision(10) << tmp << endl;
    }
    if (chk.checkBX && hasBOp_ && initialized_) {
      Teuchos::RCP<MV> BX = MVT::Clone(*this->X_,this->blockSize_);
      OPT::Apply(*BOp_,*this->X_,*BX);
      tmp = Utils::errorEquality(*BX, *BX_);
      os << " >> Error in BX == B*X    :    " << scientific << setprecision(10) << tmp << endl;
    }

    // Eta and friends
    if (chk.checkEta && initialized_) {
      tmp = orthman_->orthogError(*X_,*eta_);
      os << " >> Error in X^H B Eta == 0 :  " << scientific << setprecision(10) << tmp << endl;
      for (Array_size_type i=0; i<auxVecs_.size(); i++) {
        tmp = orthman_->orthogError(*eta_,*auxVecs_[i]);
        os << " >> Error in Eta^H B Q[" << i << "]==0 : " << scientific << setprecision(10) << tmp << endl;
      }
    }
    if (chk.checkAEta && !isSkinny_ && initialized_) {
      tmp = Utils::errorEquality(*eta_, *Aeta_, AOp_);
      os << " >> Error in AEta == A*Eta    :    " << scientific << setprecision(10) << tmp << endl;
    }
    if (chk.checkBEta && !isSkinny_ && hasBOp_ && initialized_) {
      tmp = Utils::errorEquality(*eta_, *Beta_, BOp_);
      os << " >> Error in BEta == B*Eta    :    " << scientific << setprecision(10) << tmp << endl;
    }

    // R: this is not B-orthogonality, but standard euclidean orthogonality
    if (chk.checkR && initialized_) {
      Teuchos::SerialDenseMatrix<int,ScalarType> xTx(blockSize_,blockSize_);
      MVT::MvTransMv(ONE,*X_,*R_,xTx);
      tmp = xTx.normFrobenius();
      os << " >> Error in X^H R == 0   :    " << scientific << setprecision(10) << tmp << endl;
    }

    // BR: this is B-orthogonality: this is only valid inside and immediately after solveTRSubproblem 
    // only check if B != I, otherwise, it is equivalent to the above test
    if (chk.checkBR && hasBOp_ && initialized_) {
      Teuchos::SerialDenseMatrix<int,ScalarType> xTx(blockSize_,blockSize_);
      MVT::MvTransMv(ONE,*BX_,*R_,xTx);
      tmp = xTx.normFrobenius();
      os << " >> Error in X^H B R == 0 :    " << scientific << setprecision(10) << tmp << endl;
    }

    // Z: Z is preconditioned R, should be on tangent plane
    if (chk.checkZ && initialized_) {
      Teuchos::SerialDenseMatrix<int,ScalarType> xTx(blockSize_,blockSize_);
      MVT::MvTransMv(ONE,*BX_,*Z_,xTx);
      tmp = xTx.normFrobenius();
      os << " >> Error in X^H B Z == 0 :    " << scientific << setprecision(10) << tmp << endl;
    }

    // Q
    if (chk.checkQ) {
      for (Array_size_type i=0; i<auxVecs_.size(); i++) {
        tmp = orthman_->orthonormError(*auxVecs_[i]);
        os << " >> Error in Q[" << i << "]^H B Q[" << i << "]==I: " << scientific << setprecision(10) << tmp << endl;
        for (Array_size_type j=i+1; j<auxVecs_.size(); j++) {
          tmp = orthman_->orthogError(*auxVecs_[i],*auxVecs_[j]);
          os << " >> Error in Q[" << i << "]^H B Q[" << j << "]==0: " << scientific << setprecision(10) << tmp << endl;
        }
      }
    }
    os << endl;
    return os.str();
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Print the current status of the solver
  template <class ScalarType, class MV, class OP>
  void 
  RTRBase<ScalarType,MV,OP>::currentStatus(std::ostream &os) 
  {
    using std::setprecision;
    using std::scientific;
    using std::setw;
    using std::endl;

    os <<endl;
    os <<"================================================================================" << endl;
    os << endl;
    os <<"                              RTRBase Solver Status" << endl;
    os << endl;
    os <<"The solver is "<<(initialized_ ? "initialized." : "not initialized.") << endl;
    os <<"The number of iterations performed is " << iter_       << endl;
    os <<"The current block size is             " << blockSize_  << endl;
    os <<"The number of auxiliary vectors is    " << numAuxVecs_ << endl;
    os <<"The number of operations A*x    is " << counterAOp_   << endl;
    os <<"The number of operations B*x    is " << counterBOp_    << endl;
    os <<"The number of operations Prec*x is " << counterPrec_ << endl;
    os <<"The most recent rho was " << scientific << setprecision(10) << rho_ << endl;
    os <<"The current value of f(x) is " << scientific << setprecision(10) << fx_ << endl;

    if (initialized_) {
      os << endl;
      os <<"CURRENT EIGENVALUE ESTIMATES             "<<endl;
      os << setw(20) << "Eigenvalue" 
         << setw(20) << "Residual(B)"
         << setw(20) << "Residual(2)"
         << endl;
      os <<"--------------------------------------------------------------------------------"<<endl;
      for (int i=0; i<blockSize_; i++) {
        os << scientific << setprecision(10) << setw(20) << theta_[i];
        if (Rnorms_current_) os << scientific << setprecision(10) << setw(20) << Rnorms_[i];
        else os << scientific << setprecision(10) << setw(20) << "not current";
        if (R2norms_current_) os << scientific << setprecision(10) << setw(20) << R2norms_[i];
        else os << scientific << setprecision(10) << setw(20) << "not current";
        os << endl;
      }
    }
    os <<"================================================================================" << endl;
    os << endl;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Inner product 1
  template <class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
  RTRBase<ScalarType,MV,OP>::ginner(const MV &xi) const 
  {
    std::vector<MagnitudeType> d(MVT::GetNumberVecs(xi));
    MVT::MvNorm(xi,d);
    MagnitudeType ret = 0;
    for (vecMTiter i=d.begin(); i != d.end(); ++i) {
      ret += (*i)*(*i);
    }
    return ret;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Inner product 2
  template <class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
  RTRBase<ScalarType,MV,OP>::ginner(const MV &xi, const MV &zeta) const 
  {
    std::vector<ScalarType> d(MVT::GetNumberVecs(xi));
    MVT::MvDot(xi,zeta,d);
    return SCT::real(std::accumulate(d.begin(),d.end(),SCT::zero()));
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Inner product 1 without trace accumulation
  template <class ScalarType, class MV, class OP>
  void RTRBase<ScalarType,MV,OP>::ginnersep(
      const MV &xi, 
      std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &d) const 
  {
    MVT::MvNorm(xi,d);
    for (vecMTiter i=d.begin(); i != d.end(); ++i) {
      (*i) = (*i)*(*i);
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Inner product 2 without trace accumulation
  template <class ScalarType, class MV, class OP>
  void RTRBase<ScalarType,MV,OP>::ginnersep(
      const MV &xi, const MV &zeta, 
      std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &d) const 
  {
    std::vector<ScalarType> dC(MVT::GetNumberVecs(xi));
    MVT::MvDot(xi,zeta,dC);
    vecMTiter iR=d.begin(); 
    vecSTiter iS=dC.begin();
    for (; iR != d.end(); iR++, iS++) {
      (*iR) = SCT::real(*iS);
    }
  }

  template <class ScalarType, class MV, class OP>
  Teuchos::Array<Teuchos::RCP<const MV> > RTRBase<ScalarType,MV,OP>::getAuxVecs() const {
    return auxVecs_;
  }

  template <class ScalarType, class MV, class OP>
  int RTRBase<ScalarType,MV,OP>::getBlockSize() const { 
    return(blockSize_); 
  }

  template <class ScalarType, class MV, class OP>
  const Eigenproblem<ScalarType,MV,OP>& RTRBase<ScalarType,MV,OP>::getProblem() const { 
    return(*problem_); 
  }

  template <class ScalarType, class MV, class OP>
  int RTRBase<ScalarType,MV,OP>::getMaxSubspaceDim() const {
    return blockSize_;
  }

  template <class ScalarType, class MV, class OP>
  int RTRBase<ScalarType,MV,OP>::getCurSubspaceDim() const 
  {
    if (!initialized_) return 0;
    return nevLocal_;
  }

  template <class ScalarType, class MV, class OP>
  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> 
  RTRBase<ScalarType,MV,OP>::getRitzRes2Norms() 
  {
    std::vector<MagnitudeType> ret = ritz2norms_;
    ret.resize(nevLocal_);
    return ret;
  }

  template <class ScalarType, class MV, class OP>
  std::vector<Value<ScalarType> > 
  RTRBase<ScalarType,MV,OP>::getRitzValues() 
  {
    std::vector<Value<ScalarType> > ret(nevLocal_);
    for (int i=0; i<nevLocal_; i++) {
      ret[i].realpart = theta_[i];
      ret[i].imagpart = ZERO;
    }
    return ret;
  }

  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<const MV> 
  RTRBase<ScalarType,MV,OP>::getRitzVectors() 
  {
    return X_;
  }

  template <class ScalarType, class MV, class OP>
  void RTRBase<ScalarType,MV,OP>::resetNumIters() 
  { 
    iter_=0; 
  }

  template <class ScalarType, class MV, class OP>
  int RTRBase<ScalarType,MV,OP>::getNumIters() const 
  { 
    return iter_; 
  }

  template <class ScalarType, class MV, class OP>
  RTRState<ScalarType,MV> RTRBase<ScalarType,MV,OP>::getState() const 
  {
    RTRState<ScalarType,MV> state;
    state.X = X_;
    state.AX = AX_;
    if (hasBOp_) {
      state.BX = BX_;
    }
    else {
      state.BX = Teuchos::null;
    }
    state.rho = rho_;
    state.R = R_;
    state.T = Teuchos::rcp(new std::vector<MagnitudeType>(theta_));
    return state;
  }

  template <class ScalarType, class MV, class OP>
  bool RTRBase<ScalarType,MV,OP>::isInitialized() const 
  { 
    return initialized_; 
  }

  template <class ScalarType, class MV, class OP>
  std::vector<int> RTRBase<ScalarType,MV,OP>::getRitzIndex() 
  {
    std::vector<int> ret(nevLocal_,0);
    return ret;
  }


} // end Anasazi namespace

#endif // ANASAZI_RTRBASE_HPP
