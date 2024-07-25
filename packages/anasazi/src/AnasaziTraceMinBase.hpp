// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// TODO: Modify code so R is not computed unless needed

/*! \file AnasaziTraceMinBase.hpp
  \brief Abstract base class for trace minimization eigensolvers
*/

#ifndef ANASAZI_TRACEMIN_BASE_HPP
#define ANASAZI_TRACEMIN_BASE_HPP

#include "AnasaziBasicSort.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEigensolver.hpp"
#include "AnasaziMatOrthoManager.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziSaddleOperator.hpp"
#include "AnasaziSaddleContainer.hpp"
#include "AnasaziSolverUtils.hpp"
#include "AnasaziTraceMinRitzOp.hpp"
#include "AnasaziTraceMinTypes.hpp"
#include "AnasaziTypes.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_TimeMonitor.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace Anasazi {
/**
 * @namespace Experimental
 * Namespace for new Anasazi features that are not ready for public release,
 * but are ready for evaluation by friendly expert users.
 *
 * \warning Expect header files, classes, functions, and other interfaces to change or disappear.
 * Anything in this namespace is under active development and evaluation. Documentation may be
 * sparse or not exist yet. If you understand these caveats and accept them, please feel free to
 * take a look inside and try things out.
 */
namespace Experimental {

  //! @name TraceMinBase Structures
  //@{

  /** \brief Structure to contain pointers to TraceMinBase state variables.
   *
   * This struct is utilized by TraceMinBase::initialize() and TraceMinBase::getState().
   */
  template <class ScalarType, class MV>
  struct TraceMinBaseState {
    //! \brief The current dimension of the solver.
    int curDim;
    /*! \brief The current basis.
     *
     * V has TraceMinBase::getMaxSubspaceDim() vectors, but only the first \c curDim are valid.
     */
    RCP<const MV> V;
    //! The image of V under K
    RCP<const MV> KV;
    //! The image of V under M, or Teuchos::null if M was not specified
    RCP<const MV> MopV;
    //! The current eigenvectors.
    RCP<const MV> X;
    //! The image of the current eigenvectors under K.
    RCP<const MV> KX;
    //! The image of the current eigenvectors under M, or Teuchos::null if M was not specified.
    RCP<const MV> MX;
    //! The current residual vectors
    RCP<const MV> R;
    //! The current Ritz values. This vector is a copy of the internal data.
    RCP<const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > T;
    /*! \brief The current projected K matrix.
     *
     * KK is of order TraceMinBase::getMaxSubspaceDim(), but only the principal submatrix of order \c curDim is meaningful. It is Hermitian in memory.
     *
     */
    RCP<const Teuchos::SerialDenseMatrix<int,ScalarType> > KK;
    //! The current Ritz vectors.
    RCP<const Teuchos::SerialDenseMatrix<int,ScalarType> > RV;
    //! Whether V has been projected and orthonormalized already
    bool isOrtho;
    //! Number of unconverged eigenvalues
    int NEV;
    //! Largest safe shift
    ScalarType largestSafeShift;
    //! Current Ritz shifts
    RCP< const std::vector<ScalarType> > ritzShifts;
    TraceMinBaseState() : curDim(0), V(Teuchos::null), KV(Teuchos::null), MopV(Teuchos::null),
                          X(Teuchos::null), KX(Teuchos::null), MX(Teuchos::null), R(Teuchos::null),
                          T(Teuchos::null), KK(Teuchos::null), RV(Teuchos::null), isOrtho(false),
                          NEV(0), largestSafeShift(Teuchos::ScalarTraits<ScalarType>::zero()),
                          ritzShifts(Teuchos::null) {}
  };

  //@}

  //! @name TraceMinBase Exceptions
  //@{

  /** \brief TraceMinBaseInitFailure is thrown when the TraceMinBase solver is unable to
   * generate an initial iterate in the TraceMinBase::initialize() routine.
   *
   * This exception is thrown from the TraceMinBase::initialize() method, which is
   * called by the user or from the TraceMinBase::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this exception is thrown,
   * TraceMinBase::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the solver.
   *
   */
  class TraceMinBaseInitFailure : public AnasaziError {public:
    TraceMinBaseInitFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  /** \brief TraceMinBaseOrthoFailure is thrown when the orthogonalization manager is
   * unable to orthogonalize the vectors in the current basis.
   */
  class TraceMinBaseOrthoFailure : public AnasaziError {public:
    TraceMinBaseOrthoFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  //@}

  /*! \class TraceMinBase

      \brief This is an abstract base class for the trace minimization eigensolvers.

      For more information, please see Anasazi::TraceMin (with constant subspace dimension)
      and Anasazi::TraceMinDavidson (with expanding subspaces)

      \ingroup anasazi_solver_framework

      \author Alicia Klinvex
  */

  template <class ScalarType, class MV, class OP>
  class TraceMinBase : public Eigensolver<ScalarType,MV,OP> {
  public:
    //! @name Constructor/Destructor
    //@{

    /*! \brief %TraceMinBase constructor with eigenproblem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the eigensolver, in addition
     * to a parameter list of options for the eigensolver. These options include the following:
     *   - \c "Saddle Solver Type" - a \c string specifying how to solve the saddle point problem arising at each iteration.
     *        Options are "Projected Krylov", "Schur Complement", and "Block Diagonal Preconditioned Minres". Default: "Projected Krylov"
     *         - \c "Projected Krylov": Uses projected-minres to solve the problem.
     *         - \c "Schur Complement": Explicitly forms the (inexact) Schur complement using minres.
     *         - \c "Block Diagonal Preconditioned Minres": Uses a block preconditioner on the entire saddle point problem.  For more information, please see "Overview of Anasazi and its newest eigensolver, TraceMin" on the main Anasazi page.
     *        We recommend using "Projected Krylov" in the absence of preconditioning.  If you want to use a preconditioner, "Block Diagonal Preconditioned Minres" is recommended.
     *        "Schur Complement" mainly exists for special use cases.
     *   - Ritz shift parameters
     *      - \c "When To Shift" - a \c string specifying when Ritz shifts should be performed. Options are "Never", "After Trace Levels", and "Always". Default: "Always"
     *         - \c "Never": Do not perform Ritz shifts.  This option produces guaranteed convergence but converges linearly.  Not recommended.
     *         - \c "After Trace Levels": Do not perform Ritz shifts until the trace of \f$X^TKX\f$ has stagnated (i.e. the relative change in trace has become small).
     *              The \c MagnitudeType specifying how small the relative change in trace must become may be provided via the parameter \c "Trace Threshold", whose default value is 0.02.
     *         - \c "Always": Always attempt to use Ritz shifts.
     *      - \c "How To Choose Shift" - a \c string specifying how to choose the Ritz shifts (assuming Ritz shifts are being used).
     *           Options are "Largest Converged", "Adjusted Ritz Values", and "Ritz Values". Default: "Adjusted Ritz Values"
     *         - \c "Largest Converged": Ritz shifts are chosen to be the largest converged eigenvalue.  Until an eigenvalue converges, the Ritz shifts are all 0.
     *         - \c "Adjusted Ritz Values": Ritz shifts are chosen based on the Ritz values and their associated residuals in such a way as to guarantee global convergence.
     *              This method is described in "The trace minimization method for the symmetric generalized eigenvalue problem."
     *         - \c "Ritz Values": Ritz shifts are chosen to equal the Ritz values.  This does NOT guarantee global convergence.
     *      - \c "Use Multiple Shifts" - a \c bool specifying whether to use one or many Ritz shifts (assuming shifting is enabled). Default: true
     *
     * Anasazi's trace minimization solvers are still in development, and we plan to add additional features in the future including additional saddle point solvers.
     */
    TraceMinBase( const RCP<Eigenproblem<ScalarType,MV,OP> >    &problem,
                   const RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
                   const RCP<OutputManager<ScalarType> >         &printer,
                   const RCP<StatusTest<ScalarType,MV,OP> >      &tester,
                   const RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
                   Teuchos::ParameterList &params
                 );

    //! %Anasazi::TraceMinBase destructor.
    virtual ~TraceMinBase();
    //@}


    //! @name Solver methods
    //@{

    /*! \brief This method performs trace minimization iterations until the status
     * test indicates the need to stop or an error occurs (in which case, an
     * appropriate exception is thrown).
     *
     * iterate() will first determine whether the solver is initialized; if
     * not, it will call initialize(). After
     * initialization, the solver performs TraceMin iterations until the
     * status test evaluates as ::Passed, at which point the method returns to
     * the caller.
     *
     * The trace minimization iteration proceeds as follows:
     * -# Solve the saddle point problem to obtain Delta
     * -# Add Delta to the basis
     *    In TraceMin, this is done by computing V := X - Delta, then projecting and normalizing.
     *    In TraceMinDavidson, this is done by computing V := [V Delta], then projecting and normalizing
     * -# Compute the Ritz pairs.
     * -# Update the residual.
     *
     * The status test is queried at the beginning of the iteration.
     *
     * Possible exceptions thrown include std::invalid_argument or
     * one of the TraceMinBase-specific exceptions.
     */
    void iterate();

    void harmonicIterate();

    /*! \brief Initialize the solver to an iterate, optionally providing the
     * other members of the state.
     *
     * The %TraceMinBase eigensolver contains a certain amount of state,
     * including the current Krylov basis, the current eigenvectors,
     * the current residual, etc. (see getState())
     *
     * initialize() gives the user the opportunity to manually set these,
     * although this must be done with caution, as the validity of the
     * user input will not be checked.
     *
     * \post
     * <li>isInitialized() == \c true (see post-conditions of isInitialize())
     *
     * The user has the option of specifying any component of the state using
     * initialize(). However, these arguments are assumed to match the
     * post-conditions specified under isInitialized(). Any component of the
     * state (i.e., KX) not given to initialize() will be generated.
     *
     * Note, for any pointer in \c newstate which directly points to the multivectors in
     * the solver, the data is not copied.
     */
    void initialize(TraceMinBaseState<ScalarType,MV>& newstate);

    void harmonicInitialize(TraceMinBaseState<ScalarType,MV> newstate);

    /*! \brief Initialize the solver with the initial vectors from the eigenproblem
     *  or random data.
     */
    void initialize();

    /*! \brief Indicates whether the solver has been initialized or not.
     *
     * \return bool indicating the state of the solver.
     * \post
     * If isInitialized() == \c true:
     *    - getCurSubspaceDim() > 0 and is a multiple of getBlockSize()
     *    - the first getCurSubspaceDim() vectors of V are orthogonal to auxiliary vectors and have orthonormal columns
     *    - the principal submatrix of order getCurSubspaceDim() of KK contains the projected eigenproblem matrix
     *    - X contains the Ritz vectors with respect to the current Krylov basis
     *    - T contains the Ritz values with respect to the current Krylov basis
     *    - KX == Op*X
     *    - MX == M*X if M != Teuchos::null\n
     *      Otherwise, MX == Teuchos::null
     *    - R contains the residual vectors with respect to X
     */
    bool isInitialized() const;

    /*! \brief Get access to the current state of the eigensolver.
     *
     * The data is only valid if isInitialized() == \c true.
     *
     * \returns A TraceMinBaseState object containing const pointers to the current
     * solver state. Note, these are direct pointers to the multivectors; they are not
     * pointers to views of the multivectors.
     */
    TraceMinBaseState<ScalarType,MV> getState() const;

    //@}


    //! @name Status methods
    //@{

    //! \brief Get the current iteration count.
    int getNumIters() const;

    //! \brief Reset the iteration count.
    void resetNumIters();

    /*! \brief Get access to the current Ritz vectors.

        \return A multivector with getBlockSize() vectors containing
        the sorted Ritz vectors corresponding to the most significant Ritz values.
        The i-th vector of the return corresponds to the i-th Ritz vector; there is no need to use
        getRitzIndex().
     */
    RCP<const MV> getRitzVectors();

    /*! \brief Get the Ritz values for the previous iteration.
     *
     *  \return A vector of length getCurSubspaceDim() containing the Ritz values from the
     *  previous projected eigensolve.
     */
    std::vector<Value<ScalarType> > getRitzValues();


    /*! \brief Get the index used for extracting individual Ritz vectors from getRitzVectors().
     *
     * Because the trace minimization methods are a Hermitian solvers, all Ritz values are real
     * and all Ritz vectors can be represented in a single column of a multivector. Therefore,
     * getRitzIndex() is not needed when using the output from getRitzVectors().
     *
     * \return An \c int vector of size getCurSubspaceDim() composed of zeros.
     */
    std::vector<int> getRitzIndex();


    /*! \brief Get the current residual norms, computing the norms if they are not up-to-date with the current residual vectors.
     *
     *  \return A vector of length getCurSubspaceDim() containing the norms of the
     *  residuals, with respect to the orthogonalization manager's norm() method.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getResNorms();


    /*! \brief Get the current residual 2-norms, computing the norms if they are not up-to-date with the current residual vectors.
     *
     *  \return A vector of length getCurSubspaceDim() containing the 2-norms of the
     *  current residuals.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRes2Norms();


    /*! \brief Get the 2-norms of the residuals.
     *
     * The Ritz residuals are not defined for trace minimization iterations. Hence, this method returns the
     * 2-norms of the direct residuals, and is equivalent to calling getRes2Norms().
     *
     *  \return A vector of length getBlockSize() containing the 2-norms of the direct residuals.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRitzRes2Norms();

    /*! \brief Get the dimension of the search subspace used to generate the current eigenvectors and eigenvalues.
     *
     *  \return An integer specifying the rank of the Krylov subspace currently in use by the eigensolver. If isInitialized() == \c false,
     *  the return is 0. Otherwise, it will be some strictly positive multiple of getBlockSize().
     */
    int getCurSubspaceDim() const;

    //! Get the maximum dimension allocated for the search subspace. For the trace minimization methods, this always returns numBlocks*blockSize.
    int getMaxSubspaceDim() const;

    //@}


    //! @name Accessor routines from Eigensolver
    //@{

    //! Set a new StatusTest for the solver.
    void setStatusTest(RCP<StatusTest<ScalarType,MV,OP> > test);

    //! Get the current StatusTest used by the solver.
    RCP<StatusTest<ScalarType,MV,OP> > getStatusTest() const;

    //! Get a constant reference to the eigenvalue problem.
    const Eigenproblem<ScalarType,MV,OP>& getProblem() const;

    /*! \brief Set the blocksize.
     *
     * This method is required to support the interface provided by Eigensolver. However, the preferred method
     * of setting the allocated size for the TraceMinBase eigensolver is setSize(). In fact, setBlockSize()
     * simply calls setSize(), maintaining the current number of blocks.
     *
     * The block size determines the number of Ritz vectors and values that are computed on each iteration, thereby
     * determining the increase in the subspace dimension at each iteration.
     */
    void setBlockSize(int blockSize);

    //! Get the blocksize used by the iterative solver.
    int getBlockSize() const;

    /*! \brief Set the auxiliary vectors for the solver.
     *
     *  Auxiliary vectors are ones that you want your eigenvectors to be
     *  held orthogonal to.  One example of where you may want to use this
     *  is in the computation of the Fiedler vector, where you would likely
     *  want to project against the vector of all 1s.
     *
     *  Because the current basis V cannot be assumed
     *  orthogonal to the new auxiliary vectors, a call to setAuxVecs() will
     *  reset the solver to the uninitialized state. This happens only in the
     *  case where the new auxiliary vectors have a combined dimension of
     *  greater than zero.
     *
     *  In order to preserve the current state, the user will need to extract
     *  it from the solver using getState(), orthogonalize it against the
     *  new auxiliary vectors, and reinitialize using initialize().
     */
    void setAuxVecs(const Teuchos::Array<RCP<const MV> > &auxvecs);

    //! Get the auxiliary vectors for the solver.
    Teuchos::Array<RCP<const MV> > getAuxVecs() const;

    //@}

    //! @name BlockBase-specific accessor routines
    //@{

    /*! \brief Set the blocksize and number of blocks to be used by the
     * iterative solver in solving this eigenproblem.
     *
     *  Changing either the block size or the number of blocks will reset the
     *  solver to an uninitialized state.
     */
    void setSize(int blockSize, int numBlocks);

    //! @name Output methods
    //@{

    //! This method requests that the solver print out its current status to the given output stream.
    void currentStatus(std::ostream &os);

    //@}

  protected:
    //
    // Convenience typedefs
    //
    typedef SolverUtils<ScalarType,MV,OP>                 Utils;
    typedef MultiVecTraits<ScalarType,MV>                 MVT;
    typedef OperatorTraits<ScalarType,MV,OP>              OPT;
    typedef Teuchos::ScalarTraits<ScalarType>             SCT;
    typedef typename SCT::magnitudeType                   MagnitudeType;
    typedef TraceMinRitzOp<ScalarType,MV,OP>              tracemin_ritz_op_type;
    typedef SaddleContainer<ScalarType,MV>                saddle_container_type;
    typedef SaddleOperator<ScalarType,MV,tracemin_ritz_op_type>  saddle_op_type;
    const MagnitudeType ONE;
    const MagnitudeType ZERO;
    const MagnitudeType NANVAL;
    //
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    const RCP<Eigenproblem<ScalarType,MV,OP> >     problem_;
    const RCP<SortManager<MagnitudeType> >         sm_;
    const RCP<OutputManager<ScalarType> >          om_;
    RCP<StatusTest<ScalarType,MV,OP> >             tester_;
    const RCP<MatOrthoManager<ScalarType,MV,OP> >  orthman_;
    //
    // Information obtained from the eigenproblem
    //
    RCP<const OP> Op_;
    RCP<const OP> MOp_;
    RCP<const OP> Prec_;
    bool hasM_;
    //
    // Internal timers
    // TODO: Fix the timers
    //
    RCP<Teuchos::Time> timerOp_, timerMOp_, timerSaddle_, timerSortEval_, timerDS_,
                                timerLocal_, timerCompRes_, timerOrtho_, timerInit_;
    //
    // Internal structs
    // TODO: Fix the checks
    //
    struct CheckList {
      bool checkV, checkX, checkMX,
           checkKX, checkQ, checkKK;
      CheckList() : checkV(false),checkX(false),
                    checkMX(false),checkKX(false),
                    checkQ(false),checkKK(false) {};
    };
    //
    // Internal methods
    //
    std::string accuracyCheck(const CheckList &chk, const std::string &where) const;

    //
    // Counters
    //
    int count_ApplyOp_, count_ApplyM_;

    //
    // Algorithmic parameters.
    //
    // blockSize_ is the solver block size; it controls the number of vectors added to the basis on each iteration.
    int blockSize_;
    // numBlocks_ is the size of the allocated space for the basis, in blocks.
    int numBlocks_;

    //
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;
    //
    // curDim_ reflects how much of the current basis is valid
    // NOTE: 0 <= curDim_ <= blockSize_*numBlocks_
    // this also tells us how many of the values in theta_ are valid Ritz values
    int curDim_;
    //
    // State Multivecs
    // V is the current basis
    // X is the current set of eigenvectors
    // R is the residual
    RCP<MV> X_, KX_, MX_, KV_, MV_, R_, V_;
    //
    // Projected matrices
    //
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > KK_, ritzVecs_;
    //
    // auxiliary vectors
    Teuchos::Array<RCP<const MV> > auxVecs_, MauxVecs_;
    int numAuxVecs_;
    //
    // Number of iterations that have been performed.
    int iter_;
    //
    // Current eigenvalues, residual norms
    std::vector<MagnitudeType> theta_, Rnorms_, R2norms_;
    //
    // are the residual norms current with the residual?
    bool Rnorms_current_, R2norms_current_;

    // TraceMin specific variables
    RCP<tracemin_ritz_op_type> ritzOp_;    // Operator which incorporates Ritz shifts
    enum SaddleSolType saddleSolType_;                          // Type of saddle point solver being used
    bool previouslyLeveled_;                                    // True if the trace already leveled off
    MagnitudeType previousTrace_;                               // Trace of X'KX at the previous iteration
    bool posSafeToShift_, negSafeToShift_;                      // Whether there are multiple clusters
    MagnitudeType largestSafeShift_;                            // The largest shift considered to be safe - generally the biggest converged eigenvalue
    int NEV_;                                                   // Number of eigenvalues we seek - used in computation of trace
    std::vector<ScalarType> ritzShifts_;                        // Current Ritz shifts

    // This is only used if we're using the Schur complement solver
    RCP<MV> Z_;

    // User provided TraceMin parameters
    enum WhenToShiftType whenToShift_;                          // What triggers a Ritz shift
    enum HowToShiftType howToShift_;                            // Strategy for choosing the Ritz shifts
    bool useMultipleShifts_;                                    // Whether to use one Ritz shift or many
    bool considerClusters_;                                     // Don't shift if there is only one cluster
    bool projectAllVecs_;                                       // Use all vectors in projected minres, or just 1
    bool projectLockedVecs_;                                    // Project locked vectors too in minres; does nothing if projectAllVecs = false
    bool computeAllRes_;                                        // Compute all residuals, or just blocksize ones - helps make Ritz shifts safer
    bool useRHSR_;                                              // Use -R as the RHS of projected minres rather than AX
    bool useHarmonic_;
    MagnitudeType traceThresh_;
    MagnitudeType alpha_;

    // Variables that are only used if we're shifting when the residual gets small
    // TODO: These are currently unused
    std::string shiftNorm_;                                     // Measure 2-norm or M-norm of residual
    MagnitudeType shiftThresh_;                                 // How small must the residual be?
    bool useRelShiftThresh_;                                    // Are we scaling the shift threshold by the eigenvalues?

    // TraceMin specific functions
    // Returns the trace of KK = X'KX
    ScalarType getTrace() const;
    // Returns true if the change in trace is very small (or was very small at one point)
    // TODO: Check whether I want to return true if the trace WAS small
    bool traceLeveled();
    // Get the residuals of each cluster of eigenvalues
    // TODO: Figure out whether I want to use these for all eigenvalues or just the first
    std::vector<ScalarType> getClusterResids();
    // Computes the Ritz shifts, which is not a trivial task
    // TODO: Make it safer for indefinite matrices maybe?
    void computeRitzShifts(const std::vector<ScalarType>& clusterResids);
    // Compute the tolerance for the inner solves
    // TODO: Come back to this and make sure it does what I want it to
    std::vector<ScalarType> computeTol();
    // Solves a saddle point problem
    void solveSaddlePointProblem(RCP<MV> Delta);
    // Solves a saddle point problem by using unpreconditioned projected minres
    void solveSaddleProj(RCP<MV> Delta) const;
    // Solves a saddle point problem by using preconditioned projected...Krylov...something
    // TODO: Fix this.  Replace minres with gmres or something.
    void solveSaddleProjPrec(RCP<MV> Delta) const;
    // Solves a saddle point problem by explicitly forming the inexact Schur complement
    void solveSaddleSchur (RCP<MV> Delta) const;
    // Solves a saddle point problem with a block diagonal preconditioner
    void solveSaddleBDPrec (RCP<MV> Delta) const;
    // Solves a saddle point problem with a Hermitian/non-Hermitian splitting preconditioner
    void solveSaddleHSSPrec (RCP<MV> Delta) const;
    // Computes KK = X'KX
    void computeKK();
    // Computes the eigenpairs of KK
    void computeRitzPairs();
    // Computes X = V evecs
    void computeX();
    // Updates KX and MX without a matvec by MAGIC (or by multiplying KV and MV by evecs)
    void updateKXMX();
    // Updates the residual
    void updateResidual();
    // Adds Delta to the basis
    // TraceMin and TraceMinDavidson each do this differently, which is why this is pure virtual
    virtual void addToBasis(const RCP<const MV> Delta) =0;
    virtual void harmonicAddToBasis(const RCP<const MV> Delta) =0;
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Implementations
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  // TODO: Allow the users to supply their own linear solver
  // TODO: Add additional checking for logic errors (like trying to use gmres with multiple shifts)
  template <class ScalarType, class MV, class OP>
  TraceMinBase<ScalarType,MV,OP>::TraceMinBase(
        const RCP<Eigenproblem<ScalarType,MV,OP> >    &problem,
        const RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
        const RCP<OutputManager<ScalarType> >         &printer,
        const RCP<StatusTest<ScalarType,MV,OP> >      &tester,
        const RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
        Teuchos::ParameterList &params
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
    // timers, counters
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    timerOp_(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBase::Operation Op*x")),
    timerMOp_(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBase::Operation M*x")),
    timerSaddle_(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBase::Solving saddle point problem")),
    timerSortEval_(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBase::Sorting eigenvalues")),
    timerDS_(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBase::Direct solve")),
    timerLocal_(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBase::Local update")),
    timerCompRes_(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBase::Computing residuals")),
    timerOrtho_(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBase::Orthogonalization")),
    timerInit_(Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinBase::Initialization")),
#endif
    count_ApplyOp_(0),
    count_ApplyM_(0),
    // internal data
    blockSize_(0),
    numBlocks_(0),
    initialized_(false),
    curDim_(0),
    auxVecs_( Teuchos::Array<RCP<const MV> >(0) ),
    MauxVecs_( Teuchos::Array<RCP<const MV> >(0) ),
    numAuxVecs_(0),
    iter_(0),
    Rnorms_current_(false),
    R2norms_current_(false),
    // TraceMin variables
    previouslyLeveled_(false),
    previousTrace_(ZERO),
    posSafeToShift_(false),
    negSafeToShift_(false),
    largestSafeShift_(ZERO),
    Z_(Teuchos::null)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(problem_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::TraceMinBase::constructor: user passed null problem pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(sm_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::TraceMinBase::constructor: user passed null sort manager pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(om_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::TraceMinBase::constructor: user passed null output manager pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(tester_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::TraceMinBase::constructor: user passed null status test pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(orthman_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::TraceMinBase::constructor: user passed null orthogonalization manager pointer.");
    TEUCHOS_TEST_FOR_EXCEPTION(problem_->isHermitian() == false, std::invalid_argument,
                       "Anasazi::TraceMinBase::constructor: problem is not hermitian.");

    // get the problem operators
    Op_   = problem_->getOperator();
    MOp_  = problem_->getM();
    Prec_ = problem_->getPrec();
    hasM_ = (MOp_ != Teuchos::null);

    // Set the saddle point solver parameters
    saddleSolType_ = params.get("Saddle Solver Type", PROJECTED_KRYLOV_SOLVER);
    TEUCHOS_TEST_FOR_EXCEPTION(saddleSolType_ != PROJECTED_KRYLOV_SOLVER && saddleSolType_ != SCHUR_COMPLEMENT_SOLVER && saddleSolType_ != BD_PREC_MINRES && saddleSolType_ != HSS_PREC_GMRES, std::invalid_argument,
           "Anasazi::TraceMin::constructor: Invalid value for \"Saddle Solver Type\"; valid options are PROJECTED_KRYLOV_SOLVER, SCHUR_COMPLEMENT_SOLVER, and BD_PREC_MINRES.");

    // Set the Ritz shift parameters
    whenToShift_ = params.get("When To Shift", ALWAYS_SHIFT);
    TEUCHOS_TEST_FOR_EXCEPTION(whenToShift_ != NEVER_SHIFT && whenToShift_ != SHIFT_WHEN_TRACE_LEVELS && whenToShift_ != SHIFT_WHEN_RESID_SMALL && whenToShift_ != ALWAYS_SHIFT, std::invalid_argument,
           "Anasazi::TraceMin::constructor: Invalid value for \"When To Shift\"; valid options are \"NEVER_SHIFT\", \"SHIFT_WHEN_TRACE_LEVELS\", \"SHIFT_WHEN_RESID_SMALL\", and \"ALWAYS_SHIFT\".");

    traceThresh_ = params.get("Trace Threshold", 2e-2);
    TEUCHOS_TEST_FOR_EXCEPTION(traceThresh_ < 0, std::invalid_argument,
           "Anasazi::TraceMin::constructor: Invalid value for \"Trace Threshold\"; Must be positive.");

    howToShift_ = params.get("How To Choose Shift", ADJUSTED_RITZ_SHIFT);
    TEUCHOS_TEST_FOR_EXCEPTION(howToShift_ != LARGEST_CONVERGED_SHIFT && howToShift_ != ADJUSTED_RITZ_SHIFT && howToShift_ != RITZ_VALUES_SHIFT && howToShift_ != EXPERIMENTAL_SHIFT, std::invalid_argument,
           "Anasazi::TraceMin::constructor: Invalid value for \"How To Choose Shift\"; valid options are \"LARGEST_CONVERGED_SHIFT\", \"ADJUSTED_RITZ_SHIFT\", \"RITZ_VALUES_SHIFT\".");

    useMultipleShifts_ = params.get("Use Multiple Shifts", true);

    shiftThresh_ = params.get("Shift Threshold", 1e-2);
    useRelShiftThresh_ = params.get("Relative Shift Threshold", true);
    shiftNorm_ = params.get("Shift Norm", "2");
    TEUCHOS_TEST_FOR_EXCEPTION(shiftNorm_ != "2" && shiftNorm_ != "M", std::invalid_argument,
           "Anasazi::TraceMin::constructor: Invalid value for \"Shift Norm\"; valid options are \"2\", \"M\".");

    considerClusters_ = params.get("Consider Clusters", true);

    projectAllVecs_ = params.get("Project All Vectors", true);
    projectLockedVecs_ = params.get("Project Locked Vectors", true);
    useRHSR_ = params.get("Use Residual as RHS", true);
    useHarmonic_ = params.get("Use Harmonic Ritz Values", false);
    computeAllRes_ = params.get("Compute All Residuals", true);

    // set the block size and allocate data
    int bs = params.get("Block Size", problem_->getNEV());
    int nb = params.get("Num Blocks", 1);
    setSize(bs,nb);

    NEV_ = problem_->getNEV();

    // Create the Ritz shift operator
    ritzOp_ = rcp (new tracemin_ritz_op_type (Op_, MOp_, Prec_));

    // Set the maximum number of inner iterations
    const int innerMaxIts = params.get ("Maximum Krylov Iterations", 200);
    ritzOp_->setMaxIts (innerMaxIts);

    alpha_ = params.get ("HSS: alpha", ONE);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Destructor
  template <class ScalarType, class MV, class OP>
  TraceMinBase<ScalarType,MV,OP>::~TraceMinBase() {}


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size
  // This simply calls setSize(), modifying the block size while retaining the number of blocks.
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::setBlockSize (int blockSize)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(blockSize < 1, std::invalid_argument, "Anasazi::TraceMinBase::setSize(blocksize,numblocks): blocksize must be strictly positive.");
    setSize(blockSize,numBlocks_);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Return the current auxiliary vectors
  template <class ScalarType, class MV, class OP>
  Teuchos::Array<RCP<const MV> > TraceMinBase<ScalarType,MV,OP>::getAuxVecs() const {
    return auxVecs_;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return the current block size
  template <class ScalarType, class MV, class OP>
  int TraceMinBase<ScalarType,MV,OP>::getBlockSize() const {
    return(blockSize_);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return eigenproblem
  template <class ScalarType, class MV, class OP>
  const Eigenproblem<ScalarType,MV,OP>& TraceMinBase<ScalarType,MV,OP>::getProblem() const {
    return(*problem_);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return max subspace dim
  template <class ScalarType, class MV, class OP>
  int TraceMinBase<ScalarType,MV,OP>::getMaxSubspaceDim() const {
    return blockSize_*numBlocks_;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return current subspace dim
  template <class ScalarType, class MV, class OP>
  int TraceMinBase<ScalarType,MV,OP>::getCurSubspaceDim() const {
    if (!initialized_) return 0;
    return curDim_;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return ritz residual 2-norms
  template <class ScalarType, class MV, class OP>
  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
  TraceMinBase<ScalarType,MV,OP>::getRitzRes2Norms() {
    return getRes2Norms();
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return ritz index
  template <class ScalarType, class MV, class OP>
  std::vector<int> TraceMinBase<ScalarType,MV,OP>::getRitzIndex() {
    std::vector<int> ret(curDim_,0);
    return ret;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return ritz values
  template <class ScalarType, class MV, class OP>
  std::vector<Value<ScalarType> > TraceMinBase<ScalarType,MV,OP>::getRitzValues() {
    std::vector<Value<ScalarType> > ret(curDim_);
    for (int i=0; i<curDim_; ++i) {
      ret[i].realpart = theta_[i];
      ret[i].imagpart = ZERO;
    }
    return ret;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return pointer to ritz vectors
  template <class ScalarType, class MV, class OP>
  RCP<const MV> TraceMinBase<ScalarType,MV,OP>::getRitzVectors() {
    return X_;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // reset number of iterations
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::resetNumIters() {
    iter_=0;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return number of iterations
  template <class ScalarType, class MV, class OP>
  int TraceMinBase<ScalarType,MV,OP>::getNumIters() const {
    return(iter_);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // return state pointers
  template <class ScalarType, class MV, class OP>
  TraceMinBaseState<ScalarType,MV> TraceMinBase<ScalarType,MV,OP>::getState() const {
    TraceMinBaseState<ScalarType,MV> state;
    state.curDim = curDim_;
    state.V = V_;
    state.X = X_;
    if (Op_ != Teuchos::null) {
      state.KV = KV_;
      state.KX = KX_;
    }
    else {
      state.KV = Teuchos::null;
      state.KX = Teuchos::null;
    }
    if (hasM_) {
      state.MopV = MV_;
      state.MX = MX_;
    }
    else {
      state.MopV = Teuchos::null;
      state.MX = Teuchos::null;
    }
    state.R = R_;
    state.KK = KK_;
    state.RV = ritzVecs_;
    if (curDim_ > 0) {
      state.T = rcp(new std::vector<MagnitudeType>(theta_.begin(),theta_.begin()+curDim_) );
    } else {
      state.T = rcp(new std::vector<MagnitudeType>(0));
    }
    state.ritzShifts = rcp(new std::vector<MagnitudeType>(ritzShifts_.begin(),ritzShifts_.begin()+blockSize_) );
    state.isOrtho = true;
    state.largestSafeShift = largestSafeShift_;

    return state;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Perform TraceMinBase iterations until the StatusTest tells us to stop.
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::iterate ()
  {
    if(useHarmonic_)
    {
      harmonicIterate();
      return;
    }

    //
    // Initialize solver state
    if (initialized_ == false) {
      initialize();
    }

    // as a data member, this would be redundant and require synchronization with
    // blockSize_ and numBlocks_; we'll use a constant here.
    const int searchDim = blockSize_*numBlocks_;

    // Update obtained from solving saddle point problem
    RCP<MV> Delta = MVT::Clone(*X_,blockSize_);

    ////////////////////////////////////////////////////////////////
    // iterate until the status test tells us to stop.
    // also break if our basis is full
    while (tester_->checkStatus(this) != Passed && (numBlocks_ == 1 || curDim_ < searchDim)) {

      // Print information on current iteration
      if (om_->isVerbosity(Debug)) {
        currentStatus( om_->stream(Debug) );
      }
      else if (om_->isVerbosity(IterationDetails)) {
        currentStatus( om_->stream(IterationDetails) );
      }

      ++iter_;

      // Solve the saddle point problem
      solveSaddlePointProblem(Delta);

      // Insert Delta at the end of V
      addToBasis(Delta);

      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkV = true;
        om_->print( Debug, accuracyCheck(chk, ": after addToBasis(Delta)") );
      }

      // Compute KK := V'KV
      computeKK();

      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkKK = true;
        om_->print( Debug, accuracyCheck(chk, ": after computeKK()") );
      }

      // Compute the Ritz pairs
      computeRitzPairs();

      // Compute X := V RitzVecs
      computeX();

      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkX = true;
        om_->print( Debug, accuracyCheck(chk, ": after computeX()") );
      }

      // Compute KX := KV RitzVecs and MX := MV RitzVecs (if necessary)
      updateKXMX();

      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkKX = true;
        chk.checkMX = true;
        om_->print( Debug, accuracyCheck(chk, ": after updateKXMX()") );
      }

      // Update the residual vectors
      updateResidual();
    } // end while (statusTest == false)

  } // end of iterate()



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Perform TraceMinBase iterations until the StatusTest tells us to stop.
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::harmonicIterate ()
  {
    //
    // Initialize solver state
    if (initialized_ == false) {
      initialize();
    }

    // as a data member, this would be redundant and require synchronization with
    // blockSize_ and numBlocks_; we'll use a constant here.
    const int searchDim = blockSize_*numBlocks_;

    // Update obtained from solving saddle point problem
    RCP<MV> Delta = MVT::Clone(*X_,blockSize_);

    ////////////////////////////////////////////////////////////////
    // iterate until the status test tells us to stop.
    // also break if our basis is full
    while (tester_->checkStatus(this) != Passed && (numBlocks_ == 1 || curDim_ < searchDim)) {

      // Print information on current iteration
      if (om_->isVerbosity(Debug)) {
        currentStatus( om_->stream(Debug) );
      }
      else if (om_->isVerbosity(IterationDetails)) {
        currentStatus( om_->stream(IterationDetails) );
      }

      ++iter_;

      // Solve the saddle point problem
      solveSaddlePointProblem(Delta);

      // Insert Delta at the end of V
      harmonicAddToBasis(Delta);

      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkV = true;
        om_->print( Debug, accuracyCheck(chk, ": after addToBasis(Delta)") );
      }

      // Compute KK := V'KV
      computeKK();

      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkKK = true;
        om_->print( Debug, accuracyCheck(chk, ": after computeKK()") );
      }

      // Compute the Ritz pairs
      computeRitzPairs();

      // Compute X := V RitzVecs
      computeX();

      // Get norm of each vector in X
      int nvecs;
      if(computeAllRes_)
        nvecs = curDim_;
      else
        nvecs = blockSize_;
      std::vector<int> dimind(nvecs);
      for(int i=0; i<nvecs; i++)
        dimind[i] = i;
      RCP<MV> lclX = MVT::CloneViewNonConst(*X_,dimind);
      std::vector<ScalarType> normvec(nvecs);
      orthman_->normMat(*lclX,normvec);

      // Scale X
      for(int i=0; i<nvecs; i++)
        normvec[i] = ONE/normvec[i];
      MVT::MvScale(*lclX,normvec);

      // Scale eigenvalues
      for(int i=0; i<nvecs; i++)
      {
        theta_[i] = theta_[i] * normvec[i] * normvec[i];
      }

      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkX = true;
        om_->print( Debug, accuracyCheck(chk, ": after computeX()") );
      }

      // Compute KX := KV RitzVecs and MX := MV RitzVecs (if necessary)
      updateKXMX();

      // Scale KX and MX
      if(Op_ != Teuchos::null)
      {
        RCP<MV> lclKX = MVT::CloneViewNonConst(*KX_,dimind);
        MVT::MvScale(*lclKX,normvec);
      }
      if(hasM_)
      {
        RCP<MV> lclMX = MVT::CloneViewNonConst(*MX_,dimind);
        MVT::MvScale(*lclMX,normvec);
      }

      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkKX = true;
        chk.checkMX = true;
        om_->print( Debug, accuracyCheck(chk, ": after updateKXMX()") );
      }

      // Update the residual vectors
      updateResidual();
    } // end while (statusTest == false)

  } // end of harmonicIterate()


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // initialize the solver with default state
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::initialize()
  {
    TraceMinBaseState<ScalarType,MV> empty;
    initialize(empty);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  /* Initialize the state of the solver
   *
   * POST-CONDITIONS:
   *
   * V_ is orthonormal, orthogonal to auxVecs_, for first curDim_ vectors
   * theta_ contains Ritz w.r.t. V_(1:curDim_)
   * X is Ritz vectors w.r.t. V_(1:curDim_)
   * KX = Op*X
   * MX = M*X if hasM_
   * R = KX - MX*diag(theta_)
   *
   */
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::initialize(TraceMinBaseState<ScalarType,MV>& newstate)
  {
    // TODO: Check for bad input
    // NOTE: memory has been allocated by setBlockSize(). Use setBlock below; do not Clone
    // NOTE: Overall time spent in this routine is counted to timerInit_; portions will also be counted towards other primitives

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor inittimer( *timerInit_ );
#endif

    previouslyLeveled_ = false;

    if(useHarmonic_)
    {
      harmonicInitialize(newstate);
      return;
    }

    std::vector<int> bsind(blockSize_);
    for (int i=0; i<blockSize_; ++i) bsind[i] = i;

    // in TraceMinBase, V is primary
    // the order of dependence follows like so.
    // --init->               V,KK
    //    --ritz analysis->   theta,X
    //       --op apply->     KX,MX
    //          --compute->   R
    //
    // if the user specifies all data for a level, we will accept it.
    // otherwise, we will generate the whole level, and all subsequent levels.
    //
    // the data members are ordered based on dependence, and the levels are
    // partitioned according to the amount of work required to produce the
    // items in a level.
    //
    // inconsistent multivectors widths and lengths will not be tolerated, and
    // will be treated with exceptions.
    //
    // for multivector pointers in newstate which point directly (as opposed to indirectly, via a view) to
    // multivectors in the solver, the copy will not be affected.

    // set up V and KK: get them from newstate if user specified them
    // otherwise, set them manually
    RCP<MV> lclV, lclKV, lclMV;
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > lclKK, lclRV;

    // If V was supplied by the user...
    if (newstate.V != Teuchos::null) {
      om_->stream(Debug) << "Copying V from the user\n";

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.V) != MVT::GetGlobalLength(*V_), std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Vector length of V not correct." );
      TEUCHOS_TEST_FOR_EXCEPTION( newstate.curDim < blockSize_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Rank of new state must be at least blockSize().");
      TEUCHOS_TEST_FOR_EXCEPTION( newstate.curDim > blockSize_*numBlocks_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Rank of new state must be less than getMaxSubspaceDim().");

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.V) < newstate.curDim, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Multivector for basis in new state must be as large as specified state rank.");

      curDim_ = newstate.curDim;
      // pick an integral amount
      curDim_ = (int)(curDim_ / blockSize_)*blockSize_;

      TEUCHOS_TEST_FOR_EXCEPTION( curDim_ != newstate.curDim, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Rank of new state must be a multiple of getBlockSize().");

      // put data in V
      std::vector<int> nevind(curDim_);
      for (int i=0; i<curDim_; ++i) nevind[i] = i;
      if (newstate.V != V_) {
        if (curDim_ < MVT::GetNumberVecs(*newstate.V)) {
          newstate.V = MVT::CloneView(*newstate.V,nevind);
        }
        MVT::SetBlock(*newstate.V,nevind,*V_);
      }
      lclV = MVT::CloneViewNonConst(*V_,nevind);
    } // end if user specified V
    else
    {
      // get vectors from problem or generate something, projectAndNormalize
      RCP<const MV> ivec = problem_->getInitVec();
      TEUCHOS_TEST_FOR_EXCEPTION(ivec == Teuchos::null,std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Eigenproblem did not specify initial vectors to clone from.");
      // clear newstate so we won't use any data from it below
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;
      newstate.KK      = Teuchos::null;
      newstate.KV      = Teuchos::null;
      newstate.MopV    = Teuchos::null;
      newstate.isOrtho = false;

      // If the user did not specify a current dimension, use that of the initial vectors they provided
      if(newstate.curDim > 0)
        curDim_ = newstate.curDim;
      else
        curDim_ = MVT::GetNumberVecs(*ivec);

      // pick the largest multiple of blockSize_
      curDim_ = (int)(curDim_ / blockSize_)*blockSize_;
      if (curDim_ > blockSize_*numBlocks_) {
        // user specified too many vectors... truncate
        // this produces a full subspace, but that is okay
        curDim_ = blockSize_*numBlocks_;
      }
      bool userand = false;
      if (curDim_ == 0) {
        // we need at least blockSize_ vectors
        // use a random multivec: ignore everything from InitVec
        userand = true;
        curDim_ = blockSize_;
      }

      std::vector<int> nevind(curDim_);
      for (int i=0; i<curDim_; ++i) nevind[i] = i;

      // get a pointer into V
      // lclV has curDim vectors
      //
      // get pointer to first curDim vectors in V_
      lclV = MVT::CloneViewNonConst(*V_,nevind);
      if (userand)
      {
        // generate random vector data
        MVT::MvRandom(*lclV);
      }
      else
      {
        if(newstate.curDim > 0)
        {
          if(MVT::GetNumberVecs(*newstate.V) > curDim_) {
            RCP<const MV> helperMV = MVT::CloneView(*newstate.V,nevind);
            MVT::SetBlock(*helperMV,nevind,*lclV);
          }
          // assign ivec to first part of lclV
          MVT::SetBlock(*newstate.V,nevind,*lclV);
        }
        else
        {
          if(MVT::GetNumberVecs(*ivec) > curDim_) {
            ivec = MVT::CloneView(*ivec,nevind);
          }
          // assign ivec to first part of lclV
          MVT::SetBlock(*ivec,nevind,*lclV);
        }
      }
    } // end if user did not specify V

    // A vector of relevant indices
    std::vector<int> dimind(curDim_);
    for (int i=0; i<curDim_; ++i) dimind[i] = i;

    // Compute MV if necessary
    if(hasM_ && newstate.MopV == Teuchos::null)
    {
      om_->stream(Debug) << "Computing MV\n";

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerMOp_ );
#endif
      count_ApplyM_+= curDim_;

      newstate.isOrtho = false;
      // Get a pointer to the relevant parts of MV
      lclMV = MVT::CloneViewNonConst(*MV_,dimind);
      OPT::Apply(*MOp_,*lclV,*lclMV);
    }
    // Copy MV if necessary
    else if(hasM_)
    {
      om_->stream(Debug) << "Copying MV\n";

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.MopV) != MVT::GetGlobalLength(*MV_), std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Vector length of MV not correct." );

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.MopV) < curDim_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Number of vectors in MV not correct.");

      if(newstate.MopV != MV_) {
        if(curDim_ < MVT::GetNumberVecs(*newstate.MopV)) {
          newstate.MopV = MVT::CloneView(*newstate.MopV,dimind);
        }
        MVT::SetBlock(*newstate.MopV,dimind,*MV_);
      }
      lclMV = MVT::CloneViewNonConst(*MV_,dimind);
    }
    // There is no M, so there is no MV
    else
    {
      om_->stream(Debug) << "There is no MV\n";

      lclMV = lclV;
      MV_ = V_;
    }

    // Project and normalize so that V'V = I and Q'V=0
    if(newstate.isOrtho == false)
    {
      om_->stream(Debug) << "Project and normalize\n";

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
#endif

      // These things are now invalid
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;
      newstate.KK      = Teuchos::null;
      newstate.KV      = Teuchos::null;

      int rank;
      if(auxVecs_.size() > 0)
      {
        rank = orthman_->projectAndNormalizeMat(*lclV, auxVecs_,
               Teuchos::tuple(RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
               Teuchos::null, lclMV, MauxVecs_);
      }
      else
      {
        rank = orthman_->normalizeMat(*lclV,Teuchos::null,lclMV);
      }

      TEUCHOS_TEST_FOR_EXCEPTION(rank != curDim_,TraceMinBaseInitFailure,
             "Anasazi::TraceMinBase::initialize(): Couldn't generate initial basis of full rank.");
    }

    // Compute KV
    if(Op_ != Teuchos::null && newstate.KV == Teuchos::null)
    {
      om_->stream(Debug) << "Computing KV\n";

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerOp_ );
#endif
      count_ApplyOp_+= curDim_;

      // These things are now invalid
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;
      newstate.KK      = Teuchos::null;
      newstate.KV      = Teuchos::null;

      lclKV = MVT::CloneViewNonConst(*KV_,dimind);
      OPT::Apply(*Op_,*lclV,*lclKV);
    }
    // Copy KV
    else if(Op_ != Teuchos::null)
    {
      om_->stream(Debug) << "Copying MV\n";

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.KV) != MVT::GetGlobalLength(*KV_), std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Vector length of MV not correct." );

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.KV) < curDim_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Number of vectors in KV not correct.");

      if (newstate.KV != KV_) {
        if (curDim_ < MVT::GetNumberVecs(*newstate.KV)) {
          newstate.KV = MVT::CloneView(*newstate.KV,dimind);
        }
        MVT::SetBlock(*newstate.KV,dimind,*KV_);
      }
      lclKV = MVT::CloneViewNonConst(*KV_,dimind);
    }
    else
    {
      om_->stream(Debug) << "There is no KV\n";

      lclKV = lclV;
      KV_ = V_;
    }

    // Compute KK
    if(newstate.KK == Teuchos::null)
    {
      om_->stream(Debug) << "Computing KK\n";

      // These things are now invalid
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;

      // Note: there used to be a bug here.
      // We can't just call computeKK because it only computes the new parts of KK

      // Get a pointer to the part of KK we're interested in
      lclKK = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*KK_,curDim_,curDim_) );

      // KK := V'KV
      MVT::MvTransMv(ONE,*lclV,*lclKV,*lclKK);
    }
    // Copy KK
    else
    {
      om_->stream(Debug) << "Copying KK\n";

      // check size of KK
      TEUCHOS_TEST_FOR_EXCEPTION( newstate.KK->numRows() < curDim_ || newstate.KK->numCols() < curDim_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Projected matrix in new state must be as large as specified state rank.");

      // put data into KK_
      lclKK = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*KK_,curDim_,curDim_) );
      if (newstate.KK != KK_) {
        if (newstate.KK->numRows() > curDim_ || newstate.KK->numCols() > curDim_) {
          newstate.KK = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*newstate.KK,curDim_,curDim_) );
        }
        lclKK->assign(*newstate.KK);
      }
    }

    // Compute Ritz pairs
    if(newstate.T == Teuchos::null || newstate.RV == Teuchos::null)
    {
      om_->stream(Debug) << "Computing Ritz pairs\n";

      // These things are now invalid
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;

      computeRitzPairs();
    }
    // Copy Ritz pairs
    else
    {
      om_->stream(Debug) << "Copying Ritz pairs\n";

      TEUCHOS_TEST_FOR_EXCEPTION((signed int)(newstate.T->size()) != curDim_,
             std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of T must be consistent with dimension of V.");

      TEUCHOS_TEST_FOR_EXCEPTION( newstate.RV->numRows() < curDim_ || newstate.RV->numCols() < curDim_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Ritz vectors in new state must be as large as specified state rank.");

      std::copy(newstate.T->begin(),newstate.T->end(),theta_.begin());

      lclRV = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*ritzVecs_,curDim_,curDim_) );
      if (newstate.RV != ritzVecs_) {
        if (newstate.RV->numRows() > curDim_ || newstate.RV->numCols() > curDim_) {
          newstate.RV = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*newstate.RV,curDim_,curDim_) );
        }
        lclRV->assign(*newstate.RV);
      }
    }

    // Compute X
    if(newstate.X == Teuchos::null)
    {
      om_->stream(Debug) << "Computing X\n";

      // These things are now invalid
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;

      computeX();
    }
    // Copy X
    else
    {
      om_->stream(Debug) << "Copying X\n";

      if(computeAllRes_ == false)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.X) != blockSize_ || MVT::GetGlobalLength(*newstate.X) != MVT::GetGlobalLength(*X_),
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of X must be consistent with block size and length of V.");

        if(newstate.X != X_) {
          MVT::SetBlock(*newstate.X,bsind,*X_);
        }
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.X) != curDim_ || MVT::GetGlobalLength(*newstate.X) != MVT::GetGlobalLength(*X_),
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of X must be consistent with current dimension and length of V.");

        if(newstate.X != X_) {
          MVT::SetBlock(*newstate.X,dimind,*X_);
        }
      }
    }

    // Compute KX and MX if necessary
    // TODO: These technically should be separate; it won't matter much in terms of running time though
    if((Op_ != Teuchos::null && newstate.KX == Teuchos::null) || (hasM_ && newstate.MX == Teuchos::null))
    {
      om_->stream(Debug) << "Computing KX and MX\n";

      // These things are now invalid
      newstate.R = Teuchos::null;

      updateKXMX();
    }
    // Copy KX and MX if necessary
    else
    {
      om_->stream(Debug) << "Copying KX and MX\n";
      if(Op_ != Teuchos::null)
      {
        if(computeAllRes_ == false)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.KX) != blockSize_ || MVT::GetGlobalLength(*newstate.KX) != MVT::GetGlobalLength(*KX_),
                 std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of KX must be consistent with block size and length of V.");

          if(newstate.KX != KX_) {
            MVT::SetBlock(*newstate.KX,bsind,*KX_);
          }
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.KX) != curDim_ || MVT::GetGlobalLength(*newstate.KX) != MVT::GetGlobalLength(*KX_),
                 std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of KX must be consistent with current dimension and length of V.");

          if (newstate.KX != KX_) {
            MVT::SetBlock(*newstate.KX,dimind,*KX_);
          }
        }
      }

      if(hasM_)
      {
        if(computeAllRes_ == false)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.MX) != blockSize_ || MVT::GetGlobalLength(*newstate.MX) != MVT::GetGlobalLength(*MX_),
                 std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of MX must be consistent with block size and length of V.");

          if (newstate.MX != MX_) {
            MVT::SetBlock(*newstate.MX,bsind,*MX_);
          }
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.MX) != curDim_ || MVT::GetGlobalLength(*newstate.MX) != MVT::GetGlobalLength(*MX_),
                 std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of MX must be consistent with current dimension and length of V.");

          if (newstate.MX != MX_) {
            MVT::SetBlock(*newstate.MX,dimind,*MX_);
          }
        }
      }
    }

    // Compute R
    if(newstate.R == Teuchos::null)
    {
      om_->stream(Debug) << "Computing R\n";

      updateResidual();
    }
    // Copy R
    else
    {
      om_->stream(Debug) << "Copying R\n";

      if(computeAllRes_ == false)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*X_),
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): vector length of newstate.R not correct." );
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.R) != blockSize_,
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): newstate.R must have at least block size vectors." );
        if (newstate.R != R_) {
          MVT::SetBlock(*newstate.R,bsind,*R_);
        }
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*X_),
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): vector length of newstate.R not correct." );
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.R) != curDim_,
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): newstate.R must have at least curDim vectors." );
        if (newstate.R != R_) {
          MVT::SetBlock(*newstate.R,dimind,*R_);
        }
      }
    }

    // R has been updated; mark the norms as out-of-date
    Rnorms_current_ = false;
    R2norms_current_ = false;

    // Set the largest safe shift
    largestSafeShift_ = newstate.largestSafeShift;

    // Copy over the Ritz shifts
    if(newstate.ritzShifts != Teuchos::null)
    {
      om_->stream(Debug) << "Copying Ritz shifts\n";
      std::copy(newstate.ritzShifts->begin(),newstate.ritzShifts->end(),ritzShifts_.begin());
    }
    else
    {
      om_->stream(Debug) << "Setting Ritz shifts to 0\n";
      for(size_t i=0; i<ritzShifts_.size(); i++)
        ritzShifts_[i] = ZERO;
    }

    for(size_t i=0; i<ritzShifts_.size(); i++)
      om_->stream(Debug) << "Ritz shifts[" << i << "] = " << ritzShifts_[i] << std::endl;

    // finally, we are initialized
    initialized_ = true;

    if (om_->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkV = true;
      chk.checkX = true;
      chk.checkKX = true;
      chk.checkMX = true;
      chk.checkQ = true;
      chk.checkKK = true;
      om_->print( Debug, accuracyCheck(chk, ": after initialize()") );
    }

    // Print information on current status
    if (om_->isVerbosity(Debug)) {
      currentStatus( om_->stream(Debug) );
    }
    else if (om_->isVerbosity(IterationDetails)) {
      currentStatus( om_->stream(IterationDetails) );
    }
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  /* Initialize the state of the solver
   *
   * POST-CONDITIONS:
   *
   * V_ is orthonormal, orthogonal to auxVecs_, for first curDim_ vectors
   * theta_ contains Ritz w.r.t. V_(1:curDim_)
   * X is Ritz vectors w.r.t. V_(1:curDim_)
   * KX = Op*X
   * MX = M*X if hasM_
   * R = KX - MX*diag(theta_)
   *
   */
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::harmonicInitialize(TraceMinBaseState<ScalarType,MV> newstate)
  {
    // TODO: Check for bad input
    // NOTE: memory has been allocated by setBlockSize(). Use setBlock below; do not Clone
    // NOTE: Overall time spent in this routine is counted to timerInit_; portions will also be counted towards other primitives

    std::vector<int> bsind(blockSize_);
    for (int i=0; i<blockSize_; ++i) bsind[i] = i;

    // in TraceMinBase, V is primary
    // the order of dependence follows like so.
    // --init->               V,KK
    //    --ritz analysis->   theta,X
    //       --op apply->     KX,MX
    //          --compute->   R
    //
    // if the user specifies all data for a level, we will accept it.
    // otherwise, we will generate the whole level, and all subsequent levels.
    //
    // the data members are ordered based on dependence, and the levels are
    // partitioned according to the amount of work required to produce the
    // items in a level.
    //
    // inconsistent multivectors widths and lengths will not be tolerated, and
    // will be treated with exceptions.
    //
    // for multivector pointers in newstate which point directly (as opposed to indirectly, via a view) to
    // multivectors in the solver, the copy will not be affected.

    // set up V and KK: get them from newstate if user specified them
    // otherwise, set them manually
    RCP<MV> lclV, lclKV, lclMV;
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > lclKK, lclRV;

    // If V was supplied by the user...
    if (newstate.V != Teuchos::null) {
      om_->stream(Debug) << "Copying V from the user\n";

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.V) != MVT::GetGlobalLength(*V_), std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Vector length of V not correct." );
      TEUCHOS_TEST_FOR_EXCEPTION( newstate.curDim < blockSize_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Rank of new state must be at least blockSize().");
      TEUCHOS_TEST_FOR_EXCEPTION( newstate.curDim > blockSize_*numBlocks_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Rank of new state must be less than getMaxSubspaceDim().");

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.V) < newstate.curDim, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Multivector for basis in new state must be as large as specified state rank.");

      curDim_ = newstate.curDim;
      // pick an integral amount
      curDim_ = (int)(curDim_ / blockSize_)*blockSize_;

      TEUCHOS_TEST_FOR_EXCEPTION( curDim_ != newstate.curDim, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Rank of new state must be a multiple of getBlockSize().");

      // put data in V
      std::vector<int> nevind(curDim_);
      for (int i=0; i<curDim_; ++i) nevind[i] = i;
      if (newstate.V != V_) {
        if (curDim_ < MVT::GetNumberVecs(*newstate.V)) {
          newstate.V = MVT::CloneView(*newstate.V,nevind);
        }
        MVT::SetBlock(*newstate.V,nevind,*V_);
      }
      lclV = MVT::CloneViewNonConst(*V_,nevind);
    } // end if user specified V
    else
    {
      // get vectors from problem or generate something, projectAndNormalize
      RCP<const MV> ivec = problem_->getInitVec();
      TEUCHOS_TEST_FOR_EXCEPTION(ivec == Teuchos::null,std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Eigenproblem did not specify initial vectors to clone from.");
      // clear newstate so we won't use any data from it below
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;
      newstate.KK      = Teuchos::null;
      newstate.KV      = Teuchos::null;
      newstate.MopV    = Teuchos::null;
      newstate.isOrtho = false;

      // If the user did not specify a current dimension, use that of the initial vectors they provided
      if(newstate.curDim > 0)
        curDim_ = newstate.curDim;
      else
        curDim_ = MVT::GetNumberVecs(*ivec);

      // pick the largest multiple of blockSize_
      curDim_ = (int)(curDim_ / blockSize_)*blockSize_;
      if (curDim_ > blockSize_*numBlocks_) {
        // user specified too many vectors... truncate
        // this produces a full subspace, but that is okay
        curDim_ = blockSize_*numBlocks_;
      }
      bool userand = false;
      if (curDim_ == 0) {
        // we need at least blockSize_ vectors
        // use a random multivec: ignore everything from InitVec
        userand = true;
        curDim_ = blockSize_;
      }

      std::vector<int> nevind(curDim_);
      for (int i=0; i<curDim_; ++i) nevind[i] = i;

      // get a pointer into V
      // lclV has curDim vectors
      //
      // get pointer to first curDim vectors in V_
      lclV = MVT::CloneViewNonConst(*V_,nevind);
      if (userand)
      {
        // generate random vector data
        MVT::MvRandom(*lclV);
      }
      else
      {
        if(newstate.curDim > 0)
        {
          if(MVT::GetNumberVecs(*newstate.V) > curDim_) {
            RCP<const MV> helperMV = MVT::CloneView(*newstate.V,nevind);
            MVT::SetBlock(*helperMV,nevind,*lclV);
          }
          // assign ivec to first part of lclV
          MVT::SetBlock(*newstate.V,nevind,*lclV);
        }
        else
        {
          if(MVT::GetNumberVecs(*ivec) > curDim_) {
            ivec = MVT::CloneView(*ivec,nevind);
          }
          // assign ivec to first part of lclV
          MVT::SetBlock(*ivec,nevind,*lclV);
        }
      }
    } // end if user did not specify V

    // Nuke everything from orbit
    // This is a temporary measure due to a bug in the code that I have not found yet
    // It adds a minimal amount of work
    newstate.X       = Teuchos::null;
    newstate.MX      = Teuchos::null;
    newstate.KX      = Teuchos::null;
    newstate.R       = Teuchos::null;
    newstate.T       = Teuchos::null;
    newstate.RV      = Teuchos::null;
    newstate.KK      = Teuchos::null;
    newstate.KV      = Teuchos::null;
    newstate.MopV    = Teuchos::null;
    newstate.isOrtho = false;

    // A vector of relevant indices
    std::vector<int> dimind(curDim_);
    for (int i=0; i<curDim_; ++i) dimind[i] = i;

    // Project the auxVecs out of V
    if(auxVecs_.size() > 0)
      orthman_->projectMat(*lclV,auxVecs_);

    // Compute KV
    if(Op_ != Teuchos::null && newstate.KV == Teuchos::null)
    {
      om_->stream(Debug) << "Computing KV\n";

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerOp_ );
#endif
      count_ApplyOp_+= curDim_;

      // These things are now invalid
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;
      newstate.KK      = Teuchos::null;

      lclKV = MVT::CloneViewNonConst(*KV_,dimind);
      OPT::Apply(*Op_,*lclV,*lclKV);
    }
    // Copy KV
    else if(Op_ != Teuchos::null)
    {
      om_->stream(Debug) << "Copying KV\n";

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.KV) != MVT::GetGlobalLength(*KV_), std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Vector length of KV not correct." );

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.KV) < curDim_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Number of vectors in KV not correct.");

      if (newstate.KV != KV_) {
        if (curDim_ < MVT::GetNumberVecs(*newstate.KV)) {
          newstate.KV = MVT::CloneView(*newstate.KV,dimind);
        }
        MVT::SetBlock(*newstate.KV,dimind,*KV_);
      }
      lclKV = MVT::CloneViewNonConst(*KV_,dimind);
    }
    else
    {
      om_->stream(Debug) << "There is no KV\n";

      lclKV = lclV;
      KV_ = V_;
    }



    // Project and normalize so that V'V = I and Q'V=0
    if(newstate.isOrtho == false)
    {
      om_->stream(Debug) << "Project and normalize\n";

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
#endif

      // These things are now invalid
      newstate.MopV    = Teuchos::null;
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;
      newstate.KK      = Teuchos::null;

      // Normalize lclKV
      RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > gamma = rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(curDim_,curDim_));
      int rank = orthman_->normalizeMat(*lclKV,gamma);

      // lclV = lclV/gamma
      Teuchos::SerialDenseSolver<int,ScalarType> SDsolver;
      SDsolver.setMatrix(gamma);
      SDsolver.invert();
      RCP<MV> tempMV = MVT::CloneCopy(*lclV);
      MVT::MvTimesMatAddMv(ONE,*tempMV,*gamma,ZERO,*lclV);

      TEUCHOS_TEST_FOR_EXCEPTION(rank != curDim_,TraceMinBaseInitFailure,
             "Anasazi::TraceMinBase::initialize(): Couldn't generate initial basis of full rank.");
    }

    // Compute MV if necessary
    if(hasM_ && newstate.MopV == Teuchos::null)
    {
      om_->stream(Debug) << "Computing MV\n";

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerMOp_ );
#endif
      count_ApplyM_+= curDim_;

      // Get a pointer to the relevant parts of MV
      lclMV = MVT::CloneViewNonConst(*MV_,dimind);
      OPT::Apply(*MOp_,*lclV,*lclMV);
    }
    // Copy MV if necessary
    else if(hasM_)
    {
      om_->stream(Debug) << "Copying MV\n";

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.MopV) != MVT::GetGlobalLength(*MV_), std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Vector length of MV not correct." );

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.MopV) < curDim_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Number of vectors in MV not correct.");

      if(newstate.MopV != MV_) {
        if(curDim_ < MVT::GetNumberVecs(*newstate.MopV)) {
          newstate.MopV = MVT::CloneView(*newstate.MopV,dimind);
        }
        MVT::SetBlock(*newstate.MopV,dimind,*MV_);
      }
      lclMV = MVT::CloneViewNonConst(*MV_,dimind);
    }
    // There is no M, so there is no MV
    else
    {
      om_->stream(Debug) << "There is no MV\n";

      lclMV = lclV;
      MV_ = V_;
    }

    // Compute KK
    if(newstate.KK == Teuchos::null)
    {
      om_->stream(Debug) << "Computing KK\n";

      // These things are now invalid
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;

      // Note: there used to be a bug here.
      // We can't just call computeKK because it only computes the new parts of KK

      // Get a pointer to the part of KK we're interested in
      lclKK = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*KK_,curDim_,curDim_) );

      // KK := V'KV
      MVT::MvTransMv(ONE,*lclV,*lclKV,*lclKK);
    }
    // Copy KK
    else
    {
      om_->stream(Debug) << "Copying KK\n";

      // check size of KK
      TEUCHOS_TEST_FOR_EXCEPTION( newstate.KK->numRows() < curDim_ || newstate.KK->numCols() < curDim_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Projected matrix in new state must be as large as specified state rank.");

      // put data into KK_
      lclKK = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*KK_,curDim_,curDim_) );
      if (newstate.KK != KK_) {
        if (newstate.KK->numRows() > curDim_ || newstate.KK->numCols() > curDim_) {
          newstate.KK = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*newstate.KK,curDim_,curDim_) );
        }
        lclKK->assign(*newstate.KK);
      }
    }

    // Compute Ritz pairs
    if(newstate.T == Teuchos::null || newstate.RV == Teuchos::null)
    {
      om_->stream(Debug) << "Computing Ritz pairs\n";

      // These things are now invalid
      newstate.X       = Teuchos::null;
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;
      newstate.T       = Teuchos::null;
      newstate.RV      = Teuchos::null;

      computeRitzPairs();
    }
    // Copy Ritz pairs
    else
    {
      om_->stream(Debug) << "Copying Ritz pairs\n";

      TEUCHOS_TEST_FOR_EXCEPTION((signed int)(newstate.T->size()) != curDim_,
             std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of T must be consistent with dimension of V.");

      TEUCHOS_TEST_FOR_EXCEPTION( newstate.RV->numRows() < curDim_ || newstate.RV->numCols() < curDim_, std::invalid_argument,
             "Anasazi::TraceMinBase::initialize(newstate): Ritz vectors in new state must be as large as specified state rank.");

      std::copy(newstate.T->begin(),newstate.T->end(),theta_.begin());

      lclRV = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*ritzVecs_,curDim_,curDim_) );
      if (newstate.RV != ritzVecs_) {
        if (newstate.RV->numRows() > curDim_ || newstate.RV->numCols() > curDim_) {
          newstate.RV = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*newstate.RV,curDim_,curDim_) );
        }
        lclRV->assign(*newstate.RV);
      }
    }

    // Compute X
    if(newstate.X == Teuchos::null)
    {
      om_->stream(Debug) << "Computing X\n";

      // These things are now invalid
      newstate.MX      = Teuchos::null;
      newstate.KX      = Teuchos::null;
      newstate.R       = Teuchos::null;

      computeX();
    }
    // Copy X
    else
    {
      om_->stream(Debug) << "Copying X\n";

      if(computeAllRes_ == false)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.X) != blockSize_ || MVT::GetGlobalLength(*newstate.X) != MVT::GetGlobalLength(*X_),
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of X must be consistent with block size and length of V.");

        if(newstate.X != X_) {
          MVT::SetBlock(*newstate.X,bsind,*X_);
        }
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.X) != curDim_ || MVT::GetGlobalLength(*newstate.X) != MVT::GetGlobalLength(*X_),
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of X must be consistent with current dimension and length of V.");

        if(newstate.X != X_) {
          MVT::SetBlock(*newstate.X,dimind,*X_);
        }
      }
    }

    // Compute KX and MX if necessary
    // TODO: These technically should be separate; it won't matter much in terms of running time though
    if((Op_ != Teuchos::null && newstate.KX == Teuchos::null) || (hasM_ && newstate.MX == Teuchos::null))
    {
      om_->stream(Debug) << "Computing KX and MX\n";

      // These things are now invalid
      newstate.R = Teuchos::null;

      updateKXMX();
    }
    // Copy KX and MX if necessary
    else
    {
      om_->stream(Debug) << "Copying KX and MX\n";
      if(Op_ != Teuchos::null)
      {
        if(computeAllRes_ == false)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.KX) != blockSize_ || MVT::GetGlobalLength(*newstate.KX) != MVT::GetGlobalLength(*KX_),
                 std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of KX must be consistent with block size and length of V.");

          if(newstate.KX != KX_) {
            MVT::SetBlock(*newstate.KX,bsind,*KX_);
          }
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.KX) != curDim_ || MVT::GetGlobalLength(*newstate.KX) != MVT::GetGlobalLength(*KX_),
                 std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of KX must be consistent with current dimension and length of V.");

          if (newstate.KX != KX_) {
            MVT::SetBlock(*newstate.KX,dimind,*KX_);
          }
        }
      }

      if(hasM_)
      {
        if(computeAllRes_ == false)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.MX) != blockSize_ || MVT::GetGlobalLength(*newstate.MX) != MVT::GetGlobalLength(*MX_),
                 std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of MX must be consistent with block size and length of V.");

          if (newstate.MX != MX_) {
            MVT::SetBlock(*newstate.MX,bsind,*MX_);
          }
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.MX) != curDim_ || MVT::GetGlobalLength(*newstate.MX) != MVT::GetGlobalLength(*MX_),
                 std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): Size of MX must be consistent with current dimension and length of V.");

          if (newstate.MX != MX_) {
            MVT::SetBlock(*newstate.MX,dimind,*MX_);
          }
        }
      }
    }

    // Scale X so each vector is of length 1
    {
      // Get norm of each vector in X
      const int nvecs = computeAllRes_ ? curDim_ : blockSize_;
      Teuchos::Range1D dimind2 (0, nvecs-1);
      RCP<MV> lclX = MVT::CloneViewNonConst(*X_, dimind2);
      std::vector<ScalarType> normvec(nvecs);
      orthman_->normMat(*lclX,normvec);

      // Scale X, KX, and MX accordingly
      for (int i = 0; i < nvecs; ++i) {
        normvec[i] = ONE / normvec[i];
      }
      MVT::MvScale (*lclX, normvec);
      if (Op_ != Teuchos::null) {
        RCP<MV> lclKX = MVT::CloneViewNonConst (*KX_, dimind2);
        MVT::MvScale (*lclKX, normvec);
      }
      if (hasM_) {
        RCP<MV> lclMX = MVT::CloneViewNonConst (*MX_, dimind2);
        MVT::MvScale (*lclMX, normvec);
      }

      // Scale eigenvalues
      for (int i = 0; i < nvecs; ++i) {
        theta_[i] = theta_[i] * normvec[i] * normvec[i];
      }
    }

    // Compute R
    if(newstate.R == Teuchos::null)
    {
      om_->stream(Debug) << "Computing R\n";

      updateResidual();
    }
    // Copy R
    else
    {
      om_->stream(Debug) << "Copying R\n";

      if(computeAllRes_ == false)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*X_),
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): vector length of newstate.R not correct." );
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.R) != blockSize_,
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): newstate.R must have at least block size vectors." );
        if (newstate.R != R_) {
          MVT::SetBlock(*newstate.R,bsind,*R_);
        }
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*X_),
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): vector length of newstate.R not correct." );
        TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*newstate.R) != curDim_,
               std::invalid_argument, "Anasazi::TraceMinBase::initialize(newstate): newstate.R must have at least curDim vectors." );
        if (newstate.R != R_) {
          MVT::SetBlock(*newstate.R,dimind,*R_);
        }
      }
    }

    // R has been updated; mark the norms as out-of-date
    Rnorms_current_ = false;
    R2norms_current_ = false;

    // Set the largest safe shift
    largestSafeShift_ = newstate.largestSafeShift;

    // Copy over the Ritz shifts
    if(newstate.ritzShifts != Teuchos::null)
    {
      om_->stream(Debug) << "Copying Ritz shifts\n";
      std::copy(newstate.ritzShifts->begin(),newstate.ritzShifts->end(),ritzShifts_.begin());
    }
    else
    {
      om_->stream(Debug) << "Setting Ritz shifts to 0\n";
      for(size_t i=0; i<ritzShifts_.size(); i++)
        ritzShifts_[i] = ZERO;
    }

    for(size_t i=0; i<ritzShifts_.size(); i++)
      om_->stream(Debug) << "Ritz shifts[" << i << "] = " << ritzShifts_[i] << std::endl;

    // finally, we are initialized
    initialized_ = true;

    if (om_->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkV = true;
      chk.checkX = true;
      chk.checkKX = true;
      chk.checkMX = true;
      chk.checkQ = true;
      chk.checkKK = true;
      om_->print( Debug, accuracyCheck(chk, ": after initialize()") );
    }

    // Print information on current status
    if (om_->isVerbosity(Debug)) {
      currentStatus( om_->stream(Debug) );
    }
    else if (om_->isVerbosity(IterationDetails)) {
      currentStatus( om_->stream(IterationDetails) );
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Return initialized state
  template <class ScalarType, class MV, class OP>
  bool TraceMinBase<ScalarType,MV,OP>::isInitialized() const { return initialized_; }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::setSize (int blockSize, int numBlocks)
  {
    // This routine only allocates space; it doesn't not perform any computation
    // any change in size will invalidate the state of the solver.

    TEUCHOS_TEST_FOR_EXCEPTION(blockSize < 1, std::invalid_argument, "Anasazi::TraceMinBase::setSize(blocksize,numblocks): blocksize must be strictly positive.");

    if (blockSize == blockSize_ && numBlocks == numBlocks_) {
      // do nothing
      return;
    }

    blockSize_ = blockSize;
    numBlocks_ = numBlocks;

    RCP<const MV> tmp;
    // grab some Multivector to Clone
    // in practice, getInitVec() should always provide this, but it is possible to use a
    // Eigenproblem with nothing in getInitVec() by manually initializing with initialize();
    // in case of that strange scenario, we will try to Clone from X_ first, then resort to getInitVec()
    if (X_ != Teuchos::null) { // this is equivalent to blockSize_ > 0
      tmp = X_;
    }
    else {
      tmp = problem_->getInitVec();
      TEUCHOS_TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::invalid_argument,
             "Anasazi::TraceMinBase::setSize(): eigenproblem did not specify initial vectors to clone from.");
    }

    TEUCHOS_TEST_FOR_EXCEPTION(numAuxVecs_+blockSize*static_cast<ptrdiff_t>(numBlocks) > MVT::GetGlobalLength(*tmp),std::invalid_argument,
           "Anasazi::TraceMinBase::setSize(): max subspace dimension and auxilliary subspace too large.  Potentially impossible orthogonality constraints.");

    // New subspace dimension
    int newsd = blockSize_*numBlocks_;

    //////////////////////////////////
    // blockSize dependent
    //
    ritzShifts_.resize(blockSize_,ZERO);
    if(computeAllRes_ == false)
    {
      // grow/allocate vectors
      Rnorms_.resize(blockSize_,NANVAL);
      R2norms_.resize(blockSize_,NANVAL);
      //
      // clone multivectors off of tmp
      //
      // free current allocation first, to make room for new allocation
      X_ = Teuchos::null;
      KX_ = Teuchos::null;
      MX_ = Teuchos::null;
      R_ = Teuchos::null;
      V_ = Teuchos::null;
      KV_ = Teuchos::null;
      MV_ = Teuchos::null;

      om_->print(Debug," >> Allocating X_\n");
      X_ = MVT::Clone(*tmp,blockSize_);
      if(Op_ != Teuchos::null) {
        om_->print(Debug," >> Allocating KX_\n");
        KX_ = MVT::Clone(*tmp,blockSize_);
      }
      else {
        KX_ = X_;
      }
      if (hasM_) {
        om_->print(Debug," >> Allocating MX_\n");
        MX_ = MVT::Clone(*tmp,blockSize_);
      }
      else {
        MX_ = X_;
      }
      om_->print(Debug," >> Allocating R_\n");
      R_ = MVT::Clone(*tmp,blockSize_);
    }
    else
    {
      // grow/allocate vectors
      Rnorms_.resize(newsd,NANVAL);
      R2norms_.resize(newsd,NANVAL);
      //
      // clone multivectors off of tmp
      //
      // free current allocation first, to make room for new allocation
      X_ = Teuchos::null;
      KX_ = Teuchos::null;
      MX_ = Teuchos::null;
      R_ = Teuchos::null;
      V_ = Teuchos::null;
      KV_ = Teuchos::null;
      MV_ = Teuchos::null;

      om_->print(Debug," >> Allocating X_\n");
      X_ = MVT::Clone(*tmp,newsd);
      if(Op_ != Teuchos::null) {
        om_->print(Debug," >> Allocating KX_\n");
        KX_ = MVT::Clone(*tmp,newsd);
      }
      else {
        KX_ = X_;
      }
      if (hasM_) {
        om_->print(Debug," >> Allocating MX_\n");
        MX_ = MVT::Clone(*tmp,newsd);
      }
      else {
        MX_ = X_;
      }
      om_->print(Debug," >> Allocating R_\n");
      R_ = MVT::Clone(*tmp,newsd);
    }

    //////////////////////////////////
    // blockSize*numBlocks dependent
    //
    theta_.resize(newsd,NANVAL);
    om_->print(Debug," >> Allocating V_\n");
    V_ = MVT::Clone(*tmp,newsd);
    KK_ = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );
    ritzVecs_ = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );

    if(Op_ != Teuchos::null) {
      om_->print(Debug," >> Allocating KV_\n");
      KV_ = MVT::Clone(*tmp,newsd);
    }
    else {
      KV_ = V_;
    }
    if (hasM_) {
      om_->print(Debug," >> Allocating MV_\n");
      MV_ = MVT::Clone(*tmp,newsd);
    }
    else {
      MV_ = V_;
    }

    om_->print(Debug," >> done allocating.\n");

    initialized_ = false;
    curDim_ = 0;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the auxiliary vectors
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::setAuxVecs(const Teuchos::Array<RCP<const MV> > &auxvecs) {
    typedef typename Teuchos::Array<RCP<const MV> >::iterator tarcpmv;

    // set new auxiliary vectors
    auxVecs_ = auxvecs;

    if(hasM_)
      MauxVecs_.resize(0);
    numAuxVecs_ = 0;

    for (tarcpmv i=auxVecs_.begin(); i != auxVecs_.end(); ++i) {
      numAuxVecs_ += MVT::GetNumberVecs(**i);

      if(hasM_)
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerMOp_ );
#endif
        count_ApplyM_+= MVT::GetNumberVecs(**i);

        RCP<MV> helperMV = MVT::Clone(**i,MVT::GetNumberVecs(**i));
        OPT::Apply(*MOp_,**i,*helperMV);
        MauxVecs_.push_back(helperMV);
      }
    }

    // If the solver has been initialized, V is not necessarily orthogonal to new auxiliary vectors
    if (numAuxVecs_ > 0 && initialized_) {
      initialized_ = false;
    }

    if (om_->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkQ = true;
      om_->print( Debug, accuracyCheck(chk, ": after setAuxVecs()") );
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // compute/return residual M-norms
  template <class ScalarType, class MV, class OP>
  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
  TraceMinBase<ScalarType,MV,OP>::getResNorms() {
    if (Rnorms_current_ == false) {
      // Update the residual norms
      if(computeAllRes_)
      {
        std::vector<int> curind(curDim_);
        for(int i=0; i<curDim_; i++)
          curind[i] = i;

        RCP<const MV> locR = MVT::CloneView(*R_,curind);
        std::vector<ScalarType> locNorms(curDim_);
        orthman_->norm(*locR,locNorms);

        for(int i=0; i<curDim_; i++)
          Rnorms_[i] = locNorms[i];
        for(int i=curDim_+1; i<blockSize_*numBlocks_; i++)
          Rnorms_[i] = NANVAL;

        Rnorms_current_ = true;
        locNorms.resize(blockSize_);
        return locNorms;
      }
      else
        orthman_->norm(*R_,Rnorms_);
      Rnorms_current_ = true;
    }
    else if(computeAllRes_)
    {
      std::vector<ScalarType> locNorms(blockSize_);
      for(int i=0; i<blockSize_; i++)
        locNorms[i] = Rnorms_[i];
      return locNorms;
    }

    return Rnorms_;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // compute/return residual 2-norms
  template <class ScalarType, class MV, class OP>
  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
  TraceMinBase<ScalarType,MV,OP>::getRes2Norms() {
    if (R2norms_current_ == false) {
      // Update the residual 2-norms
      if(computeAllRes_)
      {
        std::vector<int> curind(curDim_);
        for(int i=0; i<curDim_; i++)
          curind[i] = i;

        RCP<const MV> locR = MVT::CloneView(*R_,curind);
        std::vector<ScalarType> locNorms(curDim_);
        MVT::MvNorm(*locR,locNorms);

        for(int i=0; i<curDim_; i++)
        {
          R2norms_[i] = locNorms[i];
        }
        for(int i=curDim_+1; i<blockSize_*numBlocks_; i++)
          R2norms_[i] = NANVAL;

        R2norms_current_ = true;
        locNorms.resize(blockSize_);
        return locNorms;
      }
      else
        MVT::MvNorm(*R_,R2norms_);
      R2norms_current_ = true;
    }
    else if(computeAllRes_)
    {
      std::vector<ScalarType> locNorms(blockSize_);
      for(int i=0; i<blockSize_; i++)
        locNorms[i] = R2norms_[i];
      return locNorms;
    }

    return R2norms_;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set a new StatusTest for the solver.
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::setStatusTest(RCP<StatusTest<ScalarType,MV,OP> > test) {
    TEUCHOS_TEST_FOR_EXCEPTION(test == Teuchos::null,std::invalid_argument,
        "Anasazi::TraceMinBase::setStatusTest() was passed a null StatusTest.");
    tester_ = test;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the current StatusTest used by the solver.
  template <class ScalarType, class MV, class OP>
  RCP<StatusTest<ScalarType,MV,OP> > TraceMinBase<ScalarType,MV,OP>::getStatusTest() const {
    return tester_;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Print the current status of the solver
  template <class ScalarType, class MV, class OP>
  void
  TraceMinBase<ScalarType,MV,OP>::currentStatus(std::ostream &os)
  {
    using std::endl;

    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);
    os <<endl;
    os <<"================================================================================" << endl;
    os << endl;
    os <<"                          TraceMinBase Solver Status" << endl;
    os << endl;
    os <<"The solver is "<<(initialized_ ? "initialized." : "not initialized.") << endl;
    os <<"The number of iterations performed is " <<iter_<<endl;
    os <<"The block size is         " << blockSize_<<endl;
    os <<"The number of blocks is   " << numBlocks_<<endl;
    os <<"The current basis size is " << curDim_<<endl;
    os <<"The number of auxiliary vectors is "<< numAuxVecs_ << endl;
    os <<"The number of operations Op*x   is "<<count_ApplyOp_<<endl;
    os <<"The number of operations M*x    is "<<count_ApplyM_<<endl;

    os.setf(std::ios_base::right, std::ios_base::adjustfield);

    if (initialized_) {
      os << endl;
      os <<"CURRENT EIGENVALUE ESTIMATES             "<<endl;
      os << std::setw(20) << "Eigenvalue"
         << std::setw(20) << "Residual(M)"
         << std::setw(20) << "Residual(2)"
         << endl;
      os <<"--------------------------------------------------------------------------------"<<endl;
      for (int i=0; i<blockSize_; ++i) {
        os << std::setw(20) << theta_[i];
        if (Rnorms_current_) os << std::setw(20) << Rnorms_[i];
        else os << std::setw(20) << "not current";
        if (R2norms_current_) os << std::setw(20) << R2norms_[i];
        else os << std::setw(20) << "not current";
        os << endl;
      }
    }
    os <<"================================================================================" << endl;
    os << endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  template <class ScalarType, class MV, class OP>
  ScalarType TraceMinBase<ScalarType,MV,OP>::getTrace() const
  {
    ScalarType currentTrace = ZERO;

    for(int i=0; i < blockSize_; i++)
      currentTrace += theta_[i];

    return currentTrace;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  template <class ScalarType, class MV, class OP>
  bool TraceMinBase<ScalarType,MV,OP>::traceLeveled()
  {
    ScalarType ratioOfChange = traceThresh_;

    if(previouslyLeveled_)
    {
      om_->stream(Debug) << "The trace already leveled, so we're not going to check it again\n";
      return true;
    }

    ScalarType currentTrace = getTrace();

    om_->stream(Debug) << "The current trace is " << currentTrace << std::endl;

    // Compute the ratio of the change
    // We seek the point where the trace has leveled off
    // It should be reasonably safe to shift at this point
    if(previousTrace_ != ZERO)
    {
      om_->stream(Debug) << "The previous trace was " << previousTrace_ << std::endl;

      ratioOfChange = std::abs(previousTrace_-currentTrace)/std::abs(previousTrace_);
      om_->stream(Debug) << "The ratio of change is " << ratioOfChange << std::endl;
    }

    previousTrace_ = currentTrace;

    if(ratioOfChange < traceThresh_)
    {
      previouslyLeveled_ = true;
      return true;
    }
    return false;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the residual of each CLUSTER of eigenvalues
  // This is important for selecting the Ritz shifts
  template <class ScalarType, class MV, class OP>
  std::vector<ScalarType> TraceMinBase<ScalarType,MV,OP>::getClusterResids()
  {
    int nvecs;

    if(computeAllRes_)
      nvecs = curDim_;
    else
      nvecs = blockSize_;

    getRes2Norms();

    std::vector<ScalarType> clusterResids(nvecs);
    std::vector<int> clusterIndices;
    if(considerClusters_)
    {
      for(int i=0; i < nvecs; i++)
      {
        // test for cluster
        if(clusterIndices.empty() || (theta_[i-1] + R2norms_[i-1] >= theta_[i] - R2norms_[i]))
        {
          // Add to cluster
          if(!clusterIndices.empty())  om_->stream(Debug) << theta_[i-1] << " is in a cluster with " << theta_[i] << " because " << theta_[i-1] + R2norms_[i-1] << " >= " << theta_[i] - R2norms_[i] << std::endl;
          clusterIndices.push_back(i);
        }
        // Cluster completed
        else
        {
          om_->stream(Debug) << theta_[i-1] << " is NOT in a cluster with " << theta_[i] << " because " << theta_[i-1] + R2norms_[i-1] << " < " << theta_[i] - R2norms_[i] << std::endl;
          ScalarType totalRes = ZERO;
          for(size_t j=0; j < clusterIndices.size(); j++)
            totalRes += (R2norms_[clusterIndices[j]]*R2norms_[clusterIndices[j]]);

          // If the smallest magnitude value of this sign is in a cluster with the
          // largest magnitude cluster of this sign, it is not safe for the smallest
          // eigenvalue to use a shift
          if(theta_[clusterIndices[0]] < 0 && theta_[i] < 0)
            negSafeToShift_ = true;
          else if(theta_[clusterIndices[0]] > 0 && theta_[i] > 0)
            posSafeToShift_ = true;

          for(size_t j=0; j < clusterIndices.size(); j++)
            clusterResids[clusterIndices[j]] = sqrt(totalRes);

          clusterIndices.clear();
          clusterIndices.push_back(i);
        }
      }

      // Handle last cluster
      ScalarType totalRes = ZERO;
      for(size_t j=0; j < clusterIndices.size(); j++)
        totalRes += R2norms_[clusterIndices[j]];
      for(size_t j=0; j < clusterIndices.size(); j++)
        clusterResids[clusterIndices[j]] = totalRes;
    }
    else
    {
      for(int j=0; j < nvecs; j++)
        clusterResids[j] = R2norms_[j];
    }

    return clusterResids;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the Ritz shifts based on the Ritz values and residuals
  // A note on shifting: if the matrix is indefinite, you NEED to use a large block size
  // TODO: resids[i] on its own is unsafe for the generalized EVP
  //       See "A Parallel Implementation of the Trace Minimization Eigensolver"
  //       by Eloy Romero and Jose E. Roman
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::computeRitzShifts(const std::vector<ScalarType>& clusterResids)
  {
    std::vector<ScalarType> thetaMag(theta_);
    bool traceHasLeveled = traceLeveled();

    // Compute the magnitude of the eigenvalues
    for(int i=0; i<blockSize_; i++)
      thetaMag[i] = std::abs(thetaMag[i]);

    // Test whether it is safe to shift
    // TODO: Add residual shift option
    // Note: If you choose single shift with an indefinite matrix, you're gonna have a bad time...
    // Note: This is not safe for indefinite matrices, and I don't even know that it CAN be made safe.
    //       There just isn't any theoretical reason it should work.
    // TODO: It feels like this should be different for TraceMin-Base; we should be able to use all eigenvalues in the current subspace to determine whether we have a cluster
    if(whenToShift_ == ALWAYS_SHIFT || (whenToShift_ == SHIFT_WHEN_TRACE_LEVELS && traceHasLeveled))
    {
      // Set the shift to the largest safe shift
      if(howToShift_ == LARGEST_CONVERGED_SHIFT)
      {
        for(int i=0; i<blockSize_; i++)
          ritzShifts_[i] = largestSafeShift_;
      }
      // Set the shifts to the Ritz values
      else if(howToShift_ == RITZ_VALUES_SHIFT)
      {
        ritzShifts_[0] = theta_[0];

        // If we're using mulitple shifts, set them to EACH Ritz value.
        // Otherwise, only use the smallest
        if(useMultipleShifts_)
        {
          for(int i=1; i<blockSize_; i++)
            ritzShifts_[i] = theta_[i];
        }
        else
        {
          for(int i=1; i<blockSize_; i++)
            ritzShifts_[i] = ritzShifts_[0];
        }
      }
      else if(howToShift_ == EXPERIMENTAL_SHIFT)
      {
        ritzShifts_[0] = std::max(largestSafeShift_,theta_[0]-clusterResids[0]);
        for(int i=1; i<blockSize_; i++)
        {
          ritzShifts_[i] = std::max(ritzShifts_[i-1],theta_[i]-clusterResids[i]);
        }
      }
      // Use Dr. Sameh's original shifting strategy
      else if(howToShift_ == ADJUSTED_RITZ_SHIFT)
      {
        om_->stream(Debug) << "\nSeeking a shift for theta[0]=" << thetaMag[0] << std::endl;

        // This is my adjustment.  If all eigenvalues are in a single cluster, it's probably a bad idea to shift the smallest one.
        // If all of your eigenvalues are in one cluster, it's either way to early to shift or your subspace is too small
        if((theta_[0] > 0 && posSafeToShift_) || (theta_[0] < 0 && negSafeToShift_) || considerClusters_ == false)
        {
          // Initialize with a conservative shift, either the biggest safe shift or the eigenvalue adjusted by its cluster's residual
          ritzShifts_[0] = std::max(largestSafeShift_,thetaMag[0]-clusterResids[0]);

          om_->stream(Debug) << "Initializing with a conservative shift, either the most positive converged eigenvalue ("
                             << largestSafeShift_ << ") or the eigenvalue adjusted by the residual (" << thetaMag[0] << "-"
                             << clusterResids[0] << ").\n";

          // If this eigenvalue is NOT in a cluster, do an aggressive shift
          if(R2norms_[0] == clusterResids[0])
          {
            ritzShifts_[0] = thetaMag[0];
            om_->stream(Debug) << "Since this eigenvalue is NOT in a cluster, we can use the eigenvalue itself as a shift: ritzShifts[0]=" << ritzShifts_[0] << std::endl;
          }
          else
            om_->stream(Debug) << "This eigenvalue is in a cluster, so it would not be safe to use the eigenvalue itself as a shift\n";
        }
        else
        {
          if(largestSafeShift_ > std::abs(ritzShifts_[0]))
          {
            om_->stream(Debug) << "Initializing with a conservative shift...the most positive converged eigenvalue: " << largestSafeShift_ << std::endl;
            ritzShifts_[0] = largestSafeShift_;
          }
          else
            om_->stream(Debug) << "Using the previous value of ritzShifts[0]=" << ritzShifts_[0];

        }

        om_->stream(Debug) << "ritzShifts[0]=" << ritzShifts_[0] << std::endl;

        if(useMultipleShifts_)
        {
          /////////////////////////////////////////////////////////////////////////////////////////
          // Compute shifts for other eigenvalues
          for(int i=1; i < blockSize_; i++)
          {
            om_->stream(Debug) << "\nSeeking a shift for theta[" << i << "]=" << thetaMag[i] << std::endl;

            // If the previous shift was aggressive and we are not in a cluster, do an aggressive shift
            if(ritzShifts_[i-1] == thetaMag[i-1] && i < blockSize_-1 && thetaMag[i] < thetaMag[i+1] - clusterResids[i+1])
            {
              ritzShifts_[i] = thetaMag[i];
              om_->stream(Debug) << "Using an aggressive shift: ritzShifts_[" << i << "]=" << ritzShifts_[i] << std::endl;
            }
            else
            {
              if(ritzShifts_[0] > std::abs(ritzShifts_[i]))
              {
                om_->stream(Debug) << "It was unsafe to use the aggressive shift.  Choose the shift used by theta[0]="
                                   << thetaMag[0] << ": ritzShifts[0]=" << ritzShifts_[0] << std::endl;

                // Choose a conservative shift, that of the smallest positive eigenvalue
                ritzShifts_[i] = ritzShifts_[0];
              }
              else
                om_->stream(Debug) << "It was unsafe to use the aggressive shift.  We will use the shift from the previous iteration: " << ritzShifts_[i] << std::endl;

              om_->stream(Debug) << "Check whether any less conservative shifts would work (such as the biggest eigenvalue outside of the cluster, namely theta[ell] < "
                                 << thetaMag[i] << "-" << clusterResids[i] << " (" << thetaMag[i] - clusterResids[i] << ")\n";

              // If possible, choose a less conservative shift, that of the biggest eigenvalue outside of the cluster
              for(int ell=0; ell < i; ell++)
              {
                if(thetaMag[ell] < thetaMag[i] - clusterResids[i])
                {
                  ritzShifts_[i] = thetaMag[ell];
                  om_->stream(Debug) << "ritzShifts_[" << i << "]=" << ritzShifts_[i] << " is valid\n";
                }
                else
                  break;
              }
            } // end else

            om_->stream(Debug) << "ritzShifts[" << i << "]=" << ritzShifts_[i] << std::endl;
          } // end for
        } // end if(useMultipleShifts_)
        else
        {
          for(int i=1; i<blockSize_; i++)
            ritzShifts_[i] = ritzShifts_[0];
        }
      } // end if(howToShift_ == "Adjusted Ritz Values")
    } // end if(whenToShift_ == "Always" || (whenToShift_ == "After Trace Levels" && traceHasLeveled))

    // Set the correct sign
    for(int i=0; i<blockSize_; i++)
    {
      if(theta_[i] < 0)
        ritzShifts_[i] = -abs(ritzShifts_[i]);
      else
        ritzShifts_[i] = abs(ritzShifts_[i]);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  template <class ScalarType, class MV, class OP>
  std::vector<ScalarType> TraceMinBase<ScalarType,MV,OP>::computeTol()
  {
    ScalarType temp1;
    std::vector<ScalarType> tolerances(blockSize_);

    for(int i=0; i < blockSize_-1; i++)
    {
      if(std::abs(theta_[blockSize_-1]) != std::abs(ritzShifts_[i]))
        temp1 = std::abs(theta_[i]-ritzShifts_[i])/std::abs(std::abs(theta_[blockSize_-1])-std::abs(ritzShifts_[i]));
      else
        temp1 = ZERO;

      // TODO: The min and max tolerances should not be hard coded
      //       Neither should the maximum number of iterations
      tolerances[i] = std::min(temp1*temp1,0.5);
    }

    if(blockSize_ > 1)
      tolerances[blockSize_-1] = tolerances[blockSize_-2];

    return tolerances;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::solveSaddlePointProblem(RCP<MV> Delta)
  {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor lcltimer( *timerSaddle_ );
#endif

    // This case can arise when looking for the largest eigenpairs
    if(Op_ == Teuchos::null)
    {
      // dense solver
      Teuchos::SerialDenseSolver<int,ScalarType> My_Solver;

      // Schur complement
      RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > lclS = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(blockSize_,blockSize_) );

      if(computeAllRes_)
      {
        // Get the valid indices of X
        std::vector<int> curind(blockSize_);
        for(int i=0; i<blockSize_; i++)
          curind[i] = i;

        // Get a view of MX
        RCP<const MV> lclMX = MVT::CloneView(*MX_, curind);

        // form S = X' M^2 X
        MVT::MvTransMv(ONE,*lclMX,*lclMX,*lclS);

        // compute the inverse of S
        My_Solver.setMatrix(lclS);
        My_Solver.invert();

        // Delta = X - MX/S
        RCP<const MV> lclX = MVT::CloneView(*X_, curind);
        MVT::Assign(*lclX,*Delta);
        MVT::MvTimesMatAddMv( -ONE, *lclMX, *lclS, ONE, *Delta);
      }
      else
      {
        // form S = X' M^2 X
        MVT::MvTransMv(ONE,*MX_,*MX_,*lclS);

        // compute the inverse of S
        My_Solver.setMatrix(lclS);
        My_Solver.invert();

        // Delta = X - MX/S
        MVT::Assign(*X_,*Delta);
        MVT::MvTimesMatAddMv( -ONE, *MX_, *lclS, ONE, *Delta);
      }
    }
    else
    {
      std::vector<int> order(curDim_);
      std::vector<ScalarType> tempvec(blockSize_);
//      RCP<BasicSort<MagnitudeType> > sorter = rcp( new BasicSort<MagnitudeType>("SR") );

      // Stores the residual of each CLUSTER of eigenvalues
      std::vector<ScalarType> clusterResids;

/*      // Sort the eigenvalues in ascending order for the Ritz shift selection
      sorter->sort(theta_, Teuchos::rcpFromRef(order), curDim_);   // don't catch exception

      // Apply the same ordering to the residual norms
      getRes2Norms();
      for (int i=0; i<blockSize_; i++)
        tempvec[i] = R2norms_[order[i]];
      R2norms_ = tempvec;*/

      // Compute the residual of each CLUSTER of eigenvalues
      // This is important for selecting the Ritz shifts
      clusterResids = getClusterResids();

/*      // Sort the eigenvalues based on what the user wanted
      sm_->sort(theta_, Teuchos::rcpFromRef(order), blockSize_);

      // Apply the same ordering to the residual norms and cluster residuals
      for (int i=0; i<blockSize_; i++)
        tempvec[i] = R2norms_[order[i]];
      R2norms_ = tempvec;

      for (int i=0; i<blockSize_; i++)
        tempvec[i] = clusterResids[order[i]];
      clusterResids = tempvec;*/

      // Compute the Ritz shifts
      computeRitzShifts(clusterResids);

      // Compute the tolerances for the inner solves
      std::vector<ScalarType> tolerances = computeTol();

      for(int i=0; i<blockSize_; i++)
      {
        om_->stream(IterationDetails) << "Choosing Ritz shifts...theta[" << i << "]="
            << theta_[i] << ", resids[" << i << "]=" << R2norms_[i] << ", clusterResids[" << i << "]=" << clusterResids[i]
            << ", ritzShifts[" << i << "]=" << ritzShifts_[i] << ", and tol[" << i << "]=" << tolerances[i] << std::endl;
      }

      // Set the Ritz shifts for the solver
      ritzOp_->setRitzShifts(ritzShifts_);

      // Set the inner stopping tolerance
      // This uses the Ritz values to determine when to stop
      ritzOp_->setInnerTol(tolerances);

      // Solve the saddle point problem
      if(saddleSolType_ == PROJECTED_KRYLOV_SOLVER)
      {
        if(Prec_ != Teuchos::null)
          solveSaddleProjPrec(Delta);
        else
          solveSaddleProj(Delta);
      }
      else if(saddleSolType_ == SCHUR_COMPLEMENT_SOLVER)
      {
        if(Z_ == Teuchos::null || MVT::GetNumberVecs(*Z_) != blockSize_)
        {
          // We do NOT want Z to be 0, because that could result in stagnation
          // I know it's tempting to take out the MvRandom, but seriously, don't do it.
          Z_ = MVT::Clone(*X_,blockSize_);
          MVT::MvRandom(*Z_);
        }
        solveSaddleSchur(Delta);
      }
      else if(saddleSolType_ == BD_PREC_MINRES)
      {
        solveSaddleBDPrec(Delta);
//        Delta->describe(*(Teuchos::VerboseObjectBase::getDefaultOStream()),Teuchos::VERB_EXTREME);
      }
      else if(saddleSolType_ == HSS_PREC_GMRES)
      {
        solveSaddleHSSPrec(Delta);
      }
      else
        std::cout << "Invalid saddle solver type\n";
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Solve the saddle point problem using projected minres
  // TODO: We should be able to choose KX or -R as RHS.
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::solveSaddleProj (RCP<MV> Delta) const
  {
    RCP<TraceMinProjRitzOp<ScalarType,MV,OP> > projOp;

    if(computeAllRes_)
    {
      // Get the valid indices of X
      std::vector<int> curind(blockSize_);
      for(int i=0; i<blockSize_; i++)
        curind[i] = i;

      RCP<const MV> locMX = MVT::CloneView(*MX_, curind);

      // We should really project out the converged eigenvectors too
      if(projectAllVecs_)
      {
        if(projectLockedVecs_ && numAuxVecs_ > 0)
          projOp = rcp( new TraceMinProjRitzOp<ScalarType,MV,OP>(ritzOp_,locMX,orthman_,auxVecs_) );
        else
          projOp = rcp( new TraceMinProjRitzOp<ScalarType,MV,OP>(ritzOp_,locMX,orthman_) );
      }
      else
        projOp = rcp( new TraceMinProjRitzOp<ScalarType,MV,OP>(ritzOp_,locMX) );

      // Remember, Delta0 must equal 0
      // This ensures B-orthogonality between Delta and X
      MVT::MvInit(*Delta);

      if(useRHSR_)
      {
        RCP<const MV> locR = MVT::CloneView(*R_, curind);
        projOp->ApplyInverse(*locR, *Delta);
      }
      else
      {
        RCP<const MV> locKX = MVT::CloneView(*KX_, curind);
        projOp->ApplyInverse(*locKX, *Delta);
      }
    }
    else
    {
      // We should really project out the converged eigenvectors too
      if(projectAllVecs_)
      {
        if(projectLockedVecs_ && numAuxVecs_ > 0)
          projOp = rcp( new TraceMinProjRitzOp<ScalarType,MV,OP>(ritzOp_,MX_,orthman_,auxVecs_) );
        else
          projOp = rcp( new TraceMinProjRitzOp<ScalarType,MV,OP>(ritzOp_,MX_,orthman_) );
      }
      else
        projOp = rcp( new TraceMinProjRitzOp<ScalarType,MV,OP>(ritzOp_,MX_) );

      // Remember, Delta0 must equal 0
      // This ensures B-orthogonality between Delta and X
      MVT::MvInit(*Delta);

      if(useRHSR_) {
        projOp->ApplyInverse(*R_, *Delta);
      }
      else {
        projOp->ApplyInverse(*KX_, *Delta);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TODO: Fix preconditioning
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::solveSaddleProjPrec (RCP<MV> Delta) const
  {
    // If we don't have Belos installed, we can't use TraceMinProjRitzOpWithPrec
    // Of course, this problem will be detected in the constructor and an exception will be thrown
    // This is only here to make sure the code compiles properly
#ifdef HAVE_ANASAZI_BELOS
    RCP<TraceMinProjRitzOpWithPrec<ScalarType,MV,OP> > projOp;

    if(computeAllRes_)
    {
          int dimension;
          if(projectAllVecs_)
            dimension = curDim_;
          else
            dimension = blockSize_;

      // Get the valid indices of X
      std::vector<int> curind(dimension);
      for(int i=0; i<dimension; i++)
        curind[i] = i;

      RCP<const MV> locMX = MVT::CloneView(*MX_, curind);

      // We should really project out the converged eigenvectors too
      if(projectAllVecs_)
      {
        if(projectLockedVecs_ && numAuxVecs_ > 0)
          projOp = rcp( new TraceMinProjRitzOpWithPrec<ScalarType,MV,OP>(ritzOp_,locMX,orthman_,auxVecs_) );
        else
          projOp = rcp( new TraceMinProjRitzOpWithPrec<ScalarType,MV,OP>(ritzOp_,locMX,orthman_) );
      }
      else
        projOp = rcp( new TraceMinProjRitzOpWithPrec<ScalarType,MV,OP>(ritzOp_,locMX) );

      // Remember, Delta0 must equal 0
      // This ensures B-orthogonality between Delta and X
      MVT::MvInit(*Delta);

          std::vector<int> dimind(blockSize_);
      for(int i=0; i<blockSize_; i++)
        dimind[i] = i;

      if(useRHSR_)
      {
        RCP<const MV> locR = MVT::CloneView(*R_, dimind);
        projOp->ApplyInverse(*locR, *Delta);
        MVT::MvScale(*Delta, -ONE);
      }
      else
      {
        RCP<const MV> locKX = MVT::CloneView(*KX_, dimind);
        projOp->ApplyInverse(*locKX, *Delta);
      }
    }
    else
    {
      // We should really project out the converged eigenvectors too
      if(projectAllVecs_)
      {
        if(projectLockedVecs_ && numAuxVecs_ > 0)
          projOp = rcp( new TraceMinProjRitzOpWithPrec<ScalarType,MV,OP>(ritzOp_,MX_,orthman_,auxVecs_) );
        else
          projOp = rcp( new TraceMinProjRitzOpWithPrec<ScalarType,MV,OP>(ritzOp_,MX_,orthman_) );
      }
      else
        projOp = rcp( new TraceMinProjRitzOpWithPrec<ScalarType,MV,OP>(ritzOp_,MX_) );

      // Remember, Delta0 must equal 0
      // This ensures B-orthogonality between Delta and X
      MVT::MvInit(*Delta);

      if(useRHSR_)
      {
        projOp->ApplyInverse(*R_, *Delta);
        MVT::MvScale(*Delta,-ONE);
      }
      else
        projOp->ApplyInverse(*KX_, *Delta);
    }
#endif
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TODO: We can hold the Schur complement constant in later iterations
  // TODO: Make sure we're using the preconditioner correctly
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::solveSaddleSchur (RCP<MV> Delta) const
  {
    // dense solver
    Teuchos::SerialDenseSolver<int,ScalarType> My_Solver;

    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > lclL;
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > lclS = rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(blockSize_,blockSize_) );

    if(computeAllRes_)
    {
      // Get the valid indices of X
      std::vector<int> curind(blockSize_);
      for(int i=0; i<blockSize_; i++)
        curind[i] = i;

      // Z = K \ MX
      // Why would I have wanted to set Z <- X?  I want to leave Z's previous value alone...
      RCP<const MV> lclMX = MVT::CloneView(*MX_, curind);

#ifdef USE_APPLY_INVERSE
      Op_->ApplyInverse(*lclMX,*Z_);
#else
      ritzOp_->ApplyInverse(*lclMX,*Z_);
#endif

      // form S = X' M Z
      MVT::MvTransMv(ONE,*Z_,*lclMX,*lclS);

      // solve S L = I
      My_Solver.setMatrix(lclS);
      My_Solver.invert();
      lclL = lclS;

      // Delta = X - Z L
      RCP<const MV> lclX = MVT::CloneView(*X_, curind);
      MVT::Assign(*lclX,*Delta);
      MVT::MvTimesMatAddMv( -ONE, *Z_, *lclL, ONE, *Delta);
    }
    else
    {
      // Z = K \ MX
#ifdef USE_APPLY_INVERSE
      Op_->ApplyInverse(*MX_,*Z_);
#else
      ritzOp_->ApplyInverse(*MX_,*Z_);
#endif

      // form S = X' M Z
      MVT::MvTransMv(ONE,*Z_,*MX_,*lclS);

      // solve S L = I
      My_Solver.setMatrix(lclS);
      My_Solver.invert();
      lclL = lclS;

      // Delta = X - Z L
      MVT::Assign(*X_,*Delta);
      MVT::MvTimesMatAddMv( -ONE, *Z_, *lclL, ONE, *Delta);
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TODO: We can hold the Schur complement constant in later iterations
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::solveSaddleBDPrec (RCP<MV> Delta) const
  {
    RCP<MV> locKX, locMX;
    if(computeAllRes_)
    {
      std::vector<int> curind(blockSize_);
      for(int i=0; i<blockSize_; i++)
        curind[i] = i;
      locKX = MVT::CloneViewNonConst(*KX_, curind);
      locMX = MVT::CloneViewNonConst(*MX_, curind);
    }
    else
    {
      locKX = KX_;
      locMX = MX_;
    }

    // Create the operator [A BX; X'B 0]
    RCP<saddle_op_type> sadOp = rcp(new saddle_op_type(ritzOp_,locMX));

    // Create the RHS [AX; 0]
    RCP<saddle_container_type> sadRHS = rcp(new saddle_container_type(locKX));

//    locKX->describe(*(Teuchos::VerboseObjectBase::getDefaultOStream()),Teuchos::VERB_EXTREME);

//    locMX->describe(*(Teuchos::VerboseObjectBase::getDefaultOStream()),Teuchos::VERB_EXTREME);

    // Create the solution vector [Delta; L]
    MVT::MvInit(*Delta);
    RCP<saddle_container_type> sadSol = rcp(new saddle_container_type(Delta));

    // Create a minres solver
    RCP<PseudoBlockMinres<ScalarType,saddle_container_type,saddle_op_type > > sadSolver;
    if(Prec_ != Teuchos::null)
    {
      RCP<saddle_op_type> sadPrec = rcp(new saddle_op_type(ritzOp_->getPrec(),locMX,BD_PREC));
      sadSolver = rcp(new PseudoBlockMinres<ScalarType,saddle_container_type,saddle_op_type>(sadOp, sadPrec));
    }
    else {
      sadSolver = rcp(new PseudoBlockMinres<ScalarType,saddle_container_type,saddle_op_type>(sadOp));
    }

    // Set the tolerance for the minres solver
    std::vector<ScalarType> tol;
    ritzOp_->getInnerTol(tol);
    sadSolver->setTol(tol);

    // Set the maximum number of iterations
    sadSolver->setMaxIter(ritzOp_->getMaxIts());

    // Set the solution vector
    sadSolver->setSol(sadSol);

    // Set the RHS
    sadSolver->setRHS(sadRHS);

    // Solve the saddle point problem
    sadSolver->solve();
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TODO: We can hold the Schur complement constant in later iterations
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::solveSaddleHSSPrec (RCP<MV> Delta) const
  {
#ifdef HAVE_ANASAZI_BELOS
    typedef Belos::LinearProblem<ScalarType,saddle_container_type,saddle_op_type>           LP;
    typedef Belos::PseudoBlockGmresSolMgr<ScalarType,saddle_container_type,saddle_op_type>  GmresSolMgr;

    RCP<MV> locKX, locMX;
    if(computeAllRes_)
    {
      std::vector<int> curind(blockSize_);
      for(int i=0; i<blockSize_; i++)
        curind[i] = i;
      locKX = MVT::CloneViewNonConst(*KX_, curind);
      locMX = MVT::CloneViewNonConst(*MX_, curind);
    }
    else
    {
      locKX = KX_;
      locMX = MX_;
    }

    // Create the operator [A BX; X'B 0]
    RCP<saddle_op_type> sadOp = rcp(new saddle_op_type(ritzOp_,locMX,NONSYM));

    // Create the RHS [AX; 0]
    RCP<saddle_container_type> sadRHS = rcp(new saddle_container_type(locKX));

    // Create the solution vector [Delta; L]
    MVT::MvInit(*Delta);
    RCP<saddle_container_type> sadSol = rcp(new saddle_container_type(Delta));

    // Create a parameter list for the gmres solver
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

    // Set the tolerance for the gmres solver
    std::vector<ScalarType> tol;
    ritzOp_->getInnerTol(tol);
    pl->set("Convergence Tolerance",tol[0]);

    // Set the maximum number of iterations
    // TODO: Come back to this
    pl->set("Maximum Iterations", ritzOp_->getMaxIts());
    pl->set("Num Blocks", ritzOp_->getMaxIts());

    // Set the block size
    // TODO: Come back to this
        // TODO: This breaks the code right now, presumably because of a MVT cloneview issue.
    pl->set("Block Size", blockSize_);

    // Set the verbosity of gmres
//    pl->set("Verbosity", Belos::IterationDetails + Belos::StatusTestDetails + Belos::Debug);
//    pl->set("Output Frequency", 1);

    // Create the linear problem
    RCP<LP> problem = rcp(new LP(sadOp,sadSol,sadRHS));

    // Set the preconditioner
    if(Prec_ != Teuchos::null)
    {
      RCP<saddle_op_type> sadPrec = rcp(new saddle_op_type(ritzOp_->getPrec(),locMX,HSS_PREC,alpha_));
      problem->setLeftPrec(sadPrec);
    }

    // Set the problem
    problem->setProblem();

    // Create a minres solver
    RCP<GmresSolMgr> sadSolver = rcp(new GmresSolMgr(problem,pl)) ;

    // Solve the saddle point problem
    sadSolver->solve();
#else
    std::cout << "No Belos.  This is bad\n";
#endif
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute KK := V'KV
  // We only compute the NEW elements
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::computeKK()
  {
    // Get the valid indices of V
    std::vector<int> curind(curDim_);
    for(int i=0; i<curDim_; i++)
      curind[i] = i;

    // Get a pointer to the valid parts of V
    RCP<const MV> lclV = MVT::CloneView(*V_,curind);

    // Get the valid indices of KV
    curind.resize(blockSize_);
    for(int i=0; i<blockSize_; i++)
      curind[i] = curDim_-blockSize_+i;
    RCP<const MV> lclKV = MVT::CloneView(*KV_,curind);

    // Get a pointer to the valid part of KK
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > lclKK =
        rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*KK_,curDim_,blockSize_,0,curDim_-blockSize_) );

    // KK := V'KV
    MVT::MvTransMv(ONE,*lclV,*lclKV,*lclKK);

    // We only constructed the upper triangular part of the matrix, but that's okay because KK is symmetric!
    for(int r=0; r<curDim_; r++)
    {
      for(int c=0; c<r; c++)
      {
        (*KK_)(r,c) = (*KK_)(c,r);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the eigenpairs of KK, i.e. the Ritz pairs
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::computeRitzPairs()
  {
    // Get a pointer to the valid part of KK
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > lclKK =
        rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*KK_,curDim_,curDim_) );

    // Get a pointer to the valid part of ritzVecs
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > lclRV =
        rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*ritzVecs_,curDim_,curDim_) );

    // Compute Ritz pairs from KK
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerDS_ );
#endif
      int rank = curDim_;
      Utils::directSolver(curDim_, *lclKK, Teuchos::null, *lclRV, theta_, rank, 10);
      // we want all ritz values back
      // TODO: This probably should not be an ortho failure
      TEUCHOS_TEST_FOR_EXCEPTION(rank != curDim_,TraceMinBaseOrthoFailure,
             "Anasazi::TraceMinBase::computeRitzPairs(): Failed to compute all eigenpairs of KK.");
    }

    // Sort ritz pairs
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *timerSortEval_ );
#endif

      std::vector<int> order(curDim_);
      //
      // sort the first curDim_ values in theta_
      if(useHarmonic_)
      {
        Anasazi::BasicSort<ScalarType> sm;
        sm.sort(theta_, Teuchos::rcpFromRef(order), curDim_);
      }
      else
      {
        sm_->sort(theta_, Teuchos::rcpFromRef(order), curDim_);   // don't catch exception
      }
      //
      // apply the same ordering to the primitive ritz vectors
      Utils::permuteVectors(order,*lclRV);
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute X := V evecs
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::computeX()
  {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor lcltimer( *timerLocal_ );
#endif

    // Get the valid indices of V
    std::vector<int> curind(curDim_);
    for(int i=0; i<curDim_; i++)
      curind[i] = i;

    // Get a pointer to the valid parts of V
    RCP<const MV> lclV = MVT::CloneView(*V_,curind);

    if(computeAllRes_)
    {
      // Capture the relevant eigenvectors
      RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > relevantEvecs =
          rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*ritzVecs_,curDim_,curDim_) );

      // X <- lclV*S
      RCP<MV> lclX = MVT::CloneViewNonConst(*X_,curind);
      MVT::MvTimesMatAddMv( ONE, *lclV, *relevantEvecs, ZERO, *lclX );
    }
    else
    {
      // Capture the relevant eigenvectors
      RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > relevantEvecs =
          rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*ritzVecs_,curDim_,blockSize_) );

      // X <- lclV*S
      MVT::MvTimesMatAddMv( ONE, *lclV, *relevantEvecs, ZERO, *X_ );
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute KX := KV evecs and (if necessary) MX := MV evecs
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::updateKXMX()
  {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor lcltimer( *timerLocal_ );
#endif

    // Get the valid indices of V
    std::vector<int> curind(curDim_);
    for(int i=0; i<curDim_; i++)
      curind[i] = i;

    // Get pointers to the valid parts of V, KV, and MV (if necessary)
    RCP<const MV> lclV = MVT::CloneView(*V_,curind);
    RCP<const MV> lclKV = MVT::CloneView(*KV_,curind);

    if(computeAllRes_)
    {
      // Capture the relevant eigenvectors
      RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > relevantEvecs =
          rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*ritzVecs_,curDim_,curDim_) );

      // Update KX and MX
      RCP<MV> lclKX = MVT::CloneViewNonConst(*KX_,curind);
      MVT::MvTimesMatAddMv( ONE, *lclKV, *relevantEvecs, ZERO, *lclKX );
      if(hasM_)
      {
        RCP<const MV> lclMV = MVT::CloneView(*MV_,curind);
        RCP<MV> lclMX = MVT::CloneViewNonConst(*MX_,curind);
        MVT::MvTimesMatAddMv( ONE, *lclMV, *relevantEvecs, ZERO, *lclMX );
      }
    }
    else
    {
      // Capture the relevant eigenvectors
      RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > relevantEvecs =
          rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*ritzVecs_,curDim_,blockSize_) );

      // Update KX and MX
      MVT::MvTimesMatAddMv( ONE, *lclKV, *relevantEvecs, ZERO, *KX_ );
      if(hasM_)
      {
        RCP<const MV> lclMV = MVT::CloneView(*MV_,curind);
        MVT::MvTimesMatAddMv( ONE, *lclMV, *relevantEvecs, ZERO, *MX_ );
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Update the residual R := KX - MX*T
  template <class ScalarType, class MV, class OP>
  void TraceMinBase<ScalarType,MV,OP>::updateResidual () {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor lcltimer( *timerCompRes_ );
#endif

    if(computeAllRes_)
    {
      // Get the valid indices of X
      std::vector<int> curind(curDim_);
      for(int i=0; i<curDim_; i++)
        curind[i] = i;

      // Holds MX*T
      RCP<MV> MXT = MVT::CloneCopy(*MX_,curind);

      // Holds the relevant part of theta
      std::vector<ScalarType> locTheta(curDim_);
      for(int i=0; i<curDim_; i++)
        locTheta[i] = theta_[i];

      // Compute MX*T
      MVT::MvScale(*MXT,locTheta);

      // form R <- KX - MX*T
      RCP<const MV> locKX = MVT::CloneView(*KX_,curind);
      RCP<MV> locR = MVT::CloneViewNonConst(*R_,curind);
      MVT::MvAddMv(ONE,*locKX,-ONE,*MXT,*locR);
    }
    else
    {
      // Holds MX*T
      RCP<MV> MXT = MVT::CloneCopy(*MX_);

      // Holds the relevant part of theta
      std::vector<ScalarType> locTheta(blockSize_);
      for(int i=0; i<blockSize_; i++)
        locTheta[i] = theta_[i];

      // Compute MX*T
      MVT::MvScale(*MXT,locTheta);

      // form R <- KX - MX*T
      MVT::MvAddMv(ONE,*KX_,-ONE,*MXT,*R_);
    }

    // R has been updated; mark the norms as out-of-date
    Rnorms_current_ = false;
    R2norms_current_ = false;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check accuracy, orthogonality, and other debugging stuff
  //
  // bools specify which tests we want to run (instead of running more than we actually care about)
  //
  // we don't bother checking the following because they are computed explicitly:
  //    H == Prec*R
  //   KH == K*H
  //
  //
  // checkV : V orthonormal
  //          orthogonal to auxvecs
  // checkX : X orthonormal
  //          orthogonal to auxvecs
  // checkMX: check MX == M*X
  // checkKX: check KX == K*X
  // checkH : H orthonormal
  //          orthogonal to V and H and auxvecs
  // checkMH: check MH == M*H
  // checkR : check R orthogonal to X
  // checkQ : check that auxiliary vectors are actually orthonormal
  // checkKK: check that KK is symmetric in memory
  //
  // TODO:
  //  add checkTheta
  //
  template <class ScalarType, class MV, class OP>
  std::string TraceMinBase<ScalarType,MV,OP>::accuracyCheck( const CheckList &chk, const std::string &where ) const
  {
    using std::endl;

    std::stringstream os;
    os.precision(2);
    os.setf(std::ios::scientific, std::ios::floatfield);

    os << " Debugging checks: iteration " << iter_ << where << endl;

    // V and friends
    std::vector<int> lclind(curDim_);
    for (int i=0; i<curDim_; ++i) lclind[i] = i;
    RCP<const MV> lclV;
    if (initialized_) {
      lclV = MVT::CloneView(*V_,lclind);
    }
    if (chk.checkV && initialized_) {
      MagnitudeType err = orthman_->orthonormError(*lclV);
      os << " >> Error in V^H M V == I  : " << err << endl;
      for (Array_size_type i=0; i<auxVecs_.size(); ++i) {
        err = orthman_->orthogError(*lclV,*auxVecs_[i]);
        os << " >> Error in V^H M Q[" << i << "] == 0 : " << err << endl;
      }
    }

    // X and friends
    RCP<const MV> lclX;
    if(initialized_)
    {
      if(computeAllRes_)
        lclX = MVT::CloneView(*X_,lclind);
      else
        lclX = X_;
    }

    if (chk.checkX && initialized_) {
      MagnitudeType err = orthman_->orthonormError(*lclX);
      os << " >> Error in X^H M X == I  : " << err << endl;
      for (Array_size_type i=0; i<auxVecs_.size(); ++i) {
        err = orthman_->orthogError(*lclX,*auxVecs_[i]);
        os << " >> Error in X^H M Q[" << i << "] == 0 : " << err << endl;
      }
    }
    if (chk.checkMX && hasM_ && initialized_) {
      RCP<const MV> lclMX;
      if(computeAllRes_)
        lclMX = MVT::CloneView(*MX_,lclind);
      else
        lclMX = MX_;

      MagnitudeType err = Utils::errorEquality(*lclX, *lclMX, MOp_);
      os << " >> Error in MX == M*X     : " << err << endl;
    }
    if (Op_ != Teuchos::null && chk.checkKX && initialized_) {
      RCP<const MV> lclKX;
      if(computeAllRes_)
        lclKX = MVT::CloneView(*KX_,lclind);
      else
        lclKX = KX_;

      MagnitudeType err = Utils::errorEquality(*lclX, *lclKX, Op_);
      os << " >> Error in KX == K*X     : " << err << endl;
    }

    // KK
    if (chk.checkKK && initialized_) {
      Teuchos::SerialDenseMatrix<int,ScalarType> curKK(curDim_,curDim_);
      if(Op_ != Teuchos::null) {
        RCP<MV> lclKV = MVT::Clone(*V_,curDim_);
        OPT::Apply(*Op_,*lclV,*lclKV);
        MVT::MvTransMv(ONE,*lclV,*lclKV,curKK);
      }
      else {
        MVT::MvTransMv(ONE,*lclV,*lclV,curKK);
      }
      Teuchos::SerialDenseMatrix<int,ScalarType> subKK(Teuchos::View,*KK_,curDim_,curDim_);
      curKK -= subKK;
      os << " >> Error in V^H K V == KK : " << curKK.normFrobenius() << endl;

      Teuchos::SerialDenseMatrix<int,ScalarType> SDMerr(curDim_,curDim_);
      for (int j=0; j<curDim_; ++j) {
        for (int i=0; i<curDim_; ++i) {
          SDMerr(i,j) = subKK(i,j) - SCT::conjugate(subKK(j,i));
        }
      }
      os << " >> Error in KK - KK^H == 0 : " << SDMerr.normFrobenius() << endl;
    }

    // Q
    if (chk.checkQ) {
      for (Array_size_type i=0; i<auxVecs_.size(); ++i) {
        MagnitudeType err = orthman_->orthonormError(*auxVecs_[i]);
        os << " >> Error in Q[" << i << "]^H M Q[" << i << "] == I : " << err << endl;
        for (Array_size_type j=i+1; j<auxVecs_.size(); ++j) {
          err = orthman_->orthogError(*auxVecs_[i],*auxVecs_[j]);
          os << " >> Error in Q[" << i << "]^H M Q[" << j << "] == 0 : " << err << endl;
        }
      }
    }

    os << endl;

    return os.str();
  }

}} // End of namespace Anasazi

#endif

// End of file AnasaziTraceMinBase.hpp
