// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file BelosBlockGCRODRSolMgr.hpp
/// \brief A solver manager for the Block GCRO-DR (Block Recycling GMRES) linear solver.
///
#ifndef BELOS_BLOCK_GCRODR_SOLMGR_HPP
#define BELOS_BLOCK_GCRODR_SOLMGR_HPP

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosOrthoManagerFactory.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"
#include "BelosGmresIteration.hpp"
#include "BelosBlockGCRODRIter.hpp"
#include "BelosBlockGmresIter.hpp"
#include "BelosBlockFGmresIter.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"
#include "BelosOutputManager.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_as.hpp"
#ifdef BELOS_TEUCHOS_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif // BELOS_TEUCHOS_TIME_MONITOR

// MLP: Add links to examples here

namespace Belos{

  //! @name BlockGCRODRSolMgr Exceptions
  //@{

  /// \class BlockGCRODRSolMgrLinearProblemFailure
  /// \brief Thrown when the linear problem was not set up correctly.
  ///
  /// This exception is thrown if \c setProblem() is not called before \c solve() is called.
  class BlockGCRODRSolMgrLinearProblemFailure : public BelosError {
  public:
    BlockGCRODRSolMgrLinearProblemFailure(const std::string& what_arg) : BelosError(what_arg) {}
  };

  /// \class BlockGCRODRSolMgrOrthoFailure
  /// \brief Thrown when the solution manager was unable to orthogonalize the basis vectors.
  ///
  /// This exception is thrown from the \c solve() method if the orthogonalization manager is unable to generate orthonormal columns from the initial basis vectors.
  class BlockGCRODRSolMgrOrthoFailure : public BelosError {
  public:
    BlockGCRODRSolMgrOrthoFailure(const std::string& what_arg) : BelosError(what_arg) {}
  };

  /// \class BlockGCRODRSolMgrLAPACKFailure
  /// \brief Thrown when an LAPACK call fails.
  ///
  /// This exception is thrown if an LAPACK call return an error code.
  class BlockGCRODRSolMgrLAPACKFailure : public BelosError {
  public:
    BlockGCRODRSolMgrLAPACKFailure(const std::string& what_arg) : BelosError(what_arg) {}
  };

  /// \class BlockGCRODRSolMgrRecyclingFailure
  /// \brief Thrown if any problem occurs in using or creating the recycle subspace.
  ///
  /// This exception is thrown from the BlockGCRODRSolMgr::solve() method.
  class BlockGCRODRSolMgrRecyclingFailure : public BelosError {
  public:
    BlockGCRODRSolMgrRecyclingFailure(const std::string& what_arg) : BelosError(what_arg) {}
  };

  //@}

/// \class BlockGCRODRSolMgr
/// \brief A solver manager for the Block GCRO-DR (Block Recycling GMRES) linear solver.
/// \ingroup belos_solver_framework
/// \author Kirk M. Soodhalter and Michael L. Parks
///
/// GCRO-DR (also called Recycling GMRES) is a variant of GMRES that can more efficiently solve sequences linear systems.
/// It does so by "recycling" Krylov basis information from previous solves.
///
/// The original GCRO-DR algorithm can only solve one right-hand side at a time.  Block GCRO-DR extends GCRO-DR so that it can solve
/// multiple right-hand sides at a time; thus, it can solve sequences of block systems.
///

template<class ScalarType, class MV, class OP>
class BlockGCRODRSolMgr : public SolverManager<ScalarType, MV, OP> {
private:

  typedef MultiVecTraits<ScalarType,MV> MVT;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> MT;
  typedef Teuchos::ScalarTraits<MagnitudeType> SMT;
  typedef OrthoManagerFactory<ScalarType, MV, OP> ortho_factory_type;
  typedef Teuchos::SerialDenseMatrix<int,ScalarType> SDM;
  typedef Teuchos::SerialDenseVector<int,ScalarType> SDV;

public:
  //! @name Constructors/Destructor
  //@{

  /// \brief Default constructor.
  ///
  /// This constructor sets up the solver with default parameters.
  /// The linear problem must be passed in using \c setProblem() before \c solve() is called on this object.
  /// The solver values can be changed using \c setParameters().
  BlockGCRODRSolMgr();

  /*! \brief Basic constructor for GCRODRSolMgr.
   *
   * This constructor accepts the LinearProblem to be solved in
   * addition to a parameter list of options for the solver manager.
   * Some of the more important options include the following:
   * - "Num Blocks": an \c int specifying the number of blocks allocated for the Krylov basis. Default: 50.
   * - "Block Size": an \c int specifying the number of right hand sides being solved at a time.
   * - "Num Recycled Blocks": an \c int specifying the number of blocks allocated for the Krylov basis. Default: 5.
   * - "Maximum Iterations": an \c int specifying the maximum number of iterations the underlying solver is allowed to perform. Default: 5000.
   * - "Maximum Restarts": an \c int specifying the maximum number of restarts the underlying solver is allowed to perform. Default: 100.
   * - "Orthogonalization": an \c std::string specifying the desired orthogonalization. Currently supported values: "DGKS", "ICGS", "IMGS", and "TSQR" (if Belos was built with TSQR support). Default: "ICGS".
   * - "Orthogonalization Parameters": a ParameterList or RCP<(const) ParameterList> of parameters specific to the type of orthogonalization used. Defaults are set automatically.
   * - "Verbosity": a sum of MsgType specifying the verbosity. Default: Belos::Errors.
   * - "Output Style": a OutputType specifying the style of output. Default: Belos::General.
   * - "Convergence Tolerance": a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: 1e-8.
   *
   * Other supported options:

   * - "Output Frequency": an int specifying how often (in terms of number of iterations) convergence information should be output to the output stream. Default: -1 (never output convergence information).
   * - "Output Stream": a reference-counted pointer to the output stream where all solver output is sent. Default stream is std::cout (stdout, in C terms). For stderr, supply Teuchos::rcp(&std::cerr, false).
   * - "Implicit Residual Scaling": the type of scaling used in the implicit residual convergence test. Default: "Norm of Preconditioned Initial Residual".
   * - "Explicit Residual Scaling": the type of scaling used in the explicit residual convergence test. Default: "Norm of Initial Residual".
   * - "Timer Label": the string to use as a prefix for the timer labels. Default: "Belos"
   * - "Orthogonalization Constant": a \c MagnitudeType corresponding to the "depTol" parameter of DGKS orthogonalization. Ignored unless DGKS orthogonalization is used. DGKS decides the default value.
   */
  BlockGCRODRSolMgr (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                     const Teuchos::RCP<Teuchos::ParameterList> &pl);

  //! Destructor.
  virtual ~BlockGCRODRSolMgr() {};
  //@}

  /** \name Implementation of the Teuchos::Describable interface */
  //@{

  //! A description of the Block GCRODR solver manager.
  std::string description() const;

  //@}


  //! @name Accessor methods
  //@{

  //! Get current linear problem being solved for in this object.
  const LinearProblem<ScalarType,MV,OP>& getProblem() const {
    return *problem_;
  }

  //! Get a parameter list containing the valid parameters for this object.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //! Get a parameter list containing the current parameters for this object.
  Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const { return params_; }

  //! Get the residual for the most recent call to solve().
  MagnitudeType achievedTol() const {
    return achievedTol_;
  }

  //! Get the iteration count for the most recent call to \c solve().
  int getNumIters() const { return numIters_; }

  //! Whether a loss of accuracy was detected during the most recent solve.
  bool isLOADetected() const { return loaDetected_; }

  //@}

  //! @name Set methods
  //@{

  /// \brief Set the linear problem to solve on the next call to \c solve().
  void setProblem (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >& problem) {
    TEUCHOS_TEST_FOR_EXCEPTION(problem.is_null(), std::invalid_argument,
      "Belos::BlockGCRODRSolMgr::setProblem: The input LinearProblem cannot be null.");

    // Check state of problem object before proceeding
    if (! problem->isProblemSet()) {
      const bool success = problem->setProblem();
      TEUCHOS_TEST_FOR_EXCEPTION(success, std::runtime_error,
        "Belos::BlockGCRODRSolMgr::setProblem: Calling the input LinearProblem's setProblem() method failed.  This likely means that the "
        "LinearProblem has a missing (null) matrix A, solution vector X, or right-hand side vector B.  Please set these items in the LinearProblem and try again.");
    }

    problem_ = problem;
  }

  //! Set the parameters the solver should use to solve the linear problem.
  void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params );

  //@}

  //! @name Reset methods
  //@{

  /*! \brief Performs a reset of the solver manager specified by the \c ResetType.
   *
   * This informs the solver manager that the solver should prepare
   * for the next call to solve by resetting certain elements of the
   * iterative solver strategy.
   */
  void reset (const ResetType type) {
    if ((type & Belos::Problem) && ! problem_.is_null ())
      problem_->setProblem();
    else if (type & Belos::RecycleSubspace)
      keff = 0;
  }

  //@}

  //! @name Solver application methods
  //@{

  /// \brief Solve the current linear problem.
  ///
  /// This method performs possibly repeated calls to the underlying linear solver's \c iterate() routine until the problem has been
  /// solved (as decided by the solver manager) or the solver manager decides to quit.
  ///
  /// This method calls BlockGCRODRIter::iterate(), which will return either because a specially constructed status test evaluates to
  /// ::Passed or an exception is thrown.
  ///
  /// A return from BlockGCRODRIter::iterate() signifies one of the
  /// following scenarios:
  /// - the maximum number of restarts has been exceeded. In this scenario, the current solutions to the linear system will be
  ///   placed in the linear problem and return ::Unconverged.
  /// - global convergence has been met. In this case, the current solutions to the linear system will be placed in the linear
  ///   problem and the solver manager will return ::Converged.
  ///
  // \returns ::ReturnType specifying:
  /// - ::Converged: the linear problem was solved to the specification required by the solver manager.
  /// - ::Unconverged: the linear problem was not solved to the specification desired by the solver manager.
  ReturnType solve();

  //@}

private:

  // Called by all constructors; Contains init instructions common to all constructors
  void init();

  // Initialize solver state storage
  void initializeStateStorage();

  // Recycling Methods
  // Appending Function name by:
  //  "Kryl" indicates it is specialized for building a recycle space after an
  //         initial run of Block GMRES which generates an initial unaugmented block Krylov subspace
  //  "AugKryl" indicates  it is specialized for building a recycle space from the augmented Krylov subspace

  // Functions which control the building of a recycle space
  void buildRecycleSpaceKryl(int& keff, Teuchos::RCP<BlockGmresIter<ScalarType,MV,OP> > block_gmres_iter);
  void buildRecycleSpaceAugKryl(Teuchos::RCP<BlockGCRODRIter<ScalarType,MV,OP> > gcrodr_iter);

  // Recycling with Harmonic Ritz Vectors
  // Computes harmonic eigenpairs of projected matrix created during the priming solve.
  // The return value is the number of vectors needed to be stored, recycledBlocks or recycledBlocks+1.

  // HH is the projected problem from the initial cycle of Gmres, it is (at least) of dimension (numBlocks+1)*blockSize x numBlocks.
  // PP contains the harmonic eigenvectors corresponding to the recycledBlocks eigenvalues of smallest magnitude.
  int getHarmonicVecsKryl (int m, const SDM& HH, SDM& PP);

  // HH is the total block projected problem from the GCRO-DR algorithm, it is (at least) of dimension keff+(numBlocks+1)*blockSize x keff+numBlocksm.
  // VV is the Krylov vectors from the projected GMRES algorithm, which has (at least) (numBlocks+1)*blockSize vectors.
  // PP contains the harmonic eigenvectors corresponding to the recycledBlocks eigenvalues of smallest magnitude.
  int getHarmonicVecsAugKryl (int keff,
                              int m,
                              const SDM& HH,
                              const Teuchos::RCP<const MV>& VV,
                              SDM& PP);

  // Sort list of n floating-point numbers and return permutation vector
  void sort (std::vector<MagnitudeType>& dlist, int n, std::vector<int>& iperm);

  // Lapack interface
  Teuchos::LAPACK<int,ScalarType> lapack;

  //! The current linear problem to solve.
  Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;

  //Output Manager
  Teuchos::RCP<OutputManager<ScalarType> > printer_;
  Teuchos::RCP<std::ostream> outputStream_;

  //Status Test
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > sTest_;
  Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxIterTest_;
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > convTest_;
  Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > expConvTest_, impConvTest_;
  Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest_;

  //! Factory for creating MatOrthoManager subclass instances.
  ortho_factory_type orthoFactory_;

  /// \brief Orthogonalization manager.
  ///
  /// This is created by the OrthoManagerFactory instance.  The
  /// pointer may be invalidated if this solver's parameters are
  /// changed.
  Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > ortho_;

  //! This solver's current parameter list.
  Teuchos::RCP<Teuchos::ParameterList> params_;

  /// \brief Default parameter list.
  ///
  /// This is declared "mutable" so that it can be created on demand
  /// by \c getValidParameters().  The \c SolverManager interface
  /// forces \c getValidParameters() to be a const method.
  ///
  /// Using the "mutable" keyword is ugly, but the previous caching
  /// solution involved static method data in getValidParameters().
  /// That solution made creating multiple instances of this solver
  /// in different threads not safe.
  mutable Teuchos::RCP<const Teuchos::ParameterList> defaultParams_;

  //Default Solver Values
  static const bool adaptiveBlockSize_default_;
  static const std::string recycleMethod_default_;

  //Current Solver Values
  MagnitudeType convTol_, orthoKappa_, achievedTol_;
  int blockSize_, maxRestarts_, maxIters_, numIters_;
  int verbosity_, outputStyle_, outputFreq_;
  bool adaptiveBlockSize_;
  std::string orthoType_, recycleMethod_;
  std::string impResScale_, expResScale_;
  std::string label_;

  /////////////////////////////////////////////////////////////////////////
  // Solver State Storage
  /////////////////////////////////////////////////////////////////////////
  //
  // The number of blocks and recycle blocks (m and k, respectively)
  int numBlocks_, recycledBlocks_;
  // Current size of recycled subspace
  int keff;
  //
  // Residual Vector
  Teuchos::RCP<MV> R_;
  //
  // Search Space
  Teuchos::RCP<MV> V_;
  //
  // Recycle subspace and its image under action of the operator
  Teuchos::RCP<MV> U_, C_;
  //
  // Updated recycle Space and its image under action of the operator
  Teuchos::RCP<MV> U1_, C1_;
  //
  // Storage used in constructing recycle space
  Teuchos::RCP<SDM > G_;
  Teuchos::RCP<SDM > H_;
  Teuchos::RCP<SDM > B_;
  Teuchos::RCP<SDM > PP_;
  Teuchos::RCP<SDM > HP_;
  std::vector<ScalarType> tau_;
  std::vector<ScalarType> work_;
  Teuchos::RCP<SDM > F_;
  std::vector<int> ipiv_;

  //! Timer for \c solve().
  Teuchos::RCP<Teuchos::Time> timerSolve_;

  //! Whether \c setParameters() successfully finished setting parameters.
  bool isSet_;

  //! Whether a loss of accuracy was detected during the solve.
  bool loaDetected_;

  //! Whether we have generated or regenerated a recycle space yet this solve.
  bool builtRecycleSpace_;
};

  //
  // Set default solver values
  //
  template<class ScalarType, class MV, class OP>
  const bool BlockGCRODRSolMgr<ScalarType,MV,OP>::adaptiveBlockSize_default_ = true;

  template<class ScalarType, class MV, class OP>
  const std::string BlockGCRODRSolMgr<ScalarType,MV,OP>::recycleMethod_default_ = "harmvecs";

  //
  // Method definitions
  //

  template<class ScalarType, class MV, class OP>
  BlockGCRODRSolMgr<ScalarType,MV,OP>::BlockGCRODRSolMgr() {
    init();
  }

  //Basic Constructor
  template<class ScalarType, class MV, class OP>
  BlockGCRODRSolMgr<ScalarType,MV,OP>::
  BlockGCRODRSolMgr(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                    const Teuchos::RCP<Teuchos::ParameterList> &pl ) {
    // Initialize local pointers to null, and initialize local
    // variables to default values.
    init ();

    TEUCHOS_TEST_FOR_EXCEPTION(problem.is_null(), std::invalid_argument,
      "Belos::BlockGCRODR constructor: The solver manager's constructor needs "
      "the linear problem argument 'problem' to be nonnull.");

    problem_ = problem;

    // Set the parameters using the list that was passed in.  If null, we defer initialization until a non-null list is set (by the
    // client calling setParameters(), or by calling solve() -- in either case, a null parameter list indicates that default parameters should be used).
    if (! pl.is_null())
      setParameters (pl);
  }

  template<class ScalarType, class MV, class OP>
  void BlockGCRODRSolMgr<ScalarType,MV,OP>::init() {
    adaptiveBlockSize_ = adaptiveBlockSize_default_;
    recycleMethod_ = recycleMethod_default_;
    isSet_ = false;
    loaDetected_ = false;
    builtRecycleSpace_ = false;
    keff = 0;//Effective Size of Recycle Space
    //The following are all RCP smart pointers to indicated matrices/vectors.
    //Some MATLAB notation used in comments.
    R_ = Teuchos::null;//The Block Residual
    V_ = Teuchos::null;//Block Arnoldi Vectors
    U_ = Teuchos::null;//Recycle Space
    C_ = Teuchos::null;//Image of U Under Action of Operator
    U1_ = Teuchos::null;//Newly Computed Recycle Space
    C1_ = Teuchos::null;//Image of Newly Computed Recycle Space
    PP_ = Teuchos::null;//Coordinates of New Recyc. Vectors in Augmented Space
    HP_ = Teuchos::null;//H_*PP_ or G_*PP_
    G_ = Teuchos::null;//G_ such that A*[U V(:,1:numBlocks_*blockSize_)] = [C V_]*G_
    F_ = Teuchos::null;//Upper Triangular portion of QR factorization for HP_
    H_ = Teuchos::null;//H_ such that A*V(:,1:numBlocks_*blockSize_) = V_*H_ + C_*C_'*V_
    B_ = Teuchos::null;//B_ = C_'*V_

    //THIS BLOCK OF CODE IS COMMENTED OUT AND PLACED ELSEWHERE IN THE CODE
/*
    //WE TEMPORARILY INITIALIZE STATUS TESTS HERE FOR TESTING PURPOSES, BUT THEY SHOULD BE
    //INITIALIZED SOMEWHERE ELSE, LIKE THE setParameters() FUNCTION

    //INSTANTIATE AND INITIALIZE TEST OBJECTS AS NEEDED
    if (maxIterTest_.is_null()){
      maxIterTest_ = rcp (new StatusTestMaxIters<ScalarType,MV,OP> (maxIters_));
    }
    //maxIterTest_->setMaxIters (maxIters_);

    //INSTANTIATE THE PRINTER
    if (printer_.is_null()) {
      printer_ = Teuchos::rcp (new OutputManager<ScalarType> (verbosity_, outputStream_));
    }

    if (ortho_.is_null()) // || changedOrthoType) %%%In other codes, this is also triggered if orthogonalization type changed {
      // Create orthogonalization manager.  This requires that the OutputManager (printer_) already be initialized.
      Teuchos::RCP<const Teuchos::ParameterList> orthoParams = orthoFactory_.getDefaultParameters (orthoType_);
      ortho_ = orthoFactory_.makeMatOrthoManager (orthoType_, Teuchos::null, printer_, label_, orthoParams);
    }

    // Convenience typedefs
    typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
    typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;

    if (impConvTest_.is_null()) {
      impConvTest_ = rcp (new StatusTestResNorm_t (convTol_));
      impConvTest_->defineScaleForm (convertStringToScaleType (impResScale_), Belos::TwoNorm);
      impConvTest_->setShowMaxResNormOnly( true );
    }

    if (expConvTest_.is_null()) {
      expConvTest_ = rcp (new StatusTestResNorm_t (convTol_));
      expConvTest_->defineResForm (StatusTestResNorm_t::Explicit, Belos::TwoNorm);
      expConvTest_->defineScaleForm (convertStringToScaleType (expResScale_), Belos::TwoNorm);
      expConvTest_->setShowMaxResNormOnly( true );
    }

    if (convTest_.is_null()) {
      convTest_ = rcp (new StatusTestCombo_t (StatusTestCombo_t::SEQ,  impConvTest_, expConvTest_));
    }

    sTest_ = rcp (new StatusTestCombo_t (StatusTestCombo_t::OR, maxIterTest_, convTest_));

    StatusTestOutputFactory<ScalarType,MV,OP> stoFactory (outputStyle_);
    outputTest_ = stoFactory.create (printer_, sTest_, outputFreq_,
    Passed+Failed+Undefined);

 */

  }

  //  This method requires the solver manager to return a string that describes itself.
  template<class ScalarType, class MV, class OP>
  std::string BlockGCRODRSolMgr<ScalarType,MV,OP>::description() const {
    std::ostringstream oss;
    oss << "Belos::BlockGCRODRSolMgr<" << SCT::name() << ", ...>";
    oss << "{";
    oss << "Ortho Type='"<<orthoType_ ;
    oss << ", Num Blocks=" <<numBlocks_;
    oss << ", Block Size=" <<blockSize_;
    oss << ", Num Recycle Blocks=" << recycledBlocks_;
    oss << ", Max Restarts=" << maxRestarts_;
    oss << "}";
    return oss.str();
  }

   template<class ScalarType, class MV, class OP>
   Teuchos::RCP<const Teuchos::ParameterList>
   BlockGCRODRSolMgr<ScalarType,MV,OP>::getValidParameters() const {
     using Teuchos::ParameterList;
     using Teuchos::parameterList;
     using Teuchos::RCP;
     using Teuchos::rcpFromRef;

     if (defaultParams_.is_null()) {
       RCP<ParameterList> pl = parameterList ();

       const MagnitudeType convTol = SMT::squareroot (SCT::magnitude (SCT::eps()));
       const int maxRestarts = 1000;
       const int maxIters = 5000;
       const int blockSize = 2;
       const int numBlocks = 100;
       const int numRecycledBlocks = 25;
       const int verbosity = Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails;
       const Belos::OutputType outputStyle = Belos::General;
       const int outputFreq = 1;
       RCP<std::ostream> outputStream = rcpFromRef (std::cout);
       const std::string impResScale ("Norm of Preconditioned Initial Residual");
       const std::string expResScale ("Norm of Initial Residual");
       const std::string timerLabel ("Belos");
       const std::string orthoType ("ICGS");
       RCP<const ParameterList> orthoParams = orthoFactory_.getDefaultParameters (orthoType);
       //const MagnitudeType orthoKappa = SCT::magnitude (SCT::eps());
       //
       // mfh 31 Dec 2011: Negative values of orthoKappa are ignored.
       // Setting this to a negative value by default ensures that
       // this parameter is only _not_ ignored if the user
       // specifically sets a valid value.
       const MagnitudeType orthoKappa = -SMT::one();

       // Set all the valid parameters and their default values.
       pl->set ("Convergence Tolerance", convTol,
                "The tolerance that the solver needs to achieve in order for "
                "the linear system(s) to be declared converged.  The meaning "
                "of this tolerance depends on the convergence test details.");
       pl->set("Maximum Restarts", maxRestarts,
               "The maximum number of restart cycles (not counting the first) "
               "allowed for each set of right-hand sides solved.");
       pl->set("Maximum Iterations", maxIters,
               "The maximum number of iterations allowed for each set of "
               "right-hand sides solved.");
       pl->set("Block Size", blockSize,
               "How many linear systems to solve at once.");
       pl->set("Num Blocks", numBlocks,
               "The maximum number of blocks allowed in the Krylov subspace "
               "for each set of right-hand sides solved.  This includes the "
               "initial block.  Total Krylov basis storage, not counting the "
               "recycled subspace, is \"Num Blocks\" times \"Block Size\".");
       pl->set("Num Recycled Blocks", numRecycledBlocks,
               "The maximum number of vectors in the recycled subspace." );
       pl->set("Verbosity", verbosity,
               "What type(s) of solver information should be written "
               "to the output stream.");
       pl->set("Output Style", outputStyle,
               "What style is used for the solver information to write "
               "to the output stream.");
       pl->set("Output Frequency", outputFreq,
               "How often convergence information should be written "
               "to the output stream.");
       pl->set("Output Stream", outputStream,
               "A reference-counted pointer to the output stream where all "
               "solver output is sent.");
       pl->set("Implicit Residual Scaling", impResScale,
               "The type of scaling used in the implicit residual convergence test.");
       pl->set("Explicit Residual Scaling", expResScale,
               "The type of scaling used in the explicit residual convergence test.");
       pl->set("Timer Label", timerLabel,
               "The string to use as a prefix for the timer labels.");
       pl->set("Orthogonalization", orthoType,
               "The orthogonalization method to use.  Valid options: " +
               orthoFactory_.validNamesString());
       pl->set ("Orthogonalization Parameters", *orthoParams,
                "Sublist of parameters specific to the orthogonalization method to use.");
       pl->set("Orthogonalization Constant", orthoKappa,
               "When using DGKS orthogonalization: the \"depTol\" constant, used "
               "to determine whether another step of classical Gram-Schmidt is "
               "necessary.  Otherwise ignored.  Nonpositive values are ignored.");
       defaultParams_ = pl;
     }
     return defaultParams_;
   }

   template<class ScalarType, class MV, class OP>
   void
   BlockGCRODRSolMgr<ScalarType,MV,OP>::
   setParameters (const Teuchos::RCP<Teuchos::ParameterList> &params) {
     using Teuchos::isParameterType;
     using Teuchos::getParameter;
     using Teuchos::null;
     using Teuchos::ParameterList;
     using Teuchos::parameterList;
     using Teuchos::RCP;
     using Teuchos::rcp;
     using Teuchos::rcp_dynamic_cast;
     using Teuchos::rcpFromRef;
     using Teuchos::Exceptions::InvalidParameter;
     using Teuchos::Exceptions::InvalidParameterName;
     using Teuchos::Exceptions::InvalidParameterType;

     // The default parameter list contains all parameters that BlockGCRODRSolMgr understands, and none that it doesn't
     RCP<const ParameterList> defaultParams = getValidParameters();

     if (params.is_null()) {
       // Users really should supply a non-null input ParameterList,
       // so that we can fill it in.  However, if they really did give
       // us a null list, be nice and use default parameters.  In this
       // case, we don't have to do any validation.
       params_ = parameterList (*defaultParams);
     }
     else {
       // Validate the user's parameter list and fill in any missing parameters with default values.
       //
       // Currently, validation is quite strict.  The following line
       // will throw an exception for misspelled or extra parameters.
       // However, it will _not_ throw an exception if there are
       // missing parameters: these will be filled in with their
       // default values.  There is additional discussion of other
       // validation strategies in the comments of this function for
       // Belos::GCRODRSolMgr.
       params->validateParametersAndSetDefaults (*defaultParams);
       // No side effects on the solver until after validation passes.
       params_ = params;
     }

     //
     // Validate and set parameters relating to the Krylov subspace
     // dimensions and the number of blocks.
     //
     // We try to read and validate all these parameters' values
     // before setting them in the solver.  This is because the
     // validity of each may depend on the others' values.  Our goal
     // is to avoid incorrect side effects, so we don't set the values
     // in the solver until we know they are right.
     //
     {
       const int maxRestarts = params_->get<int> ("Maximum Restarts");
       TEUCHOS_TEST_FOR_EXCEPTION(maxRestarts <= 0, std::invalid_argument,
         "Belos::BlockGCRODRSolMgr: The \"Maximum Restarts\" parameter "
         "must be nonnegative, but you specified a negative value of "
         << maxRestarts << ".");

       const int maxIters = params_->get<int> ("Maximum Iterations");
       TEUCHOS_TEST_FOR_EXCEPTION(maxIters <= 0, std::invalid_argument,
         "Belos::BlockGCRODRSolMgr: The \"Maximum Iterations\" parameter "
         "must be positive, but you specified a nonpositive value of "
         << maxIters << ".");

       const int numBlocks = params_->get<int> ("Num Blocks");
       TEUCHOS_TEST_FOR_EXCEPTION(numBlocks <= 0, std::invalid_argument,
         "Belos::BlockGCRODRSolMgr: The \"Num Blocks\" parameter must be "
         "positive, but you specified a nonpositive value of " << numBlocks
         << ".");

       const int blockSize = params_->get<int> ("Block Size");
       TEUCHOS_TEST_FOR_EXCEPTION(blockSize <= 0, std::invalid_argument,
         "Belos::BlockGCRODRSolMgr: The \"Block Size\" parameter must be "
         "positive, but you specified a nonpositive value of " << blockSize
         << ".");

       const int recycledBlocks = params_->get<int> ("Num Recycled Blocks");
       TEUCHOS_TEST_FOR_EXCEPTION(recycledBlocks <= 0, std::invalid_argument,
         "Belos::BlockGCRODRSolMgr: The \"Num Recycled Blocks\" parameter must "
         "be positive, but you specified a nonpositive value of "
         << recycledBlocks << ".");
       TEUCHOS_TEST_FOR_EXCEPTION(recycledBlocks >= numBlocks,
         std::invalid_argument, "Belos::BlockGCRODRSolMgr: The \"Num Recycled "
         "Blocks\" parameter must be less than the \"Num Blocks\" parameter, "
         "but you specified \"Num Recycled Blocks\" = " << recycledBlocks
         << " and \"Num Blocks\" = " << numBlocks << ".");

       // Now that we've validated the various dimension-related
       // parameters, we can set them in the solvers.
       maxRestarts_ = maxRestarts;
       maxIters_ = maxIters;
       numBlocks_ = numBlocks;
       blockSize_ = blockSize;
       recycledBlocks_ = recycledBlocks;
     }

     // Check to see if the timer label changed.  If it did, update it in
     // the parameter list, and create a new timer with that label (if
     // Belos was compiled with timers enabled).
     {
       std::string tempLabel = params_->get<std::string> ("Timer Label");
       const bool labelChanged = (tempLabel != label_);

#ifdef BELOS_TEUCHOS_TIME_MONITOR
       std::string solveLabel = label_ + ": BlockGCRODRSolMgr total solve time";
       if (timerSolve_.is_null()) {
         // We haven't created a timer yet.
         timerSolve_ = Teuchos::TimeMonitor::getNewCounter (solveLabel);
       } else if (labelChanged) {
         // We created a timer before with a different label.  In that
         // case, clear the old timer and create a new timer with the
         // new label.  Clearing old timers prevents them from piling
         // up, since they are stored statically, possibly without
         // checking for duplicates.
         Teuchos::TimeMonitor::clearCounter (solveLabel);
         timerSolve_ = Teuchos::TimeMonitor::getNewCounter (solveLabel);
       }
#endif // BELOS_TEUCHOS_TIME_MONITOR
     }

     // Check for a change in verbosity level.  Just in case, we allow
     // the type of this parameter to be either int or MsgType (an
     // enum).
     if (params_->isParameter ("Verbosity")) {
       if (isParameterType<int> (*params_, "Verbosity")) {
         verbosity_ = params_->get<int> ("Verbosity");
       }
       else {
         verbosity_ = (int) getParameter<MsgType> (*params_, "Verbosity");
       }
     }

     // Check for a change in output style.  Just in case, we allow
     // the type of this parameter to be either int or OutputType (an
     // enum).
     if (params_->isParameter ("Output Style")) {
       if (isParameterType<int> (*params_, "Output Style")) {
         outputStyle_ = params_->get<int> ("Output Style");
       }
       else {
         outputStyle_ = (int) getParameter<OutputType> (*params_, "Output Style");
       }

       // Currently, the only way we can update the output style is to
       // create a new output test.  This is because we must first
       // instantiate a new StatusTestOutputFactory with the new
       // style, and then use the factory to make a new output test.
       // Thus, we set outputTest_ to null for now, so we remember
       // later to create a new one.
       outputTest_ = null;
     }

     // Get the output stream for the output manager.
     //
     // It has been pointed out (mfh 28 Feb 2011 in GCRODRSolMgr code)
     // that including an RCP<std::ostream> parameter makes it
     // impossible to serialize the parameter list.  If you serialize
     // the list and read it back in from the serialized
     // representation, you might not even get a valid
     // RCP<std::ostream>, let alone the same output stream that you
     // serialized.
     //
     // However, existing Belos users depend on setting the output
     // stream in the parameter list.  We retain this behavior for
     // backwards compatibility.
     //
     // In case the output stream can't be read back in, we default to
     // stdout (std::cout), just to ensure reasonable behavior.
     if (params_->isParameter ("Output Stream")) {
       try {
         outputStream_ = getParameter<RCP<std::ostream> > (*params_, "Output Stream");
       }
       catch (InvalidParameter&) {
         outputStream_ = rcpFromRef (std::cout);
       }
       // We assume that a null output stream indicates that the user
       // doesn't want to print anything, so we replace it with a
       // "black hole" stream that prints nothing sent to it.  (We
       // can't use outputStream_ = Teuchos::null, since the output
       // manager assumes that operator<< works on its output stream.)
       if (outputStream_.is_null()) {
         outputStream_ = rcp (new Teuchos::oblackholestream);
       }
     }

     outputFreq_ = params_->get<int> ("Output Frequency");
     if (verbosity_ & Belos::StatusTestDetails) {
       // Update parameter in our output status test.
       if (! outputTest_.is_null()) {
         outputTest_->setOutputFrequency (outputFreq_);
       }
     }

     // Create output manager if we need to, using the verbosity level
     // and output stream that we fetched above.  We do this here because
     // instantiating an OrthoManager using OrthoManagerFactory requires
     // a valid OutputManager.
     if (printer_.is_null()) {
       printer_ = rcp (new OutputManager<ScalarType> (verbosity_, outputStream_));
     } else {
       printer_->setVerbosity (verbosity_);
       printer_->setOStream (outputStream_);
     }

     // Get the orthogonalization manager name ("Orthogonalization").
     //
     // Getting default values for the orthogonalization manager
     // parameters ("Orthogonalization Parameters") requires knowing the
     // orthogonalization manager name.  Save it for later, and also
     // record whether it's different than before.

     //bool changedOrthoType = false; // silence warning about not being used
     if (params_->isParameter ("Orthogonalization")) {
       const std::string& tempOrthoType =
         params_->get<std::string> ("Orthogonalization");
       // Ensure that the specified orthogonalization type is valid.
       if (! orthoFactory_.isValidName (tempOrthoType)) {
         std::ostringstream os;
         os << "Belos::BlockGCRODRSolMgr: Invalid orthogonalization name \""
            << tempOrthoType << "\".  The following are valid options "
            << "for the \"Orthogonalization\" name parameter: ";
         orthoFactory_.printValidNames (os);
         TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
       }
       if (tempOrthoType != orthoType_) {
         //changedOrthoType = true;
         orthoType_ = tempOrthoType;
       }
     }

     // Get any parameters for the orthogonalization
     // ("Orthogonalization Parameters").  If not supplied, the
     // orthogonalization manager factory will supply default values.
     //
     // NOTE (mfh 12 Jan 2011, summer 2011, 31 Dec 2011) For backwards
     // compatibility with other Belos GMRES-type solvers, if params
     // has an "Orthogonalization Constant" parameter and the DGKS
     // orthogonalization manager is to be used, the value of this
     // parameter will override DGKS's "depTol" parameter.
     //
     // Users must supply the orthogonalization manager parameters as
     // a sublist (supplying it as an RCP<ParameterList> would make
     // the resulting parameter list not serializable).
     RCP<ParameterList> orthoParams = sublist (params, "Orthogonalization Parameters", true);
     TEUCHOS_TEST_FOR_EXCEPTION(orthoParams.is_null(), std::logic_error,
                                "Failed to get orthogonalization parameters.  "
                                "Please report this bug to the Belos developers.");

     // Check if the desired orthogonalization method changed, or if
     // the orthogonalization manager has not yet been instantiated.
     // If either is the case, instantiate a new MatOrthoManager
     // subclass instance corresponding to the desired
     // orthogonalization method.  We've already fetched the
     // orthogonalization method name (orthoType_) and its parameters
     // (orthoParams) above.
     //
     // As suggested (by mfh 12 Jan 2011 in Belos::GCRODRSolMgr) In
     // order to ensure parameter changes get propagated to the
     // orthomanager we simply reinstantiate the OrthoManager every
     // time, whether or not the orthogonalization method name or
     // parameters have changed.  This is not efficient for some
     // orthogonalization managers that keep a lot of internal state.
     // A more general way to fix this inefficiency would be to
     //
     // 1. make each orthogonalization implement
     //    Teuchos::ParameterListAcceptor, and
     //
     // 2. make each orthogonalization implement a way to set the
     //    OutputManager instance and timer label.
     //
     // Create orthogonalization manager.  This requires that the
     // OutputManager (printer_) already be initialized.
     ortho_ = orthoFactory_.makeMatOrthoManager (orthoType_, null, printer_,
                                                 label_, orthoParams);

     // The DGKS orthogonalization accepts a "Orthogonalization
     // Constant" parameter (also called kappa in the code, but not in
     // the parameter list).  If its value is provided in the given
     // parameter list and its value is positive, use it.  Ignore
     // negative values.
     //
     // The "Orthogonalization Constant" parameter is retained only
     // for backwards compatibility.

     //bool gotValidOrthoKappa = false; // silence warning about not being used
     if (params_->isParameter ("Orthogonalization Constant")) {
       const MagnitudeType orthoKappa =
         params_->get<MagnitudeType> ("Orthogonalization Constant");
       if (orthoKappa > 0) {
         orthoKappa_ = orthoKappa;
         //gotValidOrthoKappa = true;
         // Only DGKS currently accepts this parameter.
         if (orthoType_ == "DGKS" && ! ortho_.is_null()) {
           typedef DGKSOrthoManager<ScalarType, MV, OP> ortho_man_type;
           // This cast should always succeed; it's a bug otherwise.
           // (If the cast fails, then orthoType_ doesn't correspond
           // to the OrthoManager subclass instance that we think we
           // have, so we initialized the wrong subclass somehow.)
           rcp_dynamic_cast<ortho_man_type>(ortho_)->setDepTol (orthoKappa_);
         }
       }
     }

     // Convergence
     typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
     typedef Belos::StatusTestGenResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;

     // Check for convergence tolerance
     convTol_ = params_->get<MagnitudeType> ("Convergence Tolerance");
     if (! impConvTest_.is_null()) {
       impConvTest_->setTolerance (convTol_);
     }
     if (! expConvTest_.is_null()) {
       expConvTest_->setTolerance (convTol_);
     }

     // Check for a change in scaling, if so we need to build new residual tests.
     if (params_->isParameter ("Implicit Residual Scaling")) {
       std::string tempImpResScale =
         getParameter<std::string> (*params_, "Implicit Residual Scaling");

       // Only update the scaling if it's different.
       if (impResScale_ != tempImpResScale) {
         ScaleType impResScaleType = convertStringToScaleType (tempImpResScale);
         impResScale_ = tempImpResScale;

         if (! impConvTest_.is_null()) {
           try {
             impConvTest_->defineScaleForm (impResScaleType, Belos::TwoNorm);
           }
           catch (StatusTestError&) {
             // Delete the convergence test so it gets constructed again.
             impConvTest_ = null;
             convTest_ = null;
           }
         }
       }
     }

     if (params_->isParameter("Explicit Residual Scaling")) {
       std::string tempExpResScale =
         getParameter<std::string> (*params_, "Explicit Residual Scaling");

       // Only update the scaling if it's different.
       if (expResScale_ != tempExpResScale) {
         ScaleType expResScaleType = convertStringToScaleType (tempExpResScale);
         expResScale_ = tempExpResScale;

         if (! expConvTest_.is_null()) {
           try {
             expConvTest_->defineScaleForm (expResScaleType, Belos::TwoNorm);
           }
           catch (StatusTestError&) {
             // Delete the convergence test so it gets constructed again.
             expConvTest_ = null;
             convTest_ = null;
           }
         }
       }
     }

     //
     // Create iteration stopping criteria ("status tests") if we need
     // to, by combining three different stopping criteria.
     //
     // First, construct maximum-number-of-iterations stopping criterion.
     if (maxIterTest_.is_null()) {
       maxIterTest_ = rcp (new StatusTestMaxIters<ScalarType,MV,OP> (maxIters_));
     } else {
       maxIterTest_->setMaxIters (maxIters_);
     }

     // Implicit residual test, using the native residual to determine if
     // convergence was achieved.
     if (impConvTest_.is_null()) {
       impConvTest_ = rcp (new StatusTestResNorm_t (convTol_));
       impConvTest_->defineScaleForm (convertStringToScaleType (impResScale_),
                                      Belos::TwoNorm);
     }

     // Explicit residual test once the native residual is below the tolerance
     if (expConvTest_.is_null()) {
       expConvTest_ = rcp (new StatusTestResNorm_t (convTol_));
       expConvTest_->defineResForm (StatusTestResNorm_t::Explicit, Belos::TwoNorm);
       expConvTest_->defineScaleForm (convertStringToScaleType (expResScale_),
                                      Belos::TwoNorm);
     }
     // Convergence test first tests the implicit residual, then the
     // explicit residual if the implicit residual test passes.
     if (convTest_.is_null()) {
       convTest_ = rcp (new StatusTestCombo_t (StatusTestCombo_t::SEQ,
                                               impConvTest_,
                                               expConvTest_));
     }
     // Construct the complete iteration stopping criterion:
     //
     // "Stop iterating if the maximum number of iterations has been
     // reached, or if the convergence test passes."
     sTest_ = rcp (new StatusTestCombo_t (StatusTestCombo_t::OR,
                                          maxIterTest_, convTest_));
     // Create the status test output class.
     // This class manages and formats the output from the status test.
     StatusTestOutputFactory<ScalarType,MV,OP> stoFactory (outputStyle_);
     outputTest_ = stoFactory.create (printer_, sTest_, outputFreq_,
                                      Passed+Failed+Undefined);

     // Set the solver string for the output test
     std::string solverDesc = "Block GCRODR ";
     outputTest_->setSolverDesc (solverDesc);

     // Inform the solver manager that the current parameters were set.
     isSet_ = true;
   }

  // initializeStateStorage.
  template<class ScalarType, class MV, class OP>
  void
  BlockGCRODRSolMgr<ScalarType,MV,OP>::initializeStateStorage()
  {

    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

    // Check if there is any multivector to clone from.
    Teuchos::RCP<const MV> rhsMV = problem_->getRHS();

    //The Dimension of the Krylov Subspace
    int KrylSpaDim = (numBlocks_ - 1) * blockSize_;

    //Number of columns in [U_ V_(:,1:KrylSpaDim)]
    int augSpaDim = KrylSpaDim + recycledBlocks_ + 1;// + 1 is for possible extra recycle vector

    //Number of columns in [C_ V_]
    int augSpaImgDim = KrylSpaDim + blockSize_ + recycledBlocks_+1;

    //TEMPORARY SKELETON DEFINITION OF THIS FUNCTION TO GET THINGS WORKING
    //NOT EVERYTHING IS INITIALIZE CORRECTLY YET.

    //INITIALIZE RECYCLE SPACE VARIABLES HERE

    //WE DO NOT ALLOCATE V HERE IN THIS SOLVERMANAGER. If no recycle space exists, then block_gmres_iter
    //will allocated V for us.  If a recycle space already exists, then we will allocate V after updating the
    //recycle space for the current problem.
    // If the block Krylov subspace has not been initialized before, generate it using the RHS from lp_.
    /*if (V_ == Teuchos::null) {
      V_ = MVT::Clone( *rhsMV, (numBlocks_+1)*blockSize_ );
      }
      else{
      // Generate V_ by cloning itself ONLY if more space is needed.
      if (MVT::GetNumberVecs(*V_) < numBlocks_+1) {
      Teuchos::RCP<const MV> tmp = V_;
      V_ = MVT::Clone( *tmp, numBlocks_+1 );
      }
      }*/

    //INTITIALIZE SPACE FOR THE NEWLY COMPUTED RECYCLE SPACE VARIABLES HERE

    if (U_ == Teuchos::null) {
      U_ = MVT::Clone( *rhsMV, recycledBlocks_+1 );
    }
    else {
      // Generate U_ by cloning itself ONLY if more space is needed.
      if (MVT::GetNumberVecs(*U_) < recycledBlocks_+1) {
        Teuchos::RCP<const MV> tmp = U_;
        U_ = MVT::Clone( *tmp, recycledBlocks_+1 );
      }
    }

    // If the subspace has not been initialized before, generate it using the RHS from lp_.
    if (C_ == Teuchos::null) {
      C_ = MVT::Clone( *rhsMV, recycledBlocks_+1 );
    }
    else {
      // Generate C_ by cloning itself ONLY if more space is needed.
      if (MVT::GetNumberVecs(*C_) < recycledBlocks_+1) {
        Teuchos::RCP<const MV> tmp = C_;
        C_ = MVT::Clone( *tmp, recycledBlocks_+1 );
      }
    }

    // If the subspace has not been initialized before, generate it using the RHS from lp_.
    if (U1_ == Teuchos::null) {
      U1_ = MVT::Clone( *rhsMV, recycledBlocks_+1 );
    }
    else {
      // Generate U1_ by cloning itself ONLY if more space is needed.
      if (MVT::GetNumberVecs(*U1_) < recycledBlocks_+1) {
        Teuchos::RCP<const MV> tmp = U1_;
        U1_ = MVT::Clone( *tmp, recycledBlocks_+1 );
      }
    }

    // If the subspace has not been initialized before, generate it using the RHS from lp_.
    if (C1_ == Teuchos::null) {
      C1_ = MVT::Clone( *rhsMV, recycledBlocks_+1 );
    }
    else {
      // Generate C1_ by cloning itself ONLY if more space is needed.
      if (MVT::GetNumberVecs(*U1_) < recycledBlocks_+1) {
        Teuchos::RCP<const MV> tmp = C1_;
        C1_ = MVT::Clone( *tmp, recycledBlocks_+1 );
      }
    }

    // Generate R_ only if it doesn't exist
    if (R_ == Teuchos::null){
      R_ = MVT::Clone( *rhsMV, blockSize_ );
    }

    //INITIALIZE SOME WORK VARIABLES

    // Generate G_ only if it doesn't exist, otherwise resize it.
    if (G_ == Teuchos::null){
      G_ = Teuchos::rcp( new SDM( augSpaImgDim, augSpaDim ) );
    }
    else{
      if ( (G_->numRows() != augSpaImgDim) || (G_->numCols() != augSpaDim) )
        {
          G_->reshape( augSpaImgDim, augSpaDim );
        }
      G_->putScalar(zero);
    }

    // Generate H_ only if it doesn't exist by pointing it to a view of G_.
    if (H_ == Teuchos::null){
      H_ = Teuchos::rcp (new SDM ( Teuchos::View, *G_, KrylSpaDim + blockSize_, KrylSpaDim, recycledBlocks_+1 ,recycledBlocks_+1 ) );
    }

    // Generate F_ only if it doesn't exist, otherwise resize it.
    if (F_ == Teuchos::null){
      F_ = Teuchos::rcp( new SDM( recycledBlocks_+1, recycledBlocks_+1 ) );
    }
    else {
      if ( (F_->numRows() != recycledBlocks_+1) || (F_->numCols() != recycledBlocks_+1) ){
        F_->reshape( recycledBlocks_+1, recycledBlocks_+1 );
      }
    }
    F_->putScalar(zero);

    // Generate PP_ only if it doesn't exist, otherwise resize it.
    if (PP_ == Teuchos::null){
      PP_ = Teuchos::rcp( new SDM( augSpaImgDim, recycledBlocks_+1 ) );
    }
    else{
      if ( (PP_->numRows() != augSpaImgDim) || (PP_->numCols() != recycledBlocks_+1) ){
        PP_->reshape( augSpaImgDim, recycledBlocks_+1 );
      }
    }

    // Generate HP_ only if it doesn't exist, otherwise resize it.
    if (HP_ == Teuchos::null)
      HP_ = Teuchos::rcp( new SDM( augSpaImgDim, augSpaDim ) );
    else{
      if ( (HP_->numRows() != augSpaImgDim) || (HP_->numCols() != augSpaDim) ){
        HP_->reshape( augSpaImgDim, augSpaDim );
      }
    }

    // Size of tau_ will change during computation, so just be sure it starts with appropriate size
    tau_.resize(recycledBlocks_+1);

    // Size of work_ will change during computation, so just be sure it starts with appropriate size
    work_.resize(recycledBlocks_+1);

    // Size of ipiv_ will change during computation, so just be sure it starts with appropriate size
    ipiv_.resize(recycledBlocks_+1);

  }

template<class ScalarType, class MV, class OP>
void BlockGCRODRSolMgr<ScalarType,MV,OP>::buildRecycleSpaceKryl(int& keff, Teuchos::RCP<BlockGmresIter<ScalarType,MV,OP> > block_gmres_iter){

  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  int p = block_gmres_iter->getState().curDim; //Dimension of the Krylov space generated
  std::vector<int> index(keff);//we use this to index certain columns of U, C, and V to
  //get views into pieces of these matrices.

  //GET CORRECT PIECE OF MATRIX H APPROPRIATE TO SIZE OF KRYLOV SUBSPACE
  SDM HH(Teuchos::Copy, *H_, p+blockSize_, p);
  if(recycledBlocks_ >= p + blockSize_){//keep whole block Krylov subspace
    //IF THIS HAS HAPPENED, THIS MEANS WE CONVERGED DURING THIS CYCLE
    //THEREFORE, WE DO NOT CONSTRUCT C = A*U;
    index.resize(p);
    for (int ii=0; ii < p; ++ii) index[ii] = ii;
    Teuchos::RCP<MV> Utmp  = MVT::CloneViewNonConst( *U_, index );
    MVT::SetBlock(*V_, index, *Utmp);
    keff = p;
  }
  else{ //use a subspace selection method to get recycle space
    int info = 0;
    Teuchos::RCP<SDM > PPtmp = Teuchos::rcp (new SDM ( Teuchos::View, *PP_, p, recycledBlocks_+1 ) );
    if(recycleMethod_ == "harmvecs"){
      keff = getHarmonicVecsKryl(p, HH, *PPtmp);
      printer_->stream(Debug) << "keff = " << keff << std::endl;
    }
// Hereafter, only keff columns of PP are needed
PPtmp = Teuchos::rcp (new SDM ( Teuchos::View, *PP_, p, keff ) );
// Now get views into C, U, V
index.resize(keff);
for (int ii=0; ii<keff; ++ii) index[ii] = ii;
Teuchos::RCP<MV> Ctmp  = MVT::CloneViewNonConst( *C_, index );
Teuchos::RCP<MV> Utmp  = MVT::CloneViewNonConst( *U_, index );
Teuchos::RCP<MV> U1tmp = MVT::CloneViewNonConst( *U1_, index );
index.resize(p);
for (int ii=0; ii < p; ++ii) index[ii] = ii;
Teuchos::RCP<const MV> Vtmp = MVT::CloneView( *V_, index );

// Form U (the subspace to recycle)
// U = newstate.V(:,1:p) * PP;
MVT::MvTimesMatAddMv( one, *Vtmp, *PPtmp, zero, *U1tmp );

// Step #1: Form HP = H*P
SDM HPtmp( Teuchos::View, *HP_, p+blockSize_, keff );
HPtmp.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, *H_, *PPtmp, zero );
// Step #1.5: Perform workspace size query for QR factorization of HP
int lwork = -1;
tau_.resize(keff);
lapack.GEQRF(HPtmp.numRows(),HPtmp.numCols(),HPtmp.values(),HPtmp.stride(),&tau_[0],&work_[0],lwork,&info);
TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure, "Belos::BlockGCRODRSolMgr::solve(): LAPACK _GEQRF failed to compute a workspace size.");

// Step #2: Compute QR factorization of HP
 lwork = static_cast<int> (Teuchos::ScalarTraits<ScalarType>::magnitude (work_[0]));
work_.resize(lwork);
lapack.GEQRF(HPtmp.numRows(),HPtmp.numCols(),HPtmp.values(),HPtmp.stride(),&tau_[0],&work_[0],lwork,&info);
TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure,  "Belos::BlockGCRODRSolMgr::solve(): LAPACK _GEQRF failed to compute a QR factorization.");

// Step #3: Explicitly construct Q and R factors
// The upper triangular part of HP is copied into R and HP becomes Q.
SDM Rtmp( Teuchos::View, *F_, keff, keff );
for(int ii=0;ii<keff;ii++) { for(int jj=ii;jj<keff;jj++) Rtmp(ii,jj) = HPtmp(ii,jj); }
//lapack.ORGQR(HPtmp.numRows(),HPtmp.numCols(),HPtmp.numCols(),HPtmp.values(),HPtmp.stride(),&tau_[0],&work_[0],lwork,&info);
lapack.UNGQR(HPtmp.numRows(),HPtmp.numCols(),HPtmp.numCols(),HPtmp.values(),HPtmp.stride(),&tau_[0],&work_[0],lwork,&info);
TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure, "Belos::BlockGCRODRSolMgr::solve(): LAPACK _UNGQR failed to construct the Q factor.");
                // Now we have [Q,R] = qr(H*P)

                // Now compute C = V(:,1:p+blockSize_) * Q
                index.resize( p+blockSize_ );
                for (int ii=0; ii < (p+blockSize_); ++ii) { index[ii] = ii; }
                Vtmp = MVT::CloneView( *V_, index ); // need new view into V (p+blockSize_ vectors now; needed p above)
                MVT::MvTimesMatAddMv( one, *Vtmp, HPtmp, zero, *Ctmp );

                // Finally, compute U = U*R^{-1}.
                // This unfortuntely requires me to form R^{-1} explicitly and execute U = U * R^{-1}, as
                // backsolve capabilities don't exist in the Belos::MultiVec class

// Step #1: First, compute LU factorization of R
ipiv_.resize(Rtmp.numRows());
lapack.GETRF(Rtmp.numRows(),Rtmp.numCols(),Rtmp.values(),Rtmp.stride(),&ipiv_[0],&info);
TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure, "Belos::GCRODRSolMgr::solve(): LAPACK _GETRF failed to compute an LU factorization.");
// Step #2: Form inv(R)
lwork = Rtmp.numRows();
work_.resize(lwork);
lapack.GETRI(Rtmp.numRows(),Rtmp.values(),Rtmp.stride(),&ipiv_[0],&work_[0],lwork,&info);
//TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure,  "Belos::GCRODRSolMgr::solve(): LAPACK _GETRI failed to invert triangular matrix.");
// Step #3: Let U = U * R^{-1}
MVT::MvTimesMatAddMv( one, *U1tmp, Rtmp, zero, *Utmp );

} //end else from if(recycledBlocks_ >= p + 1)
return;
} // end buildRecycleSpaceKryl defnition

template<class ScalarType, class MV, class OP>
void BlockGCRODRSolMgr<ScalarType,MV,OP>::buildRecycleSpaceAugKryl(Teuchos::RCP<BlockGCRODRIter<ScalarType,MV,OP> > block_gcrodr_iter){
  const MagnitudeType one = Teuchos::ScalarTraits<ScalarType>::one();
  const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  std::vector<MagnitudeType> d(keff);
  std::vector<ScalarType> dscalar(keff);
  std::vector<int> index(numBlocks_+1);

  // Get the state
  BlockGCRODRIterState<ScalarType,MV> oldState = block_gcrodr_iter->getState();
  int p = oldState.curDim;

  // insufficient new information to update recycle space
  if (p<1) return;

  if(recycledBlocks_ >= keff + p){ //we add new Krylov vectors to existing recycle space
     // If this has happened, we have converged during this cycle. Therefore, do not construct C = A*C;
     index.resize(p);
     for (int ii=0; ii < p; ++ii) { index[ii] = keff+ii; } //get a view after current reycle vectors
     Teuchos::RCP<MV> Utmp  = MVT::CloneViewNonConst( *U_, index );
     for (int ii=0; ii < p; ++ii) { index[ii] =ii; }
     MVT::SetBlock(*V_, index, *Utmp);
     keff += p;
  }

  // Take the norm of the recycled vectors
  {
    index.resize(keff);
    for (int ii=0; ii<keff; ++ii) { index[ii] = ii; }
    Teuchos::RCP<MV> Utmp  = MVT::CloneViewNonConst( *U_, index );
    d.resize(keff);
    dscalar.resize(keff);
    MVT::MvNorm( *Utmp, d );
    for (int i=0; i<keff; ++i) {
      d[i] = one / d[i];
      dscalar[i] = (ScalarType)d[i];
    }
    MVT::MvScale( *Utmp, dscalar );
  }

  // Get view into current "full" upper Hessnburg matrix
  // note that p describes the dimension of the iter+1 block Krylov space so we have to adjust how we use it to index Gtmp
  Teuchos::RCP<SDM> Gtmp = Teuchos::rcp( new SDM( Teuchos::View, *G_, p+keff, p+keff-blockSize_ ) );

  // Insert D into the leading keff x keff  block of G
  for (int i=0; i<keff; ++i)
    (*Gtmp)(i,i) = d[i];

  // Compute the harmoic Ritz pairs for the generalized eigenproblem
  // getHarmonicVecsKryl assumes PP has recycledBlocks_+1 columns available
  // See previous block of comments for why we subtract p-blockSize_
  int keff_new;
  {
    SDM PPtmp( Teuchos::View, *PP_, p+keff-blockSize_, recycledBlocks_+1 );
    keff_new = getHarmonicVecsAugKryl( keff, p-blockSize_, *Gtmp, oldState.V, PPtmp );
  }

  // Code to form new U, C
  // U = [U V(:,1:p)] * P; (in two steps)

  // U(:,1:keff) = matmul(U(:,1:keff_old),PP(1:keff_old,1:keff)) (step 1)
  Teuchos::RCP<MV> U1tmp;
  {
    index.resize( keff );
    for (int ii=0; ii<keff; ++ii) { index[ii] = ii; }
    Teuchos::RCP<const MV> Utmp  = MVT::CloneView( *U_,  index );
    index.resize( keff_new );
    for (int ii=0; ii<keff_new; ++ii) { index[ii] = ii; }
    U1tmp  = MVT::CloneViewNonConst( *U1_,  index );
    SDM PPtmp( Teuchos::View, *PP_, keff, keff_new );
    MVT::MvTimesMatAddMv( one, *Utmp, PPtmp, zero, *U1tmp );
  }

  // U(:,1:keff) = U(:,1:keff) + matmul(V(:,1:m-k),PP(keff_old+1:m-k+keff_old,1:keff)) (step 2)
  {
    index.resize(p-blockSize_);
    for (int ii=0; ii < p-blockSize_; ii++) { index[ii] = ii; }
    Teuchos::RCP<const MV> Vtmp = MVT::CloneView( *V_, index );
    SDM PPtmp( Teuchos::View, *PP_, p-blockSize_, keff_new, keff );
    MVT::MvTimesMatAddMv( one, *Vtmp, PPtmp, one, *U1tmp );
  }

  // Form GP = G*P
  SDM HPtmp( Teuchos::View, *HP_, p+keff, keff_new );
  {
    SDM PPtmp( Teuchos::View, *PP_, p-blockSize_+keff, keff_new );
    HPtmp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,one,*Gtmp,PPtmp,zero);
  }

  // Workspace size query for QR factorization HP= QF (the worksize will be placed in work_[0])
  int info = 0, lwork = -1;
  tau_.resize(keff_new);
  lapack.GEQRF(HPtmp.numRows(),HPtmp.numCols(),HPtmp.values(),HPtmp.stride(),&tau_[0],&work_[0],lwork,&info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0,BlockGCRODRSolMgrLAPACKFailure,"Belos::BlockGCRODRSolMgr::solve(): LAPACK _GEQRF failed to compute a workspace size.");

  lwork = static_cast<int> (Teuchos::ScalarTraits<ScalarType>::real (work_[0]));
  work_.resize(lwork);
  lapack.GEQRF(HPtmp.numRows(),HPtmp.numCols(),HPtmp.values(),HPtmp.stride(),&tau_[0],&work_[0],lwork,&info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0,BlockGCRODRSolMgrLAPACKFailure,"Belos::BlockGCRODRSolMgr::solve(): LAPACK _GEQRF failed to compute a QR factorization.");

  // Explicitly construct Q and F factors
  // NOTE:  The upper triangular part of HP is copied into F and HP becomes Q.
  SDM Ftmp( Teuchos::View, *F_, keff_new, keff_new );
  for(int i=0;i<keff_new;i++) { for(int j=i;j<keff_new;j++) Ftmp(i,j) = HPtmp(i,j); }
  //lapack.ORGQR(HPtmp.numRows(),HPtmp.numCols(),HPtmp.numCols(),HPtmp.values(),HPtmp.stride(),&tau_[0],&work_[0],lwork,&info);
  lapack.UNGQR(HPtmp.numRows(),HPtmp.numCols(),HPtmp.numCols(),HPtmp.values(),HPtmp.stride(),&tau_[0],&work_[0],lwork,&info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0,BlockGCRODRSolMgrLAPACKFailure,"Belos::BlockGCRODRSolMgr::solve(): LAPACK _UNGQR failed to construct the Q factor.");

  // Form orthonormalized C and adjust U accordingly so that C = A*U
  // C = [C V] * Q;

  // C(:,1:keff) = matmul(C(:,1:keff_old),QQ(1:keff_old,1:keff))
  {
    Teuchos::RCP<MV> C1tmp;
    {
      index.resize(keff);
      for (int i=0; i < keff; i++) { index[i] = i; }
      Teuchos::RCP<const MV> Ctmp  = MVT::CloneView( *C_,  index );
      index.resize(keff_new);
      for (int i=0; i < keff_new; i++) { index[i] = i; }
      C1tmp  = MVT::CloneViewNonConst( *C1_,  index );
      SDM PPtmp( Teuchos::View, *HP_, keff, keff_new );
      MVT::MvTimesMatAddMv( one, *Ctmp, PPtmp, zero, *C1tmp );
    }
    // Now compute C += V(:,1:p+1) * Q
    {
      index.resize( p );
      for (int i=0; i < p; ++i) { index[i] = i; }
      Teuchos::RCP<const MV> Vtmp = MVT::CloneView( *V_, index );
      SDM PPtmp( Teuchos::View, *HP_, p, keff_new, keff, 0 );
      MVT::MvTimesMatAddMv( one, *Vtmp, PPtmp, one, *C1tmp );
    }
  }

  // C_ = C1_; (via a swap)
  std::swap(C_, C1_);

  // Finally, compute U_ = U_*R^{-1}
  // First, compute LU factorization of R
  ipiv_.resize(Ftmp.numRows());
  lapack.GETRF(Ftmp.numRows(),Ftmp.numCols(),Ftmp.values(),Ftmp.stride(),&ipiv_[0],&info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0,BlockGCRODRSolMgrLAPACKFailure,"Belos::BlockGCRODRSolMgr::solve(): LAPACK _GETRF failed to compute an LU factorization.");

  // Now, form inv(R)
  lwork = Ftmp.numRows();
  work_.resize(lwork);
  lapack.GETRI(Ftmp.numRows(),Ftmp.values(),Ftmp.stride(),&ipiv_[0],&work_[0],lwork,&info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure,"Belos::BlockGCRODRSolMgr::solve(): LAPACK _GETRI failed to compute an LU factorization.");

  {
    index.resize(keff_new);
    for (int i=0; i < keff_new; i++) { index[i] = i; }
    Teuchos::RCP<MV> Utmp  = MVT::CloneViewNonConst( *U_,  index );
    MVT::MvTimesMatAddMv( one, *U1tmp, Ftmp, zero, *Utmp );
  }

  // Set the current number of recycled blocks and subspace
  // dimension with the Block GCRO-DR iteration.
  if (keff != keff_new) {
    keff = keff_new;
    block_gcrodr_iter->setSize( keff, numBlocks_ );
    // Important to zero this out before next cyle
    SDM b1( Teuchos::View, *G_, recycledBlocks_+2, 1, 0, recycledBlocks_ );
    b1.putScalar(zero);
  }

} //end buildRecycleSpaceAugKryl definition

template<class ScalarType, class MV, class OP>
int BlockGCRODRSolMgr<ScalarType,MV,OP>::getHarmonicVecsAugKryl(int keff, int m, const SDM& GG, const Teuchos::RCP<const MV>& VV, SDM& PP){
  int i, j;
  int m2 = GG.numCols();
  bool xtraVec = false;
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
  std::vector<int> index;

  // Real and imaginary eigenvalue components
  std::vector<MagnitudeType> wr(m2), wi(m2);

  // Magnitude of harmonic Ritz values
  std::vector<MagnitudeType> w(m2);

  // Real and imaginary (right) eigenvectors; Don't zero out matrix when constructing
  SDM vr(m2,m2,false);

  // Sorted order of harmonic Ritz values
  std::vector<int> iperm(m2);

  // Set flag indicating recycle space has been generated this solve
  builtRecycleSpace_ = true;

  // Form matrices for generalized eigenproblem

  // B = G' * G; Don't zero out matrix when constructing
  SDM B(m2,m2,false);
  B.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,one,GG,GG,zero);

  // A_tmp = | C'*U        0 |
  //         | V_{m+1}'*U  I |
  SDM A_tmp( keff+m+blockSize_, keff+m );


  // A_tmp(1:keff,1:keff) = C' * U;
  index.resize(keff);
  for (int i=0; i<keff; ++i) { index[i] = i; }
  Teuchos::RCP<const MV> Ctmp  = MVT::CloneView( *C_, index );
  Teuchos::RCP<const MV> Utmp  = MVT::CloneView( *U_, index );
  SDM A11( Teuchos::View, A_tmp, keff, keff );
  MVT::MvTransMv( one, *Ctmp, *Utmp, A11 );

  // A_tmp(keff+1:m-k+keff+1,1:keff) = V' * U;
  SDM A21( Teuchos::View, A_tmp, m+blockSize_, keff, keff );
  index.resize(m+blockSize_);
  for (i=0; i < m+blockSize_; i++) { index[i] = i; }
  Teuchos::RCP<const MV> Vp = MVT::CloneView( *VV, index );
  MVT::MvTransMv( one, *Vp, *Utmp, A21 );

  // A_tmp(keff+1:m-k+keff,keff+1:m-k+keff) = eye(m-k);
  for( i=keff; i<keff+m; i++)
    A_tmp(i,i) = one;

  // A = G' * A_tmp;
  SDM A( m2, A_tmp.numCols() );
  A.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, one, GG, A_tmp, zero );

  // Compute k smallest harmonic Ritz pairs
  // SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB,
  //                   ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO,
  //                   IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE,
  //                   RCONDV, WORK, LWORK, IWORK, BWORK, INFO )
  // MLP: 'SCALING' in DGGEVX generates incorrect eigenvalues. Therefore, only permuting
  char balanc='P', jobvl='N', jobvr='V', sense='N';
  int ld = A.numRows();
  int lwork = 6*ld;
  int ldvl = ld, ldvr = ld;
  int info = 0,ilo = 0,ihi = 0;
  MagnitudeType abnrm = 0.0, bbnrm = 0.0;
  ScalarType *vl = 0; // This is never referenced by dggevx if jobvl == 'N'
  std::vector<ScalarType> beta(ld);
  std::vector<ScalarType> work(lwork);
  std::vector<MagnitudeType> rwork(lwork);
  std::vector<MagnitudeType> lscale(ld), rscale(ld);
  std::vector<MagnitudeType> rconde(ld), rcondv(ld);
  std::vector<int> iwork(ld+6);
  int *bwork = 0; // If sense == 'N', bwork is never referenced
  //lapack.GGEVX(balanc, jobvl, jobvr, sense, ld, A.values(), ld, B.values(), ld, &wr[0], &wi[0],
  //             &beta[0], vl, ldvl, vr.values(), ldvr, &ilo, &ihi, &lscale[0], &rscale[0],
  //             &abnrm, &bbnrm, &rconde[0], &rcondv[0], &work[0], lwork, &iwork[0], bwork, &info);
  lapack.GGEVX(balanc, jobvl, jobvr, sense, ld, A.values(), ld, B.values(), ld, &wr[0], &wi[0],
               &beta[0], vl, ldvl, vr.values(), ldvr, &ilo, &ihi, &lscale[0], &rscale[0],
               &abnrm, &bbnrm, &rconde[0], &rcondv[0], &work[0], lwork, &rwork[0],
               &iwork[0], bwork, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure, "Belos::BlockGCRODRSolMgr::solve(): LAPACK GGEVX failed to compute eigensolutions.");

  // Construct magnitude of each harmonic Ritz value
  // NOTE : Forming alpha/beta *should* be okay here, given assumptions on construction of matrix pencil above
  for( i=0; i<ld; i++ ) // Construct magnitude of each harmonic Ritz value
    w[i] = Teuchos::ScalarTraits<MagnitudeType>::squareroot( wr[i]*wr[i] + wi[i]*wi[i] ) / std::abs( beta[i] );

  this->sort(w,ld,iperm);

  bool scalarTypeIsComplex = Teuchos::ScalarTraits<ScalarType>::isComplex;

  // Select recycledBlocks_ smallest eigenvectors
  for( i=0; i<recycledBlocks_; i++ )
    for( j=0; j<ld; j++ )
      PP(j,i) = vr(j,iperm[ld-recycledBlocks_+i]);

  if(scalarTypeIsComplex==false) {

    // Determine exact size for PP (i.e., determine if we need to store an additional vector)
    if (wi[iperm[ld-recycledBlocks_]] != 0.0) {
      int countImag = 0;
      for ( i=ld-recycledBlocks_; i<ld; i++ )
        if (wi[iperm[i]] != 0.0) countImag++;
      // Check to see if this count is even or odd:
      if (countImag % 2) xtraVec = true;
    }

    if (xtraVec) { // we need to store one more vector
      if (wi[iperm[ld-recycledBlocks_]] > 0.0) { // I picked the "real" component
        for( j=0; j<ld; j++ )                  // so get the "imag" component
          PP(j,recycledBlocks_) = vr(j,iperm[ld-recycledBlocks_]+1);
      }
      else {                  // I picked the "imag" component
        for( j=0; j<ld; j++ ) // so get the "real" component
          PP(j,recycledBlocks_) = vr(j,iperm[ld-recycledBlocks_]-1);
      }
    }

  }

  // Return whether we needed to store an additional vector
  if (xtraVec)
    return recycledBlocks_+1;
  else
    return recycledBlocks_;

} //end getHarmonicVecsAugKryl definition

template<class ScalarType, class MV, class OP>
int BlockGCRODRSolMgr<ScalarType,MV,OP>::getHarmonicVecsKryl(int m, const SDM& HH, SDM& PP){
  bool xtraVec = false;
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  // Real and imaginary eigenvalue components
  std::vector<MagnitudeType> wr(m), wi(m);

  // Real and imaginary (right) eigenvectors; Don't zero out matrix when constructing
  SDM vr(m,m,false);

  // Magnitude of harmonic Ritz values
  std::vector<MagnitudeType> w(m);

  // Sorted order of harmonic Ritz values, also used for DGEEV
  std::vector<int> iperm(m);

  // Output info
  int info = 0;

  // Set flag indicating recycle space has been generated this solve
  builtRecycleSpace_ = true;

  // Solve linear system:  H_m^{-H}*E_m where E_m is the last blockSize_ columns of the identity matrix
  SDM HHt( HH, Teuchos::TRANS );
  Teuchos::RCP<SDM> harmRitzMatrix = Teuchos::rcp( new SDM( m, blockSize_));

  //Initialize harmRitzMatrix as E_m
  for(int i=0; i<=blockSize_-1; i++) (*harmRitzMatrix)[blockSize_-1-i][harmRitzMatrix->numRows()-1-i] = 1;

  //compute harmRitzMatrix <- H_m^{-H}*E_m
  lapack.GESV(m, blockSize_, HHt.values(), HHt.stride(), &iperm[0], harmRitzMatrix->values(), harmRitzMatrix->stride(), &info);

  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure, "Belos::BlockGCRODRSolMgr::solve(): LAPACK GESV failed to compute a solution.");
  // Compute H_m + H_m^{-H}*E_m*H_lbl^{H}*H_lbl
  // H_lbl is bottom-right block of H_, which is a blockSize_ x blockSize_ matrix

  Teuchos::SerialDenseMatrix<int, ScalarType> H_lbl(Teuchos::View, HH, blockSize_, blockSize_, (HH).numRows()-blockSize_, (HH).numCols()-blockSize_ );
  Teuchos::SerialDenseMatrix<int, ScalarType> H_lbl_t( H_lbl, Teuchos::TRANS );

  { //So that HTemp will fall out of scope

    // HH_lbl_t <- H_lbl_t*H_lbl
    Teuchos::RCP<SDM> Htemp = Teuchos::null;
    Htemp = Teuchos::rcp(new SDM(H_lbl_t.numRows(), H_lbl_t.numCols()));
    Htemp -> multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, H_lbl_t, H_lbl, zero);
    H_lbl_t.assign(*Htemp);
    // harmRitzMatrix <- harmRitzMatrix*HH_lbl_t
    Htemp = Teuchos::rcp(new SDM(harmRitzMatrix -> numRows(), harmRitzMatrix -> numCols()));
    Htemp -> multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, *harmRitzMatrix, H_lbl_t, zero);
    harmRitzMatrix -> assign(*Htemp);

    // We need to add harmRitzMatrix to the last blockSize_ columns of HH and store the result harmRitzMatrix
    int harmColIndex, HHColIndex;
    Htemp = Teuchos::rcp(new SDM(Teuchos::Copy,HH,HH.numRows()-blockSize_,HH.numCols()));
    for(int i = 0; i<blockSize_; i++){
      harmColIndex = harmRitzMatrix -> numCols() - i -1;
      HHColIndex = m-i-1;
      for(int j=0; j<m; j++) (*Htemp)[HHColIndex][j] += (*harmRitzMatrix)[harmColIndex][j];
    }
    harmRitzMatrix = Htemp;
  }

  // Revise to do query for optimal workspace first
  // Create simple storage for the left eigenvectors, which we don't care about.

  // Size of workspace and workspace for DGEEV
  int lwork = -1;
  std::vector<ScalarType> work(1);
  std::vector<MagnitudeType> rwork(2*m);

  const int ldvl = 1;
  ScalarType* vl = 0;

  // First query GEEV for the optimal workspace size
  lapack.GEEV('N', 'V', m, harmRitzMatrix->values(), harmRitzMatrix->stride(), &wr[0], &wi[0],
              vl, ldvl, vr.values(), vr.stride(), &work[0], lwork, &rwork[0], &info);

  lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ScalarType>::real (work[0])));
  work.resize( lwork );

  lapack.GEEV('N', 'V', m, harmRitzMatrix->values(), harmRitzMatrix->stride(), &wr[0], &wi[0],
              vl, ldvl, vr.values(), vr.stride(), &work[0], lwork, &rwork[0], &info);

  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure,"Belos::BlockGCRODRSolMgr::solve(): LAPACK GEEV failed to compute eigensolutions.");

  // Construct magnitude of each harmonic Ritz value
  for( int i=0; i<m; ++i ) w[i] = Teuchos::ScalarTraits<MagnitudeType>::squareroot( wr[i]*wr[i] + wi[i]*wi[i] );

  this->sort(w, m, iperm);

  bool scalarTypeIsComplex = Teuchos::ScalarTraits<ScalarType>::isComplex;

  // Select recycledBlocks_ smallest eigenvectors
  for( int i=0; i<recycledBlocks_; ++i )
    for(int j=0; j<m; j++ )
      PP(j,i) = vr(j,iperm[i]);

  if(scalarTypeIsComplex==false) {

    // Determine exact size for PP (i.e., determine if we need to store an additional vector)
    if (wi[iperm[recycledBlocks_-1]] != 0.0) {
      int countImag = 0;
      for (int i=0; i<recycledBlocks_; ++i )
        if (wi[iperm[i]] != 0.0) countImag++;
      // Check to see if this count is even or odd:
      if (countImag % 2) xtraVec = true;
    }

    if (xtraVec) { // we need to store one more vector
      if (wi[iperm[recycledBlocks_-1]] > 0.0) { // I picked the "real" component
        for(int j=0; j<m; ++j )               // so get the "imag" component
          PP(j,recycledBlocks_) = vr(j,iperm[recycledBlocks_-1]+1);
      }
      else{                         // I picked the "imag" component
        for(int j=0; j<m; ++j )     // so get the "real" component
          PP(j,recycledBlocks_) = vr(j,iperm[recycledBlocks_-1]-1);
      }
    }

  }

  // Return whether we needed to store an additional vector
  if (xtraVec) {
    printer_->stream(Debug) << "Recycled " << recycledBlocks_+1 << " vectors" << std::endl;
    return recycledBlocks_+1;
  }
  else {
    printer_->stream(Debug) << "Recycled " << recycledBlocks_ << " vectors" << std::endl;
    return recycledBlocks_;
  }
} //end getHarmonicVecsKryl

// This method sorts list of n floating-point numbers and return permutation vector
template<class ScalarType, class MV, class OP>
void BlockGCRODRSolMgr<ScalarType,MV,OP>::sort(std::vector<MagnitudeType>& dlist, int n, std::vector<int>& iperm) {
  int l, r, j, i, flag;
  int    RR2;
  MagnitudeType dRR, dK;

  // Initialize the permutation vector.
  for(j=0;j<n;j++)
    iperm[j] = j;

  if (n <= 1) return;

  l    = n / 2 + 1;
  r    = n - 1;
  l    = l - 1;
  dRR  = dlist[l - 1];
  dK   = dlist[l - 1];

  RR2 = iperm[l - 1];
  while (r != 0) {
    j = l;
    flag = 1;
    while (flag == 1) {
      i = j;
      j = j + j;
      if (j > r + 1)
        flag = 0;
      else {
        if (j < r + 1)
          if (dlist[j] > dlist[j - 1]) j = j + 1;
        if (dlist[j - 1] > dK) {
          dlist[i - 1] = dlist[j - 1];
          iperm[i - 1] = iperm[j - 1];
        }
        else {
          flag = 0;
        }
      }
    }
    dlist[i - 1] = dRR;
    iperm[i - 1] = RR2;
    if (l == 1) {
      dRR  = dlist [r];
      RR2 = iperm[r];
      dK = dlist[r];
      dlist[r] = dlist[0];
      iperm[r] = iperm[0];
      r = r - 1;
    }
    else {
      l   = l - 1;
      dRR  = dlist[l - 1];
      RR2  = iperm[l - 1];
      dK   = dlist[l - 1];
    }
  }
  dlist[0] = dRR;
  iperm[0] = RR2;
} //end sort() definition

template<class ScalarType, class MV, class OP>
ReturnType BlockGCRODRSolMgr<ScalarType,MV,OP>::solve() {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;

  // MLP: NEED TO ADD CHECK IF PARAMETERS ARE SET LATER

  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
  std::vector<int> index(numBlocks_+1);

  // MLP: EXCEPTION TESTING SHOULD GO HERE

  //THE FOLLOWING BLOCK OF CODE IS TO INFORM THE PROBLEM HOW MANY RIGHT HAND
  //SIDES WILL BE SOLVED.  IF THE CHOSEN BLOCK SIZE IS THE SAME AS THE NUMBER
  //OF RIGHT HAND SIDES IN THE PROBLEM THEN WE PASS THE PROBLEM
  //INDICES FOR ALL RIGHT HAND SIDES.  IF BLOCK SIZE IS GREATER THAN
  //THE NUMBER OF RIGHT HAND SIDES AND THE USER HAS ADAPTIVE BLOCK SIZING
  //TURNED OFF, THEN THE PROBLEM WILL GENERATE RANDOM RIGHT HAND SIDES SO WE HAVE AS MANY
  //RIGHT HAND SIDES AS BLOCK SIZE INDICATES.  THIS MAY NEED TO BE CHANGED
  //LATER TO ACCOMADATE SMARTER CHOOSING OF FICTITIOUS RIGHT HAND SIDES.
  //THIS CODE IS DIRECTLY LIFTED FROM BelosBlockGmresSolMgr.hpp.

  // Create indices for the linear systems to be solved.
  int startPtr = 0;
  int numRHS2Solve = MVT::GetNumberVecs( *(problem_->getRHS()) );
  int numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;

  //currIdx holds indices to all right-hand sides to be solved
  //or -1's to indicate that random right-hand sides should be generated
  std::vector<int> currIdx;

  if ( adaptiveBlockSize_ ) {
    blockSize_ = numCurrRHS;
    currIdx.resize( numCurrRHS  );
    for (int i=0; i<numCurrRHS; ++i)
      currIdx[i] = startPtr+i;
  }
  else {
    currIdx.resize( blockSize_ );
    for (int i=0; i<numCurrRHS; ++i)
      currIdx[i] = startPtr+i;
    for (int i=numCurrRHS; i<blockSize_; ++i)
      currIdx[i] = -1;
  }

  // Inform the linear problem of the linear systems to solve/generate.
  problem_->setLSIndex( currIdx );

  //ADD ERROR CHECKING TO MAKE SURE SIZE OF BLOCK KRYLOV SUBSPACE NOT LARGER THAN dim
  //ptrdiff_t dim = MVT::GetGlobalLength( *(problem_->getRHS()) );

  // reset loss of orthogonality flag
  loaDetected_ = false;

  // Assume convergence is achieved, then let any failed convergence set this to false.
  bool isConverged = true;

  // Initialize storage for all state variables
  initializeStateStorage();

  // Parameter list MLP: WE WILL NEED TO ASSIGN SOME PARAMETER VALUES LATER
  Teuchos::ParameterList plist;

  while (numRHS2Solve > 0){ //This loops through each block of RHS's being solved
    // **************************************************************************************************************************
    // Begin initial solver preparations.  Either update U,C for new operator or generate an initial space using a cycle of GMRES
    // **************************************************************************************************************************
    //int prime_iterations; // silence warning about not being used
    if(keff > 0){//If there is already a subspace to recycle, then use it
      // Update U, C for the new operator

//      TEUCHOS_TEST_FOR_EXCEPTION(keff < recycledBlocks_,BlockGCRODRSolMgrRecyclingFailure,
//                                 "Belos::BlockGCRODRSolMgr::solve(): Requested size of recycled subspace is not consistent with the current recycle subspace.");
                                 printer_->stream(Debug) << " Now solving RHS index " << currIdx[0] << " using recycled subspace of dimension " << keff << std::endl << std::endl;

      // Compute image of U_ under the new operator
      index.resize(keff);
      for (int ii=0; ii<keff; ++ii) { index[ii] = ii; }
      RCP<const MV> Utmp  = MVT::CloneView( *U_, index );
      RCP<MV> Ctmp  = MVT::CloneViewNonConst( *C_, index );
      problem_->apply( *Utmp, *Ctmp );

      RCP<MV> U1tmp = MVT::CloneViewNonConst( *U1_, index );

      // Orthogonalize this block
      // Get a matrix to hold the orthonormalization coefficients.
      SDM Ftmp( Teuchos::View, *F_, keff, keff );
      int rank = ortho_->normalize(*Ctmp, rcp(&Ftmp,false));
      // Throw an error if we could not orthogonalize this block
      TEUCHOS_TEST_FOR_EXCEPTION(rank != keff,BlockGCRODRSolMgrOrthoFailure,"Belos::BlockGCRODRSolMgr::solve(): Failed to compute orthonormal basis for initial recycled subspace.");

      // U_ = U_*F^{-1}
      // First, compute LU factorization of R
      int info = 0;
      ipiv_.resize(Ftmp.numRows());
      lapack.GETRF(Ftmp.numRows(),Ftmp.numCols(),Ftmp.values(),Ftmp.stride(),&ipiv_[0],&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure,"Belos::BlockGCRODRSolMgr::solve(): LAPACK _GETRF failed to compute an LU factorization.");

      // Now, form inv(F)
      int lwork = Ftmp.numRows();
      work_.resize(lwork);
      lapack.GETRI(Ftmp.numRows(),Ftmp.values(),Ftmp.stride(),&ipiv_[0],&work_[0],lwork,&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, BlockGCRODRSolMgrLAPACKFailure,"Belos::BlockGCRODRSolMgr::solve(): LAPACK _GETRI failed to invert triangular matrix.");

      // U_ = U1_; (via a swap)
      MVT::MvTimesMatAddMv( one, *Utmp, Ftmp, zero, *U1tmp );
      std::swap(U_, U1_);

      // Must reinitialize after swap
      index.resize(keff);
      for (int ii=0; ii<keff; ++ii) { index[ii] = ii; }
      Ctmp  = MVT::CloneViewNonConst( *C_, index );
      Utmp  = MVT::CloneView( *U_, index );

      // Compute C_'*R_
      SDM Ctr(keff,blockSize_);
      problem_->computeCurrPrecResVec( &*R_ );
      MVT::MvTransMv( one, *Ctmp, *R_, Ctr );

      // Update solution ( x += U_*C_'*R_ )
      RCP<MV> update = MVT::Clone( *problem_->getCurrLHSVec(), blockSize_ );
      MVT::MvInit( *update, 0.0 );
      MVT::MvTimesMatAddMv( one, *Utmp, Ctr, one, *update );
      problem_->updateSolution( update, true );

      // Update residual norm ( r -= C_*C_'*R_ )
      MVT::MvTimesMatAddMv( -one, *Ctmp, Ctr, one, *R_ );

      // We recycled space from previous call
      //prime_iterations = 0;

      // Since we have a previous recycle space, we do not need block_gmres_iter
      // and we must allocate V ourselves. See comments in initialize() routine for
      // further explanation.
      if (V_ == Teuchos::null) {
        // Check if there is a multivector to clone from.
        Teuchos::RCP<const MV> rhsMV = problem_->getRHS();
        V_ = MVT::Clone( *rhsMV, (numBlocks_+1)*blockSize_ );
      }
      else{
        // Generate V_ by cloning itself ONLY if more space is needed.
        if (MVT::GetNumberVecs(*V_) < (numBlocks_+1)*blockSize_ ) {
          Teuchos::RCP<const MV> tmp = V_;
          V_ = MVT::Clone( *tmp, (numBlocks_+1)*blockSize_ );
        }
      }
    } //end if(keff > 0)
    else{ //if there was no subspace to recycle, then build one
      printer_->stream(Debug) << " No recycled subspace available for RHS index " << std::endl << std::endl;

      Teuchos::ParameterList primeList;
      //we set this as numBlocks-1 because that is the Krylov dimension we want
      //so that the size of the Hessenberg created matches that of what we were expecting.
      primeList.set("Num Blocks",numBlocks_-1);
      primeList.set("Block Size",blockSize_);
      primeList.set("Recycled Blocks",0);
      primeList.set("Keep Hessenberg",true);
      primeList.set("Initialize Hessenberg",true);

      ptrdiff_t dim = MVT::GetGlobalLength( *(problem_->getRHS()) );
      if (blockSize_*static_cast<ptrdiff_t>(numBlocks_) > dim) {//if user has selected a total subspace dimension larger than system dimension
        ptrdiff_t tmpNumBlocks = 0;
        if (blockSize_ == 1)
          tmpNumBlocks = dim / blockSize_;  // Allow for a good breakdown.
        else{
          tmpNumBlocks = ( dim - blockSize_) / blockSize_;  // Allow for restarting.
          printer_->stream(Warnings) << "Belos::BlockGCRODRSolMgr::solve():  Warning! Requested Krylov subspace dimension is larger than operator dimension!"
                                     << std::endl << "The maximum number of blocks allowed for the Krylov subspace will be adjusted to " << tmpNumBlocks << std::endl;
          primeList.set("Num Blocks",Teuchos::as<int>(tmpNumBlocks));
        }
      }
      else{
        primeList.set("Num Blocks",numBlocks_-1);
      }
      //Create Block GMRES iteration object to perform one cycle of GMRES
      Teuchos::RCP<BlockGmresIter<ScalarType,MV,OP> > block_gmres_iter;
      block_gmres_iter = Teuchos::rcp( new BlockGmresIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,primeList) );

      // MLP: ADD LOGIC TO DEAL WITH USER ASKING TO GENERATE A LARGER SPACE THAN dim AS IN HEIDI'S BlockGmresSolMgr CODE (DIDN'T WE ALREADY DO THIS SOMEWHERE?)
      block_gmres_iter->setSize( blockSize_, numBlocks_-1 );

      //Create the first block in the current BLOCK Krylov basis (using residual)
      Teuchos::RCP<MV> V_0;
      if (currIdx[blockSize_-1] == -1) {//if augmented with random RHS's
        V_0 = MVT::Clone( *(problem_->getInitPrecResVec()), blockSize_ );
        problem_->computeCurrPrecResVec( &*V_0 );
      }
      else {
        V_0 = MVT::CloneCopy( *(problem_->getInitPrecResVec()), currIdx );
      }

      // Get a matrix to hold the orthonormalization coefficients.
      Teuchos::RCP<SDM > z_0 = Teuchos::rcp( new SDM( blockSize_, blockSize_ ) );

      // Orthonormalize the new V_0
      int rank = ortho_->normalize( *V_0, z_0 );
      (void) rank; // silence compiler warning for unused variable
      // MLP: ADD EXCEPTION IF INITIAL BLOCK IS RANK DEFFICIENT

      // Set the new state and initialize the iteration.
      GmresIterationState<ScalarType,MV> newstate;
      newstate.V = V_0;
      newstate.z = z_0;
      newstate.curDim = 0;
      block_gmres_iter->initializeGmres(newstate);

      bool primeConverged = false;

      try {
        printer_->stream(Debug) << " Preparing to Iterate!!!!" << std::endl << std::endl;
        block_gmres_iter->iterate();

        // *********************
        // check for convergence
        // *********************
        if ( convTest_->getStatus() == Passed ) {
          printer_->stream(Debug) << "We converged during the prime the pump stage" << std::endl << std::endl;
          primeConverged = !(expConvTest_->getLOADetected());
          if ( expConvTest_->getLOADetected() ) {
            // we don't have convergence
            loaDetected_ = true;
            printer_->stream(Warnings) << "Belos::BlockGCRODRSolMgr::solve(): Warning! Solver has experienced a loss of accuracy!" << std::endl;
          }
        }
        // *******************************************
        // check for maximum iterations of block GMRES
        // *******************************************
        else if( maxIterTest_->getStatus() == Passed ) {
          // we don't have convergence
          primeConverged = false;
        }
        // ****************************************************************
        // We need to recycle and continue, print a message indicating this
        // ****************************************************************
        else{ //debug statements
          printer_->stream(Debug) << " We did not converge on priming cycle of Block GMRES.  Therefore we recycle and restart. " << std::endl << std::endl;
        }
      } //end try
      catch (const GmresIterationOrthoFailure &e) {
        // If the block size is not one, it's not considered a lucky breakdown.
        if (blockSize_ != 1) {
          printer_->stream(Errors) << "Error! Caught std::exception in BlockGmresIter::iterate() at iteration "
                                   << block_gmres_iter->getNumIters() << std::endl
                                   << e.what() << std::endl;
         if (convTest_->getStatus() != Passed)
           primeConverged = false;
        }
        else {
          // If the block size is one, try to recover the most recent least-squares solution
          block_gmres_iter->updateLSQR( block_gmres_iter->getCurSubspaceDim() );
          // Check to see if the most recent least-squares solution yielded convergence.
          sTest_->checkStatus( &*block_gmres_iter );
          if (convTest_->getStatus() != Passed)
            isConverged = false;
        }
      } // end catch (const GmresIterationOrthoFailure &e)
      catch (const StatusTestNaNError& e) {
        // A NaN was detected in the solver.  Set the solution to zero and return unconverged.
        achievedTol_ = MT::one();
        Teuchos::RCP<MV> X = problem_->getLHS();
        MVT::MvInit( *X, SCT::zero() );
        printer_->stream(Warnings) << "Belos::BlockGCRODRSolMgr::solve(): Warning! NaN has been detected!" 
                                     << std::endl;
        return Unconverged;
      }
      catch (const std::exception &e) {
        printer_->stream(Errors) << "Error! Caught std::exception in BlockGmresIter::iterate() at iteration "
                                 << block_gmres_iter->getNumIters() << std::endl
                                 << e.what() << std::endl;
        throw;
      }

      // Record number of iterations in generating initial recycle spacec
      //prime_iterations = block_gmres_iter->getNumIters();//instantiated here because it is not needed outside of else{} scope;  we'll see if this is true or not

      // Update the linear problem.
      RCP<MV> update = block_gmres_iter->getCurrentUpdate();
      problem_->updateSolution( update, true );

      // Update the block residual

      problem_->computeCurrPrecResVec( &*R_ );

      // Get the state.
      newstate = block_gmres_iter->getState();
      int p = newstate.curDim;
      (void) p; // silence compiler warning for unused variable
      H_->assign(*(newstate.H));//IS THIS A GOOD IDEA?  I DID IT SO MEMBERS WOULD HAVE ACCESS TO H.
      // Get new Krylov vectors, store in V_

      // Since the block_gmres_iter returns the state
      // with const RCP's we need to cast newstate.V as
      // a non const RCP.  This is okay because block_gmres_iter
      // is about to fall out of scope, so we are simply
      // rescuing *newstate.V from being destroyed so we can
      // use it for future block recycled GMRES cycles
      V_ = rcp_const_cast<MV>(newstate.V);
      newstate.V.release();
      //COMPUTE NEW RECYCLE SPACE SOMEHOW
      buildRecycleSpaceKryl(keff, block_gmres_iter);
      printer_->stream(Debug) << "Generated recycled subspace using RHS index " << currIdx[0] << " of dimension " << keff << std::endl << std::endl;

      // Return to outer loop if the priming solve
      // converged, set the next linear system.
      if (primeConverged) {
        /* MLP:  POSSIBLY INCORRECT CODE WE ARE REPLACING *********************************
        // Inform the linear problem that we are
        // finished with this block linear system.
        problem_->setCurrLS();

        // Update indices for the linear systems to be solved.
        // KMS: Fix the numRHS2Solve; Copy from Heidi's BlockGmres code
        numRHS2Solve -= 1;
        printer_->stream(Debug) << numRHS2Solve << " Right hand sides left to solve" << std::endl;
        if ( numRHS2Solve > 0 ) {
          currIdx[0]++;
          // Set the next indices
          problem_->setLSIndex( currIdx );
        }
        else {
          currIdx.resize( numRHS2Solve );
        }
        *****************************************************************************/

        // Inform the linear problem that we are finished with this block linear system.
        problem_->setCurrLS();

        // Update indices for the linear systems to be solved.
        startPtr += numCurrRHS;
        numRHS2Solve -= numCurrRHS;
        if ( numRHS2Solve > 0 ) {
          numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;
          if ( adaptiveBlockSize_ ) {
            blockSize_ = numCurrRHS;
            currIdx.resize( numCurrRHS  );
            for (int i=0; i<numCurrRHS; ++i) currIdx[i] = startPtr+i;
          }
          else {
            currIdx.resize( blockSize_ );
            for (int i=0; i<numCurrRHS; ++i) currIdx[i] = startPtr+i;
            for (int i=numCurrRHS; i<blockSize_; ++i) currIdx[i] = -1;
          }
          // Set the next indices.
          problem_->setLSIndex( currIdx );
        }
        else {
          currIdx.resize( numRHS2Solve );
        }
        continue;//Begin solving the next block of RHS's
      }  //end if (primeConverged)
    } //end else [if(keff > 0)]

    // *****************************************
    // Initial subspace update/construction done
    // Start cycles of recycled GMRES
    // *****************************************
    Teuchos::ParameterList blockgcrodrList;
    blockgcrodrList.set("Num Blocks",numBlocks_);
    blockgcrodrList.set("Block Size",blockSize_);
    blockgcrodrList.set("Recycled Blocks",keff);

    Teuchos::RCP<BlockGCRODRIter<ScalarType,MV,OP> > block_gcrodr_iter;
    block_gcrodr_iter = Teuchos::rcp( new BlockGCRODRIter<ScalarType,MV,OP>(problem_,printer_,outputTest_,ortho_,blockgcrodrList) );
    BlockGCRODRIterState<ScalarType,MV> newstate;

    index.resize( blockSize_ );
    for(int ii = 0; ii < blockSize_; ii++) index[ii] = ii;
    Teuchos::RCP<MV> V0 =  MVT::CloneViewNonConst( *V_,  index );

    // MLP: MVT::SetBlock(*R_,index,*V0);
    MVT::Assign(*R_,*V0);

    index.resize(keff);//resize to get appropriate recycle space vectors
    for(int i=0; i < keff; i++){ index[i] = i;};
    B_ = rcp(new SDM(Teuchos::View, *G_, keff, numBlocks_*blockSize_, 0, keff));
    H_ = rcp(new SDM(Teuchos::View, *G_, (numBlocks_-1)*blockSize_ + blockSize_, (numBlocks_-1)*blockSize_, keff ,keff ));

    newstate.V = V_;
    newstate.B= B_;
    newstate.U = MVT::CloneViewNonConst(*U_, index);
    newstate.C = MVT::CloneViewNonConst(*C_, index);
    newstate.H = H_;
    newstate.curDim = blockSize_;
    block_gcrodr_iter -> initialize(newstate);

    int numRestarts = 0;

    while(1){//Each execution of this loop is a cycle of block GCRODR
      try{
        block_gcrodr_iter -> iterate();

        // ***********************
        // Check Convergence First
        // ***********************
        if( convTest_->getStatus() == Passed ) {
          // we have convergence
          break;//from while(1)
        } //end if converged

        // ***********************************
        // Check if maximum iterations reached
        // ***********************************
        else if(maxIterTest_->getStatus() == Passed ){
          // no convergence; hit maxit
          isConverged = false;
          break; // from while(1)
        } //end elseif reached maxit

        // **********************************************
        // Check if subspace full; do we need to restart?
        // **********************************************
        else if (block_gcrodr_iter->getCurSubspaceDim() == block_gcrodr_iter->getMaxSubspaceDim()){

          // Update recycled space even if we have reached max number of restarts

          // Update linear problem
          Teuchos::RCP<MV> update = block_gcrodr_iter->getCurrentUpdate();
          problem_->updateSolution(update, true);
          buildRecycleSpaceAugKryl(block_gcrodr_iter);

          printer_->stream(Debug) << " Generated new recycled subspace using RHS index " << currIdx[0] << " of dimension " << keff << std::endl << std::endl;
          // NOTE: If we have hit the maximum number of restarts, then we will quit
          if(numRestarts >= maxRestarts_) {
            isConverged = false;
            break; //from while(1)
          } //end if max restarts

          numRestarts++;

          printer_ -> stream(Debug) << " Performing restart number " << numRestarts << " of " << maxRestarts_ << std::endl << std::endl;

          // Create the restart vector (first block in the current Krylov basis)
          problem_->computeCurrPrecResVec( &*R_ );
          index.resize( blockSize_ );
          for (int ii=0; ii<blockSize_; ++ii) index[ii] = ii;
          Teuchos::RCP<MV> V0 =  MVT::CloneViewNonConst( *V_,  index );
          MVT::SetBlock(*R_,index,*V0);

          // Set the new state and initialize the solver.
          BlockGCRODRIterState<ScalarType,MV> restartState;
          index.resize( numBlocks_*blockSize_ );
          for (int ii=0; ii<(numBlocks_*blockSize_); ++ii) index[ii] = ii;
          restartState.V  = MVT::CloneViewNonConst( *V_,  index );
          index.resize( keff );
          for (int ii=0; ii<keff; ++ii) index[ii] = ii;
          restartState.U  = MVT::CloneViewNonConst( *U_,  index );
          restartState.C  = MVT::CloneViewNonConst( *C_,  index );
          B_ = rcp(new SDM(Teuchos::View, *G_, keff,                  (numBlocks_-1)*blockSize_,     0, keff));
          H_ = rcp(new SDM(Teuchos::View, *G_, numBlocks_*blockSize_, (numBlocks_-1)*blockSize_, keff ,keff ));
          restartState.B = B_;
          restartState.H = H_;
          restartState.curDim = blockSize_;
          block_gcrodr_iter->initialize(restartState);

        } //end else if need to restart

        // ****************************************************************
        // We returned from iterate(), but none of our status tests passed.
        // Something is wrong, and it is probably our fault.
        // ****************************************************************
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Belos::BlockGCRODRSolMgr::solve(): Invalid return from BlockGCRODRIter::iterate().");
        } //end else (no status test passed)

      }//end try
      catch(const BlockGCRODRIterOrthoFailure &e){ //there was an exception
        // Try to recover the most recent least-squares solution
        block_gcrodr_iter->updateLSQR( block_gcrodr_iter->getCurSubspaceDim() );
        // Check to see if the most recent least-squares solution yielded convergence.
        sTest_->checkStatus( &*block_gcrodr_iter );
        if (convTest_->getStatus() != Passed) isConverged = false;
        break;
      }  // end catch orthogonalization failure
      catch(const std::exception &e){
        printer_->stream(Errors) << "Error! Caught exception in BlockGCRODRIter::iterate() at iteration "
                                 << block_gcrodr_iter->getNumIters() << std::endl
                                 << e.what() << std::endl;
        throw;
      } //end catch standard exception
    } //end while(1)

    // Compute the current solution.
    // Update the linear problem.
    Teuchos::RCP<MV> update = block_gcrodr_iter->getCurrentUpdate();
    problem_->updateSolution( update, true );

    /* MLP:  POSSIBLY INCORRECT CODE WE ARE REPLACING *********************************
    // Inform the linear problem that we are finished with this block linear system.
    problem_->setCurrLS();

    // Update indices for the linear systems to be solved.
    numRHS2Solve -= 1;
    if ( numRHS2Solve > 0 ) {
      currIdx[0]++;
      // Set the next indices.
      problem_->setLSIndex( currIdx );
    }
    else {
      currIdx.resize( numRHS2Solve );
    }
    ******************************************************************************/

    // Inform the linear problem that we are finished with this block linear system.
    problem_->setCurrLS();

    // Update indices for the linear systems to be solved.
    startPtr += numCurrRHS;
    numRHS2Solve -= numCurrRHS;
    if ( numRHS2Solve > 0 ) {
      numCurrRHS = ( numRHS2Solve < blockSize_) ? numRHS2Solve : blockSize_;
      if ( adaptiveBlockSize_ ) {
        blockSize_ = numCurrRHS;
        currIdx.resize( numCurrRHS  );
        for (int i=0; i<numCurrRHS; ++i) currIdx[i] = startPtr+i;
      }
      else {
        currIdx.resize( blockSize_ );
        for (int i=0; i<numCurrRHS; ++i) currIdx[i] = startPtr+i;
        for (int i=numCurrRHS; i<blockSize_; ++i) currIdx[i] = -1;
      }
      // Set the next indices.
      problem_->setLSIndex( currIdx );
    }
    else {
      currIdx.resize( numRHS2Solve );
    }

    // If we didn't build a recycle space this solve but ran at least k iterations, force build of new recycle space
    if (!builtRecycleSpace_) {
      buildRecycleSpaceAugKryl(block_gcrodr_iter);
      printer_->stream(Debug) << " Generated new recycled subspace using RHS index " << currIdx[0] << " of dimension " << keff << std::endl << std::endl;
    }
  }//end while (numRHS2Solve > 0)

  // print final summary
  sTest_->print( printer_->stream(FinalSummary) );

  // print timing information
  #ifdef BELOS_TEUCHOS_TIME_MONITOR
  // Calling summarize() can be expensive, so don't call unless the user wants to print out timing details.
  // summarize() will do all the work even if it's passed a "black hole" output stream.
  if (verbosity_ & TimingDetails) Teuchos::TimeMonitor::summarize( printer_->stream(TimingDetails) );
  #endif
  // get iteration information for this solve
  numIters_ = maxIterTest_->getNumIters();

  // get residual information for this solve
  const std::vector<MagnitudeType>* pTestValues = NULL;
  pTestValues = impConvTest_->getTestValue();
  TEUCHOS_TEST_FOR_EXCEPTION(pTestValues == NULL, std::logic_error,
                             "Belos::BlockGCRODRSolMgr::solve(): The implicit convergence test's "
                             "getTestValue() method returned NULL.  Please report this bug to the "
                             "Belos developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(pTestValues->size() < 1, std::logic_error,
                             "Belos::BlockGCRODRSolMgr::solve(): The implicit convergence test's "
                             "getTestValue() method returned a vector of length zero.  Please report "
                             "this bug to the Belos developers.");
  // FIXME (mfh 12 Dec 2011) Does pTestValues really contain the
  // achieved tolerances for all vectors in the current solve(), or
  // just for the vectors from the last deflation?
  achievedTol_ = *std::max_element (pTestValues->begin(), pTestValues->end());

  if (!isConverged) return Unconverged; // return from BlockGCRODRSolMgr::solve()
    return Converged; // return from BlockGCRODRSolMgr::solve()
} //end solve()

} //End Belos Namespace

#endif /* BELOS_BLOCK_GCRODR_SOLMGR_HPP */
