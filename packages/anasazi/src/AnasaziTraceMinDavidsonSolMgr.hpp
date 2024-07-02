// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziTraceMinDavidsonSolMgr.hpp
 *  \brief The Anasazi::TraceMinDavidsonSolMgr provides a solver manager for the TraceMinDavidson eigensolver with expanding subspaces.
 *
 *  For TraceMin with a constant subspace dimension, please see Anasazi::TraceMinSolMgr.
*/

#ifndef ANASAZI_TRACEMIN_DAVIDSON_SOLMGR_HPP
#define ANASAZI_TRACEMIN_DAVIDSON_SOLMGR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTraceMinBaseSolMgr.hpp"
#include "AnasaziTraceMinDavidson.hpp"

/** \example TraceMinDavidsonGeneralizedEx.cpp
    This is an example of how to use the TraceMinDavidsonSolMgr solver manager to solve a generalized eigenvalue problem, using Tpetra data stuctures.
*/

/** \example TraceMinDavidsonLaplacianEx.cpp
    This is an example of how to use the TraceMinDavidsonSolMgr solver manager to compute the Fiedler vector, using Tpetra data stuctures and an Ifpack2 preconditioner.
*/

/** \example TraceMinDavidsonSpecTransEx.cpp
    This is an example of how to use the TraceMinDavidsonSolMgr solver manager to compute the largest eigenpairs of a matrix via an invisible spectral transformation, using Tpetra data stuctures.
*/

/** \example TraceMinDavidsonUserOpEx.cpp
    This is an example of how to use the TraceMinDavidsonSolMgr solver manager to solve a standard eigenvalue problem where the matrix is not explicitly available, using Tpetra data stuctures.
*/



namespace Anasazi {
namespace Experimental {

/*! \class TraceMinDavidsonSolMgr
 *
 *  \brief The Anasazi::TraceMinDavidsonSolMgr provides a flexible solver manager over the TraceMinDavidson eigensolver.
 *
 * This solver manager implements a hard-locking mechanism, whereby eigenpairs designated to be locked are moved from the eigensolver and placed in
 * auxilliary storage. The eigensolver is then restarted and continues to iterate, orthogonal to the locked eigenvectors.
 *
 * The solver manager provides to the solver a StatusTestCombo object constructed as follows:<br>
 *    &nbsp;&nbsp;&nbsp;<tt>combo = globaltest OR lockingtest OR debugtest</tt><br>
 * where
 *    - \c globaltest terminates computation when global convergence has been detected.<br>
 *      It is encapsulated in a StatusTestWithOrdering object, to ensure that computation is terminated
 *      only after the most significant eigenvalues/eigenvectors have met the convergence criteria.<br>
 *      If not specified via setGlobalStatusTest(), \c globaltest is a StatusTestResNorm object which tests the
 *      2-norms of the direct residuals relative to the Ritz values.
 *    - \c lockingtest halts TraceMinDavidson::iterate() in order to deflate converged eigenpairs for locking.<br>
 *      It will query the underlying TraceMinDavidson eigensolver to determine when eigenvectors should be locked.<br>
 *      If not specified via setLockingStatusTest(), \c lockingtest is a StatusTestResNorm object.
 *    - \c debugtest allows a user to specify additional monitoring of the iteration, encapsulated in a StatusTest object<br>
 *      If not specified via setDebugStatusTest(), \c debugtest is ignored.<br>
 *      In most cases, it should return ::Failed; if it returns ::Passed, solve() will throw an AnasaziError exception.
 *
 * Additionally, the solver manager will terminate and restart solve() when the subspace is full,
 * and it will terminate solve and return after a specified number of restarts.
 *
 * Much of this behavior is controlled via parameters and options passed to the
 * solver manager. For more information, see TraceMinDavidsonSolMgr().

 \ingroup anasazi_solver_framework

 \author Alicia Klinvex
 */

template<class ScalarType, class MV, class OP>
class TraceMinDavidsonSolMgr : public TraceMinBaseSolMgr<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  public:

  //! @name Constructors
  //@{

  /*! \brief Basic constructor for TraceMinDavidsonSolMgr.
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager.
   * Since this class inherits from TraceMinBaseSolMgr, it accepts the same options as TraceMinBaseSolMgr(), with a few additions:
   *   - \c "Use Harmonic Ritz Values" - a \c bool specifying whether to compute the harmonic Ritz values.  Unless you're computing interior eigenpairs, this should be false. Default: false
   *   - \c "Block Size" - an \c int specifying the block size to be used by the underlying solver.
   *                       If the eigenvalues are clustered, you may want to use a larger block size to capture the correct multiplicity.
   *                       Default: 1
   *   - \c "Num Blocks" - an \c int specifying the maximum number of blocks in the subspace. After we compute this many blocks, we will restart.
   *   - \c "Num Restart Blocks" - an \c int specifying how many blocks we keep when restarting.  We will extract the most significant
   *                               blocks from the previous restart and use them to reinitialize TraceMinDavidson. Default: NEV / block size
   *   - \c "Maximum Restarts" - an \c int specifying the maximum number of restarts to be performed. Default: 50
   *
   * A special note about Jacobi-Davidson:\n
   * Since Jacobi-Davidson is a special case of TraceMin-Davidson, it does not have its own solver class.
   * If you want to run Jacobi-Davidson, all you need to do is set the following options:
   *   - \c "When To Shift" = \c "Always" (which is the default)
   *   - \c "How To Choose Shift" = \c "Ritz Values" (not recommended)
   *
   * Choosing the Ritz shifts to be equal to the Ritz values destroys TraceMin-Davidson's guaranteed global convergence,
   * so it is recommended that you consider TraceMin-Davidson (which is identical in every other way)
   * with its default adjusted Ritz shifts instead.
   */
  TraceMinDavidsonSolMgr( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
                          Teuchos::ParameterList &pl );
  //@}

  private:
    int maxRestarts_;

    // Returns true if the subspace is full
    bool needToRestart(const Teuchos::RCP< TraceMinBase<ScalarType,MV,OP> > solver)
    {
      return (solver->getCurSubspaceDim() == solver->getMaxSubspaceDim());
    };

    // Performs a restart by reinitializing TraceMinDavidson with the most significant part of the basis
    bool performRestart(int &numRestarts, Teuchos::RCP< TraceMinBase<ScalarType,MV,OP> > solver);

    // Returns a TraceMinDavidson solver
    Teuchos::RCP< TraceMinBase<ScalarType,MV,OP> > createSolver(
            const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
            const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &outputtest,
            const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
            Teuchos::ParameterList &plist)
    {
      return Teuchos::rcp( new TraceMinDavidson<ScalarType,MV,OP>(this->problem_,sorter,this->printer_,outputtest,ortho,plist) );
    };
};


//---------------------------------------------------------------------------//
// Prevent instantiation on complex scalar type
// FIXME: this really is just a current flaw in the implementation, TraceMin
// *should* work for Hermitian matrices
//---------------------------------------------------------------------------//
template <class MagnitudeType, class MV, class OP>
class TraceMinDavidsonSolMgr<std::complex<MagnitudeType>,MV,OP>
{
  public:

    typedef std::complex<MagnitudeType> ScalarType;
    TraceMinDavidsonSolMgr(
            const RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
            Teuchos::ParameterList &pl )
    {
        // Provide a compile error when attempting to instantiate on complex type
        MagnitudeType::this_class_is_missing_a_specialization();
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// Basic constructor for TraceMinDavidsonSolMgr
template<class ScalarType, class MV, class OP>
TraceMinDavidsonSolMgr<ScalarType,MV,OP>::TraceMinDavidsonSolMgr( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem, Teuchos::ParameterList &pl ) :
  TraceMinBaseSolMgr<ScalarType,MV,OP>(problem,pl)
{
  // TODO: Come back tot these exceptions and make the descriptions better.
  maxRestarts_ = pl.get("Maximum Restarts", 50);
  TEUCHOS_TEST_FOR_EXCEPTION(maxRestarts_ <= 0, std::invalid_argument,
         "Anasazi::TraceMinDavidsonSolMgr::constructor(): \"Maximum Restarts\" must be strictly positive.");

  this->useHarmonic_ = pl.get("Use Harmonic Ritz Values", false);

  TEUCHOS_TEST_FOR_EXCEPTION(this->useHarmonic_ && problem->getM() != Teuchos::null, std::invalid_argument, "Anasazi::TraceMinDavidsonSolMgr::constructor(): Harmonic Ritz values do not currently work with generalized eigenvalue problems.  Please disable \"Use Harmonic Ritz Values\".");

  // block size: default is 1
  this->blockSize_ = pl.get("Block Size", 1);
  TEUCHOS_TEST_FOR_EXCEPTION(this->blockSize_ <= 0, std::invalid_argument,
         "Anasazi::TraceMinDavidsonSolMgr::constructor(): \"Block Size\" must be strictly positive.");

  this->numRestartBlocks_ = (int)std::ceil(this->problem_->getNEV() / this->blockSize_);
  this->numRestartBlocks_ = pl.get("Num Restart Blocks", this->numRestartBlocks_);
  TEUCHOS_TEST_FOR_EXCEPTION(this->numRestartBlocks_ <= 0, std::invalid_argument,
         "Anasazi::TraceMinDavidsonSolMgr::constructor(): \"Num Restart Blocks\" must be strictly positive.");

  this->numBlocks_ = pl.get("Num Blocks", 3*this->numRestartBlocks_);
  TEUCHOS_TEST_FOR_EXCEPTION(this->numBlocks_ <= 1, std::invalid_argument,
         "Anasazi::TraceMinDavidsonSolMgr::constructor(): \"Num Blocks\" must be greater than 1.  If you only wish to use one block, please use TraceMinSolMgr instead of TraceMinDavidsonSolMgr.");

  TEUCHOS_TEST_FOR_EXCEPTION(this->numRestartBlocks_ >= this->numBlocks_, std::invalid_argument,
         "Anasazi::TraceMinDavidsonSolMgr::constructor(): \"Num Blocks\" must be strictly greater than \"Num Restart Blocks\".");

  std::stringstream ss;
  ss << "Anasazi::TraceMinDavidsonSolMgr::constructor(): Potentially impossible orthogonality requests. Reduce basis size (" << static_cast<ptrdiff_t>(this->numBlocks_)*this->blockSize_ << ") or locking size (" << this->maxLocked_ << ") because " << static_cast<ptrdiff_t>(this->numBlocks_) << "*" << this->blockSize_ << " + " << this->maxLocked_ << " > " << MVT::GetGlobalLength(*this->problem_->getInitVec()) << ".";
  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<ptrdiff_t>(this->numBlocks_)*this->blockSize_ + this->maxLocked_ > MVT::GetGlobalLength(*this->problem_->getInitVec()),
         std::invalid_argument, ss.str());

  TEUCHOS_TEST_FOR_EXCEPTION(this->maxLocked_ + this->blockSize_ < this->problem_->getNEV(), std::invalid_argument,
         "Anasazi::TraceMinDavidsonSolMgr: Not enough storage space for requested number of eigenpairs.");
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Performs a restart by reinitializing TraceMinDavidson with the most significant part of the basis
template <class ScalarType, class MV, class OP>
bool TraceMinDavidsonSolMgr<ScalarType,MV,OP>::performRestart(int &numRestarts, Teuchos::RCP< TraceMinBase<ScalarType,MV,OP> > solver)
{
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor restimer(*this->_timerRestarting);
#endif

  if ( numRestarts >= maxRestarts_ ) {
    return false; // break from while(1){tm_solver->iterate()}
  }
  numRestarts++;

  this->printer_->stream(IterationDetails) << " Performing restart number " << numRestarts << " of " << maxRestarts_ << std::endl << std::endl;

  TraceMinBaseState<ScalarType,MV> oldstate = solver->getState();
  TraceMinBaseState<ScalarType,MV> newstate;
  int newdim = this->numRestartBlocks_*this->blockSize_;
  std::vector<int> indToCopy(newdim);
  for(int i=0; i<newdim; i++) indToCopy[i] = i;

  // Copy the relevant parts of the old state to the new one
  // This may involve computing parts of X
  if(this->useHarmonic_)
  {
    newstate.V = MVT::CloneView(*solver->getRitzVectors(),indToCopy);
    newstate.curDim = newdim;

  }
  else
  {
    this->copyPartOfState (oldstate, newstate, indToCopy);
  }

  // send the new state to the solver
  newstate.NEV = oldstate.NEV;
  solver->initialize(newstate);

  return true;
}


}} // end Anasazi namespace

#endif /* ANASAZI_TraceMinDavidson_SOLMGR_HPP */
