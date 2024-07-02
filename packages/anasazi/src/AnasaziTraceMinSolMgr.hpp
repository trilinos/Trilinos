// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_TRACEMIN_SOLMGR_HPP
#define ANASAZI_TRACEMIN_SOLMGR_HPP

/*! \file AnasaziTraceMinSolMgr.hpp
 *  \brief The Anasazi::TraceMinSolMgr provides a solver manager for the TraceMin eigensolver with a constant subspace dimension.
 *
 *  For TraceMin with expanding subspaces, please see Anasazi::TraceMinDavidsonSolMgr.
 */

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverUtils.hpp"

#include "AnasaziTraceMin.hpp"
#include "AnasaziTraceMinBaseSolMgr.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestWithOrdering.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_TimeMonitor.hpp"
#ifdef TEUCHOS_DEBUG
#  include <Teuchos_FancyOStream.hpp>
#endif
#ifdef HAVE_MPI
#include <mpi.h>
#endif


namespace Anasazi {
namespace Experimental {

template<class ScalarType, class MV, class OP>

/*! \class TraceMinSolMgr
 *
 *  \brief The Anasazi::TraceMinSolMgr provides a flexible solver manager over the TraceMin eigensolver.
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
 *    - \c lockingtest halts TraceMin::iterate() in order to deflate converged eigenpairs for locking.<br>
 *      It will query the underlying TraceMin eigensolver to determine when eigenvectors should be locked.<br>
 *      If not specified via setLockingStatusTest(), \c lockingtest is a StatusTestResNorm object.
 *    - \c debugtest allows a user to specify additional monitoring of the iteration, encapsulated in a StatusTest object<br>
 *      If not specified via setDebugStatusTest(), \c debugtest is ignored.<br> 
 *      In most cases, it should return ::Failed; if it returns ::Passed, solve() will throw an AnasaziError exception.
 *
 * Additionally, the solver manager will terminate solve() after a specified number of iterations.
 * 
 * Much of this behavior is controlled via parameters and options passed to the
 * solver manager. For more information, see TraceMinSolMgr().

 \ingroup anasazi_solver_framework

 \author Alicia Klinvex
 */

class TraceMinSolMgr : public TraceMinBaseSolMgr<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:

  //! @name Constructors
  //@{ 

  /*! \brief Basic constructor for TraceMinSolMgr. 
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager. 
   * Since this class inherits from TraceMinBaseSolMgr, it accepts the same options as TraceMinBaseSolMgr(), with a few additions:
   *   - \c "Block Size" - an \c int specifying the block size to be used by the underlying solver. 
   *                       A larger block size means more work per iteration, but it may also decrease the number of
   *                       iterations required. Default: 2*problem->getNEV()
   *   - \c "Maximum Iterations" - an \c int specifying the maximum number of TraceMin iterations to be performed. Default: 100
   */
  TraceMinSolMgr( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );
  //@}

  private:
  
  int maxits_;

  // Test whether we have exceeded the maximum number of iterations
  bool exceededMaxIter() { return (this->iter_ >= maxits_); };

  // TraceMin does not restart, so this will always return false
  bool needToRestart(const Teuchos::RCP< TraceMinBase<ScalarType,MV,OP> > solver) { return false; };

  // TraceMin does not restart, so this will throw an exception
  bool performRestart(int &numRestarts, Teuchos::RCP< TraceMinBase<ScalarType,MV,OP> > solver)
  { TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Anasazi::TraceMinSolMgr::performRestart(): TraceMin does not perform restarts!"); };

  // Returns a new TraceMin solver object
  Teuchos::RCP< TraceMinBase<ScalarType,MV,OP> > createSolver( 
            const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
            const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &outputtest,
            const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
            Teuchos::ParameterList &plist
          );
};


//---------------------------------------------------------------------------//
// Prevent instantiation on complex scalar type
// FIXME: this really is just a current flaw in the implementation, TraceMin
// *should* work for Hermitian matrices
//---------------------------------------------------------------------------//
template <class MagnitudeType, class MV, class OP>
class TraceMinSolMgr<std::complex<MagnitudeType>,MV,OP>
{
  public:

    typedef std::complex<MagnitudeType> ScalarType;
    TraceMinSolMgr(
            const RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
            Teuchos::ParameterList &pl )
    {
        // Provide a compile error when attempting to instantiate on complex type
        MagnitudeType::this_class_is_missing_a_specialization();
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor - accepts maximum iterations in addition to the other parameters of the abstract base class
template<class ScalarType, class MV, class OP>
TraceMinSolMgr<ScalarType,MV,OP>::TraceMinSolMgr( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem, Teuchos::ParameterList &pl ) :
      TraceMinBaseSolMgr<ScalarType,MV,OP>(problem,pl)
{
  // Get the maximum number of iterations
  maxits_ = pl.get("Maximum Iterations", 100);
  TEUCHOS_TEST_FOR_EXCEPTION(maxits_ < 1, std::invalid_argument, "Anasazi::TraceMinSolMgr::constructor(): \"Maximum Iterations\" must be strictly positive.");

  // block size: default is 2* nev()
  // TODO: Find out minimum value
  this->blockSize_ = pl.get("Block Size",2*this->problem_->getNEV());
  TEUCHOS_TEST_FOR_EXCEPTION(this->blockSize_ < this->problem_->getNEV(), std::invalid_argument,
         "Anasazi::TraceMinSolMgr::constructor(): \"Block Size\" must be greater than or equal to the number of desired eigenpairs.");

  this->useHarmonic_ = pl.get("Use Harmonic Ritz Values", false);
  TEUCHOS_TEST_FOR_EXCEPTION(this->useHarmonic_, std::invalid_argument,
         "Anasazi::TraceMinSolMgr::constructor(): Please disable the harmonic Ritz values.  It doesn't make sense to use them with TraceMin, which does not use expanding subspaces.  Perhaps you wanted TraceMin-Davidson?");

  // TraceMin does not restart, so the number of blocks and number of restart blocks will always be 1
  this->numBlocks_ = 1;
  this->numRestartBlocks_ = 1;

  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<ptrdiff_t>(this->numBlocks_)*this->blockSize_ + this->maxLocked_ > MVT::GetGlobalLength(*this->problem_->getInitVec()),
         std::invalid_argument,
         "Anasazi::TraceMinSolMgr::constructor(): Potentially impossible orthogonality requests. Reduce basis size or locking size.");

  TEUCHOS_TEST_FOR_EXCEPTION(this->maxLocked_ + this->blockSize_ < this->problem_->getNEV(), std::invalid_argument,
         "Anasazi::TraceMinDavidsonSolMgr: Not enough storage space for requested number of eigenpairs.");
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Returns a new TraceMin solver object
template <class ScalarType, class MV, class OP>
Teuchos::RCP< TraceMinBase<ScalarType,MV,OP> > TraceMinSolMgr<ScalarType,MV,OP>::createSolver( 
            const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
            const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &outputtest,
            const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
            Teuchos::ParameterList &plist
          )
{
  return Teuchos::rcp( new TraceMin<ScalarType,MV,OP>(this->problem_,sorter,this->printer_,outputtest,ortho,plist) );
}


}} // end Anasazi namespace

#endif /* ANASAZI_TRACEMIN_SOLMGR_HPP */
