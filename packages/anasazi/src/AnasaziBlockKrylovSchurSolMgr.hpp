// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP
#define ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP

/*! \file AnasaziBlockKrylovSchurSolMgr.hpp
 *  \brief The Anasazi::BlockKrylovSchurSolMgr provides a solver manager for the BlockKrylovSchur eigensolver.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziSolverUtils.hpp"

#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestWithOrdering.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_TimeMonitor.hpp"


/** \example BlockKrylovSchur/BlockKrylovSchurEpetraEx.cpp
    This is an example of how to use the Anasazi::BlockKrylovSchurSolMgr solver manager, using Epetra data structures.
*/

/** \example BlockKrylovSchur/BlockKrylovSchurEpetraExGenAmesos.cpp
    This is an example of how to use the Anasazi::BlockKrylovSchurSolMgr solver manager to solve a generalized eigenvalue problem, using Epetra data stuctures and the Amesos solver package.
*/

/** \example BlockKrylovSchur/BlockKrylovSchurEpetraExGenAztecOO.cpp
    This is an example of how to use the Anasazi::BlockKrylovSchurSolMgr solver manager to solve a generalized eigenvalue problem, using Epetra data stuctures and the AztecOO solver package.
*/

/** \example BlockKrylovSchur/BlockKrylovSchurEpetraExGenBelos.cpp
    This is an example of how to use the Anasazi::BlockKrylovSchurSolMgr solver manager to solve a generalized eigenvalue problem, using Epetra data stuctures and the Belos solver package.
*/

/** \example BlockKrylovSchur/BlockKrylovSchurEpetraExSVD.cpp
    This is an example of how to use the Anasazi::BlockKrylovSchurSolMgr solver manager to compute an SVD, using Epetra data structures.
*/

namespace Anasazi {


/*! \class BlockKrylovSchurSolMgr
 *
 *  \brief The Anasazi::BlockKrylovSchurSolMgr provides a flexible solver manager over the BlockKrylovSchur eigensolver.
 *
 * The solver manager provides to the solver a StatusTestCombo object constructed as follows:<br>
 *    &nbsp;&nbsp;&nbsp;<tt>combo = globaltest OR debugtest</tt><br>
 * where
 *    - \c globaltest terminates computation when global convergence has been detected.<br>
 *      It is encapsulated in a StatusTestWithOrdering object, to ensure that computation is terminated
 *      only after the most significant eigenvalues/eigenvectors have met the convergence criteria.<br>
 *      If not specified via setGlobalStatusTest(), this test is a StatusTestResNorm object which tests the
 *      2-norms of the Ritz residuals relative to the Ritz values.
 *    - \c debugtest allows a user to specify additional monitoring of the iteration, encapsulated in a StatusTest object<br>
 *      If not specified via setDebugStatusTest(), \c debugtest is ignored.<br> 
 *      In most cases, it should return ::Failed; if it returns ::Passed, solve() will throw an AnasaziError exception.
 *
 * Additionally, the solver manager will terminate solve() after a specified number of restarts.
 * 
 * Much of this behavior is controlled via parameters and options passed to the
 * solver manager. For more information, see BlockKrylovSchurSolMgr().

 \ingroup anasazi_solver_framework

 \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
 */

template<class ScalarType, class MV, class OP>
class BlockKrylovSchurSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef MultiVecTraitsExt<ScalarType,MV> MVText;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:

  //! @name Constructors/Destructor
  //@{ 

  /*! \brief Basic constructor for BlockKrylovSchurSolMgr.
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - Solver parameters
   *      - \c "Which" - a \c string specifying the desired eigenvalues: SM, LM, SR or LR. Default: "LM"
   *      - \c "Block Size" - a \c int specifying the block size to be used by the underlying block Krylov-Schur solver. Default: 1
   *      - \c "Num Blocks" - a \c int specifying the number of blocks allocated for the Krylov basis. Default: 3*nev
   *      - \c "Extra NEV Blocks" - a \c int specifying the number of extra blocks the solver should keep in addition to those
             required to compute the number of eigenvalues requested.  Default: 0
   *      - \c "Maximum Restarts" - a \c int specifying the maximum number of restarts the underlying solver is allowed to perform. Default: 20
   *      - \c "Orthogonalization" - a \c string specifying the desired orthogonalization:  DGKS and SVQB. Default: "SVQB"
   *      - \c "Verbosity" - a sum of MsgType specifying the verbosity. Default: Anasazi::Errors
   *   - Convergence parameters
   *      - \c "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision.
   *      - \c "Relative Convergence Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding convergence. Default: true
   */
  BlockKrylovSchurSolMgr( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~BlockKrylovSchurSolMgr() {};
  //@}
  
  //! @name Accessor methods
  //@{ 

  //! Return the eigenvalue problem.
  const Eigenproblem<ScalarType,MV,OP>& getProblem() const {
    return *_problem;
  }

  //! Get the iteration count for the most recent call to \c solve().
  int getNumIters() const {
    return _numIters;
  }

  /*! \brief Return the Ritz values from the most recent solve.
   */
  std::vector<Value<ScalarType> > getRitzValues() const {
    std::vector<Value<ScalarType> > ret( _ritzValues );
    return ret;
  }

  /*! \brief Return the timers for this object. 
   *
   * The timers are ordered as follows:
   *   - time spent in solve() routine
   *   - time spent restarting
   */
   Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const {
     return Teuchos::tuple(_timerSolve, _timerRestarting);
   }

  //@}

  //! @name Solver application methods
  //@{ 
    
  /*! \brief This method performs possibly repeated calls to the underlying eigensolver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
   * quit.
   *
   * This method calls BlockKrylovSchur::iterate(), which will return either because a specially constructed status test evaluates to ::Passed
   * or an exception is thrown.
   *
   * A return from BlockKrylovSchur::iterate() signifies one of the following scenarios:
   *    - the maximum number of restarts has been exceeded. In this scenario, the solver manager will place\n
   *      all converged eigenpairs into the eigenproblem and return ::Unconverged.
   *    - global convergence has been met. In this case, the most significant NEV eigenpairs in the solver and locked storage  \n
   *      have met the convergence criterion. (Here, NEV refers to the number of eigenpairs requested by the Eigenproblem.)    \n
   *      In this scenario, the solver manager will return ::Converged.
   *
   * \returns ::ReturnType specifying:
   *     - ::Converged: the eigenproblem was solved to the specification required by the solver manager.
   *     - ::Unconverged: the eigenproblem was not solved to the specification desired by the solver manager.
  */
  ReturnType solve();

  //! Set the status test defining global convergence.
  void setGlobalStatusTest(const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &global);

  //! Get the status test defining global convergence.
  const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > & getGlobalStatusTest() const;

  //! Set the status test for debugging.
  void setDebugStatusTest(const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &debug);

  //! Get the status test for debugging.
  const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > & getDebugStatusTest() const;

  //@}

  private:
  Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > _problem;
  Teuchos::RCP<SortManager<MagnitudeType> > _sort;

  std::string _whch, _ortho; 
  MagnitudeType _ortho_kappa;

  MagnitudeType _convtol;
  int _maxRestarts;
  bool _relconvtol,_conjSplit;
  int _blockSize, _numBlocks, _stepSize, _nevBlocks, _xtra_nevBlocks;
  int _numIters;
  int _verbosity;
  bool _inSituRestart;

  std::vector<Value<ScalarType> > _ritzValues;

  Teuchos::RCP<Teuchos::Time> _timerSolve, _timerRestarting;

  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > globalTest_;
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > debugTest_;

  int _printNum;
};


// Constructor
template<class ScalarType, class MV, class OP>
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::BlockKrylovSchurSolMgr( 
        const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl ) : 
  _problem(problem),
  _whch("LM"),
  _ortho("SVQB"),
  _ortho_kappa(-1.0),
  _convtol(0),
  _maxRestarts(20),
  _relconvtol(true),
  _conjSplit(false),
  _blockSize(0),
  _numBlocks(0),
  _stepSize(0),
  _nevBlocks(0),
  _xtra_nevBlocks(0),
  _numIters(0),
  _verbosity(Anasazi::Errors),
  _inSituRestart(false),
  _timerSolve(Teuchos::TimeMonitor::getNewTimer("Anasazi: BlockKrylovSchurSolMgr::solve()")),
  _timerRestarting(Teuchos::TimeMonitor::getNewTimer("Anasazi: BlockKrylovSchurSolMgr restarting")),
  _printNum(-1)
{
  TEUCHOS_TEST_FOR_EXCEPTION(_problem == Teuchos::null,               std::invalid_argument, "Problem not given to solver manager.");
  TEUCHOS_TEST_FOR_EXCEPTION(!_problem->isProblemSet(),               std::invalid_argument, "Problem not set.");
  TEUCHOS_TEST_FOR_EXCEPTION(_problem->getInitVec() == Teuchos::null, std::invalid_argument, "Problem does not contain initial vectors to clone from.");

  const int nev = _problem->getNEV();

  // convergence tolerance
  _convtol = pl.get("Convergence Tolerance",MT::prec());
  _relconvtol = pl.get("Relative Convergence Tolerance",_relconvtol);
  
  // maximum number of restarts
  _maxRestarts = pl.get("Maximum Restarts",_maxRestarts);

  // block size: default is 1
  _blockSize = pl.get("Block Size",1);
  TEUCHOS_TEST_FOR_EXCEPTION(_blockSize <= 0, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Block Size\" must be strictly positive.");

  // set the number of blocks we need to save to compute the nev eigenvalues of interest.
  _xtra_nevBlocks = pl.get("Extra NEV Blocks",0);
  if (nev%_blockSize) {
    _nevBlocks = nev/_blockSize + _xtra_nevBlocks + 1;
  } else {
    _nevBlocks = nev/_blockSize + _xtra_nevBlocks;
  }

  _numBlocks = pl.get("Num Blocks",3*_nevBlocks);
  TEUCHOS_TEST_FOR_EXCEPTION(_numBlocks <= _nevBlocks, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Num Blocks\" must be strictly positive and large enough to compute the requested eigenvalues.");

  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<ptrdiff_t>(_numBlocks)*_blockSize > MVText::GetGlobalLength(*_problem->getInitVec()),
                     std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: Potentially impossible orthogonality requests. Reduce basis size.");
  
  // step size: the default is _maxRestarts*_numBlocks, so that Ritz values are only computed every restart.
  if (_maxRestarts) {
    _stepSize = pl.get("Step Size", (_maxRestarts+1)*(_numBlocks+1));
  } else {
    _stepSize = pl.get("Step Size", _numBlocks+1);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(_stepSize < 1, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Step Size\" must be strictly positive.");

  // get the sort manager
  if (pl.isParameter("Sort Manager")) {
    _sort = Teuchos::getParameter<Teuchos::RCP<Anasazi::SortManager<MagnitudeType> > >(pl,"Sort Manager");
  } else {
    // which values to solve for
    _whch = pl.get("Which",_whch);
    TEUCHOS_TEST_FOR_EXCEPTION(_whch != "SM" && _whch != "LM" && _whch != "SR" && _whch != "LR" && _whch != "SI" && _whch != "LI",
                       std::invalid_argument, "Invalid sorting string.");
    _sort = Teuchos::rcp( new BasicSort<MagnitudeType>(_whch) );
  }

  // which orthogonalization to use
  _ortho = pl.get("Orthogonalization",_ortho);
  if (_ortho != "DGKS" && _ortho != "SVQB") {
    _ortho = "SVQB";
  }

  // which orthogonalization constant to use
  _ortho_kappa = pl.get("Orthogonalization Constant",_ortho_kappa);

  // verbosity level
  if (pl.isParameter("Verbosity")) {
    if (Teuchos::isParameterType<int>(pl,"Verbosity")) {
      _verbosity = pl.get("Verbosity", _verbosity);
    } else {
      _verbosity = (int)Teuchos::getParameter<Anasazi::MsgType>(pl,"Verbosity");
    }
  }

  // restarting technique: V*Q or applyHouse(V,H,tau)
  if (pl.isParameter("In Situ Restarting")) {
    if (Teuchos::isParameterType<bool>(pl,"In Situ Restarting")) {
      _inSituRestart = pl.get("In Situ Restarting",_inSituRestart);
    } else {
      _inSituRestart = ( Teuchos::getParameter<int>(pl,"In Situ Restarting") != 0 );
    }
  }

  _printNum = pl.get<int>("Print Number of Ritz Values",-1);
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType 
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::solve() {

  const int nev = _problem->getNEV();
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  Teuchos::BLAS<int,ScalarType> blas;
  Teuchos::LAPACK<int,ScalarType> lapack;
  typedef SolverUtils<ScalarType,MV,OP> msutils;

  //////////////////////////////////////////////////////////////////////////////////////
  // Output manager
  Teuchos::RCP<BasicOutputManager<ScalarType> > printer = Teuchos::rcp( new BasicOutputManager<ScalarType>(_verbosity) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Status tests
  //
  // convergence
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > convtest;
  if (globalTest_ == Teuchos::null) {
    convtest = Teuchos::rcp( new StatusTestResNorm<ScalarType,MV,OP>(_convtol,nev,StatusTestResNorm<ScalarType,MV,OP>::RITZRES_2NORM,_relconvtol) );
  }
  else {
    convtest = globalTest_;
  }
  Teuchos::RCP<StatusTestWithOrdering<ScalarType,MV,OP> > ordertest 
    = Teuchos::rcp( new StatusTestWithOrdering<ScalarType,MV,OP>(convtest,_sort,nev) );
  // for a non-short-circuited OR test, the order doesn't matter
  Teuchos::Array<Teuchos::RCP<StatusTest<ScalarType,MV,OP> > > alltests;
  alltests.push_back(ordertest);

  if (debugTest_ != Teuchos::null) alltests.push_back(debugTest_);

  Teuchos::RCP<StatusTestCombo<ScalarType,MV,OP> > combotest
    = Teuchos::rcp( new StatusTestCombo<ScalarType,MV,OP>( StatusTestCombo<ScalarType,MV,OP>::OR, alltests) );
  // printing StatusTest
  Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputtest;
  if ( printer->isVerbosity(Debug) ) {
    outputtest = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer,combotest,1,Passed+Failed+Undefined ) );
  }
  else {
    outputtest = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer,combotest,1,Passed ) );
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Orthomanager
  Teuchos::RCP<OrthoManager<ScalarType,MV> > ortho; 
  if (_ortho=="SVQB") {
    ortho = Teuchos::rcp( new SVQBOrthoManager<ScalarType,MV,OP>(_problem->getM()) );
  } else if (_ortho=="DGKS") {
    if (_ortho_kappa <= 0) {
      ortho = Teuchos::rcp( new BasicOrthoManager<ScalarType,MV,OP>(_problem->getM()) );
    }
    else {
      ortho = Teuchos::rcp( new BasicOrthoManager<ScalarType,MV,OP>(_problem->getM(),_ortho_kappa) );
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(_ortho!="SVQB"&&_ortho!="DGKS",std::logic_error,"Anasazi::BlockKrylovSchurSolMgr::solve(): Invalid orthogonalization type.");
  }
  
  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",_blockSize);
  plist.set("Num Blocks",_numBlocks);
  plist.set("Step Size",_stepSize);
  if (_printNum == -1) {
    plist.set("Print Number of Ritz Values",_nevBlocks*_blockSize);
  }
  else {
    plist.set("Print Number of Ritz Values",_printNum);
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockKrylovSchur solver
  Teuchos::RCP<BlockKrylovSchur<ScalarType,MV,OP> > bks_solver 
    = Teuchos::rcp( new BlockKrylovSchur<ScalarType,MV,OP>(_problem,_sort,printer,outputtest,ortho,plist) );
  // set any auxiliary vectors defined in the problem
  Teuchos::RCP< const MV > probauxvecs = _problem->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    bks_solver->setAuxVecs( Teuchos::tuple< Teuchos::RCP<const MV> >(probauxvecs) );
  }

  // Create workspace for the Krylov basis generated during a restart
  // Need at most (_nevBlocks*_blockSize+1) for the updated factorization and another block for the current factorization residual block (F).
  //  ---> (_nevBlocks*_blockSize+1) + _blockSize
  // If Hermitian, this becomes _nevBlocks*_blockSize + _blockSize
  // we only need this if there is the possibility of restarting, ex situ
  Teuchos::RCP<MV> workMV;
  if (_maxRestarts > 0) {
    if (_inSituRestart==true) {
      // still need one work vector for applyHouse()
      workMV = MVT::Clone( *_problem->getInitVec(), 1 );
    }
    else { // inSituRestart == false
      if (_problem->isHermitian()) {
        workMV = MVT::Clone( *_problem->getInitVec(), _nevBlocks*_blockSize + _blockSize );
      } else {
        workMV = MVT::Clone( *_problem->getInitVec(), _nevBlocks*_blockSize+1 + _blockSize );
      }
    }
  } else {
    workMV = Teuchos::null;
  }

  // go ahead and initialize the solution to nothing in case we throw an exception
  Eigensolution<ScalarType,MV> sol;
  sol.numVecs = 0;
  _problem->setSolution(sol);

  int numRestarts = 0;
  int cur_nevBlocks = 0;

  // enter solve() iterations
  {
    Teuchos::TimeMonitor slvtimer(*_timerSolve);
  
    // tell bks_solver to iterate
    while (1) {
      try {
        bks_solver->iterate();
    
        ////////////////////////////////////////////////////////////////////////////////////
        //
        // check convergence first
        //
        ////////////////////////////////////////////////////////////////////////////////////
        if (ordertest->getStatus() == Passed ) {
          // we have convergence
          // ordertest->whichVecs() tells us which vectors from solver state are the ones we want
          // ordertest->howMany() will tell us how many
          break;
        }
        ////////////////////////////////////////////////////////////////////////////////////
        //
        // check for restarting, i.e. the subspace is full
        //
        ////////////////////////////////////////////////////////////////////////////////////
        // this is for the Hermitian case, or non-Hermitian conjugate split situation.
        // --> for the Hermitian case the current subspace dimension needs to match the maximum subspace dimension
        // --> for the non-Hermitian case:
        //     --> if a conjugate pair was detected in the previous restart then the current subspace dimension needs to match the
        //         maximum subspace dimension (the BKS solver keeps one extra vector if the problem is non-Hermitian).
        //     --> if a conjugate pair was not detected in the previous restart then the current subspace dimension will be one less
        //         than the maximum subspace dimension.
        else if ( (bks_solver->getCurSubspaceDim() == bks_solver->getMaxSubspaceDim()) ||
                  (!_problem->isHermitian() && !_conjSplit && (bks_solver->getCurSubspaceDim()+1 == bks_solver->getMaxSubspaceDim())) ) {
  
          // Update the Schur form of the projected eigenproblem, then sort it.
          if (!bks_solver->isSchurCurrent()) {
            bks_solver->computeSchurForm( true );

            // Check for convergence, just in case we wait for every restart to check
            outputtest->checkStatus( &*bks_solver );  
          }

          // Don't bother to restart if we've converged or reached the maximum number of restarts
          if ( numRestarts >= _maxRestarts || ordertest->getStatus() == Passed) {
            break; // break from while(1){bks_solver->iterate()}
          }

          // Start restarting timer and increment counter 
          Teuchos::TimeMonitor restimer(*_timerRestarting);
          numRestarts++;
  
          printer->stream(Debug) << " Performing restart number " << numRestarts << " of " << _maxRestarts << std::endl << std::endl;
  
          // Get the most current Ritz values before we continue.
          _ritzValues = bks_solver->getRitzValues();

          // Get the state.
          BlockKrylovSchurState<ScalarType,MV> oldState = bks_solver->getState();

          // Get the current dimension of the factorization
          int curDim = oldState.curDim;

          // Determine if the storage for the nev eigenvalues of interest splits a complex conjugate pair.
          std::vector<int> ritzIndex = bks_solver->getRitzIndex();
          if (ritzIndex[_nevBlocks*_blockSize-1]==1) {
            _conjSplit = true;
            cur_nevBlocks = _nevBlocks*_blockSize+1;
          } else {
            _conjSplit = false;
            cur_nevBlocks = _nevBlocks*_blockSize;
          }

          // Update the Krylov-Schur decomposition

          // Get a view of the Schur vectors of interest.
          Teuchos::SerialDenseMatrix<int,ScalarType> Qnev(Teuchos::View, *(oldState.Q), curDim, cur_nevBlocks);

          // Get a view of the current Krylov basis.
          std::vector<int> curind( curDim );
          for (int i=0; i<curDim; i++) { curind[i] = i; }
          Teuchos::RCP<const MV> basistemp = MVT::CloneView( *(oldState.V), curind );

          // Compute the new Krylov basis: Vnew = V*Qnev
          // 
          // this will occur ex situ in workspace allocated for this purpose (tmpMV)
          // or in situ in the solver's memory space.
          //
          // we will also set a pointer for the location that the current factorization residual block (F),
          // currently located after the current basis in oldstate.V, will be moved to
          //
          Teuchos::RCP<MV> newF;
          if (_inSituRestart) {
            //
            // get non-const pointer to solver's basis so we can work in situ
            Teuchos::RCP<MV> solverbasis = Teuchos::rcp_const_cast<MV>(oldState.V);
            Teuchos::SerialDenseMatrix<int,ScalarType> copyQnev(Qnev);
            // 
            // perform Householder QR of copyQnev = Q [D;0], where D is unit diag. We will want D below.
            std::vector<ScalarType> tau(cur_nevBlocks), work(cur_nevBlocks);
            int info;
            lapack.GEQRF(curDim,cur_nevBlocks,copyQnev.values(),copyQnev.stride(),&tau[0],&work[0],work.size(),&info);
            TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
                               "Anasazi::BlockDavidsonSolMgr::solve(): error calling GEQRF during restarting.");
            // we need to get the diagonal of D
            std::vector<ScalarType> d(cur_nevBlocks);
            for (int j=0; j<copyQnev.numCols(); j++) {
              d[j] = copyQnev(j,j);
            }
            if (printer->isVerbosity(Debug)) {
              Teuchos::SerialDenseMatrix<int,ScalarType> R(Teuchos::Copy,copyQnev,cur_nevBlocks,cur_nevBlocks);
              for (int j=0; j<R.numCols(); j++) {
                R(j,j) = SCT::magnitude(R(j,j)) - 1.0;
                for (int i=j+1; i<R.numRows(); i++) {
                  R(i,j) = zero;
                }
              }
              printer->stream(Debug) << "||Triangular factor of Su - I||: " << R.normFrobenius() << std::endl;
            }
            // 
            // perform implicit V*Qnev
            // this actually performs V*[Qnev Qtrunc*M] = [newV truncV], for some unitary M
            // we are interested in only the first cur_nevBlocks vectors of the result
            curind.resize(curDim);
            for (int i=0; i<curDim; i++) curind[i] = i;
            {
              Teuchos::RCP<MV> oldV = MVT::CloneViewNonConst(*solverbasis,curind);
              msutils::applyHouse(cur_nevBlocks,*oldV,copyQnev,tau,workMV);
            }
            // multiply newV*D
            // get pointer to new basis
            curind.resize(cur_nevBlocks);
            for (int i=0; i<cur_nevBlocks; i++) { curind[i] = i; }
            {
              Teuchos::RCP<MV> newV = MVT::CloneViewNonConst( *solverbasis, curind );
              MVT::MvScale(*newV,d);
            }
            // get pointer to new location for F
            curind.resize(_blockSize);
            for (int i=0; i<_blockSize; i++) { curind[i] = cur_nevBlocks + i; }
            newF = MVT::CloneViewNonConst( *solverbasis, curind );
          }
          else {
            // get pointer to first part of work space
            curind.resize(cur_nevBlocks);
            for (int i=0; i<cur_nevBlocks; i++) { curind[i] = i; }
            Teuchos::RCP<MV> tmp_newV = MVT::CloneViewNonConst(*workMV, curind );
            // perform V*Qnev
            MVT::MvTimesMatAddMv( one, *basistemp, Qnev, zero, *tmp_newV );
            tmp_newV = Teuchos::null;
            // get pointer to new location for F
            curind.resize(_blockSize);
            for (int i=0; i<_blockSize; i++) { curind[i] = cur_nevBlocks + i; }
            newF = MVT::CloneViewNonConst( *workMV, curind );
          }

          // Move the current factorization residual block (F) to the last block of newV.
          curind.resize(_blockSize);
          for (int i=0; i<_blockSize; i++) { curind[i] = curDim + i; }
          Teuchos::RCP<const MV> oldF = MVT::CloneView( *(oldState.V), curind );
          for (int i=0; i<_blockSize; i++) { curind[i] = i; }
          MVT::SetBlock( *oldF, curind, *newF );
          newF = Teuchos::null;

          // Update the Krylov-Schur quasi-triangular matrix.
          //
          // Create storage for the new Schur matrix of the Krylov-Schur factorization
          // Copy over the current quasi-triangular factorization of oldState.H which is stored in oldState.S.
          Teuchos::SerialDenseMatrix<int,ScalarType> oldS(Teuchos::View, *(oldState.S), cur_nevBlocks+_blockSize, cur_nevBlocks);
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > newH = 
            Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( oldS ) );
          //
          // Get a view of the B block of the current factorization
          Teuchos::SerialDenseMatrix<int,ScalarType> oldB(Teuchos::View, *(oldState.H), _blockSize, _blockSize, curDim, curDim-_blockSize);
          //
          // Get a view of the a block row of the Schur vectors.
          Teuchos::SerialDenseMatrix<int,ScalarType> subQ(Teuchos::View, *(oldState.Q), _blockSize, cur_nevBlocks, curDim-_blockSize);
          //
          // Get a view of the new B block of the updated Krylov-Schur factorization
          Teuchos::SerialDenseMatrix<int,ScalarType> newB(Teuchos::View, *newH,  _blockSize, cur_nevBlocks, cur_nevBlocks);
          //
          // Compute the new B block.
          blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, _blockSize, cur_nevBlocks, _blockSize, one, 
                     oldB.values(), oldB.stride(), subQ.values(), subQ.stride(), zero, newB.values(), newB.stride() );


          //
          // Set the new state and initialize the solver.
          BlockKrylovSchurState<ScalarType,MV> newstate;
          if (_inSituRestart) {
            newstate.V = oldState.V;
          } else {
            newstate.V = workMV;
          }
          newstate.H = newH;
          newstate.curDim = cur_nevBlocks;
          bks_solver->initialize(newstate);
  
        } // end of restarting
        ////////////////////////////////////////////////////////////////////////////////////
        //
        // we returned from iterate(), but none of our status tests Passed.
        // something is wrong, and it is probably our fault.
        //
        ////////////////////////////////////////////////////////////////////////////////////
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::BlockKrylovSchurSolMgr::solve(): Invalid return from bks_solver::iterate().");
        }
      }
      catch (const AnasaziError &err) {
        printer->stream(Errors) 
          << "Anasazi::BlockKrylovSchurSolMgr::solve() caught unexpected exception from Anasazi::BlockKrylovSchur::iterate() at iteration " << bks_solver->getNumIters() << std::endl
          << err.what() << std::endl
          << "Anasazi::BlockKrylovSchurSolMgr::solve() returning Unconverged with no solutions." << std::endl;
        return Unconverged;
      }
    }

    //
    // free temporary space
    workMV = Teuchos::null;

    // Get the most current Ritz values before we return
    _ritzValues = bks_solver->getRitzValues();

    sol.numVecs = ordertest->howMany();
    printer->stream(Debug) << "ordertest->howMany() : " << sol.numVecs << std::endl;
    std::vector<int> whichVecs = ordertest->whichVecs();

    // Place any converged eigenpairs in the solution container.
    if (sol.numVecs > 0) {

      // Next determine if there is a conjugate pair on the boundary and resize.
      std::vector<int> tmpIndex = bks_solver->getRitzIndex();
      for (int i=0; i<(int)_ritzValues.size(); ++i) {
        printer->stream(Debug) << _ritzValues[i].realpart << " + i " << _ritzValues[i].imagpart << ", Index = " << tmpIndex[i] << std::endl;
      }
      printer->stream(Debug) << "Number of converged eigenpairs (before) = " << sol.numVecs << std::endl;
      for (int i=0; i<sol.numVecs; ++i) {
        printer->stream(Debug) << "whichVecs[" << i << "] = " << whichVecs[i] << ", tmpIndex[" << whichVecs[i] << "] = " << tmpIndex[whichVecs[i]] << std::endl;
      }
      if (tmpIndex[whichVecs[sol.numVecs-1]]==1) {
        printer->stream(Debug) << "There is a conjugate pair on the boundary, resizing sol.numVecs" << std::endl;
        whichVecs.push_back(whichVecs[sol.numVecs-1]+1);
	sol.numVecs++;
        for (int i=0; i<sol.numVecs; ++i) {
          printer->stream(Debug) << "whichVecs[" << i << "] = " << whichVecs[i] << ", tmpIndex[" << whichVecs[i] << "] = " << tmpIndex[whichVecs[i]] << std::endl;
        }
      }

      bool keepMore = false;
      int numEvecs = sol.numVecs;
      printer->stream(Debug) << "Number of converged eigenpairs (after) = " << sol.numVecs << std::endl;
      printer->stream(Debug) << "whichVecs[sol.numVecs-1] > sol.numVecs-1 : " << whichVecs[sol.numVecs-1] << " > " << sol.numVecs-1 << std::endl;
      if (whichVecs[sol.numVecs-1] > (sol.numVecs-1)) {
	keepMore = true;
	numEvecs = whichVecs[sol.numVecs-1]+1;  // Add 1 to fix zero-based indexing
        printer->stream(Debug) << "keepMore = true; numEvecs = " << numEvecs << std::endl;
      }

      // Next set the number of Ritz vectors that the iteration must compute and compute them.
      bks_solver->setNumRitzVectors(numEvecs);
      bks_solver->computeRitzVectors();

      // If the leading Ritz pairs are the converged ones, get the information 
      // from the iteration to the solution container. Otherwise copy the necessary
      // information using 'whichVecs'.
      if (!keepMore) {
	sol.index = bks_solver->getRitzIndex();
	sol.Evals = bks_solver->getRitzValues();
	sol.Evecs = MVT::CloneCopy( *(bks_solver->getRitzVectors()) );
      }

      // Resize based on the number of solutions being returned and set the number of Ritz
      // vectors for the iteration to compute.
      sol.Evals.resize(sol.numVecs);
      sol.index.resize(sol.numVecs);
 
      // If the converged Ritz pairs are not the leading ones, copy over the information directly.      
      if (keepMore) {
	std::vector<Anasazi::Value<ScalarType> > tmpEvals = bks_solver->getRitzValues();
	for (int vec_i=0; vec_i<sol.numVecs; ++vec_i) {
	  sol.index[vec_i] = tmpIndex[whichVecs[vec_i]];
	  sol.Evals[vec_i] = tmpEvals[whichVecs[vec_i]];
	}
	sol.Evecs = MVT::CloneCopy( *(bks_solver->getRitzVectors()), whichVecs );
      }

      // Set the solution space to be the Ritz vectors at this time.
      sol.Espace = sol.Evecs;
    } 
  }

  // print final summary
  bks_solver->currentStatus(printer->stream(FinalSummary));

  // print timing information
  Teuchos::TimeMonitor::summarize(printer->stream(TimingDetails));

  _problem->setSolution(sol);
  printer->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << std::endl;

  // get the number of iterations performed during this solve.
  _numIters = bks_solver->getNumIters();

  if (sol.numVecs < nev) {
    return Unconverged; // return from BlockKrylovSchurSolMgr::solve() 
  }
  return Converged; // return from BlockKrylovSchurSolMgr::solve() 
}


template <class ScalarType, class MV, class OP>
void 
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::setGlobalStatusTest(
    const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &global) 
{
  globalTest_ = global;
}

template <class ScalarType, class MV, class OP>
const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > & 
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::getGlobalStatusTest() const 
{
  return globalTest_;
}

template <class ScalarType, class MV, class OP>
void 
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::setDebugStatusTest(
    const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > &debug)
{
  debugTest_ = debug;
}

template <class ScalarType, class MV, class OP>
const Teuchos::RCP< StatusTest<ScalarType,MV,OP> > & 
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::getDebugStatusTest() const
{
  return debugTest_;
}

} // end Anasazi namespace

#endif /* ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP */
