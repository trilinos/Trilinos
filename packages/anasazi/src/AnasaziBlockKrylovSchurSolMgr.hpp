
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
 *  \brief The Anasazi::BlockKrylovSchurSolMgr provides a powerful solver manager for the BlockKrylovSchur eigensolver.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziModalSolverUtils.hpp"

#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestOrderedResNorm.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "Teuchos_BLAS.hpp"

/*! \class Anasazi::BlockKrylovSchurSolMgr
 *
 *  \brief The Anasazi::BlockKrylovSchurSolMgr provides a powerful and fully-featured solver manager over the BlockKrylovSchur eigensolver.
 *
 * This solver manager implements a hard-locking mechanism, whereby eigenpairs designated to be locked are moved from the eigensolver and placed in
 * auxiliary storage. The eigensolver is then restarted and continues to iterate, always orthogonal to the locked eigenvectors.

 \ingroup anasazi_solvermanagers

 \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
 */

namespace Anasazi {

template<class ScalarType, class MV, class OP>
class BlockKrylovSchurSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
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
   *   - "Which" - a \c string specifying the desired eigenvalues: SM, LM, SR or LR. Default: "LM"
   *   - "Block Size" - a \c int specifying the block size to be used by the underlying LOBPCG solver. Default: problem->getNEV()
   *   - "Num Blocks" - a \c int specifying the number of blocks allocated for the Krylov basis. Default: 2
   *   - "Maximum Restarts" - a \c int specifying the maximum number of restarts the underlying solver is allowed to perform. Default: 20
   *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Anasazi::Errors
   *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision.
   *   - "Relative Convergence Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding convergence. Default: true
   *   - "Use Locking" - a \c bool specifying whether the algorithm should employ locking of converged eigenpairs. Default: false
   *   - "Max Locked" - a \c int specifying the maximum number of eigenpairs to be locked. Default: problem->getNEV()
   *   - "Locking Quorum" - a \c int specifying the number of eigenpairs that must meet the locking criteria before locking actually occurs. Default: 1
   *   - "Locking Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide locking. Default: 0.1*convergence tolerance
   *   - "Relative Locking Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding locking. Default: true
   */
  BlockKrylovSchurSolMgr( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~BlockKrylovSchurSolMgr() {};
  //@}
  
  //! @name Accessor methods
  //@{ 

  Eigenproblem<ScalarType,MV,OP>& getProblem() const {
    return *_problem;
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
   *    - the locking conditions have been met. In this scenario, some of the current eigenpairs will be removed\n
   *      from the eigensolver and placed into auxiliary storage. The eigensolver will be restarted with the remaining part of the Krylov subspace\n
   *      and some random information to replace the removed subspace.
   *    - global convergence has been met. In this case, the most significant NEV eigenpairs in the solver and locked storage  \n
   *      have met the convergence criterion. (Here, NEV refers to the number of eigenpairs requested by the Eigenproblem.)    \n
   *      In this scenario, the solver manager will return ::Converged.
   *
   * \returns ::ReturnType specifying:
   * <ul>
   *     - ::Converged: the eigenproblem was solved to the specification required by the solver manager.
   *     - ::Unconverged: the eigenproblem was not solved to the specification desired by the solver manager.
   * </ul>
  */
  ReturnType solve();
  //@}

  private:
  Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem;
  Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > _sort;

  string _whch; 

  MagnitudeType _convtol;
  int _maxRestarts;
  bool _relconvtol,_conjSplit;
  int _blockSize, _numBlocks, _stepSize;
  int _verbosity;
};


// Constructor
template<class ScalarType, class MV, class OP>
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::BlockKrylovSchurSolMgr( 
        const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl ) : 
  _problem(problem),
  _whch("LM"),
  _convtol(0),
  _maxRestarts(20),
  _relconvtol(true),
  _conjSplit(false),
  _blockSize(0),
  _numBlocks(0),
  _stepSize(0),
  _verbosity(Anasazi::Errors)
{
  TEST_FOR_EXCEPTION(_problem == Teuchos::null,               std::invalid_argument, "Problem not given to solver manager.");
  TEST_FOR_EXCEPTION(!_problem->isProblemSet(),               std::invalid_argument, "Problem not set.");
  TEST_FOR_EXCEPTION(_problem->getInitVec() == Teuchos::null, std::invalid_argument, "Problem does not contain initial vectors to clone from.");

  // which values to solve for
  _whch = pl.get("Which",_whch);
  if (_whch != "SM" && _whch != "LM" && _whch != "SR" && _whch != "LR") {
    _whch = "LM";
  }

  // convergence tolerance
  _convtol = pl.get("Convergence Tolerance",MT::prec());
  _relconvtol = pl.get("Relative Convergence Tolerance",_relconvtol);
  
  // maximum number of restarts
  _maxRestarts = pl.get("Maximum Restarts",_maxRestarts);

  // block size: default is nev()
  _blockSize = pl.get("Block Size",_problem->getNEV());
  TEST_FOR_EXCEPTION(_blockSize <= 0, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Block Size\" must be strictly positive.");
  _numBlocks = pl.get("Num Blocks",2);
  TEST_FOR_EXCEPTION(_numBlocks <= 1, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Num Blocks\" must be strictly positive.");

  TEST_FOR_EXCEPTION(_numBlocks*_blockSize > MVT::GetVecLength(*_problem->getInitVec()),
                     std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: Potentially impossible orthogonality requests. Reduce basis size or locking size.");

  // step size: the default is _maxRestarts*_numBlocks, so that Ritz values are only computed every restart.
  _stepSize = pl.get("Step Size", _maxRestarts*_numBlocks);
  TEST_FOR_EXCEPTION(_stepSize <= 1, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Step Size\" must be strictly positive.");

  // get the sort manager
  if (pl.isParameter("Sort Manager")) {
    _sort = pl.get<Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > >("Sort Manager");
  } else {
    _sort = Teuchos::rcp( new BasicSort<ScalarType,MV,OP>(_whch) );
  }

  // verbosity level
  _verbosity = pl.get("Verbosity", _verbosity);

}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType 
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::solve() {

  const int nev = _problem->getNEV();

  //////////////////////////////////////////////////////////////////////////////////////
  // Output manager
  Teuchos::RefCountPtr<BasicOutputManager<ScalarType> > printer = Teuchos::rcp( new BasicOutputManager<ScalarType>(_verbosity) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Status tests
  //
  // convergence
  Teuchos::RefCountPtr<StatusTestOrderedResNorm<ScalarType,MV,OP> > convtest 
      = Teuchos::rcp( new StatusTestOrderedResNorm<ScalarType,MV,OP>(_sort,_convtol,nev,false,_relconvtol) );

  // printing StatusTest
  Teuchos::RefCountPtr<StatusTestOutput<ScalarType,MV,OP> > outputtest
    = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer,convtest,1,Passed ) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Orthomanager
  Teuchos::RefCountPtr<SVQBOrthoManager<ScalarType,MV,OP> > ortho 
    = Teuchos::rcp( new SVQBOrthoManager<ScalarType,MV,OP>(_problem->getM()) );

  // utils
  ModalSolverUtils<ScalarType,MV,OP> msutils(printer);

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",_blockSize);
  plist.set("Num Blocks",_numBlocks);
  plist.set("Step Size",_stepSize);

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockKrylovSchur solver
  Teuchos::RefCountPtr<BlockKrylovSchur<ScalarType,MV,OP> > bks_solver 
    = Teuchos::rcp( new BlockKrylovSchur<ScalarType,MV,OP>(_problem,_sort,printer,outputtest,ortho,plist) );
  // set any auxiliary vectors defined in the problem
  Teuchos::RefCountPtr< const MV > probauxvecs = _problem->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    bks_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(probauxvecs) );
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Storage
  // workMV will be used when restarting because of 
  // a) full basis, in which case we need 2*blockSize, for X and KX
  // b) locking, in which case we will need as many vectors as in the current basis, 
  //    a number which will be always <= (numblocks-1)*blocksize
  //    [note: this is because we never lock with curdim == numblocks*blocksize]
  Teuchos::RefCountPtr<MV> workMV;
  workMV = MVT::Clone(*_problem->getInitVec(),2*_blockSize);

  // go ahead and initialize the solution to nothing in case we throw an exception
  Eigensolution<ScalarType,MV> sol;
  sol.numVecs = 0;
  _problem->setSolution(sol);

  int numRestarts = 0;

  // tell bks_solver to iterate
  while (1) {
    try {
      bks_solver->iterate();

      // check convergence first
      if (convtest->getStatus() == Passed ) {
        // we have convergence
        // convtest->whichVecs() tells us which vectors from solver state are the ones we want
        // convtest->howMany() will tell us how many
        break;
      }
      // check for restarting before locking: if we need to lock, it will happen after the restart
      else if ( bks_solver->getCurSubspaceDim() == bks_solver->getMaxSubspaceDim() ) {

        if ( numRestarts >= _maxRestarts ) {
          break; // break from while(1){bks_solver->iterate()}
        }
        numRestarts++;

        printer->stream(Debug) << " Performing restart number " << numRestarts << " of " << _maxRestarts << endl << endl;

        // the solver has filled its basis. 
        std::vector<int> b1ind(_blockSize), b2ind(_blockSize);
        for (int i=0;i<_blockSize;i++) {
          b1ind[i] = i;
          b2ind[i] = _blockSize+i;
        }

	if (bks_solver->isSchurCurrent())
	  cout << "The Schur factorization is current!" << endl;
	else 
	  cout << "The Schur factorization is not current!" << endl;

        BlockKrylovSchurState<ScalarType,MV> newstate;
        //newstate.V = newV;
        bks_solver->initialize(newstate);

      }
      else {
        TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::BlockKrylovSchurSolMgr::solve(): Invalid return from bks_solver::iterate().");
      }
    }
    catch (std::exception e) {
      printer->stream(Errors) << "Error! Caught exception in BlockKrylovSchur::iterate() at iteration " << bks_solver->getNumIters() << endl 
                              << e.what() << endl;
      throw;
    }
  }

  // clear temp space
  workMV = Teuchos::null;

  sol.numVecs = convtest->howMany();
  if (sol.numVecs > 0) {
    sol.Evecs = MVT::Clone(*_problem->getInitVec(),sol.numVecs);
    sol.Espace = sol.Evecs;
    sol.Evals.resize(sol.numVecs);

    // copy them into the solution
    std::vector<int> which = convtest->whichVecs();
    // indices between [0,blockSize) refer to vectors/values in the solver
    // indices between [blockSize,blocksize+numlocked) refer to locked vectors/values
    // everything has already been ordered by the solver; we just have to partition the two references
    std::vector<int> inlocked(0), insolver(0);
    for (unsigned int i=0; i<which.size(); i++) {
      if (which[i] < nev) {
        insolver.push_back(which[i]);
      }
      else {
        // sanity check
        TEST_FOR_EXCEPTION(which[i] > nev,std::logic_error,"Anasazi::BlockKrylovSchurSolMgr::solve(): indexing mistake.");
        inlocked.push_back(which[i] - _blockSize);
      }
    }

    TEST_FOR_EXCEPTION(insolver.size() + inlocked.size() != (unsigned int)sol.numVecs,std::logic_error,"Anasazi::BlockKrylovSchurSolMgr::solve(): indexing mistake.");

    // set the vecs,vals in the solution
    if (insolver.size() > 0) {
      // set vecs
      int lclnum = insolver.size();
      std::vector<int> tosol(lclnum);
      for (int i=0; i<lclnum; i++) tosol[i] = i;
      Teuchos::RefCountPtr<const MV> v = MVT::CloneView(*bks_solver->getRitzVectors(),insolver);
      MVT::SetBlock(*v,tosol,*sol.Evecs);
      // set vals
      std::vector<MagnitudeType> fromsolver = bks_solver->getRitzValues();
      for (unsigned int i=0; i<insolver.size(); i++) {
        sol.Evals[i] = fromsolver[insolver[i]];
      }
    }

    // setup sol.index, remembering that all eigenvalues are real so that index = {0,...,0}
    sol.index.resize(sol.numVecs,0);

    // sort the eigenvalues and permute the eigenvectors appropriately
    {
      std::vector<int> order(sol.numVecs);
      std::vector<ScalarType> vals_st(sol.numVecs);
      std::copy(sol.Evals.begin(),sol.Evals.end(),vals_st.begin());
      _sort->sort( NULL, sol.numVecs, &vals_st[0], &order );
      for (int i=0; i<sol.numVecs; i++) {
        sol.Evals[i] = SCT::real( vals_st[i] );
      }
      // now permute the eigenvectors according to order
      msutils.permuteVectors(sol.numVecs,order,*sol.Evecs);
    }
  }

  // print final summary
  bks_solver->currentStatus(printer->stream(FinalSummary));

  // print timing information
  Teuchos::TimeMonitor::summarize(printer->stream(TimingDetails));

  _problem->setSolution(sol);
  printer->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << endl;

  if (sol.numVecs < nev) {
    return Unconverged; // return from BlockKrylovSchurSolMgr::solve() 
  }
  return Converged; // return from BlockKrylovSchurSolMgr::solve() 
}


} // end Anasazi namespace

#endif /* ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP */
