
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

#ifndef ANASAZI_BLOCKDAVIDSON_SOLMGR_HPP
#define ANASAZI_BLOCKDAVIDSON_SOLMGR_HPP

/*! \file AnasaziBlockDavidsonSolMgr.hpp
 *  \brief The Anasazi::BlockDavidsonSolMgr provides a powerful solver manager for the BlockDavidson eigensolver.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziModalSolverUtils.hpp"

#include "AnasaziBlockDavidson.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestOrderedResNorm.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "Teuchos_BLAS.hpp"

/*! \class Anasazi::BlockDavidsonSolMgr
 *
 *  \brief The Anasazi::BlockDavidsonSolMgr provides a powerful and fully-featured solver manager over the BlockDavidson eigensolver.
 *
 *  Features include: FINISH
 *  <ul>
 *  <li>
 *  </ul>
 *
 */

namespace Anasazi {

template<class ScalarType, class MV, class OP>
class BlockDavidsonSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:

  //@{ \name Constructors/Destructor.

  //! Default Constructor.
  BlockDavidsonSolMgr() {};

  //! Basic Constructor.
  BlockDavidsonSolMgr( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~BlockDavidsonSolMgr() {};
  //@}
  
  //@{ \name Accessor methods

  Eigenproblem<ScalarType,MV,OP>& getProblem() const {
    return *_problem;
  }

  //@}

  //@{ \name Solver application methods.
    
  /*! \brief This method performs possibly repeated calls to the underlying eigensolver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
   * quit.
   *
   * In the case of BlockDavidsonSolMgr, FINISH
   *
   * \returns ReturnType specifying:
   * <ul>
   *    <li>Converged: the eigenproblem was solved to the specification required by the solver manager.
   *    <li>Unconverged: the eigenproblem was not solved to the specification desired by the solver manager.
   * </ul>
  */
  ReturnType solve();
  //@}

  private:
  Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem;

  string _whch; 

  MagnitudeType _convtol, _locktol;
  int _maxRestarts;
  bool _useLocking;
  bool _relconvtol, _rellocktol;
  int _blockSize, _numBlocks;
  int _maxLocked;
  int _verbosity;
  int _lockQuorum;
};


// Constructor
template<class ScalarType, class MV, class OP>
BlockDavidsonSolMgr<ScalarType,MV,OP>::BlockDavidsonSolMgr( 
        const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl ) : 
  _problem(problem),
  _whch("SR"),
  _convtol(MT::prec()),
  _locktol(MT::prec()/10),
  _maxRestarts(20),
  _useLocking(false),
  _relconvtol(true),
  _rellocktol(true),
  _blockSize(0),
  _numBlocks(0),
  _maxLocked(0),
  _verbosity(Anasazi::Errors),
  _lockQuorum(1)
{
  TEST_FOR_EXCEPTION(_problem == Teuchos::null, AnasaziError, "Problem not given to solver manager.");
  TEST_FOR_EXCEPTION(!_problem->isProblemSet(), AnasaziError, "Problem not set.");
  TEST_FOR_EXCEPTION(!_problem->isHermitian(), AnasaziError, "Problem not symmetric.");

  // which values to solve for
  _whch = pl.get("Which",_whch);
  if (_whch != "SM" && _whch != "LM" && _whch != "SR" && _whch != "LR") {
    _whch = "SR";
  }

  // convergence tolerance
  _convtol = pl.get("Convergence Tolerance",_convtol);
  _relconvtol = pl.get("Relative Convergence Tolerance",_relconvtol);
  
  // locking tolerance
  _useLocking = pl.get("Use Locking",_useLocking);
  _rellocktol = pl.get("Relative Locking Tolerance",_rellocktol);
  _locktol = pl.get("Locking Tolerance",_locktol);

  // maximum number of restarts
  _maxRestarts = pl.get("Maximum Restarts",_maxRestarts);

  // block size: default is nev()
  _blockSize = pl.get("Block Size",_problem->getNEV());
  TEST_FOR_EXCEPTION(_blockSize <= 0, std::invalid_argument,
                     "Anasazi::BlockDavidsonSolMgr: \"Block Size\" must be strictly positive.");
  _numBlocks = pl.get("Num Blocks",1);
  TEST_FOR_EXCEPTION(_numBlocks <= 0, std::invalid_argument,
                     "Anasazi::BlockDavidsonSolMgr: \"Num Blocks\" must be strictly positive.");

  // max locked: default is nev(), must satisfy _maxLocked + _blockSize >= nev
  if (_useLocking) {
    _maxLocked = pl.get("Max Locked",_problem->getNEV());
  }
  else {
    _maxLocked = 0;
  }
  if (_maxLocked == 0) {
    _useLocking = false;
  }
  TEST_FOR_EXCEPTION(_maxLocked < 0, std::invalid_argument,
                     "Anasazi::BlockDavidsonSolMgr: \"Max Locked\" must be positive.");
  TEST_FOR_EXCEPTION(_maxLocked + _blockSize < _problem->getNEV(), 
                     std::logic_error,
                     "Anasazi::BlockDavidsonSolMgr: Not enough storage space for requested number of eigenpairs.");
  TEST_FOR_EXCEPTION(_numBlocks*_blockSize + _maxLocked > MVT::GetVecLength(*_problem->getInitVec()),
                     std::logic_error,
                     "Anasazi::BlockDavidsonSolMgr: Potentially impossible orthogonality requests. Reduce basis size or locking size.");

  if (_useLocking) {
    _lockQuorum = pl.get("Locking Quorum",_lockQuorum);
    TEST_FOR_EXCEPTION(_lockQuorum <= 0,
                       std::logic_error,
                       "Anasazi::BlockDavidsonSolMgr: \"Locking Quorum\" must be strictly positive.");
  }

  // verbosity level
  _verbosity = pl.get("Verbosity", _verbosity);

}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType 
BlockDavidsonSolMgr<ScalarType,MV,OP>::solve() {

  const int nev = _problem->getNEV();

  //////////////////////////////////////////////////////////////////////////////////////
  // Sort manager
  Teuchos::RefCountPtr<BasicSort<ScalarType,MV,OP> > sorter = Teuchos::rcp( new BasicSort<ScalarType,MV,OP>(_whch) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Output manager
  Teuchos::RefCountPtr<BasicOutputManager<ScalarType> > printer = Teuchos::rcp( new BasicOutputManager<ScalarType>(_verbosity) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Status tests
  //
  // convergence
  Teuchos::RefCountPtr<StatusTestOrderedResNorm<ScalarType,MV,OP> > convtest 
      = Teuchos::rcp( new StatusTestOrderedResNorm<ScalarType,MV,OP>(sorter,_convtol,nev,false,_relconvtol) );
  // locking
  Teuchos::RefCountPtr<StatusTestResNorm<ScalarType,MV,OP> > locktest;
  if (_useLocking) {
    locktest = Teuchos::rcp( new StatusTestResNorm<ScalarType,MV,OP>(_locktol,_lockQuorum,false,_rellocktol) );
  }
  // combo class
  Teuchos::Array<Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > > alltests;
  // for an OR test, the order doesn't matter
  alltests.push_back(convtest);
  if (locktest != Teuchos::null)   alltests.push_back(locktest);
  // combo: convergence || locking 
  Teuchos::RefCountPtr<StatusTestCombo<ScalarType,MV,OP> > combotest
    = Teuchos::rcp( new StatusTestCombo<ScalarType,MV,OP>( StatusTestCombo<ScalarType,MV,OP>::OR, alltests) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Orthomanager
  Teuchos::RefCountPtr<SVQBOrthoManager<ScalarType,MV,OP> > ortho 
    = Teuchos::rcp( new SVQBOrthoManager<ScalarType,MV,OP>(_problem->getM()) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",_blockSize);
  plist.set("Num Blocks",_numBlocks);

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockDavidson solver
  Teuchos::RefCountPtr<BlockDavidson<ScalarType,MV,OP> > bd_solver 
    = Teuchos::rcp( new BlockDavidson<ScalarType,MV,OP>(_problem,sorter,printer,combotest,ortho,plist) );
  // set any auxilliary vectors defined in the problem
  Teuchos::RefCountPtr< const MV > probauxvecs = _problem->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    bd_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(probauxvecs) );
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Storage
  int numlocked = 0;
  Teuchos::RefCountPtr<MV> lockvecs;
  if (_maxLocked > 0) {
    lockvecs = MVT::Clone(*_problem->getInitVec(),_maxLocked);
  }
  std::vector<MagnitudeType> lockvals;

  // go ahead and initialize the solution to nothing in case we throw an exception
  Eigensolution<ScalarType,MV> sol;
  sol.numVecs = 0;
  _problem->setSolution(sol);

  int numRestarts = 0;

  // tell bd_solver to iterate
  while (1) {
    try {
      bd_solver->iterate();

      // check convergence first
      if (convtest->getStatus() == Passed ) {
        // we have convergence or not
        // convtest->whichVecs() tells us which vectors from lockvecs and solver->getEvecs() are the ones we want
        // convtest->howMany() will tell us how many
        break;
      }
      // check locking if we didn't converge
      else if (locktest != Teuchos::null && locktest->getStatus() == Passed) {

        // remove the locked vectors,values from bd_solver: put them in newvecs, newvals
        int numnew = locktest->howMany();
        TEST_FOR_EXCEPTION(numnew <= 0,std::logic_error,"Anasazi::BlockDavidsonSolMgr::solve(): status test mistake.");
        // get the indices
        std::vector<int> newind = locktest->whichVecs();

        // don't lock more than _maxLocked; we didn't allocate enough space.
        if (numlocked + numnew > _maxLocked) {
          numnew = _maxLocked - numlocked;
          newind.resize(numnew);
        }

        {
          // debug printing
          printer->print(Debug,"Locking vectors: ");
          for (unsigned int i=0; i<newind.size(); i++) {printer->stream(Debug) << " " << newind[i];}
          printer->print(Debug,"\n");
        }
        std::vector<MagnitudeType> newvals(numnew);
        Teuchos::RefCountPtr<const MV> newvecs;
        {
          // get the vectors
          newvecs = MVT::CloneView(*bd_solver->getEvecs(),newind);
          // get the values
          std::vector<MagnitudeType> allvals = bd_solver->getEigenvalues();
          for (int i=0; i<numnew; i++) {
            newvals[i] = allvals[newind[i]];
          }
        }
        // put newvecs into lockvecs
        {
          std::vector<int> indlock(numnew);
          for (int i=0; i<numnew; i++) indlock[i] = numlocked+i;
          MVT::SetBlock(*newvecs,indlock,*lockvecs);
          newvecs = Teuchos::null;
        }
        // put newvals into lockvals
        lockvals.insert(lockvals.end(),newvals.begin(),newvals.end());
        numlocked += numnew;
        // add locked vecs as aux vecs, along with aux vecs from problem
        {
          std::vector<int> indlock(numlocked);
          for (int i=0; i<numlocked; i++) indlock[i] = i;
          Teuchos::RefCountPtr<const MV> curlocked = MVT::CloneView(*lockvecs,indlock);
          if (probauxvecs != Teuchos::null) {
            bd_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(probauxvecs,curlocked) );
          }
          else {
            bd_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(curlocked) );
          }
        }
        // add locked vals to convtest
        convtest->setAuxVals(lockvals);

        //
        // restart the solver
        //
        // finish
        // we have to get the projected stiffness matrix from the solver, compute the projected eigenvectors and sort them
        // then all of the ones that don't correspond to the ones we wish to lock get orthogonalized using QR factorization
        // then we generate new (slightly smaller basis) and augment this basis with random information to preserve rank 
        // (which must be a muliple of blockSize) and hand this to the solver

        Teuchos::BLAS<int,ScalarType> blas;
        {
          const ScalarType ONE = SCT::one();
          const ScalarType ZERO = SCT::zero();

          BlockDavidsonState<ScalarType,MV> state = bd_solver->getState();
          // don't need the following
          state.X = Teuchos::null;
          state.KX = Teuchos::null;
          state.MX = Teuchos::null;
          state.R = Teuchos::null;
          state.H = Teuchos::null;
          state.T = Teuchos::null;
          //
          // get current size of basis, build index, and get a view of the current basis and project stiffness matrix
          int curdim = state.curDim;
          std::vector<int> curind(curdim);
          for (int i=0; i<curdim; i++) curind[i] = i;
          Teuchos::RefCountPtr<const MV> curV = MVT::CloneView(*state.V,curind);
          Teuchos::RefCountPtr<const Teuchos::SerialDenseMatrix<int,ScalarType> > curKK;
          curKK = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*state.KK,curdim,curdim) );
          //
          // compute eigenvectors of the projected stiffness matrix
          ModalSolverUtils<ScalarType,MV,OP> msutils(printer);
          Teuchos::SerialDenseMatrix<int,ScalarType> S(curdim,curdim);
          std::vector<MagnitudeType> theta(curdim);
          int rank = curdim;
          msutils.directSolver(curdim,*curKK,0,&S,&theta,&rank,10);
          TEST_FOR_EXCEPTION(rank != curdim,std::logic_error,"Anasazi::BlockDavidsonSolMgr::solve(): direct solve did not compute all eigenvectors."); // this should never happen
          // 
          // sort the eigenvalues (so that we can order the eigenvectors): finish
          //
          // select the non-locked eigenvectors
          std::vector<int> unlockind(curdim-numnew);
          set_difference(curind.begin(),curind.end(),newind.begin(),newind.end(),unlockind.begin());
          Teuchos::SerialDenseMatrix<int,ScalarType> Sunlocked(curdim,curdim-numnew);
          for (int i=0; i<curdim-numnew; i++) {
            blas.COPY(curdim, S[unlockind[i]], 1, Sunlocked[i], 1);
          }
          // 
          // compute the qr factorization: finish
          //
          // allocate space for the new basis, compute it
          Teuchos::RefCountPtr<MV> newV = MVT::Clone(*curV,curdim);
          Teuchos::RefCountPtr<MV> genV;
          { 
            std::vector<int> genind(curdim-numnew);
            for (int i=0; i<curdim-numnew; i++) genind[i] = i;
            genV = MVT::CloneView(*newV,genind);
            MVT::MvTimesMatAddMv(ONE,*curV,Sunlocked,ZERO,*genV);
          }
          //
          // allocate space for the new projected stiffness matrix, compute leading block
          Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > newKK;
          newKK = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(curdim,curdim) );
          {
            Teuchos::SerialDenseMatrix<int,ScalarType> tmpKK(curdim,curdim-numnew),
                                                       newKK11(Teuchos::View,*newKK,curdim-numnew,curdim-numnew),
                                                       curKKsym(*curKK);
            for (int j=0; j<curdim; j++) {
              for (int i=j+1; i<curdim; i++) {
                curKKsym(i,j) = curKKsym(j,i);
              }
            }
            tmpKK.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,curKKsym,Sunlocked,ZERO);
            newKK11.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ONE,Sunlocked,tmpKK,ZERO);
          }
          //
          // generate random data to fill the rest of the basis back to curdim
          Teuchos::RefCountPtr<MV> augV;
          {
            std::vector<int> augind(curdim-numnew);
            for (int i=0; i<numnew; i++) augind[i] = curdim-numnew+i;
            augV = MVT::CloneView(*newV,augind);
            MVT::MvRandom(*augV);
          }
          // 
          // orthogonalize it against auxvecs and the current basis
          {
            Teuchos::Array<Teuchos::RefCountPtr<const MV> > against = bd_solver->getAuxVecs();
            against.push_back(genV);
            ortho->projectAndNormalize(*augV,Teuchos::null,
                                       Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null),Teuchos::null,
                                       against);

          }
          //
          // project the stiffness matrix on the new part
          {
            Teuchos::RefCountPtr<MV> augKV = MVT::Clone(*augV,curdim-numnew);
            OPT::Apply(*_problem->getOperator(),*augV,*augKV);
            Teuchos::SerialDenseMatrix<int,ScalarType> KK12(Teuchos::View,*newKK,curdim-numnew,numnew,0,curdim-numnew),
                                                       KK22(Teuchos::View,*newKK,numnew,numnew,curdim-numnew,curdim-numnew);
            MVT::MvTransMv(ONE,*genV,*augKV,KK12);
            MVT::MvTransMv(ONE,*augV,*augKV,KK22);
          }
          //
          // pass all of this to the solver
          genV = augV = Teuchos::null;
          state.V = newV;
          state.KK = newKK;
          state.curDim = curdim;
          bd_solver->initialize(state);
        }

        if (numlocked == _maxLocked) {
          // disabled locking now by setting quorum to unreachable number
          locktest->setQuorum(_blockSize*_numBlocks+1);
        }

      }
      // didn't lock, didn't converge; it must be time to restart
      // however, if numBlocks == 1, we never restart, and solver should not have quit.
      else if ( bd_solver->getCurSubspaceDim() == bd_solver->getMaxSubspaceDim() ) {

        if ( numRestarts >= _maxRestarts ) {
          break; // break from while(1){bd_solver->iterate}
        }
        numRestarts++;

        printer->stream(Debug) << " Performing restart number " << numRestarts << " of " << _maxRestarts << endl << endl;

        // the solver has filled its basis. 
        // the current eigenvectors will be used to restart the basis.
        Teuchos::RefCountPtr<MV> newV, newKV, newMV;
        { 
          // we want curstate (and its pointers) to go out of scope and be deleted
          BlockDavidsonState<ScalarType,MV> curstate = bd_solver->getState();
          newV = MVT::CloneCopy(*curstate.X);
        }

        if (_problem->getM() != Teuchos::null) {
          newMV = MVT::Clone(*newV,_blockSize);
          OPT::Apply(*_problem->getM(),*newV,*newMV);
        }
        else {
          newMV = Teuchos::null;
        }

        // send this basis to the orthomanager to ensure orthonormality
        Teuchos::Array<Teuchos::RefCountPtr<const MV> > curauxvecs = bd_solver->getAuxVecs();
        Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
        ortho->projectAndNormalize(*newV,newMV,dummy,Teuchos::null,curauxvecs);

        // compute K*newV
        newKV = MVT::Clone(*newV,_blockSize);
        OPT::Apply(*_problem->getOperator(),*newV,*newKV);

        // compute projected stiffness matrix
        Teuchos::RefCountPtr< Teuchos::SerialDenseMatrix<int,ScalarType> > 
            newKK = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(_blockSize,_blockSize) );
        MVT::MvTransMv(SCT::one(),*newV,*newKV,*newKK);

        // initialize() will do the rest
        BlockDavidsonState<ScalarType,MV> newstate;
        newstate.V = newV;
        newstate.KK = newKK;
        bd_solver->initialize(newstate);
      }
      else {
        TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::BlockDavidsonSolMgr::solve(): Invalid return from bd_solver::iterate().");
      }
    }
    catch (std::exception e) {
      printer->stream(Errors) << "Error! Caught exception in BlockDavidson::iterate() at iteration " << bd_solver->getNumIters() << endl 
                              << e.what() << endl;
      throw;
    }
  }

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
      if (which[i] < _blockSize) {
        insolver.push_back(which[i]);
      }
      else {
        // sanity check
        TEST_FOR_EXCEPTION(which[i] >= numlocked+_blockSize,std::logic_error,"Anasazi::BlockDavidsonSolMgr::solve(): indexing mistake.");
        inlocked.push_back(which[i] - _blockSize);
      }
    }

    TEST_FOR_EXCEPTION(insolver.size() + inlocked.size() != (unsigned int)sol.numVecs,std::logic_error,"Anasazi::BlockDavidsonSolMgr::solve(): indexing mistake.");

    // set the vecs,vals in the solution
    if (insolver.size() > 0) {
      // set vecs
      int lclnum = insolver.size();
      std::vector<int> tosol(lclnum);
      for (int i=0; i<lclnum; i++) tosol[i] = i;
      Teuchos::RefCountPtr<const MV> v = MVT::CloneView(*bd_solver->getEvecs(),insolver);
      MVT::SetBlock(*v,tosol,*sol.Evecs);
      // set vals
      std::vector<MagnitudeType> fromsolver = bd_solver->getEigenvalues();
      for (unsigned int i=0; i<insolver.size(); i++) {
        sol.Evals[i] = fromsolver[insolver[i]];
      }
    }

    // get the vecs,vals from locked storage
    if (inlocked.size() > 0) {
      int solnum = insolver.size();
      // set vecs
      int lclnum = inlocked.size();
      std::vector<int> tosol(lclnum);
      for (int i=0; i<lclnum; i++) tosol[i] = solnum + i;
      Teuchos::RefCountPtr<const MV> v = MVT::CloneView(*lockvecs,inlocked);
      MVT::SetBlock(*v,tosol,*sol.Evecs);
      // set vals
      for (unsigned int i=0; i<inlocked.size(); i++) {
        sol.Evals[i+solnum] = lockvals[inlocked[i]];
      }
    }

    // setup sol.index, remembering that all eigenvalues are real so that index = {0,...,0}
    sol.index.resize(sol.numVecs,0);
  }

  // print final summary
  bd_solver->currentStatus(printer->stream(FinalSummary));

  // print timing information
  Teuchos::TimeMonitor::summarize(printer->stream(TimingDetails));

  _problem->setSolution(sol);
  printer->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << endl;

  if (sol.numVecs < nev) {
    return Unconverged; // return from BlockDavidsonSolMgr::solve() 
  }
  return Converged; // return from BlockDavidsonSolMgr::solve() 
}


} // end Anasazi namespace

#endif /* ANASAZI_BLOCKDAVIDSON_SOLMGR_HPP */
