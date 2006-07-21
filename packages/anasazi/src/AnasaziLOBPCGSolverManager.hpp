
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

#ifndef ANASAZI_LOBPCG_SOLVERMANAGER_HPP
#define ANASAZI_LOBPCG_SOLVERMANAGER_HPP

/*! \file AnasaziLOBPCGSolverManager.hpp
 *  \brief The Anasazi::LOBPCGSolverManager provides a powerful solver manager for the LOBPCG eigensolver.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"

#include "AnasaziLOBPCG.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestOrderedResNorm.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziBasicOutputManager.hpp"

/*! \class Anasazi::LOBPCGSolverManager
 *
 *  \brief The Anasazi::LOBPCGSolverManager provides a powerful and fully-featured solver manager over the LOBPCG eigensolver.
 *
 *  Features include: FINISH
 *  <ul>
 *  <li>
 *  </ul>
 *
 */

namespace Anasazi {

template<class ScalarType, class MV, class OP>
class LOBPCGSolverManager : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:

  //@{ \name Constructors/Destructor.

  //! Default Constructor.
  LOBPCGSolverManager() {};

  //! Basic Constructor.
  LOBPCGSolverManager( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~LOBPCGSolverManager() {};
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
   * In the case of LOBPCGSolverManager, FINISH
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
  int _maxIters;
  bool _useLocking;
  bool _relconvtol, _rellocktol;
  int _blockSize;
  bool _fullOrtho;
  int _maxLocked;
  int _verbosity;
  int _lockQuorum;
};


// Constructor
template<class ScalarType, class MV, class OP>
LOBPCGSolverManager<ScalarType,MV,OP>::LOBPCGSolverManager( 
        const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl ) : 
  _problem(problem),
  _whch("SR"),
  _convtol(MT::prec()),
  _locktol(MT::prec()/10),
  _maxIters(100),
  _useLocking(false),
  _relconvtol(true),
  _rellocktol(true),
  _blockSize(0),
  _fullOrtho(true),
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

  // maximum number of iterations
  _maxIters = pl.get("Maximum Iterations",_maxIters);

  // block size: default is nev()
  _blockSize = pl.get("Block Size",_problem->getNEV());
  TEST_FOR_EXCEPTION(_blockSize <= 0, std::invalid_argument,
                     "Anasazi::LOBPCGSolverManager: \"Block Size\" must be strictly positive.");

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
                     "Anasazi::LOBPCGSolverManager: \"Max Locked\" must be positive.");
  TEST_FOR_EXCEPTION(_maxLocked + _blockSize < _problem->getNEV(), 
                     std::logic_error,
                     "Anasazi::LOBPCGSolverManager: Not enough storage space for requested number of eigenpairs.");

  if (_useLocking) {
    _lockQuorum = pl.get("Locking Quorum",_lockQuorum);
    TEST_FOR_EXCEPTION(_lockQuorum <= 0,
                       std::logic_error,
                       "Anasazi::LOBPCGSolverManager: \"Locking Quorum\" must be strictly positive.");
  }

  // full orthogonalization: default true
  _fullOrtho = pl.get("Full Ortho",_fullOrtho);

  // verbosity level
  _verbosity = pl.get("Verbosity", _verbosity);
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType 
LOBPCGSolverManager<ScalarType,MV,OP>::solve() {

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
  // maximum number of iterations: optional test
  Teuchos::RefCountPtr<StatusTestMaxIters<ScalarType,MV,OP> > maxtest;
  if (_maxIters > 0) {
    maxtest = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>(_maxIters) );
  }
  // convergence
  Teuchos::RefCountPtr<StatusTestOrderedResNorm<ScalarType,MV,OP> > convtest 
      = Teuchos::rcp( new StatusTestOrderedResNorm<ScalarType,MV,OP>(sorter,_convtol,nev,false,_relconvtol) );
  // locking
  Teuchos::RefCountPtr<StatusTestResNorm<ScalarType,MV,OP> > locktest;
  if (_useLocking) {
    locktest = Teuchos::rcp( new StatusTestResNorm<ScalarType,MV,OP>(_locktol,_lockQuorum,false,_rellocktol) );
  }
  Teuchos::Array<Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > > alltests;
  // for an OR test, the order doesn't matter
  alltests.push_back(convtest);
  if (maxtest != Teuchos::null) alltests.push_back(maxtest);
  if (locktest != Teuchos::null)   alltests.push_back(locktest);
  // combo: convergence || locking || max iters
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
  plist.set("Full Ortho",_fullOrtho);

  //////////////////////////////////////////////////////////////////////////////////////
  // LOBPCG solver
  Teuchos::RefCountPtr<LOBPCG<ScalarType,MV,OP> > lobpcg_solver 
    = Teuchos::rcp( new LOBPCG<ScalarType,MV,OP>(_problem,sorter,printer,combotest,ortho,plist) );
  // set any auxilliary vectors defined in the problem
  Teuchos::RefCountPtr< const MV > probauxvecs = _problem->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    lobpcg_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(probauxvecs) );
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

  // tell the lobpcg_solver to iterate
  while (1) {
    try {
      lobpcg_solver->iterate();

      // check convergence first
      if (convtest->getStatus() == Passed || (maxtest != Teuchos::null && maxtest->getStatus() == Passed) ) {
        // we have convergence or not
        // convtest->whichVecs() tells us which vectors from lockvecs and solver->getEvecs() are the ones we want
        // convtest->howMany() will tell us how many
        break;
      }
      // check locking if we didn't converge
      else if (locktest != Teuchos::null && locktest->getStatus() == Passed) {

        // remove the locked vectors,values from lobpcg_solver: put them in newvecs, newvals
        int numnew = locktest->howMany();
        TEST_FOR_EXCEPTION(numnew <= 0,std::logic_error,"Anasazi::LOBPCGSolverManager::solve(): status test mistake.");
        // get the indices
        std::vector<int> indnew = locktest->whichVecs();
        {
          // debug printing
          printer->print(Debug,"Locking vectors: ");
          for (unsigned int i=0; i<indnew.size(); i++) {printer->stream(Debug) << " " << indnew[i];}
          printer->print(Debug,"\n");
        }
        std::vector<MagnitudeType> newvals(numnew);
        Teuchos::RefCountPtr<const MV> newvecs;
        {
          // work in a local scope, to hide the variabes needed for extracting this info
          // get the vectors
          newvecs = MVT::CloneView(*lobpcg_solver->getEvecs(),indnew);
          // get the values
          std::vector<MagnitudeType> allvals = lobpcg_solver->getEigenvalues();
          for (int i=0; i<numnew; i++) {
            newvals[i] = allvals[indnew[i]];
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
            lobpcg_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(probauxvecs,curlocked) );
          }
          else {
            lobpcg_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(curlocked) );
          }
        }
        // add locked vals to convtest
        convtest->setAuxVals(lockvals);
        // fill out the empty state in the solver
        {
          LOBPCGState<ScalarType,MV> state = lobpcg_solver->getState();
          // don't need the following: KX, KP, R, T, H, KH, MH
          // if hasP(), then ortho it against new aux vecs (and maybe against X); otherwise, it is invalid, so wipe it
          state.R = Teuchos::null;
          state.T = Teuchos::null;
          state.KX = Teuchos::null;
          state.KP = Teuchos::null;
          state.H = Teuchos::null;
          state.KH = Teuchos::null;
          state.MH = Teuchos::null;
          // make copies of the X state and clear the old X state
          Teuchos::RefCountPtr<MV> newX = MVT::CloneCopy(*state.X);
          state.X = Teuchos::null;
          Teuchos::RefCountPtr<MV> newMX, newnewMX;
          // get select part, set to random
          Teuchos::RefCountPtr<MV> newnewX = MVT::CloneView(*newX,indnew);
          MVT::MvRandom(*newnewX);
          if (state.MX != Teuchos::null) {
            newMX = MVT::CloneCopy(*state.MX);
            state.MX = Teuchos::null;
            newnewMX = MVT::CloneView(*newMX,indnew);
            // apply M to randomly generated new bit
            OPT::Apply(*_problem->getM(),*newnewX,*newnewMX);
          }
          // clear the views
          newnewX = Teuchos::null;
          newnewMX = Teuchos::null;

          Teuchos::Array<Teuchos::RefCountPtr<const MV> > curauxvecs = lobpcg_solver->getAuxVecs();
          Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
          // ortho X against the aux vectors
          ortho->projectAndNormalize(*newX,newMX,dummy,Teuchos::null,curauxvecs);
          // save it
          state.X = newX;
          state.MX = newMX;

          if (lobpcg_solver->hasP()) {
            Teuchos::RefCountPtr<MV> newMP, newP = MVT::CloneCopy(*state.P);
            if (state.MP != Teuchos::null) {
              newMP = MVT::CloneCopy(*state.MP);
            }
            if (_fullOrtho) {
              // ortho P against aux vectors and new X
              curauxvecs.push_back(newX);
              ortho->projectAndNormalize(*newP,newMP,dummy,Teuchos::null,curauxvecs);
            }
            else {
              // ortho P against aux vectors
              ortho->projectAndNormalize(*newP,newMP,dummy,Teuchos::null,curauxvecs);
            }
            state.P = newP;
            state.MP = newMP;
          }
          else {
            state.P = Teuchos::null;
            state.MP = Teuchos::null;
          }
          // set the new state
          lobpcg_solver->initialize(state);
        }

        if (numlocked == _maxLocked) {
          // disabled locking now
          locktest->setQuorum(_blockSize+1);
        }
      }
      else {
        TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::LOBPCGSolverManager::solve(): Invalid return from lobpcg_solver::iterate().");
      }
    }
    catch (LOBPCGRitzFailure re) {
      if (_fullOrtho) {
        // if we are already using full orthogonalization, there isn't much we can do here. 
        // the most recent information in the status tests is still valid, and can be used to extract/return the 
        // eigenpairs that have converged.
        break; // while(1)
      }
      printer->stream(Errors) << "Error! Caught LOBPCGRitzFailure at iteration " << lobpcg_solver->getNumIters() << endl
                              << "Full orthogonalization is off; will try to recover." << endl;
      // otherwise, get the current "basis" from the solver, orthonormalize it, do a rayleigh-ritz, and restart with the ritz vectors
      // if there aren't enough, break and quit with what we have
      LOBPCGState<ScalarType,MV> curstate = lobpcg_solver->getState();
      Teuchos::RefCountPtr<MV> restart, Krestart, Mrestart;
      int localsize = lobpcg_solver->hasP() ? 3*_blockSize : 2*_blockSize;
      bool hasM = _problem->getM() != Teuchos::null;
      restart  = MVT::Clone(*curstate.X, localsize);
      std::vector<int> ind(_blockSize);
      // set restart = [X H P] and Krestart = [X H P]
      // put X into [0 , blockSize)
      for (int i=0; i < _blockSize; i++) ind[i] = i;
      MVT::SetBlock(*curstate.X,ind,*restart);
      // put H into [_blockSize , 2*blockSize)
      for (int i=0; i < _blockSize; i++) ind[i] = _blockSize+i;
      MVT::SetBlock(*curstate.H,ind,*restart);
      // optionally, put P into [2*blockSize,3*blockSize)
      if (localsize == 3*_blockSize) {
        for (int i=0; i < _blockSize; i++) ind[i] = 2*_blockSize+i;
        MVT::SetBlock(*curstate.P,ind,*restart);
      }
      // project against auxvecs and locked vecs, and orthonormalize the basis
      if (hasM) {
        Mrestart = MVT::Clone(*restart, localsize);
        OPT::Apply(*_problem->getM(),*restart,*Mrestart);
      }
      else {
        Mrestart = restart;
      }
      Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
      Teuchos::Array<Teuchos::RefCountPtr<const MV> > Q;
      {
        if (numlocked > 0) {
          std::vector<int> indlock(numlocked);
          for (int i=0; i<numlocked; i++) indlock[i] = i;
          Teuchos::RefCountPtr<const MV> curlocked = MVT::CloneView(*lockvecs,indlock);
          Q.push_back(curlocked);
        }
        if (probauxvecs != Teuchos::null) {
          Q.push_back(probauxvecs);
        }
      }
      int rank = ortho->projectAndNormalize(*restart,Mrestart,dummy,Teuchos::null,Q);
      if (rank < _blockSize) {
        // quit
        printer->stream(Errors) << "Error! Recovered basis only rank " << rank << ". Block size is " << _blockSize << ".\n"
                                << "Recovery failed." << endl;
        throw;
      }
      // reduce multivec size if necessary
      if (rank < localsize) {
        localsize = rank;
        ind.resize(localsize);
        for (int i=0; i<localsize; i++) ind[i] = i;
        // grab the first part of restart,Mrestart
        restart = MVT::CloneView(*restart,ind);
        if (hasM) {
          Mrestart = MVT::CloneView(*Mrestart,ind);
        }
        else {
          Mrestart = restart;
        }
      }
      // project the matrices
      Krestart = MVT::Clone(*restart, localsize);
      OPT::Apply(*_problem->getOperator(),*restart,*Krestart);
      Teuchos::SerialDenseMatrix<int,ScalarType> KK(localsize,localsize), MM(localsize,localsize), S(localsize,localsize);
      std::vector<MagnitudeType> theta(localsize);
      // KK = restart^H K restart
      MVT::MvTransMv(1.0,*restart,*Krestart,KK);
      // MM = restart^H M restart
      MVT::MvTransMv(1.0,*restart,*Mrestart,MM);
      rank = localsize;
      ModalSolverUtils<ScalarType,MV,OP> msutils(printer);
      msutils.directSolver(localsize,KK,&MM,&S,&theta,&rank,1);
      if (rank < _blockSize) {
        printer->stream(Errors) << "Error! Recovered basis of rank " << rank << " produced only " << rank << "ritz vectors.\n"
                                << "Block size is " << _blockSize << ".\n"
                                << "Recovery failed." << endl;
        throw;
      }
      theta.resize(rank);
      // sort the ritz values using the sort manager
      {
        Teuchos::BLAS<int,ScalarType> blas;
        std::vector<int> order(rank);
        std::vector<ScalarType> theta_st(theta.size());
        std::copy(theta.begin(),theta.end(),theta_st.begin());
        sorter->sort( lobpcg_solver.get(), rank, &(theta_st[0]), &order );   // don't catch exception
        
        //  Reorder theta according to sorting results from theta_st
        std::vector<MagnitudeType> theta_copy(theta);
        for (int i=0; i<rank; i++) {
          theta[i] = theta_copy[order[i]];
        }
        // Sort the primitive ritz vectors
        Teuchos::SerialDenseMatrix<int,ScalarType> copyS( S );
        for (int i=0; i<rank; i++) {
          blas.COPY(localsize, copyS[order[i]], 1, S[i], 1);
        }
      }
      //
      Teuchos::SerialDenseMatrix<int,ScalarType> S1(Teuchos::View,S,localsize,_blockSize);
      // compute the ritz vectors
      LOBPCGState<ScalarType,MV> newstate;
      Teuchos::RefCountPtr<MV> newX, newKX, newMX;
      // X
      newX = MVT::Clone(*restart,_blockSize);
      MVT::MvTimesMatAddMv(1.0,*restart,S1,0.0,*newX);
      newstate.X = newX;
      // KX
      newKX = MVT::Clone(*restart,_blockSize);
      MVT::MvTimesMatAddMv(1.0,*Krestart,S1,0.0,*newKX);
      newstate.KX = newKX;
      // MX
      if (hasM) {
        newMX = MVT::Clone(*restart,_blockSize);
        MVT::MvTimesMatAddMv(1.0,*Mrestart,S1,0.0,*newMX);
        newstate.MX = newMX;
      }
      theta.resize(_blockSize);
      newstate.T = Teuchos::rcp( new std::vector<MagnitudeType>(theta) );
      // clear all of this memory before we initialize
      curstate = LOBPCGState<ScalarType,MV>();
      restart = Teuchos::null;
      Krestart = Teuchos::null;
      Mrestart = Teuchos::null;
      newX = Teuchos::null;
      newMX = Teuchos::null;
      newKX = Teuchos::null;
      // initialize
      lobpcg_solver->initialize(newstate);
    }
    // don't catch any other exceptions
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
        TEST_FOR_EXCEPTION(which[i] >= numlocked+_blockSize,std::logic_error,"Anasazi::LOBPCGSolverManager::solve(): indexing mistake.");
        inlocked.push_back(which[i] - _blockSize);
      }
    }

    TEST_FOR_EXCEPTION(insolver.size() + inlocked.size() != (unsigned int)sol.numVecs,std::logic_error,"Anasazi::LOBPCGSolverManager::solve(): indexing mistake.");

    // set the vecs,vals in the solution
    if (insolver.size() > 0) {
      // set vecs
      int lclnum = insolver.size();
      std::vector<int> tosol(lclnum);
      for (int i=0; i<lclnum; i++) tosol[i] = i;
      Teuchos::RefCountPtr<const MV> v = MVT::CloneView(*lobpcg_solver->getEvecs(),insolver);
      MVT::SetBlock(*v,tosol,*sol.Evecs);
      // set vals
      std::vector<MagnitudeType> fromsolver = lobpcg_solver->getEigenvalues();
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

  _problem->setSolution(sol);
  printer->stream(Debug) << "Returning " << sol.numVecs << " to eigenproblem." << endl;

  if (sol.numVecs < nev) {
    return Unconverged; // return from LOBPCGSolverManager::solve() 
  }
  return Converged; // return from LOBPCGSolverManager::solve() 
}


} // end Anasazi namespace

#endif /* ANASAZI_LOBPCG_SOLVERMANAGER_HPP */
