
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

#ifndef ANASAZI_LOBPCG_SOLMGR_HPP
#define ANASAZI_LOBPCG_SOLMGR_HPP

/*! \file AnasaziLOBPCGSolMgr.hpp
 *  \brief The Anasazi::LOBPCGSolMgr provides a powerful solver manager for the LOBPCG eigensolver.
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
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziBasicOutputManager.hpp"

/*! \class Anasazi::LOBPCGSolMgr
 *
 *  \brief The Anasazi::LOBPCGSolMgr provides a powerful and fully-featured solver manager over the LOBPCG eigensolver.
 *
 * This solver manager exists to provide a flexible manager over the Anasazi::LOBPCG eigensolver intended for general use. Features 
 * provided by this solver manager include:
 *   - locking of converged eigenpairs
 *   - global convergence on only the significant eigenpairs (instead of any eigenpairs with low residual)
 *   - recovery from Anasazi::LOBPCGRitzFailure when full orthogonalization is disabled
 *
 * These options are all driven by a list of parameters handed to the solver manager at construction. For more information, see Anasazi::LOBPCGSolMgr::LOBPCGSolMgr().

 \ingroup anasazi_solvermanagers

 \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
 */


namespace Anasazi {

template<class ScalarType, class MV, class OP>
class LOBPCGSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:

  //! @name Constructors/Destructor
  //@{ 

  /*! \brief Basic constructor for LOBPCGSolMgr.
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - \c "Which" - a \c string specifying the desired eigenvalues: SM, LM, SR or LR. Default: "SR"
   *   - \c "Block Size" - a \c int specifying the block size to be used by the underlying LOBPCG solver. Default: problem->getNEV()
   *   - \c "Full Ortho" - a \c bool specifying whether the underlying solver should employ the full orthogonalization scheme. Default: true
   *   - \c "Recover" - a \c bool specifying whether the solver manager should attempt to recover in the case of a LOBPCGRitzFailure when full orthogonalization is disabled. Default: true
   *   - \c "Maximum Iterations" - a \c int specifying the maximum number of iterations the underlying solver is allowed to perform. Default: 100
   *   - \c "Verbosity" - a sum of MsgType specifying the verbosity. Default: Anasazi::Errors
   *   - \c "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision.
   *   - \c "Relative Convergence Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding convergence. Default: true
   *   - \c "Use Locking" - a \c bool specifying whether the algorithm should employ locking of converged eigenpairs. Default: false
   *   - \c "Max Locked" - a \c int specifying the maximum number of eigenpairs to be locked. Default: problem->getNEV()
   *   - \c "Locking Quorum" - a \c int specifying the number of eigenpairs that must meet the locking criteria before locking actually occurs. Default: 1
   *   - \c "Locking Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide locking. Default: 0.1*convergence tolerance
   *   - \c "Relative Locking Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding locking. Default: true
   */
  LOBPCGSolMgr( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~LOBPCGSolMgr() {};
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
   * This method calls LOBPCG::iterate(), which will return either because a specially constructed status test evaluates to ::Passed
   * or an exception is thrown.
   *
   * A return from LOBPCG::iterate() signifies one of the following scenarios:
   *    - the maximum number of iterations has been exceeded. In this scenario, the solver manager will place\n
   *      all converged eigenpairs into the eigenproblem and return ::Unconverged.
   *    - the locking conditions have been met. In this scenario, some of the current eigenpairs will be removed\n
   *      from the eigensolver and placed into auxiliary storage. The eigensolver will be restarted with the remaining\n
   *      eigenpairs and some random information to replace the removed eigenpairs.
   *    - global convergence has been met. In this case, the most significant NEV eigenpairs in the solver and locked storage  \n
   *      have met the convergence criterion. (Here, NEV refers to the number of eigenpairs requested by the Eigenproblem.)    \n
   *      In this scenario, the solver manager will return ::Converged.
   *    - an LOBPCGRitzFailure exception has been thrown. If full orthogonalization is enabled and recovery from this exception\n
   *      is requested, the solver manager will attempt to recover from this exception by gathering the current eigenvectors,  \n
   *      preconditioned residual, and search directions from the eigensolver, orthogonormalizing the basis composed of these  \n
   *      three, projecting the eigenproblem, and restarting the eigensolver with the solution of the project eigenproblem. Any \n
   *      additional failure that occurs during this recovery effort will result in the eigensolver returning ::Unconverged.
   *
   * \returns ::ReturnType specifying:
   *    - ::Converged: the eigenproblem was solved to the specification required by the solver manager.
   *    - ::Unconverged: the eigenproblem was not solved to the specification desired by the solver manager
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
  bool _recover;
};


// Constructor
template<class ScalarType, class MV, class OP>
LOBPCGSolMgr<ScalarType,MV,OP>::LOBPCGSolMgr( 
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
  _lockQuorum(1),
  _recover(true)
{
  TEST_FOR_EXCEPTION(_problem == Teuchos::null,              std::invalid_argument, "Problem not given to solver manager.");
  TEST_FOR_EXCEPTION(!_problem->isProblemSet(),              std::invalid_argument, "Problem not set.");
  TEST_FOR_EXCEPTION(!_problem->isHermitian(),               std::invalid_argument, "Problem not symmetric.");
  TEST_FOR_EXCEPTION(_problem->getInitVec() == Teuchos::null,std::invalid_argument, "Problem does not contain initial vectors to clone from.");


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
                     "Anasazi::LOBPCGSolMgr: \"Block Size\" must be strictly positive.");

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
                     "Anasazi::LOBPCGSolMgr: \"Max Locked\" must be positive.");
  TEST_FOR_EXCEPTION(_maxLocked + _blockSize < _problem->getNEV(), 
                     std::invalid_argument,
                     "Anasazi::LOBPCGSolMgr: Not enough storage space for requested number of eigenpairs.");

  if (_useLocking) {
    _lockQuorum = pl.get("Locking Quorum",_lockQuorum);
    TEST_FOR_EXCEPTION(_lockQuorum <= 0,
                       std::invalid_argument,
                       "Anasazi::LOBPCGSolMgr: \"Locking Quorum\" must be strictly positive.");
  }

  // full orthogonalization: default true
  _fullOrtho = pl.get("Full Ortho",_fullOrtho);

  // verbosity level
  _verbosity = pl.get("Verbosity", _verbosity);

  // recover from LOBPCGRitzFailure
  _recover = pl.get("Recover",_recover);
}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType 
LOBPCGSolMgr<ScalarType,MV,OP>::solve() {

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
  // printing StatusTest
  Teuchos::RefCountPtr<StatusTestOutput<ScalarType,MV,OP> > outputtest
    = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer,combotest,1,Passed ) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Orthomanager
  Teuchos::RefCountPtr<SVQBOrthoManager<ScalarType,MV,OP> > ortho 
    = Teuchos::rcp( new SVQBOrthoManager<ScalarType,MV,OP>(_problem->getM()) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",_blockSize);
  plist.set("Full Ortho",_fullOrtho);

  // utils
  ModalSolverUtils<ScalarType,MV,OP> msutils(printer);

  //////////////////////////////////////////////////////////////////////////////////////
  // LOBPCG solver
  Teuchos::RefCountPtr<LOBPCG<ScalarType,MV,OP> > lobpcg_solver 
    = Teuchos::rcp( new LOBPCG<ScalarType,MV,OP>(_problem,sorter,printer,outputtest,ortho,plist) );
  // set any auxiliary vectors defined in the problem
  Teuchos::RefCountPtr< const MV > probauxvecs = _problem->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    lobpcg_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(probauxvecs) );
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Storage
  // 
  // lockvecs will contain eigenvectors that have been determined "locked" by the status test
  int numlocked = 0;
  Teuchos::RefCountPtr<MV> lockvecs;
  if (_useLocking) {
    lockvecs = MVT::Clone(*_problem->getInitVec(),_maxLocked);
  }
  std::vector<MagnitudeType> lockvals;
  // workMV will be used as work space for LOBPCGRitzFailure recovery and locking
  // it will be partitioned in these cases as follows:
  // for LOBPCGRitzFailure recovery:
  // workMV = [X H P OpX OpH OpP], where OpX OpH OpP will be used for K and M
  // total size: 2*3*blocksize
  // for locking
  // workMV = [X P MX MP], with MX,MP needing storage only if hasM==true
  // total size: 2*blocksize or 4*blocksize
  Teuchos::RefCountPtr<MV> workMV;
  if (_fullOrtho == false && _recover == true) {
    workMV = MVT::Clone(*_problem->getInitVec(),2*3*_blockSize);
  }
  else if (_useLocking) {
    if (_problem->getM() != Teuchos::null) {
      workMV = MVT::Clone(*_problem->getInitVec(),4*_blockSize);
    }
    else {
      workMV = MVT::Clone(*_problem->getInitVec(),2*_blockSize);
    }
  }

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
        // convtest->whichVecs() tells us which vectors from lockvecs and solver->getRitzVectors() are the ones we want
        // convtest->howMany() will tell us how many
        break;
      }
      // check locking if we didn't converge
      else if (locktest != Teuchos::null && locktest->getStatus() == Passed) {

        // remove the locked vectors,values from lobpcg_solver: put them in newvecs, newvals
        int numnew = locktest->howMany();
        TEST_FOR_EXCEPTION(numnew <= 0,std::logic_error,"Anasazi::LOBPCGSolMgr::solve(): status test mistake.");
        // get the indices
        std::vector<int> indnew = locktest->whichVecs();

        // don't lock more than _maxLocked; we didn't allocate enough space.
        if (numlocked + numnew > _maxLocked) {
          numnew = _maxLocked - numlocked;
          indnew.resize(numnew);
        }

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
          newvecs = MVT::CloneView(*lobpcg_solver->getRitzVectors(),indnew);
          // get the values
          std::vector<MagnitudeType> allvals = lobpcg_solver->getRitzValues();
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
          Teuchos::RefCountPtr<MV> newstateX, newstateMX, newstateP, newstateMP;
          //
          // workMV will be partitioned as follows: workMV = [X P MX MP], 
          //
          // make a copy of the current X,MX state
          std::vector<int> bsind(_blockSize); 
          for (int i=0; i<_blockSize; i++) bsind[i] = i;
          newstateX = MVT::CloneView(*workMV,bsind);
          MVT::SetBlock(*state.X,bsind,*newstateX);

          if (state.MX != Teuchos::null) {
            std::vector<int> block3(_blockSize);
            for (int i=0; i<_blockSize; i++) block3[i] = 2*_blockSize+i;
            newstateMX = MVT::CloneView(*workMV,block3);
            MVT::SetBlock(*state.MX,bsind,*newstateMX);
          }
          //
          // get select part, set to random, apply M
          {
            Teuchos::RefCountPtr<MV> newX = MVT::CloneView(*newstateX,indnew);
            MVT::MvRandom(*newX);

            if (newstateMX != Teuchos::null) {
              Teuchos::RefCountPtr<MV> newMX = MVT::CloneView(*newstateMX,indnew);
              OPT::Apply(*_problem->getM(),*newX,*newMX);
            }
          }

          Teuchos::Array<Teuchos::RefCountPtr<const MV> > curauxvecs = lobpcg_solver->getAuxVecs();
          Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
          // ortho X against the aux vectors
          ortho->projectAndNormalize(*newstateX,newstateMX,dummy,Teuchos::null,curauxvecs);

          if (lobpcg_solver->hasP()) {
            //
            // get P and optionally MP, orthogonalize against X and auxiliary vectors
            std::vector<int> block2(_blockSize);
            for (int i=0; i<_blockSize; i++) block2[i] = _blockSize+i;
            newstateP = MVT::CloneView(*workMV,block2);
            MVT::SetBlock(*state.P,bsind,*newstateP);

            if (state.MP != Teuchos::null) {
              std::vector<int> block4(_blockSize);
              for (int i=0; i<_blockSize; i++) block4[i] = 3*_blockSize+i;
              newstateMP = MVT::CloneView(*workMV,block4);
              MVT::SetBlock(*state.MP,bsind,*newstateMP);
            }

            if (_fullOrtho) {
              // ortho P against the new aux vectors and new X
              curauxvecs.push_back(newstateX);
              ortho->projectAndNormalize(*newstateP,newstateMP,dummy,Teuchos::null,curauxvecs);
            }
            else {
              // ortho P against the new aux vectors
              ortho->projectAndNormalize(*newstateP,newstateMP,dummy,Teuchos::null,curauxvecs);
            }
          }
          // set the new state
          LOBPCGState<ScalarType,MV> newstate;
          newstate.X  = newstateX;
          newstate.MX = newstateMX;
          newstate.P  = newstateP;
          newstate.MP = newstateMP;
          lobpcg_solver->initialize(newstate);
        }

        if (numlocked == _maxLocked) {
          // disabled locking now
          locktest->setQuorum(_blockSize+1);
        }
      }
      else {
        TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::LOBPCGSolMgr::solve(): Invalid return from lobpcg_solver::iterate().");
      }
    }
    catch (LOBPCGRitzFailure re) {
      if (_fullOrtho==true || _recover==false) {
        // if we are already using full orthogonalization, there isn't much we can do here. 
        // the most recent information in the status tests is still valid, and can be used to extract/return the 
        // eigenpairs that have converged.
        break; // while(1)
      }
      printer->stream(Errors) << "Error! Caught LOBPCGRitzFailure at iteration " << lobpcg_solver->getNumIters() << endl
                              << "Full orthogonalization is off; will try to recover." << endl;
      // otherwise, get the current "basis" from the solver, orthonormalize it, do a rayleigh-ritz, and restart with the ritz vectors
      // if there aren't enough, break and quit with what we have
      //
      // workMV = [X H P OpX OpH OpP], where OpX OpH OpP will be used for K and M
      LOBPCGState<ScalarType,MV> curstate = lobpcg_solver->getState();
      Teuchos::RefCountPtr<MV> restart, Krestart, Mrestart;
      int localsize = lobpcg_solver->hasP() ? 3*_blockSize : 2*_blockSize;
      bool hasM = _problem->getM() != Teuchos::null;
      {
        std::vector<int> recind(localsize);
        for (int i=0; i<localsize; i++) recind[i] = i;
        restart = MVT::CloneView(*workMV,recind);
      }
      {
        std::vector<int> recind(localsize);
        for (int i=0; i<localsize; i++) recind[i] = localsize+i;
        Krestart = MVT::CloneView(*workMV,recind);
      }
      if (hasM) {
        Mrestart = Krestart;
      }
      else {
        Mrestart = restart;
      }
      //
      // set restart = [X H P] and Mrestart = M*[X H P]
      //
      // put X into [0 , blockSize)
      {
        std::vector<int> blk1(_blockSize);
        for (int i=0; i < _blockSize; i++) blk1[i] = i;
        MVT::SetBlock(*curstate.X,blk1,*restart);

        // put MX into [0 , blockSize)
        if (hasM) {
          MVT::SetBlock(*curstate.MX,blk1,*Mrestart);
        }
      }
      //
      // put H into [_blockSize , 2*blockSize)
      {
        std::vector<int> blk2(_blockSize);
        for (int i=0; i < _blockSize; i++) blk2[i] = _blockSize+i;
        MVT::SetBlock(*curstate.H,blk2,*restart);

        // put MX into [0 , blockSize)
        if (hasM) {
          MVT::SetBlock(*curstate.MH,blk2,*Mrestart);
        }
      }
      // optionally, put P into [2*blockSize,3*blockSize)
      if (localsize == 3*_blockSize) {
        std::vector<int> blk3(_blockSize);
        for (int i=0; i < _blockSize; i++) blk3[i] = 2*_blockSize+i;
        MVT::SetBlock(*curstate.P,blk3,*restart);

        // put MX into [0 , blockSize)
        if (hasM) {
          MVT::SetBlock(*curstate.MP,blk3,*Mrestart);
        }
      }
      // project against auxvecs and locked vecs, and orthonormalize the basis
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
        break;
      }
      // reduce multivec size if necessary
      if (rank < localsize) {
        localsize = rank;
        std::vector<int> redind(localsize);
        for (int i=0; i<localsize; i++) redind[i] = i;
        // grab the first part of restart,Krestart
        restart = MVT::CloneView(*restart,redind);
        Krestart = MVT::CloneView(*Krestart,redind);
        if (hasM) {
          Mrestart = Krestart;
        }
        else {
          Mrestart = restart;
        }
      }
      Teuchos::SerialDenseMatrix<int,ScalarType> KK(localsize,localsize), MM(localsize,localsize), S(localsize,localsize);
      std::vector<MagnitudeType> theta(localsize);
      // project the matrices
      //
      // MM = restart^H M restart
      MVT::MvTransMv(1.0,*restart,*Mrestart,MM);
      // 
      // compute Krestart = K*restart
      OPT::Apply(*_problem->getOperator(),*restart,*Krestart);
      //
      // KK = restart^H K restart
      MVT::MvTransMv(1.0,*restart,*Krestart,KK);
      rank = localsize;
      msutils.directSolver(localsize,KK,&MM,&S,&theta,&rank,1);
      if (rank < _blockSize) {
        printer->stream(Errors) << "Error! Recovered basis of rank " << rank << " produced only " << rank << "ritz vectors.\n"
                                << "Block size is " << _blockSize << ".\n"
                                << "Recovery failed." << endl;
        break;
      }
      theta.resize(rank);
      //
      // sort the ritz values using the sort manager
      {
        Teuchos::BLAS<int,ScalarType> blas;
        std::vector<int> order(rank);
        // make a ScalarType copy for sorting
        std::vector<ScalarType> theta_st(theta.size());
        std::copy(theta.begin(),theta.end(),theta_st.begin());
        // sort
        sorter->sort( lobpcg_solver.get(), rank, &(theta_st[0]), &order );   // don't catch exception
        // put back into theta
        for (int i=0; i<rank; i++) {
          theta[i] = SCT::real(theta_st[i]);
        }
        // Sort the primitive ritz vectors
        Teuchos::SerialDenseMatrix<int,ScalarType> curS(Teuchos::View,S,rank,rank);
        msutils.permuteVectors(order,curS);
      }
      //
      Teuchos::SerialDenseMatrix<int,ScalarType> S1(Teuchos::View,S,localsize,_blockSize);
      //
      // compute the ritz vectors: store them in Krestart
      LOBPCGState<ScalarType,MV> newstate;
      Teuchos::RefCountPtr<MV> newX; 
      {
        std::vector<int> bsind(_blockSize);
        for (int i=0; i<_blockSize; i++) bsind[i] = i;
        newX = MVT::CloneView(*Krestart,bsind);
      }
      MVT::MvTimesMatAddMv(1.0,*restart,S1,0.0,*newX);
      // send X and theta into the solver
      newstate.X = newX;
      theta.resize(_blockSize);
      newstate.T = Teuchos::rcp( new std::vector<MagnitudeType>(theta) );
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
        TEST_FOR_EXCEPTION(which[i] >= numlocked+_blockSize,std::logic_error,"Anasazi::LOBPCGSolMgr::solve(): indexing mistake.");
        inlocked.push_back(which[i] - _blockSize);
      }
    }

    TEST_FOR_EXCEPTION(insolver.size() + inlocked.size() != (unsigned int)sol.numVecs,std::logic_error,"Anasazi::LOBPCGSolMgr::solve(): indexing mistake.");

    // set the vecs,vals in the solution
    if (insolver.size() > 0) {
      // set vecs
      int lclnum = insolver.size();
      std::vector<int> tosol(lclnum);
      for (int i=0; i<lclnum; i++) tosol[i] = i;
      Teuchos::RefCountPtr<const MV> v = MVT::CloneView(*lobpcg_solver->getRitzVectors(),insolver);
      MVT::SetBlock(*v,tosol,*sol.Evecs);
      // set vals
      std::vector<MagnitudeType> fromsolver = lobpcg_solver->getRitzValues();
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

    // sort the eigenvalues and permute the eigenvectors appropriately
    {
      std::vector<int> order(sol.numVecs);
      std::vector<ScalarType> vals_st(sol.numVecs);
      std::copy(sol.Evals.begin(),sol.Evals.end(),vals_st.begin());
      sorter->sort( NULL, sol.numVecs, &vals_st[0], &order );
      for (int i=0; i<sol.numVecs; i++) {
        sol.Evals[i] = SCT::real( vals_st[i] );
      }
      // now permute the eigenvectors according to order
      msutils.permuteVectors(sol.numVecs,order,*sol.Evecs);
    }
  }

  // print final summary
  lobpcg_solver->currentStatus(printer->stream(FinalSummary));

  // print timing information
  Teuchos::TimeMonitor::summarize(printer->stream(TimingDetails));

  _problem->setSolution(sol);
  printer->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << endl;

  if (sol.numVecs < nev) return Unconverged; // return from LOBPCGSolMgr::solve() 
  return Converged; // return from LOBPCGSolMgr::solve() 
}


} // end Anasazi namespace

#endif /* ANASAZI_LOBPCG_SOLMGR_HPP */
