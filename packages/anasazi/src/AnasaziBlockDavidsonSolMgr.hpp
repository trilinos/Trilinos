
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
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "Teuchos_BLAS.hpp"

/*! \class Anasazi::BlockDavidsonSolMgr
 *
 *  \brief The Anasazi::BlockDavidsonSolMgr provides a powerful and fully-featured solver manager over the BlockDavidson eigensolver.
 *
 * This solver mangaer implements a hard-locking mechanism, whereby eigenpairs designated to be locked are moved from the eigensolver and placed in
 * auxiliary storage. The eigensolver is then restarted and continues to iterate, always orthogonal to the locked eigenvectors.

 \ingroup anasazi_solvermanagers

 \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
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

  //! @name Constructors/Destructor
  //@{ 

  /*! \brief Basic constructor for BlockDavidsonSolMgr.
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - "Which" - a \c string specifying the desired eigenvalues: SM, LM, SR or LR. Default: "SR"
   *   - "Block Size" - a \c int specifying the block size to be used by the underlying block Davidson solver. Default: problem->getNEV()
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
  BlockDavidsonSolMgr( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~BlockDavidsonSolMgr() {};
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
   * This method calls BlockDavidson::iterate(), which will return either because a specially constructed status test evaluates to ::Passed
   * or an exception is thrown.
   *
   * A return from BlockDavidson::iterate() signifies one of the following scenarios:
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
  _convtol(0),
  _locktol(0),
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
  TEST_FOR_EXCEPTION(_problem == Teuchos::null,               std::invalid_argument, "Problem not given to solver manager.");
  TEST_FOR_EXCEPTION(!_problem->isProblemSet(),               std::invalid_argument, "Problem not set.");
  TEST_FOR_EXCEPTION(!_problem->isHermitian(),                std::invalid_argument, "Problem not symmetric.");
  TEST_FOR_EXCEPTION(_problem->getInitVec() == Teuchos::null, std::invalid_argument, "Problem does not contain initial vectors to clone from.");

  // which values to solve for
  _whch = pl.get("Which",_whch);
  TEST_FOR_EXCEPTION(_whch != "SM" && _whch != "LM" && _whch != "SR" && _whch != "LR",std::invalid_argument, "Invalid sorting string.");

  // convergence tolerance
  _convtol = pl.get("Convergence Tolerance",MT::prec());
  _relconvtol = pl.get("Relative Convergence Tolerance",_relconvtol);
  
  // locking tolerance
  _useLocking = pl.get("Use Locking",_useLocking);
  _rellocktol = pl.get("Relative Locking Tolerance",_rellocktol);
  _locktol = pl.get("Locking Tolerance",_convtol/10.0);

  // maximum number of restarts
  _maxRestarts = pl.get("Maximum Restarts",_maxRestarts);

  // block size: default is nev()
  _blockSize = pl.get("Block Size",_problem->getNEV());
  TEST_FOR_EXCEPTION(_blockSize <= 0, std::invalid_argument,
                     "Anasazi::BlockDavidsonSolMgr: \"Block Size\" must be strictly positive.");
  _numBlocks = pl.get("Num Blocks",2);
  TEST_FOR_EXCEPTION(_numBlocks <= 1, std::invalid_argument,
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
                     std::invalid_argument,
                     "Anasazi::BlockDavidsonSolMgr: Not enough storage space for requested number of eigenpairs.");
  TEST_FOR_EXCEPTION(_numBlocks*_blockSize + _maxLocked > MVT::GetVecLength(*_problem->getInitVec()),
                     std::invalid_argument,
                     "Anasazi::BlockDavidsonSolMgr: Potentially impossible orthogonality requests. Reduce basis size or locking size.");

  if (_useLocking) {
    _lockQuorum = pl.get("Locking Quorum",_lockQuorum);
    TEST_FOR_EXCEPTION(_lockQuorum <= 0,
                       std::invalid_argument,
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
      = Teuchos::rcp( new StatusTestOrderedResNorm<ScalarType,MV,OP>(sorter,_convtol,nev,StatusTestOrderedResNorm<ScalarType,MV,OP>::RES_ORTH,_relconvtol) );
  // locking
  Teuchos::RefCountPtr<StatusTestResNorm<ScalarType,MV,OP> > locktest;
  if (_useLocking) {
    locktest = Teuchos::rcp( new StatusTestResNorm<ScalarType,MV,OP>(_locktol,_lockQuorum,StatusTestResNorm<ScalarType,MV,OP>::RES_ORTH,_rellocktol) );
  }
  // combo class
  Teuchos::Array<Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > > alltests;
  // for an OR test, the order doesn't matter
  alltests.push_back(convtest);
  if (locktest != Teuchos::null)   alltests.push_back(locktest);
  // combo: convergence || locking 
  Teuchos::RefCountPtr<StatusTestCombo<ScalarType,MV,OP> > combotest
    = Teuchos::rcp( new StatusTestCombo<ScalarType,MV,OP>( StatusTestCombo<ScalarType,MV,OP>::OR, alltests) );
  // printing StatusTest
  Teuchos::RefCountPtr<StatusTestOutput<ScalarType,MV,OP> > outputtest
    = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer,combotest,1,Passed ) );

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

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockDavidson solver
  Teuchos::RefCountPtr<BlockDavidson<ScalarType,MV,OP> > bd_solver 
    = Teuchos::rcp( new BlockDavidson<ScalarType,MV,OP>(_problem,sorter,printer,outputtest,ortho,plist) );
  // set any auxiliary vectors defined in the problem
  Teuchos::RefCountPtr< const MV > probauxvecs = _problem->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    bd_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(probauxvecs) );
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Storage
  // for locked vectors
  int numlocked = 0;
  Teuchos::RefCountPtr<MV> lockvecs;
  // lockvecs is used to hold the locked eigenvectors, as well as for temporary storage when locking
  // when locking, we will lock some number of vectors numnew, where numnew <= maxlocked - curlocked
  // we will allocated numnew random vectors, which will go into workMV (see below)
  // we will also need numnew storage for the image of these random vectors under A and M; 
  // columns [curlocked+1,curlocked+numnew] will be used for this storage
  if (_maxLocked > 0) {
    lockvecs = MVT::Clone(*_problem->getInitVec(),_maxLocked);
  }
  std::vector<MagnitudeType> lockvals;
  // workMV will be used when restarting because of 
  // a) full basis, in which case we need 2*blockSize, for X and KX
  // b) locking, in which case we will need as many vectors as in the current basis, 
  //    a number which will be always <= (numblocks-1)*blocksize
  //    [note: this is because we never lock with curdim == numblocks*blocksize]
  Teuchos::RefCountPtr<MV> workMV;
  if (_useLocking) {
    workMV = MVT::Clone(*_problem->getInitVec(),(_numBlocks-1)*_blockSize);
  }
  else {
    workMV = MVT::Clone(*_problem->getInitVec(),2*_blockSize);
  }

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
        // we have convergence
        // convtest->whichVecs() tells us which vectors from lockvecs and solver state are the ones we want
        // convtest->howMany() will tell us how many
        break;
      }
      // check for restarting before locking: if we need to lock, it will happen after the restart
      else if ( bd_solver->getCurSubspaceDim() == bd_solver->getMaxSubspaceDim() ) {

        if ( numRestarts >= _maxRestarts ) {
          break; // break from while(1){bd_solver->iterate()}
        }
        numRestarts++;

        printer->stream(Debug) << " Performing restart number " << numRestarts << " of " << _maxRestarts << endl << endl;

        // the solver has filled its basis. 
        // the current eigenvectors will be used to restart the basis.
        std::vector<int> b1ind(_blockSize), b2ind(_blockSize);
        for (int i=0;i<_blockSize;i++) {
          b1ind[i] = i;
          b2ind[i] = _blockSize+i;
        }

        // these will be pointers into workMV
        Teuchos::RefCountPtr<MV> newV, newKV, newMV;
        { 
          newV = MVT::CloneView(*workMV,b1ind);
          MVT::SetBlock(*bd_solver->getRitzVectors(),b1ind,*newV);
        }

        if (_problem->getM() != Teuchos::null) {
          newMV = MVT::CloneView(*workMV,b2ind);
          OPT::Apply(*_problem->getM(),*newV,*newMV);
        }
        else {
          newMV = Teuchos::null;
        }

        // send this basis to the orthomanager to ensure orthonormality
        // we don't want anything in our Krylov basis that hasn't been sent through the ortho manager
        Teuchos::Array<Teuchos::RefCountPtr<const MV> > curauxvecs = bd_solver->getAuxVecs();
        Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
        ortho->projectAndNormalize(*newV,newMV,dummy,Teuchos::null,curauxvecs);
        // don't need newMV anymore, and the storage it points to will be needed for newKV
        newMV = Teuchos::null;

        // compute K*newV
        newKV = MVT::CloneView(*workMV,b2ind);
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

        // don't need either of these anymore
        newKV = Teuchos::null;
        newV = Teuchos::null;
      }
      // check locking if we didn't converge or restart
      else if (locktest != Teuchos::null && locktest->getStatus() == Passed) {

        // get number,indices of vectors to be locked
        int numnew = locktest->howMany();
        TEST_FOR_EXCEPTION(numnew <= 0,std::logic_error,"Anasazi::BlockDavidsonSolMgr::solve(): status test mistake.");
        std::vector<int> newind = locktest->whichVecs();

        // don't lock more than _maxLocked; we didn't allocate enough space.
        if (numlocked + numnew > _maxLocked) {
          numnew = _maxLocked - numlocked;
          // just use the first of them
          newind.resize(numnew);
        }

        {
          // debug printing
          printer->print(Debug,"Locking vectors: ");
          for (unsigned int i=0; i<newind.size(); i++) {printer->stream(Debug) << " " << newind[i];}
          printer->print(Debug,"\n");
        }

        BlockDavidsonState<ScalarType,MV> state = bd_solver->getState();
        //
        // get current size of basis, build index, and get a view of the current basis and projected stiffness matrix
        int curdim = state.curDim;

        // workMV will be partitioned as follows:
        // workMV = [genV augV ...]
        // genV will be of size curdim - numnew, and contain the generated basis
        // augV will be of size numnew, and contain random directions to make up for
        //      the lost space
        //
        // lockvecs will be partitioned as follows:
        // lockvecs = [curlocked augTmp ...]
        // augTmp will be used for the storage of M*augV and K*augV
        // later, the locked vectors (storeg inside the eigensolver and referenced via 
        // an MV view) will be moved into lockvecs on top of augTmp when the space is no
        // longer needed.
        Teuchos::RefCountPtr<const MV> newvecs, curlocked;
        Teuchos::RefCountPtr<MV> genV, augV, augTmp;
        Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > newKK;
        newKK = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(curdim,curdim) );
        // newvecs,newvals will contain the to-be-locked eigenpairs
        std::vector<MagnitudeType> newvals(numnew);
        {
          // setup genV
          std::vector<int> genind(curdim-numnew);
          for (int i=0; i<curdim-numnew; i++) genind[i] = i;
          genV = MVT::CloneView(*workMV,genind);
          // setup augV
          std::vector<int> augind(numnew);
          for (int i=0; i<numnew; i++) augind[i] = curdim-numnew+i;
          augV = MVT::CloneView(*workMV,augind);
          // setup curlocked
          if (numlocked > 0) {
            std::vector<int> lockind(numlocked);
            for (int i=0; i<numlocked; i++) lockind[i] = i;
            curlocked = MVT::CloneView(*lockvecs,lockind);
          }
          else {
            curlocked = Teuchos::null;
          }
          // setup augTmp
          std::vector<int> augtmpind(numnew); 
          for (int i=0; i<numnew; i++) augtmpind[i] = numlocked+i;
          augTmp = MVT::CloneView(*lockvecs,augtmpind);
          // setup newvecs
          newvecs = MVT::CloneView(*bd_solver->getRitzVectors(),newind);
          // setup newvals
          std::vector<MagnitudeType> allvals = bd_solver->getRitzValues();
          for (int i=0; i<numnew; i++) {
            newvals[i] = allvals[newind[i]];
          }
        }

        //
        // restart the solver
        //
        // we have to get the projected stiffness matrix from the solver, compute the projected eigenvectors and sort them
        // then all of the ones that don't correspond to the ones we wish to lock get orthogonalized using QR factorization
        // then we generate new, slightly smaller basis (in genV) and augment this basis with random information (in augV)
        // to preserve rank (which must be a muliple of blockSize)
        Teuchos::BLAS<int,ScalarType> blas;
        {
          const ScalarType ONE = SCT::one();
          const ScalarType ZERO = SCT::zero();

          std::vector<int> curind(curdim);
          for (int i=0; i<curdim; i++) curind[i] = i;
          Teuchos::RefCountPtr<const MV> curV = MVT::CloneView(*state.V,curind);
          Teuchos::RefCountPtr<const Teuchos::SerialDenseMatrix<int,ScalarType> > curKK;
          curKK = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*state.KK,curdim,curdim) );
          //
          // compute eigenvectors of the projected stiffness matrix
          Teuchos::SerialDenseMatrix<int,ScalarType> S(curdim,curdim);
          std::vector<MagnitudeType> theta(curdim);
          int rank = curdim;
          msutils.directSolver(curdim,*curKK,0,&S,&theta,&rank,10);
          TEST_FOR_EXCEPTION(rank != curdim,std::logic_error,"Anasazi::BlockDavidsonSolMgr::solve(): direct solve did not compute all eigenvectors."); // this should never happen
          // 
          // sort the eigenvalues (so that we can order the eigenvectors)
          {
            std::vector<int> order(curdim);
            // make a ScalarType copy of theta
            sorter->sort(bd_solver.get(),curdim,theta,&order);
            //
            // apply the same ordering to the primitive ritz vectors
            msutils.permuteVectors(order,S);
          }
          // select the non-locked eigenvectors
          std::vector<int> unlockind(curdim-numnew);
          set_difference(curind.begin(),curind.end(),newind.begin(),newind.end(),unlockind.begin());
          Teuchos::SerialDenseMatrix<int,ScalarType> Sunlocked(curdim,curdim-numnew);
          for (int i=0; i<curdim-numnew; i++) {
            blas.COPY(curdim, S[unlockind[i]], 1, Sunlocked[i], 1);
          }
          // 
          // compute the qr factorization
          {
            Teuchos::LAPACK<int,ScalarType> lapack;
            int info, lwork = (curdim*numnew)*lapack.ILAENV( 1, "geqrf", "", curdim, curdim-numnew );
            std::vector<ScalarType> tau(curdim-numnew), work(lwork);
            lapack.GEQRF(curdim,curdim-numnew,Sunlocked.values(),Sunlocked.stride(),&tau[0],&work[0],lwork,&info);
            TEST_FOR_EXCEPTION(info != 0,std::logic_error,"Anasazi::BlockDavidsonSolMgr::solve(): Error calling GEQRF in locking code.");
            lapack.UNGQR(curdim,curdim-numnew,curdim-numnew,Sunlocked.values(),Sunlocked.stride(),&tau[0],&work[0],lwork,&info);
            TEST_FOR_EXCEPTION(info != 0,std::logic_error,"Anasazi::BlockDavidsonSolMgr::solve(): Error calling UNGQR in locking code.");

            if (printer->isVerbosity(Debug)) {
              Teuchos::SerialDenseMatrix<int,ScalarType> StS(curdim-numnew,curdim-numnew);
              int info = StS.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ONE,Sunlocked,Sunlocked,ZERO);
              TEST_FOR_EXCEPTION(info != 0, std::logic_error, "Anasazi::BlockDavidsonSolMgr::solve(): Input error to SerialDenseMatrix::multiply.");
              for (int i=0; i<curdim-numnew; i++) {
                StS(i,i) -= ONE;
              }
              printer->stream(Debug) << "Locking: Error in Snew^T Snew == I : " << StS.normFrobenius() << endl;
            }
          }

          // compute the first part of the new basis: genV = curV*Sunlocked
          MVT::MvTimesMatAddMv(ONE,*curV,Sunlocked,ZERO,*genV);
          //
          // compute leading block of the new projected stiffness matrix
          {
            Teuchos::SerialDenseMatrix<int,ScalarType> tmpKK(curdim,curdim-numnew),
                                                       newKK11(Teuchos::View,*newKK,curdim-numnew,curdim-numnew),
                                                       curKKsym(*curKK);
            for (int j=0; j<curdim; j++) {
              for (int i=j+1; i<curdim; i++) {
                curKKsym(i,j) = curKKsym(j,i);
              }
            }
            int info = tmpKK.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,curKKsym,Sunlocked,ZERO);
            TEST_FOR_EXCEPTION(info != 0, std::logic_error, "Anasazi::BlockDavidsonSolMgr::solve(): Input error to SerialDenseMatrix::multiply.");
            info = newKK11.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ONE,Sunlocked,tmpKK,ZERO);
            TEST_FOR_EXCEPTION(info != 0, std::logic_error, "Anasazi::BlockDavidsonSolMgr::solve(): Input error to SerialDenseMatrix::multiply.");
          }
          //
          // generate random data to fill the rest of the basis back to curdim
          MVT::MvRandom(*augV);
          // 
          // orthogonalize it against auxvecs, the current basis, and all locked vectors
          {
            Teuchos::Array<Teuchos::RefCountPtr<const MV> > against;
            Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
            if (probauxvecs != Teuchos::null) against.push_back(probauxvecs);
            if (curlocked != Teuchos::null)   against.push_back(curlocked);
            against.push_back(newvecs);
            against.push_back(genV);
            ortho->projectAndNormalize(*augV,Teuchos::null,dummy,Teuchos::null,against);
          }
          //
          // project the stiffness matrix on the new part
          {
            OPT::Apply(*_problem->getOperator(),*augV,*augTmp);
            Teuchos::SerialDenseMatrix<int,ScalarType> KK12(Teuchos::View,*newKK,curdim-numnew,numnew,0,curdim-numnew),
                                                       KK22(Teuchos::View,*newKK,numnew,numnew,curdim-numnew,curdim-numnew);
            MVT::MvTransMv(ONE,*genV,*augTmp,KK12);
            MVT::MvTransMv(ONE,*augV,*augTmp,KK22);
          }
        }

        // done with pointers
        augV = genV = augTmp = Teuchos::null;
        curlocked = Teuchos::null;

        // put newvecs into lockvecs, newvals into lockvals
        {
          lockvals.insert(lockvals.end(),newvals.begin(),newvals.end());

          std::vector<int> indlock(numnew);
          for (int i=0; i<numnew; i++) indlock[i] = numlocked+i;
          MVT::SetBlock(*newvecs,indlock,*lockvecs);
          newvecs = Teuchos::null;

          numlocked += numnew;
          std::vector<int> curind(numlocked);
          for (int i=0; i<numlocked; i++) curind[i] = i;
          curlocked = MVT::CloneView(*lockvecs,curind);
        }
        // add locked vecs as aux vecs, along with aux vecs from problem
        // add lockvals to convtest
        {
          convtest->setAuxVals(lockvals);

          Teuchos::Array< Teuchos::RefCountPtr<const MV> > aux;
          if (probauxvecs != Teuchos::null) aux.push_back(probauxvecs);
          aux.push_back(curlocked);
          bd_solver->setAuxVecs(aux);
        }

        // get pointer to new basis, KK
        {
          std::vector<int> curind(curdim);
          for (int i=0; i<curdim; i++) curind[i] = i;
          state.V = MVT::CloneView(*workMV,curind);
          state.KK = newKK;
          state.curDim = curdim;
          // clear the rest
          state.X = Teuchos::null;
          state.KX = Teuchos::null;
          state.MX = Teuchos::null;
          state.R = Teuchos::null;
          state.H = Teuchos::null;
          state.T = Teuchos::null;
          //
          // pass new state to the solver
          bd_solver->initialize(state);
        }

        if (numlocked == _maxLocked) {
          // disabled locking now by setting quorum to unreachable number
          locktest->setQuorum(_blockSize+1);
        }
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
      Teuchos::RefCountPtr<const MV> v = MVT::CloneView(*bd_solver->getRitzVectors(),insolver);
      MVT::SetBlock(*v,tosol,*sol.Evecs);
      // set vals
      std::vector<MagnitudeType> fromsolver = bd_solver->getRitzValues();
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
      sorter->sort(bd_solver.get(), sol.numVecs, sol.Evals, &order );
      // now permute the eigenvectors according to order
      msutils.permuteVectors(sol.numVecs,order,*sol.Evecs);
    }
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
