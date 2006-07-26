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

/*! \file AnasaziBlockDavidson.hpp
  \brief Implementation of the block Davidson method
*/

#ifndef ANASAZI_BLOCK_DAVIDSON_HPP
#define ANASAZI_BLOCK_DAVIDSON_HPP

#include "AnasaziEigensolver.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziModalSolverUtils.hpp"
#include "AnasaziSortManager.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!     \class Anasazi::BlockDavidson

        \brief This class implements the block Davidson method, an iterative
        method for solving symmetric eigenvalue problems.

        \author Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {
  
  template <class ScalarType, class MV, class OP>
  class BlockDavidson : public Eigensolver<ScalarType,MV,OP> { 
  public:
    //@{ \name Constructor/Destructor.
    
    //! %Anasazi::BlockDavidson constructor.
    BlockDavidson( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
                   const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sm,
                   const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
                   Teuchos::ParameterList &pl
                   );
    
    //! %Anasazi::BlockDavidson destructor.
    virtual ~BlockDavidson() {};
    //@}
    
    //@{ \name Solver application methods.
    
    /*! \brief This method uses iterate to compute approximate solutions to the
      original problem.  It may return without converging if it has taken the
      maximum number of iterations or numerical breakdown is observed.
    */
    ReturnType solve();
    //@}
    
    //@{ \name Solver status methods.
    
    //! Get the current iteration count.
    int GetNumIters() const { return(_iter); };
    
    //! Get the current restart count of the iteration method.
    /*! Some eigensolvers can perform restarts (i.e.Arnoldi) to reduce memory
      and orthogonalization costs.  For other eigensolvers, like LOBPCG or block Davidson,
      this is not a valid stopping criteria.
    */
    int GetNumRestarts() const { return(_numRestarts); };

    //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
    int GetBlockSize() const { return(_blockSize); }
    
    /*! \brief Get a constant reference to the current linear problem, 
      which may include a current solution.
    */
    Eigenproblem<ScalarType,MV,OP>& GetEigenproblem() const { return(*_problem); };
    
    //@}
    
    //@{ \name Output methods.
    
    //! This method requests that the solver print out its current status to screen.
    void currentStatus();

    //@}
  private:
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    typedef typename std::vector<ScalarType>::iterator STiter;
    typedef typename std::vector<MagnitudeType>::iterator MTiter;

    /*! \brief These methods will not be defined.
     */
    BlockDavidson(const BlockDavidson<ScalarType,MV,OP> &method);
    BlockDavidson<ScalarType,MV,OP>& operator=(const BlockDavidson<ScalarType,MV,OP> &method);
    //
    // Internal methods
    //
    void accuracyCheck(const MV *X, const MV *MX, const MV *Q) const; 
    //
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem; 
    Teuchos::RefCountPtr<OutputManager<ScalarType> > _om; 
    Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > _sm; 
    Teuchos::ParameterList _pl;
    //
    // Information obtained from the eigenproblem
    //
    Teuchos::RefCountPtr<OP> _Op;
    Teuchos::RefCountPtr<OP> _MOp;
    Teuchos::RefCountPtr<OP> _Prec;
    Teuchos::RefCountPtr<MV> _evecs;
    Teuchos::RefCountPtr<std::vector<ScalarType> > _evals;
    const int _nev;  
    //
    // Internal data.
    //
    const int _maxIter, _blockSize;
    const MagnitudeType _residual_tolerance;
    int _numBlocks, _numRestarts, _iter, _dimSearch, _knownEV;
    //
    // Internal utilities class required by eigensolver.
    //
    ModalSolverUtils<ScalarType,MV,OP> _MSUtils;
    std::vector<MagnitudeType> _theta, _normR, _resids;
    //
    // Output stream from the output manager
    //
    std::ostream& _os;
    //
    // Internal timers
    //
    bool _restartTimers;
    Teuchos::RefCountPtr<Teuchos::Time> _timerOp, _timerMOp, _timerPrec,
                                        _timerSortEval, _timerDS,
                                        _timerOrtho, _timerTotal;
    //
    // Counters
    //
    int _count_ApplyOp, _count_ApplyM, _count_ApplyPrec;

  };
  //
  // Implementation
  //
  template <class ScalarType, class MV, class OP>
  BlockDavidson<ScalarType,MV,OP>::BlockDavidson(const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
                                                 const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sm,
                                                 const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
                                                 Teuchos::ParameterList &pl
                                                 ):
    _problem(problem), 
    _om(om),
    _sm(sm),
    _pl(pl),
    _Op(_problem->GetOperator()),
    _MOp(_problem->GetM()),
    _Prec(_problem->GetPrec()),
    _evecs(_problem->GetEvecs()), 
    _evals(problem->GetEvals()), 
    _nev(problem->GetNEV()), 
    _maxIter(_pl.get("Max Iters", 300)),
    _blockSize(_pl.get("Block Size", 1)),
    _residual_tolerance(_pl.get("Tol", (MagnitudeType)(1.0e-6) )),
    _numBlocks(_pl.get("Max Blocks", 25)), 
    _numRestarts(0), 
    _iter(0), 
    _dimSearch(_blockSize*_numBlocks),    
    _knownEV(0),
    _MSUtils(om),
    _os(_om->GetOStream()),
    _restartTimers(_pl.get("Restart Timers", false)),
    _timerOp(Teuchos::TimeMonitor::getNewTimer("Operation Op*x")),
    _timerMOp(Teuchos::TimeMonitor::getNewTimer("Operation M*x")),
    _timerPrec(Teuchos::TimeMonitor::getNewTimer("Operation Prec*x")),
    _timerSortEval(Teuchos::TimeMonitor::getNewTimer("Sorting eigenvalues")),
    _timerDS(Teuchos::TimeMonitor::getNewTimer("Direct solve")),
    _timerOrtho(Teuchos::TimeMonitor::getNewTimer("Orthogonalization")),
    _timerTotal(Teuchos::TimeMonitor::getNewTimer("Total time")),
    _count_ApplyOp(0),
    _count_ApplyM(0),
    _count_ApplyPrec(0)
  {     
    //
    // Define dense local matrices and arrays
    //
    // theta = Storage for local eigenvalues     (size: dimSearch)
    // normR = Storage for the norm of residuals (size: blockSize)
    // resids = Storage for the norm of the computed eigenvalues
    //
    _theta.resize( _dimSearch );
    _normR.resize( _blockSize );
    _resids.resize( _nev );
    //    
  }

  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::currentStatus() 
  {
    int i;
    if (_om->doPrint()) {
      _os.setf(ios::scientific, ios::floatfield);
      _os.precision(6);
      _os <<endl;
      _os <<"******************* CURRENT STATUS *******************"<<endl;
      _os <<"The number of iterations performed is " <<_iter<<endl;
      _os <<"The number of restarts performed is "<<_numRestarts<<endl;
      _os <<"The block size is "<<_blockSize<<endl;
      _os <<"The number of eigenvalues requested is "<<_nev<<endl;
      _os <<"The number of eigenvalues computed is "<<_knownEV<<endl;
      _os <<"The requested residual tolerance is "<<_residual_tolerance<<endl;
      _os <<"The number of operations Op*x   is "<<_count_ApplyOp<<endl;
      _os <<"The number of operations M*x    is "<<_count_ApplyM<<endl;
      _os <<"The number of operations Prec*x is "<<_count_ApplyPrec<<endl;
      _os <<endl;
      _os <<"COMPUTED EIGENVALUES                 "<<endl;
      _os << std::setw(16) << std::right << "Eigenvalue" 
          << std::setw(16) << std::right << "Ritz Residual"
          << endl;
      _os <<"------------------------------------------------------"<<endl;
      if ( _knownEV > 0 ) {
        for (i=0; i<_knownEV; i++) {
          _os << std::setw(16) << std::right << (*_evals)[i]
              << std::setw(16) << std::right << _resids[i]
              << endl;
        }
      } 
      else {
        _os <<"[none computed]"<<endl;
      }
      _os <<endl;
      _os <<"CURRENT EIGENVALUE ESTIMATES             "<<endl;
      _os << std::setw(16) << std::right << "Ritz value" 
          << std::setw(16) << std::right << "Residual"
          << endl;
      _os <<"------------------------------------------------------"<<endl;
      if ( _iter > 0 ) {
        for (i=0; i<_blockSize; i++) {
          _os << std::setw(16) << std::right << _theta[i] 
              << std::setw(16) << std::right << _normR[i]
              << endl;
        }
      } 
      else {
        _os <<"[none computed]" << endl;
      }
      _os << "******************************************************" << endl;  
      _os << endl; 
    }
  }
  
  template <class ScalarType, class MV, class OP>
  ReturnType BlockDavidson<ScalarType,MV,OP>::solve () 
  {
    Teuchos::TimeMonitor LocalTimer(*_timerTotal,_restartTimers);

    if ( _restartTimers ) {
      _timerOp->reset();
      _timerMOp->reset();
      _timerPrec->reset();
      _timerSortEval->reset();
      _timerDS->reset();
      _timerOrtho->reset();
      _count_ApplyOp = 0;
      _count_ApplyM = 0;
      _count_ApplyPrec = 0;
    }

    //
    // Check the Anasazi::Eigenproblem was set by user, if not, return failed.
    //
    if ( !_problem->IsProblemSet() ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error ))
        _os << "ERROR : Anasazi::Eigenproblem was not set, call Anasazi::Eigenproblem::SetProblem() before calling solve"<< endl;
      return Failed;
    }
    //
    // Check the Anasazi::Eigenproblem is symmetric, if not, return failed.
    //
    if ( !_problem->IsSymmetric() ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error ))
        _os << "ERROR : Anasazi::Eigenproblem is not symmetric" << endl;
      return Failed;
    }
    //
    // Retrieve the initial vector and operator information from the Anasazi::Eigenproblem.
    //
    Teuchos::RefCountPtr<MV> iVec = _problem->GetInitVec();
    //
    if ( iVec.get() == 0 ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : Initial vector is not specified, set initial vector in eigenproblem "<<endl;
      return Failed;
    }

    int dim = MVT::GetVecLength( *iVec );
    //
    // Check that the maximum number of blocks for the eigensolver is a positive number
    //    
    if ( _numBlocks<=0 ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : numBlocks = "<< _numBlocks <<" [ should be positive number ] " << endl;
      return Failed;
    } 
    //
    // Check that the maximum number of iterations is a positive number
    //    
    if ( _maxIter<=0 ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : maxIter = "<< _maxIter <<" [ should be positive number ] " << endl;
      return Failed;
    } 
    //
    // Check that the search subspace is larger than the number of eigenvalues requested
    //
    if ( _numBlocks*_blockSize < _nev ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : Search space dimension (numBlocks*blockSize) = "<< _numBlocks*_blockSize 
            << " [ should be greater than "<< _nev << " ] " << endl;
      return Failed;
    } 
    //
    // If the search subspace dimension is the same size as the number of requested eigenvalues,
    // then we must be computing all of them.
    //
    if ( (_numBlocks*_blockSize == _nev) && (_nev != dim) ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error ))
        _os << "ERROR : Search space dimension (numBlocks*blockSize) = "<< _numBlocks*_blockSize 
            << " [ should be greater than "<< _nev << " ] " << endl;
      return Failed;
    }
    //
    // If the search subspace dimension is larger than the dimension of the operator, reset
    // the maximum number of blocks accordingly.
    //    
    if (_numBlocks*_blockSize > dim ) {
      if (_om->isVerbosityAndPrint( Anasazi::Warning ))
        _os << "WARNING : Search space dimension (numBlocks*blockSize) = "<< _numBlocks*_blockSize 
            <<" [ should not be greater than " << dim << " ] " << endl;
      
      // Set the maximum number of blocks in the factorization below the dimension of the space.
      _numBlocks = dim / _blockSize;
      
      if (_om->isVerbosityAndPrint( Anasazi::Warning ))
        _os << "WARNING : numBlocks reset to "<< _numBlocks << endl;
    }
    //
    // Reinitialize internal data and pointers, prepare for solve
    //
    _numRestarts = 0; 
    _iter = 0; 
    _knownEV = 0;
    //
    // Necessary variables
    //
    int i, j;
    int info, nb, lwork;
    int bStart = 0, offSet = 0;
    int nFound = _blockSize;
    bool reStart = false;
    bool criticalExit = false;
    ReturnType ret;
    Teuchos::RefCountPtr<MV> Xcurrent, Xprev, Xtotal, Xnext;
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::LAPACK<int,ScalarType> lapack;
    std::vector<int> _order;
    //
    // Define local block vectors
    //
    // MX = Working vectors (storing M*X if M is specified, else pointing to X)
    // KX = Working vectors (storing K*X)
    //
    // R = Residuals
    //
    Teuchos::RefCountPtr<MV> X, MX, KX, R;
    X = MVT::Clone( *iVec, _dimSearch + _blockSize );
    KX = MVT::Clone( *iVec, _blockSize );
    R = MVT::Clone( *iVec, _blockSize );
    //
    // Check to see if there is a mass matrix, so we know how much space is required.
    // [ If there isn't a mass matrix we can use a view of X.]
    //
    if (_MOp.get()) {
      MX = MVT::Clone( *iVec, _blockSize );
    }
    else {
      std::vector<int> index( _blockSize );
      for (int i=0; i<_blockSize; i++) {
        index[i] = i;
      }
      MX = MVT::CloneView( *X, index );
    }
    //
    // KK = Local stiffness matrix               (size: dimSearch x dimSearch)
    //
    // S = Local eigenvectors                    (size: dimSearch x dimSearch)
    //
    Teuchos::SerialDenseMatrix<int,ScalarType> KK( _dimSearch, _dimSearch ), S( _dimSearch, _dimSearch );
    // work vector for GEQRF and ORGQR
    std::vector<ScalarType> tau( _dimSearch );

    //
    // Initialize the workspace.
    //
    {
      // view vector into the basis
      Teuchos::RefCountPtr<MV> tmpXinit;
      // index vector
      std::vector<int> index;

      //
      // Determine how many init vectors were specified by the user
      //
      int numIVecs = MVT::GetNumberVecs( *iVec );
      if (numIVecs > _blockSize) {
        numIVecs = _blockSize;
      }

      //
      // Get a view into the basis
      //
      index.resize( numIVecs );
      for (i=0; i<numIVecs; i++) {
        index[i] = i;
      }
      tmpXinit = MVT::CloneView( *X, index );

      //
      // Copy the first numIVecs of the initial vectors into the first
      // numIVecs vectors of the basis (any additional vectors in iVec are ignored)
      //
      MVT::SetBlock( *iVec, index, *tmpXinit );

      //
      // Augment the initial vectors with random vectors if necessary
      //
      int leftOver = _blockSize - numIVecs;
      if (leftOver > 0) {
        index.resize(leftOver);
        for (i=0; i<leftOver; i++) {
          index[i] = numIVecs + i;
        }
        tmpXinit = MVT::CloneView( *X, index );
        MVT::MvRandom( *tmpXinit );
      }
    }

    //
    // Determine the maximum number of blocks for this factorization.
    // ( NOTE:  This will be the _numBlocks since we don't know about already converged vectors here )
    int maxBlock = _numBlocks;
    //
    // Compute the required workspace for this factorization.
    //
    lwork = lapack.ILAENV( 1, "geqrf", "", maxBlock*_blockSize, maxBlock*_blockSize );
    lwork *= _blockSize;
    std::vector<ScalarType> work( lwork );
    //
    while ( _iter <= _maxIter ) {
       
      for (nb = bStart; nb < maxBlock ; nb++) 
      {
        // 
        // Increment iteration counter
        //
        _iter++;
        if ( _iter > _maxIter )
          break;
        //
        // Get a view of the current vectors being worked on.
        // The offset is dependent on the number of iterations been taken so far (nb)
        // and the number of known eigenvalues (_knownEV).
        //
        int localSize = nb*_blockSize;
        std::vector<int> index( _blockSize );
        for (i=0; i < _blockSize; i++)
          index[i] = localSize + _knownEV + i;
        //
        Xcurrent = MVT::CloneView( *X, index );
        //
        if (_knownEV + localSize > 0) {
          index.resize( _knownEV + localSize );
          for (i=0; i < _knownEV + localSize; i++)
            index[i] = i;
          
          Xprev = MVT::CloneView( *X, index );
        }
        //
        // Apply the mass matrix.
        //        
        if (_MOp.get()) {
          _timerMOp->start();
          ret = OPT::Apply( *_MOp, *Xcurrent, *MX );
          _timerMOp->stop();
          _count_ApplyM += MVT::GetNumberVecs( *Xcurrent );
          if (ret != Ok) {
            if (_om->isVerbosityAndPrint(Error)) {
              _os << "ERROR : Applying M operator in BlockReduction" << endl;
            }
            criticalExit = true;
            break; // break out of for(nb) loop
          }
        }
        //
        // Orthonormalize Xcurrent against the known eigenvectors and previous vectors.
        //
        if (nb == bStart) {
          if (nFound > 0) {
            if (_knownEV == 0) {
              _timerOrtho->start();
              info = _MSUtils.massOrthonormalize( *Xcurrent, *MX, _MOp.get(), *Xcurrent, nFound, 2 );
              _timerOrtho->stop();
            }
            else {
              _timerOrtho->start();
              info = _MSUtils.massOrthonormalize( *Xcurrent, *MX, _MOp.get(), *Xprev,    nFound, 0 );
              _timerOrtho->stop();
            }
          }
          nFound = 0;
        } 
        else {
          _timerOrtho->start();
          info = _MSUtils.massOrthonormalize( *Xcurrent, *MX, _MOp.get(), *Xprev, _blockSize, 0 );
          _timerOrtho->stop();
        }
        //
        // Exit the code if there has been a problem.
        //
        if (info < 0) 
          return Failed;
        //
        // Check orthogonality of X ( if required )
        //
        if (_om->isVerbosity( OrthoDetails ) ) {
          if (localSize > 0) {
            accuracyCheck( Xcurrent.get(), MX.get(), Xprev.get() );
          }
          else {
            accuracyCheck( Xcurrent.get(), MX.get(), NULL );
          }
        }
        //
        // Apply the stiffness matrix.
        //
        _timerOp->start();
        ret = OPT::Apply( *_Op, *Xcurrent, *KX );
        _timerOp->stop();
        _count_ApplyOp += MVT::GetNumberVecs( *Xcurrent );
        if (ret != Ok) {
          if (_om->isVerbosityAndPrint(Error)) {
            _os << "ERROR : Applying Op operator in BlockReduction" << endl;
          }
          criticalExit = true;
          break; // break out of for(nb) loop
        }
        //
        // Update the local stiffness matrix ( Xtotal^T * K * Xcurrent where Xtotal = [Xprev Xcurrent] )
        // Note:  Only the upper half of the matrix is stored in KK
        //
        index.resize( localSize + _blockSize );
        for (i=0; i < localSize + _blockSize; i++) {
          index[i] = _knownEV + i;
        }
        Xtotal = MVT::CloneView( *X, index );
        Teuchos::SerialDenseMatrix<int,ScalarType> subKK( Teuchos::View, KK, localSize+_blockSize, _blockSize, 0, localSize );
        MVT::MvTransMv( one, *Xtotal, *KX, subKK );
        //
        // Perform spectral decomposition
        //
        int nevLocal = localSize+_blockSize;
        _timerDS->start();
        info = _MSUtils.directSolver(localSize+_blockSize, KK, 0, &S, &_theta, &nevLocal, 10);
        _timerDS->stop();
        //
        // Exit the code if there has been a problem.
        //
        if (info != 0) {
          //
          // Detect critical failure
          // 
          if (info<0) {
            criticalExit = true;
            break;
          }
          //
          // Restart:  Spectral decomposition failed.
          // Reinitialize all counters and randomize the starting block
          //
          reStart = true;
          _numRestarts++;
          index.resize( _blockSize );
          for (i=0; i<_blockSize; i++) {
            index[i] = _knownEV + i;
          }
          Teuchos::RefCountPtr<MV> Xinit = MVT::CloneView( *X, index );
          MVT::MvRandom( *Xinit );
          nFound = _blockSize;
          bStart = 0;
          break;
        }
        //
        //---------------------------------------------------
        // Sort the ritz values using the sort manager
        //---------------------------------------------------
        // The sort manager is templated on ScalarType
        // Make a ScalarType copy of _theta for sorting
        std::vector<ScalarType> _theta_st(_theta.size());
        _timerSortEval->start();
        std::copy<MTiter,STiter>(_theta.begin(),_theta.begin()+localSize+_blockSize,_theta_st.begin());
        _order.resize(localSize+_blockSize);
        ret = _sm->sort( this, localSize+_blockSize, &(_theta_st[0]), &_order );
        //reorder _theta according to sorting results from _theta_st
        std::vector<MagnitudeType> _theta_copy(_theta);
        for (i=0; i<localSize+_blockSize; i++) {
          _theta[i] = _theta_copy[_order[i]];
        }
        _timerSortEval->stop();
        if (ret != Ok) {
          if (_om->isVerbosityAndPrint(Error)) {
            _os << "ERROR : Sorting in solve()" << endl;
          }
          criticalExit = true;
          break; // break out of for(nb) loop
        }
        // Sort the primitive ritz vectors
        // We need the first _blockSize vectors ordered to generate the next
        // columns immediately below, as well as when/if we restart.
        Teuchos::SerialDenseMatrix<int,ScalarType> copyS( S );
        for (i=0; i<localSize+_blockSize; i++) {
          blas.COPY(localSize+_blockSize, copyS[_order[i]], 1, S[i], 1);
        }
        // Create a view matrix of the first _blockSize vectors
        Teuchos::SerialDenseMatrix<int,ScalarType> subS( Teuchos::View, S, localSize+_blockSize, _blockSize );
        //
        //---------------------------------------------------
        // Update the search space 
        //---------------------------------------------------
        // KX = Xtotal * S where S is the eigenvectors of the projected problem.
        //
        MVT::MvTimesMatAddMv( one, *Xtotal, subS, zero, *KX );
        //
        // Apply the mass matrix for the next block
        // 
        if (_MOp.get()) {
          _timerMOp->start();
          ret = OPT::Apply( *_MOp, *KX, *MX );
          _timerMOp->stop();
          _count_ApplyM += MVT::GetNumberVecs( *KX );
          if (ret != Ok) {
            if (_om->isVerbosityAndPrint(Error)) {
              _os << "ERROR : Applying M operator in BlockReduction" << endl;
            }
            criticalExit = true;
            break; // break out of for(nb) loop
          }
        }
        //
        // Apply the stiffness matrix for the next block
        //
        _timerOp->start();
        ret = OPT::Apply( *_Op, *KX, *R );
        _timerOp->stop();
        _count_ApplyOp += MVT::GetNumberVecs( *KX );
        if (ret != Ok) {
          if (_om->isVerbosityAndPrint(Error)) {
            _os << "ERROR : Applying Op operator in BlockReduction" << endl;
          }
          criticalExit = true;
          break; // break out of for(nb) loop
        }
        //
        // Compute the residual :
        // R = KX - MX*diag(theta)
        // 
        Teuchos::SerialDenseMatrix<int,ScalarType> D(_blockSize, _blockSize);
        for (i=0; i<_blockSize; i++ ) {
          D(i,i) = -_theta[i];
        }
        //
        if (_MOp.get()) {
          MVT::MvTimesMatAddMv( one, *MX, D, one, *R );
        }
        else {
          MVT::MvTimesMatAddMv( one, *KX, D, one, *R );
        }
        MVT::MvNorm( *R, &_normR );
        //
        // Scale the norms of residuals with the eigenvalues and check for converged eigenvectors.
        //
        nFound = 0;
        for (j=0; j<_blockSize; j++) {
          // Scale eigenvalues if _theta is non-zero.
          if ( _theta[j] != zero ) {
            _normR[j] /= SCT::magnitude(_theta[j]);
          }
          // Check for convergence
          if (_normR[j] < _residual_tolerance) {
            nFound ++;          
          }
        }
        // Print information on current iteration
        if (_om->isVerbosityAndPrint( IterationDetails )) {
          _os << " Iteration " << _iter << " - Number of converged eigenvectors ";
          _os << _knownEV + nFound << endl;
        } 
        
        if (_om->isVerbosityAndPrint( IterationDetails )) {
          _os << endl;
          cout.precision(2);
          cout.setf(ios::scientific, ios::floatfield);
          for (i=0; i<_blockSize; ++i) {
            _os << " Iteration " << _iter << " - Scaled Norm of Residual " << i;
            _os << " = " << _normR[i] << endl;
          }
          _os << endl;
          cout.precision(2);
          for (i=0; i<localSize + _blockSize; ++i) {
            _os << " Iteration "<< _iter << " - Ritz eigenvalue " << i;
            cout.setf((fabs(_theta[i]) < 0.01) ? ios::scientific : ios::fixed, ios::floatfield);
            _os << " = " << _theta[i] << endl;
          }
          _os << endl;
        }
        //
        // Exit the loop to treat the converged eigenvectors
        //
        if (nFound > 0) {
          nb += 1;
          offSet = 0;
          break;
        }
        //
        // Apply the preconditioner on the residuals
        //
        if (maxBlock == 1) {
          if (_Prec.get()) {
            _timerPrec->start();
            ret = OPT::Apply( *_Prec, *R, *Xcurrent );
            _timerPrec->stop();
            _count_ApplyPrec += MVT::GetNumberVecs( *R );
            if (ret != Ok) {
              if (_om->isVerbosityAndPrint(Error)) {
                _os << "ERROR : Applying Prec operator in BlockReduction" << endl;
              }
              criticalExit = true;
              break; // break out of for(nb) loop
            }
          }
          else
            MVT::MvAddMv( one, *R, zero, *R, *Xcurrent );
          //
          // Update the preconditioned residual 
          //
          MVT::MvAddMv( one, *KX, -one, *Xcurrent, *Xcurrent );
          break;
        } // if (maxBlock == 1)

        if (nb == maxBlock - 1) {
          nb++;
          break;
        }
        //
        // Prepare next block in factorization
        //
        index.resize( _blockSize );
        for( i=0; i<_blockSize; i++) {
          index[i] = _knownEV + localSize + _blockSize + i;
        }
        Xnext = MVT::CloneView( *X, index );
        if (_Prec.get()) {
          _timerPrec->start();
          ret = OPT::Apply( *_Prec, *R, *Xnext );
          _timerPrec->stop();
          _count_ApplyPrec += MVT::GetNumberVecs( *R );
          if (ret != Ok) {
            if (_om->isVerbosityAndPrint(Error)) {
              _os << "ERROR : Applying Prec operator in BlockReduction" << endl;
            }
            criticalExit = true;
            break; // break out of for(nb) loop
          }
        }
        else {
          MVT::MvAddMv( one, *R, zero, *R, *Xnext );
        }
      } // for (nb = bStart; nb < maxBlock; nb++)

      //
      // Check if there is any reason we need to skip the rest of the while loop.
      //
      if (_iter > _maxIter) {
        break;
      }

      // Experienced unrecoverable error
      if ( criticalExit ) {
        break;
      }

      if (reStart == true) {
        reStart = false;
        continue;
      }

      //
      // Store the final converged eigenvectors
      //
      if (_knownEV + nFound >= _nev) {
        std::vector<int> index(1);
        for (j=0; j<_blockSize; j++) {
            if (_normR[j] < _residual_tolerance) {
              index[0] = j;
              Teuchos::RefCountPtr<MV> tmpKX = MVT::CloneView( *KX, index );
              index[0] = _knownEV;
              MVT::SetBlock( *tmpKX, index, *_evecs );
              (*_evals)[_knownEV] = _theta[j];
              _resids[_knownEV] = _normR[j];
              _knownEV++;
            }
            if (_knownEV == _nev) {
              break;
            }
        }
        break;
      } // if (_knownEV + nFound >= _nev)      
      //
      // Store the converged eigenvectors and define the restarting block
      // in the particular case of 1 block ( maxBlock == 1 ).
      //
      if (maxBlock == 1) {
        if (nFound > 0) {
          std::vector<int> index(1);
          int tmp_ptr = _knownEV + nFound;
          nFound = 0;
          for (j=0; j<_blockSize; j++) {
            //
            // Get a view of the current prospective eigenvector.
            //
            index[0] = j;
            Teuchos::RefCountPtr<MV> tmpKX = MVT::CloneView( *KX, index );
            if (_normR[j] < _residual_tolerance) {
              index[0] = _knownEV;
              MVT::SetBlock( *tmpKX, index, *X );              
              MVT::SetBlock( *tmpKX, index, *_evecs );              
              (*_evals)[_knownEV] = _theta[j];
              _resids[_knownEV] = _normR[j];
              _knownEV++;
              nFound++;              
            }
            else {
              index[0] = tmp_ptr + (j-nFound);
              MVT::SetBlock( *tmpKX, index, *X );
            }
          } // for (j=0; j<_blockSize; j++)
          //
          index.resize( nFound );
          for (i=0; i<nFound; i++) {
            index[i] = _knownEV + _blockSize - nFound + i;
          }
          Xnext = MVT::CloneView( *X, index );
          MVT::MvRandom( *Xnext );
        }
        else {
          nFound = _blockSize;
        }
        continue;
      } // if (maxBlock==1)
      //
      // Define the restarting block when maxBlock > 1
      //
      if (nFound > 0) {
        //
        // Copy the converged eigenvalues and residuals.
        //
        int count=0;
        std::vector<int> index(nFound);
        for (j=0; j<_blockSize; j++) {
          if (_normR[j] < _residual_tolerance) {
            index[count] = j;
            (*_evals)[_knownEV + count] = _theta[j];
            _resids[_knownEV + count] = _normR[j];
            count++;
          }
        }
        Teuchos::RefCountPtr<MV> tmpKX = MVT::CloneView( *KX, index );
        for (j=0; j<nFound; j++) {
          index[j] = _knownEV + j;
        }
        MVT::SetBlock( *tmpKX, index, *_evecs );

        int firstIndex = _blockSize;
        //
        // Find the first index where the unconverged eigenvectors start.
        //
        for (j=0; j<_blockSize; j++) {
          if (_normR[j] >= _residual_tolerance) {
            firstIndex = j;
            break;
          }
        }
        //
        // If the firstIndex comes before the number found, then we need
        // to move the converged eigenvectors to the front of the spectral
        // transformation.
        //
        std::vector<ScalarType> tmp_swap_vec(nb*_blockSize);
        while (firstIndex < nFound) {
          for (j=firstIndex; j<_blockSize; j++) {
            if (_normR[j] < _residual_tolerance) {
              //
              // Swap and j-th and firstIndex-th position
              //
              blas.COPY(nb*_blockSize, S[ j ], 1, &tmp_swap_vec[0], 1);
              blas.COPY(nb*_blockSize, S[ firstIndex ], 1, S[ j ], 1 );
              blas.COPY(nb*_blockSize, &tmp_swap_vec[0], 1, S[ firstIndex ], 1 );
              // Swap _theta
              std::swap<MagnitudeType>(_theta[j],_theta[firstIndex]);
              // Swap _normR
              std::swap<MagnitudeType>(_normR[j],_normR[firstIndex]);
              break;
            }
          } 
          for (j=0; j<_blockSize; j++) {
            if (_normR[j] >= _residual_tolerance) {
              firstIndex = j;
              break;
            }
          }
        } // while (firstIndex < nFound)
      } // if (nFound > 0)
      //
      // Define the restarting size
      //
      bStart = ((nb - offSet) > 2 ) ? (nb - offSet)/2 : 0;
      //
      // Define the restarting space and local stiffness matrix
      //
      KK.putScalar( zero );
      for (j=0; j<bStart*_blockSize; j++) {
        KK(j,j) = _theta[j + nFound];
      }
      //
      // Form the restarting space
      //
      int oldCol = nb*_blockSize;
      int newCol = nFound + (bStart+1)*_blockSize;
      newCol = (newCol > oldCol) ? oldCol : newCol;
      std::vector<int> index( oldCol );
      lapack.GEQRF(oldCol, newCol, S.values(), S.stride(), &tau[0], &work[0], lwork, &info);
      lapack.UNGQR(oldCol, newCol, newCol, S.values(), S.stride(), &tau[0], &work[0], lwork, &info);      
      for (i=0; i<oldCol; i++) {
        index[i] = _knownEV + i;
      }
      Teuchos::RefCountPtr<MV> oldX = MVT::CloneView( *X, index );
      index.resize( newCol );
      for (i=0; i<newCol; i++) {
        index[i] = _knownEV + i; 
      }
      Teuchos::RefCountPtr<MV> newX = MVT::CloneView( *X, index );
      Teuchos::RefCountPtr<MV> temp_newX = MVT::Clone( *X, newCol );
      Teuchos::SerialDenseMatrix<int,ScalarType> _Sview( Teuchos::View, S, oldCol, newCol );
      MVT::MvTimesMatAddMv( one, *oldX, _Sview, zero, *temp_newX );
      MVT::MvAddMv( one, *temp_newX, zero, *temp_newX, *newX ); 
      
      if (nFound == 0) {
        offSet++;
      }
      
      _knownEV += nFound;
      maxBlock = (_dimSearch/_blockSize) - (_knownEV/_blockSize);
      //
      // Put random vectors if the Rayleigh Ritz vectors are not enough
      // 
      newCol = nFound + (bStart+1)*_blockSize;
      if (newCol > oldCol) {
        index.resize( nFound );
        for (i=0; i<nFound; i++) {
          index[i] = _knownEV + _blockSize - nFound + i;
        }
        Xnext = MVT::CloneView( *X, index );
        MVT::MvRandom( *Xnext );
        continue;
      }
      //
      // Reset nFound counter.
      // 
      nFound = 0;
      //
    } // while (_iter < _maxIter)
    //
    // Sort the computed eigenvalues, eigenvectors, and residuals
    //
    if ((info==0) && (_knownEV > 0)) {
      // Sort the eigenvalues
      _order.resize(_knownEV);
      _timerSortEval->start();
      ret = _sm->sort( this, _knownEV, &(*_evals)[0], &_order );
      _timerSortEval->stop();
      if (ret != Ok) {
        if (_om->isVerbosityAndPrint(Error)) {
          _os << "ERROR : Sorting in solve()"
              << endl;
        }
        criticalExit = true;
      }
      else {
        // use _order to permute _evecs and _resids
        _MSUtils.permuteVectors(_knownEV,_order,*_evecs,&_resids);
      }
    }

    //
    // Print out a final summary if necessary
    //
    if (_om->isVerbosity( FinalSummary )) {
      currentStatus();
    }

    // Print timing details 
    _timerTotal->stop();
    if (_om->isVerbosity( Anasazi::TimingDetails )) {
      if (_om->doPrint())
        _os <<"********************TIMING DETAILS********************"<<endl;
      Teuchos::TimeMonitor::summarize( _os );
      if (_om->doPrint())
        _os <<"******************************************************"<<endl;
    }

    if (_knownEV == _nev) {
      return Ok;
    }
    else {
      return Unconverged;
    }

  } // end solve()

  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::accuracyCheck(const MV *X, const MV *MX, const MV *Q) const 
    {
      cout.precision(2);
      cout.setf(ios::scientific, ios::floatfield);
      ScalarType tmp;
      
      if (X) {
        if (_MOp.get()) {
          if (MX) {
            tmp = _MSUtils.errorEquality(X, MX, _MOp.get());
            if (_om->doPrint())
              _os << " >> Difference between MX and M*X = " << tmp << endl;
          }
          tmp = _MSUtils.errorOrthonormality(X, _MOp.get());
          if (_om->doPrint())
            _os << " >> Error in X^T M X - I = " << tmp << endl;
        }
        else {
          tmp = _MSUtils.errorOrthonormality(X, 0);
          if (_om->doPrint())
            _os << " >> Error in X^T X - I = " << tmp << endl;
        }
      }
      
      if (Q == 0) {
	_os << endl;     
        return;
      }      

      if (_MOp.get()) {
        tmp = _MSUtils.errorOrthonormality(Q, _MOp.get());
        if (_om->doPrint())
          _os << " >> Error in Q^T M Q - I = " << tmp << endl;
        if (X) {
          tmp = _MSUtils.errorOrthogonality(Q, X, _MOp.get());
          if (_om->doPrint())
            _os << " >> Orthogonality Q^T M X up to " << tmp << endl;
        }
      }
      else {
        tmp = _MSUtils.errorOrthonormality(Q, 0);
        if (_om->doPrint())
          _os << " >> Error in Q^T Q - I = " << tmp << endl;
        if (X) {
          tmp = _MSUtils.errorOrthogonality(Q, X, 0);
          if (_om->doPrint())
            _os << " >> Orthogonality Q^T X up to " << tmp << endl;
        }
      } 
      _os << endl;     
    }
  
  } // End of namespace Anasazi

#endif

// End of file AnasaziBlockDavidson.hpp
