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

/*! \file AnasaziBlockKrylovSchur.hpp
  \brief Implementation of the block Krylov-Schur method
*/

#ifndef ANASAZI_BLOCK_KRYLOV_SCHUR_HPP
#define ANASAZI_BLOCK_KRYLOV_SCHUR_HPP

#include "AnasaziEigensolver.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziSortManager.hpp"
#include "AnasaziOutputManager.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

#if defined(HAVE_COMPLEX)
#define ANSZI_CPLX_CLASS std::complex
#elif  defined(HAVE_COMPLEX_H)
#define ANSZI_CPLX_CLASS ::complex
#endif

/*!  \class Anasazi::BlockKrylovSchur

  \brief This class implements the Restarted Block Krylov Schur Method,
  an iterative method for solving eigenvalue problems.
  
  This method is a block version of the method presented by G.W. Stewart 
  in "A Krylov-Schur Algorithm for Large Eigenproblems", 
  SIAM J. Matrix Anal. Appl., Vol 23(2001), No. 3, pp. 601-614.
  
  \author Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {

  template <class ScalarType, class MV, class OP>
  class BlockKrylovSchur : public Eigensolver<ScalarType,MV,OP> { 
  public:
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
  
    //@{ \name Constructor/Destructor.
    
    //! %Anasazi::BlockKrylovSchur constructor.
    BlockKrylovSchur( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
                      const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sm,
                      const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
                      Teuchos::ParameterList &pl
                      );
    
    //! %Anasazi::BlockKrylovSchur destructor.
    virtual ~BlockKrylovSchur() {};
    //@}
    
    //@{ \name Solver application methods.
    
    /*! \brief This method uses iterate to compute approximate solutions to the
      original problem.  It may return without converging if it has taken the
      maximum number of iterations or numerical breakdown is observed.
    */
    ReturnType solve();
    //@}
    
    //@{ \name Solver status methods.
    
    //! This method returns the computed Ritz values.
    Teuchos::RefCountPtr<const std::vector<ScalarType> > GetRitzValues() const { return(_ritzvalues); };
    
    //! This method returns the Ritz residuals for the computed eigenpairs.
    Teuchos::RefCountPtr<const std::vector<MagnitudeType> > GetRitzResiduals() const { return(_ritzresiduals); };
    
    //! Get the current iteration count.
    int GetNumIters() const { return(_iter); };
    
    //! Get the current restart count of the iteration method.
    /*! Some eigensolvers can perform restarts (i.e.Arnoldi) to reduce memory
      and orthogonalization costs.  For other eigensolvers that don't
      perform restarts (i.e. LOBPCG), this is not a valid stopping criteria.
    */
    int GetNumRestarts() const { return(_numRestarts); };
  
    //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
    int GetBlockSize() const { return(_blockSize); }
    
    //! Get the total length of the Krylov-Schur factorization.
    /*! This number will be the product of the length and the blocksize given by the user.
     */
    int GetKrylovFactLength() const { return(_totallength); };

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
    /*! \brief This method performs a given number of steps of the block Krylov-Schur method
     */
    void iterate( const int steps );
    
    /*! \brief These methods will not be defined.
     */
    BlockKrylovSchur(const BlockKrylovSchur<ScalarType,MV,OP> &method);
    BlockKrylovSchur<ScalarType,MV,OP>& operator=(const BlockKrylovSchur<ScalarType,MV,OP> &method);
    
    /*! \brief Internal methods
     */
    void QRFactorization( MV&, Teuchos::SerialDenseMatrix<int,ScalarType>& );
    void ComputeSchurForm( const bool apply );
    void SortSchurForm( Teuchos::SerialDenseMatrix<int,ScalarType>& H, 
                        Teuchos::SerialDenseMatrix<int,ScalarType>& Q );
    int BlockReduction();
    void BlkOrth( MV& Vec_in, const int j );
    void BlkOrthSing( MV& Vec_in, const int j );
    void ComputeEvecs();
    void Restart();
    void SetBlkTols();
    void CheckBlkArnRed( const int j );
    void CheckSchurVecs( const int j ); 
    
    Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem; 
    Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > _sm; 
    Teuchos::RefCountPtr<OutputManager<ScalarType> > _om; 
    Teuchos::ParameterList &_pl;
    //
    // Information obtained from the eigenproblem
    //
    Teuchos::RefCountPtr<OP> _Op;
    Teuchos::RefCountPtr<OP> _MOp;
    Teuchos::RefCountPtr<MV> _evecs;
    Teuchos::RefCountPtr<std::vector<ScalarType> > _evals;
    const int _nev;  

    int _maxBlocks;
    const int _restarts, _blockSize, _stepSize;
    const MagnitudeType _residual_tolerance;
    int _numRestarts, _iter, _jstart, _jend, _nevblock, _totallength;
    int _offset, _maxoffset;
    bool _isdecompcurrent, _isevecscurrent, _exit_flg, _dep_flg, _error_flg;
    MagnitudeType _schurerror, _dep_tol, _blk_tol, _sing_tol, _def_tol;
    std::vector<int> _order;
    Teuchos::RefCountPtr<MV> _basisvecs;
    Teuchos::SerialDenseMatrix<int,ScalarType> _hessmatrix;
    Teuchos::RefCountPtr<std::vector<ScalarType> > _ritzvalues;
    Teuchos::RefCountPtr<std::vector<MagnitudeType> > _ritzresiduals;
    //
    // Output stream from the output manager
    //
    std::ostream& _os;
    //
    // Internal timers
    //
    bool _restartTimers;
    Teuchos::RefCountPtr<Teuchos::Time> _timerOp, _timerSortEval,
                                        _timerCompSF, _timerSortSF,
                                        _timerCompEvec, _timerQRFact, 
                                        _timerOrtho, _timerTotal;
    //
    // Counters
    //
    int _count_ApplyOp;
  };
  
  //----------------------------------------------------------------------------------------
  // Implementation
  //----------------------------------------------------------------------------------------
  template <class ScalarType, class MV, class OP>
  BlockKrylovSchur<ScalarType,MV,OP>::BlockKrylovSchur(
      const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
      const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sm,
      const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
      Teuchos::ParameterList &pl
      ): 
    _problem(problem), 
    _sm(sm),
    _om(om),
    _pl(pl),
    _Op(_problem->GetOperator()),
    _MOp(_problem->GetM()),
    _evecs(_problem->GetEvecs()), 
    _evals(_problem->GetEvals()), 
    _nev(_problem->GetNEV()), 
    _maxBlocks(_pl.get("Max Blocks", 25)), 
    _restarts(_pl.get("Max Restarts", 0)),
    _blockSize(_pl.get("Block Size", 1 )),
    _stepSize(_pl.get("Step Size", _maxBlocks*(_restarts+1)*_blockSize)),
    _residual_tolerance(_pl.get("Tol", (MagnitudeType)(1.0e-6) )),
    _numRestarts(0), 
    _iter(0), 
    _jstart(0), 
    _jend(0), 
    _nevblock(0),
    _totallength(0),
    _offset(0),
    _maxoffset(0),
    _isdecompcurrent(false),
    _isevecscurrent(false),
    _exit_flg(false),
    _dep_flg(false),
    _error_flg(false),
    _schurerror(1.0), 
    _dep_tol(1.0), 
    _blk_tol(1.0),
    _sing_tol(1.0),
    _def_tol(1.0),
    _os(_om->GetOStream()),
    _restartTimers(_pl.get("Restart Timers", false)),
    _timerOp(Teuchos::TimeMonitor::getNewTimer("Operation Op*x")),
    _timerSortEval(Teuchos::TimeMonitor::getNewTimer("Sorting eigenvalues")),
    _timerCompSF(Teuchos::TimeMonitor::getNewTimer("Computing Schur form")),
    _timerSortSF(Teuchos::TimeMonitor::getNewTimer("Sorting Schur form")),
    _timerCompEvec(Teuchos::TimeMonitor::getNewTimer("Computing eigenvectors")),
    _timerQRFact(Teuchos::TimeMonitor::getNewTimer("QR factorization")),
    _timerOrtho(Teuchos::TimeMonitor::getNewTimer("Orthogonalization")),
    _timerTotal(Teuchos::TimeMonitor::getNewTimer("Total time")),
    _count_ApplyOp(0)
  {     
  }

  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  
  template <class ScalarType, class MV, class OP>
  void 
  BlockKrylovSchur<ScalarType,MV,OP>::SetBlkTols() {
    typedef typename Teuchos::ScalarTraits<MagnitudeType> MGT;
    const MagnitudeType two = 2.0;
    const MagnitudeType eps = SCT::eps();
    _blk_tol = 10.0*MGT::squareroot(eps);
    _sing_tol = 10.0 * eps;
    _dep_tol = 1.0/MGT::squareroot(two);
    _def_tol = eps;
  }

  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  
  template <class ScalarType, class MV, class OP>
  void 
  BlockKrylovSchur<ScalarType,MV,OP>::currentStatus() {
    int i;
    if (_om->doPrint()) {
      _os.setf(ios::scientific, ios::floatfield);
      _os.precision(6);
      _os <<endl;
      _os <<"******************* CURRENT STATUS *******************"<<endl;
      _os <<"The number of iterations performed is " <<_iter<<endl;
      _os <<"The number of restarts performed is " 
          << (_numRestarts>_restarts ? _numRestarts-1 : _numRestarts) 
          << " of " << _restarts << endl;
      _os <<"The current block size is "<<_blockSize<<endl;
      _os <<"The number of eigenvalues requested is "<<_nev<<endl;
      _os <<"The requested residual tolerance is "<<_residual_tolerance<<endl;        
      _os <<"The error for the partial Schur decomposition is "<< _schurerror <<endl;
      _os <<"The number of operations Op*x is "<<_count_ApplyOp<<endl;
      //
      //  Determine status of solver and output information correctly.
      //
      if ( _schurerror < _residual_tolerance ) {
        _os <<endl;
        _os <<"COMPUTED EIGENVALUES                  "<<endl;
      } 
      else {
        if (_exit_flg) {
          _os <<"ERROR: Complete orthogonal basis could not be computed"<<endl;
        }
        else if (_error_flg) {
          _os << "ERROR: Encountered unrecoverable error" << endl;
        }
        _os <<endl;
        _os <<"CURRENT EIGENVALUE ESTIMATES              "<<endl;
      }
      //
      //  Print out current computed eigenvalues.  If we don't have all the requested
      //  eigenvalues yet, print out the ones we have.
      //
      int _nevtemp = _nev;
      if (_jstart < _nevblock) { _nevtemp = _jstart*_blockSize; }
      //
      if (_problem->IsSymmetric() || Teuchos::ScalarTraits<ScalarType>::isComplex) {
        _os << std::setw(16) << std::right << "Eigenvalue" 
            << std::setw(16) << std::right << "Ritz Residual"
            << endl;
        _os <<"------------------------------------------------------"<<endl;
        if ( _nevtemp == 0 ) {
          _os <<"[none computed]"<<endl;
        } else {
          for (i=0; i<_nevtemp; i++) {
            _os << std::setw(16) << std::right << (*_evals)[i]
                << std::setw(16) << std::right << (*_ritzresiduals)[i]
                << endl;
          }
        }
      } 
      else {
        _os << std::setw(16) << std::right << "Real Part"
            << std::setw(16) << std::right << "Imag Part" 
            << std::setw(16) << std::right << "Ritz Residual"
            << endl;
        _os <<"------------------------------------------------------"<<endl;
        if ( _nevtemp == 0 ) {
          _os <<"[none computed]"<<endl;
        } else {
          for (i=0; i<_nevtemp; i++) {
            _os << std::setw(16) << std::right << (*_evals)[i]
                << std::setw(16) << std::right << (*_evals)[_nev+i]
                << std::setw(16) << std::right << (*_ritzresiduals)[i]<<endl;
          }
        }
      }
      _os <<"******************************************************"<<endl << endl;;
    }        
  }

  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  
  template <class ScalarType, class MV, class OP>
  ReturnType 
  BlockKrylovSchur<ScalarType,MV,OP>::solve () 
  {
    Teuchos::TimeMonitor LocalTimer(*_timerTotal,_restartTimers);

    if ( _restartTimers ) {
      _timerOp->reset();
      _timerCompSF->reset();
      _timerSortSF->reset();
      _timerSortEval->reset();
      _timerCompEvec->reset();
      _timerQRFact ->reset();
      _timerOrtho->reset();
      _count_ApplyOp = 0;
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
    // Retrieve the initial vector and operator information from the Anasazi::Eigenproblem.
    //
    Teuchos::RefCountPtr<MV> ivec = _problem->GetInitVec();
    //
    if ( ivec.get() == 0 ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : Initial vector is not specified, set initial vector in eigenproblem "<<endl;
      return Failed;
    }

    int dim = MVT::GetVecLength( *ivec );
    //
    // Check that the maximum number of blocks for the eigensolver is a positive number
    //    
    if ( _maxBlocks<=0 ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : maxBlocks = "<< _maxBlocks <<" [ should be positive number ] " << endl;
      return Failed;
    } 

    if ( _maxBlocks*_blockSize < _nev ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : Krylov subspace dimension (maxBlocks*blockSize) = "<< _maxBlocks*_blockSize 
            << " [ should be greater than "<< _nev << " ] " << endl;
      return Failed;
    } 
    //
    // Check that the Krylov subspace is larger than the number of eigenvalues requested
    //
    if ( _maxBlocks*_blockSize < _nev ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : Krylov subspace dimension (maxBlocks*blockSize) = "<< _maxBlocks*_blockSize 
            << " [ should be greater than "<< _nev << " ] " << endl;
      return Failed;
    } 
    //
    // If the Krylov subspace dimension is the same size as the number of requested eigenvalues,
    // then we must be computing all of them.
    //
    if ( (_maxBlocks*_blockSize == _nev) && (_nev != dim) ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error ))
        _os << "ERROR : Krylov subspace dimension (maxBlocks*blockSize) = "<< _maxBlocks*_blockSize 
            << " [ should be greater than "<< _nev << " ] " << endl;
      return Failed;
    }
    //
    // If the Krylov subspace dimension is larger than the dimension of the operator, reset
    // the maximum number of blocks accordingly.
    //    
    if (_maxBlocks*_blockSize > dim ) {
      if (_om->isVerbosityAndPrint( Anasazi::Warning ))
        _os << "WARNING : Krylov subspace dimension (maxBlocks*blockSize) = "<< _maxBlocks*_blockSize 
            <<" [ should not be greater than " << dim << " ] " << endl;
      
      // Set the maximum number of blocks in the factorization below the dimension of the space.
      _maxBlocks = dim / _blockSize;
      
      if (_om->isVerbosityAndPrint( Anasazi::Warning ))
        _os << "WARNING : maxBlocks reset to "<< _maxBlocks << endl;
    }
    
    if ( (_maxBlocks*_blockSize == _nev) && (_nev != dim) ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error ))
        _os << "ERROR : Krylov subspace dimension (maxBlocks*blockSize) = "<< _maxBlocks*_blockSize 
            << " [ should be greater than "<< _nev << " ] " << endl;
      return Failed;
    }
    //
    // Check that the step size for the eigensolver is a positive number.
    //
    if ( _stepSize<=0 ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : stepSize = "<< _stepSize <<" [ should be positive number ] " << endl;
      return Failed;
    } 

    //
    // Determine _nevblock : how many blocks it will take to contain the _nev eigenvalues/vectors
    // NOTE: An additional block is kept if _nev is a perfect multiple of _blockSize because of the
    // potential presence of complex eigenvalue pairs.  Additional blocks can be retained, up to
    // _maxoffset if the block ends with one eigenvalue of a complex conjugate pair.
    //
    _nevblock = _nev/_blockSize + 1;
    //
    // Use alternate settings if the number of eigenvalues requested is equal to the dimension
    if ( _nev == dim ) {
      _nevblock = _nev/_blockSize;  
      if (_nev%_blockSize) 
        _nevblock++;    
    }
    _maxoffset = (_maxBlocks-_nevblock)/2;
    _totallength = _blockSize*_maxBlocks;
    //
    // Reinitialize internal data and pointers, prepare for solve
    // 
    _numRestarts=0; _iter=0; _jstart=0; _jend=0; 
    _isdecompcurrent=false; _isevecscurrent=false; 
    _exit_flg=false; _dep_flg=false; _error_flg=false;
    _schurerror=1.0; 
    //
    // Make room for the Arnoldi vectors and F.
    //
    _basisvecs = MVT::Clone( *ivec, (_maxBlocks+1)*_blockSize );
    //
    // Create the rectangular Hessenberg matrix
    //
    _hessmatrix.shape((_maxBlocks+1)*_blockSize, _maxBlocks*_blockSize); 
    //
    // Create the vectors for eigenvalues and their residual errors and
    // initialize them.
    //
    _ritzvalues = Teuchos::rcp(new std::vector<ScalarType>(2*_totallength));
    _ritzresiduals = Teuchos::rcp(new std::vector<MagnitudeType>(_totallength)); 
    _order.resize( _totallength );
    const MagnitudeType one = 1.0;
    const ScalarType zero = 0.0;
    for (int i=0; i< _totallength; i++) {
      (*_ritzvalues)[i] = zero; (*_ritzvalues)[_totallength+i] = zero;
      (*_ritzresiduals)[i] = one;
    }
    //
    // Set the tolerances for block orthogonality
    //
    SetBlkTols();  
    //
    // Right now the solver will just go the remaining iterations, but this design will allow
    // for checking of the residuals every so many iterations, independent of restarts.
    //
    while (_schurerror > _residual_tolerance && _numRestarts <= _restarts 
           && !_exit_flg && !_error_flg) {
      iterate( _stepSize );
    }
    //
    // Compute the current approximate eigenvectors before returning.
    //
    if (!_error_flg) {
      Teuchos::TimeMonitor CompTimer(*_timerCompEvec);
      ComputeEvecs();    
    }
    //
    // Print out a final summary before returning.
    //
    if (_om->isVerbosity( FinalSummary )) {
      currentStatus();
    }

    // Print timing details 
    _timerTotal->stop();

    // Reset format that will be used to print the summary
    Teuchos::TimeMonitor::format().setPageWidth(54);

    if (_om->isVerbosity( Anasazi::TimingDetails )) {
      if (_om->doPrint())
        _os <<"********************TIMING DETAILS********************"<<endl;
      Teuchos::TimeMonitor::summarize( _os );
      if (_om->doPrint())
        _os <<"******************************************************"<<endl;
    }

    // Return Failed if we did not meet the specified tolerance
    if (_schurerror > _residual_tolerance) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : Specified residual tolerance ( "<< _residual_tolerance<<" ) was not reached"<< endl;
      return Unconverged;
    }

    // Output warning if numerical breakdown was detected in orthogonalization
    if (_exit_flg) {
      if (_om->isVerbosityAndPrint( Anasazi::Warning )) 
        _os << "WARNING : Numerical breakdown detected in Anasazi::BlockKrylovSchur"<< endl;
    }
    
    if (_error_flg) {
      return Failed;
    }
    
    return Ok;
  }

  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  
  template <class ScalarType, class MV, class OP>
  void 
  BlockKrylovSchur<ScalarType,MV,OP>::iterate(const int steps) {
    int i=0;
    int tempsteps = steps;
    int blk_red_steps = 0;
    //
    // If this is the first steps of Block Krylov Schur, initialize the first block of _basisvecs
    //
    if (!_iter) {
      int numIVecs = MVT::GetNumberVecs( *(_problem->GetInitVec()) );
      if (numIVecs > _blockSize)
        numIVecs = _blockSize;
      std::vector<int> index( numIVecs );
      for (i=0; i<numIVecs; i++) {
        index[i] = i;
      }
      //
      // Copy the first _blockSize of the initial vectors into the first _blockSize
      // of _basisvecs, any additional vectors will be ignored.
      //
      MVT::SetBlock( *(_problem->GetInitVec()), index, *_basisvecs );
      //
      // Augment the initial vectors with random vectors if necessary
      //
      numIVecs = _blockSize - MVT::GetNumberVecs( *(_problem->GetInitVec()) );
      if (numIVecs > 0) {
        index.resize(numIVecs);
        for (i=0; i<numIVecs; i++)
          index[i] = _blockSize - numIVecs + i;
        Teuchos::RefCountPtr<MV> tmpIVec = MVT::CloneView( *_basisvecs, index );
        MVT::MvRandom( *tmpIVec );
      }  
      //
      // Orthogonalize the first block of vectors.
      //      
      index.resize(_blockSize);
        for (i=0; i<_blockSize; i++)
          index[i] = i;
      Teuchos::RefCountPtr<MV> U_vec = MVT::CloneView( *_basisvecs, index );
      Teuchos::SerialDenseMatrix<int,ScalarType> G10( _blockSize,_blockSize );
      QRFactorization( *U_vec, G10 );
    }
    //
    // Leave the iteration method now if the orthogonal subspace can't be extended.
    //
    if (_exit_flg) { return; }        
    //
    // Now we go the number of steps requested by the user.  This may cause
    // a restart or hit the number of maximum iterations (restarts).  
    //
    while(tempsteps > 0 && _numRestarts <= _restarts 
          && !_exit_flg && !_error_flg && _schurerror > _residual_tolerance) {
      _isevecscurrent = false;
      // If we don't need to restart, just get it over with and return.
      if (_jstart+tempsteps < _maxBlocks) {
        _jend = _jstart+tempsteps;
        blk_red_steps = BlockReduction();
        if (_error_flg)  return;
        //
        // Move the pointer and update the iteration count.
        //
        _iter += blk_red_steps;
        tempsteps -= blk_red_steps;
        _jstart += blk_red_steps;
        ComputeSchurForm( false );
        _isdecompcurrent = false;
        // Output current information if necessary
        if (_om->isVerbosity( IterationDetails )) {
          currentStatus();
        }
      }
      // Finish off this factorization and restart.
      else {
        _jend = _maxBlocks;
        blk_red_steps = BlockReduction();
        if (_error_flg)  return;
        //
        // Move the pointer and update the iteration count.
        //
        _iter += blk_red_steps;
        tempsteps -= blk_red_steps;
        _jstart += blk_red_steps;
        //
        //  Compute the Schur factorization and prepare for a restart.  Don't
        //  compute restart if at end of iterations.
        //
        if (_numRestarts < _restarts) {
          ComputeSchurForm( true );  
          Restart();  
          _isdecompcurrent = true;
          _numRestarts++;
        } else {
          ComputeSchurForm( false );
          _numRestarts++;
          _isdecompcurrent = false;
        }
        //
        // Output current information if necessary
        if (_om->isVerbosity( IterationDetails )) {
          currentStatus();
        }
      }
    }
  }


  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  
  template <class ScalarType, class MV, class OP>
  int 
  BlockKrylovSchur<ScalarType,MV,OP>::BlockReduction () {

    int i,j;
    ReturnType ret;
    std::vector<int> index( _blockSize );
    Teuchos::RefCountPtr<MV> U_vec, F_vec;
    
    for ( j = _jstart; j < _jend; j++ ) {
      //
      // Associate the j-th block of _basisvecs with U_vec.
      //
      for ( i=0; i<_blockSize; i++ ) {
        index[i] = j*_blockSize+i;
      }
      U_vec = MVT::CloneView( *_basisvecs, index );
      //
      // Associate (j+1)-st block of ArnoldiVecs with F_vec.
      //
      //for ( i=0; i<_blockSize; i++ ) {
      //        index[i] = (j+1)*_blockSize+i;
      //}
      F_vec = MVT::Clone( *_basisvecs, _blockSize );
      //F_vec = MVT::CloneView( *_basisvecs, index );
      //
      //  Compute F_vec = OP * U_vec
      //
      {
	Teuchos::TimeMonitor OpTimer(*_timerOp);
	ret = OPT::Apply( *_Op, *U_vec, *F_vec ); 
      }
      _count_ApplyOp += MVT::GetNumberVecs(*U_vec);
      if (ret != Ok) {
        // Apply() failed. Return error code.
        if (_om->isVerbosityAndPrint(Error)) {
          _os << "ERROR : Applying Op operator in BlockReduction" << endl;
        }
        _error_flg = true;
        return -1;
      }
      //
      // Use previous dependency information to decide which orthogonalization
      // method to use for the new block.  The global variable _dep_flg tells us
      // if we've had problems with orthogonality before.  If no problems have
      // been detected before we will use standard block orthogonalization.
      // However, if there are problems with this, we will use a more stringent
      // approach.
      //
      if (!_dep_flg) {
        BlkOrth( *F_vec, j );
      }
      //
      // If any block dependency was detected previously, then the more stringent
      // orthogonalization will be used.  If this method can't resolve the
      // dependency, then the _exit_flg will be set indicating that we can't proceed
      // any further.
      //
      if (_dep_flg) {
        BlkOrthSing( *F_vec, j );
      }
      //
      // If we cannot go any further with the factorization, then we need to exit
      // this method.
      //
      if (_exit_flg) { break; }
      if (_error_flg) { return -1; }
    }
    //
    // Return the number of steps in the block reduction that were accomplished.
    //
    return (j-_jstart);

  } // end BlockReduction()
  
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  template <class ScalarType, class MV, class OP>
  void 
  BlockKrylovSchur<ScalarType,MV,OP>::BlkOrth( MV& Vec_in, const int j ) {
    //
    // Orthogonalization is first done between the new block of
    // vectors and all previous blocks, then the vectors within the
    // new block are orthogonalized.

    // Local timer
    Teuchos::TimeMonitor LocalTimer(*_timerOrtho);

    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    const int max_num_orth = 2;
    int i, row_offset, col_offset;
    std::vector<int> index( _blockSize );
    std::vector<MagnitudeType> norm1( _blockSize );
    std::vector<MagnitudeType> norm2( _blockSize );
    //
    // Associate (j+1)-st block of ArnoldiVecs with F_vec.
    //
    for ( i=0; i<_blockSize; i++ ) {
      index[i] = (j+1)*_blockSize+i;
    }
    Teuchos::RefCountPtr<MV> F_vec = MVT::CloneView( *_basisvecs, index );
    MVT::MvAddMv( one, Vec_in, zero, Vec_in, *F_vec );
    //
    // Zero out the full block column of the Hessenberg matrix
    // even though we're only going to set the coefficients in
    // rows [0:(j+1)*_blockSize-1]
    //
    int num_prev = (j+1)*_blockSize;
    //
    Teuchos::SerialDenseMatrix<int,ScalarType> tmp_hess( Teuchos::View, _hessmatrix, num_prev, _blockSize, 0, j*_blockSize );
    tmp_hess.putScalar( zero );
    //
    // Grab all previous Arnoldi vectors
    //
    index.resize( num_prev );
    for (i=0; i<num_prev; i++){
      index[i] = i;
    }
    Teuchos::RefCountPtr<MV> V_prev = MVT::CloneView( *_basisvecs, index );
    //
    // Create a matrix to store the product trans(V_prev)*B*F_vec
    //
    Teuchos::SerialDenseMatrix<int,ScalarType> dense_mat(num_prev, _blockSize );
    //
    MVT::MvNorm( *F_vec, &norm1 );
    //
    // Check the norm of the candidate block of vectors to make sure they're
    // not zero.  [ This might happen if the matrix is the zero matrix ]
    //
    for (i=0; i<_blockSize; i++) {
      if (norm1[i] == zero) {
        _dep_flg = true;
        if (_om->isVerbosityAndPrint( Debug )){
          _os << "Col " << num_prev+i << " is the zero vector" << endl;
          _os << endl;
        }
      }          
    }
    //
    // Perform two steps of block classical Gram-Schmidt so that
    // F_vec is B-orthogonal to the columns of V_prev.
    //
    for ( int num_orth=0; num_orth<max_num_orth; num_orth++ ) {
      //
      // Compute trans(V_prev)*B*F_vec and store in the j'th diagonal
      // block of the Hessenberg matrix
      //
      _problem->InnerProd( *V_prev, *F_vec, dense_mat );
      //
      // Update the orthogonalization coefficients for the j-th block
      // column of the Hessenberg matrix.
      //
      tmp_hess += dense_mat;
      //
      // F_vec <- F_vec - V(0:(j+1)*block-1,:) * H(0:num_prev-1,j:num_prev-1)
      //
      MVT::MvTimesMatAddMv( -one, *V_prev, dense_mat, one, *F_vec );
      //
    } // end for num_orth=0;...)
    //
    MVT::MvNorm( *F_vec, &norm2 );
    //
    // Check to make sure the new block of Arnoldi vectors are
    // not dependent on previous Arnoldi vectors.  
    //
    for (i=0; i<_blockSize; i++){
      if (norm2[i] < norm1[i] * _blk_tol) {
        _dep_flg = true;
        if (_om->isVerbosityAndPrint( Warning )) {
          _os << "WARNING : Col " << num_prev+i << " is dependent on previous "
               << "Arnoldi vectors in V_prev" << endl;
          _os << endl;
        }
      }
    } // end for (i=0;...)
    //
    if (_om->isVerbosity( OrthoDetails )) {
      if (_om->doPrint()) {
        _os << "Checking Orthogonality after BlkOrth()"
            << " Iteration: " << j << endl<<endl;
      }
      CheckBlkArnRed(j);
      // This may have triggered an error flag, in which case we
      // should quit.
      if (_error_flg) {
        return;
      }
    }
    //
    // If dependencies have not already been detected, compute
    // the QR factorization of the next block. Otherwise,
    // this block of Arnoldi vectors will be re-computed via and
    // implementation of A. Ruhe's block Arnoldi.
    //
    if (!_dep_flg) {
      //
      // Compute the QR factorization of F_vec
      //
      row_offset = (j+1)*_blockSize; col_offset = j*_blockSize;
      Teuchos::SerialDenseMatrix<int,ScalarType> sub_blockSize_hess(Teuchos::View, _hessmatrix, _blockSize, 
                                                                    _blockSize, row_offset, col_offset);
      
      QRFactorization( *F_vec, sub_blockSize_hess );
    }
    //
  }  // end BlkOrth()
  
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  template <class ScalarType, class MV, class OP>
  void 
  BlockKrylovSchur<ScalarType,MV,OP>::BlkOrthSing( MV& Vec_in, const int j ) {
    //
    // This is a variant of A. Ruhe's block Arnoldi
    // The orthogonalization of the vectors F_vec is done
    // one at a time. If a dependency is detected, a random
    // vector is added and orthogonalized against all previous
    // Arnoldi vectors.
    
    // Local timer
    Teuchos::TimeMonitor LocalTimer(*_timerOrtho);

    const int IntOne = 1;
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    int i, num_prev;
    std::vector<int> index( _blockSize );
    Teuchos::SerialDenseVector<int,ScalarType> dense_vec;
    std::vector<MagnitudeType> norm1(IntOne);
    std::vector<MagnitudeType> norm2(IntOne);
    //
    // Place the candidate vectors Vec_in into the (j+1)-st block of ArnoldiVecs.
    //
    for ( i=0; i<_blockSize; i++ ) {
      index[i] = (j+1)*_blockSize+i;
    }
    MVT::SetBlock( Vec_in, index, *_basisvecs ); 
    //
    // Zero out the full block column of the Hessenberg matrix
    //
    int n_row = _hessmatrix.numRows();
    //
    Teuchos::SerialDenseMatrix<int,ScalarType> tmp_hess( Teuchos::View, _hessmatrix, n_row, _blockSize, 0, j*_blockSize );
    tmp_hess.putScalar( zero );
    //
    Teuchos::RefCountPtr<MV>  q_vec, Q_vec, tptr;
    tptr = MVT::Clone( *_basisvecs, IntOne ); 
    //
    // Start a loop to orthogonalize each of the _blockSize
    // columns of the (j+1)-st block of _basisvecs against all 
    // the others.
    //
    std::vector<int> index2(IntOne);
    //
    for (int iter=0; iter<_blockSize; iter++){
      num_prev = (j+1)*_blockSize + iter; // number of previous _basisvecs
      dense_vec.size(num_prev);
      Teuchos::SerialDenseMatrix<int,ScalarType> tmp_hess_col( Teuchos::View, _hessmatrix, num_prev, 1, 0, j*_blockSize+iter );
      //
      // Grab the next column of _basisvecs
      //
      index2[0] = num_prev;
      q_vec = MVT::CloneView( *_basisvecs, index2 ); 
      //
      // Grab all previous columns of _basisvecs
      //
      index.resize(num_prev);
      for (i=0; i<num_prev; i++){
        index[i] = i;
      }
      Q_vec = MVT::CloneView( *_basisvecs, index ); 
      //
      // Create matrix to store product trans(Q_vec)*B*q_vec
      //
      // Do one step of classical Gram-Schmidt B-orthogonalization
      // with a 2nd correction step if needed.
      //
      MVT::MvNorm( *q_vec, &norm1 );
      //
      // Leave if this is the zero vector, there is no more we can do here.
      //
      if (norm1[0] == zero) { 
        if (_om->isVerbosityAndPrint( Debug )) {
          _os << "Column " << num_prev << " of _basisvecs is the zero vector" 
               << endl<<endl;
        }
        _exit_flg = true; 
        break; 
      }
      //
      // Compute trans(Q_vec)*B*q_vec
      //
      _problem->InnerProd( *Q_vec, *q_vec, dense_vec );
      //
      // Sum results [0:num_prev-1] into column (num_prev-_blockSize)
      // of the Hessenberg matrix
      //
      tmp_hess_col += dense_vec;
      //
      // Compute q_vec<- q_vec - Q_vec * dense_vec
      //
      MVT::MvTimesMatAddMv( -one, *Q_vec, dense_vec, one, *q_vec );
      //
      MVT::MvNorm( *q_vec, &norm2 );
      //
      if (norm2[0] < norm1[0] * _dep_tol) {
        //
        // Repeat process with newly computed q_vec
        //
        // Compute trans(Q_vec)*q_vec
        //
        _problem->InnerProd( *Q_vec, *q_vec, dense_vec );
        //
        // Sum results [0:num_prev-1] into column (num_prev-_blockSize)
        // of the Hessenberg matrix
        //
        tmp_hess_col += dense_vec;
              //
        // Compute q_vec<- q_vec - Q_vec * dense_vec
        //
        MVT::MvTimesMatAddMv( -one, *Q_vec, dense_vec, one, *q_vec );
        //
        MVT::MvNorm( *q_vec, &norm2 );
      }
      //
      // Check for linear dependence
      //
      if (norm2[0] < norm1[0] * _sing_tol) {
        if (_om->isVerbosityAndPrint( Warning )) {
          _os << "WARNING : Column " << num_prev << " of _basisvecs is dependent" 
              << endl<<endl;
        }
        //
        // Create a random vector and orthogonalize it against all
        // previous cols of _basisvecs
        // We could try adding a random unit vector instead -- not
        // sure if this would make any difference.
        //
        MVT::MvRandom( *tptr );
        MVT::MvNorm( *tptr, &norm1 );
        //
        // This code  is automatically doing 2 steps of B-orthogonalization
        // after adding a random vector. We could do one step of
        // orthogonalization with a correction step if needed.
        //
        for (int num_orth=0; num_orth<2; num_orth++){
          _problem->InnerProd( *Q_vec, *tptr, dense_vec );
          // Note that we don't change the entries of the
          // Hessenberg matrix when we orthogonalize a
          // random vector
          MVT::MvTimesMatAddMv( -one, *Q_vec, dense_vec, one, *tptr );
        }
        //
        MVT::MvNorm( *tptr, &norm2 );
        //
        if (norm2[0] > norm1[0] * _sing_tol){
          // Copy vector into the current column of _basisvecs
          MVT::MvAddMv( one, *tptr, zero, *tptr, *q_vec );
          MVT::MvNorm( *q_vec, &norm2 );
          //
          // Normalize the new q_vec
          //
          ScalarType rjj = one/norm2[0];
          MVT::MvAddMv( rjj, *q_vec, zero, *q_vec, *q_vec );
          //
          // Enter a zero in the [(j+1)*_blockSize + iter] row in the
          // [(j*_blockSize + iter] column of the Hessenberg matrix
          //
          _hessmatrix((j+1)*_blockSize+iter, j*_blockSize+iter) = zero;
        }
        else {
          // Can't produce a new orthonormal basis vector
          // Clean up and exit this block Arnoldi factorization!
          _exit_flg = true;
          return;
        }
      }
      else {
        //
        // Normalize the new q_vec
        //
        ScalarType rjj = one/norm2[0];
        MVT::MvAddMv( rjj, *q_vec, zero, *q_vec, *q_vec );
        //
        // Enter norm of q_vec to the [(j+1)*_blockSize + iter] row
        // in the [(j*_blockSize + iter] column of the Hessenberg matrix
        //
        _hessmatrix((j+1)*_blockSize+iter, j*_blockSize+iter) = norm2[0];
      } // end else ...
    } // end for (i=0;...)
    //
    if (_om->isVerbosity( OrthoDetails )) {
      if (_om->doPrint()) {
        _os << "Checking Orthogonality after BlkOrthSing()"
            << " Iteration: " << j << endl<<endl;
      }
      CheckBlkArnRed(j);
      // This may have triggered an error flag, in which case we
      // should quit.
      if (_error_flg) {
        return;
      }
    }
  } // end BlkOrthSing()

  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  template <class ScalarType, class MV, class OP>
  void 
  BlockKrylovSchur<ScalarType,MV,OP>::QRFactorization (MV& VecIn, 
						       Teuchos::SerialDenseMatrix<int,ScalarType>& R) {
    // local timer
    Teuchos::TimeMonitor LocalTimer(*_timerQRFact);

    int i,j;
    int nb = MVT::GetNumberVecs( VecIn ); assert (nb == _blockSize);
    const int IntOne=1;
    std::vector<int> index;
    std::vector<int> index2(IntOne);
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    bool addvec = false; 
    //
    std::vector<MagnitudeType> norm1(IntOne), norm2(IntOne);
    Teuchos::RefCountPtr<MV> qj, Qj, tptr;
    tptr = MVT::Clone( *_basisvecs, IntOne ); 
    //
    // Zero out the array that will contain the Fourier coefficients.
    //
    R.putScalar( zero );
    //
    // Start the loop to orthogonalize the nb columns of VecIn.
    //
    for ( j=0; j<nb; j++ ) {
      //
      // Grab the j-th column of VecIn (the first column is indexed to 
      // be the zero-th one).
      //
      index2[0] = j;
      qj = MVT::CloneView( VecIn, index2 );
      //
      // If we are beyong the 1st column, orthogonalize against the previous
      // vectors in the current block.
      //
      if ( j ) {
        Teuchos::SerialDenseMatrix<int,ScalarType> R_col( Teuchos::View, R, j, 1, 0, j );
        index.resize( j );
        for ( i=0; i<j; i++ ) {
          index[i] = i;
        }
        //
        // Grab the first j columns of VecIn (that are now an orthogonal
        // basis for first j columns of the entering VecIn).
        //
        Qj = MVT::CloneView( VecIn, index );
        Teuchos::SerialDenseVector<int,ScalarType> rj(j);
        _problem->MvNorm( *qj, &norm1 );
        //
        // Do one step of classical Gram-Schmidt orthogonalization
        // with a second correction step if needed
        //
        // Determine the Fourier coefficients for B-orthogonalizing column
        // j of VecIn against columns 0:j-1 of VecIn. In other words,
        // result = trans(Qj)*B*qj.
        //
        _problem->InnerProd( *Qj, *qj, rj );
        //
        // Sum results[0:j-1] into column j of R.
        //
        R_col += rj;
        //
        // Compute qj <- qj - Qj * rj.
        //
        MVT::MvTimesMatAddMv( -one, *Qj, rj, one, *qj );
        //
        _problem->MvNorm( *qj, &norm2 );
        //
        if (norm2[0] < norm1[0] * _dep_tol){
          //
          // Repeat process with newly computed qj
          //
          _problem->InnerProd( *Qj, *qj, rj );
          //
          // Sum results[0:j-1] into column j of R.
          //
          R += rj;
          //
          // Compute qj <- qj - Qj * rj.
          //
          MVT::MvTimesMatAddMv( -one, *Qj, rj, one, *qj );
          //
          _problem->MvNorm( *qj, &norm2 );
        }
        //
        // Check for dependencies
        //
        if (_iter) {
          // This is not the 1st block. A looser tolerance is used to
          // determine dependencies. If a dependency is detected, a flag
          // is set so we can back out this method and out of BlkOrth.
          // The method BlkOrthSing is used to construct the new block
          // of orthonormal basis vectors one at a time. If a dependency
          // is detected within this method, a random vector is added
          // and orthogonalized against all previous basis vectors.
          //
          if (norm2[0] < norm1[0] * _blk_tol) {
            if (_om->isVerbosityAndPrint( Debug )) {
              _os << "Column " << j << " of current block is dependent"<<endl;
            }
            _dep_flg = true;
            return;
          }
        }
        else {
          // This is the 1st block of basis vectors.
          // Use a tighter tolerance to determine dependencies, because
          // if a dependency is detected we will be adding a random
          // vector and orthogonalizing it against previous vectors
          // in the 1st block
          //
          if (norm2[0] < norm1[0] * _sing_tol) {
            // The 1st block of vectors are dependent
            // Add a random vector and orthogonalize it against
            // previous vectors in block.
            //
            addvec = true;
            Teuchos::SerialDenseVector<int,ScalarType> tj(j);
            //
            MVT::MvRandom( *tptr );
            _problem->MvNorm( *tptr, &norm1 );
            //
            for (int num_orth=0; num_orth<2; num_orth++){
              _problem->InnerProd( *Qj, *tptr, tj );
              MVT::MvTimesMatAddMv( -one, *Qj, tj, one, *tptr );
            }
            _problem->MvNorm( *tptr, &norm2 );
            //
            if (norm2[0] > norm1[0] * _sing_tol){
              // Copy vector into current column of _basisvecs
              MVT::MvAddMv( one, *tptr, zero, *tptr, *qj );
            }
            else {
              _exit_flg = true;
              return;
            }
          }
        } // if (_iter) ...
      } // if (j) ...
      //
      // If we have not exited, compute the norm of column j of
      // VecIn (qj), then normalize qj to make it into a unit vector
      //
      std::vector<MagnitudeType> normq(IntOne);
      _problem->MvNorm( *qj, &normq );
      //
      ScalarType rjj = one / normq[0];
      MVT::MvAddMv ( rjj, *qj, zero, *qj, *qj );
      //
      if (addvec){
        // We've added a random vector, so
        // enter a zero in j'th diagonal element of R
        R(j,j) = zero;
      }
      else {
        R(j,j) = normq[0];
      }
    } // for (j=0; j<nb; j++) ...
    //
  }
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::ComputeEvecs() {
    //
    Teuchos::TimeMonitor LocalTimer(*_timerCompEvec);
    int i=0,j=0,k=0;
    const int IntOne=1;
    int n=_jstart*_blockSize, info=0;
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::LAPACK<int,ScalarType> lapack;
    Teuchos::LAPACK<int,MagnitudeType> lapack_mag;
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::RefCountPtr<MV> basistemp;
    Teuchos::SerialDenseMatrix<int,ScalarType> Q(n,n);
    std::vector<int> index( n );
    //
    // Initialize the eigenvectors to zero.
    //
    MVT::MvInit( *_evecs, 0.0 );
    //
    //  Set the index array.
    //
    for (k=0; k<n; k++) { index[k] = k; }
    //
    //  By some chance the number of eigenvalues requested may be more than
    //  the size of the factorization.  To prevent a seg fault in creating the
    //  eigenvectors, determine if this is the current situation.
    //
    int curr_nev = _nev;
    if ( n < _nev ) { curr_nev = n; }
    //
    //  Get a view into the current Hessenberg matrix.
    //
    Teuchos::SerialDenseMatrix<int,ScalarType> Hj(Teuchos::Copy, _hessmatrix, n, n, i, j);
    //
    //  If the Krylov-Schur decomposition is not current, compute residuals
    //  like we are going to restart to update decomposition.
    //
    if (!_isdecompcurrent) { 
      //
      // Update the Krylov basis using the Schur form.  
      // Take into account that deflated blocks need not be updated.
      //        
      SortSchurForm( Hj, Q );
      if (_error_flg) {
        return;
      }
      basistemp = MVT::Clone( *_basisvecs, n );
      Teuchos::RefCountPtr<MV> basistemp2 = MVT::CloneView( *_basisvecs, index );
      MVT::MvTimesMatAddMv ( one, *basistemp2, Q, zero, *basistemp );
    } 
    else {
      //
      // We can aquire the Ritz vectors from the current decomposition.
      //
      basistemp = MVT::CloneCopy( *_basisvecs, index );
    }
    //
    // Check the Schur form.
    //
    if (_om->isVerbosity( OrthoDetails )) {
      if (_om->doPrint()) {
        _os << "Checking Schur Factorization"
            << " Iteration: " << _jstart << endl<<endl;
      }
      CheckSchurVecs( _jstart );
    }
    //
    //  If the operator is symmetric, then the Ritz vectors are the eigenvectors.
    //  So, copy the Ritz vectors.  Else, we need to compute the eigenvectors of the
    //  Schur form to compute the eigenvectors of the non-symmetric operator.
    //
    if (_problem->IsSymmetric()) {
      index.resize( curr_nev );
      for (k=0; k<curr_nev; k++)
        index[k] = k;
      MVT::SetBlock( *basistemp, index, *_evecs );
    } 
    else {  
      //
      //  Now compute the eigenvectors of the Schur form
      //  Reset the dense matrix and compute the eigenvalues of the Schur form.
      //

      // Allocate the work space. This space will be used below for calls to:
      // * TREVC (requires 3*N for real, 2*N for complex) 
      int lwork = 3*n;
      std::vector<ScalarType> work( lwork );
      std::vector<MagnitudeType> rwork( n );
      char side = 'R';
      int mm; 
      const int ldvl = 1;
      ScalarType vl[ ldvl ];
      lapack.TREVC( side, n, Hj.values(), Hj.stride(), vl, ldvl,
                    Q.values(), Q.stride(), n, &mm, &work[0], &rwork[0], &info );
      assert(info==0);
      //
      //  Convert back to approximate eigenvectors of the operator and compute their norms.
      //
      std::vector<MagnitudeType> evecnrm( n );
      Teuchos::RefCountPtr<MV> evecstemp = MVT::Clone( *_basisvecs, n );
      MVT::MvTimesMatAddMv( one, *basistemp, Q, zero, *evecstemp );
      MVT::MvNorm( *evecstemp, &evecnrm );
      //
      // Sort the eigenvectors.
      //
      int conjprs=0;
      std::vector<int> index2(IntOne), indexi;
      Teuchos::RefCountPtr<MV> evecstempr, evecr1;
      ScalarType t_evecnrm;
      i = 0;

      while ( i < curr_nev ) {        
        if ((*_ritzvalues)[_totallength+i] != zero) {
          t_evecnrm = one/lapack_mag.LAPY2(evecnrm[i],evecnrm[i+1]);
          index2[0] = i;

          // Copy the real part of the eigenvector.  Scale by square-root of 2 to normalize the vector.
          evecstempr = MVT::CloneView( *evecstemp, index2 );
          evecr1 = MVT::CloneView( *_evecs, index2 );
          MVT::MvAddMv( t_evecnrm, *evecstempr, zero, *evecstempr, *evecr1 );
          //
          // Make sure we don't set the eigenvalue for a extra conjugate pair, we may not have room.
          //
          if ( i+1 < _nev ) {
            index2[0] = i+1;
            evecr1 = MVT::CloneView( *_evecs, index2 );
            MVT::MvAddMv( t_evecnrm, *evecstempr, zero, *evecstempr, *evecr1 );
          }
          // Note where imaginary part of eigenvector is.
          indexi.push_back(i+1);
          
          // Increment counters.
          conjprs++;
          i = i+2;
        } else {
          // Copy the real part of the eigenvector, scale to be norm one.
          // We don't have to do anything for the imaginary
          // part since we initialized the vectors to zero.
          index2[0] = i;
          evecstempr = MVT::CloneView( *evecstemp, index2 );
          evecr1 = MVT::CloneView( *_evecs, index2 );
          MVT::MvAddMv( one/evecnrm[i], *evecstempr, zero, *evecstempr, *evecr1 );
          // Increment counter.
          i++;
        }
      }
      // Set the imaginary part of the eigenvectors if conjugate pairs exist.
      // If the last eigenvector has a split conjugate pair, don't set negative imaginary
      // part.
      if (conjprs) {        
        Teuchos::RefCountPtr<MV> evecstempi, eveci1;
        std::vector<int> indexi_pnev( IntOne );
        //
        // There is storage for an extra eigenvector.  
        // So, when the last eigenvalues is the first of a conjugate pair, that eigenvector will be computed.
        //
        for (i=0; i<conjprs; i++) {
          indexi_pnev[0] = indexi[i] + _nev;
          index2[0] = indexi[i];          
          t_evecnrm = one/lapack_mag.LAPY2(evecnrm[indexi[i]],evecnrm[indexi[i]-1]);
          evecstempi = MVT::CloneView( *evecstemp, index2 ); 

          if ( indexi_pnev[0] < 2*_nev ) {
            eveci1 = MVT::CloneView( *_evecs, indexi_pnev );
            MVT::MvAddMv( t_evecnrm*Teuchos::ScalarTraits<ScalarType>::magnitude((*_ritzvalues)[_totallength+indexi[i]])/(*_ritzvalues)[_totallength+indexi[i]],
                          *evecstempi, zero, *evecstempi, *eveci1 );
          }
          // Change index and set non-conjugate part of imag eigenvector.
          indexi[i]--;
          indexi_pnev[0]--;
          eveci1 = MVT::CloneView( *_evecs, indexi_pnev );
          MVT::MvAddMv( t_evecnrm*Teuchos::ScalarTraits<ScalarType>::magnitude((*_ritzvalues)[_totallength+indexi[i]])/(*_ritzvalues)[_totallength+indexi[i]],
                        *evecstempi, zero, *evecstempi, *eveci1 );
        }
      }
    }
    // Set flag to indicate that the eigenvectors are current.
    _isevecscurrent = true;
  }
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::ComputeSchurForm( const bool apply )
  {
    // local timer
    Teuchos::TimeMonitor LocalTimer(*_timerCompSF);

    int m = _jstart*_blockSize, n=_jstart*_blockSize;
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::BLAS<int,ScalarType> blas; 
    Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > Hj;
    Teuchos::SerialDenseMatrix<int,ScalarType> Q( n, n );

    if (apply) {
      Hj = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::View, _hessmatrix, m, n ) );
    } else {        
      Hj = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::Copy, _hessmatrix, m, n ) );
    }
    //
    SortSchurForm( *Hj, Q );
    if (_error_flg) {
      return;
    }
    //
    if (_nevblock <= _jstart || apply ) {
      //
      // Necessary variables.
      //
      int i=0,j=0       ;
      int mm1 = (_jstart-1)*_blockSize;
      int _nevtemp, numimag;
      //
      // Determine new offset depending upon placement of conjugate pairs.        
      // ( if we are restarting, determine the new offset ) 
      if (apply) {
        _offset = _maxoffset;
        for (i=0; i<_maxoffset; i++) {
          numimag = 0;
          for (j=0; j<(_nevblock+i)*_blockSize; j++) { 
            if ((*_ritzvalues)[_totallength+j]!=zero) { numimag++; }; 
          }
          if (!(numimag % 2)) { _offset = i; break; }
        }
      }
      _nevtemp = n;
      if (_jstart > _nevblock+_offset)
        _nevtemp = (_nevblock+_offset)*_blockSize;
      //
      Teuchos::SerialDenseMatrix<int,ScalarType> sub_block_hess(Teuchos::View, _hessmatrix, _blockSize, _blockSize, m, mm1);
      Teuchos::SerialDenseMatrix<int,ScalarType> sub_block_q(Teuchos::View, Q, _blockSize, _nevtemp, mm1 );
      Teuchos::SerialDenseMatrix<int,ScalarType> sub_block_b( _blockSize, _nevtemp );
      blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, _blockSize, _nevtemp, _blockSize, one, 
                 sub_block_hess.values(), sub_block_hess.stride(), sub_block_q.values(), 
                 sub_block_q.stride(), zero, sub_block_b.values(), _blockSize );
      //
      //---------------------------------------------------
      // Compute Schur decomposition error
      //
      // The residual for the Schur decomposition A(VQ) = (VQ)T + FB_m^TQ
      // where HQ = QT is || FB_m^TQ || = || H_{m+1,m}*B_m^TQ ||.
      //
      // We are only interested in the partial Krylov-Schur decomposition corresponding
      // to the _nev eigenvalues of interest or the _nevblock*_blockSize number of
      // eigenvalues we're keeping.
      // NOTE:  The Schur error is not updated if the Schur decomposition is
      //        not large enough to compute _nev eigenvalues, else we could accidently
      //        satisfy a condition for convergence.
      //---------------------------------------------------
      //
      if (_nevblock <= _jstart ) {
        //
        Teuchos::SerialDenseMatrix<int,ScalarType> sub_block_b2(Teuchos::View, sub_block_b, _blockSize, _nev);
        _schurerror = sub_block_b2.normFrobenius();
      }
      if (apply) {
        //
        //  We are going to restart, so update the Krylov-Schur decomposition.
        //
        // Update the Krylov-Schur basis.  
        //        
        std::vector<int> index(_nevtemp);
        for (i = 0; i < _nevtemp; i++ )
          index[i] = i;
        Teuchos::SerialDenseMatrix<int,ScalarType> Qnev(Teuchos::View, Q, n, _nevtemp);
        Teuchos::RefCountPtr<MV> basistemp = MVT::CloneView( *_basisvecs, index );
        index.resize( n );
        for (i = 0; i < n; i++ )
          index[i] = i;
        Teuchos::RefCountPtr<MV> basistemp2 = MVT::CloneCopy( *_basisvecs, index );        
        MVT::MvTimesMatAddMv ( one, *basistemp2, Qnev, zero, *basistemp );
        //
        // Update the Krylov-Schur quasi-triangular matrix.
        //
        Teuchos::SerialDenseMatrix<int,ScalarType> Hjp1(Teuchos::View, _hessmatrix, _blockSize, _nevtemp, _nevtemp );
        Hjp1.assign( sub_block_b );
      }
    }
  }
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::SortSchurForm( Teuchos::SerialDenseMatrix<int,ScalarType>& H, 
                                                          Teuchos::SerialDenseMatrix<int,ScalarType>& Q ) {
    // Local timer
    Teuchos::TimeMonitor LocalTimer(*_timerSortSF);

    ReturnType ret;
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    Teuchos::LAPACK<int,ScalarType> lapack; 
    Teuchos::LAPACK<int,MagnitudeType> lapack_mag;
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::BLAS<int,MagnitudeType> blas_mag;
    int i, j, info=0;
    int n = H.numRows(), ldh = H.stride(), ldq = Q.stride(); 
    int m = H.numRows(), mm1 = H.numRows() - _blockSize;
    ScalarType* ptr_h = H.values();
    ScalarType* ptr_q = Q.values();
    //
    //  If the operator is symmetric, analyze the block tridiagonal matrix
    //  and enforce symmetry.
    //
    if (_problem->IsSymmetric()) {
      if (_numRestarts > 0 && _restarts!=0) {
        //
        // The method has been restarted, so more caution must be used in
        // imposing symmetry.
        //
        for(j=_nevblock*_blockSize; j<n; j++) {
          for(i=0; i<j; i++) {
            H( i, j ) = H( j, i );
          }
        }
      } else {
        //
        // We haven't restarted, so just enforce symmetry throughout Hj.
        //
        for( j=0; j<n; j++ ) {
          for( i=0; i<j; i++ ) {
            H( i, j ) = H( j, i );
          }
        }
      }
    }
    //
    //---------------------------------------------------
    // Compute the current eigenvalue estimates
    // ---> Use driver GEES to first reduce to upper Hessenberg 
    //         form and then compute Schur form, outputting eigenvalues
    //---------------------------------------------------
    //
    // Allocate the work space. This space will be used below for calls to:
    // * GEES  (requires 3*N for real, 2*N for complex)
    // * TREVC (requires 3*N for real, 2*N for complex) 
    // * TREXC (requires N for real, none for complex)
    // Furthermore, both TREVC and GEES require a real array of length N (for
    // complex datatypes)
    int lwork = 3*n;
    std::vector<ScalarType> work( lwork );
    std::vector<MagnitudeType> rwork( n );
    std::vector<int> bwork( n );
    int sdim = 0; 
    char jobvs = 'V';
    lapack.GEES( jobvs, n, ptr_h, ldh, &sdim, &(*_ritzvalues)[0],
                 &(*_ritzvalues)[_totallength], ptr_q, ldq, &work[0], lwork, &rwork[0], &bwork[0], &info );
    assert(info==0);
    //
    //---------------------------------------------------
    // Compute the current Ritz residuals for ALL the eigenvalues estimates (Ritz values)
    //           || Ax - x\theta || = || FB_m^Ts || 
    //                              = || V_m+1*H_{m+1,m}*B_m^T*s ||
    //                              = || H_{m+1,m}*B_m^T*s ||
    //
    // where V_m is the current Krylov-Schur basis and x = V_m*s
    // NOTE: This means that s = e_i if the problem is symmetric, else the eigenvectors
    //       of the Schur form need to be computed.
    //
    // First compute H_{m+1,m}*B_m^T, then determine what 's' is.
    //---------------------------------------------------
    //
    // H_{m+1,m}
    Teuchos::SerialDenseMatrix<int,ScalarType> sub_block_hess(Teuchos::View, _hessmatrix, _blockSize, _blockSize, m, mm1);
    //
    // Last block rows of Q since the previous B_m is E_m (the last m-block of canonical basis vectors)
    Teuchos::SerialDenseMatrix<int,ScalarType> sub_block_q(Teuchos::View, Q, _blockSize, n, mm1 );
    //
    // Compute H_{m+1,m}*B_m^T
    Teuchos::SerialDenseMatrix<int,ScalarType> sub_block_b( _blockSize, n );
    blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, _blockSize, n, _blockSize, one, 
               sub_block_hess.values(), sub_block_hess.stride(), sub_block_q.values(), 
               sub_block_q.stride(), zero, sub_block_b.values(), _blockSize );
    //
    // Determine what 's' is and compute Ritz residuals.
    //
    ScalarType* b_ptr = sub_block_b.values();
    if (_problem->IsSymmetric()) {
      //
      // 's' is the i-th canonical basis vector.
      //
      for (i=0; i<n ; i++) {
        (*_ritzresiduals)[i] = blas.NRM2(_blockSize, b_ptr + i*_blockSize, 1);
      }   
    } else {
      //
      //  's' is the eigenvector of the block upper triangular, Schur matrix.
      //
      char side = 'R';
      int mm;
      const int ldvl = 1;
      ScalarType vl[ ldvl ];
      Teuchos::SerialDenseMatrix<int,ScalarType> Q_temp( n, n );
      Teuchos::SerialDenseMatrix<int,ScalarType> S( _blockSize, n );
      lapack.TREVC( side, n, H.values(), H.stride(), vl, ldvl,
                    Q_temp.values(), Q_temp.stride(), n, &mm, &work[0], &rwork[0], &info );
      assert(info==0);
      //
      // Scale the eigenvectors so that their euclidean norms are all one.
      // ( conjugate pairs get normalized by the sqrt(2) )
      //
      ScalarType temp;
      ScalarType* qt_ptr = Q_temp.values();
      i = 0;
      while( i < n ) {
        if ( (*_ritzvalues)[_totallength+i] != zero ) {
          temp = lapack_mag.LAPY2( blas.NRM2( n, qt_ptr+i*n, 1 ), blas.NRM2( n, qt_ptr+(i+1)*n, 1 ) );
          blas.SCAL( n, one/temp, qt_ptr+i*n, 1 );
          blas.SCAL( n, one/temp, qt_ptr+(i+1)*n, 1 );
          i = i+2;
        } else {
          temp = blas.NRM2( n, qt_ptr+i*n, 1 );
          blas.SCAL( n, one/temp, qt_ptr+i*n, 1 );
          i++;
        }
      }
      //
      // Compute H_{m+1,m}*B_m^T*S where the i-th column of S is 's' for the i-th Ritz-value
      //
      blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, _blockSize, n, n, one, 
                 sub_block_b.values(), sub_block_b.stride(), Q_temp.values(), 
                 Q_temp.stride(), zero, S.values(), S.stride() );
      ScalarType* s_ptr = S.values();
      i = 0;
      while( i < n ) {
        if ( (*_ritzvalues)[_totallength+i] != zero ) {
          (*_ritzresiduals)[i] = lapack_mag.LAPY2( blas.NRM2(_blockSize, s_ptr + i*_blockSize, 1),
                                            blas.NRM2(_blockSize, s_ptr + (i+1)*_blockSize, 1) );
          (*_ritzresiduals)[i+1] = (*_ritzresiduals)[i];
          i = i+2;
        } else {
          (*_ritzresiduals)[i] = blas.NRM2(_blockSize, s_ptr + i*_blockSize, 1);
          i++;
        }
      }
    }
    //
    //---------------------------------------------------
    // Sort the eigenvalues
    //---------------------------------------------------
    //
    if (_problem->IsSymmetric()) {
      if (_om->isVerbosityAndPrint( Debug )) {
        _os <<endl<<"Iteration ( "<< _iter <<" ) : [ Ritz Values Before Sorting (from SortManager) ]"<<endl;
        for(i=0; i<n; i++) {
          _os << (*_ritzvalues)[i] << endl;
        }
        _os <<endl;
      }
      //
      
      {
	Teuchos::TimeMonitor SortTimer(*_timerSortEval);
	ret = _sm->sort( this, n, &(*_ritzvalues)[0], &_order );
      }
      if (ret != Ok) {
        if (_om->isVerbosityAndPrint(Error)) {
          _os << "ERROR : Sorting in SortSchurForm"
              << endl;
        }
        _error_flg = true;
        return;
      }
      //
      if (_om->isVerbosityAndPrint( Debug )) {
        _os <<endl<<"Iteration ( "<< _iter <<" ) : [ Ritz Values After Sorting (from SortManager) ]"<<endl;
        for(i=0; i<n; i++)
          _os << (*_ritzvalues)[i] << "\t" << _order[i] << endl;
        _os <<endl;
      }
      //
    } else {
      if (_om->isVerbosityAndPrint( Debug )) {
        _os <<endl<<"Iteration ( "<< _iter <<" ) : [ Ritz Values Before Sorting (from SortManager) ]"<<endl;
        for(i=0; i<n; i++)
          _os << (*_ritzvalues)[i] << "\t" << (*_ritzvalues)[_totallength+i] << endl;
        _os <<endl;
      }
      //
      {
	Teuchos::TimeMonitor SortTimer(*_timerSortEval);
	ret = _sm->sort( this, n, &(*_ritzvalues)[0], &(*_ritzvalues)[_totallength], &_order );
      }
      if (ret != Ok) {
        if (_om->isVerbosityAndPrint(Error)) {
          _os << "Error sorting in SortSchurForm!"
              << endl;
        }
        _error_flg = true;
        return;
      }
      //
      if (_om->isVerbosityAndPrint( Debug )) {
        _os <<endl<<"Iteration ( "<< _iter <<" ) : [ Ritz Values After Sorting (from SortManager) ]"<<endl;
        for(i=0; i<n; i++)
          _os << (*_ritzvalues)[i] << "\t" << (*_ritzvalues)[_totallength+i] << "\t" << _order[i] << endl;
        _os <<endl;
      }
    }
    //
    // Re-sort _ritzresiduals based on _order
    //
    std::vector<MagnitudeType> ritz2( n );
    for (i=0; i<n; i++) { ritz2[i] = (*_ritzresiduals)[ _order[i] ]; }
    blas_mag.COPY( n, &ritz2[0], 1, &(*_ritzresiduals)[0], 1 );
    //
    // Copy the nev eigenvalues into the proper vectors
    // NOTE:  If we don't have nev Ritz values, then only n are copied
    //
    ( n > _nev ? blas.COPY( _nev, &(*_ritzvalues)[0], 1, &(*_evals)[0], 1 ) : blas.COPY( n, &(*_ritzvalues)[0], 1, &(*_evals)[0], 1 ) );
    if (!_problem->IsSymmetric() )
      ( n > _nev ? blas.COPY( _nev, &(*_ritzvalues)[_totallength], 1, &(*_evals)[_nev], 1 ) : blas.COPY( n, &(*_ritzvalues)[_totallength], 1, &(*_evals)[_nev], 1 ) );
    //
    //---------------------------------------------------
    // Reorder real Schur factorization, remember to add one to the indices for the
    // fortran call and determine offset.  The offset is necessary since the TREXC
    // method reorders in a nonsymmetric fashion, thus we use the reordering in
    // a stack-like fashion.  Also take into account conjugate pairs, which may mess
    // up the reordering, since the pair is moved if one of the pair is moved.
    //---------------------------------------------------
    //
    if (_om->isVerbosityAndPrint( Debug )) {
      _os <<endl<<"Iteration ( "<< _iter <<" ) : [ Before sorting the Schur form ( H ) ]"<<endl;
      H.print(_os);          
      _os <<endl;
    }
    //
    int _nevtemp = 0;
    char compq = 'V';
    std::vector<int> offset2( n );
    std::vector<int> order2( n );
    i = 0; 
    while (i < n) {
      if ((*_ritzvalues)[_totallength+i] != zero) {
        offset2[_nevtemp] = 0;
        for (j=i; j<n; j++) {
          if (_order[j] > _order[i]) { offset2[_nevtemp]++; }
        }
        order2[_nevtemp] = _order[i];
        i = i+2;
      } else {
        offset2[_nevtemp] = 0;
        for (j=i; j<n; j++) {
          if (_order[j] > _order[i]) { offset2[_nevtemp]++; }
        }
        order2[_nevtemp] = _order[i];
        i++;
      }
      _nevtemp++;
    }
    for (i=_nevtemp-1; i>=0; i--) {
      lapack.TREXC( compq, n, ptr_h, ldh, ptr_q, ldq, order2[i]+1+offset2[i], 
                    1, &work[0], &info );
      assert(info==0);
    }
    //
    if (_om->isVerbosityAndPrint( Debug )) {
      _os <<endl<<"Iteration ( "<< _iter <<" ) : [ After sorting and reordering the Schur form ( H ) ]"<<endl;
      H.print(_os);          
      _os <<endl;
    }
    //
    // Determine largest off diagonal element of Schur matrix for symmetric case.
    // (this is being disabled for the time being, since it's not being used)
    /*
    ScalarType _maxsymmelem = zero;
    if (_problem->IsSymmetric()) {
      for(j=0; j<n; j++){
        for(i=0; i<j; i++) {
          if(Teuchos::ScalarTraits<ScalarType>::magnitude(H(i, j))>_maxsymmelem) { _maxsymmelem = H(i, j); }
        }
      }
    }
    */
  }
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::Restart() {
    //  This method assumes the ComputeResiduals has been called before it
    //  to compute the Schur vectors and residuals.  This information is used to 
    //  restart the factorization.
    //
    int i;
    int _nevtemp = (_nevblock+_offset)*_blockSize;
    std::vector<int> index( _blockSize );
    //
    //  Move the F_vec block to the _jstart+1 position.        
    //
    for (i=0; i<_blockSize; i++) {
      index[i] = _jstart*_blockSize + i;
    }
    Teuchos::RefCountPtr<MV> F_vec = MVT::CloneCopy( *_basisvecs, index );
    for (i=0; i<_blockSize; i++) {
      index[i] = _nevtemp + i;
    }
    MVT::SetBlock( *F_vec, index, *_basisvecs );
    //
    //  Reset the pointer.
    //
    _jstart = _nevblock+_offset; 
  }
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::CheckSchurVecs( const int j ) {
    //
    // Check the difference between the projection of A with the Schur vectors and the Schur matrix.
    // 
    int i, n = j*_blockSize;
    ReturnType ret;
    std::vector<int> index( n );
    for( i=0; i<n; i++ ) { index[i] = i; } 
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    Teuchos::SerialDenseMatrix<int,ScalarType> Hj( Teuchos::View, _hessmatrix, n, n );
    Teuchos::SerialDenseMatrix<int,ScalarType> SchurProj( n, n );
    Teuchos::RefCountPtr<MV> Z = MVT::CloneView( *_basisvecs, index );
    Teuchos::RefCountPtr<MV> basistemp = MVT::Clone( *_basisvecs, n );
    {
      Teuchos::TimeMonitor OpTimer(*_timerOp);
    ret = OPT::Apply( *_Op, *Z, *basistemp );
    }
    _count_ApplyOp += MVT::GetNumberVecs(*Z);
    if (ret != Ok) {
      if (_om->isVerbosityAndPrint(Error)) {
        _os << "ERROR : Applying Op operator in CheckSchurVecs"
            << endl;
      }
      _error_flg = true;
      return;
    }
    MVT::MvTransMv( one, *Z, *basistemp, SchurProj );
    SchurProj.scale( -one );
    SchurProj += Hj;
    if (_om->doPrint()) {
      _os << endl 
          << "Error in Schur Projection ( || (VQ)^T*A*(VQ) - S || ) at restart " 
          << _numRestarts << " is "<< SchurProj.normFrobenius() << endl;
    }
  }
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  
  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::CheckBlkArnRed( const int j ) {
    int i,k,m=(j+1)*_blockSize;
    std::vector<int> index( _blockSize );
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    ReturnType ret;       
    
    for ( i=0; i<_blockSize; i++ ) {
      index[i] = m+i;
    }
    Teuchos::RefCountPtr<MV> F_vec = MVT::CloneView( *_basisvecs, index );
    
    std::vector<MagnitudeType> ptr_norms( m );
    ScalarType sum = zero;
    
    MVT::MvNorm( *F_vec, &ptr_norms );
    for ( i=0; i<_blockSize; i++ ) {
      sum += ptr_norms[i];
    }
    
    index.resize( m );
    for ( i=0; i<m; i++ ) {
      index[i] = i;  
    }
    Teuchos::RefCountPtr<MV> Vj = MVT::CloneView( *_basisvecs, index );

    if (_om->doPrint()) {
      _os << endl;
      _os << "********Block Arnoldi Iteration "<<j+1<<"******** " << endl;
      _os << endl;
    }    

    Teuchos::SerialDenseMatrix<int,ScalarType> VTV(m,m);
    _problem->InnerProd( *Vj, *Vj, VTV );
    ScalarType* ptr=VTV.values();
    ScalarType column_sum;
    
    for (k=0; k<m; k++) {
      column_sum=zero;
      for (i=0; i<m; i++) {
        if (i==k) {
          ptr[i] -= one;
        }
        column_sum += ptr[i];
      }

      if (_om->doPrint()) {
        _os << endl 
            <<" V^T*B*V-I " << "for column " << k << " is " 
            << Teuchos::ScalarTraits<ScalarType>::magnitude(column_sum) 
            << endl;
      }
      ptr += m;
    }
    
    Teuchos::SerialDenseMatrix<int,ScalarType> E(m,_blockSize);
    
    _problem->InnerProd( *Vj, *F_vec, E );
    ScalarType* ptr_Ej=E.values();
    
    for (k=0;k<_blockSize;k++) {
      column_sum=zero;
      for (i=0; i<m; i++) {
        column_sum += ptr_Ej[i];
      }
      ptr_Ej += m;
      if (ptr_norms[k]) column_sum = column_sum/ptr_norms[k];

      if (_om->doPrint()) {
        _os << endl 
            << " B-Orthogonality with F " << "for column " << k 
            << " is " << Teuchos::ScalarTraits<ScalarType>::magnitude(column_sum) 
            << endl;
      }
    }
    
    Teuchos::RefCountPtr<MV> AVj = MVT::Clone( *_basisvecs, m ); 
    {
      Teuchos::TimeMonitor OpTimer(*_timerOp);
      ret = OPT::Apply( *_Op, *Vj, *AVj );
    }
    _count_ApplyOp += MVT::GetNumberVecs( *Vj );
    if (ret != Ok) {
      // Apply() failed. Return error code.
      if (_om->isVerbosityAndPrint(Error)) {
        _os << "Error applying Op operator in CheckBlkArnRed!"
            << endl;
      }
      _error_flg = true;
      return;
    }
    Teuchos::SerialDenseMatrix<int,ScalarType> Hj(Teuchos::View, _hessmatrix, m, m);
    MVT::MvTimesMatAddMv( -one, *Vj, Hj, one, *AVj );
    index.resize( _blockSize );
    for ( i=0; i<_blockSize; i++ ) {  
      index[i] = j*_blockSize+i;
    }
    
    Teuchos::RefCountPtr<MV> Fj = MVT::CloneView( *AVj, index );
    MVT::MvAddMv(-one, *F_vec, one, *Fj, *Fj);
    
    MVT::MvNorm( *AVj, &ptr_norms );
    
    for ( i=0; i<m; i++ ) { 
      if (_om->doPrint()) {
        _os << endl
            << " Arnoldi relation " << "for column " << i << " is " << ptr_norms[i] 
            << endl;  
      }
    }
  }

} // End of namespace Anasazi

#endif

// End of file AnasaziBlockKrylovSchur.hpp



