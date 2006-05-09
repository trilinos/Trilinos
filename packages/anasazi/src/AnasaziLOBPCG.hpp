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

/*! \file AnasaziLOBPCG.hpp
  \brief Implementation of the locally-optimal block preconditioned conjugate gradient method
*/

#ifndef ANASAZI_LOBPCG_HPP
#define ANASAZI_LOBPCG_HPP

#include "AnasaziEigensolver.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziSortManager.hpp"
#include "AnasaziModalSolverUtils.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!        \class Anasazi::LOBPCG

        \brief This class implements the Locally-Optimal Block Preconditioned
        Conjugate Gradient method for solving symmetric eigenvalue problems.

        This implementation is a modification of the one found in :
        A. Knyazev, "Toward the optimal preconditioned eigensolver:
        Locally optimal block preconditioner conjugate gradient method",
        SIAM J. Sci. Comput., vol 23, n 2, pp. 517-541

        \author Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {

  template <class ScalarType, class MV, class OP>
  class LOBPCG : public Eigensolver<ScalarType,MV,OP> { 
  public:
    //@{ \name Constructor/Destructor.

    //! %Anasazi::LOBPCG constructor.
    LOBPCG( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
            const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sm,
            const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
            Teuchos::ParameterList &pl 
            );

    //! %Anasazi::LOBPCG destructor.
    virtual ~LOBPCG() {};
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
    
    //! Get the solvers native residuals for the current eigenpairs. 
    /*! This is not be the same as the true residuals for most solvers. Sometimes the native
      residuals are not in multivector form, so the norm type is solver dependent.  
      
      \note
      <ol>
      <li> If the native residual is in multivector form then a non-null pointer will be
      returned, else the normvec will be populated with the current residual norms. 
      <li> If the native residual is returned in multivector form, the memory is managed
      by the calling routine.
      </ol>
    */
    Teuchos::RefCountPtr<const MV> GetNativeResiduals( ScalarType* normvec ) const { return Teuchos::null; };
    
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
    //
    // Internal methods
    //
    void accuracyCheck(const MV *X, const MV *MX,
                       const MV *R, const MV *Q,
                       const MV *H, const MV *P) const;
    //
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem;
    const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > _sm;
    const Teuchos::RefCountPtr<OutputManager<ScalarType> > _om;
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
    int _numRestarts, _iter, _knownEV, _nevLocal;
    std::vector<MagnitudeType> _theta;
    std::vector<MagnitudeType> _normR, _resids;
    bool _error_flg;
    //
    // Internal utilities class required by eigensolver.
    //
    ModalSolverUtils<ScalarType,MV,OP> _MSUtils;
    //
    // Output stream from the output manager
    //
    std::ostream& _os;
    //
    // Internal timers
    //
    bool _restartTimers;
    Teuchos::RefCountPtr<Teuchos::Time> _timerOp, _timerMOp, _timerPrec,
                                        _timerSortEval, 
                                        _timerLocalProj, _timerDS,
                                        _timerLocalUpdate, _timerCompRes,
                                        _timerOrtho, 
                                        _timerRestart, _timerTotal;
    //
    // Counters
    //
    int _count_ApplyOp, _count_ApplyM, _count_ApplyPrec;
  };
  //
  // Implementation
  //
  template <class ScalarType, class MV, class OP>
  LOBPCG<ScalarType,MV,OP>::LOBPCG(const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
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
    _Prec(_problem->GetPrec()),
    _evecs(_problem->GetEvecs()), 
    _evals(problem->GetEvals()), 
    _nev(problem->GetNEV()), 
    _maxIter(_pl.get("Max Iters", 300)),
    _blockSize(_pl.get("Block Size", 1)),
    _residual_tolerance(_pl.get("Tol", (MagnitudeType)(1.0e-6) )),
    _numRestarts(0), 
    _iter(0), 
    _knownEV(0),
    _nevLocal(0),
    _error_flg(false),
    _MSUtils(om),
    _os(_om->GetOStream()),
    _restartTimers(_pl.get("Restart Timers",false)),
    _timerOp(Teuchos::TimeMonitor::getNewTimer("Operation Op*x")),
    _timerMOp(Teuchos::TimeMonitor::getNewTimer("Operation M*x")),
    _timerPrec(Teuchos::TimeMonitor::getNewTimer("Operation Prec*x")),
    _timerSortEval(Teuchos::TimeMonitor::getNewTimer("Sorting eigenvalues")),
    _timerLocalProj(Teuchos::TimeMonitor::getNewTimer("Local projection")),
    _timerDS(Teuchos::TimeMonitor::getNewTimer("Direct solve")),
    _timerLocalUpdate(Teuchos::TimeMonitor::getNewTimer("Local update")),
    _timerCompRes(Teuchos::TimeMonitor::getNewTimer("Computing residuals")),
    _timerOrtho(Teuchos::TimeMonitor::getNewTimer("Orthogonalization")),
    _timerRestart(Teuchos::TimeMonitor::getNewTimer("Restarting")),
    _timerTotal(Teuchos::TimeMonitor::getNewTimer("Total time")),
    _count_ApplyOp(0),
    _count_ApplyM(0),
    _count_ApplyPrec(0)
  {     
  }

  template <class ScalarType, class MV, class OP>
  void 
  LOBPCG<ScalarType,MV,OP>::currentStatus() 
  {
    int i;
    if (_om->doPrint()) {
      cout.setf(ios::scientific, ios::floatfield);  
      cout.precision(6);
      _os <<endl;
      _os <<"******************* CURRENT STATUS *******************"<<endl;
      _os <<"The number of iterations performed is " <<_iter<<endl;
      _os <<"The number of restarts performed is "<<_numRestarts<<endl;
      _os <<"The current block size is "<<_blockSize<<endl;
      _os <<"The number of eigenvalues requested is "<<_nev<<endl;
      _os <<"The number of eigenvalues computed is "<<_knownEV<<endl;
      _os <<"The requested residual tolerance is "<<_residual_tolerance<<endl;
      _os <<"The number of operations Op*x   is "<<_count_ApplyOp<<endl;
      _os <<"The number of operations M*x    is "<<_count_ApplyM<<endl;
      _os <<"The number of operations Prec*x is "<<_count_ApplyPrec<<endl;
      _os << endl;
      _os <<"COMPUTED EIGENVALUES                 "<<endl;
      _os.setf(ios_base::right, ios_base::adjustfield);
      _os << std::setw(16) << "Eigenvalue" 
          << std::setw(16) << "Ritz Residual"
          << endl;
      _os <<"------------------------------------------------------"<<endl;
      if ( _knownEV > 0 ) {
        for (i=0; i<_knownEV; i++) {
          _os << std::setw(16) << (*_evals)[i] 
              << std::setw(16) << _resids[i] 
              << endl;
        }
      } 
      else {
        _os <<"[none computed]"<<endl;
      }
      _os <<endl;
      _os <<"CURRENT EIGENVALUE ESTIMATES             "<<endl;
      _os << std::setw(16) << "Ritz value" 
          << std::setw(16) << "Residual"
          << endl;
      _os <<"------------------------------------------------------"<<endl;
      if ( _iter > 0 || _nevLocal > 0 ) {
        for (i=0; i<_blockSize; i++) {
          _os << std::setw(16) << _theta[i] 
              << std::setw(16) << _normR[i] 
              << endl;
        }
      } 
      else {
        _os <<"[none computed]"<<endl;
      }
      _os << "******************************************************"<<endl;  
      _os << endl;
    }
  }

  template <class ScalarType, class MV, class OP>
  ReturnType LOBPCG<ScalarType,MV,OP>::solve () 
  {
    Teuchos::TimeMonitor LocalTimer(*_timerTotal,_restartTimers);

    if ( _restartTimers ) {
      _timerOp->reset();
      _timerMOp->reset();
      _timerPrec->reset();
      _timerSortEval->reset();
      _timerLocalProj->reset();
      _timerDS->reset();
      _timerLocalUpdate->reset();
      _timerCompRes->reset();
      _timerRestart->reset();
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
    // Check that the maximum number of iterations is a positive number
    //    
    if ( _maxIter<=0 ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error )) 
        _os << "ERROR : maxIter = "<< _maxIter <<" [ should be positive number ] " << endl;
      return Failed;
    } 
    //
    // If the search subspace dimension is larger than the dimension of the operator, reset
    // the maximum number of blocks accordingly.
    //    
    if ( _blockSize > dim ) {
      if (_om->isVerbosityAndPrint( Anasazi::Error ))
        _os << "ERROR : Search space dimension (blockSize) = "<< _blockSize 
            <<" [ should not be greater than " << dim << " ] " << endl;
      return Failed;
    }
    //
    // Reinitialize internal data and pointers, prepare for solve
    //
    _numRestarts = 0; 
    _iter = 0; 
    _knownEV = 0;
    _error_flg = false;
    //
    // Necessary variables
    //
    int i, j;
    int info = 0, leftOver;
    bool reStart = false;
    bool criticalExit = false;
    bool haveMass = (_MOp.get()!=0);
    ReturnType ret;
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::LAPACK<int,ScalarType> lapack;
    std::vector<int> _order;
    //
    // Internal storage for eigensolver
    //
    // Working vectors
    Teuchos::RefCountPtr<MV> X, MX, KX;
    Teuchos::RefCountPtr<MV> R;
    //
    X = MVT::Clone( *iVec, _blockSize );
    KX = MVT::Clone( *iVec, _blockSize );
    R = MVT::Clone( *iVec, _blockSize );
    //
    // Initialize the workspace.
    //
    {
      int numIVecs = MVT::GetNumberVecs( *iVec );
      if (numIVecs > _blockSize)
        numIVecs = _blockSize;
      std::vector<int> index( numIVecs );
      for (i=0; i<numIVecs; i++) {
        index[i] = i;
      }
      //
      // Copy the first numIVecs of the initial vectors into the first
      // numIVecs X (any additional vectors in iVec are ignored)
      //
      MVT::SetBlock( *iVec, index, *X );
      //
      // Augment the initial vectors with random vectors if necessary
      //
      int leftOver = _blockSize - numIVecs;
      if (leftOver) {
        index.resize(leftOver);
        for (i=0; i<leftOver; i++)
          index[i] = numIVecs + i;
        Teuchos::RefCountPtr<MV> tmpIVec = MVT::CloneView( *X, index );
        MVT::MvRandom( *tmpIVec );
      }
    }

    //
    // Check to see if there is a mass matrix, so we know how much space is required.
    // [ If there isn't a mass matrix we can use a view of X.]
    //
    if (haveMass) {
      MX = MVT::Clone( *iVec, _blockSize );
    }
    else {
      std::vector<int> index( _blockSize );
      for (int i=0; i<_blockSize; i++)
        index[i] = i;
      MX = MVT::CloneView( *X, index );
    }
    //
    // Preconditioned search space working vectors.
    Teuchos::RefCountPtr<MV> H, MH, KH;
    //
    H = MVT::Clone( *iVec, _blockSize );
    KH = MVT::Clone( *iVec, _blockSize );
    //
    if (haveMass) {
      MH = MVT::Clone( *iVec, _blockSize );
    }
    else {
      std::vector<int> index( _blockSize );
      for (int i=0; i<_blockSize; i++)
        index[i] = i;
      MH = MVT::CloneView( *H, index );
    }
    //
    // Search direction working vectors.
    Teuchos::RefCountPtr<MV> P, MP, KP;
    //
    P = MVT::Clone( *iVec, _blockSize );
    KP = MVT::Clone( *iVec, _blockSize );
    //
    if (haveMass) {
      MP = MVT::Clone( *iVec, _blockSize );
    }
    else {
      std::vector<int> index( _blockSize );
      for (int i=0; i<_blockSize; i++)
        index[i] = i;
      MP = MVT::CloneView( *P, index );
    }
    //
    // theta = Storage for local eigenvalues     (size: 3*_blockSize)
    // normR = Storage for the norm of residuals (size: blockSize)
    // resids = Storage for the residuals of the converged eigenvalues (size: _nev)
    //
    _theta.resize( 3*_blockSize );
    _normR.resize( _blockSize );
    _resids.resize( _nev );

    // Define dense local matrices and arrays
    //
    // KK = Local stiffness matrix               (size: 3*_blockSize x 3*_blockSize)
    //
    // MM = Local mass matrix                    (size: 3*_blockSize x 3*_blockSize)
    //
    // S = Local eigenvectors                    (size: 3*_blockSize x 3*_blockSize)
    //
    Teuchos::SerialDenseMatrix<int,ScalarType> KK, MM, S;
    //
    KK.shape( 3*_blockSize, 3*_blockSize );
    MM.shape( 3*_blockSize, 3*_blockSize );
    S.shape( 3*_blockSize, 3*_blockSize );        
    //
    // Miscellaneous definitions
    int localSize = 0;
    const int twoBlocks = 2*_blockSize;
    const int threeBlocks = 3*_blockSize;
    int _nFound = _blockSize;
    //    
    for ( _iter=0; _iter <= _maxIter; _iter++ ) {
      
      if ( (_iter==0) || (reStart == true) ) {
        
        reStart = false;
        localSize = _blockSize;
        
        if (_nFound > 0) {
          
          std::vector<int> index( _nFound );
          for (i=0; i<_nFound; i++) {
            index[i] = _blockSize-_nFound + i;
          }
          Teuchos::RefCountPtr<MV> X2 = MVT::CloneView( *X, index );
          Teuchos::RefCountPtr<MV> MX2 = MVT::CloneView( *MX, index );
          Teuchos::RefCountPtr<MV> KX2 = MVT::CloneView( *KX, index );
          
          // Apply the mass matrix to X
          if (haveMass) {
	    {
	      Teuchos::TimeMonitor MOpTimer( *_timerMOp );
	      ret = OPT::Apply( *_MOp, *X2, *MX2 );
	    }
            _count_ApplyM += MVT::GetNumberVecs( *X2 );
            if (ret != Ok) {
              if (_om->isVerbosityAndPrint(Error)) {
                _os << "ERROR : Applying M operator in BlockReduction" << endl;
              }
              _error_flg = true;
              break; // break out of for(_iter) loop
            }
          }
          
          if (_knownEV > 0) {
            
            // Orthonormalize X against the known eigenvectors with Gram-Schmidt
            // Note: Use R as a temporary work space
            index.resize( _knownEV );
            for (i=0; i<_knownEV; i++) {
              index[i] = i;
            }
            Teuchos::RefCountPtr<MV> copyQ = MVT::CloneView( *_evecs, index );
            
	    {
	      Teuchos::TimeMonitor OrthoTimer( *_timerOrtho );
	      info = _MSUtils.massOrthonormalize( *X, *MX, _MOp.get(), *copyQ, _nFound, 0 );
	    }
            
            // Exit the code if the orthogonalization did not succeed
            if (info < 0) {
              info = -10;
              //return info;
            }
          } // if (knownEV > 0) 
          
          // Apply the stiffness matrix to X
	  {
	    Teuchos::TimeMonitor OpTimer( *_timerOp );
	    ret = OPT::Apply( *_Op, *X2, *KX2 );
	  }
          _count_ApplyOp += MVT::GetNumberVecs( *X2 );
          if (ret != Ok) {
            if (_om->isVerbosityAndPrint(Error)) {
              _os << "ERROR : Applying Op operator in BlockReduction" << endl;
            }
            _error_flg = true;
            break; // break out of for(_iter) loop
          }
          
        } // if (_nFound > 0)
        
      } // if ( (_iter==0) || reStart == true)

      else {
        
        // Apply the preconditioner on the residuals
        if (_Prec.get()) {
	  {
	    Teuchos::TimeMonitor PrecTimer( *_timerPrec );
	    ret = OPT::Apply( *_Prec, *R, *H );
	  }
          _count_ApplyPrec += MVT::GetNumberVecs( *R );
          if (ret != Ok) {
            if (_om->isVerbosityAndPrint(Error)) {
              _os << "ERROR : Applying Prec operator in BlockReduction" << endl;
            }
            _error_flg = true;
            break; // break out of for(_iter) loop
          }
        }
        else {
          MVT::MvAddMv( one, *R, zero, *R, *H );
        }

        // Apply the mass matrix on H
        if (haveMass) {
	  {
	    Teuchos::TimeMonitor MOpTimer( *_timerMOp );
	    ret = OPT::Apply( *_MOp, *H, *MH);
	  }
          _count_ApplyM += MVT::GetNumberVecs( *H );
          if (ret != Ok) {
            if (_om->isVerbosityAndPrint(Error)) {
              _os << "ERROR : Applying M operator in BlockReduction" << endl;
            }
            _error_flg = true;
            break; // break out of for(_iter) loop
          }
        }
        
        if (_knownEV > 0) {

          // Orthogonalize H against the known eigenvectors
          // Note: Use R as a temporary work space
          std::vector<int> index( _knownEV );
          for (i=0; i<_knownEV; i++) {
            index[i] = i;
          }
          Teuchos::RefCountPtr<MV> copyQ = MVT::CloneView( *_evecs, index );

	  {
	    Teuchos::TimeMonitor OrthoTimer( *_timerOrtho );
	    _MSUtils.massOrthonormalize( *H, *MH, _MOp.get(), *copyQ, _blockSize, 1);
	  }
        } // if (knownEV > 0)
        
        // Apply the stiffness matrix to H
	{
	  Teuchos::TimeMonitor OpTimer( *_timerOp );
	  ret = OPT::Apply( *_Op, *H, *KH);
	}
        _count_ApplyOp += MVT::GetNumberVecs( *H );
        if (ret != Ok) {
          if (_om->isVerbosityAndPrint(Error)) {
            _os << "ERROR : Applying Op operator in BlockReduction" << endl;
          }
          _error_flg = true;
          break; // break out of for(_iter) loop
        }
        
        if (localSize == _blockSize)
          localSize += _blockSize;

      } // if ( (_iter==0) || (reStart==true))

      // Form "local" mass and stiffness matrices
      {
	Teuchos::TimeMonitor LocalProjTimer( *_timerLocalProj );
	Teuchos::SerialDenseMatrix<int,ScalarType> KK11( Teuchos::View, KK, _blockSize, _blockSize );
	MVT::MvTransMv( one, *X, *KX, KK11 );
	
	Teuchos::SerialDenseMatrix<int,ScalarType> MM11( Teuchos::View, MM, _blockSize, _blockSize );
	MVT::MvTransMv( one, *X, *MX, MM11 );
      
	if (localSize > _blockSize) {
	  Teuchos::SerialDenseMatrix<int,ScalarType> KK12( Teuchos::View, KK, _blockSize, _blockSize, 
							   0, _blockSize );
	  MVT::MvTransMv( one, *X, *KH, KK12 );
	  
	  Teuchos::SerialDenseMatrix<int,ScalarType> KK22( Teuchos::View, KK, _blockSize, _blockSize, 
							   _blockSize, _blockSize );
	  MVT::MvTransMv( one, *H, *KH, KK22 );
	  
	  Teuchos::SerialDenseMatrix<int,ScalarType> MM12( Teuchos::View, MM, _blockSize, _blockSize, 
							   0, _blockSize );
	  MVT::MvTransMv( one, *X, *MH, MM12 );
	  
	  Teuchos::SerialDenseMatrix<int,ScalarType> MM22( Teuchos::View, MM, _blockSize, _blockSize, 
							   _blockSize, _blockSize );
	  MVT::MvTransMv( one, *H, *MH, MM22 );
	  
	  if (localSize > twoBlocks) {
	    
	    Teuchos::SerialDenseMatrix<int,ScalarType> KK13( Teuchos::View, KK, _blockSize, _blockSize, 
							     0, twoBlocks );
	    MVT::MvTransMv( one, *X, *KP, KK13 );
	    
	    Teuchos::SerialDenseMatrix<int,ScalarType> KK23( Teuchos::View, KK, _blockSize, _blockSize, 
							     _blockSize, twoBlocks );
	    MVT::MvTransMv( one, *H, *KP, KK23 );
	    
	    Teuchos::SerialDenseMatrix<int,ScalarType> KK33( Teuchos::View, KK, _blockSize, _blockSize, 
							     twoBlocks, twoBlocks );
	    MVT::MvTransMv( one, *P, *KP, KK33 );
	    
	    Teuchos::SerialDenseMatrix<int,ScalarType> MM13( Teuchos::View, MM, _blockSize, _blockSize, 
							     0, twoBlocks );
	    MVT::MvTransMv( one, *X, *MP, MM13 );
	    
	    Teuchos::SerialDenseMatrix<int,ScalarType> MM23( Teuchos::View, MM, _blockSize, _blockSize, 
							     _blockSize, twoBlocks );
	    MVT::MvTransMv( one, *H, *MP, MM23 );
	    
	    Teuchos::SerialDenseMatrix<int,ScalarType> MM33( Teuchos::View, MM, _blockSize, _blockSize, 
							     twoBlocks, twoBlocks );
	    MVT::MvTransMv( one, *P, *MP, MM33 );
	    
	  } // if (localSize > twoBlocks)
	  
	} // if (localSize > blockSize)

      } // end timing block

      // Perform a spectral decomposition
      _nevLocal = localSize;
      {
	Teuchos::TimeMonitor DSTimer( *_timerDS );
	info = _MSUtils.directSolver(localSize, KK, &MM, &S, &_theta, &_nevLocal, 
				     (_blockSize == 1) ? 1 : 0);
      }
      
      if (info < 0) {
        // Stop when spectral decomposition has a critical failure
        criticalExit = true;
        break;
      } // if (info < 0)

      // Check for restarting
      if (_nevLocal < _blockSize) {
        if (   _om->isVerbosityAndPrint( IterationDetails ) 
             ||_om->isVerbosityAndPrint( Debug )            ) {
          _os << " Iteration " << _iter;
          _os << "- Failure for spectral decomposition - RESTART with new random search\n";
        }
        if (_blockSize == 1) {
          MVT::MvRandom( *X );
          _nFound = _blockSize;
        }
        else {
          std::vector<int> index( _blockSize-1 );
          for (i=0; i<_blockSize-1; i++)
            index[i] = i+1;
          Teuchos::RefCountPtr<MV> Xinit = MVT::CloneView( *X, index );
          MVT::MvRandom( *Xinit );
          _nFound = _blockSize - 1;
        } // if (_blockSize == 1)
        reStart = true;
        _numRestarts += 1;
        info = 0;
        continue;
      } // if (_nevLocal < _blockSize)
    

      // We can reduce the size of the local problem if the directSolver detects rank deficiency 
      if (_nevLocal == _blockSize) {
        localSize = _blockSize;
      }
      else if (_nevLocal <= twoBlocks) {
        localSize = twoBlocks;
      }

      //
      //---------------------------------------------------
      // Sort the ritz values using the sort manager
      //---------------------------------------------------
      // The sort manager is templated on ScalarType
      // Make a ScalarType copy of _theta for sorting
      std::vector<ScalarType> _theta_st(_theta.size());
      {
	Teuchos::TimeMonitor SortTimer( *_timerSortEval );
	std::copy<MTiter,STiter>(_theta.begin(),_theta.begin()+_nevLocal,_theta_st.begin());
	_order.resize(_nevLocal);
	ret = _sm->sort( this, _nevLocal, &(_theta_st[0]), &_order );

	//  Reorder _theta according to sorting results from _theta_st
	std::vector<MagnitudeType> _theta_copy(_theta);
	for (i=0; i<_nevLocal; i++) {
	  _theta[i] = _theta_copy[_order[i]];
	}
      }
      if (ret != Ok) {
        if (_om->isVerbosityAndPrint(Error)) {
          _os << "ERROR : Sorting in solve()" << endl;
        }
        _error_flg = true;
        break;
      }
      // Sort the primitive ritz vectors
      // We need the first _blockSize vectors ordered to generate the next
      // columns immediately below, as well as later, when/if we restart.
      Teuchos::SerialDenseMatrix<int,ScalarType> copyS( S );
      for (i=0; i<_nevLocal; i++) {
        blas.COPY(_nevLocal, copyS[_order[i]], 1, S[i], 1);
      }
      
      // Compute the residuals
      Teuchos::SerialDenseMatrix<int,ScalarType> S11( Teuchos::View, S, _blockSize, _blockSize );
      {
	Teuchos::TimeMonitor CompResTimer( *_timerCompRes );
	MVT::MvTimesMatAddMv( one, *KX, S11, zero, *R );
	
	if (localSize >= twoBlocks) {
	  Teuchos::SerialDenseMatrix<int,ScalarType> S21( Teuchos::View, S, _blockSize, _blockSize, _blockSize );
	  MVT::MvTimesMatAddMv( one, *KH, S21, one, *R );
	  
	  if (localSize == threeBlocks) {
	    Teuchos::SerialDenseMatrix<int,ScalarType> S31( Teuchos::View, S, _blockSize, _blockSize, twoBlocks  );          
	    MVT::MvTimesMatAddMv( one, *KP, S31, one, *R );
	    
	  } // if (localSize == threeBlocks)
	} // if (localSize >= twoBlocks )
	
	// Replace S with S*Lambda
	for (j = 0; j < _blockSize; ++j) {
	  blas.SCAL(localSize, _theta[j], S[j], 1);
	}
	
	MVT::MvTimesMatAddMv( -one, *MX, S11, one, *R );
	
	if (localSize >= twoBlocks) {
	  Teuchos::SerialDenseMatrix<int,ScalarType> S21( Teuchos::View, S, _blockSize, _blockSize, _blockSize );
	  MVT::MvTimesMatAddMv( -one, *MH, S21, one, *R );
	  
	  if (localSize == threeBlocks) {
	    Teuchos::SerialDenseMatrix<int,ScalarType> S31( Teuchos::View, S, _blockSize, _blockSize, twoBlocks  );          
	    MVT::MvTimesMatAddMv( -one, *MP, S31, one, *R );
	  }
	} // if (localSize >= twoBlocks)
	
	// Restore S from S*Lambda back to S
	for (j = 0; j < _blockSize; ++j) {
	  blas.SCAL(localSize, one/_theta[j], S[j], 1);
	}

      } // end timing block

      // Compute the norms of the residuals
      MVT::MvNorm( *R, &_normR );
      
      // Scale the norms of residuals with the eigenvalues
      // Count the converged eigenvectors
      _nFound = 0;
      for (j = 0; j < _blockSize; ++j) {
        // Scale eigenvalues if _theta is non-zero.
        if ( _theta[j] != zero ) {
          _normR[j] /= SCT::magnitude(_theta[j]);
        }
        // Check for convergence
        if (_normR[j] < _residual_tolerance) {
          _nFound ++;
        }
      }
      
      // Store the residual history
      /*      if (localVerbose > 2) {
              memcpy(resHistory + historyCount*blockSize, normR, blockSize*sizeof(double));
              historyCount += 1;
              }
      */      
      
      // Print information on current iteration
      if (_om->isVerbosity( IterationDetails )) {
        currentStatus();
      }
      
      if (_nFound == 0) {
        // Update the spaces
        // Note: Use R as a temporary work space
        // Note: S11 was previously defined above

	{
	  Teuchos::TimeMonitor LocalUpdateTimer( *_timerLocalUpdate );
	  
	  MVT::MvAddMv( one, *X, zero, *X, *R );        
	  MVT::MvTimesMatAddMv( one, *R, S11, zero, *X );
	  
	  MVT::MvAddMv( one, *KX, zero, *KX, *R );
	  MVT::MvTimesMatAddMv( one, *R, S11, zero, *KX );
	  
	  if (haveMass) {
	    MVT::MvAddMv( one, *MX, zero, *MX, *R );
	    MVT::MvTimesMatAddMv( one, *R, S11, zero, *MX );
	  }
	  
	  if (localSize == twoBlocks) {
	    Teuchos::SerialDenseMatrix<int,ScalarType> S21( Teuchos::View, S, _blockSize, _blockSize, _blockSize );
	    MVT::MvTimesMatAddMv( one, *H, S21, zero, *P );
	    MVT::MvAddMv( one, *P, one, *X, *X );
	    
	    MVT::MvTimesMatAddMv( one, *KH, S21, zero, *KP );
	    MVT::MvAddMv( one, *KP, one, *KX, *KX );
	    
	    if (haveMass) {
	      MVT::MvTimesMatAddMv( one, *MH, S21, zero, *MP );
	      MVT::MvAddMv( one, *MP, one, *MX, *MX );
	    }
	  } // if (localSize == twoBlocks)
	  else if (localSize == threeBlocks) {          
	    Teuchos::SerialDenseMatrix<int,ScalarType> S21( Teuchos::View, S, _blockSize, _blockSize, _blockSize );          
	    Teuchos::SerialDenseMatrix<int,ScalarType> S31( Teuchos::View, S, _blockSize, _blockSize, twoBlocks  );          
	    MVT::MvAddMv( one, *P, zero, *P, *R );
	    MVT::MvTimesMatAddMv( one, *R, S31, zero, *P );
	    MVT::MvTimesMatAddMv( one, *H, S21, one, *P );
	    MVT::MvAddMv( one, *P, one, *X, *X );
	    
	    MVT::MvAddMv( one, *KP, zero, *KP, *R );
	    MVT::MvTimesMatAddMv( one, *R, S31, zero, *KP );
	    MVT::MvTimesMatAddMv( one, *KH, S21, one, *KP );
	    MVT::MvAddMv( one, *KP, one, *KX, *KX );
	    
	    if (haveMass) {
	      MVT::MvAddMv( one, *MP, zero, *MP, *R );
	      MVT::MvTimesMatAddMv( one, *R, S31, zero, *MP );
	      MVT::MvTimesMatAddMv( one, *MH, S21, one, *MP );
	      MVT::MvAddMv( one, *MP, one, *MX, *MX );
	    }
	  } // if (localSize == threeBlocks)
	  
	} // end timing block
	
        // Compute the new residuals
	{
	  Teuchos::TimeMonitor CompResTimer( *_timerCompRes );
	  MVT::MvAddMv( one, *KX, zero, *KX, *R );
	  Teuchos::SerialDenseMatrix<int,ScalarType> T( _blockSize, _blockSize );
	  for (j = 0; j < _blockSize; ++j) {
	    T(j,j) = _theta[j];
	  }
	  MVT::MvTimesMatAddMv( -one, *MX, T, one, *R );
	}
	
        // When required, monitor some orthogonalities
        if (_om->isVerbosity( OrthoDetails )) {
          if (_knownEV == 0) {
            accuracyCheck(X.get(), MX.get(), R.get(), 0, 0, 0);
          }
          else {
            std::vector<int> index2( _knownEV );
            for (j=0; j<_knownEV; j++)
              index2[j] = j;
            Teuchos::RefCountPtr<MV> copyQ = MVT::CloneView( *_evecs, index2 );
            accuracyCheck(X.get(), MX.get(), R.get(), copyQ.get(), (localSize>_blockSize) ? H.get() : 0,
                          (localSize>twoBlocks) ? P.get() : 0);
          }
        } // if (isVerbosity( OrthoDetails ))
        
        if (localSize < threeBlocks)
          localSize += _blockSize;
        continue;
      } // if (_nFound == 0)      


      // Order the Ritz eigenvectors by putting the converged vectors at the beginning
      int firstIndex = _blockSize;
      for (j = 0; j < _blockSize; j++) {
        if (_normR[j] >= _residual_tolerance) {
          firstIndex = j;
          break;
        }
      } // for (j = 0; j < _blockSize; j++)
      //
      // If the firstIndex comes before the number found, then we need
      // to move the converged eigenvectors to the front of the spectral
      // transformation.
      //

      while (firstIndex < _nFound) {
        for (j = firstIndex; j < _blockSize; j++) {
          if (_normR[j] < _residual_tolerance) {
            //
            // Swap and j-th and firstIndex-th position
            //
            std::swap_ranges<ScalarType*,ScalarType*>(S[j],S[j]+localSize,S[firstIndex]);

            // Swap _theta
            std::swap<MagnitudeType>(_theta[j],_theta[firstIndex]);
            
            // Swap _normR
            std::swap<MagnitudeType>(_normR[j],_normR[firstIndex]);
            break;
          }
        }
        for (j = 0; j < _blockSize; j++) {
          if (_normR[j] >= _residual_tolerance) {
            firstIndex = j;
            break;
          }
        } // for (j = 0; j < _blockSize; j++)
      } // while (firstIndex < _nFound)      


      //
      // Copy the converged eigenvalues and residuals.
      //
      leftOver = 0;
      if (_nFound > 0) {
        if ( _knownEV + _nFound  > _nev ) {
          leftOver = _knownEV + _nFound - _nev;
        }
        // Copy first _nFound-leftOver elements in _theta to range
        // [_knownEV,_knownEV+_nFound-leftOver) in *_evals
        std::copy<MTiter,STiter>(_theta.begin(),_theta.begin()+_nFound-leftOver,
                                 _evals->begin()+_knownEV);
        // Copy first _nFound-leftOver elements in _normR to range
        // [_knownEV,_knownEV+_nFound-leftOver) in _resids
        std::copy<MTiter,MTiter>(_normR.begin(),_normR.begin()+_nFound-leftOver,
                                 _resids.begin()+_knownEV);
      }


      //
      // Compute and store the converged eigenvectors 
      //
      if (_nFound > 0) {

        // Get a view of the eigenvector storage where the new eigenvectors will be put.
        std::vector<int> index2( _nFound - leftOver );
        for(j=0; j<_nFound-leftOver; j++)
          index2[j] = _knownEV + j;
        Teuchos::RefCountPtr<MV> EVconv = MVT::CloneView( *_evecs, index2 );
        
        Teuchos::SerialDenseMatrix<int,ScalarType> S11conv( Teuchos::View, S, _blockSize, _nFound-leftOver );
        MVT::MvTimesMatAddMv( one, *X, S11conv, zero, *EVconv );
        
        if (localSize >= twoBlocks) {
          Teuchos::SerialDenseMatrix<int,ScalarType> S21( Teuchos::View, S, _blockSize, _nFound-leftOver, _blockSize );
          MVT::MvTimesMatAddMv( one, *H, S21, one, *EVconv );
          
          if (localSize == threeBlocks) {
            Teuchos::SerialDenseMatrix<int,ScalarType> S31( Teuchos::View, S, _blockSize, _nFound-leftOver, twoBlocks  );
            MVT::MvTimesMatAddMv( one, *P, S31, one, *EVconv );
          }
        }	

        //
        // Sort the computed eigenvalues, eigenvectors, and residuals
        //
        _order.resize(_knownEV+_nFound-leftOver);
	{
	  Teuchos::TimeMonitor SortTimer( *_timerSortEval );
	  ret = _sm->sort( this, _knownEV+_nFound-leftOver, &(*_evals)[0], &_order);
	}
        if (ret != Ok) {
          if (_om->isVerbosityAndPrint(Error)) {
            _os << "ERROR : Sorting in solve()" << endl;
          }
          _error_flg = true;
          break;
        }
        // use _order to permute _evecs and _resids
        _MSUtils.permuteVectors(_knownEV+_nFound-leftOver,_order,*_evecs,&_resids);
        
        // Increment number of known eigenpairs.
        _knownEV += (_nFound-leftOver);
        
        // We don't need to define restarting vectors, so break out of this loop.
        if ( _knownEV >= _nev )
          break;
      }

      // Define the restarting vectors
      {
	Teuchos::TimeMonitor RestartTimer( *_timerRestart );
	leftOver = (_nevLocal < _blockSize + _nFound) ? _nevLocal - _nFound : _blockSize;
	
	// Grab a view of the restarting vectors (a subview of X)
	std::vector<int> index2( leftOver );
	for (j=0; j<leftOver; j++) {
	  index2[j] = j;
	}
	Teuchos::RefCountPtr<MV>  Xtmp = MVT::CloneView( *X, index2 );
	Teuchos::RefCountPtr<MV> KXtmp = MVT::CloneView(*KX, index2 );
	Teuchos::RefCountPtr<MV> MXtmp = MVT::CloneView(*MX, index2 );
	
	Teuchos::SerialDenseMatrix<int,ScalarType> S11new( Teuchos::View, S, _blockSize, leftOver, 0, _nFound );
	
	MVT::MvAddMv( one, *X, zero, *X, *R );
	MVT::MvTimesMatAddMv( one, *R, S11new, zero, *Xtmp );
	
	MVT::MvAddMv( one, *KX, zero, *KX, *R );
	MVT::MvTimesMatAddMv( one, *R, S11new, zero, *KXtmp );
	
	if (haveMass) {
	  MVT::MvAddMv( one, *MX, zero, *MX, *R );
	  MVT::MvTimesMatAddMv( one, *R, S11new, zero, *MXtmp );
	}
	
	if (localSize >= twoBlocks) {
	  Teuchos::SerialDenseMatrix<int,ScalarType> S21new( Teuchos::View, S, _blockSize, leftOver, _blockSize, _nFound );        
	  MVT::MvTimesMatAddMv( one, *H, S21new, one, *Xtmp );
	  MVT::MvTimesMatAddMv( one, *KH, S21new, one, *KXtmp );
	  if (haveMass) {
	    MVT::MvTimesMatAddMv( one, *MH, S21new, one, *MXtmp );
	  }
	  
	  if (localSize >= threeBlocks) {
	    Teuchos::SerialDenseMatrix<int,ScalarType> S31new( Teuchos::View, S, _blockSize, leftOver, twoBlocks, _nFound );
	    MVT::MvTimesMatAddMv( one, *P, S31new, one, *Xtmp );
	    MVT::MvTimesMatAddMv( one, *KP, S31new, one, *KXtmp );
	    if (haveMass) {
	      MVT::MvTimesMatAddMv( one, *MP, S31new, one, *MXtmp );
	    }
	  }
	}
	if (_nevLocal < _blockSize + _nFound) {
	  // Put new random vectors at the end of the block
	  index2.resize( _blockSize - leftOver );
	  for (j=0; j<_blockSize-leftOver; j++)
	    index2[j] = leftOver + j;
	  Xtmp = MVT::CloneView( *X, index2 );
	  MVT::MvRandom( *Xtmp );
	}
	else {
	  _nFound = 0;
	}
	reStart = true;

      } // end timing block
	
    } // for (_iter = 0; _iter < _maxIter; ++_iter)

    //
    // Print out a final summary if necessary
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

    if (info != 0) {
      return Failed;
    }
    else if (_knownEV != _nev) {
      return Unconverged;
    }
    return Ok;
  }
  
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::accuracyCheck(const MV *X, const MV *MX,
                                                      const MV *R, const MV *Q,
                                                      const MV *H, const MV *P) const 
  {    
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    bool haveMass = (_MOp.get()!=0);
    MagnitudeType tmp;
    
    _os << " Checking Orthogonality after Iteration : "<< _iter << endl;
    if (X) {
      if (haveMass) {
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
    
    if ((R) && (X)) {
      tmp = _MSUtils.errorOrthogonality(X, R);
      if (_om->doPrint())
        _os << " >> Orthogonality X^T R up to " << tmp << endl;
    }
    
    if (Q == 0) {
      _os << endl;
      return;
    }    

    if (haveMass) {
      tmp = _MSUtils.errorOrthonormality(Q, _MOp.get());
      if (_om->doPrint())
        _os << " >> Error in Q^T M Q - I = " << tmp << endl;
      if (X) {
        tmp = _MSUtils.errorOrthogonality(Q, X, _MOp.get());
        if (_om->doPrint())
          _os << " >> Orthogonality Q^T M X up to " << tmp << endl;
      }
      if (H) {
        tmp = _MSUtils.errorOrthogonality(Q, H, _MOp.get());
        if (_om->doPrint())
          _os << " >> Orthogonality Q^T M H up to " << tmp << endl;
      }
      if (P) {
        tmp = _MSUtils.errorOrthogonality(Q, P, _MOp.get());
        if (_om->doPrint())
          _os << " >> Orthogonality Q^T M P up to " << tmp << endl;
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
      if (H) {
        tmp = _MSUtils.errorOrthogonality(Q, H, 0);
        if (_om->doPrint())
          _os << " >> Orthogonality Q^T H up to " << tmp << endl;
      }
      if (P) {
        tmp = _MSUtils.errorOrthogonality(Q, P, 0);
        if (_om->doPrint())
          _os << " >> Orthogonality Q^T P up to " << tmp << endl;
      }
    }
    _os << endl;
  }      
  
} // end Anasazi namespace

#endif // ANASAZI_LOBPCG_HPP
