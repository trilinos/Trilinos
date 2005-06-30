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
#include "AnasaziModalSolverUtils.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"

/*!	\class Anasazi::LOBPCG

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
    // Internal methods
    //
  void accuracyCheck(const MV *X, const MV *MX,
		     const MV *R, const MV *Q,
		     const MV *H, const MV *P) const;
    //
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem; 
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
    const ScalarType _residual_tolerance;
    int _numRestarts, _iter, _knownEV, _nevLocal;
    std::vector<ScalarType> _theta, _normR, _resids;
    //
    // Internal utilities class required by eigensolver.
    //
    ModalSolverUtils<ScalarType,MV,OP> _MSUtils;
    //
    // Output stream from the output manager
    std::ostream& _os;

    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
  };
  //
  // Implementation
  //
  template <class ScalarType, class MV, class OP>
  LOBPCG<ScalarType,MV,OP>::LOBPCG(const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
				   const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
				   Teuchos::ParameterList &pl
				   ):
    _problem(problem), 
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
    _residual_tolerance(_pl.get("Tol", 1.0e-6)),
    _numRestarts(0), 
    _iter(0), 
    _knownEV(0),
    _nevLocal(0),
    _MSUtils(om),
    _os(_om->GetOStream())
  {     
  }

  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::currentStatus() 
  {
    int i;
    if (_om->doPrint()) {
      cout.setf(ios::scientific, ios::floatfield);  
      cout.precision(6);
      _os <<endl;
      _os <<"********************CURRENT STATUS********************"<<endl;
      _os <<"Iterations :\t"<<_iter<<endl;
      _os <<"Restarts :\t"<<_numRestarts<<endl;
      _os <<"Block Size :\t"<<_blockSize<<endl;
      _os <<"Requested Eigenvalues : "<<_nev<<endl;
      _os <<"Computed Eigenvalues : "<<_knownEV<<endl;
      _os <<"Residual Tolerance : "<<_residual_tolerance<<endl;
      _os <<"------------------------------------------------------"<<endl;
      _os <<"Computed Eigenvalues: "<<endl;
      _os <<"------------------------------------------------------"<<endl;
      _os <<"Eigenvalue\tResidual"<<endl;
      _os <<"------------------------------------------------------"<<endl;
      if ( _knownEV > 0 ) {
	for (i=0; i<_knownEV; i++)
	  _os <<(*_evals)[i]<<"\t"<<_resids[i]<<endl;
      } else {
	_os <<"[none computed]"<<endl;
      }
      _os <<"------------------------------------------------------"<<endl;
      _os <<"Current Eigenvalue Estimates (Ritz Values): "<<endl;
      _os <<"------------------------------------------------------"<<endl;
      _os <<"Ritz Value\tResidual"<<endl;
      _os <<"------------------------------------------------------"<<endl;
      if ( _iter > 0 || _nevLocal > 0 ) {
	for (i=0; i<_blockSize; i++)
	  _os <<_theta[i]<<"\t"<<_normR[i]<<endl;
      } else {
	_os <<"[none computed]"<<endl;
      }
      _os << endl;
    }
  }
  
  template <class ScalarType, class MV, class OP>
  ReturnType LOBPCG<ScalarType,MV,OP>::solve () 
  {
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
    // Reinitialize internal data and pointers, preparse for solve
    //
    _numRestarts = 0; 
    _iter = 0; 
    _knownEV = 0;
    //
    // Necessary variables
    //
    int i, j;
    int info, nb, leftOver;
    int bStart = 0, offSet = 0;
    bool reStart = false;
    bool criticalExit = false;
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::LAPACK<int,ScalarType> lapack;
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
    MVT::MvRandom( *X );
    //
    // Check to see if there is a mass matrix, so we know how much space is required.
    // [ If there isn't a mass matrix we can use a view of X.]
    //
    if (_MOp.get())
      MX = MVT::Clone( *iVec, _blockSize );
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
    if (_MOp.get())
      MH = MVT::Clone( *iVec, _blockSize );
    else {
      std::vector<int> index( _blockSize );
      for (int i=0; i<_blockSize; i++)
	index[i] = i;
      MH = MVT::CloneView( *X, index );
    }
    //
    // Search direction working vectors.
    Teuchos::RefCountPtr<MV> P, MP, KP;
    //
    P = MVT::Clone( *iVec, _blockSize );
    KP = MVT::Clone( *iVec, _blockSize );
    //
    if (_MOp.get())
      MP = MVT::Clone( *iVec, _blockSize );
    else {
      std::vector<int> index( _blockSize );
      for (int i=0; i<_blockSize; i++)
	index[i] = i;
      MP = MVT::CloneView( *X, index );
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
    // Initialize the workspace.
    //
    MVT::MvRandom( *X );
    //
    // Miscellaneous definitions
    int localSize;
    int twoBlocks = 2*_blockSize;
    int threeBlocks = 3*_blockSize;
    int _nFound = _blockSize;
    //    
    for ( _iter=0; _iter <= _maxIter; _iter++ ) {
      
      if ( (_iter==0) || (reStart == true) ) {
	
	reStart = false;
	localSize = _blockSize;
	
	if (_nFound > 0) {
	  
	  std::vector<int> index( _nFound );
	  for (i=0; i<_nFound; i++)
	    index[i] = _blockSize-_nFound + i;
	  Teuchos::RefCountPtr<MV> X2 = MVT::CloneView( *X, index );
	  Teuchos::RefCountPtr<MV> MX2 = MVT::CloneView( *MX, index );
	  Teuchos::RefCountPtr<MV> KX2 = MVT::CloneView( *KX, index );
	  
	  // Apply the mass matrix to X
	  //timeMassOp -= MyWatch.WallTime();
	  if (_MOp.get())
	    OPT::Apply( *_MOp, *X2, *MX2 );
	  //timeMassOp += MyWatch.WallTime();
	  //massOp += _nFound;
	  
	  if (_knownEV > 0) {
	    
	    // Orthonormalize X against the known eigenvectors with Gram-Schmidt
	    // Note: Use R as a temporary work space
	    index.resize( _knownEV );
	    for (i=0; i<_knownEV; i++)
	      index[i] = i;
	    Teuchos::RefCountPtr<MV> copyQ = MVT::CloneView( *_evecs, index );
	    
	    //timeOrtho -= MyWatch.WallTime();
	    info = _MSUtils.massOrthonormalize( *X, *MX, _MOp.get(), *copyQ, _nFound, 0 );
	    //timeOrtho += MyWatch.WallTime();
	    
	    // Exit the code if the orthogonalization did not succeed
	    if (info < 0) {
	      info = -10;
	      //return info;
	    }
	  } // if (knownEV > 0) 
	  
	  // Apply the stiffness matrix to X
	  //        timeStifOp -= MyWatch.WallTime();
	  OPT::Apply( *_Op, *X2, *KX2 );
	  //timeStifOp += MyWatch.WallTime();
	  //stifOp += _nFound;
	  
	} // if (_nFound > 0)
	
      } // if ( (_iter==0) || reStart == true)

      else {
	
	// Apply the preconditioner on the residuals
	if (_Prec.get()) {
	  //timePrecOp -= MyWatch.WallTime();
	  OPT::Apply( *_Prec, *R, *H );
	  //timePrecOp += MyWatch.WallTime();
	  //precOp += blockSize;
	}
	else {
	  MVT::MvAddMv( one, *R, zero, *R, *H );
	}

	// Apply the mass matrix on H
	// timeMassOp -= MyWatch.WallTime();
	if (_MOp.get())
	  OPT::Apply( *_MOp, *H, *MH);
	//timeMassOp += MyWatch.WallTime();
	//massOp += blockSize;
	
	if (_knownEV > 0) {

	  // Orthogonalize H against the known eigenvectors
	  // Note: Use R as a temporary work space
	  std::vector<int> index( _knownEV );
	  for (i=0; i<_knownEV; i++)
	    index[i] = i;
	  Teuchos::RefCountPtr<MV> copyQ = MVT::CloneView( *_evecs, index );

	  //timeOrtho -= MyWatch.WallTime();
	  _MSUtils.massOrthonormalize( *H, *MH, _MOp.get(), *copyQ, _blockSize, 1);
	  //timeOrtho += MyWatch.WallTime();

	} // if (knownEV > 0)
	
	// Apply the stiffness matrix to H
	// timeStifOp -= MyWatch.WallTime();
	OPT::Apply( *_Op, *H, *KH);
	//timeStifOp += MyWatch.WallTime();
	//stifOp += blockSize;
	
	if (localSize == _blockSize)
	  localSize += _blockSize;

      } // if ( (_iter==0) || (reStart==true))

      // Form "local" mass and stiffness matrices
      //timeLocalProj -= MyWatch.WallTime();
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
      //timeLocalProj += MyWatch.WallTime();
          
      // Perform a spectral decomposition
      //timeLocalSolve -= MyWatch.WallTime();
      _nevLocal = localSize;
      info = _MSUtils.directSolver(localSize, KK, &MM, &S, &_theta, &_nevLocal, 
				   (_blockSize == 1) ? 1 : 0);
      //timeLocalSolve += MyWatch.WallTime();
      
      if (info < 0) {
	// Stop when spectral decomposition has a critical failure
        break;
      } // if (info < 0)
      
      // Check for restarting
      if ((_theta[0] < 0.0) || (_nevLocal < _blockSize)) {
	if (_om->isVerbosityAndPrint( IterationDetails ) ) {
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
      } // if ((theta[0] < 0.0) || (_nevLocal < _blockSize))
    

      if ((localSize == twoBlocks) && (_nevLocal == _blockSize)) {
	//for (j = 0; j < _nevLocal; ++j) 
	//memcpy(S + j*_blockSize, S + j*twoBlocks, _blockSize*sizeof(double)); 
	_os << "localSize == twoBlocks && _nevLocal == _blockSize"<<endl;
	localSize = _blockSize;
      }
      
      if ((localSize == threeBlocks) && (_nevLocal <= twoBlocks)) {
	//for (j = 0; j < _nevLocal; ++j) 
	// memcpy(S + j*twoBlocks, S + j*threeBlocks, twoBlocks*sizeof(double)); 
	_os << "localSize == threeBlocks && _nevLocal <= twoBlocks"<<endl;
	localSize = twoBlocks;
      }
      
      // Compute the residuals
      //timeResidual -= MyWatch.WallTime();
      Teuchos::SerialDenseMatrix<int,ScalarType> S11( Teuchos::View, S, _blockSize, _blockSize );
      MVT::MvTimesMatAddMv( 1.0, *KX, S11, 0.0, *R );
      
      if (localSize >= twoBlocks) {
	Teuchos::SerialDenseMatrix<int,ScalarType> S21( Teuchos::View, S, _blockSize, _blockSize, _blockSize );
	MVT::MvTimesMatAddMv( 1.0, *KH, S21, 1.0, *R );
	
	if (localSize == threeBlocks) {
	  Teuchos::SerialDenseMatrix<int,ScalarType> S31( Teuchos::View, S, _blockSize, _blockSize, twoBlocks  );	  
	  MVT::MvTimesMatAddMv( 1.0, *KP, S31, 1.0, *R );

	} // if (localSize == threeBlocks)
      } // if (localSize >= twoBlocks )
      
      for (j = 0; j < _blockSize; ++j)
	blas.SCAL(localSize, _theta[j], S[j], 1);
      
      MVT::MvTimesMatAddMv( -one, *MX, S11, one, *R );
      
      if (localSize >= twoBlocks) {
	Teuchos::SerialDenseMatrix<int,ScalarType> S21( Teuchos::View, S, _blockSize, _blockSize, _blockSize );
	MVT::MvTimesMatAddMv( -one, *MH, S21, one, *R );
	
	if (localSize == threeBlocks) {
	  Teuchos::SerialDenseMatrix<int,ScalarType> S31( Teuchos::View, S, _blockSize, _blockSize, twoBlocks  );	  
	  MVT::MvTimesMatAddMv( -one, *MP, S31, one, *R );
	}
      } // if (localSize >= twoBlocks)
            
      for (j = 0; j < _blockSize; ++j)
	blas.SCAL(localSize, one/_theta[j], S[j], 1);

      //timeResidual += MyWatch.WallTime();
      
      // Compute the norms of the residuals
      //timeNorm -= MyWatch.WallTime();
      _problem->MvNorm( *R, &_normR );
      
      // Scale the norms of residuals with the eigenvalues
      // Count the converged eigenvectors
      _nFound = 0;
      for (j = 0; j < _blockSize; ++j) {
	_normR[j] = (_theta[j] == 0.0) ? _normR[j] : _normR[j]/_theta[j];
	if (_normR[j] < _residual_tolerance) 
	  _nFound ++;
      }
   
      //timeNorm += MyWatch.WallTime();
      
      // Store the residual history
      /*      if (localVerbose > 2) {
	      memcpy(resHistory + historyCount*blockSize, normR, blockSize*sizeof(double));
	      historyCount += 1;
	      }
      */      
      
      // Print information on current iteration
      if (_om->isVerbosityAndPrint( IterationDetails )) 
	currentStatus();

      
      if (_nFound == 0) {
	// Update the spaces
	// Note: Use R as a temporary work space
	// Note: S11 was previously defined above
	//timeLocalUpdate -= MyWatch.WallTime();
	
	MVT::MvAddMv( one, *X, zero, *X, *R );	
	MVT::MvTimesMatAddMv( one, *R, S11, zero, *X );
	
	MVT::MvAddMv( one, *KX, zero, *KX, *R );
	MVT::MvTimesMatAddMv( one, *R, S11, zero, *KX );
	
	if (_MOp.get()) {
	  MVT::MvAddMv( one, *MX, zero, *MX, *R );
	  MVT::MvTimesMatAddMv( one, *R, S11, zero, *MX );
	}
	
	if (localSize == twoBlocks) {
	  Teuchos::SerialDenseMatrix<int,ScalarType> S21( Teuchos::View, S, _blockSize, _blockSize, _blockSize );
	  MVT::MvTimesMatAddMv( one, *H, S21, zero, *P );
	  MVT::MvAddMv( one, *P, one, *X, *X );
	  
	  MVT::MvTimesMatAddMv( one, *KH, S21, zero, *KP );
	  MVT::MvAddMv( one, *KP, one, *KX, *KX );
	  
	  if (_MOp.get()) {
	    MVT::MvTimesMatAddMv( one, *MH, S21, zero, *MP );
	    MVT::MvAddMv( one, *MP, one, *MX, *MX );
	  }
	} // if (localSize == twoBlocks)
	
	if (localSize == threeBlocks) {	  
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
	  
	  if (_MOp.get()) {
	    MVT::MvAddMv( one, *MP, zero, *MP, *R );
	    MVT::MvTimesMatAddMv( one, *R, S31, zero, *MP );
	    MVT::MvTimesMatAddMv( one, *MH, S21, one, *MP );
	    MVT::MvAddMv( one, *MP, one, *MX, *MX );
	  }
	} // if (localSize == threeBlocks)
	
	//timeLocalUpdate += MyWatch.WallTime();
	
	// Compute the new residuals
	//timeResidual -= MyWatch.WallTime();
	MVT::MvAddMv( one, *KX, zero, *KX, *R );
	Teuchos::SerialDenseMatrix<int,ScalarType> T( _blockSize, _blockSize );
	for (j = 0; j < _blockSize; ++j) 
	  T(j,j) = _theta[j];
	MVT::MvTimesMatAddMv( -one, *MX, T, one, *R );
	//timeResidual += MyWatch.WallTime();
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
      for (j = 0; j < _blockSize; ++j) {
	if (_normR[j] >= _residual_tolerance) {
	  firstIndex = j;
	  break;
	}
      } // for (j = 0; j < _blockSize; ++j)
      //
      // If the firstIndex comes before the number found, then we need
      // to move the converged eigenvectors to the front of the spectral
      // transformation.
      //
      ScalarType tmp_swap;
      std::vector<ScalarType> tmp_swap_vec(localSize);
      while (firstIndex < _nFound) {
	for (j = firstIndex; j < _blockSize; ++j) {
	  if (_normR[j] < _residual_tolerance) {
	    //
	    // Swap and j-th and firstIndex-th position
	    //
	    blas.COPY(localSize, S[ j ], 1, &tmp_swap_vec[0], 1);
	    blas.COPY(localSize, S[ firstIndex ], 1, S[ j ], 1 );
	    blas.COPY(localSize, &tmp_swap_vec[0], 1, S[ firstIndex ], 1 );

	    // Swap _theta
	    tmp_swap = _theta[j];
	    _theta[j] = _theta[firstIndex];
	    _theta[firstIndex] = tmp_swap;

	    // Swap _normR
	    tmp_swap = _normR[j];
	    _normR[j] = _normR[firstIndex];
	    _normR[firstIndex] = tmp_swap;
	    break;
	  }
	}
	for (j = 0; j < _blockSize; ++j) {
	  if (_normR[j] >= _residual_tolerance) {
	    firstIndex = j;
	    break;
	  }
	} // for (j = 0; j < _blockSize; ++j)
      } // while (firstIndex < _nFound)      
      //
      // Copy the converged eigenvalues and residuals.
      //
      leftOver = 0;
      if (_nFound > 0) {
	if ( _knownEV + _nFound  > _nev ) 
	  leftOver = _knownEV + _nFound -_nev;
	blas.COPY( _nFound-leftOver, &_theta[0], 1, &(*_evals)[_knownEV], 1 ); 
	blas.COPY( _nFound-leftOver, &_normR[0], 1, &_resids[_knownEV], 1 );
      }
      //
      // Compute the current eigenvector estimates
      //
      MVT::MvTimesMatAddMv( one, *X, S11, zero, *R );
      
      if (localSize >= twoBlocks) {
	Teuchos::SerialDenseMatrix<int,ScalarType> S21( Teuchos::View, S, _blockSize, _blockSize, _blockSize );	  
	MVT::MvTimesMatAddMv( one, *H, S21, one, *R );
	
	if (localSize == threeBlocks) {
	  Teuchos::SerialDenseMatrix<int,ScalarType> S31( Teuchos::View, S, _blockSize, _blockSize, twoBlocks  );
	  MVT::MvTimesMatAddMv( one, *P, S31, one, *R );
	}
      }
      
      // Store the converged eigenvectors 
      if (_nFound) {

	leftOver = 0;
	if ( _knownEV + _nFound  > _nev)
	  leftOver = _knownEV + _nFound - _nev;
	
	// Get a view of the converged eigenvectors
	std::vector<int> index2( _nFound - leftOver );
	for(j=0; j<_nFound-leftOver; j++)
	  index2[j] = j;
	Teuchos::RefCountPtr<MV> Rconv = MVT::CloneView( *R, index2 );
	
	// Get a view of the eigenvector storage where the new eigenvectors will be put.
	for(j=0; j<_nFound-leftOver; j++)
	  index2[j] = _knownEV + j;
	Teuchos::RefCountPtr<MV> EVconv = MVT::CloneView( *_evecs, index2 );
	
	// Put the converged eigenvectors in the storage
	MVT::MvAddMv( one, *Rconv, zero, *Rconv, *EVconv );

	// Sort the eigenpairs
	//timePostProce -= MyWatch.WallTime();
	if ((info==0) && (_knownEV > 0))
	  _MSUtils.sortScalars_Vectors(_knownEV, &(*_evals)[0], _evecs.get(), &_resids);     
	//timePostProce += MyWatch.WallTime();

	// Increment number of known eigenpairs.
	_knownEV += (_nFound-leftOver);

      }

      // We don't need to define restarting vectors, so break out of this loop.
      if ( _knownEV >= _nev )
	break;

      // Define the restarting vectors
      //timeRestart -= MyWatch.WallTime();
      leftOver = (_nevLocal < _blockSize + _nFound) ? _nevLocal - _nFound : _blockSize;
      Teuchos::SerialDenseMatrix<int,ScalarType> S11new( Teuchos::View, S, _blockSize, leftOver, 0, _nFound );
      
      MVT::MvAddMv( one, *X, zero, *X, *R );
      MVT::MvTimesMatAddMv( one, *R, S11new, zero, *X );
      
      MVT::MvAddMv( one, *KX, zero, *KX, *R );
      MVT::MvTimesMatAddMv( one, *R, S11new, zero, *KX );
      
      if (_MOp.get()) {
	MVT::MvAddMv( one, *MX, zero, *MX, *R );
	MVT::MvTimesMatAddMv( one, *R, S11new, zero, *MX );
      }

      if (localSize >= twoBlocks) {
	Teuchos::SerialDenseMatrix<int,ScalarType> S21new( Teuchos::View, S, _blockSize, leftOver, _blockSize, _nFound );	
	MVT::MvTimesMatAddMv( one, *H, S21new, one, *X );
	MVT::MvTimesMatAddMv( one, *KH, S21new, one, *KX );
	if (_MOp.get()) 
	  MVT::MvTimesMatAddMv( one, *MH, S21new, one, *MX );
	
	if (localSize >= threeBlocks) {
	  Teuchos::SerialDenseMatrix<int,ScalarType> S31new( Teuchos::View, S, _blockSize, leftOver, twoBlocks, _nFound );
	  MVT::MvTimesMatAddMv( one, *P, S31new, one, *X );
	  MVT::MvTimesMatAddMv( one, *KP, S31new, one, *KX );
	  if (_MOp.get())
	    MVT::MvTimesMatAddMv( one, *MP, S31new, one, *MX );
        }
      }
      if (_nevLocal < _blockSize + _nFound) {
	// Put new random vectors at the end of the block
	std::vector<int> index2( _blockSize - leftOver );
	for (j=0; j<_blockSize-leftOver; j++)
	  index2[j] = leftOver + j;
	Teuchos::RefCountPtr<MV> Xtmp = MVT::CloneView( *X, index2 );
	MVT::MvRandom( *Xtmp );
      }
      else {
	_nFound = 0;
      }
      reStart = true;
      //timeRestart += MyWatch.WallTime();
      
    } // for (_iter = 0; _iter < _maxIter; ++_iter)

    //timeOuterLoop += MyWatch.WallTime();
    //highMem = (highMem > currentSize()) ? highMem : currentSize();
    //
    // Print out a final summary if necessary
    //
    if (_om->isVerbosity( FinalSummary ))
      currentStatus();    

    if (info != 0)
      return Failed;

    return Ok;
  }
  
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::accuracyCheck(const MV *X, const MV *MX,
						      const MV *R, const MV *Q,
						      const MV *H, const MV *P) const 
  {    
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    ScalarType tmp;
    
    _os << " Checking Orthogonality after Iteration : "<< _iter << endl;
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
    
    if ((R) && (X)) {
      tmp = _MSUtils.errorOrthogonality(X, R);
      if (_om->doPrint())
	_os << " >> Orthogonality X^T R up to " << tmp << endl;
    }
    
    if (Q == 0)
      return;
    
    if (_MOp.get()) {
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
