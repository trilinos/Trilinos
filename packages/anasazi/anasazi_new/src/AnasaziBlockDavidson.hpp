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

#ifndef ANASAZI_BLOCK_DAVIDSON_HPP
#define ANASAZI_BLOCK_DAVIDSON_HPP

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

/*!	\class Anasazi::BlockDavidson

	\brief This class implements the Restarted Block Krylov Schur Method,
	an iterative method for solving eigenvalue problems.

	This method is a block version of the method presented by G.W. Stewart 
	in "A Krylov-Schur Algorithm for Large Eigenproblems", 
	SIAM J. Matrix Anal. Appl., Vol 28, No. 8, pp. 601-614.

	\author Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {
  
  template <class ScalarType, class MV, class OP>
  class BlockDavidson : public Eigensolver<ScalarType,MV,OP> { 
  public:
    //@{ \name Constructor/Destructor.
    
    //! %Anasazi::BlockDavidson constructor.
    BlockDavidson( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
		   const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
		   const ScalarType tol=1.0e-6,
		   const int blockSize = 1,
		   const int length=25, 
		   const int maxIter=300 
		   );
    
    //! %Anasazi::BlockDavidson destructor.
    virtual ~BlockDavidson() {};
    //@}
    
    //@{ \name Solver application methods.
    
    /*! \brief This method uses iterate to compute approximate solutions to the
      original problem.  It may return without converging if it has taken the
      maximum number of iterations or numerical breakdown is observed.
    */
    void solve();
    //@}
    
    //@{ \name Solver status methods.
    
    //! Get the current iteration count.
    int GetNumIters() const { return(_iter); };
    
    //! Get the current restart count of the iteration method.
    /*! Some eigensolvers can perform restarts (i.e.Arnoldi) to reduce memory
      and orthogonalization costs.  For other eigensolvers, like LOBPCG or block Davidson,
      this is not a valid stopping criteria.
    */
    int GetNumRestarts() const { return(_restartIter); };

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
    void accuracyCheck(const MV *X, const MV *MX, const MV *Q) const; 
    //
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem; 
    const Teuchos::RefCountPtr<OutputManager<ScalarType> > _om; 
    //
    // Information obtained from the eigenproblem
    //
    Teuchos::RefCountPtr<OP> _Op;
    Teuchos::RefCountPtr<OP> _BOp;
    Teuchos::RefCountPtr<OP> _Prec;
    Teuchos::RefCountPtr<MV> _evecs;
    Teuchos::RefCountPtr<std::vector<ScalarType> > _evals;
    const int _nev;  
    //
    // Internal data.
    //
    const int _numBlocks, _maxIter, _blockSize;
    const ScalarType _residual_tolerance;
    int _restartIter, _iter, _dimSearch, _knownEV;
    int _lwork;
    //
    // Internal utilities class required by eigensolver.
    //
    ModalSolverUtils<ScalarType,MV,OP> _MSUtils;
    //
    // Internal storage for eigensolver
    //
    Teuchos::RefCountPtr<MV> _Xvec;
    Teuchos::RefCountPtr<MV> _MXvec;
    Teuchos::RefCountPtr<MV> _KXvec; 
    Teuchos::RefCountPtr<MV> _Rvec;
    Teuchos::SerialDenseMatrix<int,ScalarType> _KKsdm, _Ssdm;
    std::vector<ScalarType> _theta, _normR, _work, _resids;

    // Output stream from the output manager
    std::ostream& _os;

    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
  };
  //
  // Implementation
  //
  template <class ScalarType, class MV, class OP>
  BlockDavidson<ScalarType,MV,OP>::BlockDavidson(const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
					    const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
					    const ScalarType tol,
					    const int blockSize,
					    const int numBlocks, 
					    const int maxIter
					    ): 
    _problem(problem), 
    _om(om),
    _Op(_problem->GetOperator()),
    _BOp(_problem->GetB()),
    _Prec(_problem->GetPrec()),
    _evecs(_problem->GetEvecs()), 
    _evals(problem->GetEvals()), 
    _nev(problem->GetNEV()), 
    _numBlocks(numBlocks), 
    _maxIter(maxIter),
    _blockSize(blockSize),
    _residual_tolerance(tol),
    _restartIter(0), 
    _iter(0), 
    _dimSearch(0),    
    _knownEV(0),
    _lwork(0),
    _MSUtils(om),
    _os(_om->GetOStream())
  {     
    //
    // Retrieve the initial vector and operator information from the Anasazi::Eigenproblem.
    //
    Teuchos::RefCountPtr<MV> iVec = _problem->GetInitVec();
    //
    // Define local block vectors
    //
    // MX = Working vectors (storing M*X if M is specified, else pointing to X)
    // KX = Working vectors (storing K*X)
    //
    // R = Residuals
    //
    _dimSearch = _blockSize*_numBlocks;
    _Xvec = MVT::Clone( *iVec, _dimSearch + _blockSize );
    _KXvec = MVT::Clone( *iVec, _blockSize );
    _Rvec = MVT::Clone( *iVec, _blockSize );
    //
    // Check to see if there is a mass matrix, so we know how much space is required.
    // [ If there isn't a mass matrix we can use a view of _Xvec.]
    //
    if (_BOp.get())
      _MXvec = MVT::Clone( *iVec, _blockSize );
    else {
      std::vector<int> index( _blockSize );
      for (int i=0; i<_blockSize; i++)
	index[i] = i;
      _MXvec = MVT::CloneView( *_Xvec, index );
    }
    //
    // Initialize the workspace.
    //
    MVT::MvRandom( *_Xvec );
    //
    // Define dense local matrices and arrays
    //
    // theta = Storage for local eigenvalues     (size: dimSearch)
    // normR = Storage for the norm of residuals (size: blockSize)
    //
    // KK = Local stiffness matrix               (size: dimSearch x dimSearch)
    //
    // S = Local eigenvectors                    (size: dimSearch x dimSearch)
    //
    // tmpKK = Local workspace                   (size: dimSearch x dimSearch)
    //
    _theta.resize( _dimSearch );
    _normR.resize( _blockSize );
    _resids.resize( _nev );
    //    
    _KKsdm.shape( _dimSearch, _dimSearch );
    _Ssdm.shape( _dimSearch, _dimSearch );
        
  }

  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::currentStatus() 
  {
    int i;
    if (_om->doPrint()) {
      cout.setf(ios::scientific, ios::floatfield);  
      cout.precision(6);
      _os <<endl;
      _os <<"********************CURRENT STATUS********************"<<endl;
      _os <<"Iterations :\t"<<_iter<<endl;
      _os <<"Restarts :\t"<<_restartIter<<endl;
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
    }
  }
  
  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::solve () 
  {
    int i, j;
    int info, nb;
    int bStart = 0, offSet = 0;
    int nFound = _blockSize;
    bool reStart = false;
    bool criticalExit = false;
    Teuchos::RefCountPtr<MV> Xcurrent, Xprev, Xtotal, Xnext;
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::LAPACK<int,ScalarType> lapack;
    //
    // Determine the maximum number of blocks for this factorization.
    // ( NOTE:  This will be the _numBlocks since we don't know about already converged vectors here )
    int maxBlock = _numBlocks;
    //
    // Compute the required workspace for this factorization.
    //
    _lwork = lapack.ILAENV( 1, "geqrf", "", maxBlock*_blockSize, maxBlock*_blockSize );
    _lwork *= _blockSize;
    _work.resize( _lwork );
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
	Xcurrent = MVT::CloneView( *_Xvec, index );
	//
	if (_knownEV + localSize > 0) {
	  index.resize( _knownEV + localSize );
	  for (i=0; i < _knownEV + localSize; i++)
	    index[i] = i;
	  
	  Xprev = MVT::CloneView( *_Xvec, index );
	}
	//
	// Apply the mass matrix.
	//	
	if (_BOp.get())
	  OPT::Apply( *_BOp, *Xcurrent, *_MXvec );
	//
	// Orthonormalize Xcurrent again the known eigenvectors and previous vectors.
	//
	if (nb == bStart) {
	  if (nFound > 0) {
	    if (_knownEV == 0) {
	      info = _MSUtils.massOrthonormalize( *Xcurrent, *_MXvec, _BOp.get(), *Xcurrent, nFound, 2 );
	    }
	    else {
	      info = _MSUtils.massOrthonormalize( *Xcurrent, *_MXvec, _BOp.get(), *Xprev, nFound, 0 );
	    }
	  }
	  nFound = 0;
	} 
	else {
	  info = _MSUtils.massOrthonormalize( *Xcurrent, *_MXvec, _BOp.get(), *Xprev, _blockSize, 0 );
	}
	//
	// Exit the code if there has been a problem.
	//
	if (info < 0) 
	  return;
	//
	// Check orthogonality of X ( if required )
	//
	if (_om->isVerbosity( Debug ) ) {
	  if (localSize > 0)
	    accuracyCheck( Xcurrent.get(), _MXvec.get(), Xprev.get() );
	  else
	    accuracyCheck( Xcurrent.get(), _MXvec.get(), NULL );
	}
	//
	// Apply the stiffness matrix.
	//
	OPT::Apply( *_Op, *Xcurrent, *_KXvec );
	//
	// Update the local stiffness matrix ( Xtotal^T * K * Xcurrent where Xtotal = [Xprev Xcurrent] )
	// Note:  Only the upper half of the matrix is stored in KK
	//
	index.resize( localSize + _blockSize );
	for (i=0; i < localSize + _blockSize; i++)
	  index[i] = _knownEV + i;
	Xtotal = MVT::CloneView( *_Xvec, index );
	Teuchos::SerialDenseMatrix<int,ScalarType> subKK( Teuchos::View, _KKsdm, localSize+_blockSize, _blockSize, 0, localSize );
	MVT::MvTransMv( one, *Xtotal, *_KXvec, subKK );
	//
	// Perform spectral decomposition
	//
	info = _MSUtils.directSolver(localSize+_blockSize, _KKsdm, 0, &_Ssdm, &_theta, localSize+_blockSize, 10);
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
	  _restartIter++;
	  index.resize( _blockSize );
	  for (i=0; i<_blockSize; i++)
	    index[i] = _knownEV + i;
	  Teuchos::RefCountPtr<MV> Xinit = MVT::CloneView( *_Xvec, index );
	  MVT::MvRandom( *Xinit );
	  nFound = _blockSize;
	  bStart = 0;
	  break;
	}
	//
	// Update the search space :
	// KX = Xtotal * S where S is the eigenvectors of the projected problem.
	//
	Teuchos::SerialDenseMatrix<int,ScalarType> subS( Teuchos::View, _Ssdm, localSize+_blockSize, _blockSize );
	MVT::MvTimesMatAddMv( one, *Xtotal, subS, zero, *_KXvec );
	//
	// Apply the mass matrix for the next block
	// 
	if (_BOp.get())
	  OPT::Apply( *_BOp, *_KXvec, *_MXvec );
	//
	// Apply the stiffness matrix for the next block
	//
	OPT::Apply( *_Op, *_KXvec, *_Rvec );
	//
	// Compute the residual :
	// R = KX - diag(theta)*MX
	// 
	Teuchos::SerialDenseMatrix<int,ScalarType> D(_blockSize, _blockSize);
	for (i=0; i<_blockSize; i++ )
	  D(i,i) = -_theta[i];
	//
	if (_BOp.get()) {
	  MVT::MvTimesMatAddMv( one, *_MXvec, D, one, *_Rvec );
	}
	else {
	  MVT::MvTimesMatAddMv( one, *_KXvec, D, one, *_Rvec );
	}
	_problem->MvNorm( *_Rvec, &_normR );
	//
	// Scale the norms of residuals with the eigenvalues and check for converged eigenvectors.
	//
	nFound = 0;
	for (j=0; j<_blockSize; j++) {
	  // Scale eigenvalues if _theta is non-zero.
	  if ( _theta[j] != zero )
	    _normR[j] /= _theta[j];
	  // Check for convergence
	  if (_normR[j] < _residual_tolerance)
	    nFound ++;	  
	}
	// Print information on current iteration
	if (_om->isVerbosityAndPrint( Debug )) {
	  _os << " Iteration " << _iter << " - Number of converged eigenvectors ";
	  _os << _knownEV + nFound << endl;
	} 
	
	if (_om->isVerbosityAndPrint( Debug )) {
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
	    OPT::Apply( *_Prec, *_Rvec, *Xcurrent );
	  }
	  else
	    MVT::MvAddMv( one, *_Rvec, zero, *_Rvec, *Xcurrent );
	  //
	  // Update the preconditioned residual 
	  //
	  MVT::MvAddMv( one, *_KXvec, -one, *Xcurrent, *Xcurrent );
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
	for( i=0; i<_blockSize; i++) 
	  index[i] = _knownEV + localSize + _blockSize + i;
	Xnext = MVT::CloneView( *_Xvec, index );
	if (_Prec.get()) {
	  OPT::Apply( *_Prec, *_Rvec, *Xnext );
	}
	else 
	  MVT::MvAddMv( one, *_Rvec, zero, *_Rvec, *Xnext );
	
      } // for (nb = bStart; nv < maxBlock; nb++)

      //
      // Check if there is any reason we need to skip the rest of the while loop.
      //
      if (_iter > _maxIter)
	break;

      if (reStart == true) {
	reStart = false;
	continue;
      }

      if (criticalExit == true)
	break;

      //
      // Store the final converged eigenvectors
      //
      if (_knownEV + nFound >= _nev) {
	std::vector<int> index(1);
	for (j=0; j<_blockSize; j++) {
	    if (_normR[j] < _residual_tolerance) {
	      index[0] = j;
	      Teuchos::RefCountPtr<MV> tmp_KXvec = MVT::CloneView( *_KXvec, index );
	      index[0] = _knownEV;
	      MVT::SetBlock( *tmp_KXvec, index, *_evecs );
	      (*_evals)[_knownEV] = _theta[j];
	      _resids[_knownEV] = _normR[j];
	      _knownEV++;
	    }
	    if (_knownEV == _nev)
	      break;
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
	    Teuchos::RefCountPtr<MV> tmp_KXvec = MVT::CloneView( *_KXvec, index );
	    if (_normR[j] < _residual_tolerance) {
	      index[0] = _knownEV;
	      MVT::SetBlock( *tmp_KXvec, index, *_Xvec );	      
	      MVT::SetBlock( *tmp_KXvec, index, *_evecs );	      
	      (*_evals)[_knownEV] = _theta[j];
	      _resids[_knownEV] = _normR[j];
	      _knownEV++;
	      nFound++;	      
	    }
	    else {
	      index[0] = tmp_ptr + (j-nFound);
	      MVT::SetBlock( *tmp_KXvec, index, *_Xvec );
	    }
	  } // for (j=0; j<_blockSize; j++)
	  //
	  index.resize( nFound );
	  for (i=0; i<nFound; i++) 
	    index[i] = _knownEV + _blockSize - nFound + i;
	  Xnext = MVT::CloneView( *_Xvec, index );
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
	ScalarType tmp_swap;
	std::vector<ScalarType> tmp_swap_vec(nb*_blockSize);
	while (firstIndex < nFound) {
	  for (j=firstIndex; j<_blockSize; j++) {
	    if (_normR[j] < _residual_tolerance) {
	      //
	      // Swap and j-th and firstIndex-th position
	      //
	      blas.COPY(nb*_blockSize, _Ssdm[ j ], 1, &tmp_swap_vec[0], 1);
	      blas.COPY(nb*_blockSize, _Ssdm[ firstIndex ], 1, _Ssdm[ j ], 1 );
	      blas.COPY(nb*_blockSize, &tmp_swap_vec[0], 1, _Ssdm[ firstIndex ], 1 );
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
	  for (j=0; j<_blockSize; j++) {
	    if (_normR[j] >= _residual_tolerance) {
	      firstIndex = j;
	      break;
	    }
	  }
	} // while (firstIndex < nFound)
	//
	// Copy the converged eigenvalues and residuals.
	//
	blas.COPY( nFound, &_theta[0], 1, &(*_evals)[_knownEV], 1 ); 
	blas.COPY( nFound, &_normR[0], 1, &_resids[_knownEV], 1 );
      } // if (nFound > 0)
      //
      // Define the restarting size
      //
      bStart = ((nb - offSet) > 2 ) ? (nb - offSet)/2 : 0;
      //
      // Define the restarting space and local stiffness matrix
      //
      _KKsdm.putScalar( zero );
      for (j=0; j<bStart*_blockSize; j++)
	_KKsdm(j,j) = _theta[j + nFound];
      //
      // Form the restarting space
      //
      int oldCol = nb*_blockSize;
      int newCol = nFound + (bStart+1)*_blockSize;
      newCol = (newCol > oldCol) ? oldCol : newCol;
      std::vector<int> index( oldCol );
      lapack.GEQRF(oldCol, newCol, _Ssdm.values(), _Ssdm.stride(), &_theta[0], &_work[0], _lwork, &info);
      lapack.ORGQR(oldCol, newCol, newCol, _Ssdm.values(), _Ssdm.stride(), &_theta[0], &_work[0], _lwork, &info);      
      for (i=0; i<oldCol; i++)
	index[i] = _knownEV + i;
      Teuchos::RefCountPtr<MV> oldX = MVT::CloneView( *_Xvec, index );
      index.resize( newCol );
      for (i=0; i<newCol; i++)
	index[i] = _knownEV + i; 
      Teuchos::RefCountPtr<MV> newX = MVT::CloneView( *_Xvec, index );
      Teuchos::RefCountPtr<MV> temp_newX = MVT::Clone( *_Xvec, newCol );
      Teuchos::SerialDenseMatrix<int,ScalarType> _Sview( Teuchos::View, _Ssdm, oldCol, newCol );
      MVT::MvTimesMatAddMv( one, *oldX, _Sview, zero, *temp_newX );
      MVT::MvAddMv( one, *temp_newX, zero, *temp_newX, *newX ); 
      
      if (nFound == 0)
	offSet++;
      
      _knownEV += nFound;
      maxBlock = (_dimSearch/_blockSize) - (_knownEV/_blockSize);
      //
      // Put random vectors if the Rayleigh Ritz vectors are not enough
      // 
      newCol = nFound + (bStart+1)*_blockSize;
      if (newCol > oldCol) {
	index.resize( nFound );
	for (i=0; i<nFound; i++)
	  index[i] = _knownEV + _blockSize - nFound + i;
	Xnext = MVT::CloneView( *_Xvec, index );
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
    // Sort the eigenvectors
    //
    if ((info==0) && (_knownEV > 0))
      _MSUtils.sortScalars_Vectors(_knownEV, &(*_evals)[0], _evecs.get(), &_resids); 
    //
    // Print out a final summary if necessary
    //
    if (_om->isVerbosity( FinalSummary ))
      currentStatus();

  } // end solve()

  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::accuracyCheck(const MV *X, const MV *MX, const MV *Q) const 
    {
      cout.precision(2);
      cout.setf(ios::scientific, ios::floatfield);
      ScalarType tmp;
      
      if (X) {
	if (_BOp.get()) {
	  if (MX) {
	    tmp = _MSUtils.errorEquality(X, MX, _BOp.get());
	    if (_om->doPrint())
	      _os << " >> Difference between MX and M*X = " << tmp << endl;
	  }
	  tmp = _MSUtils.errorOrthonormality(X, _BOp.get());
	  if (_om->doPrint())
	    _os << " >> Error in X^T M X - I = " << tmp << endl;
	}
	else {
	  tmp = _MSUtils.errorOrthonormality(X, 0);
	  if (_om->doPrint())
	    _os << " >> Error in X^T X - I = " << tmp << endl;
	}
      }
      
      if (Q == 0)
	return;
      
      if (_BOp.get()) {
	tmp = _MSUtils.errorOrthonormality(Q, _BOp.get());
	if (_om->doPrint())
	  _os << " >> Error in Q^T M Q - I = " << tmp << endl;
	if (X) {
	  tmp = _MSUtils.errorOrthogonality(Q, X, _BOp.get());
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
    }
  
  } // End of namespace Anasazi

#endif

// End of file AnasaziBlockDavidson.hpp



